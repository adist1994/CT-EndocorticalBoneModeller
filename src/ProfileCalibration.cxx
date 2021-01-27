//
// Created by rap58 on 20/10/15.
//

#include "ProfileCalibration.h"
#include "utilities.h"
#include "corticalbone.h"


namespace itk
{
    ProfileCalibration::ProfileCalibration() {

        m_lowMedian = m_highMedian = nan("1");
        m_modeIndex = 0;
        m_requireAllInsideImage = true;
        m_ProcessedTimeStamp.Modified();
    }

    bool ProfileCalibration::Process() { // process if m_startIndex is before the edgeIndex

        bool status;
        if(!IsPtInsideImage()) {
            // ensure the correct display - i.e. up-to-date
            m_lowMedian = m_highMedian = nan("1");
            status = false;
        } else {
            status = CalculateCalibrationValues();

        }
        m_ProcessedTimeStamp.Modified();
        return status;
    }

    bool ProfileCalibration::IsValid() const {

        return !isnan(m_lowMedian) && !isnan(m_highMedian) && m_modeIndex>=kMedian && m_modeIndex<=kMaximum;

    }

    bool ProfileCalibration::IsUpToDate() const {
        if(m_ProcessedTimeStamp.GetMTime() < m_profile->GetSizeModificationTime()) {
            return false;
        } else if(m_ProcessedTimeStamp.GetMTime() < m_profile->GetPositionModificationTime()) {
            return false;
        } else {
            return true;
        }
    }

    bool ProfileCalibration::AreProfilesSufficientlyInsideImage() const {
        return m_profile->GetStartIndex()==0 && IsPtInsideImage(); // start must be inside
    }

    // private
    bool ProfileCalibration::CalculateCalibrationValues() {

        int numberOfSamples = m_profile->GetNumberOfSamples();

        int startIndex, endIndex, maxIndex; ScalarType maxValue;
        m_profile->GetSampledExtentIndices(startIndex, endIndex);
        m_profile->GetSampledMax(maxIndex, maxValue);

        // calculate the indicies  for the low (1st 20%) and high median
        ScalarType lowRatio = 0.2, highRatio = 0.1;
        int nLow = (int)round(numberOfSamples*lowRatio);
        int nHigh = round(numberOfSamples*highRatio/2.0) * 2; // rounded to the nearest 2

        // check for valid values
        if(startIndex!=0 || endIndex!=numberOfSamples-1 || !IsPtInsideImage() || maxIndex<nLow || maxIndex+nHigh>numberOfSamples-1 || m_modeIndex<kMedian || m_modeIndex>kMaximum) {
            m_highMedian = m_lowMedian = nan("1");
            return false;
        }

        itkArrayType meanValues = m_profile->GetValues();

        // calculate the medians
        m_lowMedian = m_linearRegressionCalculator.CalculateArrayMedian(meanValues, 0, nLow);
        if(m_modeIndex==kMedian) {
            m_highMedian = m_linearRegressionCalculator.CalculateArrayMedian(meanValues, maxIndex - nHigh / 2, maxIndex + nHigh / 2);
        } else if(m_modeIndex==kMaximum) {
            m_highMedian = maxValue;
        }
        return true;

    }

    //void ProfileThresholder::GenerateData() {    Process(); }

    // Public Methods - Getters
    ProfileCalibration::itkArrayType ProfileCalibration::GetParameterValues() const {

        itkArrayType parameters = itkArrayType(m_numberOfParameters); // name = "Transform Parameters"
        parameters[0] = m_lowMedian; parameters[1] = m_highMedian;
        return parameters;

    }

    ProfileCalibration::itkArrayType ProfileCalibration::GetDisplayValues() const {

        // if(!IsUpToDate()) { Process(); } // non-const

        itkArrayType displayValues(m_numberOfParameters);

        if(IsUpToDate()) {
            displayValues[0] = m_lowMedian;
            displayValues[1] = m_highMedian;
        } else {
            displayValues[0] = nan("1");
            displayValues[1] = nan("1");
        }


        return displayValues;

    }

    ProfileCalibration::TableType ProfileCalibration::GetImageDataTable() const {

        int numberOfSamples = m_profile->GetNumberOfSamples();
        TableType table = TableType::New();
        table->SetNumberOfRows(numberOfSamples);

        if(IsPtInsideImage()) {
            ArrayType positions = Utilities::convertItkToVtkArray(m_profile->GetPositions(), "Positions");
            ArrayType values = Utilities::convertItkToVtkArray(m_profile->GetValues(), "Profile Mean");

            table->AddColumn(positions);
            table->AddColumn(values); // the average

        } else {
            // do nothing
        }
        return table;
    }

    ProfileCalibration::TableType ProfileCalibration::GetDisplayModelTable() const {

        TableType table = TableType::New();

        table->SetNumberOfRows(m_numberOfParameters);

        ScalarType length= m_profile->GetProfileLength();
        ScalarType startPosition = m_profile->GetStartXValue();

        // create and add the low median value line
        if(!isnan(m_lowMedian)) {
            ArrayType lowPositions = ArrayType::New();
            lowPositions->SetName("Low Median Positions");
            lowPositions->SetNumberOfValues(2);
            lowPositions->SetValue(0, startPosition);
            lowPositions->SetValue(1, startPosition + length);
            table->AddColumn(lowPositions);

            ArrayType lowValues = ArrayType::New();
            lowValues->SetName("Low Median");
            lowValues->SetNumberOfValues(2);
            lowValues->SetValue(0, m_lowMedian);
            lowValues->SetValue(1, m_lowMedian);
            table->AddColumn(lowValues);
        }

        // create and add the low median value line
        if(!isnan(m_lowMedian)) {
            ArrayType highPositions = ArrayType::New();
            highPositions->SetName("High Median Positions");
            highPositions->SetNumberOfValues(2);
            highPositions->SetValue(0, startPosition);
            highPositions->SetValue(1, startPosition + length);
            table->AddColumn(highPositions);

            ArrayType highValues = ArrayType::New();
            highValues->SetName("High Median");
            highValues->SetNumberOfValues(2);
            highValues->SetValue(0, m_highMedian);
            highValues->SetValue(1, m_highMedian);
            table->AddColumn(highValues);
        }

        return table;

    }

    ProfileCalibration::TableType ProfileCalibration::GetProcessingModelTable() const {
        TableType table = TableType::New();
        return table;
    }

    unsigned int ProfileCalibration::GetNumberOfParameterValues() const {
        return m_numberOfParameters;
    }

    std::string ProfileCalibration::GetParameterName(int index) const {

        if(index==kLowMedian) {
            return "Low Median";
        } else if(index==kHighMedian) {
            return "High Median";
        } else if(index==itk::CorticalBone::kInvalid) {
            return "N/A";
        } else {
            std::cerr<<"Error: ProfileCalibration::GetParameterName() invalid index:"<<index<<endl;
            return "Invalid Index";
        }

    }

    std::string ProfileCalibration::GetParameterNameShort(int index) const {

        if(index==kLowMedian) {
            return "Low";
        } else if(index==kHighMedian) {
            return "High";
        } else if(index==itk::CorticalBone::kInvalid) {
            return "N/A";
        } else {
            std::cerr<<"Error: ProfileCalibration::GetParameterNameShort() invalid index:"<<index<<endl;
            return "Invalid Index";
        }
    }

    unsigned int ProfileCalibration::GetNumberOfDisplayValues(bool includeImportedOptions) const {

        return GetNumberOfParameterValues();

    }

    std::string ProfileCalibration::GetDisplayName(int index) const {

        return GetParameterName(index);

    }

    std::string ProfileCalibration::GetDisplayNameShort(int index) const {

        return GetParameterNameShort(index);
    }

    ProfileCalibration::itkArray2DType ProfileCalibration::GetDisplayRanges(bool includeImportedOptions) const {

        itkArray2DType displayRanges(m_numberOfParameters, 2);

        displayRanges.put(kLowMedian,0,-200);     displayRanges.put(kLowMedian,1,200);
        displayRanges.put(kHighMedian,0,600.0);            displayRanges.put(kHighMedian,1,1600);



        return displayRanges;

    }

    int ProfileCalibration::GetModelScheme() {
        return m_modeIndex;
    }



    // Public Methods - Setters
    void  ProfileCalibration::SetModelScheme(int index) {
        if(m_modeIndex>=kMedian && m_modeIndex<=kMaximum) {
            m_modeIndex=index;
        } else {
            m_modeIndex=-1;
        }
    }


}