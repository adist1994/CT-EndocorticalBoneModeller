#include <iostream>
#include <math.h> 
#include <istream>
#include <sstream>
#include <vtkTable.h>

#include "linearTransform.h"
#include "profileClassifier.h"
#include <vtkSortDataArray.h>

namespace itk {

    // private static variable declaration
    LinearRegressionCalculator LinearTransform::m_linearRegressionCalculator = LinearRegressionCalculator();

    // Setup
    LinearTransform::LinearTransform() {

        m_positions = m_meanProfiles = itkArrayType(1); m_positions.Fill(nan("1")); m_meanProfiles.Fill(nan("1"));
        m_multipleProfiles = itkArray2DType(1,1); m_multipleProfiles.Fill(nan("1"));

        ClearThresholds();

        m_modeIndex = 0;

        m_valid = false;

        m_ProcessedTimeStamp.Modified();
    }

    void LinearTransform::SetProfilePositionInfo(int startIndex, int endIndex, int numberOfProfiles, double startPosition) {
        m_linearRegressionCalculator.SetProfilePositionInfo(startIndex, endIndex, numberOfProfiles, startPosition);
    }

    void LinearTransform::SetProfileMaxPositionInfo(int maxIndex, ScalarType maxValue, itkArrayType maxValues, itkIntArrayType maxIndices) {
        m_linearRegressionCalculator.SetProfileMaxPositionInfo(maxIndex, maxValue, maxValues, maxIndices);
    }

    void LinearTransform::SetProfileSampleInfo(ScalarType increment, int numberOfSamples, ScalarType length, int edgeIndex) {
        m_linearRegressionCalculator.SetProfileSampleInfo(increment, numberOfSamples, length, edgeIndex);
    }

    void LinearTransform::SetProfilePositions(itkArrayType positions) {

        if(m_positions != positions) {
            m_positions = positions; // todo - error check and enure the array is the correct size
            m_transformUpToDate = false;
        }
    }

    void LinearTransform::SetProfileValues(itkArrayType values) {

        if(m_meanProfiles != values) {
            m_meanProfiles = values;
            m_transformUpToDate = false;
        }
    }

    void LinearTransform::SetMultipleProfileValues(itkArray2DType multipleValues) {
        if(m_multipleProfiles != multipleValues) {
            m_multipleProfiles = multipleValues; m_transformUpToDate = false;
        }
    }

    void LinearTransform::SetClassifierLevels(ScalarType stDensity, ScalarType cbDensity, ScalarType threshold) {

        if(m_thresholdIndex <kVLowPercent || m_thresholdIndex > kMedianManual) {
            cerr<<"Warning in LinearTransform::SetClassifierLevels invalid percentID = "<< m_thresholdIndex <<endl;
            ClearThresholds();
        } else if(stDensity != m_ST || cbDensity != m_CB || threshold != m_threshold) {

            m_ST = stDensity; m_CB = cbDensity; m_threshold = threshold;
            m_transformUpToDate = false; m_thresholdsSet = true;
        }
    }

    void LinearTransform::ClearThresholds() {
        ClearThresholdValues();

        m_thresholdPValue = m_thresholdPercentage = nan("1"); //m_thresholdIndex = -1;
    }

    void LinearTransform::ClearThresholdValues() {
        m_ST = m_CB = m_threshold = nan("1");
        m_thresholdsSet = false; m_transformUpToDate = false;
    }

    void LinearTransform::SetMode(int index) {
        m_modeIndex = index;
    }

    int  LinearTransform::GetMode() {
        return m_modeIndex;
    }

    // get state
    bool LinearTransform::IsValid() const {
        return m_valid;
    }

    bool LinearTransform::IsUpToDate() const {
        return m_transformUpToDate;
    }

    bool LinearTransform::AreThresholdsSet() const {
        return m_thresholdsSet;
    }

    bool LinearTransform::EstimateThresholds() {

        ClearResults(); ClearThresholdValues(); // thresholdsSet = false; in clearThresholdValues

        // check if enough info to find edges
        ScalarType minRatio = 0.15; int minN = int(minRatio * m_linearRegressionCalculator.GetProfileSampleNumber());

        if(m_thresholdIndex < kVLowPercent && m_thresholdIndex > kHighPercent) {
            cerr<<"Error in LinearTransform::EstimateThresholds() invalid percent ID of "<< m_thresholdIndex <<endl;
        } else if (m_linearRegressionCalculator.GetProfileStartIndex() == 0 && minN < m_linearRegressionCalculator.GetProfileMaxIndex()) { // stricter 'in-bounds' requirements as this will be averaged between all mesh vertices
            m_STLine = m_linearRegressionCalculator.CalculateLinearRegressionMultiple(m_positions, m_multipleProfiles, 0, 0+minN);
            m_ST = m_STLine.meanY;
            m_CB = m_linearRegressionCalculator.CalculateArrayMedian(m_linearRegressionCalculator.GetProfileMaxValues(), 0, m_linearRegressionCalculator.GetProfileNumber()-1);

            m_threshold =  m_STLine.meanY + m_thresholdPValue * m_STLine.stdY;
            m_thresholdsSet = !isnan(m_threshold); // just in case only one

        } else { // not enough info to find thresholds

            //cerr<<"Warning in ClassifierTransform::GetThreshold() profile start or end is out of bounds"<<endl;
        }
        return m_thresholdsSet;
    }

    // Get values
    unsigned int LinearTransform::GetNumberOfParameters() const {
        return m_numberOfParameters;
    }

    unsigned int LinearTransform::GetNumberOfDisplayValues(bool importSet) const {
        if(importSet) {
            return m_numberOfDisplayValues;
        } else {
            return m_numberOfDisplayValues;
        }
    }

    unsigned int LinearTransform::GetNumberOfSamples() const {
        int numberOfSamples = m_linearRegressionCalculator.GetProfileSampleNumber();
        if(numberOfSamples==-1) {
            throw std::logic_error( "Error: ProfileProperties not yet set in the LinearTransform class" );
        }

        return numberOfSamples;
    }

    LinearTransform::ProfilePropertiesType LinearTransform::GetProfileProperties() {
        return m_linearRegressionCalculator.GetProfileProperties();
    }

    bool LinearTransform::GetClassifierThresholds(ScalarType &stDensity, ScalarType &cbDensity, ScalarType &threshold) {

        if(m_thresholdsSet) {
            stDensity = m_ST; cbDensity = m_CB; threshold = m_threshold;
        }
        return m_thresholdsSet;
    }

    int LinearTransform::GetClassifierLevelIndex() {
        return m_thresholdIndex; // should be set to -1 if not calculated
    }

    void LinearTransform::SetClassifierLevelIndex(int thresholdIndex) {
        if(thresholdIndex<kVLowPercent || thresholdIndex > kMedianManual) {
            cerr<<"Warning in LinearTransform::SetClassifierLevelsIndex invalid percentID = "<<thresholdIndex<<endl;
            ClearThresholds();
        } else if(thresholdIndex != m_thresholdIndex && thresholdIndex <= kMedianManual) {
            ClearThresholds();
            m_thresholdIndex = thresholdIndex;
            m_thresholdPercentage = LookupPercentDensity(m_thresholdIndex);
            m_thresholdPValue = LookupPValue(m_thresholdIndex);

        }
    }

    std::string LinearTransform::GetClassifierLevelName() {
        if(m_thresholdIndex ==kVLowPercent) {
            return "Very High Threshold [99.99%]";
        } else if(m_thresholdIndex ==kLowPercent) {
            return  "High Threshold [99.9%]";
        } else if(m_thresholdIndex ==kMidPercent) {
            return "Medium Threshold [99.0%]";
        } else if(m_thresholdIndex ==kHighPercent) {
            return "Low Threshold [95.0%]";
        } else if(m_thresholdIndex == kManualThreshold) {
            return "Manually Set"; // todo decide value to use for manaual
        } else if(m_thresholdIndex == kMedianMidpoint) {
            return "Median Midpoint"; // todo decide value to use for manaual
        } else if(m_thresholdIndex == kMedianMidpoint) {
            return "Median Manual"; // todo decide value to use for manaual
        } else {
            return "Invalid Threshold; index="+std::to_string(m_thresholdIndex);
        }
    }

    LinearTransform::itkArrayType LinearTransform::GetProfilePositions() const {
        return m_positions;
    }

    LinearTransform::itkArrayType LinearTransform::GetProfileValues() const {
        return m_meanProfiles;
    }

    LinearTransform::itkArray2DType LinearTransform::GetMultipleProfileValues() const {
        return m_multipleProfiles;
    }

    LinearTransform::ScalarType LinearTransform::GetErrorMean() const {
        return CalculateAbsMeanError();
    }

    //----------- private ---------------//
    bool LinearTransform::ReadyToProcess(bool thresholdsRequired) const {

        int numberOfSamples = m_linearRegressionCalculator.GetProfileSampleNumber();
        int startIndex = m_linearRegressionCalculator.GetProfileStartIndex();
        int endIndex = m_linearRegressionCalculator.GetProfileEndIndex();

        if(m_positions.Size() != numberOfSamples || m_meanProfiles.Size() != numberOfSamples || m_multipleProfiles.rows() != numberOfSamples) {
            return false;
        } else if(startIndex == -1 || m_linearRegressionCalculator.GetProfileMaxIndex() == -1) {
            return false;
        } else if( startIndex >= m_linearRegressionCalculator.GetProfileEdgeIndex() || endIndex != numberOfSamples-1) {
            return false;
        } else if(thresholdsRequired && (m_thresholdIndex >= kManualThreshold || m_thresholdIndex <= kMedianManual)) {
            return m_thresholdsSet && !isnan(m_ST) && !isnan(m_CB) && !isnan(m_threshold);
        } else if(thresholdsRequired && m_thresholdIndex < kManualThreshold && m_thresholdIndex > kMedianManual) {
            return m_thresholdsSet && !isnan(m_ST) && !isnan(m_CB) && !isnan(m_threshold) && !isnan(m_thresholdPercentage);
        } else {
            return true;
        }
    }

    LinearTransform::ScalarType LinearTransform::LookupPercentDensity(int index) {
        if(index==kVLowPercent) {
            return 0.01;
        } else if(index==kLowPercent) {
            return 0.1;
        } else if(index==kMidPercent) {
            return 1.0;
        } else if(index==kHighPercent) {
            return 5.0;
        } else if(index == kManualThreshold) {
            return nan("1"); // todo decide value to use for manaual
        } else if(index == kMedianMidpoint) {
            return nan("1"); // todo decide value to use for manaual
        } else if(index == kMedianManual) {
            return nan("1"); // todo decide value to use for manaual
        } else {
            return nan("1");
        }
    }

    LinearTransform::ScalarType LinearTransform::LookupPValue(int index) {
        if(index==kVLowPercent) {
            return 3.62;
        } else if(index==kLowPercent) {
            return 3.08;
        } else if(index==kMidPercent) {
            return 2.33;
        } else if(index==kHighPercent) {
            return 1.645;
        } else if(index == kManualThreshold) {
            return nan("1");
        } else if(index == kManualThreshold) {
            return nan("1");
        } else {
            return nan("1");
        }
    }

    LinearTransform::ScalarType LinearTransform::LookupPercentThreshold(int index) {
        if(index==kVLowPercent) { // actually v high in gui
            return 5;
        } else if(index==kLowPercent) {
            return 10;
        } else if(index==kMidPercent) {
            return 15;
        } else if(index==kHighPercent) {
            return 20;
        } else if(index == kManualThreshold) {
            return nan("1");
        } else if(index == kMedianMidpoint) {
            return nan("1");
        } else if(index == kMedianManual) {
            return nan("1");
        } else {
            return nan("1");
        }
    }

}