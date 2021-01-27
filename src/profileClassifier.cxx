// includes
#include <iostream>
#include <math.h> 
#include <istream>
#include <sstream>
#include <vtkTable.h>

#include "profileClassifier.h"
#include "corticalbone.h"
#include <vtkSortDataArray.h>

namespace itk
{
    ProfileClassifier::ProfileClassifier() {

        m_MaximumOffset = nan("1");

        m_classifierTF = ClassifierTransform::New();
        m_requireAllInsideImage = false;
        m_resampleAboutPeriostealEdge = false;
    }

    bool ProfileClassifier::SetClassifierLevel(ScalarType stDensity, ScalarType cbDensity, ScalarType threshold) {

        m_classifierTF->SetClassifierLevels(stDensity, cbDensity, threshold);

        return true;
    }

    int ProfileClassifier::GetClassifierLevelIndex() {

        return m_classifierTF->GetClassifierLevelIndex();
    }

    void ProfileClassifier::SetClassifierLevelIndex(int percentIndex) {
        m_classifierTF->SetClassifierLevelIndex(percentIndex);
    }

    std::string ProfileClassifier::GetClassifierLevelName() {

        return m_classifierTF->GetClassifierLevelName();
    }

    bool ProfileClassifier::GetClassifierThresolds(ScalarType &stDensity, ScalarType &cbDensity, ScalarType &threshold) {

        return m_classifierTF->GetClassifierThresholds(stDensity, cbDensity, threshold);
    }

    int ProfileClassifier::GetModelScheme() {
        return m_classifierTF->GetMode();
    }

    std::string ProfileClassifier::GetModelSchemeName() {
        return m_classifierTF->GetModeName();
    }

    void  ProfileClassifier::SetModelScheme(int index) {
        m_classifierTF->SetMode(index);
    }

    void ProfileClassifier::ClearThresholds() {
        m_classifierTF->ClearThresholds();
    }

    bool ProfileClassifier::Process() { // process if m_startIndex is before the edgeIndex

        if(!m_profile->IsProfileAveragingOn()) {
            cerr<<"Cannont perform profile thresholding without multiple profiles"<<endl;
            InitaliseProfileArrays();
            return false;
        }

        // set sampled values
        m_classifierTF->SetProfilePositions(m_profile->GetPositions());
        m_classifierTF->SetProfileValues(m_profile->GetValues());
        m_classifierTF->SetMultipleProfileValues(m_profile->GetMultipleValues());

        // set sampled properties
        int startIndex, stopIndex, maxIndex; ScalarType maxValue;
        itkArrayType maxValues = m_profile->GetMaxValues(); itkIntArrayType maxIndices = m_profile->GetMaxIndices();
        m_profile->GetSampledExtentIndices(startIndex, stopIndex);
        m_profile->GetSampledMax(maxIndex, maxValue);
        m_classifierTF->SetProfileMaxPositionInfo(maxIndex, maxValue, maxValues, maxIndices);
        m_classifierTF->SetProfilePositionInfo(startIndex, stopIndex, m_profile->GetNumberOfProfiles(), m_profile->GetStartXValue());
        m_classifierTF->SetProfileSampleInfo(m_profile->GetIncrement(), m_profile->GetNumberOfSamples(), m_profile->GetProfileLength(), m_profile->GetEdgeIndex());

        bool runStatus = m_classifierTF->Threshold();
        if(runStatus) {
            m_ProcessedTimeStamp.Modified();
        }
        return m_classifierTF->IsUpToDate();

    }

    bool ProfileClassifier::InitialiseDensities() { // only estimate if entirely in bounds

        if(!IsPtInsideImage(true)) { // first check to see if entirely inbounds
            return false;
        }

        // set sampled values
        m_classifierTF->SetProfilePositions(m_profile->GetPositions());
        m_classifierTF->SetProfileValues(m_profile->GetValues());
        m_classifierTF->SetMultipleProfileValues(m_profile->GetMultipleValues());

        // set sampled properties
        int startIndex, stopIndex, maxIndex; ScalarType maxValue;
        itkArrayType maxValues = m_profile->GetMaxValues(); itkIntArrayType maxIndices = m_profile->GetMaxIndices();
        m_profile->GetSampledExtentIndices(startIndex, stopIndex);
        m_profile->GetSampledMax(maxIndex, maxValue);
        m_classifierTF->SetProfileMaxPositionInfo(maxIndex, maxValue, maxValues, maxIndices);
        m_classifierTF->SetProfilePositionInfo(startIndex, stopIndex, m_profile->GetNumberOfProfiles(), m_profile->GetStartXValue());
        m_classifierTF->SetProfileSampleInfo(m_profile->GetIncrement(), m_profile->GetNumberOfSamples(), m_profile->GetProfileLength(), m_profile->GetEdgeIndex());

        bool runStatus = m_classifierTF->EstimateThresholds();

        return runStatus;
    }

    void ProfileClassifier::TurnOnProfileResampling(ScalarType maximum) {
        if(isnan(maximum) && isinf(maximum)) {
            m_resampleAboutPeriostealEdge = false;
            m_MaximumOffset = nan("1");
            cerr<<"Error in ProfileClassifier::TurnOnProfileResampling invalid maximum value"<<endl;
        } else {
            m_resampleAboutPeriostealEdge = true;
            m_MaximumOffset = maximum;
        }
        m_classifierTF->SetProfileRangeLimits(m_MaximumOffset);
    }

    void ProfileClassifier::TurnOffProfileResampling() {
        m_resampleAboutPeriostealEdge = false;
        m_MaximumOffset = nan("1");
        m_classifierTF->SetProfileRangeLimits(m_MaximumOffset);
    }

    bool ProfileClassifier::GetResamplingMode() {
        return m_resampleAboutPeriostealEdge;
    }

    bool ProfileClassifier::IsValid() const {
        return m_classifierTF->IsValid();

    }

    bool ProfileClassifier::IsNan() const {
        return false;
    }

    bool ProfileClassifier::IsUpToDate() const {
        return m_classifierTF->IsUpToDate();

    }

    bool ProfileClassifier::AreProfilesSufficientlyInsideImage() const {

        LinearTransform::ProfilePropertiesType profileProperties = LinearTransform::GetProfileProperties();

        return profileProperties.startIndex != -1 && profileProperties.startIndex < profileProperties.edgeIndex && profileProperties.endIndex == profileProperties.numberOfSamples-1;
    }

    bool ProfileClassifier::AreThresholdsSet() const {
        return m_classifierTF->AreThresholdsSet();

    }

    // private
    void ProfileClassifier::InitaliseProfileArrays() {

        // reset the mean profile display array
        itkArrayType values = itkArrayType(1); values.Fill(nan("1"));
        itkArray2DType multipleValues = itkArray2DType(1,1); values.Fill(nan("1"));
        itkArrayType positions = itkArrayType(1); positions.Fill(nan("1"));

        m_classifierTF->SetProfileValues(values);
        m_classifierTF->SetMultipleProfileValues(multipleValues);
        m_classifierTF->SetProfilePositions(positions);
    }

    //void ProfileThresholder::GenerateData() {    Process(); }

    bool ProfileClassifier::ReadyToProcess(bool thresholdsRequired) {

        return m_classifierTF->ReadyToProcess(thresholdsRequired);
    }

    // Public Methods - Getters
    ProfileClassifier::itkArrayType ProfileClassifier::GetParameterValues() const {
        return m_classifierTF->GetParameterValues();

    }

    ProfileClassifier::itkArrayType ProfileClassifier::GetDisplayValues() const {
        return m_classifierTF->GetDisplayValues();

    }

    ProfileClassifier::TableType ProfileClassifier::GetImageDataTable() const {

        TableType table;
        if(IsPtInsideImage()) {
            table = m_classifierTF->GetProfileDisplayTable();
        } else {
            table = TableType::New();
            table->SetNumberOfRows(m_profile->GetNumberOfSamples());
        }
        return table;
    }

    ProfileClassifier::TableType ProfileClassifier::GetDisplayModelTable() const {
        return m_classifierTF->GetModelDisplayTable();

    }

    ProfileClassifier::TableType ProfileClassifier::GetDisplayModelTable(vtkSmartPointer<vtkDoubleArray> parametersIn, double offset) const {
        return m_classifierTF->GetModelDisplayTable(parametersIn, offset);
    }

    ProfileClassifier::TableType ProfileClassifier::GetProcessingModelTable() const {
        return m_classifierTF->GetClassifierTable();
    }

    unsigned int ProfileClassifier::GetNumberOfParameterValues() const {
        return m_classifierTF->GetNumberOfParameters();

    }

    std::string ProfileClassifier::GetParameterName(int index) const {
        return m_classifierTF->GetParameterName(index);

    }

    std::string ProfileClassifier::GetParameterNameShort(int index) const {
        return m_classifierTF->GetParameterNameShort(index);

    }

    unsigned int ProfileClassifier::GetNumberOfDisplayValues(bool includeImportedOptions) const {
        return m_classifierTF->GetNumberOfDisplayValues(includeImportedOptions);

    }

    std::string ProfileClassifier::GetDisplayName(int index) const {
        return m_classifierTF->GetDisplayName(index);

    }

    std::string ProfileClassifier::GetDisplayNameShort(int index) const {
        return m_classifierTF->GetDisplayNameShort(index);

    }

    ProfileClassifier::itkArray2DType ProfileClassifier::GetDisplayRanges(bool includeImportedOptions) const {
        return m_classifierTF->GetDisplayRanges(includeImportedOptions);

    }

    ProfileClassifier::itkArrayType ProfileClassifier::GetClassifications() const {
        return m_classifierTF->GetClassifications();
    }

    ProfileClassifier::itkArrayType ProfileClassifier::GetPercentages() const {
        return m_classifierTF->GetPercentages();
    }

    // getters for display values
    ProfileClassifier::ScalarType ProfileClassifier::GetErrorMean() const {
        return m_classifierTF->GetErrorMean();

    }

}