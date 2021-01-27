#include <iostream>
#include <math.h> 
#include <istream>
#include <sstream>

#include "classifierTransform.h"
#include "profileClassifier.h"
#include "utilities.h"
#include <vtkSortDataArray.h>

namespace itk {
    // Setup
    ClassifierTransform::ClassifierTransform() : LinearTransform() {

        //m_profileDetails = ProfileProperties();
        m_numberOfParameters = kPorosityParameter +1;
        m_numberOfDisplayValues = kProfileNumber +1;

        m_classificationType = -1;

        m_profileRange = nan("1");

        ClearResults();
    }

    // Process events
    bool ClassifierTransform::Threshold() {

        if(!ReadyToProcess() || isnan(m_profileRange)) {
            return false;
        }

        if(!IsUpToDate()) {

            ClearResults();

            //ClassifyGlobalSingleThreshold();//Classifing();
            PercentageBasedClassifing();

            m_ProcessedTimeStamp.Modified();

            CalculateValidity();

            m_transformUpToDate = true;
        }

        return m_valid;

    }

    void ClassifierTransform::SetProfileRangeLimits(ScalarType maximum) {
        m_profileRange = maximum;
    }

    // return mode name
    std::string  ClassifierTransform::GetModeName() {

        if(m_modeIndex==kGlobal) {
            return "GlobalScheme";
        } else if(m_modeIndex==kMedian) {
            return "MedianScheme";
        } else {
            return "Undefined";
        }
    }

    // Get parameter values and names
    ClassifierTransform::itkArrayType ClassifierTransform::GetParameterValues() const { // returns curve results - GetTissueCurves
        itkArrayType parameters = itkArrayType(m_numberOfParameters); // name - "Transform Parameters"
        double increment = m_linearRegressionCalculator.GetProfileIncrement();
        double profileStartLocation = m_linearRegressionCalculator.GetProfileStartPosition();
        

        if(m_transformUpToDate) {
            parameters[kSTMean] = m_STLine.meanY;
            parameters[kCBMean] = m_CBLine.meanY;
            parameters[kDBMean] = m_DBLine.meanY;
            parameters[kTBMean] = m_TBLine.meanY;
            parameters[kECMean] = m_ECLine.meanY;
            parameters[kTotMean] = m_totLine.meanY;

            parameters[kSTSlope] = m_STLine.slope;
            parameters[kCBSlope] = m_CBLine.slope;
            parameters[kDBSlope] = m_DBLine.slope;
            parameters[kTBSlope] = m_TBLine.slope;
            parameters[kECSlope] = m_ECLine.slope;
            parameters[kTotSlope] = m_totLine.slope;

            parameters[kSTIntercept] = m_STLine.intercept;
            parameters[kCBIntercept] = m_CBLine.intercept;
            parameters[kDBIntercept] = m_DBLine.intercept;
            parameters[kTBIntercept] = m_TBLine.intercept;
            parameters[kECIntercept] = m_ECLine.intercept;
            parameters[kTotIntercept] = m_totLine.intercept;

            parameters[kSTSTD] = m_STLine.stdY;
            parameters[kCBSTD] = m_CBLine.stdY;
            parameters[kDBSTD] = m_DBLine.stdY;
            parameters[kTBSTD] = m_TBLine.stdY;
            parameters[kECSTD] = m_ECLine.stdY;
            parameters[kTotSTD] = m_totLine.stdY;

            parameters[kSTStart] = profileStartLocation + m_stIndices[0] * increment;
            parameters[kCBStart] = profileStartLocation + m_cbIndices[0] * increment;
            parameters[kDBStart] = profileStartLocation + m_dbIndices[0] * increment;
            parameters[kTBStart] = profileStartLocation + m_tbIndices[0] * increment;
            parameters[kECStart] = profileStartLocation + m_ecIndices[0] * increment;
            parameters[kTotStart] = profileStartLocation + m_totIndices[0] * increment;

            parameters[kSTEnd] = profileStartLocation + m_stIndices[1] * increment;
            parameters[kCBEnd] = profileStartLocation + m_cbIndices[1] * increment;
            parameters[kDBEnd] = profileStartLocation + m_dbIndices[1] * increment;
            parameters[kTBEnd] = profileStartLocation + m_tbIndices[1] * increment;
            parameters[kECEnd] = profileStartLocation + m_ecIndices[1] * increment;
            parameters[kTotEnd] = profileStartLocation + m_totIndices[1] * increment;

            parameters[kPorosityParameter] = m_porosity;
        } else {
            parameters.Fill(nan("1"));
        }

        return parameters;
    }

    std::string ClassifierTransform::GetParameterName(int index) const {

        if(index == kSTMean) {
            return "Soft Tissue Mean";
        } else if(index == kCBMean) {
            return "Cortical Bone Mean";
        } else if(index == kDBMean) {
            return "Dense Bone Mean";
        } else if(index == kTBMean) {
            return "Trabecular Bone Mean";
        } else if(index == kECMean) {
            return "Endocortical Region Mean";
        } else if(index == kTotMean) {
            return "Total Cortical Region Mean";
        } else if(index == kSTSlope) {
            return "Soft Tissue Slope";
        } else if(index == kCBSlope) {
            return "Cortical Bone Slope";
        } else if(index == kDBSlope) {
            return "Dense Bone Slope";
        } else if(index == kTBSlope) {
            return "Trabecular Bone Slope";
        } else if(index == kECSlope) {
            return "Endocortical Region Slope";
        } else if(index == kTotSlope) {
            return "Total Cortical Region Slope";
        } else if(index == kSTIntercept) {
            return "Soft Tissue Intercept";
        } else if(index == kCBIntercept) {
            return "Cortical Bone Intercept";
        } else if(index == kDBIntercept) {
            return "Dense Bone Intercept";
        } else if(index == kTBIntercept) {
            return "Trabecular Bone Intercept";
        } else if(index == kECIntercept) {
            return "Endocortical Region Intercept";
        } else if(index == kTotIntercept) {
            return "Total Cortex Region Intercept";
        } else if(index == kSTSTD) {
            return "Soft Tissue Standard Deviation";
        } else if(index == kCBSTD) {
            return "Cortical Bone Standard Deviation";
        } else if(index == kDBSTD) {
            return "Dense Bone Standard Deviation";
        } else if(index == kTBSTD) {
            return "Trabecular Bone Standard Deviation";
        } else if(index == kECSTD) {
            return "Endocortical Region Standard Deviation";
        } else if(index == kTotSTD) {
            return "Total Cortex Region Standard Deviation";
        } else if(index == kSTStart) {
            return "Soft Tissue Start";
        } else if(index == kSTEnd) {
            return "Soft Tissue End";
        } else if(index == kCBStart) {
            return "Cortical Bone Tissue Start";
        } else if(index == kCBEnd) {
            return "Cortical Bone End";
        } else if(index == kDBStart) {
            return "Dense Bone Tissue Start";
        } else if(index == kDBEnd) {
            return "Dense Bone End";
        } else if(index == kTBStart) {
            return "Trabecular Bone Tissue Start";
        } else if(index == kTBEnd) {
            return "Trabecular Bone End";
        } else if(index == kECStart) {
            return "Endocortical Region Start";
        } else if(index == kECEnd) {
            return "Endocortical Region End";
        } else if(index == kTotStart) {
            return "Total Cortex Region Start";
        } else if(index == kTotEnd) {
            return "Total Cortex Region End";
        } else if(index == kPorosityParameter) {
            return "Porosity";
        } else {
            std::cout<<"Error: ClassifierTransform::GetParameterName() invalid GetParameterName index.";
            return "Invalid Selection";
        }
    }

    std::string ClassifierTransform::GetParameterNameShort(int index) const {

        if(index == kSTMean) {
            return "\u03BC<sub>ST</sub>";
        } else if(index == kCBMean) {
            return "\u03BC<sub>CB</sub>";
        } else if(index == kDBMean) {
            return "\u03BC<sub>DB</sub>";
        } else if(index == kTBMean) {
            return "\u03BC<sub>TB</sub>";
        } else if(index == kECMean) {
            return "\u03BC<sub>EC</sub>";
        } else if(index == kTotMean) {
            return "\u03BC<sub>total</sub>";
        } else if(index == kSTSlope) {
            return "m<sub>ST</sub>";
        } else if(index == kCBSlope) {
            return "m<sub>CB</sub>";
        } else if(index == kDBSlope) {
            return "m<sub>DB</sub>";
        } else if(index == kTBSlope) {
            return "m<sub>TB</sub>";
        } else if(index == kECSlope) {
            return "m<sub>EC</sub>";
        } else if(index == kTotSlope) {
            return "m<sub>Total</sub>";
        } else if(index == kSTIntercept) {
            return "b<sub>ST</sub>";
        } else if(index == kCBIntercept) {
            return "b<sub>CB</sub>";
        } else if(index == kDBIntercept) {
            return "b<sub>DB</sub>";
        } else if(index == kTBIntercept) {
            return "b<sub>TB</sub>";
        } else if(index == kECIntercept) {
            return "b<sub>EC</sub>";
        } else if(index == kTotIntercept) {
            return "b<sub>Total</sub>";
        } else if(index == kSTSTD) {
            return "\u03C3<sub>ST</sub>";
        } else if(index == kCBSTD) {
            return "\u03C3<sub>CB</sub>";
        } else if(index == kDBSTD) {
            return "\u03C3<sub>DB</sub>";
        } else if(index == kTBSTD) {
            return "\u03C3<sub>TB</sub>";
        } else if(index == kECSTD) {
            return "\u03C3<sub>MC</sub>";
        } else if(index == kTotSTD) {
            return "\u03C3<sub>Total</sub>";
        } else if(index == kSTStart) {
            return "X<sub>STStart</sub>";
        } else if(index == kCBStart) {
            return "X<sub>CBStart</sub>";
        } else if(index == kDBStart) {
            return "X<sub>DBStart</sub>";
        } else if(index == kTBStart) {
            return "X<sub>TBStart</sub>";
        } else if(index == kECStart) {
            return "X<sub>ECStart</sub>";
        } else if(index == kTotStart) {
            return "X<sub>TotalStart</sub>";
        } else if(index == kSTEnd) {
            return "X<sub>STEnd</sub>";
        } else if(index == kCBEnd) {
            return "X<sub>CBEnd</sub>";
        } else if(index == kDBEnd) {
            return "X<sub>DBEnd</sub>";
        } else if(index == kTBEnd) {
            return "X<sub>TBEnd</sub>";
        } else if(index == kECEnd) {
            return "X<sub>ECEnd</sub>";
        } else if(index == kTotEnd) {
            return "X<sub>TotalEnd</sub>";
        } else if(index == kPorosityParameter) {
            return "\u03D5";
        } else {
            std::cout<<"Error: ClassifierTransform::GetParameterNameShort() invalid GetParameterNameShort index.";
            return "Invalid Selection";
        }
    }

    // Get display values and names
    ClassifierTransform::itkArrayType ClassifierTransform::GetDisplayValues() const { // returns all the underlying classifier values

        itkArrayType displayValues(m_numberOfDisplayValues);


            if (m_valid && m_transformUpToDate) {
                displayValues[kCorticalThresold] = m_CB;
                displayValues[kNotBoneThreshold] = m_ST;
                displayValues[kThreshold] = m_threshold;
                displayValues[kClassification] = m_classificationType;
                displayValues[kCorticalDensity] = m_CBLine.meanY;
                displayValues[kDenseDensity] = m_DBLine.meanY;
                displayValues[kTotalCorticalDensity] = m_totLine.meanY;
                displayValues[kNotBoneDensity] = m_STLine.meanY;
                displayValues[kTrabecularDensity] = m_TBLine.meanY;
                displayValues[kCorticalWidth] = (m_cbIndices[1]- m_cbIndices[0])*m_linearRegressionCalculator.GetProfileIncrement();
                displayValues[kEndocorticalWidth] = (m_ecIndices[1]- m_ecIndices[0])*m_linearRegressionCalculator.GetProfileIncrement();
                displayValues[kEndocorticalSlope] = m_ECLine.slope;
                displayValues[kPorosity] = m_porosity;
                displayValues[kProfileNumber] = m_linearRegressionCalculator.GetProfileNumber();
            } else {
                for (int i = 0; i < m_numberOfDisplayValues; i++) {
                    displayValues[i] = nan("1");
                }
            }
        return displayValues;
    }

    std::string ClassifierTransform::GetDisplayName(int index) const {

            if(index == kInvalid) {
                return "Blank";
            } else if(index == kCorticalThresold) {
                return "CB Threshold";
            } else if(index == kNotBoneThreshold) {
                return "NB Threshold";
            } else if(index == kThreshold) {
                return "Threshold Density";
            } else if(index == kClassification) {
                return "Classification";
            } else if(index == kCorticalDensity) {
                return "Cortical Bone Density";
            } else if(index == kDenseDensity) {
                return "Dense Bone Density";
            } else if(index == kTotalCorticalDensity) {
                return "Cortical + Endo Density";
            } else if(index == kNotBoneDensity) {
                return "Not Bone Density";
            } else if(index == kTrabecularDensity) {
                return "Trabecular Bone Density";
            } else if(index == kCorticalWidth) {
                return "Cortical Width";
            } else if(index == kEndocorticalWidth) {
                return "Endocrotical Region Width";
            } else if(index == kEndocorticalSlope) {
                return "Endocortical Region Slope";
            } else if(index == kPorosity) {
                return "Porosity";
            } else if(index == kProfileNumber) {
                return "Number of Sampled profiles";
            } else {
                std::cerr<<"Error: ClassifierTransform::GetDisplayName() invalid GetDisplayName index: "<<index<<endl;
                return "Invalid Selection";
            }
    }

    std::string ClassifierTransform::GetDisplayNameShort(int index) const {
            if(index == kInvalid) {
                return "Blank";
            } else if(index == kCorticalThresold) {
                return "thresh<sub>cb</sub>";
            } else if(index == kNotBoneThreshold) {
                return "thresh<sub>nb</sub>";
            } else if(index == kThreshold) {
                return "Y<sub>thresh</sub>";
            } else if(index == kClassification) {
                return "Classification";
            } else if(index == kCorticalDensity) {
                return "Y<sub>cb</sub>";
            } else if(index == kDenseDensity) {
                return "Y<sub>db</sub>";
            } else if(index == kTotalCorticalDensity) {
                return "Y<sub>cb+ec</sub>";
            } else if(index == kNotBoneDensity) {
                return "Y<sub>nb</sub>";
            } else if(index == kTrabecularDensity) {
                return "Y<sub>tb</sub>";
            } else if(index == kCorticalWidth) {
                return "T<sub>cb</sub>";
            } else if(index == kEndocorticalWidth) {
                return "T<sub>ec</sub>";
            } else if(index == kEndocorticalSlope) {
                return "m<sub>ec</sub>";
            } else if(index == kPorosity) {
                return "\u03D5";
            } else if(index == kProfileNumber) {
                return "Number of Sampled profiles";
            } else {
                std::cerr<<"Error: ClassifierTransform::GetDisplayName() invalid GetDisplayName index: "<<index<<endl;
                return "Invalid Selection";
            }
    }

    ClassifierTransform::itkArray2DType ClassifierTransform::GetDisplayRanges(bool includeImportedOptions) const {

        int numberOfDisplays = (includeImportedOptions) ? m_numberOfDisplayValues + 3 : m_numberOfDisplayValues;

        itkArray2DType displayRanges(numberOfDisplays, 2);

        displayRanges.put(kCorticalThresold,0,400.0);               displayRanges.put(kCorticalThresold,1,1600.0);
        displayRanges.put(kNotBoneThreshold,0,-200.0);              displayRanges.put(kNotBoneThreshold,1,200.0);
        displayRanges.put(kThreshold,0,200.0);                      displayRanges.put(kThreshold,1,600.0);
        displayRanges.put(kClassification,0,0.0);                   displayRanges.put(kClassification,1,4.0);
        displayRanges.put(kCorticalDensity, 0,400.0);               displayRanges.put(kCorticalDensity, 1,1600.0);
        displayRanges.put(kDenseDensity, 0,400.0);                  displayRanges.put(kDenseDensity, 1,1600.0);
        displayRanges.put(kTotalCorticalDensity, 0,200.0);          displayRanges.put(kTotalCorticalDensity, 1,1600.0);
        displayRanges.put(kNotBoneDensity,0,-200.0);                displayRanges.put(kNotBoneDensity,1,200.0);
        displayRanges.put(kTrabecularDensity,0,-200.0);             displayRanges.put(kTrabecularDensity,1,800.0);
        displayRanges.put(kCorticalWidth, 0,0.0);                   displayRanges.put(kCorticalWidth, 1,6.0);
        displayRanges.put(kEndocorticalWidth, 0,0.0);               displayRanges.put(kEndocorticalWidth, 1,6.0);
        displayRanges.put(kEndocorticalSlope, 0,-1500.0);           displayRanges.put(kEndocorticalSlope, 1,0.0);
        displayRanges.put(kPorosity, 0,0.0);                        displayRanges.put(kPorosity, 1,1.0);
        displayRanges.put(kProfileNumber,0,0.0);                    displayRanges.put(kProfileNumber,1,1000.0);

        return displayRanges;
    }

    ClassifierTransform::itkArrayType ClassifierTransform::GetClassifications() const {

        if(!ReadyToProcess() || !IsUpToDate()) {
            int numberOfSamples = m_linearRegressionCalculator.GetProfileSampleNumber();
            itkArrayType array = itkArrayType(numberOfSamples); array.Fill(nan("1"));
            return array;
        } else if(m_modeIndex==kMedian) {
            return m_filteredTissueValues;
        } else {
            return m_tissueValues;
        }
    }

    ClassifierTransform::itkArrayType ClassifierTransform::GetPercentages() const {
        if(!ReadyToProcess() || !IsUpToDate()) {
            int numberOfSamples = m_linearRegressionCalculator.GetProfileSampleNumber();
            itkArrayType array = itkArrayType(numberOfSamples); array.Fill(nan("1"));
            return array;
        } else {
            return m_percentArray; // % bone;
        }
    }

    // get results and display tables
    ClassifierTransform::TableType ClassifierTransform::GetProfileDisplayTable() const {
        TableType profileTable = TableType::New();
        int numberOfSamples = m_linearRegressionCalculator.GetProfileSampleNumber();

        profileTable->SetNumberOfRows(numberOfSamples);

        // add positions and mean values
        ArrayType positions = Utilities::convertItkToVtkArray(m_positions, "Positions");
        ArrayType values = Utilities::convertItkToVtkArray(m_meanProfiles, "Image Profile");
        profileTable->AddColumn(positions); profileTable->AddColumn(values); // the average


        int numberOfProfiles = m_linearRegressionCalculator.GetProfileNumber();
        for(int j=0; j<numberOfProfiles; j++) {
            ArrayType array = ArrayType::New(); array->SetNumberOfValues(numberOfSamples);
            std::stringstream converter; converter << j;
            std::string name = "Profile Values " + converter.str();
            array->SetName(name.c_str());

            for(int i=0; i<numberOfSamples; i++) {
                array->SetValue(i, m_multipleProfiles[i][j]);
            }
            profileTable->AddColumn(array);
        }

        return profileTable;
    }

    ClassifierTransform::TableType ClassifierTransform::GetModelDisplayTable() const {

        int numberOfRows = 2;
        ScalarType increment = m_linearRegressionCalculator.GetProfileIncrement();
        ScalarType profileStartLocation = m_linearRegressionCalculator.GetProfileStartPosition();
        int startIndex = m_linearRegressionCalculator.GetProfileStartIndex();

        TableType profileTable = TableType::New();
        profileTable->SetNumberOfRows(2);

        if(!m_transformUpToDate) {
            return profileTable;
        }

        // st segment
        if(m_STLine.n != -1) {

            ArrayType positions = ArrayType::New(), values = ArrayType::New();
            positions->SetName("ST Positions"); positions->SetNumberOfValues(numberOfRows);
            values->SetName("ST Line"); values->SetNumberOfValues(numberOfRows);

            ScalarType startLocation = profileStartLocation + increment * m_stIndices[0], endLocation = profileStartLocation + increment * m_stIndices[1];
            positions->SetValue(0, startLocation); values->SetValue(0, m_linearRegressionCalculator.CalculateLinearValue(m_STLine, startLocation));
            positions->SetValue(1, endLocation); values->SetValue(1, m_linearRegressionCalculator.CalculateLinearValue(m_STLine, endLocation));
            profileTable->AddColumn(positions); profileTable->AddColumn(values);
        }

        // cb segment
        if(m_CBLine.n != -1) {

            ArrayType positions = ArrayType::New(), values = ArrayType::New();
            positions->SetName("CB Positions"); positions->SetNumberOfValues(numberOfRows);
            values->SetName("CB Line"); values->SetNumberOfValues(numberOfRows);

            ScalarType startLocation = profileStartLocation + increment * m_cbIndices[0], endLocation = profileStartLocation + increment * m_cbIndices[1];
            positions->SetValue(0, startLocation); values->SetValue(0, m_linearRegressionCalculator.CalculateLinearValue(m_CBLine, startLocation));
            positions->SetValue(1, endLocation); values->SetValue(1, m_linearRegressionCalculator.CalculateLinearValue(m_CBLine, endLocation));
            profileTable->AddColumn(positions); profileTable->AddColumn(values);
        }

        // dense bone (CB minus the endocortical region)
        if(m_DBLine.n != -1) {

            ArrayType positions = ArrayType::New(), values = ArrayType::New();
            positions->SetName("DB Positions"); positions->SetNumberOfValues(numberOfRows);
            values->SetName("DB Line"); values->SetNumberOfValues(numberOfRows);

            ScalarType startLocation = profileStartLocation + increment * m_dbIndices[0], endLocation = profileStartLocation + increment * m_dbIndices[1];
            positions->SetValue(0, startLocation); values->SetValue(0, m_linearRegressionCalculator.CalculateLinearValue(m_DBLine, startLocation));
            positions->SetValue(1, endLocation); values->SetValue(1, m_linearRegressionCalculator.CalculateLinearValue(m_DBLine, endLocation));
            profileTable->AddColumn(positions); profileTable->AddColumn(values);
        }

        // tb segment
        if(m_TBLine.n != -1) {

            ArrayType positions = ArrayType::New(), values = ArrayType::New();
            positions->SetName("TB Positions"); positions->SetNumberOfValues(numberOfRows);
            values->SetName("TB Line"); values->SetNumberOfValues(numberOfRows);

            ScalarType startLocation = profileStartLocation + increment * m_tbIndices[0], endLocation = profileStartLocation + increment * m_tbIndices[1];
            positions->SetValue(0, startLocation); values->SetValue(0, m_linearRegressionCalculator.CalculateLinearValue(m_TBLine, startLocation));
            positions->SetValue(1, endLocation); values->SetValue(1, m_linearRegressionCalculator.CalculateLinearValue(m_TBLine, endLocation));
            profileTable->AddColumn(positions); profileTable->AddColumn(values);
        }

        // ec region segment
        if(m_ECLine.n != -1) {

            ArrayType positions = ArrayType::New(), values = ArrayType::New();
            positions->SetName("EC Positions"); positions->SetNumberOfValues(numberOfRows);
            values->SetName("EC Line"); values->SetNumberOfValues(numberOfRows);

            ScalarType startLocation = profileStartLocation + increment * m_ecIndices[0], endLocation = profileStartLocation + increment * m_ecIndices[1];
            positions->SetValue(0, startLocation); values->SetValue(0, m_linearRegressionCalculator.CalculateLinearValue(m_ECLine, startLocation));
            positions->SetValue(1, endLocation); values->SetValue(1, m_linearRegressionCalculator.CalculateLinearValue(m_ECLine, endLocation));
            profileTable->AddColumn(positions); profileTable->AddColumn(values);
        }

        // total bone (CB minus the endocortical region)
        if(m_totLine.n != -1) {

            ArrayType positions = ArrayType::New(), values = ArrayType::New();
            positions->SetName("Tot Positions"); positions->SetNumberOfValues(numberOfRows);
            values->SetName("Tot Line"); values->SetNumberOfValues(numberOfRows);

            ScalarType startLocation = profileStartLocation + increment * m_totIndices[0], endLocation = profileStartLocation + increment * m_totIndices[1];
            positions->SetValue(0, startLocation); values->SetValue(0, m_linearRegressionCalculator.CalculateLinearValue(m_totLine, startLocation));
            positions->SetValue(1, endLocation); values->SetValue(1, m_linearRegressionCalculator.CalculateLinearValue(m_totLine, endLocation));
            profileTable->AddColumn(positions); profileTable->AddColumn(values);
        }

        return profileTable;
    }

    ClassifierTransform::TableType ClassifierTransform::GetModelDisplayTable(ArrayType parametersIn, double offset) const {

        int numberOfRows = 2;

        TableType profileTable = TableType::New();
        profileTable->SetNumberOfRows(numberOfRows);

        ArrayType positions, values; ScalarType x1, x2, m, b;

        // st line
        x1 = parametersIn->GetValue(kSTStart); x2 = parametersIn->GetValue(kSTEnd);
        m = parametersIn->GetValue(kSTSlope); b = parametersIn->GetValue(kSTIntercept);

        positions = ArrayType::New(); positions->SetName("ST Positions");
        positions->SetNumberOfValues(numberOfRows);
        positions->SetValue(0, x1-offset); positions->SetValue(1, x2-offset);
        profileTable->AddColumn(positions);

        values = ArrayType::New(); values->SetName("ST Values");
        values->SetNumberOfValues(numberOfRows);
        values->SetValue(0, x1*m+b); values->SetValue(1, x2*m+b);
        profileTable->AddColumn(values);

        // cb line
        x1 = parametersIn->GetValue(kDBStart); x2 = parametersIn->GetValue(kDBEnd);
        m = parametersIn->GetValue(kDBSlope); b = parametersIn->GetValue(kDBIntercept);

        positions = ArrayType::New(); positions->SetName("DB Positions");
        positions->SetNumberOfValues(numberOfRows);
        positions->SetValue(0, x1-offset); positions->SetValue(1, x2-offset);
        profileTable->AddColumn(positions);

        values = ArrayType::New(); values->SetName("DB Values");
        values->SetNumberOfValues(numberOfRows);
        values->SetValue(0, x1*m+b); values->SetValue(1, x2*m+b);
        profileTable->AddColumn(values);

        // cb line
        x1 = parametersIn->GetValue(kCBStart); x2 = parametersIn->GetValue(kCBEnd);
        m = parametersIn->GetValue(kCBSlope); b = parametersIn->GetValue(kCBIntercept);

        positions = ArrayType::New(); positions->SetName("CB Positions");
        positions->SetNumberOfValues(numberOfRows);
        positions->SetValue(0, x1-offset); positions->SetValue(1, x2-offset);
        profileTable->AddColumn(positions);

        values = ArrayType::New(); values->SetName("CB Values");
        values->SetNumberOfValues(numberOfRows);
        values->SetValue(0, x1*m+b); values->SetValue(1, x2*m+b);
        profileTable->AddColumn(values);

        // tb line
        x1 = parametersIn->GetValue(kTBStart); x2 = parametersIn->GetValue(kTBEnd);
        m = parametersIn->GetValue(kTBSlope); b = parametersIn->GetValue(kTBIntercept);

        positions = ArrayType::New(); positions->SetName("TB Positions");
        positions->SetNumberOfValues(numberOfRows);
        positions->SetValue(0, x1-offset); positions->SetValue(1, x2-offset);
        profileTable->AddColumn(positions);

        values = ArrayType::New(); values->SetName("TB Values");
        values->SetNumberOfValues(numberOfRows);
        values->SetValue(0, x1*m+b); values->SetValue(1, x2*m+b);
        profileTable->AddColumn(values);

        // ec line
        x1 = parametersIn->GetValue(kECStart); x2 = parametersIn->GetValue(kECEnd);
        m = parametersIn->GetValue(kECSlope); b = parametersIn->GetValue(kECIntercept);

        positions = ArrayType::New(); positions->SetName("MC Positions");
        positions->SetNumberOfValues(numberOfRows);
        positions->SetValue(0, x1-offset); positions->SetValue(1, x2-offset);
        profileTable->AddColumn(positions);

        values = ArrayType::New(); values->SetName("MC Values");
        values->SetNumberOfValues(numberOfRows);
        values->SetValue(0, x1*m+b); values->SetValue(1, x2*m+b);
        profileTable->AddColumn(values);

        return profileTable;
    }

    ClassifierTransform::TableType ClassifierTransform::GetClassifierTable() const {

        int numberOfSamples = m_linearRegressionCalculator.GetProfileSampleNumber();

        vtkSmartPointer<vtkDoubleArray> positions = vtkSmartPointer<vtkDoubleArray>::New();
        positions->SetName("position"); positions->SetNumberOfValues(numberOfSamples*2);
        vtkSmartPointer<vtkDoubleArray> tissueClassifications = vtkSmartPointer<vtkDoubleArray>::New();
        tissueClassifications->SetName("Tissue Classifications"); tissueClassifications->SetNumberOfValues(numberOfSamples*2);
        vtkSmartPointer<vtkDoubleArray> percentages = vtkSmartPointer<vtkDoubleArray>::New();
        percentages->SetName("Bone Percentages"); percentages->SetNumberOfValues(numberOfSamples*2);


        if(m_transformUpToDate) {
            // set positions
            ScalarType increment = m_linearRegressionCalculator.GetProfileIncrement();
            ScalarType  startPosition = m_linearRegressionCalculator.GetProfileStartPosition();
            positions->SetValue(0, startPosition);
            for(int i=1; i<numberOfSamples*2-1; i++) {
                ScalarType index = ((int)((i-1)/2))+0.5;
                positions->SetValue(i, startPosition+increment*index);
            }
            positions->SetValue(numberOfSamples*2-1, m_positions[numberOfSamples-1]);
            // set classifications
            for(int i=0; i<numberOfSamples*2; i++) {
                tissueClassifications->SetValue(i, m_tissueValues[i/2]);
                percentages->SetValue(i, m_percentArray[i/2]*10); // to give a range 0 to 1000 todo - remove or standardise
            }

        } else {
            for(int i=0; i<numberOfSamples*2-2; i++) {
                positions->SetValue(i, nan("1"));
                tissueClassifications->SetValue(i, nan("1"));
                percentages->SetValue(i, nan("1"));
            }
        }

        TableType classificationTable = TableType::New();
        classificationTable->SetNumberOfRows(numberOfSamples*2-2);
        classificationTable->AddColumn(positions);
        classificationTable->AddColumn(tissueClassifications);
        classificationTable->AddColumn(percentages);

        if(m_modeIndex==kMedian && m_transformUpToDate) {
            vtkSmartPointer<vtkDoubleArray> filteredTissueClassifications = vtkSmartPointer<vtkDoubleArray>::New();
            filteredTissueClassifications->SetName("Filtered Tissue Classifications"); filteredTissueClassifications->SetNumberOfValues(numberOfSamples*2);
            for(int i=0; i<numberOfSamples*2; i++) {
                filteredTissueClassifications->SetValue(i, m_filteredTissueValues[i/2]);
            }
            classificationTable->AddColumn(filteredTissueClassifications);
        }

        return classificationTable;

    }


    //----------- private ---------------//
    // class based classifier

    void ClassifierTransform::CalculateTissueRegressionLines() { // assumes multiple profiles

        m_STLine = m_linearRegressionCalculator.CalculateLinearRegressionMultiple(m_positions, m_multipleProfiles, (int)ceil(m_stIndices[0]), (int)floor(m_stIndices[1]));
        m_CBLine = m_linearRegressionCalculator.CalculateLinearRegressionMultiple(m_positions, m_multipleProfiles, (int)ceil(m_cbIndices[0]), (int)floor(m_cbIndices[1]));
        m_TBLine = m_linearRegressionCalculator.CalculateLinearRegressionMultiple(m_positions, m_multipleProfiles, (int)ceil(m_tbIndices[0]), (int)floor(m_tbIndices[1]));
        m_ECLine = m_linearRegressionCalculator.CalculateLinearRegressionMultiple(m_positions, m_multipleProfiles, (int)ceil(m_ecIndices[0]), (int)floor(m_ecIndices[1]));
        m_DBLine = m_linearRegressionCalculator.CalculateLinearRegressionMultiple(m_positions, m_multipleProfiles, (int)ceil(m_dbIndices[0]), (int)floor(m_dbIndices[1]));
        m_totLine = m_linearRegressionCalculator.CalculateLinearRegressionMultiple(m_positions, m_multipleProfiles, (int)ceil(m_totIndices[0]), (int)floor(m_totIndices[1]));
    }

    // percentage based classifier
    bool ClassifierTransform::PercentageBasedClassifing() {
        CalculatePercentages();
        bool status = CalculateTissueSectionsStrict();
        //bool status = CalculateTissueSectionsDirect();
        CalculateTissueRegressionLines();
        return status;
    }

    bool ClassifierTransform::CalculatePercentages() { // telly the percentage of each profile exceeding the threshold at each profile location

        // get general profile details
        int numberOfSamples = m_linearRegressionCalculator.GetProfileSampleNumber();
        int startIndex = m_linearRegressionCalculator.GetProfileStartIndex();
        int endIndex = m_linearRegressionCalculator.GetProfileEndIndex();
        int numberOfProfiles = m_linearRegressionCalculator.GetProfileNumber();

        // check if worth proceeding
        if(numberOfProfiles <=1 || startIndex==kInvalid || endIndex!=numberOfSamples-1 || !m_thresholdsSet) {
            cerr<<"Warning: Terminating ClassifierTransform::ClassifyGlobalSingleThreshold() as multiple profiles required and valid start and end indices are required"<<endl;
            return false; // only run if multiple profiles have been sampled.
        }


        // cycle through the all the positions and categorise by number of profiles greater than the threshold
        for(int i=0; i<numberOfSamples; i++) {
            if(i>=startIndex && i<=endIndex) {
                int numberAboveThreshold = 0;
                for(int j=0; j<numberOfProfiles; j++) {
                    if(m_multipleProfiles[i][j] >= m_threshold) {
                        numberAboveThreshold++;
                    }
                }
                m_percentArray[i]=100.0*((ScalarType)numberAboveThreshold)/((ScalarType)numberOfProfiles); // % bone
            }
        }
        return true;
    }

    bool ClassifierTransform::CalculateTissueSectionsDirect() {
        int edgeIndex = m_linearRegressionCalculator.GetProfileEdgeIndex();
        int numberOfSamples = m_linearRegressionCalculator.GetProfileSampleNumber();
        int numberOfProfiles = m_linearRegressionCalculator.GetProfileNumber();
        double profileIncrement = m_linearRegressionCalculator.GetProfileIncrement();
        int startIndex, endIndex;

        // i. initial in bounds check - the start needs to be before -3mm
        if(profileIncrement*m_linearRegressionCalculator.GetProfileStartIndex() + m_linearRegressionCalculator.GetProfileStartPosition()>-3) {
            ClearResults();
            return false;
        }

        // 1. find the max 'percentage' in-between -3 to 3mm
        double maxPercentage = 0; int maxPercentageIndex = 0;
        startIndex = std::max((int)(edgeIndex - m_profileRange / profileIncrement), 0);
        endIndex = std::min((int)(edgeIndex + m_profileRange / profileIncrement), numberOfSamples);
        for(int i=startIndex; i<=endIndex; i++) {
            if(m_percentArray[i]>maxPercentage) { // record first instance of max
                maxPercentageIndex =i; maxPercentage=m_percentArray[i];
            }
        }

        if(maxPercentage < 25.0) { // todo - consider weather or not to try look further along the profile for the maximum
            ClearResults();
            return false;
        }

        // 2. initial estimate of cortical region
        //    a. xP - 1st % reaching 1/2 of the max %
        //    b. xE - last % below 1/2 of max % after max % has occured
        startIndex = m_linearRegressionCalculator.GetProfileStartIndex();
        endIndex = m_linearRegressionCalculator.GetProfileEndIndex();
        bool prevGreaterThanThreshold = false;
        ScalarType xPIndex=nan("1"), xEIndex=nan("1");
        double thresholdPercentage = maxPercentage*0.5;
        for(int i=startIndex; i<=endIndex; i++) {

            if(m_percentArray[i]<=thresholdPercentage && prevGreaterThanThreshold && i>=maxPercentageIndex) {
                xEIndex = (i*(m_percentArray[i-1]-thresholdPercentage)+(i-1)*(thresholdPercentage-m_percentArray[i]))/(m_percentArray[i-1]-m_percentArray[i]);
                //xEIndex=i-1;
                break; // must estimate cb end as used as such // todo linearly interploate
            }

            if(m_percentArray[i]>=thresholdPercentage && i>=1) { // check to avoid indexing out of bounds
                prevGreaterThanThreshold = true;
                if(isnan(xPIndex)) {
                    //xPIndex = i; // todo linearly interpolate
                    xPIndex = (i*(thresholdPercentage-m_percentArray[i-1])+(i-1)*(m_percentArray[i]-thresholdPercentage))/(m_percentArray[i]-m_percentArray[i-1]);
                }
            } else {
                prevGreaterThanThreshold = false;
            }
        }
        if(startIndex == -1 || isnan(xPIndex) || xPIndex<=0) {
            ClearResults();
            return false;
        } else {
            m_stIndices[0] = startIndex; m_stIndices[1] = m_cbIndices[0] = m_dbIndices[0] = m_totIndices[0] = xPIndex;
        }

        if(isnan(xEIndex)) { // if never drops back down to threshold select min through 10
            double minPercent = maxPercentage; // todo - decide weather or not to presue the last ditch attempt to get a reasonable measurement
            for(int i=maxPercentageIndex; i<=numberOfSamples*0.8; i++) {
                if(m_percentArray[i]<minPercent) { // record first instance of max
                    minPercent=m_percentArray[i]; xEIndex=i;
                }
            }
            if(isnan(xEIndex)) {
                ClearResults();
                return false;
            }
        }

        // 3. final estimate of cortical region
        //    a. xP - as above
        //    b. xEcb - 1st percentage below the max (as measured from the centre of the cortical region
        ScalarType xEcbIndex = nan("1");
        startIndex = (int) (0.5*(xPIndex+xEIndex));
        startIndex = (startIndex > maxPercentageIndex) ? startIndex : maxPercentageIndex; // ensure start after max
        thresholdPercentage=maxPercentage-100.0/numberOfProfiles;
        if(m_percentArray[startIndex]>thresholdPercentage) {
            endIndex=(int)floor(xEIndex);
            for(int i=startIndex; i<=endIndex; i++) {
                if(m_percentArray[i]<=thresholdPercentage) { // todo - consider making it look for a minimum
                    //xEcbIndex=i; break; // todo linearly interpolate
                    xEcbIndex = (i*(m_percentArray[i-1]-thresholdPercentage)+(i-1)*(thresholdPercentage-m_percentArray[i]))/(m_percentArray[i-1]-m_percentArray[i]);
                    break;
                }
            }
        } else {
            endIndex=(int)ceil(xPIndex);
            for(int i=startIndex; i>=endIndex; i--) {
                if(m_percentArray[i]>thresholdPercentage) {
                    //xEcbIndex=i+1; break; // todo linearly interpolate
                    xEcbIndex = ((i+1)*(m_percentArray[i]-thresholdPercentage)+i*(thresholdPercentage-m_percentArray[i+1]))/(m_percentArray[i]-m_percentArray[i+1]);
                    break;
                }
            }
        }
        if(isnan(xEcbIndex)) {
            ClearResults();
            return false;
        } else {
            m_dbIndices[1] = m_ecIndices[0] = xEcbIndex;
        }

        // 4.1 estimate of the end of the endocortical region
        //    a. xEcb - as above
        //    b. xEtb - 1st minimum after the percentages fall below the estimated mean tb percentage
        //    c. xE - refine to the 1/2 percent using the meanTBPercentage
        endIndex = m_linearRegressionCalculator.GetProfileEndIndex(); startIndex = ceil(xEcbIndex);
        //endIndex = (endIndex < xEIndex + (int)(10.0/profileIncrement)) ? endIndex : xEIndex + (int)(10.0/profileIncrement); // todo consider limiting amount of tb to use for estimating the mean
        double meanTBPercentage = m_linearRegressionCalculator.CalculateArrayMean(m_percentArray, ceil(xEIndex), endIndex);
        thresholdPercentage=0.5*(maxPercentage+meanTBPercentage);
        bool prevLessThanThreshold = false; prevGreaterThanThreshold = true;
        ScalarType xEtbIndex = nan("1");
        for(int i=startIndex; i<=endIndex; i++) {

            if(prevLessThanThreshold && (m_percentArray[i-1]<m_percentArray[i] || m_percentArray[i]==0) && i-2>=0) {
                //xEtbIndex = i-1; break; // todo linearly interpolate
                xEtbIndex=((i-1.5)*(m_percentArray[i-2]-m_percentArray[i-1])+(i-0.5)*(m_percentArray[i]-m_percentArray[i-1]))/(m_percentArray[i-2]+m_percentArray[i]-2.0*m_percentArray[i-1]);
                break;
            }

            if(m_percentArray[i]<=thresholdPercentage && prevGreaterThanThreshold && i>=maxPercentageIndex) { // i.e. 1/2 between the cb max and the average tb percentage
                //xEIndex=i-1; prevGreaterThanThreshold=false; // must estimate cb end as used as such todo - decide which to use // todo linearly interpolate

                xEIndex=(i*(m_percentArray[i-1]-thresholdPercentage)+(i-1)*(thresholdPercentage-m_percentArray[i]))/(m_percentArray[i-1]-m_percentArray[i]);
                prevGreaterThanThreshold=false;
            }

            prevLessThanThreshold= m_percentArray[i]<=meanTBPercentage;
        }
        if(isnan(xEtbIndex)||isnan(xEIndex)) {
            ClearResults();
            return false;
        } else {
            m_ecIndices[1] = m_totIndices[1] = m_tbIndices[0] = xEtbIndex;
            m_tbIndices[1] = endIndex; m_cbIndices[1] = xEIndex;
        }

        // 5. calculate the porosity from %'s
        startIndex = (int)ceil(m_totIndices[0]); endIndex = (int)floor(m_totIndices[1]);
        ScalarType percentageSum=0;
        for(int i=startIndex; i<=endIndex; i++) {
            percentageSum+=m_percentArray[i]; // in percent
        }
        m_porosity = (percentageSum/(endIndex+1-startIndex))/maxPercentage;

        // 6. calculate the linear regression lines of the the resulting sections
        //    a. Not Bone - xP-6 to xP
        //    b. Cortical Bone - xP to xE
        //    c. Trabecular Bone - xE to xP+12
        // preformed one level up

        //cout<<" % threshold="<<thresholdPercentage<<", tb threshold="<<meanTBPercentage<<", xP="<<xPIndex*profileIncrement-6<<", xE="<<xEIndex*profileIncrement-6<<", xEcb="<<xEcbIndex*profileIncrement-6<<", xEtb="<<xEtbIndex*profileIncrement-6<<endl;

        return !(isnan(xEtbIndex) || endIndex == -1);
    }

    bool ClassifierTransform::CalculateTissueSectionsStrict() {
        int edgeIndex = m_linearRegressionCalculator.GetProfileEdgeIndex();
        int numberOfSamples = m_linearRegressionCalculator.GetProfileSampleNumber();
        int numberOfProfiles = m_linearRegressionCalculator.GetProfileNumber();
        double profileIncrement = m_linearRegressionCalculator.GetProfileIncrement();
        int startIndex, endIndex;

        // i. initial in bounds check - the start needs to be before -3mm
        if(profileIncrement*m_linearRegressionCalculator.GetProfileStartIndex() + m_linearRegressionCalculator.GetProfileStartPosition()>-3) {
            ClearResults();
            return false;
        }

        // 1. find the max 'percentage' in-between -3 to 3mm
        double maxPercentage = 0; int maxPercentageIndex = 0;
        startIndex = std::max((int)(edgeIndex - m_profileRange / profileIncrement), 0);
        endIndex = std::min((int)(edgeIndex + m_profileRange / profileIncrement), numberOfSamples);
        for(int i=startIndex; i<=endIndex; i++) {
            if(m_percentArray[i]>maxPercentage) { // record first instance of max
                maxPercentageIndex =i; maxPercentage=m_percentArray[i];
            }
        }

        if(maxPercentage < 25.0) { // todo - consider weather or not to try look further along the profile for the maximum
            ClearResults();
            return false;
        }

        // 2. initial estimate of cortical region
        //    a. xP - 1st % reaching 1/2 of the max %
        //    b. xE - last % below 1/2 of max % after max % has occured
        startIndex = m_linearRegressionCalculator.GetProfileStartIndex();
        endIndex = m_linearRegressionCalculator.GetProfileEndIndex();
        bool prevGreaterThanThreshold = false;
        ScalarType xPIndex=nan("1"), xEIndex=nan("1");
        double thresholdPercentage = maxPercentage*0.5;
        for(int i=startIndex; i<=endIndex; i++) {

            if(m_percentArray[i]<=thresholdPercentage && prevGreaterThanThreshold && i>=maxPercentageIndex) {
                xEIndex = (i*(m_percentArray[i-1]-thresholdPercentage)+(i-1)*(thresholdPercentage-m_percentArray[i]))/(m_percentArray[i-1]-m_percentArray[i]);
                break; // must estimate cb end as used as such
            }

            if(m_percentArray[i]>=thresholdPercentage && i>=1) { // check to avoid indexing out of bounds
                prevGreaterThanThreshold = true;
                if(isnan(xPIndex)) {
                    xPIndex = (i*(thresholdPercentage-m_percentArray[i-1])+(i-1)*(m_percentArray[i]-thresholdPercentage))/(m_percentArray[i]-m_percentArray[i-1]);
                }
            } else {
                prevGreaterThanThreshold = false;
            }
        }
        if(startIndex == -1 || isnan(xPIndex) || xPIndex<=0) {
            ClearResults();
            return false;
        } else {
            m_stIndices[0] = startIndex; m_stIndices[1] = m_cbIndices[0] = m_dbIndices[0] = m_totIndices[0] = xPIndex;
        }

        if(isnan(xEIndex)) { // if never drops back down to threshold select min through 10
            double minPercent = maxPercentage; // todo - decide weather or not to presue the last ditch attempt to get a reasonable measurement
            for(int i=maxPercentageIndex; i<=numberOfSamples*0.8; i++) {
                if(m_percentArray[i]<minPercent) { // record minimum values // todo consider linear interpolation
                    minPercent=m_percentArray[i]; xEIndex=i;
                }
            }
            if(isnan(xEIndex)) {
                ClearResults();
                return false;
            }
        }

        // 3. final estimate of cortical region
        //    a. xP - as above
        //    b. xEcb - 1st percentage below the max (as measured from the centre of the cortical region
        ScalarType xEcbIndex = nan("1");
        startIndex = maxPercentageIndex;
        thresholdPercentage=maxPercentage-100.0/numberOfProfiles;
        endIndex=(int)floor(xEIndex);
        for(int i=startIndex; i<=endIndex; i++) {
            if(m_percentArray[i]<=thresholdPercentage) { // todo - consider making it look for a minimum
                xEcbIndex = (i*(m_percentArray[i-1]-thresholdPercentage)+(i-1)*(thresholdPercentage-m_percentArray[i]))/(m_percentArray[i-1]-m_percentArray[i]);
                break;
            }
        }
        if(isnan(xEcbIndex)) {
            ClearResults();
            return false;
        } else {
            m_dbIndices[1] = m_ecIndices[0] = xEcbIndex;
        }

        // 4.1 estimate of the end of the endocortical region
        //    a. xEcb - as above
        //    b. xEtb - 1st minimum after the percentages fall below the estimated mean tb percentage
        //    c. xE - refine to the 1/2 percent using the meanTBPercentage
        endIndex = m_linearRegressionCalculator.GetProfileEndIndex(); startIndex = (int)ceil(xEcbIndex);
        //endIndex = (endIndex < xEIndex + (int)(10.0/profileIncrement)) ? endIndex : xEIndex + (int)(10.0/profileIncrement); // todo consider limiting amount of tb to use for estimating the mean
        double meanTBPercentage = m_linearRegressionCalculator.CalculateArrayMean(m_percentArray, (int)ceil(xEIndex), endIndex);
        thresholdPercentage=0.5*(maxPercentage+meanTBPercentage);
        bool prevLessThanThreshold = false; prevGreaterThanThreshold = true;
        ScalarType xEtbIndex = nan("1");
        for(int i=startIndex; i<=endIndex; i++) {

            if(prevLessThanThreshold && (m_percentArray[i-1]<m_percentArray[i] || m_percentArray[i-1]==0) && i-2>=0) {
                xEtbIndex=((i-1.5)*(m_percentArray[i-2]-m_percentArray[i-1])+(i-0.5)*(m_percentArray[i]-m_percentArray[i-1]))/(m_percentArray[i-2]+m_percentArray[i]-2.0*m_percentArray[i-1]);
                break;
            }

            if(m_percentArray[i]<=thresholdPercentage && prevGreaterThanThreshold && i>=maxPercentageIndex) { // i.e. 1/2 between the cb max and the average tb percentage

                xEIndex=(i*(m_percentArray[i-1]-thresholdPercentage)+(i-1)*(thresholdPercentage-m_percentArray[i]))/(m_percentArray[i-1]-m_percentArray[i]);
                prevGreaterThanThreshold=false;
            }

            prevLessThanThreshold = (m_percentArray[i]<=meanTBPercentage);
        }
        if(isnan(xEtbIndex)||isnan(xEIndex)) {
            ClearResults();
            return false;
        } else {
            m_ecIndices[1] = m_totIndices[1] = m_tbIndices[0] = xEtbIndex;
            m_tbIndices[1] = endIndex; m_cbIndices[1] = xEIndex;
        }

        // 5. calculate the porosity from %'s
        startIndex = (int)ceil(m_totIndices[0]); endIndex = (int)floor(m_totIndices[1]);
        ScalarType percentageSum=0;
        for(int i=startIndex; i<=endIndex; i++) {
            percentageSum+=m_percentArray[i]; // in percent
        }
        m_porosity = (percentageSum/(endIndex+1-startIndex))/maxPercentage;


        // 6. calculate the linear regression lines of the the resulting sections
        //    a. Not Bone - xP-6 to xP
        //    b. Cortical Bone - xP to xE
        //    c. Trabecular Bone - xE to xP+12
        // preformed one level up

        //cout<<" % threshold="<<thresholdPercentage<<", tb threshold="<<meanTBPercentage<<", xP="<<xPIndex*profileIncrement-6<<", xE="<<xEIndex*profileIncrement-6<<", xEcb="<<xEcbIndex*profileIncrement-6<<", xEtb="<<xEtbIndex*profileIncrement-6<<endl;

        return !(isnan(xEtbIndex) || endIndex == -1);
    }

    void ClassifierTransform::MedianFilter(itkIntArrayType classArray, itkArrayType valueArray, int filterSize) {

        int filterHalfWidth = floor(filterSize/2.0);
        int startIndex = m_linearRegressionCalculator.GetProfileStartIndex();
        int endIndex = m_linearRegressionCalculator.GetProfileEndIndex();

        m_filteredTissueValues = itkArrayType(valueArray); m_filteredTissueTypes = itkIntArrayType(classArray); // todo names - "Filtered Densities" and "Filtered Tissue Types"

        if(startIndex==kInvalid || endIndex==kInvalid) {
            cerr<<"Warning in ClassifierTransform::MedianFilter() invalid start/end indices"<<endl;
            return;
        }

        int filterCount[3] = {0,0,0};
        ScalarType filterValues[3] = {m_ST, m_threshold, m_CB};

        for(int i=startIndex - filterHalfWidth; i<=endIndex; i++) {

            int addIndex = i+filterHalfWidth, removeIndex = i - (filterHalfWidth + 1); // add at fistance filterHW, remove at distance filter HW + 1
            // add value at the extreme end of filter (if it's within the range of indices with valid values)
            if(addIndex>= startIndex && addIndex<=endIndex) {
                int tissueType = classArray[addIndex];
                filterCount[tissueType]++;
            }
            // remove values now outside the filter
            if(removeIndex >= startIndex && removeIndex <= endIndex) {
                int tissueType = classArray[removeIndex];
                filterCount[tissueType]--;
            }
            //calculate and set the filter median
            if(i>=startIndex && i<=endIndex) {
                int medianType = round((kST*filterCount[kST] + kCB*filterCount[kCB] + kTB*filterCount[kTB])/(filterCount[kST] + filterCount[kTB] + filterCount[kCB]));
                m_filteredTissueValues[i] = filterValues[medianType];
                m_filteredTissueTypes[i] = medianType;
                //cerr<<"filtering cnt[0]="<<filterCount[0]<<", cnt[1]="<<filterCount[1]<<", cnt[2]="<<filterCount[2]<<"median type="<<medianType<<endl;
            }
        }

        return;
    }

    void ClassifierTransform::CalculateValidity() {

        m_valid = !isnan(m_ST) && !isnan(m_CB) && !isnan(m_threshold) && (m_ST < m_CB);
        return;
    }

    ClassifierTransform::ScalarType ClassifierTransform::CalculateAbsMeanError() const {

        if(!m_valid) {
            return nan("1");
        } else {
            return nan("1"); // AbsMeanError not calculated
        }

    }

    void ClassifierTransform::ClearResults() {


        m_classificationType = -1;

        m_stIndices[0] = m_stIndices[1] = m_cbIndices[0] = m_cbIndices[1] = m_dbIndices[0] = m_dbIndices[1] = nan("1");
        m_tbIndices[0] = m_tbIndices[1] = m_ecIndices[0] = m_ecIndices[1] = m_totIndices[0] = m_totIndices[1] = nan("1");

        m_porosity = nan("1");

        m_STLine = m_DBLine = m_totLine = m_CBLine = m_TBLine = m_ECLine = RegressionLineType();

        int numberOfSamples = m_linearRegressionCalculator.GetProfileSampleNumber();
        if(numberOfSamples == kInvalid) {
            return;
        }

        // todo names - m_tissueValues="Tissue Densities", m_tissueTypes="Tissue Classes"
        m_tissueValues = itkArrayType(numberOfSamples); m_tissueValues.Fill(nan("1"));
        m_tissueTypes = itkIntArrayType(numberOfSamples); m_tissueTypes.Fill(kInvalid);
        m_percentArray = itkArrayType(numberOfSamples); m_percentArray.Fill(nan("1"));

        m_filteredTissueValues = itkArrayType(numberOfSamples); m_filteredTissueValues.Fill(nan("1"));
        m_filteredTissueTypes = itkIntArrayType(numberOfSamples); m_filteredTissueTypes.Fill(kInvalid);

    }
}
