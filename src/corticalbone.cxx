#include <iostream>
#include <math.h> 
#include <istream>
#include <sstream>

#include "corticalbone.h"

//#include <itkArray.h>

// for defining the mesh type
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVRMLImporter.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkOBJReader.h>
#include <vtkOBJWriter.h>
#include <vtkDataSet.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkImageMapper3D.h>
#include <vtkPolyDataNormals.h>

// manually read in sw vrml file
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>

// for defining the DICOM type
#include "itkImageSeriesReader.h"
#include "itkCastImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include <itkImageToVTKImageFilter.h>

#include "vtkImageActor.h"
#include "vtkMapper.h"
#include "vtkPolyDataMapper.h"
#include "utilities.h"

// define interpolator type
#include <itkLinearInterpolateImageFunction.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <vtkPointData.h>
#include <vtkTable.h>

#include <itkMinimumMaximumImageCalculator.h>
#include <vtkLookupTable.h>

// define display types
#include <vtkSphereSource.h>
#include <vtkLineSource.h>

#include <vtkKdTreePointLocator.h>
#include <vtkMath.h>

#include <vtkColorTransferFunction.h>
#include <itkPowellOptimizer.h>

#include <vtkDelimitedTextWriter.h>
#include <vtkDelimitedTextReader.h>
#include <itkTIFFImageIO.h>
#include <itkCastImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <vtkPLYReader.h>
#include <vtkDICOMImageReader.h>
#include <itkChangeInformationImageFilter.h>


//--------------- Cortical Bone Methods --------------------//

//----------- Public Methods --------------//

namespace itk {

    //---------------- Public ---------------------//
    CorticalBone::CorticalBone() {


        verbose = false;
        // re-create all other objects each time a new project is loaded
        initialiseSettings();

    }

    void CorticalBone::setVerbose(bool verboseIn) {
        verbose=verboseIn;
    }

    bool CorticalBone::reset() {

        initialiseSettings();

        return false;
    }

    /* Setters */
    bool CorticalBone::loadMesh(std::string fileName) {

        // get extension
        std::string extn = fileName.substr(fileName.rfind('.'), fileName.size()-fileName.rfind('.'));

        mesh = vtkSmartPointer<vtkPolyData>::New();


        meshSet = false; meshMeasured = false;
        if( extn.compare(".obj") == 0 ) {
            meshSet = this->openOBJ(fileName, mesh);
        } else if( extn.compare(".wrl") == 0 ) {
            //meshSet = this->openVRML(fileName, mesh);
            meshSet = this->openVRML_SW(fileName, mesh);
        } else if( extn.compare(".stl") == 0 ) {
            meshSet = this->openSTL(fileName, mesh);

        } else if( extn.compare(".ply") == 0 ) {
            meshSet = this->openPLY(fileName, mesh);
        }

        if(!meshSet) { // exit if failed to open mesh
            return meshSet;
        }

        // check for point normals and create if missing
        if(mesh->GetPointData()->GetArray("Normals")!=NULL) {
            meshNormals = mesh->GetPointData()->GetArray("Normals");
        } else if(mesh->GetPointData()->GetNormals()!=NULL) {
            meshNormals = mesh->GetPointData()->GetNormals();
        } else { // create missing normal values
            vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
            normalGenerator->SetInputData(mesh);
            normalGenerator->ComputePointNormalsOn();
            normalGenerator->ComputeCellNormalsOff();
            normalGenerator->Update();
            mesh = normalGenerator->GetOutput();
            meshNormals = mesh->GetPointData()->GetNormals();
        }

        // setup mesh tree
        meshTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
        meshTree->SetDataSet(mesh);
        meshTree->BuildLocator();

        // try load mesh mask
        openMeshMask(fileName);

        if(meshSet && imageSet) {
            // setup model registration
            setSampling();
            setImageBuffer();
            createRegistrationMethods();
        }

        return meshSet;
    }

    bool CorticalBone::loadImage(std::string fileName) {

        // get extension
        std::string extn;
        if(fileName.rfind('.')!=std::string::npos) {
            extn = fileName.substr(fileName.rfind('.'), fileName.size()-fileName.rfind('.'));

        } else { // note for .dcm either return a .dcm file or just the directory
            extn=".dcm";
        }

        Utilities::measureTime(true);

        // get image
        bool status;
        if( extn.compare(".dcm") == 0 ) {
            if(verbose){cout<<Utilities::getTabString()<<"Opened a DICOM CT image."<<endl;}

            std::string filePath = fileName.substr(0, fileName.rfind('/')); // TODO make path separators independent of OS
            status = this->openDICOM(filePath);

        } else if( extn.compare(".mhd") == 0 ) {
            if(verbose){cout<<Utilities::getTabString()<<"Opened a RAW CT image."<<endl;}

            status = this->openRAW(fileName);
            //status = this->openQCT(fileName);

        } else if(extn.compare(".tif") == 0 ) {

            cout<<"Requires a name list when opening .tif image series"<<endl;
            status = false;
            //if(verbose){cout<<Utilities::getTabString()<<"Opened a TIF CT image."<<endl;}

            //this->openTIF(fileName);
        } else if( extn.compare(".QCT") == 0) {
            if(verbose){cout<<Utilities::getTabString()<<"Opening a QCT CT images."<<endl;}

            status = this->openQCT(fileName);

        } else {

            std::string filePath = fileName.substr(0, fileName.rfind('/')); // TODO make path separators independent of OS

            // 1st try open dcm anyway
            if(!openDICOM(filePath)) {
                if(verbose){cout<<Utilities::getTabString()<<"Unsupported image extension encountered:"<<extn<<endl;}
                cerr<<"Unsupported image extension encountered:"<<extn<<endl;
                imageSet=false; meshMeasured = false;
                return false;
            } else {
                status = true;
            }

        }

        if(status) {
            if(verbose){cout<<Utilities::getTabString()<<Utilities::getTabString()<<Utilities::measureTime(false)<<" to load the image"<<endl;}
        } else {
            if(verbose){cout<<Utilities::getTabString()<<Utilities::getTabString()<<"Error in CorticalBone::loadImage() while opening a file. No file opened."<<endl;}
            return false;
        }

        return setupImage();
    }

    bool CorticalBone::loadImage(std::vector<std::string> nameList, double zSpacing, double zOffset) {

        // get extension
        std::string extn;
        if(nameList.at(0).rfind('.')!=std::string::npos) {
            extn = nameList.at(0).substr(nameList.at(0).rfind('.'), nameList.at(0).size()-nameList.at(0).rfind('.'));
        } else {
            extn=".dcm";
        }

        Utilities::measureTime(true);

        // get image
        if(extn.compare(".tif") == 0 ) {
            if(verbose){cout<<Utilities::getTabString()<<"Opened a TIF CT image: "<<nameList.at(0).c_str()<<endl;}
            this->openTIF(nameList, zSpacing, zOffset);
        } else {
            if(verbose){cout<<Utilities::getTabString()<<"Unsupported image extension encountered for opening with prespecified name list:"<<extn<<endl;}
            imageSet=false; meshMeasured = false;
            return false;
        }

        if(verbose){cout<<Utilities::getTabString()<<Utilities::getTabString()<<Utilities::measureTime(false)<<" to load the image"<<endl;}
        cerr<<Utilities::measureTime(false)<<" to load the image"<<endl;


        return setupImage();
    }

    bool CorticalBone::setImportProfiles(std::string fileName) {

        clearPreviousResults(); parametersImported = false;

        std::string fileNameStub = fileName.substr(0, fileName.rfind('_'));
        std::string extn = ".txt";

        if(verbose){cout<<"imported paths"<<Utilities::getTabString()<<"- parameters: "<<fileNameStub+"_Parameters"+extn<<endl;}
        if(verbose){cout<<Utilities::getTabString()<<"- profiles: "<<fileNameStub+"_ImageProfiles"+extn<<endl;}

        importedParametersTable = vtkSmartPointer<vtkTable>::New();
        if (!openImportParameters(fileNameStub+"_Parameters"+extn)) { // set importedModelIndex to whatever is specified in the file
            if(verbose){cout<<Utilities::getTabString()<<"Error: while importing the parameters file."<<endl;}
            return parametersImported;
        }

        importedPositionsArray = vtkSmartPointer<vtkDoubleArray>::New();
        if(!openProfiles(fileNameStub + "_Positions" + extn, importedPositionsArray, "Positions",
                         importedModelIndex)) { // add rows to imported profile tables
            if(verbose){cout<<Utilities::getTabString()<<"Error: while importing the positions file."<<endl;}
            return parametersImported;
        }

        importedProfilesArray = vtkSmartPointer<vtkDoubleArray>::New();
        if(!openProfiles(fileNameStub + "_ImageProfiles" + extn, importedProfilesArray, "Image Profiles",
                         importedModelIndex)) { // add rows to imported profile tables
            if(verbose){cout<<Utilities::getTabString()<<"Error: while importing the image profiles file."<<endl;}
            return parametersImported;
        }

        if(importedModelIndex==kHighResClassifier) {
            importedClassificationArray = vtkSmartPointer<vtkDoubleArray>::New();
            importedPercentageArray = vtkSmartPointer<vtkDoubleArray>::New();

            std::string classificationName = fileNameStub+"_ClassificationProfiles"+extn;
            std::string percentageName = fileNameStub+"_PercentageProfiles"+extn;
            ifstream classificationFile (classificationName.c_str()), percentageFile (percentageName.c_str());
            bool status = classificationFile.is_open() && percentageFile.is_open();
            if(status) {
                status = openProfiles(classificationName, importedClassificationArray, "Classification Profiles", importedModelIndex);
            }
            if(status) {
                status = openProfiles(percentageName, importedPercentageArray, "Percentage Profiles", importedModelIndex);
            }
            if(!status) {
                vtkIdType maxI = mesh->GetNumberOfPoints(), maxJ = numberOfSamples;
                importedClassificationArray->SetNumberOfComponents(numberOfSamples);
                importedClassificationArray->SetNumberOfTuples(mesh->GetNumberOfPoints());

                importedPercentageArray->SetNumberOfComponents(numberOfSamples);
                importedPercentageArray->SetNumberOfTuples(mesh->GetNumberOfPoints());
                for(int i=0; i<maxI; i++) {
                    for(int j=0; j<maxJ; j++) {
                        importedClassificationArray->SetComponent(i,j, nan("1"));
                        importedPercentageArray->SetComponent(i,j, nan("1"));
                    }
                }
            }
        }

        // hack to set fabricated data to having an ramp model
        if(importedModelIndex == kFabricated) { importedModelIndex= kEndostealRamp; }

        clearPreviousResults(); parametersImported = true; // update the display values to include enteries for the import errors
        return parametersImported;
    }

    bool CorticalBone::setCalibration(double p0, double p1, double p2) {
        if(phantomIndex!=kNoCal && !isnan(calP0) && !isnan(calP1) && !isnan(calP2)) {
            calP0=p0; calP1=p1; calP2=p2; calibrationSet=true; phantomChanged = false;

            if(verbose){std::cout<<"Calibration (Manually set) [BMD = p0 + p1 x HU + p2 * HU<sup>2</sup>]: p0 = "<< calP0 <<", p1 = "<<
                        calP1 << ", p2 = "<< calP2 <<std::endl;}

            // todo - only if kCal set to manual? if not check for conflicts. set values must match those that would have been caluculated.
        } else {
            calP0=0.0; calP1=1.0; calP2=0.0; calibrationSet=false; phantomChanged = false;
        }
        modelProfile->SetCalibration(calP0,calP1,calP2);

        // reset results
        clearPreviousResults();

        return calibrationSet;
    }

    bool CorticalBone::setCalibrationPtGeometry(double radius, int number) {
        if(isnan(radius)||isinf(radius)||number<0) {
            ImageType::SpacingType imageSpacing = image->GetSpacing();
            cntrlRadius = 10 * (imageSpacing[0] + imageSpacing[1] + imageSpacing[2]);
            cntrlNumber = 0; calNumber = 0; return false;
        } else {
            cntrlRadius = radius; cntrlNumber =number; calNumber = number; // asumption is this is the manula ctrl mode
            return true;
        }
    }

    bool CorticalBone::getCalibrationPtGeometry(double &radius, int &number) {

        radius= cntrlRadius;
        number= (int)cntrlNumber;
        return !(isnan(cntrlRadius) || isinf(cntrlRadius) || cntrlNumber == 0);

    }

    void CorticalBone::setDisplayIndex(int index) {

        if(index==kInvalid) {
            displayIndex=kInvalid;
            return;
        }

        int maxNumbeOfDisplays = getNumberOfDisplays();

        if(index >= 0 && index < maxNumbeOfDisplays) {
            displayIndex = index;
        } else { // set blank
            std::cerr<<"Error - display index entered is out of bounds: "<<index<<std::endl;
            displayIndex = kInvalid;
        }

    }

    int CorticalBone::getDisplayIndex() {
        return displayIndex;
    }

    void CorticalBone::setPhantomIndex(int index) {

        if(index==kManualControlPtsCal) {
            cntrlNumber = calNumber = 0; cntrlRadius =nan("1"); // undefined until set
        } else if(index==kMindwaySolidCal) {
            cntrlNumber = 4; calNumber = 5;
            cntrlRadius = 6.0;
        } else if(index==kBoneDensityCal) {
            cntrlNumber = 4; calNumber = 3;
            cntrlRadius = 6.0;
        } else if(index==kEuropeanSpineCal) {
            cntrlNumber = 12; calNumber = 12;
            cntrlRadius = 2.0; // 10*(imageSpacing[0]+imageSpacing[1]+imageSpacing[2])
        } else if(index==kNoCal) {
            cntrlNumber = calNumber = 0;
            ImageType::SpacingType imageSpacing = image->GetSpacing();
            cntrlRadius = nan("1");
        } else if(index==kManualLinearCal || index==kManualQuadraticCal) {
            cntrlNumber = calNumber = 0; cntrlRadius = nan("1"); // not relevant
        }

        if(index > kNoCal && index <= kManualQuadraticCal) {
            phantomIndex = index;
            if(calibrationSet) {
                clearCalibration();
            }
            if(meshMeasured) {
                clearPreviousResults();
            }
            phantomChanged = true;

        } else if(kNoCal == phantomIndex) {
            phantomChanged=false;
        } else {
            std::cerr<<"Error - phantom index entered is out of bounds: "<<index<<std::endl;
            phantomIndex=kNoCal; phantomChanged=false;
        }
    }

    void CorticalBone::setModelIndex(int index) {

        if(index < kThreeTierRect || index > kCalibration) {
            std::cerr<<"Error - model index entered is out of bounds: "<<index<<std::endl;
            return;
        }

        modelIndex = index;
        clearPreviousResults();

        if(modelIndex<=kEndostealRamp) {
            modelProfile->TurnOffAllSamplesRequired();
            modelProfile->TurnOnMeanOnly();
            modelProfile->TurnOffMaxDetection();
        } else if(modelIndex==kCalibration) {
            modelProfile->TurnOnAllSamplesRequired();
            modelProfile->TurnOnMeanOnly();
            modelProfile->TurnOnMaxDetection();
        } else {
            modelProfile->TurnOffAllSamplesRequired();
            modelProfile->TurnOffMeanOnly();
            modelProfile->TurnOnMaxDetection();
        }
    }

    void CorticalBone::setOptimiserIndex(int index) {
        if(index >= kLMOptimiser && index <= kEvolutionaryOptimiser) {
            registrationMethod->SetOptimiserSelection(index);
            clearPreviousResults();
        } else {
            registrationMethod->SetOptimiserSelection(kLMOptimiser);
            std::cerr<<"Error - optimiser index entered is out of bounds: "<<index<<std::endl;
        }
    }

    void CorticalBone::setFittingSchemeIndex(int index) {


        if(modelIndex <= kEndostealRamp && index >= kStdAFitting && index <= kCBSmoothingDFitting) {
            fittingSchemeIndex = index;

        } else if(modelIndex == kHighResClassifier && index >= itk::ClassifierTransform::kGlobal && index <= itk::ClassifierTransform::kMedian) {
            fittingSchemeIndex = index;
            classifierMethod->SetModelScheme(fittingSchemeIndex); // make so classifier index = 0

        } else if(modelIndex == kCalibration && index >= itk::ProfileCalibration::kMedian && index <= itk::ProfileCalibration::kMaximum) {
            fittingSchemeIndex = index;
            calibrationMethod->SetModelScheme(fittingSchemeIndex);
        } else {
            std::cerr<<"Error - fitting scheme index entered is out of bounds: "<<index<<std::endl;
            fittingSchemeIndex = -1;
        }

        clearPreviousResults();
        resetCBConstraint();

    }

    double CorticalBone::getProfileLength () {
        return profileLength;
    }

    // model info
    void CorticalBone::getModelSelection(int& modelIndexIn, std::string& modelName) {
        modelIndexIn=modelIndex;
        modelName = getModelName();
    }

    void CorticalBone::getOptimiserSelection(int& optimiserIndexIn, std::string& optimiserName){
        optimiserIndexIn=registrationMethod->GetOptimiserSelection();
        optimiserName=getOptimisationName();
    }

    void CorticalBone::getSchemeSelection(int& schemeIndexIn, std::string& schemeName){
        schemeIndexIn=fittingSchemeIndex;
        if(modelIndex<=kEndostealRamp) {
            if(fittingSchemeIndex== kStdAFitting) {
                schemeName="Std A Fitting";
            } else if(fittingSchemeIndex== kStdBFitting) {
                schemeName="Std B Fitting";
            } else if(fittingSchemeIndex== kStdCFitting) {
                schemeName="Std C Fitting";
            } else if(fittingSchemeIndex== kStdDFitting) {
                schemeName="Std D Fitting";
            } else if(fittingSchemeIndex== kCBMV2AFitting) {
                schemeName="CBMV2AFitting";
            } else if(fittingSchemeIndex== kCBMV2BFitting) {
                schemeName="CBMV2BFitting";
            } else if(fittingSchemeIndex== kCBMV2CFitting) {
                schemeName="CBMV2CFitting";
            } else if(fittingSchemeIndex== kCBMV2DFitting) {
                schemeName="CBMV2DFitting";
            } else if(fittingSchemeIndex== kUnconstrainedAFitting) {
                schemeName="UnconstrainedAFitting";
            } else if(fittingSchemeIndex== kUnconstrainedBFitting) {
                schemeName="UnconstrainedBFitting";
            } else if(fittingSchemeIndex== kUnconstrainedCFitting) {
                schemeName="UnconstrainedCFitting";
            } else if(fittingSchemeIndex== kCBSmoothingAFitting) {
                schemeName="CBSmoothingAFitting";
            } else if(fittingSchemeIndex== kCBSmoothingBFitting) {
                schemeName="CBSmoothingBFitting";
            } else if(fittingSchemeIndex== kCBSmoothingCFitting) {
                schemeName="CBSmoothingCFitting";
            } else if(fittingSchemeIndex== kCBSmoothingDFitting) {
                schemeName="CBSmoothingDFitting";
            } else {
                schemeName="Invalid";
            }
        } else {
            schemeName = classifierMethod->GetModelSchemeName();
        }
    }

    // CB density Mode
    void CorticalBone::setFixedCBDensity(double fixedCBDensity, bool globalFixedDensity) {

        // calibration correction to avoid calibration twice
        //fixedCBDensity = fixedCBDensity * this->calibrationSlope + this->calibrationIntercept; - no longer necessary


        if(globalFixedDensity) {
            globalFixedCB = fixedCBDensity;
            if(globalFixedCB>maxImageValue) { // update max density if grteater than image
                registrationMethod->SetMaxBMDDensity(globalFixedCB); // todo - consider removeing this and leaving it as fixed to the calibrated 1600 value
            }
        }

        registrationMethod->FixParameter(fixedCBDensity, ModelRegistrationType::kYcbParam);

        // add value to 'flexParameterArray' - unless already added
        if(transformFlexMap[ModelRegistrationType::kYcbParam]==1.0) {

            // Note: rect and ramp can be the same as ModelRegistrationType::kCorticalDensityParam index is less than divergance

            int flexIndex=0; // get index to remove
            for(int i=0; i<ModelRegistrationType::kYcbParam; i++) {
                flexIndex += transformFlexMap[i];
            }

            transformFlexMap[ModelRegistrationType::kYcbParam]=0.0;

            // remove value
            removeParametersValue(rectFlexParameters, flexIndex);
            removeParametersValue(rectFlexScales, flexIndex);
            removeParametersValue(rampFlexScales, flexIndex);

        }

        if(globalFixedDensity) {
            clearPreviousResults();
        }

    }

    void CorticalBone::removeFixedCBDensity(bool clearGlobalFixedDensity) {

        // set map value 0, remove fixed value, insert flex value
        registrationMethod->FreeParameter(ModelRegistrationType::kYcbParam);

        // reinsert value into 'flexParameterArray' - unless not removed
        if(transformFlexMap[ModelRegistrationType::kYcbParam]==0.0) {

            // Note: rect and ramp can be the same as ModelRegistrationType::kCorticalDensityParam index is less than divergance

            int flexIndex=0; // get index to remove
            for(int i=0; i<ModelRegistrationType::kYcbParam; i++) {
                flexIndex += transformFlexMap[i];
            }

            transformFlexMap[ModelRegistrationType::kYcbParam]=1.0;

            // insert value
            insertParametersValue(rectFlexParameters, rectAllParameters[ModelRegistrationType::kYcbParam], flexIndex);
            insertParametersValue(rectFlexScales, rectAllScales[ModelRegistrationType::kYcbParam], flexIndex);
            insertParametersValue(rampFlexScales, rampAllScales[ModelRegistrationType::kYcbParam], flexIndex);

        }

        if(clearGlobalFixedDensity) {
            clearPreviousResults();
            globalFixedCB = nan("1");
            registrationMethod->SetMaxBMDDensity(maxBMDValue); // maxImageValue
        }

    }

    void CorticalBone::turnOnFWHMMode() {
        registrationMethod->TurnOnFWHMMode();

        setFixedCBDensity(registrationMethod->GetMaxDensity());
        clearPreviousResults();
    }

    void CorticalBone::turnOffFWHMMode() {
        // set registration methods to normal mode
        registrationMethod->TurnOffFWHMMode();

        clearPreviousResults();
    }

    void CorticalBone::turnOffFixedSigma() {
        modelProfile->TurnOffGlobalSigma();
        removeFixedSigma();
        clearPreviousResults();
    }

    void CorticalBone::turnOnFixedSigma(double x, double y, double z){
        double sigma[Dimension] = {x, y, z};
        modelProfile->TurnOnGlobalSigma(sigma);
        clearPreviousResults();
    }

    // profile averaging
    void CorticalBone::turnOffProfileAveraging() {
        modelProfile->TurnOffMultipleProfiles();
        clearPreviousResults();
    }

    void CorticalBone::turnOnProfileAveraging(double x, double y, double z) {
        double hrFWHMValues[Dimension] = {x, y, z};
        modelProfile->TurnOnMultipleProfiles(hrFWHMValues);
        clearPreviousResults();
    }

    bool CorticalBone::getProfileAveragingValues(double &x, double &y, double &z) {
        return modelProfile->GetRadius(x, y, z);
    }

    // fixed sample number
    void CorticalBone::removeFixedSampleMode() {
        sampleNumberFixed = false;
        if(meshSet==true && imageSet==true) {
            setSampling();
            modelProfile->SetNumberOfSamples(numberOfSamples);

            clearPreviousResults();
        }
    }

    bool CorticalBone::setFixedSampleNumber(int sampleNumber) {
        sampleNumberFixed = true;
        numberOfSamples = sampleNumber;
        modelProfile->SetNumberOfSamples(numberOfSamples);

        clearImageProfilesArray();

        clearPreviousResults();

        if(verbose){cout<<Utilities::getTabString()<<"Number of Samples Set Manually: "<<numberOfSamples<<endl;}

        return true; // return false if invalid sampleNumber
    }

    int CorticalBone::getSampleNumber() {
        return numberOfSamples;
    }

    // smoothing value
    void CorticalBone::setSmoothingRadius(double smoothingValue) {
        smoothingRadius = smoothingValue;

        clearPreviousResults();

    }

    void CorticalBone::turnOffSmoothingMode() {

        smoothingSet = false;
        smoothingRadius = nan("1");

        clearPreviousResults();
    }

    void CorticalBone::turnOnSmoothingMode() {
        smoothingSet = true;

        clearPreviousResults();
    }

    double CorticalBone::getSmoothingValue() {
        return smoothingRadius;
    }


    // import profiles
    void CorticalBone::removeImportedProfile() {
        parametersImported = false;
        importedProfilesArray = NULL; importedPositionsArray = NULL;
        importedClassificationArray = NULL;
        importedPercentageArray = NULL; importedParametersTable = NULL;

        clearPreviousResults();

    }



    /* Savers */
    bool CorticalBone::saveMesh(std::string fileName) {
        bool status = false;

        // get extension
        std::string extn = fileName.substr(fileName.rfind('.'), fileName.size()-fileName.rfind('.'));

        if( extn.compare(".obj") == 0 ) {
            status = this->saveOBJ(fileName, mesh);
        } else if( extn.compare(".wrl") == 0 ) {
            status = this->saveVRML(fileName, mesh);
        } else if( extn.compare(".stl") == 0 ) {
            status = this->saveSTL(fileName, mesh);
        }

        return status;

    }

    bool CorticalBone::saveValueArrays(std::string baseFileName) {

        if(!meshMeasured) {
            cout<<"Warning in CorticalBone::saveValueArrays() no results to save."<<endl;
            return false;
        }

        // get extension
        std::string extn = ".txt"; std::string objExtn = ".obj"; std::string wrlExtn = ".wrl"; std::string projExtn = ".xml";
        std::string prameterLabel = "_Parameters"+extn;
        std::string displayLabel = "_Displays"+extn;
        std::string imageProfileLabel = "_ImageProfiles"+extn, imagePositionsLabel = "_Positions"+extn;
        std::string periostealMeshLabel = "_PeriostealMesh";

        // trim off any extns or parameter/display/imageprofile/periostealmesh labels
        baseFileName = baseFileName.substr(0, baseFileName.rfind(prameterLabel));
        baseFileName = baseFileName.substr(0, baseFileName.rfind(displayLabel));
        baseFileName = baseFileName.substr(0, baseFileName.rfind(imageProfileLabel));
        baseFileName = baseFileName.substr(0, baseFileName.rfind(periostealMeshLabel));
        baseFileName = baseFileName.substr(0, baseFileName.rfind(extn)); // as mac auto adds .txt extn
        baseFileName = baseFileName.substr(0, baseFileName.rfind(projExtn)); // just incase called with project name

        saveTable(baseFileName+prameterLabel, modelParametersTable, "Parameter Values");
        saveTable(baseFileName+displayLabel, modelDisplayTable, "Display Values");
        saveProfiles(baseFileName + imageProfileLabel, imageProfilesArray, "Image Profiles");
        saveProfiles(baseFileName + imagePositionsLabel, imagePositionsArray, "Positions");
        if(modelIndex!=kCalibration) {
            saveOBJ(baseFileName+periostealMeshLabel+objExtn, periostealMesh);
            saveVRML(baseFileName+periostealMeshLabel+wrlExtn, periostealMesh);
        }

        if(modelIndex==kHighResClassifier) {
            std::string classificationLabel = "_ClassificationProfiles"+extn, percentageLabel = "_PercentageProfiles"+extn;
            saveProfiles(baseFileName + classificationLabel, imageClassificationArray, "Classification Profiles");
            saveProfiles(baseFileName + percentageLabel, imagePercentageArray, "Percentage Profiles");
        }

        return true;
    }

    bool CorticalBone::saveCalibration(std::string baseFileName) {


        if(phantomIndex==kManualLinearCal || phantomIndex==kManualQuadraticCal) {
            //cout<<"Warning in CorticalBone::saveCalibration() no control points to save."<<endl;
            return false;
        } else if(!calibrationSet) {
            //cout<<"Warning in CorticalBone::saveCalibration() no results to save."<<endl;
            return false;
        }

        // get extension
        std::string extn = ".txt"; std::string projExtn = ".xml";
        std::string prameterLabel = "_Parameters"+extn;
        std::string displayLabel = "_Displays"+extn;
        std::string imageProfileLabel = "_ImageProfiles"+extn, imagePositionsLabel = "_Positions"+extn;
        std::string periostealMeshLabel = "_PeriostealMesh";

        // trim off any extns or parameter/display/imageprofile/periostealmesh labels
        baseFileName = baseFileName.substr(0, baseFileName.rfind(prameterLabel));
        baseFileName = baseFileName.substr(0, baseFileName.rfind(displayLabel));
        baseFileName = baseFileName.substr(0, baseFileName.rfind(imageProfileLabel));
        baseFileName = baseFileName.substr(0, baseFileName.rfind(periostealMeshLabel));
        baseFileName = baseFileName.substr(0, baseFileName.rfind(extn)); // as mac auto adds .txt extn
        baseFileName = baseFileName.substr(0, baseFileName.rfind(projExtn)); // just incase called with project name

        // get extension
        std::string calibrationLabel = "_Calibration";
        std::string calibrationHeader1, calibrationHeader2;


        // create array of data to store
        vtkSmartPointer<vtkDoubleArray> calibrationArray = vtkSmartPointer<vtkDoubleArray>::New();
        if(phantomIndex==kManualControlPtsCal) {
            calibrationHeader1 ="Manual Calibration Control Points";
            calibrationHeader2="x, y, z, HU";
            calibrationArray->SetNumberOfComponents(4);
            calibrationArray->SetNumberOfTuples(cntrlNumber);
            for(vtkIdType i=0; i<cntrlNumber; i++) {
                calibrationArray->SetComponent(i, 0, cntrlPts->GetComponent(i, 0));
                calibrationArray->SetComponent(i, 1, cntrlPts->GetComponent(i, 1));
                calibrationArray->SetComponent(i, 2, cntrlPts->GetComponent(i, 2));
                calibrationArray->SetComponent(i, 3, cntrlVals->GetValue(i));
            }
        } else {
            calibrationHeader1 ="Manual Calibration Control Points";
            calibrationHeader2="HU, BMD";
            int number = calBMDVals->GetNumberOfTuples();

            calibrationArray->SetNumberOfComponents(2);
            calibrationArray->SetNumberOfTuples(number);
            for(vtkIdType i=0; i<number; i++) {
                calibrationArray->SetComponent(i, 0, calHUVals->GetValue(i));
                calibrationArray->SetComponent(i, 1, calBMDVals->GetValue(i));
            }
        }


        // construct header information
        std::string delimiter = ",";
        calibrationHeader1.append(delimiter).append(getPhantomName(phantomIndex));

        saveArray(baseFileName + calibrationLabel + extn, calibrationArray, calibrationHeader1, calibrationHeader2);

        return true;
    }

    /* Loaders */
    bool CorticalBone::loadValueArrays(std::string fileName) {

        clearPreviousResults(); // setup parameters

        std::string fileNameStub = fileName.substr(0, fileName.rfind('.'));
        std::string extn = ".txt"; std::string meshExtn = ".obj";

        if (!openParameters(fileNameStub+"_Parameters"+extn)) {
            clearPreviousResults();
            return false;
        }

        if (!openDisplays(fileNameStub+"_Displays"+extn)) {
            clearPreviousResults();
            return false;
        }

        if (!openProfiles(fileNameStub + "_Positions" + extn, imagePositionsArray, "Positions", modelIndex)) {
            clearPreviousResults();
            return false;
        }

        if (!openProfiles(fileNameStub + "_ImageProfiles" + extn, imageProfilesArray, "Image Profiles", modelIndex)) {
            clearPreviousResults();
            return false;
        }

        if(modelIndex==kHighResClassifier) { // todo simplify
            std::string classificationName = fileNameStub+"_ClassificationProfiles"+extn;
            std::string percentageName = fileNameStub+"_PercentageProfiles"+extn;
            ifstream classificationFile (classificationName.c_str()), percentageFile (percentageName.c_str());
            bool status = classificationFile.is_open() && percentageFile.is_open();
            if(status) {
                status = openProfiles(classificationName, imageClassificationArray, "Classification Profiles",
                                      modelIndex);
            }
            if(status) {
                status = openProfiles(percentageName, imageClassificationArray, "Percentage Profiles", modelIndex);
            }
            if(!status) {
                vtkIdType maxI = mesh->GetNumberOfPoints(), maxJ = numberOfSamples;
                imageClassificationArray = vtkSmartPointer<vtkDoubleArray>::New();
                imageClassificationArray->SetNumberOfComponents(numberOfSamples);
                imageClassificationArray->SetNumberOfTuples(mesh->GetNumberOfPoints());

                imagePercentageArray = vtkSmartPointer<vtkDoubleArray>::New();
                imagePercentageArray->SetNumberOfComponents(numberOfSamples);
                imagePercentageArray->SetNumberOfTuples(mesh->GetNumberOfPoints());
                for(int i=0; i<maxI; i++) {
                    for(int j=0; j<maxJ; j++) {
                        imageClassificationArray->SetComponent(i,j, nan("1"));
                        imagePercentageArray->SetComponent(i,j, nan("1"));
                    }
                }
            }
        }

        modelFitStatusArray.Fill(kValid); // todo tidy up

        maskResultArrays(); // sets modelFitStatusArray[i]==kMasked; todo tidy up

        // set the validity array
        int maxI = modelDisplayTable->GetNumberOfRows();
        for(int i=0; i<maxI; i++) {
            if(!(meshMask[i]==kMasked) && isnan(modelDisplayTable->GetValue(i,0).ToDouble())) { // check if out of bounds / invalid

                bool outOfBouunds = true;
                for(int j=0; j<numberOfSamples; j++) {
                    if(!isnan(imageProfilesArray->GetComponent(i, j))) {
                        outOfBouunds=false; break;
                    }
                }
                if(outOfBouunds) {
                    modelFitStatusArray[i] = kOutOfBounds;
                    outOfBoundsCount++;
                } else {
                    modelFitStatusArray[i] = kInvalid;
                    nanCount++;
                }
            }

        }
        meshMeasured = true;
        
        // try load periosteal mesh - generate if missing
        ifstream profilesFile (fileName.c_str());
        if(profilesFile.good()) { // open if it already exists
            periostealMesh = vtkSmartPointer<vtkPolyData>::New();
            openOBJ(fileNameStub+"_PeriostealMesh"+meshExtn, periostealMesh); // todo - consider requiring a saved mesh in the future
        } else { // generate a new one if it doesn't
            calculatePeriostealMesh();
        }

        loadedResults = meshMeasured;

        return meshMeasured;
    }

    /* Getters */
    vtkSmartPointer<vtkPolyData> CorticalBone::getMesh() {
        return mesh;
    }

    vtkSmartPointer<vtkPolyData> CorticalBone::getPeriostealMesh() {
        if(meshMeasured) {
            return periostealMesh;
        } else {
            cout<<"Warning in CorticalBone::getPeriostealMesh() periostealMesh doesn't exist, so mesh returned instead."<<endl;
            return mesh;
        }
    }

    ImageType::Pointer CorticalBone::getDICOM() {
        return image;
    }

    vtkSmartPointer<vtkTable> CorticalBone::getImageProfileTable(vtkIdType pointId) {

        vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

        if(modelFitStatusArray[pointId]==kOutOfBounds || modelFitStatusArray[pointId]==kMasked) { // empty table
            return table;
        }

        if(modelIndex<=kEndostealRamp) { // a fitter model
            table = registrationMethod->GetImageDataTable(); // positions are offset for some reason
        } else if(modelIndex==kHighResClassifier) { // a classifier / thresholder model
            table = classifierMethod->GetImageDataTable();
        } else if(modelIndex==kCalibration) {
            table = calibrationMethod->GetImageDataTable();
        }

        return table;
    }

    vtkSmartPointer<vtkTable> CorticalBone::getImportedProfileTable(vtkIdType pointId) {

        if(!parametersImported || modelIndex==kCalibration) { // out of bounds - return a empty table
            vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
            return table;
        }

        // create table
        vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
        table->SetNumberOfRows(numberOfSamples);

        // create position and image profiles to add to the table
        vtkSmartPointer<vtkDoubleArray> values = vtkSmartPointer<vtkDoubleArray>::New();
        values->SetNumberOfValues(numberOfSamples); values->SetName("Imported Image");

        vtkSmartPointer<vtkDoubleArray> positions = vtkSmartPointer<vtkDoubleArray>::New();
        positions->SetName("Positions"); positions->SetNumberOfValues(numberOfSamples);

        double offset = calculateImportOffset(pointId);
        for(int i=0; i<numberOfSamples; i++) {
            positions->SetValue(i, importedPositionsArray->GetComponent(pointId, i)-offset);
            values->SetValue(i, importedProfilesArray->GetComponent(pointId, i));
        }

        // add to table
        table->AddColumn(positions); table->AddColumn(values);

        if (importedModelIndex==kHighResClassifier) { // a classifier / thresholder model
            vtkSmartPointer<vtkDoubleArray> classifications = vtkSmartPointer<vtkDoubleArray>::New();
            classifications->SetNumberOfValues(numberOfSamples); classifications->SetName("Imported Classifications");
            vtkSmartPointer<vtkDoubleArray> percentages = vtkSmartPointer<vtkDoubleArray>::New();
            percentages->SetNumberOfValues(numberOfSamples); percentages->SetName("Imported Percentages");
            for(int i=0; i<numberOfSamples; i++) {
                classifications->SetValue(i, importedClassificationArray->GetComponent(pointId, i));
                percentages->SetValue(i, importedPercentageArray->GetComponent(pointId, i)*10);
            }
            table->AddColumn(classifications); table->AddColumn(percentages);
        }

        return table;
    }

    vtkSmartPointer<vtkTable> CorticalBone::getDisplayModelTable(vtkIdType pointId) {

        vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

        if(modelFitStatusArray[pointId]==kOutOfBounds || modelFitStatusArray[pointId]==kMasked) { // empty table
            return table;
        }

        if(modelIndex<=kEndostealRamp) { // a fitter model
            table = registrationMethod->GetDisplayModelTable();
        } else if (modelIndex==kHighResClassifier) { // a classifier / thresholder model
            table = classifierMethod->GetDisplayModelTable();
        } else if(modelIndex==kCalibration) {
            table = calibrationMethod->GetDisplayModelTable();
        }

        //writeModelParameterToOutput(pointId, modelParametersTable, modelIndex); // output model values only display also updated. todo - why does it just give zero and nan values?
        return table;
    }

    vtkSmartPointer<vtkTable> CorticalBone::getWeightsTable(vtkIdType pointId) {

        vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

        if(modelFitStatusArray[pointId]==kOutOfBounds || modelFitStatusArray[pointId]==kMasked) { // empty table
            return table;
        }

        if(modelIndex<=kEndostealRamp) { // a fitter model
            table = registrationMethod->GetWeightTable();
        }

        //writeModelParameterToOutput(pointId, modelParametersTable, modelIndex); // output model values only display also updated. todo - why does it just give zero and nan values?
        return table;
    }

    vtkSmartPointer<vtkTable> CorticalBone::getImportedDisplayModelTable(vtkIdType pointId) {

        vtkSmartPointer<vtkTable> table;

        if(!parametersImported || modelIndex==kCalibration) { // out of bounds - return a empty table
            table = vtkSmartPointer<vtkTable>::New();
            return table;
        }

        int n = importedParametersTable->GetNumberOfColumns();

        double offset = calculateImportOffset(pointId);
        if(importedModelIndex<=kEndostealRamp) { // a fitter model

            ParametersType parameters(n);
            for(int i=0;i<n;i++) {
                parameters[i] = importedParametersTable->GetValue(pointId, i).ToDouble();
            }

            table = registrationMethod->GetDisplayModelTable(parameters, importedModelIndex, offset);

        } else if(importedModelIndex==kHighResClassifier) { // a classifier / threshold model

            vtkSmartPointer<vtkDoubleArray> parameters = vtkSmartPointer<vtkDoubleArray>::New();
            parameters->SetNumberOfValues(n); parameters->SetName("Imported Model parameters");
            for(int i=0; i<n; i++) {
                parameters->SetValue(i, importedParametersTable->GetValue(pointId, i).ToDouble());
            }
            table = classifierMethod->GetDisplayModelTable(parameters, offset);
        } else {
            cerr<<"Invalid model import index in CorticalBone::getImportedDisplayModelTable = "<<importedModelIndex<<endl;
            return table;
        }

        if(verbose){cout<<"Imported "; writeModelParameterToOutput(pointId, importedParametersTable, importedModelIndex);}
        return table;
    }

    vtkSmartPointer<vtkTable> CorticalBone::getProcessingModelTable(vtkIdType pointId) {

        if(modelFitStatusArray[pointId]==kOutOfBounds || modelFitStatusArray[pointId]==kMasked) { // empty table
            vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
            return table;
        }

        if(modelIndex<=kEndostealRamp) { // a fitter model
            return registrationMethod->GetProcessorModelTable();
        } else if(modelIndex==kHighResClassifier) { // a classifier / thresholder model
            return classifierMethod->GetProcessingModelTable();
        } else if(modelIndex==kCalibration) {
            return calibrationMethod->GetProcessingModelTable();
        } else { // empty table
            cerr<<"error in CorticalBone::getProcessingModelTable() invalid modelIndex="<<modelIndex<<endl;
            vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
            return table;
        }
    }

    vtkSmartPointer<vtkTable> CorticalBone::getImportedProcessingModelTable(vtkIdType pointId) {
        vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
        if(importedModelIndex<=kEndostealRamp) { // a fitter model

        } else { // a classifier / thresholder model

        }
        return table;
    }

    vtkSmartPointer<vtkTable> CorticalBone::getErrorAreaTable(vtkIdType pointId) {

        vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

//    if(parametersImported==true) {
//      
//      double offset = calculateImportOffset(pointID);
//      double sampleSpacing = profileLength/numberOfSamples; 
//      int offsetIndex = ceil(fabs(offset)/sampleSpacing), numberOfOverlap = numberOfSamples - offsetIndex;
//      
//      vtkSmartPointer<vtkDoubleArray> localModel = registrationMethod->GetUnblurredModelValues();
//      
//      vtkSmartPointer<vtkDoubleArray> overlapX = vtkSmartPointer<vtkDoubleArray>::New(); 
//      vtkSmartPointer<vtkDoubleArray> overlapMinY = vtkSmartPointer<vtkDoubleArray>::New(); 
//      vtkSmartPointer<vtkDoubleArray> overlapMaxY = vtkSmartPointer<vtkDoubleArray>::New(); 
//      overlapX->SetName("Position"); overlapMinY->SetName("Min Y"); overlapMaxY->SetName("Max Y");
//      overlapX->SetNumberOfValues(numberOfOverlap); overlapMinY->SetNumberOfValues(numberOfOverlap); overlapMaxY->SetNumberOfValues(numberOfOverlap);
//    
//      for (int i=0; i<numberOfOverlap; i++) {
//        int localIndex = (importedPositionOffset>=0) ? i : i + offsetIndex;
//        int importedIndex = (importedPositionOffset>=0) ? i + offsetIndex : i;
//        overlapX->SetValue(i, positionsArray->GetValue(localIndex));
//        double localValue = localModel->GetValue(localIndex); // get localModel value at an arbitrary location
//        
//        double importedValue = importedProfilesArray->GetComponent(pointID, importedIndex);
//        
//        
//        double minY = std::min(localValue, importedValue), maxY = std::max(localValue, importedValue);
//        overlapMinY->SetValue(i, minY); overlapMaxY->SetValue(i, maxY);
//      }
//      
//      table->SetNumberOfRows(numberOfOverlap);
//      table->AddColumn(overlapX);
//      table->AddColumn(overlapMinY);
//      table->AddColumn(overlapMaxY);
//    }

        return table;
    }

    bool CorticalBone::getCalibrationPoints(int &phantomType, vtkSmartPointer<vtkDoubleArray> &pts) {
        pts = cntrlPts; phantomType = phantomIndex;
        return calibrationSet;
    }

    bool CorticalBone::getCalibrationValues(double& p0, double& p1, double& p2) {
        if(calibrationSet == true) {
            p1 = calP1;
            p0 = calP0;
            p2 = calP2;
        } else {
            p0 = 0.0; p1 = 1.0; p2 = 0.0;
        }
        return calibrationSet;
    }

    bool CorticalBone::getCalibrationPhantom(std::string& phantomType, int& phantomIndexIn) {
        phantomType = getPhantomName(phantomIndex);
        phantomIndexIn = phantomIndex;
        return calibrationSet;
    }

    int CorticalBone::getState() {
        if(imageSet==false && meshSet==false) {
            return 0; // none open
        } else if(imageSet==false && meshSet==true) {
            return 1; // mesh open
        } else if (phantomChanged) {
            return 2; // calibration phantom changed and calibration not rerun
        } else if(imageSet==true && meshSet==false) {
            return 2; // image open
        } else if(imageSet==true && meshSet==true && meshMeasured==false) {
            return 3; // both open, ? calibration
        } else if(imageSet==true && meshSet==true && meshMeasured==true) {
            return 4; // both open, ? calibrated, measurement made
        } else {
            return -1;
        }

    }

    double CorticalBone::getFixedCB() {
        if(isnan(globalFixedCB)) {
            cerr<<"WARNING: in CorticalBone::getFixedCB() CB value is not fixed"<<endl;
            return nan("1");
        } else {
            return globalFixedCB;
        }
    }

    bool CorticalBone::getFixedSigma(double &x, double &y, double &z) {
        if(modelProfile->IsSigmaSet()) {
            double sigma[Dimension];
            modelProfile->GetGlobalSigma(sigma);
            x = sigma[0]; y = sigma[1]; z = sigma[2];
            return true;
        } else {
            return false;
        }
    }

    double CorticalBone::getPeriostealOffset() {
        return modelProfile->GetPeriostealOffset();
    }

    bool CorticalBone::calculateThresholds(int thresholdIndex) {
        if(!classifierMethod->AreThresholdsSet() || classifierMethod->GetClassifierLevelIndex() != thresholdIndex) {

            if(thresholdIndex >= ClassifierTransform::kVLowPercent && thresholdIndex <= ClassifierTransform::kHighPercent) {
                if(verbose){cout<<"<b>----- Setting Global Thresholds -----</b>"<<endl;}
                Utilities::measureTime(true);
                estimateMinBasedThresholdsOverMesh(thresholdIndex);
                if(verbose){cout<<Utilities::getTabString()<<Utilities::measureTime(false)<<" to estimate the classifier thresholds"<<std::endl;}
                if(verbose){cout<<"<b>-----   Global Thresholds Set   -----</b>"<<endl;}
            } else if (thresholdIndex == ClassifierTransform::kMedianMidpoint) {
                if(verbose){cout<<"<b>----- Setting Global Thresholds -----</b>"<<endl;}
                Utilities::measureTime(true);
                estimateMedianBasedThresholdsOverMesh(thresholdIndex);
                if(verbose){cout<<Utilities::getTabString()<<Utilities::measureTime(false)<<" to estimate the classifier thresholds"<<std::endl;}
                if(verbose){cout<<"<b>-----   Global Thresholds Set   -----</b>"<<endl;}
            } else if (thresholdIndex == ClassifierTransform::kMedianManual) {
                cerr<<"Error in CorticalBone::calculateThresholds - Median manual thresholding selected without a weight " <<endl;
            } else {
                cerr<<"Error in CorticalBone::calculateThresholds - invalid classifier threshold selection type: "<< thresholdIndex <<endl;
            }
            thresholdWeight = nan("1");
        }
        return classifierMethod->AreThresholdsSet();
    }

    bool CorticalBone::calculateThresholds(int thresholdIndex, double weight) {
        if(!classifierMethod->AreThresholdsSet() || classifierMethod->GetClassifierLevelIndex() != thresholdIndex || thresholdWeight != weight) {

            if(thresholdIndex == ClassifierTransform::kMedianManual) {
                if(verbose){cout<<"<b>----- Setting Global Thresholds -----</b>"<<endl;}
                Utilities::measureTime(true);
                thresholdWeight=weight;
                estimateMedianBasedThresholdsOverMesh(thresholdIndex);
                if(verbose){cout<<Utilities::getTabString()<<Utilities::measureTime(false)<<" to estimate the classifier thresholds"<<std::endl;}
                if(verbose){cout<<"<b>-----   Global Thresholds Set   -----</b>"<<endl;}
            } else if (thresholdIndex >= ClassifierTransform::kVLowPercent && thresholdIndex <= ClassifierTransform::kMedianMidpoint) {
                thresholdWeight = nan("1");
                cerr<<"Error in CorticalBone::calculateThresholds - weight provided when not needed: " << thresholdIndex <<endl;
            } else {
                thresholdWeight = nan("1");
                cerr<<"Error in CorticalBone::calculateThresholds - invalid classifier threshold selection type: "<< thresholdIndex <<endl;
            }
        }
        return classifierMethod->AreThresholdsSet();
    }

    void CorticalBone::setClassifierThresholdInfo(double softTissue, double corticalBone, double threshold,
                                                  int thresholdIndex) {
        if(meshSet && imageSet && thresholdIndex != ClassifierTransform::kMedianManual) {
            if(verbose){cout<<"<b>----- Manually Setting Global Thresholds -----</b>"<<endl;}
            thresholdWeight=nan("1");
            classifierMethod->SetClassifierLevelIndex(thresholdIndex);
            classifierMethod->SetClassifierLevel(softTissue, corticalBone, threshold);
            if(verbose){cout<<"<b>-----   Global Thresholds Set   -----</b>"<<endl;}
            clearPreviousResults();
        } else if(thresholdIndex != ClassifierTransform::kMedianManual) {
            cerr<<"Error in CorticalBone::setClassifierThresholdInfo - weight not included for MedianManual mode" <<endl;
        } else {
            cerr<<"Error in CorticalBone::setClassifierThresholdInfo - mesh or image not set" <<endl;
        }
    }

    void CorticalBone::setClassifierThresholdInfo(double softTissue, double corticalBone, double threshold, double weight,
                                                  int thresholdIndex) {
        if(meshSet && imageSet && thresholdIndex == ClassifierTransform::kMedianManual) {
            if(verbose){cout<<"<b>----- Manually Setting Global Thresholds -----</b>"<<endl;}
            thresholdWeight = weight;
            classifierMethod->SetClassifierLevelIndex(thresholdIndex);
            classifierMethod->SetClassifierLevel(softTissue, corticalBone, threshold);
            if(verbose){cout<<"<b>-----   Global Thresholds Set   -----</b>"<<endl;}
            clearPreviousResults();
        } else if(thresholdIndex != ClassifierTransform::kMedianManual) {
            cerr<<"Error in CorticalBone::setClassifierThresholdInfo - not MedianManual mode" <<endl;
        } else {
            cerr<<"Error in CorticalBone::setClassifierThresholdInfo - mesh or image not set" <<endl;
        }
    }

    bool CorticalBone::getClassifierInfo(double &softTissue, double &corticalBone, double &threshold, int &thresholdindex, std::string &name) {
        if(meshSet && imageSet && classifierMethod->AreThresholdsSet()) {
            thresholdindex = classifierMethod->GetClassifierLevelIndex();
            classifierMethod->GetClassifierThresolds(softTissue, corticalBone, threshold);
            name = classifierMethod->GetClassifierLevelName();
            return true;
        } else {
            return false;
        }
    }

    bool CorticalBone::getClassifierInfo(double &softTissue, double &corticalBone, double &threshold, double &weight, int &thresholdindex, std::string &name) {
        if(meshSet && imageSet && classifierMethod->AreThresholdsSet() && classifierMethod->GetClassifierLevelIndex()==ClassifierTransform::kMedianManual) {
            thresholdindex = classifierMethod->GetClassifierLevelIndex();
            classifierMethod->GetClassifierThresolds(softTissue, corticalBone, threshold);
            name = classifierMethod->GetClassifierLevelName();
            weight = thresholdWeight;
            return true;
        } else {
            return false;
        }
    }

    bool CorticalBone::getClassifierIndex(double &index) { // only true if already processed
        if(meshSet && imageSet && classifierMethod->AreThresholdsSet()) {
            index = classifierMethod->GetClassifierLevelIndex();
            return true;
        } else {
            return false;
        }
    }

    bool CorticalBone::getclassifierName(std::string &name) {
        if(meshSet && imageSet && classifierMethod->AreThresholdsSet()) {
            name = classifierMethod->GetClassifierLevelName();
            return true;
        } else {
            return false;
        }
    }

    void CorticalBone::removeThresholds() {
        if(meshSet && imageSet) {
            classifierMethod->ClearThresholds();
            clearPreviousResults();
        }
    }

    /* state */
    bool CorticalBone::isCalibrated() {
        return calibrationSet;
    }

    bool CorticalBone::isCBFixed() { // note - is CBFixed globally and smoothing / FWHM not on

        return !isnan(globalFixedCB); // so can have smoothing/unconstatined with constraint, smoothing/unconstatining just adjusts/removes the constraint in later fitting runs
    }

    bool CorticalBone::isSigmaFixed() {
        return modelProfile->IsSigmaSet();
    }

    bool CorticalBone::isSetToFWHMMode() { // TODO consider combining with 'isCBFixed' into a single 'getCBConstraintMode()' method)
        return registrationMethod->GetFWHMMode();
    }

    bool CorticalBone::isMeshSet() {
        return meshSet;
    }

    bool CorticalBone::isImageSet() {
        return imageSet;
    }

    bool CorticalBone::isMeshMeasured() {
        return meshMeasured;
    }

    bool CorticalBone::isParametersImported() {
        return parametersImported;
    }

    bool CorticalBone::isClassifierMode() {
        return (modelIndex==kHighResClassifier);
    }

    bool CorticalBone::isProfileAveragingOn() {
        return modelProfile->IsProfileAveragingOn();
    }

    bool CorticalBone::isSampleNumberFixed() {
        return sampleNumberFixed;
    }

    bool CorticalBone::isSmoothingOn() {
        return smoothingSet;
    }

    bool CorticalBone::areThresholdsSet() {
        if(meshSet && imageSet) {
            return classifierMethod->AreThresholdsSet();
        } else {
            return false;
        }
    }

    bool CorticalBone::areResultsLoaded() { // are their locally generated results as apposed to loaded or no results
        return loadedResults;
    }

    std::string CorticalBone::getDisplayName() {

        if(modelIndex <= kEndostealRamp) {
            return registrationMethod->GetDisplayName(displayIndex);
        } else if(modelIndex == kHighResClassifier) {
            return classifierMethod->GetDisplayName(displayIndex);
        } else if(modelIndex == kCalibration) {
            return calibrationMethod->GetDisplayName(displayIndex);
        } else {
            cerr<<"Error in CorticalBone::getDisplayName(); invalid model index="<<modelIndex<<endl;
            return std::string("");
        }
    }

    /* processing methods */
    vtkSmartPointer<vtkDoubleArray> CorticalBone::getControlPointValues() {
        if(cntrlNumber > 0 && calibrationSet) {
            return cntrlVals;
        } else {
            vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
            return vals;
        }
    }

    bool CorticalBone::runCalibration(vtkSmartPointer<vtkDoubleArray> pts) {
        cntrlPts = pts;

        calP0=0.0; calP1=1.0; calP2=0.0; calibrationSet=false;

        if(phantomIndex==kNoCal) {
            noPhantomCalibration();
        } else if(phantomIndex==kMindwaySolidCal) {
            mindwaysCalibration();
        } else if (phantomIndex== kBoneDensityCal) {
            BoneDensityCalibration();
        } else if(phantomIndex== kManualQuadraticCal || phantomIndex==kManualLinearCal) {
            // calibration already preformed. check not nan
            if(isnan(calP0) || isnan(calP1) || isnan(calP2)) {
                calP0=0.0; calP1=1.0; calP2=0.0; calibrationSet=false;
            } else {
                phantomChanged = false; calibrationSet = true;
                // calibrationParametersSet=true; - leave as is
            }

        } else if(phantomIndex==kEuropeanSpineCal){
            ESPCalibration();
        } else if (phantomIndex == kManualControlPtsCal) {
            ControlPointCalibration();
        } else {
            cerr<<"Error in CorticalBone::runCalibration() invalid phantomIndex."<<endl;
        }

        // now update the reg method and displays
        modelProfile->SetCalibration(calP0,calP1,calP2); // set regardless of outcome as if it fails the cal vals are set to linear

        if(meshMeasured) {
            clearPreviousResults();
        }
        return calibrationSet;
    }

    bool CorticalBone::runModellingAtPoint(vtkIdType pointId) {

        if(!imageSet || !meshSet) {
            // throw exception
            throw std::logic_error( "mesh or image not yet set." );
        } else if(phantomChanged) {
            meshMeasured=false;
            cerr<<"Warning ignore CorticalBone::runModellingAtPoint() as calibration need to be completed"<<endl;
            return meshMeasured;
        }

        if(verbose){cout<<"<b>---- Run CBM At Point: "<<pointId<<" ----</b>"<<endl;}

        updateModelProfile(pointId, false);

        if( meshMask[pointId]==kMasked ) {
            if(verbose){cout<<Utilities::getTabString()<<"Masked Point: ignore"<<endl;}
        } else if(modelIndex <= kEndostealRamp) { // i.e. not thresholding
            if(verbose){cout<<Utilities::getTabString()<<"Model: "<<getModelName()<<","<<Utilities::getTabString()<<"Optimiser: "<<getOptimisationName()<<endl;}

            if(modelFitStatusArray[pointId] == kUndefined) { // run fitting as not already run
                runModelFittingPipeline(pointId, fittingSchemeIndex, modelIndex);
                updatePeriostealMeshNode(pointId);
            } else { // look up previous results
                ParametersType scales = (modelIndex==kThreeTierRect) ? rectFlexScales : rampFlexScales;
                registrationMethod->SetModelSelection(modelIndex); registrationMethod->SetTransformParameterScales(scales);
                registrationMethod->SetCombinedParameters(getTableRow( modelParametersTable, pointId));
                registrationMethod->Initialize();
            }

            if(modelFitStatusArray[pointId] == kOutOfBounds) {
                if(verbose){std::cout<<Utilities::getTabString()<<"Model Fitting Failed: Specified Point Outside Image Bounds."<<std::endl;}
            } else if(modelFitStatusArray[pointId]==kInvalid) {
                if(verbose){std::cout<<Utilities::getTabString()<<"Model Fitting Failed: Invalid Model Parameters. "<<std::endl;}
            }

        } else if(modelIndex == kHighResClassifier) { // i.e. thresholding

            if(verbose){cout<<Utilities::getTabString()<<"Model: "<<getModelName()<<endl;}

            if(modelFitStatusArray[pointId] == kUndefined) {
                runThresholdingPipeline(pointId);
                updatePeriostealMeshNode(pointId);
            } else {
                classifierMethod->Process(); // re-run but do not re-save the results
            }

            if(modelFitStatusArray[pointId] == kOutOfBounds) {
                if(verbose){std::cout<<"Thresholding Failed: Specified Point Outside Image Bounds."<<std::endl;}
            } else if(modelFitStatusArray[pointId]==kInvalid) {
                if(verbose){std::cout<<"Thresholding Failed: Error. "<<std::endl;}
            }
        } else if(modelIndex == kCalibration) {


            if(verbose){cout<<Utilities::getTabString()<<"Calibration: "<<getModelName()<<endl;}

            if(modelFitStatusArray[pointId] == kUndefined) {
                runGenerateCalibrationValues(pointId);
            } else {
                runGenerateCalibrationValues(pointId); // re-run but do not re-save the results
            }

            if(modelFitStatusArray[pointId] == kOutOfBounds) {
                if(verbose){std::cout<<"Calibration Failed: Specified Point Outside Image Bounds."<<std::endl;}
            } else if(modelFitStatusArray[pointId]==kInvalid) {
                if(verbose){std::cout<<"Calibration Failed: Maximum profile values is in the first 20% of the profile. "<<std::endl;}
            }

        } else {
            cerr<<"Error in CorticalBone::runModellingAtPoint Invalid model index of: "<<modelIndex<<endl;
            return false;
        }

        writeModelResults(pointId);

        if(verbose){cout<<"<b>---- Finish Running CBM ----</b>"<<endl;}

        return modelFitStatusArray[pointId]==kValid;
    }

    vtkIdType CorticalBone::getClosestPoint(double point[Dimension]) {
        vtkIdType closestPtId = meshTree->FindClosestPoint(point);

        return closestPtId;
    }

    bool CorticalBone::runModellingOverMesh() {

        if(phantomChanged) {
            meshMeasured=false;
            cerr<<"Warning ignore CorticalBone::runModellingOverMesh() as calibration need to be completed"<<endl;
            return meshMeasured;
        } else if(imageSet && meshSet && !meshMeasured) {

            if(verbose){cout<<"Model: "<<getModelName();
                if(modelIndex<=kEndostealRamp){cout<<Utilities::getTabString()<<"Optimiser: "<<getOptimisationName();}
                cout<<endl;}

            Utilities::measureTime(true);

            if(modelIndex == kHighResClassifier) { // thresholding
                runClassifierOverMesh();
            } else if(modelIndex == kCalibration) {
                runCalibrationOverMesh();
            } else if(fittingSchemeIndex < kCBSmoothingAFitting) { // model fitting
                resetCBConstraint();
                runCBMOverMeshWithNoSmoothing();
            } else { // model fitting
                resetCBConstraint();
                runCBMOverMeshWithCBSmoothing();
            }
            
            //calculatePeriostealMesh();

            if(verbose){
                cout<<Utilities::getTabString()<<Utilities::measureTime(false)<<" elapsed time; "<<failureCount<<" failures (with "<<nanCount
                        << " of these Nan), and "<<outOfBoundsCount<<" points out of Bounds of a total of "<<mesh->GetNumberOfPoints()-maskedCount;
                if(maskedCount==0) {cout<<" points"<<std::endl;} else {cout<<" unmasked points."<<std::endl;}
            }


        } else if(meshMeasured && !loadedResults) {
            if(verbose){cout<<"Model: "<<getModelName();
                if(modelIndex<=kEndostealRamp){cout<<Utilities::getTabString()<<"Optimiser: "<<getOptimisationName();}
                cout<<endl;
                std::cout<<"Previously Run; "<<failureCount<<" failures (with "<<nanCount
                        << " of these Nan), and "<<outOfBoundsCount<<" points out of Bounds of a total of "<<mesh->GetNumberOfPoints()-maskedCount;
                if(maskedCount==0) {cout<<" points"<<std::endl;} else {cout<<" unmasked points."<<std::endl;}
            }
        } else if(loadedResults) {
            if(verbose){cout<<"Model: "<<getModelName();
                if(modelIndex<=kEndostealRamp){cout<<Utilities::getTabString()<<"Optimiser: "<<getOptimisationName();}
                cout<<endl;
                std::cout<<"Previously Run loaded results; "<<nanCount<<" nan values (either out of bounds or failures) of a total of "<<mesh->GetNumberOfPoints()-maskedCount;
                if(maskedCount==0) {cout<<" points"<<std::endl;} else {cout<<" unmasked points."<<std::endl;}
            }
        }
        return meshMeasured;
    }

    bool CorticalBone::runCBMOverMeshWithCBSmoothing() {

        if(!smoothingSet) {
            meshMeasured=false; cerr<<"Error in CorticalBone::runCBMOverMeshWithCBSmoothing() smoothing values not set"<<endl;
            return meshMeasured;
        }

        vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
        vtkIdType numPoints = points->GetNumberOfPoints();
        double unconstrainedCBEstimates[numPoints], unconstrainedABSErrors[numPoints], smoothedCBEstimates[numPoints];
        double unconstrainedCBPrecisions[numPoints];

        // select model parameters
        int localSchemeIndex1, localSchemeIndex2, localModelIndex;

        if(fittingSchemeIndex == kCBSmoothingAFitting) { // general - combined cortical thickness, cortical density and endocortical region  - SELECTED METHOD
            localSchemeIndex1 = kUnconstrainedAFitting; localSchemeIndex2 = kStdBFitting;
            localModelIndex = modelIndex;
        } else if(fittingSchemeIndex == kCBSmoothingBFitting) { // cortical - bias for accurate cortical measurements
            localSchemeIndex1 = kUnconstrainedAFitting; localSchemeIndex2 = kStdCFitting;
            localModelIndex = modelIndex;
        } else if(fittingSchemeIndex == kCBSmoothingCFitting) { // endo & outer cortical - bias for accurate outer and endo measurements
            localSchemeIndex1 = kUnconstrainedBFitting; localSchemeIndex2 = kStdBFitting;
            localModelIndex = modelIndex;
        } else { // old oposit model to give more accurate desities
            localSchemeIndex1 = kUnconstrainedCFitting; localSchemeIndex2 = kStdCFitting;
            localModelIndex = (modelIndex==kThreeTierRect) ? kEndostealRamp : kThreeTierRect;
        }

        // run constrained to unconstrained fitting scheme - store the optimised unconstrained CB density estimate
        for(vtkIdType pointId = 0; pointId < numPoints; ++pointId) {

            if(meshMask[pointId]==kMasked) {
                unconstrainedCBEstimates[pointId] = nan("1");
                unconstrainedABSErrors[pointId] = nan("1");
            } else {
                // perform registration
                updateModelProfile(pointId);


                runModelFittingPipeline(pointId, localSchemeIndex1, localModelIndex);
                if(fittingSchemeIndex == kCBSmoothingAFitting) {
                    unconstrainedCBEstimates[pointId] = registrationMethod->GetCorticalDensity();
                } else {
                    unconstrainedCBEstimates[pointId] = registrationMethod->GetMeanCorticalDensity(); // average for endo, CB for rect
                }
                unconstrainedABSErrors[pointId] = registrationMethod->GetErrorMean();
                //cerr<<"pointId="<<pointId<<" Covariance Matrix Sizes";
                unconstrainedCBPrecisions[pointId] = registrationMethod->GetCorticalBonePrecision();
                //std::clog<<"pointID="<<pointId<<", precision="<<unconstrainedCBPrecisions[pointId]<<", ABS Error="<<unconstrainedABSErrors[pointId]<<endl;
            }
        }


        // run over mesh smoothing the values
        double radius = smoothingRadius*3; // note inter-vertex distances are generally in the range of .9 to 4mm
        double normalScale = 1/(sqrt(2*M_PI*smoothingRadius*smoothingRadius)), normalExpScale = -1/(2*smoothingRadius*smoothingRadius);
        for(vtkIdType pointId = 0; pointId < numPoints; ++pointId) {

            // location of point to smooth
            double location[Dimension];
            meshTree->GetDataSet()->GetPoint(pointId, location);

            // calculate local smoothed density
            vtkSmartPointer<vtkIdList> nearPoints = vtkSmartPointer<vtkIdList>::New();
            meshTree->FindPointsWithinRadius(radius, location, nearPoints);
            vtkIdType numberNearPoints = nearPoints->GetNumberOfIds();

            // distance and error weights
            double smoothedCBEstimate = 0, weightsSum = 0;

            for(vtkIdType i=0; i<numberNearPoints; i++) {

                vtkIdType nearPtId = nearPoints->GetId(i);

                if(!isnan(unconstrainedCBEstimates[nearPtId]) && !isnan(unconstrainedABSErrors[nearPtId]) && !isnan(unconstrainedCBPrecisions[nearPtId]) && unconstrainedCBPrecisions[nearPtId]!=0) { // only include in estimate if valid measurement made there
                    // todo - what if precision=0. Note in bern scan 2L a pt 7216, precision is 0, so inf weighting
                    double nearLocation[Dimension];
                    meshTree->GetDataSet()->GetPoint(nearPtId, nearLocation);

                    double distance = sqrt(vtkMath::Distance2BetweenPoints(location, nearLocation));

                    double distanceWeight = normalScale * exp(distance*distance * normalExpScale);
                    double errorWeight = 1.0/unconstrainedABSErrors[nearPtId];
                    double precisionWeight = 1.0/unconstrainedCBPrecisions[nearPtId]; // related to covariance
                    //cout<<distance<<"("<<unconstrainedCBEstimates[nearPtId]<<"):   1/"<<unconstrainedABSErrors[nearPtId]<<"="<<errorWeight<<" x "<<distanceWeight<<"="<<errorWeight*distanceWeight<<";\t\t";
                    smoothedCBEstimate += unconstrainedCBEstimates[nearPtId] * distanceWeight * errorWeight * precisionWeight;
                    weightsSum += distanceWeight * errorWeight * precisionWeight;
                    //cout<<smoothedCBEstimate<<", "<<weightsSum<<endl;
                }
            }
            //cout<<smoothedCBEstimate;
            smoothedCBEstimate = (weightsSum != 0) ? (smoothedCBEstimate / weightsSum) : globalFixedCB;
            smoothedCBEstimates[pointId] = smoothedCBEstimate;
        }

        // if version D perform a second distance weighting
        if(true || fittingSchemeIndex == kCBSmoothingDFitting) {
            // increase the smoothing radius by 3x
            normalScale = 1/(sqrt(2*M_PI*radius*radius)), normalExpScale = -1/(2*radius*radius); radius *= 3;
            double doubleSmoothedCBEstimates[numPoints];
            for (vtkIdType pointId = 0; pointId < numPoints; ++pointId) {

                // location of point to smooth
                double location[Dimension]; meshTree->GetDataSet()->GetPoint(pointId, location);

                // calculate local smoothed density
                vtkSmartPointer<vtkIdList> nearPoints = vtkSmartPointer<vtkIdList>::New();
                meshTree->FindPointsWithinRadius(radius, location, nearPoints);
                vtkIdType numberNearPoints = nearPoints->GetNumberOfIds();

                // distance weighted smoothing
                double smoothedCBEstimate = 0, weightsSum = 0;
                for (vtkIdType i = 0; i < numberNearPoints; i++) {

                    vtkIdType nearPtId = nearPoints->GetId(i);

                    if (!isnan(smoothedCBEstimates[nearPtId])) { // only include if valid estimate

                        double nearLocation[Dimension]; meshTree->GetDataSet()->GetPoint(nearPtId, nearLocation);
                        double distance = sqrt(vtkMath::Distance2BetweenPoints(location, nearLocation));

                        double distanceWeight = normalScale * exp(distance * distance * normalExpScale);
                        smoothedCBEstimate += smoothedCBEstimates[nearPtId] * distanceWeight;
                        weightsSum += distanceWeight;
                        //cout<<smoothedCBEstimate<<", "<<weightsSum<<endl;
                    }
                }
                //cout<<smoothedCBEstimate;
                doubleSmoothedCBEstimates[pointId] = (weightsSum != 0) ? (smoothedCBEstimate / weightsSum) : globalFixedCB;
            }

            for (vtkIdType pointId = 0; pointId < numPoints; ++pointId) {
                smoothedCBEstimates[pointId] = doubleSmoothedCBEstimates[pointId];
            }
        }

        // reset the error counts
        modelFitStatusArray.SetSize((unsigned int)mesh->GetNumberOfPoints());
        modelFitStatusArray.Fill(kUndefined);
        nanCount = 0; failureCount = 0; outOfBoundsCount = 0;
        maskResultArrays();

        // rerun but with new CB estimates
        for(vtkIdType pointId = 0; pointId < numPoints; ++pointId) {

            if(modelFitStatusArray[pointId]!=kMasked) {
                // update the CB estimate
                if(!isnan(smoothedCBEstimates[pointId])) {
                    setFixedCBDensity(smoothedCBEstimates[pointId], false);
                } else {
                    removeFixedCBDensity(false); // todo decide wither or not to use global, unconstrained or fit-failure
                }


                // perform registration
                updateModelProfile(pointId);
                // todo what fitting scheme to use (i.e. what weightings, & what scales?)
                runModelFittingPipeline(pointId, localSchemeIndex2, modelIndex);
                updatePeriostealMeshNode(pointId);
            }
        }

        meshMeasured = true;
        return meshMeasured;

    }

    bool CorticalBone::runCBMOverMeshWithNoSmoothing() {

        vtkSmartPointer<vtkPoints> points = mesh->GetPoints();

        vtkIdType numPoints = points->GetNumberOfPoints();

        // cycle through mesh
        for(vtkIdType pointId = 0; pointId < numPoints; ++pointId) {

            // perform registration
            if(modelFitStatusArray[pointId] == kUndefined) {
                updateModelProfile(pointId);
                runModelFittingPipeline(pointId, fittingSchemeIndex, modelIndex);
                updatePeriostealMeshNode(pointId);
            }

        }

        meshMeasured = true;
        return meshMeasured;
    }

    bool CorticalBone::runClassifierOverMesh() {

        if(!classifierMethod->AreThresholdsSet()) {
            clearPreviousResults();
            cerr<<" Error in CorticalBone::runThresholdingOverMesh() thresholds not set"<<endl;
            return false;
        }

        vtkIdType numPoints = mesh->GetPoints()->GetNumberOfPoints();

        // cycle through mesh
        for(vtkIdType pointId = 0; pointId < numPoints; ++pointId) {

            //cerr<<"point ID="<<pointId<<", ";
            // perform registration
            updateModelProfile(pointId);
            if(modelFitStatusArray[pointId] == kUndefined) {
                runThresholdingPipeline(pointId);
                updatePeriostealMeshNode(pointId);
            }
        }

        meshMeasured = true;
        return meshMeasured;
    }

    bool CorticalBone::runCalibrationOverMesh() {

        vtkIdType numPoints = mesh->GetPoints()->GetNumberOfPoints();

        // cycle through mesh
        for(vtkIdType pointId = 0; pointId < numPoints; ++pointId) {

            updateModelProfile(pointId);
            if(modelFitStatusArray[pointId] == kUndefined) {
                runGenerateCalibrationValues(pointId);
            }
        }

        meshMeasured = true;
        return meshMeasured;
    }

    bool CorticalBone::estimateMinBasedThresholdsOverMesh(int percentIndex) {

        classifierMethod->SetClassifierLevelIndex(percentIndex);

        clearPreviousResults();

        vtkSmartPointer<vtkPoints> points = mesh->GetPoints();

        vtkIdType numPoints = points->GetNumberOfPoints();

        double stEstimateSum=0, cbEstimateSum=0, thresholdSum=0; int validCount = 0;

        modelProfile->TurnOnAllSamplesRequired(); // try speed up by reducing the number of calculations / stored values

        // cycle through mesh - take average - change to take median as more noise resistant
        for(vtkIdType pointId = 0; pointId < numPoints; ++pointId) {

            if(meshMask[pointId]!=kMasked && meshMask[pointId]!=kDblPeak) {

                //cerr<<"Point ID="<<pointId;
                // perform registration
                updateModelProfile(pointId);
                if (classifierMethod->InitialiseDensities()) { // an 'IsInsideImage' test built into the InitialiseDensities method
                    double st, cb, threshold, percent;
                    classifierMethod->GetClassifierThresolds(st, cb, threshold);
                    stEstimateSum += st;
                    cbEstimateSum += cb;
                    thresholdSum += threshold;

                    validCount++;
                }
            }
        }

        bool status = classifierMethod->SetClassifierLevel(stEstimateSum / validCount, cbEstimateSum / validCount, thresholdSum / validCount);
        meshMeasured = false;

        modelProfile->TurnOffAllSamplesRequired(); // todo assumes the modelIndex=kHighResClassifier. check is true

        //cerr<<"Classifier Thresholds for percent Index="<<highResolutionMethod->GetClassifierLevelName()<<", Calculated from "<<validCount<<" of "<<numPoints<<" measurement points."<<endl;
        //cerr<<"Thresholds: Soft Tissue="<<stEstimateSum / validCount<<", Cortical Bone Density="<<cbEstimateSum / validCount<<", Thresholds="<<thresholdSum / validCount<<endl;

        if(verbose){cout<<"Classifier Thresholds for percent Index="<< classifierMethod->GetClassifierLevelName()<<", Calculated from "<<validCount<<" of "<<numPoints<<" measurement points."<<endl;}
        if(verbose){cout<<"Thresholds: Soft Tissue="<<stEstimateSum / validCount<<", Cortical Bone Density="<<cbEstimateSum / validCount<<", Thresholds="<<thresholdSum / validCount<<endl;}

        if(!status) { cerr<<"Error in CorticalBone::estimateMinBasedThresholdsOverMesh()"<<endl; }

        return meshMeasured;
    }

    bool CorticalBone::estimateMedianBasedThresholdsOverMesh(int thresholdIndex) {

        classifierMethod->SetClassifierLevelIndex(thresholdIndex);

        clearPreviousResults();

        vtkSmartPointer<vtkPoints> points = mesh->GetPoints();

        vtkIdType numPoints = points->GetNumberOfPoints();

        // todo - create two arrays to house the median maxs and mins at each point
        vtkIdType upperNumberLimit = modelProfile->GetMaximumPossibleProfiles() * (numPoints - maskedCount);
        itk::Array<double> allMaximums = itk::Array<double>(upperNumberLimit); itk::Array<double> allMinimums = itk::Array<double>(upperNumberLimit);

        int count = 0;

        modelProfile->TurnOnAllSamplesRequired(); // temporarily turn on to try speed up by reducing the number of calculations / stored values
        modelProfile->TurnOnMinDetection();

        // cycle through mesh - take average - change to take median as more noise resistant
        for(vtkIdType pointId = 0; pointId < numPoints; ++pointId) {

            if(meshMask[pointId]!=kMasked) { // an 'IsInsideImage' test built into the InitialiseDensities method
                // treat kDblPeak the same as another peak won't significantly effect the max/min values
                updateModelProfile(pointId);
                if(modelProfile->IsPtInsideImage()) {
                    itk::Array<double> maximums = modelProfile->GetMaxValues();
                    itk::Array<double> minimums = modelProfile->GetMinValues();

                    int n = maximums.size();
                    for (int i = 0; i < n; i++) {
                        allMaximums[count] = maximums[i];
                        allMinimums[count] = minimums[i];
                        count++;
                    }
                }
            }
        }


        double maximum = Utilities::CalculateArrayMedian(allMaximums, 0, count-1);
        double minimum = Utilities::CalculateArrayMedian(allMinimums, 0, count-1);
        double threshold; bool status;
        if(thresholdIndex==ClassifierTransform::kMedianMidpoint) {
            threshold = (maximum + minimum) / 2;
            status = classifierMethod->SetClassifierLevel(minimum, maximum, threshold);
        } else if(thresholdIndex==ClassifierTransform::kMedianManual) {
            threshold = thresholdWeight * maximum + (1.0-thresholdWeight)*minimum;
            status = classifierMethod->SetClassifierLevel(minimum, maximum, threshold);
        } else {
            cerr<<"Invalid thresholdIndex in estimateMedianBasedThresholdsOverMesh of "<<thresholdIndex<<endl;
            clearPreviousResults();
        }

        meshMeasured = false;

        modelProfile->TurnOffAllSamplesRequired();
        modelProfile->TurnOffMinDetection();

        //cerr<<"Classifier Thresholds for percent Index="<<highResolutionMethod->GetClassifierLevelName()<<", Calculated from "<<validCount<<" of "<<numPoints<<" measurement points."<<endl;
        //cerr<<"Thresholds: Soft Tissue="<<stEstimateSum / validCount<<", Cortical Bone Density="<<cbEstimateSum / validCount<<", Thresholds="<<thresholdSum / validCount<<endl;

        if(verbose){cout<<"Classifier Median based threshold  calculated from "<<count<<" profiles from "<<numPoints-maskedCount<<" unmasked measurement points."<<endl;}
        if(verbose){cout<<"Thresholds: Not bone="<<minimum<<", Cortical Bone Density="<<maximum<<", Threshold="<<threshold<<endl;}

        if(!status) { cerr<<"Error in CorticalBone::estimateMinBasedThresholdsOverMesh()"<<endl; }

        return meshMeasured;
    }

    /* Display API */
    void CorticalBone::getMeasurementLocation(double location[Dimension]) {

        modelProfile->GetEdge(location);
    }

    vtkSmartPointer<vtkPoints> CorticalBone::getProfilePoints() {
        return modelProfile->GetProfilePoints();
    }

    void CorticalBone::getProfileStart(double start[Dimension]) {
        modelProfile->GetStart(start);
    }

    void CorticalBone::getProfileEnd(double end[Dimension]) {
        modelProfile->GetEnd(end);
    }

    void CorticalBone::getDisplayRange(double &min, double &max) {
        itk::Array2D<CoordinateType> displayRanges;

        if(modelIndex<=kEndostealRamp) {
            displayRanges = registrationMethod->GetDisplayRanges(parametersImported);
        } else if(modelIndex==kHighResClassifier) {
            displayRanges = classifierMethod->GetDisplayRanges(parametersImported);
        } else if(modelIndex==kCalibration) {
            displayRanges = calibrationMethod->GetDisplayRanges(parametersImported);
        } else {
            cerr<<"Error in CorticalBone::getDisplayRange(); Invalid model index="<<modelIndex;
            min=nan("1"); max=nan("1");
            return;
        }

        min = displayRanges.get(displayIndex,0); max = displayRanges.get(displayIndex,1);
    }

    vtkSmartPointer<vtkDoubleArray> CorticalBone::getDisplayArray() {
        return vtkDoubleArray::SafeDownCast( modelDisplayTable->GetColumn(displayIndex) );
    }

    //----------------------------------------------------------------------------//
    //------------------ Private Methods -----------------------------------------//
    //----------------------------------------------------------------------------//

    /* Display Updates */
    // parameters
    void CorticalBone::clearParametersTable() {

        int n = getNumberOfParameters();

        modelParametersTable = vtkSmartPointer<vtkTable>::New();
        for(int i=0; i<n; i++) {

            std::string parameterName = getParameterName(i);
            vtkSmartPointer<vtkDoubleArray> parameterArray = vtkSmartPointer<vtkDoubleArray>::New();
            parameterArray->SetName(parameterName.c_str());
            modelParametersTable->AddColumn(parameterArray);
        }
        modelParametersTable->SetNumberOfRows(mesh->GetNumberOfPoints());
    }

    void CorticalBone::clearMeshAttributeTable() {

        modelDisplayTable = vtkSmartPointer<vtkTable>::New();
        int n = getNumberOfDisplays();
        for(int i=0; i<n; i++) {

            std::string displayName = getDisplayName(i);

            vtkSmartPointer<vtkDoubleArray> displayArray = vtkSmartPointer<vtkDoubleArray>::New();
            displayArray->SetNumberOfComponents(1);
            displayArray->SetName(displayName.c_str());
            modelDisplayTable->AddColumn(displayArray);
        }

        modelDisplayTable->SetNumberOfRows(mesh->GetNumberOfPoints());
    }

    void CorticalBone::clearImageProfilesArray() {
        imageProfilesArray = vtkSmartPointer<vtkDoubleArray>::New();
        imagePositionsArray = vtkSmartPointer<vtkDoubleArray>::New();

        imageProfilesArray->SetNumberOfComponents(numberOfSamples);
        imageProfilesArray->SetNumberOfTuples(mesh->GetNumberOfPoints());

        imagePositionsArray->SetNumberOfComponents(numberOfSamples);
        imagePositionsArray->SetNumberOfTuples(mesh->GetNumberOfPoints());

        if(modelIndex==kHighResClassifier) {
            imageClassificationArray = vtkSmartPointer<vtkDoubleArray>::New();

            imageClassificationArray->SetNumberOfComponents(numberOfSamples);
            imageClassificationArray->SetNumberOfTuples(mesh->GetNumberOfPoints());

            imagePercentageArray = vtkSmartPointer<vtkDoubleArray>::New();

            imagePercentageArray->SetNumberOfComponents(numberOfSamples);
            imagePercentageArray->SetNumberOfTuples(mesh->GetNumberOfPoints());
        }

    }

    void CorticalBone::maskResultArrays() {

        vtkIdType numPoints = mesh->GetPoints()->GetNumberOfPoints();
        for(vtkIdType i=0; i<numPoints; i++) {
            if (meshMask[i] == kMasked) {
                modelFitStatusArray[i] = kMasked;

                for (int j = 0; j < modelParametersTable->GetNumberOfColumns(); j++) {
                    modelParametersTable->SetValue(i, j, nan("1"));
                }

                for (int j = 0; j < numberOfSamples; j++) {
                    imageProfilesArray->SetComponent(i, j, nan("1"));
                    imagePositionsArray->SetComponent(i, j, nan("1"));
                }

                for (int j = 0; j < modelDisplayTable->GetNumberOfColumns(); j++) {
                    modelDisplayTable->SetValue(i, j, nan("1"));
                }
            }
        }

        if(modelIndex==kHighResClassifier) {
            for(vtkIdType i=0; i<numPoints; i++) {
                if (meshMask[i] == kMasked) {
                    for (int j = 0; j < numberOfSamples; j++) {
                        imageClassificationArray->SetComponent(i, j, nan("1"));
                        imagePercentageArray->SetComponent(i, j, nan("1"));
                    }
                }
            }
        }
    }

    // display values
    void CorticalBone::storeDisplayValues(vtkIdType pointId) {

        itk::Array<double> displayValues;

        if(!modelProfile->IsPtInsideImage()) {
            displayValues.SetSize((unsigned int)modelDisplayTable->GetNumberOfColumns()); displayValues.Fill(nan("1"));
        } else if(modelIndex <= kEndostealRamp) {

            displayValues = registrationMethod->GetDisplayValues(); // true so nan if model is invalid

            if(parametersImported) { // get the error between the local model and the imported data - modelAtributeTable values set within method
                calculateImportErrors(pointId);
            }
        } else if(modelIndex == kHighResClassifier) { // no import values
            displayValues = classifierMethod->GetDisplayValues();
        } else if(modelIndex == kCalibration) { // no import values
            displayValues = calibrationMethod->GetDisplayValues();
        } else {
            cerr<<"Error in CorticalBone::storeParameterValues() invalid modelIndex="<<modelIndex<<endl;
            return;
        }

        int n = displayValues.size();
        for(int i=0; i<n; i++) {
            modelDisplayTable->SetValue(pointId, i, displayValues[i]);
        }

    }

    void CorticalBone::storeParameterValues(vtkIdType pointId) {

        if(modelIndex<kThreeTierRect || modelIndex>kCalibration) {
            cerr<<"Error in CorticalBone::storeParameterValues() invalid modelIndex="<<modelIndex<<endl;
            return;
        }

        itk::Array<double> parameters;
        if(!modelProfile->IsPtInsideImage()) {
            parameters.SetSize((unsigned int)modelParametersTable->GetNumberOfColumns()); parameters.Fill(nan("1"));
        } else if(modelIndex <=kEndostealRamp) {
            parameters = registrationMethod->GetCombinedParameters();
        } else if(modelIndex == kHighResClassifier) {
            parameters = classifierMethod->GetParameterValues();
        } else if(modelIndex == kCalibration) {
            parameters = calibrationMethod->GetParameterValues();
        }

        for(int i=0; i<parameters.GetSize(); i++) {
            modelParametersTable->SetValue(pointId, i, parameters[i]);
        }

        if(modelIndex <=kEndostealRamp && modelParametersTable->GetNumberOfColumns() > parameters.GetSize()) { // manually set endosteal width if missing
            modelParametersTable->SetValue(pointId, modelParametersTable->GetNumberOfColumns()-1, 0);
        }

    }

    void CorticalBone::storeSampledProfiles(vtkIdType pointId) {

        itk::Array<double> imageProfile = modelProfile->GetValues();
        itk::Array<double> positions = modelProfile->GetPositions();
        for(int i=0; i<numberOfSamples; i++) {
            imageProfilesArray->SetComponent(pointId, i, imageProfile[i]);
            imagePositionsArray->SetComponent(pointId, i, positions[i]);
        }


        if(modelIndex == kHighResClassifier) {
            itk::Array<double> classifications = classifierMethod->GetClassifications();
            itk::Array<double> percentages = classifierMethod->GetPercentages();

            for(int i=0; i<numberOfSamples; i++) {
                imageClassificationArray->SetComponent(pointId, i, classifications[i]);
                imagePercentageArray->SetComponent(pointId, i, percentages[i]);
            }
        }
    }
    
    bool CorticalBone::calculatePeriostealMesh() { // generates a mesh representing the periosteal surface
        bool status = true;
        
        if(meshMeasured && modelIndex<=kHighResClassifier) {
            periostealMesh = vtkSmartPointer<vtkPolyData>::New();
            periostealMesh->DeepCopy(mesh);
            vtkSmartPointer<vtkPoints> points = periostealMesh->GetPoints();
            
            int xPIndex = (modelIndex <= kEndostealRamp) ? (int)ModelRegistrationType::kXpParam : (int)ClassifierTransform::kCBStart;
            double initialMeshLocation = 0; //defaultPeriostealEdgeRatio * profileLength;
            
            for(vtkIdType pointId=0; pointId<periostealMesh->GetNumberOfPoints(); pointId++) { // surface normals are assumed to point out; + offset = move out
                
                // get offset from origional edge
                double offset = initialMeshLocation - modelParametersTable->GetValue(pointId, xPIndex).ToDouble(); // + offset means move out, - offset means move in
                if(!isnan(offset)) { // only update position if not an nan measure
                    // get normal at this location
                    double normal[3]; meshNormals->GetTuple(pointId, normal);
                    // get point at this location
                    double location[Dimension]; periostealMesh->GetPoint(pointId, location);
                    // calculate the new point location
                    location[0]=location[0]+normal[0]*offset; location[1]=location[1]+normal[1]*offset; location[1]=location[1]+normal[1]*offset;
                    points->SetPoint(pointId, location);
                }
            }
            periostealMesh->SetPoints(points);
            
        } else {
            periostealMesh=NULL;
            status=false;
        }
        return status;
    }

    void CorticalBone::resetPeriostealMesh() {
        if(meshSet) {
            periostealMesh = vtkSmartPointer<vtkPolyData>::New();
            periostealMesh->DeepCopy(mesh);
        }
    }

    bool CorticalBone::updatePeriostealMeshNode(vtkIdType pointId) {

        if(meshSet && modelFitStatusArray[pointId]==kValid) {

            vtkSmartPointer<vtkPoints> points = periostealMesh->GetPoints();

            int xPIndex = (modelIndex <= kEndostealRamp) ? (int) ModelRegistrationType::kXpParam
                                                         : (int) ClassifierTransform::kCBStart;
            double initialMeshLocation = 0;

            // surface normals are assumed to point out; + offset = move out

            // get offset from origional edge
            double offset = initialMeshLocation - modelParametersTable->GetValue(pointId,
                                                                                 xPIndex).ToDouble(); // + offset means move out, - offset means move in
            if (!isnan(offset)) { // only update position if not an nan measure
                // get normal at this location
                double normal[3];
                meshNormals->GetTuple(pointId, normal);
                // get point at this location
                double location[Dimension];
                periostealMesh->GetPoint(pointId, location);
                // calculate the new point location
                location[0] = location[0] + normal[0] * offset;
                location[1] = location[1] + normal[1] * offset;
                location[1] = location[1] + normal[1] * offset;
                points->SetPoint(pointId, location);
            }

            /*double periostealPoint[Dimension]; // todo - us to replace above and coorectly calculate position along a potentially curved profile
            modelProfile->GetPeriostealMeshNode(periostealPoint);
            points->SetPoint(pointId, periostealPoint);*/

            periostealMesh->SetPoints(points);
        }
        return true;
    }

    void CorticalBone::setModelToImport(vtkIdType pointId) { // todo - do I need to set / reset calibration
        // set model index
        registrationMethod->SetModelSelection(importedModelIndex);
        // remove any calibration
        modelProfile->SetCalibration(0.0,1.0,0.0); // all imported values are already expected to be calibrated
        // get parameters and scale as necessary
        ParametersType importedParameters = getTableRow(importedParametersTable, pointId);

        // set parameters and max CB density
        registrationMethod->SetMaxBMDDensity(importedParameters[ModelRegistrationType::kYcbParam]);
        registrationMethod->SetCombinedParameters(importedParameters);
    }

    void CorticalBone::calculateImportErrors(vtkIdType pointId) { // only relevant when the local model is fitted as it must be continuous along the profile

        if(modelIndex>kEndostealRamp) {
            cerr<<"CorticalBone::calculateImportErrors was entered in error for a local method that is not fitting based"<<endl;
            return;
        }

        int n = registrationMethod->GetNumberOfDisplayValues(parametersImported);

        if(!isImportedModelValid(pointId) || !registrationMethod->IsValid()) {
            modelDisplayTable->SetValue(pointId, n-3, nan("1"));
            modelDisplayTable->SetValue(pointId, n-2, nan("1"));
            modelDisplayTable->SetValue(pointId, n-1, nan("1"));
            return;
        }

        // calculate the remainder
        double increment = modelProfile->GetIncrement();
        double offset = calculateImportOffset(pointId);

        if(isnan(offset)) { // hr is out of bounds or LR is invalid
            modelDisplayTable->SetValue(pointId, ModelRegistrationType::kImportImageBias, nan("1"));
            modelDisplayTable->SetValue(pointId, ModelRegistrationType::kImportRMSError, nan("1"));
            modelDisplayTable->SetValue(pointId, ModelRegistrationType::kImportImageSTD, nan("1"));
            return;
        }

        // calculate errors - re-sample along local image to
        double imageAlignmentBiases[numberOfSamples];
        double imageAlignmentBias = 0, ModelToImageSqrError = 0;
        int numberOfValidValues=0;
        itk::Array<double> imageValues = modelProfile->GetOffsetProfileImageValues(offset);
        for(int i=0; i<numberOfSamples; i++) {

            double importValue = importedProfilesArray->GetComponent(pointId, i);

            if(!isnan(importValue)) {
                numberOfValidValues++;


                double modelValue = registrationMethod->GetModelValue(i*increment-offset);

                imageAlignmentBias += importValue - imageValues[i];
                ModelToImageSqrError += (imageValues[i] - modelValue) * (imageValues[i] - modelValue);
                imageAlignmentBiases[i] = importValue - imageValues[i];
            } else {
                imageAlignmentBiases[i]=nan("1");
            }
        }

        double meanImageAlignmentBias = imageAlignmentBias / numberOfValidValues;

        // calculate the image standard deviation
        double importVariance = 0;
        for (int i = 0; i < numberOfSamples; i++) {
            if (!isnan(imageAlignmentBiases[i])) {
                importVariance += (imageAlignmentBiases[i] - meanImageAlignmentBias) *
                                  (imageAlignmentBiases[i] - meanImageAlignmentBias) / numberOfValidValues;
            }
        }
        modelDisplayTable->SetValue(pointId, ModelRegistrationType::kImportImageBias, meanImageAlignmentBias);
        modelDisplayTable->SetValue(pointId, ModelRegistrationType::kImportRMSError, sqrt(ModelToImageSqrError / numberOfValidValues));
        modelDisplayTable->SetValue(pointId, ModelRegistrationType::kImportImageSTD, sqrt(importVariance));
    }

    double CorticalBone::calculateImportOffset(vtkIdType pointId) {
        double offset, xPLocal, xPImported;

        if(importedModelIndex <= kEndostealRamp) {
            xPImported = importedParametersTable->GetValue(pointId, ModelRegistrationType::kXpParam).ToDouble();
        } else if (importedModelIndex == kHighResClassifier) {
            if(!isnan(importedParametersTable->GetValue(pointId, ClassifierTransform::kCBStart).ToDouble())) {
                xPImported = 0.5 * (importedParametersTable->GetValue(pointId, ClassifierTransform::kSTEnd).ToDouble() + importedParametersTable->GetValue(pointId, ClassifierTransform::kCBStart).ToDouble());
            } else { // if no CB
                xPImported = importedParametersTable->GetValue(pointId, ClassifierTransform::kSTEnd).ToDouble();
            }
        } else {
            xPImported = 0;
            cerr<<"Error: In CorticalBone::calculateImportOffset incorrect 'Imported Model Index' = "<<importedModelIndex<<endl;
        }

        if(modelIndex <= kEndostealRamp) {
            xPLocal = modelParametersTable->GetValue(pointId, ModelRegistrationType::kXpParam).ToDouble();
        } else if (modelIndex == kHighResClassifier) {
            if(!isnan(modelParametersTable->GetValue(pointId, ClassifierTransform::kCBStart).ToDouble())) {
                xPLocal = 0.5 * (modelParametersTable->GetValue(pointId, ClassifierTransform::kSTEnd).ToDouble()
                                 + modelParametersTable->GetValue(pointId, ClassifierTransform::kCBStart).ToDouble());
            } else {
                xPLocal = modelParametersTable->GetValue(pointId, ClassifierTransform::kSTEnd).ToDouble();
            }
        } else {
            xPLocal = 0;
            cerr<<"Error: In CorticalBone::calculateImportOffset incorrect 'Model Index' = "<<modelIndex<<endl;
        }

        offset = xPImported - xPLocal;

        return offset;
    }

    bool CorticalBone::isImportedModelValid(vtkIdType pointId) {

        if(importedModelIndex == kHighResClassifier) { // TODO - check validity of imported threshold models
            return true;
        } else if(importedModelIndex == kCalibration) {
            return false;
        }

        // store current registration settings
        int currentModelIndex = registrationMethod->GetModelSelection();
        int localMaxDensity = registrationMethod->GetMaxDensity();
        ParametersType currentParameters = registrationMethod->GetCombinedParameters();
        double localP0, localP1, localP2; modelProfile->GetCalibration(localP0, localP1, localP2);

        // get imported parameters
        setModelToImport(pointId);

        // test validity
        bool validity = registrationMethod->IsValid();

        // restore to local state
        registrationMethod->SetModelSelection(currentModelIndex);
        modelProfile->SetCalibration(localP0, localP1, localP2);
        registrationMethod->SetMaxBMDDensity(localMaxDensity);
        registrationMethod->SetCombinedParameters(currentParameters);

        return validity;
    }

    void CorticalBone::clearPreviousResults() {

        // reset measurements & run status
        clearMeshAttributeTable();
        clearParametersTable();
        clearImageProfilesArray();
        resetPeriostealMesh();

        modelFitStatusArray.SetSize(mesh->GetNumberOfPoints());
        modelFitStatusArray.Fill(kUndefined);
        nanCount = 0; failureCount = 0; outOfBoundsCount = 0;

        // preset any masked out values
        maskResultArrays();

        // reset display state - which leads display to be reset to blank
        meshMeasured = false;
        displayIndex = kInvalid;

        loadedResults = false; // no loaded results

    }

    void CorticalBone::clearCalibration() {

        // result calibration values and state
        calibrationSet = false;

        calP1 = 0.0; calP0 = 1.0; calP2 = 0.0;
    }

    /* private methods */
    void CorticalBone::initialiseSettings() {

        // setup state
        imageSet = false; meshSet = false; meshMeasured = false; calibrationSet = false; phantomChanged = false;
        parametersImported = false; sampleNumberFixed = false; smoothingSet = false;

        loadedResults=false; // none loaded to begin with

        // set default values
        displayIndex = kInvalid; // set to blank initially
        fittingSchemeIndex = kStdAFitting;
        modelIndex = kThreeTierRect; // set to rect initially
        phantomIndex = kNoCal; // set to blank initially

        importedModelIndex = kUndefined;

        // profile initial values - set properly in setSampling()
        numberOfSamples = 500;
        profileLength = 20;//18;
        maximumOffset = 3;
        defaultPeriostealEdgeRatio = 0.3;//1.0/3.0;

        smoothingRadius = nan("1");

        thresholdWeight = nan("1");

        // create parameter and display arrays
        modelParametersTable = vtkSmartPointer<vtkTable>::New();
        modelDisplayTable = vtkSmartPointer<vtkTable>::New();

        // calibration values
        calP1 = 1.0; calP0 = calP2 = 0.0; // default to linear relationship
        cntrlRadius = nan("1"); cntrlNumber = 0; calNumber = 0;

        cntrlPts = NULL; cntrlPts = NULL; calHUVals = calBMDVals = NULL; cntrlVals = NULL;


        // fixed cb value
        globalFixedCB = nan("1");

        // free up image and mesh memory
        image = NULL; mesh = NULL; periostealMesh = NULL;
        meshTree = NULL; meshNormals = NULL;

        // reset mask count
        maskedCount=0;

    }

    /* Mesh Loading */
    bool CorticalBone::openVRML(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn){

        // import VRML
        vtkSmartPointer<vtkVRMLImporter> importer = vtkSmartPointer<vtkVRMLImporter>::New();
        importer->SetFileName ( fileName.c_str() );
        importer->Read();
        importer->Update();

        // convert to vtkDataSet
        vtkSmartPointer<vtkActorCollection> actors = importer->GetRenderer()->GetActors();
        actors->InitTraversal();
        vtkSmartPointer<vtkDataSet> pDataset = actors->GetNextActor()->GetMapper()->GetInput();

        //Convert to vtkPolyData
        vtkSmartPointer<vtkPolyData> localMesh = vtkPolyData::SafeDownCast ( pDataset );
        meshIn->DeepCopy(localMesh);
        return true;

    }

    bool CorticalBone::openVRML_SW(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn) {
        
        bool status = true;

        // create point, normal and poly objs to store file values in
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkDataArray> normals = vtkSmartPointer<vtkFloatArray>::New();
        vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();

        // create io stream
        std::string currLine;
        ifstream swVRMLFile (fileName.c_str());

        // set offset correct for sw file
        std::string offsetMessage;
        float swOffset[3] = {0,0,0};
        if(imageSet==true) {
            offsetMessage = "IMPORTANT: assumes wrl file created by Stradwin, so an offset correction will be made.";
            for(unsigned int i=0; i<2; i++) { // no z offset. hence upper limit of 1
                swOffset[i] = -1*image->GetSpacing()[i]/2;
            }
        } else {
            offsetMessage = "WARNING: wrl file opened before DICOM. No sw offset correction made.";
        }

        if (swVRMLFile.is_open())  {
            while ( getline (swVRMLFile,currLine) ) {

                // check for start of points section
                if(currLine.find("coord Coordinate") != std::string::npos) {
                    //cerr<<" Entering coordinate section"<<endl;

                    int ptCount=0;

                    // cycle through points and store
                    while ( getline (swVRMLFile,currLine) ) {

                        // check for end of points section
                        if(currLine.find("] }") != std::string::npos) {
                            break;
                        }

                        ptCount++;

                        // trim leading zeros
                        currLine = currLine.substr(currLine.find_first_not_of(' '));

                        // parse for three numbers
                        float point[3];
                        std::istringstream strStream( currLine );
                        for(int i=0; i<3;i++) {        // See the WARNING above for WHY we're doing this!
                            std::string var;
                            getline( strStream, var, ' ' );
                            //cerr << var << ", ";
                            point[i]= std::atof(var.c_str())+swOffset[i];
                        }
                        //cerr<<point[0]<<", "<<point[1]<<", "<<point[2]<<"; Pt index = "<<ptCount<<endl;
                        points->InsertNextPoint(point[0], point[1], point[2]);

                    }
                }
                    // check for start of normals section
                else if(currLine.find("normal Normal") != std::string::npos) {

                    //cerr<<" Entering normal section"<<endl;

                    normals->SetNumberOfComponents(3);
                    //normals->SetNumberOfTuples(points->GetNumberOfPoints()); do not allocate right now as insert values

                    // cycle through normals and store
                    while ( getline (swVRMLFile,currLine) ) {

                        // check for end of points section
                        if(currLine.find("] }") != std::string::npos) {
                            break;
                        }

                        // trim leading zeros
                        currLine = currLine.substr(currLine.find_first_not_of(' '));

                        // parse for three numbers
                        float normal[3];
                        std::istringstream strStream( currLine );
                        for(int i=0; i<3;i++) {        // See the WARNING above for WHY we're doing this!
                            std::string var;
                            getline( strStream, var, ' ' );
                            var=var.substr(0,var.find_first_of(','));
                            //cerr << var << ", ";
                            normal[i]= std::atof(var.c_str());
                        }
                        //cerr<<endl;
                        normals->InsertNextTuple(normal);
                        //normals->InsertNextPoint(normal[0], normal[1], normal[2]);
                    }
                }
                    // check for start of triangle section
                else if(currLine.find("coordIndex") != std::string::npos) {

                    //cerr<<" Entering index section"<<endl;
                    // cycle through indices and store
                    while ( getline (swVRMLFile,currLine) ) {

                        // check for end of points section
                        if(currLine.find("]") != std::string::npos) {
                            break;
                        }

                        // trim leading zeros
                        currLine = currLine.substr(currLine.find_first_not_of(' '));

                        // parse for three numbers
                        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                        std::istringstream strStream( currLine );
                        for(int i=0; i<3;i++) {        // See the WARNING above for WHY we're doing this!
                            std::string var;
                            getline( strStream, var, ',' );
                            //cerr << var << ", ";
                            //index[i]= std::atof(var);
                            triangle->GetPointIds()->SetId ( i, std::atof(var.c_str()) );
                        }
                        //cerr<<endl;
                        polys->InsertNextCell ( triangle );
                    }
                }
            }
            swVRMLFile.close();
        }


        // check reletive lengths of the points and normals
        if(points->GetNumberOfPoints()!=normals->GetNumberOfTuples()) {
            cerr<<"ERROR in CorticalBone::openVRML_SW() number of normals and points do not match."<<endl;
            status = false;
        }

        // Add the geometry and topology to the mesh
        meshIn->SetPoints ( points );
        meshIn->SetPolys ( polys );
        meshIn->GetPointData()->SetNormals(normals);

        if(verbose){cout<<Utilities::getTabString()<<"Opened a '.wrl' mesh file with "<<meshIn->GetPoints()->GetNumberOfPoints()<<" vertices."<<endl;}
        if(verbose){cout<<Utilities::getTabString()<<Utilities::getTabString()<<offsetMessage<<endl;}
        
        return status;
        
    }

    bool CorticalBone::openSTL(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn){

        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName(fileName.c_str());
        reader->Update();

        vtkSmartPointer<vtkPolyData> localMesh = reader->GetOutput();
        meshIn->DeepCopy(localMesh);


        if(verbose){cout<<Utilities::getTabString()<<"Opened a '.stl' mesh file with "<<meshIn->GetPoints()->GetNumberOfPoints()<<" vertices."<<endl;}
        return true;
    }

    bool CorticalBone::openPLY(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn) {

        int status = vtkPLYReader::CanReadFile(fileName.c_str());


        vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
        reader->SetFileName(fileName.c_str());
        reader->Update();

        vtkSmartPointer<vtkPolyData> localMesh = reader->GetOutput();
        meshIn->DeepCopy(localMesh);

        if(verbose){cout<<Utilities::getTabString()<<"Opened a '.ply' mesh file with "<<meshIn->GetPoints()->GetNumberOfPoints()<<" vertices."<<endl;}
        return true;

    }

    bool CorticalBone::openOBJ(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn) {

        vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
        reader->SetFileName(fileName.c_str());
        reader->Update();

        vtkSmartPointer<vtkPolyData> localMesh = reader->GetOutput();
        meshIn->DeepCopy(localMesh);

        if(verbose){cout<<Utilities::getTabString()<<"Opened a '.obj' mesh file with "<<meshIn->GetPoints()->GetNumberOfPoints()<<" vertices."<<endl;}
        
        return true;

    }

    bool CorticalBone::openMeshMask(std::string fileName) {

        // create mesh mask
        int n = mesh->GetNumberOfPoints();
        meshMask.SetSize(n); meshMask.Fill(kValid); // initially assume nothing masked out

        bool status;

        // check for and try open a sel file - preference given to sel file
        std::string maskName = fileName.substr(0, fileName.rfind('.')) + "_Mask.sel";
        status = openMaskSelFile(maskName);

        if(status==true) {
            return status;
        }

        // check for and try open a txt file
        maskName = fileName.substr(0, fileName.rfind('.')) + "_Mask.txt";
        status = openMaskTxtFile(maskName);
        return status;

    }

    bool CorticalBone::openMaskTxtFile(std::string fileName) {

        vtkIdType n = mesh->GetNumberOfPoints(); maskedCount=0;

        // check if there is a test file and try open
        ifstream maskFile (fileName.c_str());

        // read in values
        if (maskFile.good() && maskFile.is_open())  {

            std::string currLine;

            // open file and load in values
            getline(maskFile, currLine); // skip first comment line
            for (int i=0; i<n; i++) {
                getline(maskFile, currLine);

                int value = std::atoi(currLine.c_str());
                if(value==0) {
                    meshMask[i] = kMasked; // masked
                    maskedCount++;
                } else if(value==2) {
                    meshMask[i] = kDblPeak;
                }
            }
            if(verbose) {cout<<Utilities::getTabString()<<"Opened a '.txt' mesh mask file."<<endl;}
            return true;
        } else {
            return false;
        }


    }

    bool CorticalBone::openMaskSelFile(std::string fileName) { // output of wxRegSurf I think

        vtkIdType n = mesh->GetNumberOfPoints(); maskedCount=0;

        // check if there is a sel file and try open
        ifstream maskFile (fileName.c_str());
        if(maskFile.good() && maskFile.is_open()) {

            std::string currLine; getline(maskFile, currLine); // skip first comment line

            currLine = currLine.substr(currLine.find_first_not_of(' '));

            // parse single line for all numbers
            std::istringstream strStream( currLine );
            for(int i=0; i<n;i++) {
                std::string var;
                getline( strStream, var, ' ' );
                //cerr << var << ", ";
                int value=std::atoi(var.c_str());
                if(value==0) { // sel file equals zero if not selected
                    meshMask[i] = kMasked; // masked
                    maskedCount++;
                } else if(value==2) {
                    meshMask[i] = kDblPeak;
                }
            }

            if(verbose) {cout<<Utilities::getTabString()<<"Opened a '.sel' mesh mask file."<<endl;}
            return true;
        }  else {
            return false;
        }
    }

    /* CT Image Loading */
    bool CorticalBone::openDICOM(std::string filePath) {
        // instantiate the series reader
        ImageReaderType::Pointer reader = ImageReaderType::New();

        // connect reader to a DICOM v3 class = GDCMImageIO
        ImageDICOMIOType::Pointer dicomIO = ImageDICOMIOType::New();
        reader->SetImageIO( dicomIO );

        // generates the DICOM slice names given the directory
        DICOMNamesGeneratorType::Pointer nameGenerator = DICOMNamesGeneratorType::New();
        cerr<<"DICOM directory="<<filePath.c_str()<<endl;
        nameGenerator->SetInputDirectory( filePath );

        // generate file names and past to reader
        FileNamesContainer fileNames = nameGenerator->GetInputFileNames();
        reader->SetFileNames( fileNames );

        // actually go about loading / reading in the image. Termed 'invoking'
        try {
            reader->Update();

            /*// http://public.kitware.com/pipermail/vtkusers/2014-August/084854.html

            // check for image orientation patient not aligned with x, y
            vtkSmartPointer<vtkDICOMImageReader> reader2d = vtkSmartPointer<vtkDICOMImageReader>::New();
            reader2d->SetFileName(filePath.c_str()); reader2d->Update();

            float *orientation = reader2d->GetImageOrientationPatient();
            float *xdir = &orientation[0]; float *ydir = &orientation[3];
            float zdir[3]; vtkMath::Cross(xdir, ydir, zdir);

            if(xdir[kX]!=1.0 && ydir[kY]!=1.0) {
                itk::ChangeInformationImageFilter< ImageType >::Pointer filter = itk::ChangeInformationImageFilter< ImageType >::New();
                filter->SetInput( reader->GetOutput() );

                const ImageType::DirectionType direction = image->GetDirection();
                const ImageType::DirectionType newDirection;
                newDirection[0][kX]=xdir[kX];
                        = direction * rotation.GetMatrix();
                filter->SetOutputDirection( newDirection );
                filter->ChangeDirectionOn();

            }



*/


            // create an image to process
            image = reader->GetOutput();

            ImageType::DirectionType directionMatrix = image->GetDirection();

            if(directionMatrix[0][kX]!=1.0 && directionMatrix[1][kY]!=1.0) {
                directionMatrix[0][kX]=1.0; directionMatrix[0][kY]=0.0; directionMatrix[0][kZ]=0.0;
                directionMatrix[1][kX]=0.0; directionMatrix[1][kY]=1.0; directionMatrix[1][kZ]=0.0;
                directionMatrix[2][kX]=0.0; directionMatrix[2][kY]=0.0; directionMatrix[2][kZ]=1.0;
                image->SetDirection(directionMatrix);
            }

            return true;
        } catch (itk::ExceptionObject &ex) {
            std::clog << ex << std::endl;
            //return EXIT_FAILURE;
            return false;
        }
    }

    bool CorticalBone::openRAW(std::string fileName) {

        ImageReaderType::Pointer reader = ImageReaderType::New();
        reader->SetFileName(fileName);

        try {
            reader->Update();
        } catch (itk::ExceptionObject &excp) {
            std::cerr << "Exception thrown while reading RAW image file: " << excp.GetDescription() <<endl;
            return false; // error code
        }

        image = reader->GetOutput();
        return true;
    }

    bool CorticalBone::openQCT(std::string fileName) {

        // todo - read image

        /*ImageQCTIOType::Pointer io = ImageQCTIOType::New();

        ImageReaderType::Pointer reader = ImageReaderType::New();
        reader->SetFileName(fileName);
        reader->SetImageIO(io);

        try {
            reader->Update();
        } catch (itk::ExceptionObject &excp) {
            std::cerr << "Exception thrown while reading RAW image file: " << excp.GetDescription() <<endl;
            return false; // error code
        }

        image = reader->GetOutput();
        return true;*/

        ImageType::IndexType start = {0,0,0};
        ImageType::SizeType size;
        ImageType::PointType origin;
        origin.Fill(0.0);
        ImageType::SpacingType spacing;
        spacing[0]=1.0; spacing[1]=1.0; spacing[2]=1.0; //={1.0,1.0,1.0};
        PixelType qctToHU_offset=-1000;

        // read in binary file
        ifstream imageFile (fileName.c_str(), ios::out | ios::binary);

        // read in values
        if (!imageFile.good() && !imageFile.is_open())  {
            return false;
        }

        std::string fileNameOut = fileName.substr(0, fileName.rfind('.')) + ".txt";

        //-------- read in header --------//
        // check endien-ness
        bool swapEndien = false;
        vtkTypeUInt32 integer32;
        imageFile.seekg(0,ios::beg);
        imageFile.read((char*)&integer32,sizeof(integer32));
        if(integer32 !=0xface424b) {
            swapEndien=true;
        }

        // check for 1 at byte 32
        imageFile.seekg(32);
        imageFile.read((char*)&integer32,sizeof(vtkTypeUInt32));
        if(swapEndien) {integer32=Utilities::swap32(integer32);}
        if(integer32 !=1) {
            cerr<<"Failed to open QCT file. Incorrect header format. Byte 32 != 1"<<endl;
            imageFile.close(); return false;
        }

        // get number of frames
        imageFile.seekg(36);
        imageFile.read((char*)&integer32,sizeof(vtkTypeUInt32));
        if(swapEndien) {integer32=Utilities::swap32(integer32);}
        size[2]= integer32;
        if ((size[2] < 1) || (size[2] > 2048)) {
            cerr<<"Failed to open QCT file. Invalid number of slices="<<size[2]<<endl;
            imageFile.close(); return false;
        }
        //cerr<<"Number of frames="<<size[2]<<endl;

        // get the fame data locations
        itk::Array<vtkTypeUInt32> sliceMemoeryOffsets((unsigned int)size[2]);
        for (int i=0; i<size[2]; i++) {
            imageFile.read((char*)&sliceMemoeryOffsets[i],sizeof(vtkTypeUInt32));
            if(swapEndien) {sliceMemoeryOffsets[i]=Utilities::swap32(sliceMemoeryOffsets[i]);}
            //cerr<<i<<": "<<sliceMemoeryOffsets[i]<<endl;
        }

        //------- read in slice header ------// // todo - rigerously check each slice has matching headers
        bool stackedUp=true;
        int sliceIndex=0;
        // get slice size
        imageFile.seekg(sliceMemoeryOffsets[sliceIndex]);
        imageFile.read((char*)&integer32,sizeof(vtkTypeUInt32));
        if(swapEndien) {integer32=Utilities::swap32(integer32);}
        if (integer32 !=1) {
            cerr<<"Failed to open QCT file. Invalid header for slice="<<sliceIndex<<endl;
            imageFile.close(); return false;
        }
        imageFile.read((char*)&integer32,sizeof(vtkTypeUInt32));
        if(swapEndien) {integer32=Utilities::swap32(integer32);}
        if (integer32 !=0) {
            cerr<<"Failed to open QCT file. Invalid header for slice="<<sliceIndex<<endl;
            imageFile.close(); return false;
        }

        // read numbers of pixels
        imageFile.read((char*)&integer32,sizeof(vtkTypeUInt32));
        if(swapEndien) {integer32=Utilities::swap32(integer32);} size[0]= integer32;
        imageFile.read((char*)&integer32,sizeof(vtkTypeUInt32));
        if(swapEndien) {integer32=Utilities::swap32(integer32);} size[1]= integer32;
        if (size[0] < 1 || size[1] < 1) {
            cerr<<"Failed to open QCT file. Invalid number of pixels in slice header="<<sliceIndex<<endl;
            imageFile.close(); return false;
        }
        //cerr<<"Number of pixels=["<<size[0]<<","<<size[1]<<"]"<<endl;

        // read in the scales
        vtkTypeFloat32 float32;
        imageFile.seekg(sliceMemoeryOffsets[sliceIndex]+sizeof(vtkTypeUInt32)*6);
        imageFile.read((char*)&float32,sizeof(vtkTypeFloat32));
        if(swapEndien) {float32=Utilities::swap32(float32);} spacing[0]=float32;
        imageFile.read((char*)&float32,sizeof(vtkTypeFloat32));
        if(swapEndien) {float32=Utilities::swap32(float32);} spacing[1]=float32;
        if (spacing[0] < 0 || spacing[1] < 0) {
            cerr<<"Failed to open QCT file. Invalid spacing in slice header="<<sliceIndex<<endl;
            imageFile.close(); return false;
        }

        // read in z spacing
        vtkTypeFloat32 zSliceMinimum;
        imageFile.seekg(sliceMemoeryOffsets[sliceIndex]+sizeof(vtkTypeUInt32)*10); // slice 0
        imageFile.read((char*)&zSliceMinimum,sizeof(vtkTypeFloat32));
        if(swapEndien) {zSliceMinimum=Utilities::swap32(zSliceMinimum);}
        imageFile.seekg(sliceMemoeryOffsets[sliceIndex+1]+sizeof(vtkTypeUInt32)*10); // slice 1
        imageFile.read((char*)&float32,sizeof(vtkTypeFloat32));
        if(swapEndien) {float32=Utilities::swap32(float32);} spacing[2]=float32-zSliceMinimum; // /10?
        if(spacing[2]<0) {
            stackedUp=false; spacing[2]=fabs(spacing[2]);
            imageFile.seekg(sliceMemoeryOffsets[size[2]-1]+sizeof(vtkTypeUInt32)*10);
            imageFile.read((char*)&zSliceMinimum,sizeof(vtkTypeFloat32));
            if(swapEndien) {zSliceMinimum=Utilities::swap32(zSliceMinimum);}
        }
        if (spacing[0] <= 0 || spacing[1] <= 0 || spacing[2] <=0) {
            cerr<<"Failed to open QCT file. Invalid z spacing between slice headers="<<sliceIndex<<" and "<<sliceIndex+1<<endl;
            imageFile.close(); return false;
        }
        //cerr<<"Image resolution=["<<spacing[0]<<","<<spacing[1]<<","<<spacing[2]<<"]"<<endl;
        origin[2]=zSliceMinimum;
        //cerr<<"Image origin=["<<origin[0]<<","<<origin[1]<<","<<origin[2]<<"]"<<endl;



        // set up image size/dimensions
        image = ImageType::New();

        ImageType::RegionType region;
        region.SetSize( size );
        region.SetIndex( start );
        image->SetRegions( region );

        image->SetOrigin( origin );
        image->SetSpacing( spacing );

        /*// output slice header
        cerr<<" slice 0: "<<endl;
        imageFile.seekg(sliceMemoeryOffsets[0]); cerr<<std::dec;
        for (int i=0; i<6; i++) {
            imageFile.read((char*)&integer32,sizeof(vtkTypeUInt32));
            cerr<<i<<": "<< integer32 <<endl;
        }
        for (int i=0; i<12; i++) {
            imageFile.read((char*)&float32,sizeof(vtkTypeFloat32));
            cerr<<i+7<<": "<<float32<<endl;
        } cerr<<endl<<endl;
        for (int i=0; i<14; i++) {
            imageFile.read((char*)&integer32,sizeof(vtkTypeUInt32));
            cerr<<i+7<<": "<<integer32<<endl;
        } cerr<<endl<<endl;
        vtkTypeUInt16 integer16;
        for (int i=0; i<size[0]*size[1]; i++) {
            imageFile.read((char*)&integer16,sizeof(vtkTypeUInt16));
            cerr<<integer16<<", "; //cerr<<i<<": "<<integer32<<endl;
            if(i%size[0]==0) {cerr<<endl;}
        } cerr<<endl<<endl;*/


        // set up image
        image->Allocate();
        sliceIndex = stackedUp ? 0 : (int)size[2]-1;
        for(int k=0; k<size[2]; k++) {

            imageFile.clear(); // as if it reaches the end of file it will nolonger seek correctly
            imageFile.seekg(sliceMemoeryOffsets[sliceIndex]+sizeof(vtkTypeUInt32)*34, ios::beg);
            for(int j=0; j<size[0]; j++) { // todo - address: seemingly reversed order to deal with image & sw mesh mess-match
                for(int i=0; i<size[1]; i++) {
                    ImageType::IndexType pixelIndex = {i,j,k}; // {i,j,k};

                    vtkTypeUInt16 integer16; imageFile.read((char*)&integer16,sizeof(vtkTypeUInt16));
                    if(swapEndien) {integer16=Utilities::swap16(integer16);}

                    image->SetPixel(pixelIndex, (PixelType)integer16+qctToHU_offset);
                }
            }
            sliceIndex = stackedUp ? sliceIndex+1 : sliceIndex-1;
        }

        imageFile.close();
        return true;

        //( wxFileName file, long frame, long real_frame, unsigned short *data, float *pos )


        // create class based off: class RawImageIO:public ImageIOBase

        /*image = ImageType::New();

        ImageType::IndexType start;
        start[0] = 0;
        start[1] = 0;
        start[2] = 0;

        ImageType::SizeType size;
        size[0] = 200;
        size[1] = 200;
        size[2] = 200;

        ImageType::RegionType region;
        region.SetSize( size );
        region.SetIndex( start );
        image->SetRegions( region );
        image->Allocate();

         for(int i=0; i<size[0]; i++) {
            for(int j=0; j<size[1]; j++) {
                for(int k=0; k<size[2]; k++) {
                    ImageType::IndexType pixelIndex = {i,j,k};
                    image->SetPixel(pixelIndex, i*j*k);
                }
            }
        }*/

        /*    unsigned short buffer[2048];
            qct_frame_t header;
            FILE *data_file;
            data32 l, frames, i, d, *offset;
            int val, start_f, day, month, year, bday, bmonth, byear;
            size_t sz;

            // Open file and read frame number and offset locations
            if ( NULL == (data_file = fopen(fileName, "rb")) ) return 0;
            sz = fread(buffer, sizeof(char), 32, data_file); // skip first 32 bytes
            sz = fread(&l, sizeof(data32), 1, data_file); // next value should be 1
            if (l != 1) {
                fclose(data_file);
                return 0;
            }
            sz = fread(&frames, sizeof(data32), 1, data_file); // next value should be frames
            if ((frames < 1) || (frames > 2048)) {
                fclose(data_file);
                return 0;
            }
            offset = new data32[frames];
            sz = fread(offset, sizeof(data32), frames, data_file); // following values are offsets into frame + header

            // Pick out the date of birth, which is forward 140 bytes
            sz = fread(buffer, sizeof(char), 140, data_file); // skip next 300 bytes
            sz = fread(&d, sizeof(data32), 1, data_file);
            bday =   (d & 0x0000001F);
            bmonth = (d & 0x000001E0) >> 5;
            byear =  (d & 0xFFFFFE00) >> 9;

            // Also pick out the date, which is forward 300 bytes (i.e. another 156)
            sz = fread(buffer, sizeof(char), 156, data_file); // skip next 300 bytes
            sz = fread(&d, sizeof(data32), 1, data_file);
            day =   (d & 0x0000001F);
            month = (d & 0x000001E0) >> 5;
            year =  (d & 0xFFFFFE00) >> 9;

            // Keep looking for the frame header until we get to a frame which has header.l3 set to 1
            for (start_f=0; start_f<frames; start_f++) {
                fseek(data_file, offset[start_f], SEEK_SET);
                sz = fread(&header, sizeof(qct_frame_t), 1, data_file);
                if (header.l3 == 1) break;
            }

            // Check we have enough frames for requested frame
            if ((start_f+real_frame) >= frames) {
                delete[] offset;
                fclose(data_file);
                return 0;
            }

            // Possibly skip to the frame we really want
            if (real_frame > 0) {
                fseek(data_file, offset[start_f+real_frame], SEEK_SET);
                sz = fread(&header, sizeof(qct_frame_t), 1, data_file);
            }

            // Check for reasonable header values
            if ((header.l1 != 1) || (header.xpixels<0) || (header.ypixels<0) || (header.xscale<0.0) || (header.yscale<0.0)) {
                delete[] offset;
                fclose(data_file);
                return 0;
            }

            // Return now with a positive value if we are just checking for a credible QCT file
            if (real_frame < 0) {
                fclose(data_file);
                return 1;
            }

            if (data == NULL) {

                // Possibly set up all the relevant stradwin resources
                // We assume everything else has been set to default
                res->buf_dicom->write_no_changeprocs(true);
                res->buf_width->write_no_changeprocs((long)header.xpixels);
                res->buf_height->write_no_changeprocs((long)header.ypixels);
                res->pos_rec->write_no_changeprocs(true);
                res->bin_im_filename->write_no_changeprocs(file.GetFullPath());
                res->vid_xpos->write_no_changeprocs((long)0);
                res->vid_ypos->write_no_changeprocs((long)0);
                res->xscale->write_no_changeprocs(header.xscale/10.0);
                res->yscale->write_no_changeprocs(header.yscale/10.0);
                res->xtrans->write_no_changeprocs(0.0);
                res->ytrans->write_no_changeprocs(0.0);
                res->ztrans->write_no_changeprocs(0.0);
                res->azimuth->write_no_changeprocs(0.0);
                res->elevation->write_no_changeprocs(0.0);
                res->roll->write_no_changeprocs(0.0);
                res->buf_frames->write_no_changeprocs(frames-start_f);

            } else {

                // Get the data
                dicom_has_pos = true;
                dicom.valid = true;
                dicom.byte = false;
                dicom.rf = false;
                dicom.hum = 1.0;
                dicom.hub = -1000.0;
                dicom_content[frame] = 0;
                dicom.thick = header.thickness/10.0;
                dicom_pos[6*frame+0] = 0.0;
                dicom_pos[6*frame+1] = 0.0;
                dicom_pos[6*frame+2] = header.offset/10.0;
                dicom_pos[6*frame+3] = 0.0;
                dicom_pos[6*frame+4] = 0.0;
                dicom_pos[6*frame+5] = 0.0;
                dicom_xscale[frame] = header.xscale;
                dicom_yscale[frame] = header.yscale;
                if (pos != NULL) {
                    for (int n=0; n<6; n++) pos[n] = dicom_pos[6*frame+n];
                }
                if ((header.xpixels==res->buf_width->read_proposed()) && (header.ypixels==res->buf_height->read_proposed())) {
                    for (int y=0; y<header.ypixels; y++) {
                        sz = fread(data+y*header.xpixels, 2, header.xpixels, data_file);
                    }
                }
                if (frame == 0) res->alert_string->write_no_changeprocs(wxString::Format("Scanned on %02d/%02d/%04d, DOB %02d/%02d/%04d", day, month, year, bday, bmonth, byear));
            }

            delete[] offset;
            fclose(data_file);
            return header.xpixels * header.ypixels;


            */

        /* cout<<"corticalBone::openQCT() not yet implemented."<<endl;
        return false; */


    }

    void CorticalBone::openTIF(std::vector<std::string> fileNames, double zSpacing, double zOffset) {

        // define types for tiff images
        typedef unsigned short TiffPixelType;
        typedef itk::Image< TiffPixelType, Dimension > TiffImageType;
        typedef itk::ImageSeriesReader< TiffImageType > TiffImageReaderType;
        //typedef itk::CastImageFilter<TiffImageType, ImageType> CastImageFilterType;
        typedef itk::RescaleIntensityImageFilter<TiffImageType, ImageType> CastImageFilterType;

        // instantiate the series reader - Tiff
        TiffImageReaderType::Pointer reader = TiffImageReaderType::New();

        // set file names
        reader->SetFileNames( fileNames );

        // actually go about loading / reading in the image. Termed 'invoking'
        try {
            reader->Update();
        } catch (itk::ExceptionObject &ex) {
            std::clog << ex << std::endl;
            //return EXIT_FAILURE;
        }

        //// tiff image
        //TiffImageType::Pointer tiffImage = reader->GetOutput();
        //cerr<<"Tiff image loaded."<<endl;

        // image casting filter
        CastImageFilterType::Pointer filter = CastImageFilterType::New();
        filter->SetInput( reader->GetOutput() );
        filter->Update();

        // create an image to process
        image = filter->GetOutput();
        cerr<<"Tiff image converted."<<endl;

//        // manually populate image
//        image = ImageType::New();
//        image->SetRegions(tiffImage->GetLargestPossibleRegion());
//        image->Allocate();

        // set spacing
        ImageType::SpacingType imageSpacing = image->GetSpacing(); //tiffImage->GetSpacing();
        if(isnan(zSpacing)) { // if the z spacing is not specified set isotropic
            imageSpacing[2]=(imageSpacing[0]+imageSpacing[1])/2.0;
            if(imageSpacing[0]==imageSpacing[1]) { // isotropic
                cout<<Utilities::getTabString()<<"Warning: No z spacing included with the .tiff image stack so set isotropic"<<endl;
            } else {
                cout<<Utilities::getTabString()<<"Warning: No z spacing included with the .tiff image stack so set as the x and y spacing average"<<endl;
            }
        } else {
            imageSpacing[2]=zSpacing;
        }
        image->SetSpacing(imageSpacing);

        // set origin
        //ImageType::PointType imageOrigin = tiffImage->GetOrigin();
        if(!isnan(zOffset)) { // only set origin if not nan
            ImageType::PointType imageOrigin = image->GetOrigin();
            imageOrigin[2]=zOffset;
            image->SetOrigin(imageOrigin);
        } else {
            //cout<<Utilities::getTabString()<<"Warning: No z offset provided."<<endl;
        }
        //image->SetOrigin(imageOrigin);


//        // manually convert image as cast filter fails to properly convert - todo use iamge iterators instead
//        ImageRegion<Dimension>::SizeType imageSize = image->GetLargestPossibleRegion().GetSize( );
//        for(int i=0; i<imageSize[0]; i++) {
//            for(int j=0; j<imageSize[1]; j++) {
//                for(int k=0; k<imageSize[2]; k++) {
//                    ImageType::IndexType index = {{i,j,k}};
//                   // PixelType signedValue = ((PixelType)tiffImage->GetPixel(index))-SHRT_MAX/2;
//                    PixelType signedValue = (PixelType)(((long)tiffImage->GetPixel(index))-SHRT_MAX/2);
//                    //PixelType signedValue = (PixelType)(tiffImage->GetPixel(index)/2);
//                    image->SetPixel(index, signedValue);
//                }
//            }
//        }
//        cerr<<"Tiff image pixels converted."<<endl;

        return;
    }

    bool CorticalBone::setupImage() {

        imageSet=true; meshMeasured = false;


        // set up image interpolator
        Utilities::measureTime(true);
        interpolator = InterpolatorType::New();
        interpolator->SetInputImage(image);
        //m_interpolator->SetSplineOrder(3); // ideally set to order 3 but it runs like a dog. Todo - investigate performance improvements

        // create profile object - needs image and interpolator
        createModelProfile();


        // set up display
        ImageMinMaxCalculatorType::Pointer imageIntensityRange = ImageMinMaxCalculatorType::New();
        imageIntensityRange->SetImage(image);
        imageIntensityRange->ComputeMaximum();
        maxImageValue = imageIntensityRange->GetMaximum();

        if(verbose){cout<<Utilities::getTabString()<<Utilities::getTabString()<<Utilities::measureTime(false)<<" to create the image interpolator objects"<<endl;}


        ImageType::SpacingType imageSpacing = image->GetSpacing();
        cntrlRadius = 10 * (imageSpacing[0] + imageSpacing[1] + imageSpacing[2]);

        // set buffered region
        if(meshSet && imageSet) {
            setSampling();
            setImageBuffer();
            // setup model registration
            createRegistrationMethods();
        }
        return imageSet;

    }

    /* Save Methods */
    bool CorticalBone::saveVRML(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn) {

        // #include <vtkVRMLExporter.h>

        // import VRML
//        vtkSmartPointer<vtkVRMLExporter> exporter = vtkSmartPointer<vtkVRMLExporter>::New();
//        exporter->SetFileName ( fileName.c_str() ); // seems to require a render window its data input

        // manually save the wrl

        // get mesh data
        vtkSmartPointer<vtkPoints> points = meshIn->GetPoints();
        vtkSmartPointer<vtkCellArray> tris = meshIn->GetPolys();
        vtkSmartPointer<vtkDataArray> normals;
        if(meshIn->GetPointData()->GetArray("Normals")!=NULL) {
            normals = meshIn->GetPointData()->GetArray("Normals");
        } else if(meshIn->GetPointData()->GetNormals()!=NULL) {
            normals = meshIn->GetPointData()->GetNormals();
        } else { // create missing normal values
            vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
            normalGenerator->SetInputData(mesh);
            normalGenerator->ComputePointNormalsOn();
            normalGenerator->ComputeCellNormalsOff();
            normalGenerator->Update();
            meshIn = normalGenerator->GetOutput();
            normals = meshIn->GetPointData()->GetNormals();
        }
        int n = normals->GetNumberOfTuples();

        // set up file
        ofstream vrmlFile;
        vrmlFile.open (fileName.c_str());
        std::string comma = ",", space=" ";

        // write header information
        vrmlFile << "#VRML V2.0 utf8" << endl;
        vrmlFile << "# LAD sphere radius hint 5.0" << endl;
        vrmlFile << "Transform { children [ Shape {" << endl;
        vrmlFile << "appearance Appearance { material Material {" << endl;
        vrmlFile << "ambientIntensity 0.2" << endl;
        vrmlFile << "diffuseColor 1.000 1.000 0.200" << endl;
        vrmlFile << "specularColor 0.1 0.1 0.1" << endl;
        vrmlFile << "emissiveColor 0 0 0" << endl;
        vrmlFile << "shininess     0.2" << endl;
        vrmlFile << "transparency  0.200" << endl;
        vrmlFile << "} }" << endl;


         // write start of body
        vrmlFile << "geometry IndexedFaceSet {" << endl;

        // write start of coordinates - pts
        vrmlFile << "coord Coordinate { point [" << endl;

        for(int i=0; i<n; i++) {
            std::ostringstream lineStream; lineStream.precision(9);
            double point[Dimension]; points->GetPoint(i, point);

            lineStream << point[0] << space << point[1] << space << point[2];

            if(i!=n-1) {
                lineStream << comma;
            }

            vrmlFile << lineStream.str() << endl;
        }

        vrmlFile << "] }" << endl;

        // write start of normals
        vrmlFile << "normal Normal { vector [" << endl;

        for(int i=0; i<n; i++) {
            std::ostringstream lineStream; lineStream.precision(9);
            double normal[Dimension]; normals->GetTuple(i, normal);

            lineStream << normal[0] << space << normal[1] << space << normal[2];

            if(i!=n-1) {
                lineStream << comma;
            }

            vrmlFile << lineStream.str() << endl;
        }

        vrmlFile << "]" << endl;


        // write start of tris
        vrmlFile << "coordIndex [" << endl;
        vtkIdType nTris = tris->GetNumberOfCells();
        vtkIdType cellLocation = 0; // index at the start of each cell entry (i.e. all tri 0 verticies have id 0)
        for(int i=0; i<nTris; i++) {
            vtkIdType numIds; // to hold the size of the cell
            vtkIdType *tri; // to hold the ids in the cell
            tris->GetCell(cellLocation, numIds, tri);


            std::ostringstream lineStream; lineStream.precision(9);
            lineStream << tri[0] << comma << space << tri[1] << comma << space << tri[2] << comma << space << "-1";

            if(i!=n-1) {
                lineStream << comma;
            }

            vrmlFile << lineStream.str() << endl;

            cellLocation += 1 + numIds; // incrment to the start of the next cell
        }

        vrmlFile << "] }" << endl;

        // write end of body
        vrmlFile << "}" << endl;
        vrmlFile << "} ] }" << endl;


        // close file and return
        vrmlFile.close();

        return true;
    }

    bool CorticalBone::saveSTL(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn) {
        vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(meshIn);
        writer->Write();
        
        return true;
    }

    bool CorticalBone::saveOBJ(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn) {

        vtkSmartPointer<vtkOBJWriter> writer = vtkSmartPointer<vtkOBJWriter>::New();
        writer->SetInputData(meshIn);
        writer->SetFileName(fileName.c_str());
        writer->Update();

        return true;
    }

    /* Open / Save Parameter Values */
    bool CorticalBone::openImportParameters(std::string fileName) {


        // check header information
        char delimiter = ',';

        ifstream parametersFile (fileName.c_str());
        std::string header;
        if (!parametersFile.is_open())  {
            cerr<<"Error: Unable to open import file - "<<fileName.c_str()<<endl;
            return false;
        }

        // header line 1
        std::string headerLine1;
        getline (parametersFile, headerLine1);

        // check its a parameter file
        std::istringstream header1Stream( headerLine1 );
        std::string headerEntry;
        getline (header1Stream, headerEntry,delimiter); // ignore
        if(strcmp(headerEntry.c_str(), "Parameter Values")!=0) {
            cerr<<"Error: headerEntry="<<headerEntry<<endl;
            return false;
        }

        // get model type of import
        getline (header1Stream, headerEntry,delimiter);
        if(strcmp(headerEntry.c_str(), getModelName(kThreeTierRect).c_str())==0) {
            importedModelIndex = kThreeTierRect;
            if(verbose){cout<<Utilities::getTabString()<<"type: Rect Function"<<endl;}
        } else if(strcmp(headerEntry.c_str(), getModelName(kEndostealRamp).c_str())==0) {
            importedModelIndex = kEndostealRamp;
            if(verbose){cout<<Utilities::getTabString()<<"type: Ramp Function"<<endl;}
        } else if(strcmp(headerEntry.c_str(), getModelName(kFabricated).c_str())==0) {
            importedModelIndex = kFabricated;
            if(verbose){cout<<Utilities::getTabString()<<"type: Fabricated"<<endl;}
        } else if(strcmp(headerEntry.c_str(), getModelName(kHighResClassifier).c_str())==0) {
            importedModelIndex = kHighResClassifier;
            if(verbose){cout<<Utilities::getTabString()<<"type: Classifier"<<endl;}
        } else if(strcmp(headerEntry.c_str(), getModelName(kCalibration).c_str())==0) {
            importedModelIndex = kInvalid;
            cerr<<Utilities::getTabString()<<"type: Calibration Unsupported"<<endl;
            return false;
        } else {
            importedModelIndex = kInvalid;
            cerr<<Utilities::getTabString()<<"type: Unsupported"<<endl;
            return false;
        }

        // header line 2
        std::string headerLine2;
        getline (parametersFile, headerLine2);

        // get the number of parameters TODO - check it is the amount expected
        std::istringstream header2Stream( headerLine2 ); std::string columnName;
        while(std::getline(header2Stream, columnName, delimiter)) {
            vtkSmartPointer<vtkDoubleArray> columnArray = vtkSmartPointer<vtkDoubleArray>::New();
            columnArray->SetName(columnName.c_str());
            importedParametersTable->AddColumn(columnArray);
        }

        vtkIdType maxI = mesh->GetNumberOfPoints(); importedParametersTable->SetNumberOfRows(maxI);


        // read in parameter values
        bool status = readTableValues(parametersFile, importedParametersTable);
        parametersFile.close();

        return status;
    }

    bool CorticalBone::openParameters(std::string fileName) {


        char delimiter = ',';

        ifstream parametersFile (fileName.c_str());

        // check if file exists
        if(!parametersFile.good()) {
            return false; // silent as this is currently called each time a new project is opened irrespective of weither results files exist
        }

        std::string header;
        if (!parametersFile.is_open())  {
            cerr<<"Error in CorticalBone::openParameters(): Unable to open file - "<<fileName.c_str()<<endl;
            return false;
        }

        // check header line 1
        std::string headerLine1;
        getline (parametersFile, headerLine1);

        // check its a parameter file
        std::istringstream header1Stream( headerLine1 );
        std::string headerEntry;
        getline (header1Stream, headerEntry,delimiter);
        if(strcmp(headerEntry.c_str(), "Parameter Values")!=0) {
            cerr<<"Error in CorticalBone::openParameters(): headerEntry isn't 'Parameter Values', its:"<<headerEntry<<endl;
            return false;
        }

        // check  the model type matches the project
        getline (header1Stream, headerEntry,delimiter);
        if(strcmp(headerEntry.c_str(), getModelName(modelIndex).c_str())!=0) {
            cerr<<"Miss-match in CorticalBone::openParameters() between the model types specified in the parameters header file ("<<headerEntry<<") and the project model index ("<<getModelName(modelIndex)<<")"<<endl;
            return false;

        }
//        // check to make sure correct scheme
//        getline (header1Stream, headerEntry,delimiter);
//        if(strcmp(headerEntry.c_str(), getModelName(fittingSchemeIndex).c_str())!=0) {
//            cerr<<"Miss-match between the scheme types specified in the header files that of the project"<<endl;
//            return false;
//        }


        // check header line 2 - number of parameters
        std::string headerLine2;
        getline (parametersFile, headerLine2);

        std::istringstream header2Stream( headerLine2 ); std::string columnName; int columnNumber = 0;
        while(std::getline(header2Stream, columnName, delimiter)) {
            columnNumber++;
        }

        if(columnNumber!=modelParametersTable->GetNumberOfColumns()) {
            cerr<<"Miss-match in CorticalBone::openParameters() between the number of parameters specified in the header files ("<<columnNumber<<") and that of the project ("<<getNumberOfParameters()<<")"<<endl;
            return false;
        }


        // read in parameter values
        bool status = readTableValues(parametersFile, modelParametersTable);
        parametersFile.close();

        return status;
    }

    bool CorticalBone::openDisplays(std::string fileName) { // only results not imports


        // check header infomation
        char delimiter = ',';

        ifstream displaysFile(fileName.c_str());

        // check if file exists
        if(!displaysFile.good()) {
            return false; // silent as this is currently called each time a new project is opened irrespective of weither results files exist
        }

        std::string header;
        if (!displaysFile.is_open())  {
            cerr<<"Error: Unable to open file - "<<fileName.c_str()<<endl;
            return false;
        }

        // header line 1
        std::string headerLine1;
        getline (displaysFile, headerLine1);

        // check its a display file
        std::istringstream header1Stream( headerLine1 );
        std::string headerEntry;
        getline (header1Stream, headerEntry,delimiter); // ignore
        if(strcmp(headerEntry.c_str(), "Display Values")!=0) {
            cerr<<"Error: headerEntry="<<headerEntry<<endl;
            return false;
        }
        // check to make sure correct model type
        getline (header1Stream, headerEntry, delimiter);
        if(strcmp(headerEntry.c_str(), getModelName(modelIndex).c_str())!=0) {
            cerr<<"Miss-match between the model types specified in the header files ("<<headerEntry<<") and the project model index ("<<modelIndex<<")"<<endl;
            return false;
        }
//        // check to make sure correct scheme
//        getline (header1Stream, headerEntry,delimiter);
//        if(strcmp(headerEntry.c_str(), getModelName(fittingSchemeIndex).c_str())==0) {
//            cerr<<"Miss-match between the scheme types specified in the header files that of the project"<<endl;
//            return false;
//        }

        std::string headerLine2;
        getline (displaysFile, headerLine2);

        // get the number of displays and ensure it matches the selection
        std::istringstream header2Stream( headerLine2 ); std::string columnName; int columnNumber = 0;
        while(std::getline(header2Stream, columnName, delimiter)) {
            columnNumber++;
        }

        if(columnNumber!= modelDisplayTable->GetNumberOfColumns()) {
            cerr<<"Miss-match between the number of displays specified in the header files ("<<columnNumber<<") and that of the project ("<<
                                                                                                             modelDisplayTable->GetNumberOfColumns()<<")"<<endl;
            return false;
        }

        bool status = readTableValues(displaysFile, modelDisplayTable);
        displaysFile.close();

        return status;
    }

    bool CorticalBone::openProfiles(std::string fileName, vtkSmartPointer<vtkDoubleArray> &array, std::string name,
                                    int index) {
        char delimiter = ',';

        ifstream profilesFile (fileName.c_str());
        if (!profilesFile.is_open())  {
            cerr<<"Error: Unable to open file - "<<fileName.c_str()<<endl;
            return false;
        }

        // check header 1
        std::string headerLine1;
        getline (profilesFile, headerLine1);

        // check it is an image profiles file
        std::istringstream header1Stream( headerLine1 );
        std::string headerEntry;
        getline (header1Stream, headerEntry,delimiter); // ignore
        if(strcmp(headerEntry.c_str(), name.c_str())!=0) {
            cerr<<"Error: headerEntry isn't '"<<name.c_str()<<"' it is: "<<headerEntry<<endl;
            return false;
        }

        // check for the correct model type
        getline (header1Stream, headerEntry,delimiter);
        if(strcmp(headerEntry.c_str(), getModelName(index).c_str())!=0) {
            cerr<<"Miss-match between the model types specified in the '"<<name.c_str()<<"' header file ("<<headerEntry<<") and the 'Previously Defined'  ("<<getModelName(index)<<") model"<<endl;
            return false;
        }

        // check mheader line 2
        std::string headerLine2;
        getline (profilesFile, headerLine2);

        // get the number for parameters
        std::istringstream header2Stream( headerLine2 ); std::string importedNumberOfSamplesString;
        std::getline(header2Stream, importedNumberOfSamplesString, ':');
        std::getline(header2Stream, importedNumberOfSamplesString);
        int loadedNumberOfSamples = std::atoi(importedNumberOfSamplesString.c_str());
        if(loadedNumberOfSamples != numberOfSamples) {
            cerr << "Miss-match between the loaded (" << loadedNumberOfSamples << ") and local (" << numberOfSamples << ") 'Number Of Samples'" << endl;
            return false;
        }

        vtkIdType maxI = mesh->GetNumberOfPoints(), maxJ = numberOfSamples;
        array->SetNumberOfComponents(maxJ); array->SetNumberOfTuples(maxI);

        // read in the imported image profiles
        bool status = readArrayValues(profilesFile, array);

        profilesFile.close();

        return status;
    }

    bool CorticalBone::readTableValues(ifstream &file, vtkSmartPointer<vtkTable>& table) {

        char delimiter = ',';


        int maxI = table->GetNumberOfRows(), maxJ = table->GetNumberOfColumns();

        // populate the Parameters Table
        std::string line;
        for ( int i=0; i<maxI; i++ ) {

            getline(file, line);
            if(file.fail()) {
                cerr<<"Error ["<<i<<"]"<<line<<endl;
                table = NULL;
                return false;
            }

            std::istringstream dataStream(line);
            for( int j = 0; j<maxJ; j++ ) {
                std::string var;
                getline(dataStream, var, delimiter);

                if(dataStream.fail()) {
                    cerr<<"Error ["<<i<<"]["<<j<<"]"<<var<<endl;
                    table = NULL;
                    return false;
                }
                table->SetValue(i, j, std::atof(var.c_str()));
            }
        }
        if(getline(file, line)) {
            cerr<<"Error: too many lines in the selected parameters file"<<endl;
            table = NULL;
            return false;
        }
        
        return true;
    }

    bool CorticalBone::readArrayValues(ifstream &file, vtkSmartPointer<vtkDoubleArray>& array) {

        char delimiter = ',';

        vtkIdType maxI = array->GetNumberOfTuples(), maxJ = array->GetNumberOfComponents();


        // populate the Parameters Table
        std::string line;
        for ( int i=0; i<maxI; i++ ) {

            getline(file, line);
            if(file.fail()) {
                cerr<<"Error ["<<i<<"]"<<line<<endl;
                array = NULL;
                return false;
            }

            std::istringstream dataStream(line);

            for( int j = 0; j<maxJ; j++ ) {
                std::string var;
                getline(dataStream, var, delimiter);

                if(dataStream.fail()) {
                    cerr<<"Error ["<<i<<"]["<<j<<"]"<<var<<endl;
                    array = NULL;
                    return false;
                }
                array->SetComponent(i, j, std::atof(var.c_str()));
            }
        }
        if(getline(file, line)) {
            cerr<<"Error: too many lines in the selected parameters file"<<endl;
            array = NULL;
            return false;
        }

        return true;
    }

    bool CorticalBone::saveTable(std::string fileName, vtkSmartPointer<vtkTable> table, std::string description) {

        int iMax = table->GetNumberOfRows(), jMax = table->GetNumberOfColumns();

        // construct header information
        std::string delimiter = ",";
        std::string header_line1 = description;
        header_line1.append(delimiter).append(getModelName()).append(delimiter).append(getSchemeName()).append(delimiter).append(getOptimisationName()).append(delimiter).append(getPhantomName(phantomIndex)).append(delimiter).append(getCBDensityState());

        std::string header_line2;
        for(int j=0; j<jMax; j++) {
            // trim ""'s
            std::string columnName = table->GetColumnName(j);
            header_line2.append(columnName.c_str());
            if(j!=jMax-1) {
                header_line2.append(delimiter);
            }
        }

        ofstream tableFile;
        tableFile.open (fileName.c_str()); //cerr<<"Opened File: "<<fileName.c_str()<<endl;
        tableFile << header_line1 << endl;
        tableFile << header_line2 << endl;

        for(int i=0; i<iMax; i++) {
            std::ostringstream lineStream; lineStream.precision(9);
            for(int j=0; j<jMax; j++) {

                lineStream << table->GetValue(i, j);

                if(j!=jMax-1) {
                    lineStream << delimiter;
                }
            }

            tableFile << lineStream.str() << endl;
        }
        tableFile.close();
        
        return true;
    }

    bool CorticalBone::saveArray(std::string fileName, vtkSmartPointer<vtkDoubleArray> array,
                                    std::string header_line1, std::string header_line2) {

        std::string delimiter = ",";

        ofstream profilesStream;
        profilesStream.open (fileName.c_str()); cerr << "Saving File: " << fileName.c_str() << endl;
        profilesStream << header_line1 << endl;
        profilesStream << header_line2 << endl;

        for(int i=0; i<array->GetNumberOfTuples(); i++) {
            std::string line;
            for(int j=0; j<array->GetNumberOfComponents(); j++) {

                std::ostringstream imageValueStream;
                imageValueStream << std::setprecision(10) << array->GetComponent(i, j);

                line.append(imageValueStream.str()).append(delimiter);
            }

            profilesStream << line << endl;
        }
        profilesStream.close();

        return true;
    }

    bool CorticalBone::saveProfiles(std::string fileName, vtkSmartPointer<vtkDoubleArray> array,
                                    std::string description) {

        // construct header information
        std::string delimiter = ",";
        std::string header_line1 = description;
        header_line1.append(delimiter).append(getModelName()).append(delimiter).append(getSchemeName()).append(delimiter).append(getOptimisationName()).append(delimiter).append(getPhantomName(phantomIndex)).append(delimiter).append(getCBDensityState());

        std::ostringstream numberOfSamplesStream;
        numberOfSamplesStream << numberOfSamples;
        std::string header_line2 = "Number of Samples: " + numberOfSamplesStream.str();

        return saveArray(fileName, array, header_line1, header_line2);

    }

    /* Creation / Setup Methods */
    void CorticalBone::createParameters() {

        //---------- Rect Initial TF parameters ------------//
        int n = registrationMethod->GetNumberOfCombinedParameters(kThreeTierRect);
        // flexible parameters
        rectAllParameters.SetSize(n);
        rectAllParameters[0] = 0.0;     // X_P
        rectAllParameters[1] = profileLength * 1.0/3.0;     // X_E
        rectAllParameters[2] = 0;                           // Y_TB - 25
        rectAllParameters[3] = 1450; //1000;   // Y_CB - 800 or maxImageValue*0.8
        rectAllParameters[4] = 0;                           // Y_ST - -50
        rectAllParameters[5] = 0.75;                           // Sigma - 0.2

        // scale parameters
        rectAllScales.SetSize(n); // scale for max range -1 to 1
        rectAllScales[0] = 1.0/profileLength;       // Thickness
        rectAllScales[1] = 1.0/profileLength;       // Centre Position
        rectAllScales[2] = 1.0/maxBMDValue;              // Soft Tissue Density - once devided by maxImageValue
        rectAllScales[3] = 1.0/maxBMDValue;             // Cortical Bone Density
        rectAllScales[4] = 1.0/maxBMDValue;              // Trabeculae Bone Density
        rectAllScales[5] = 4.0/profileLength;           // Sigma

        // fixed parameters
        rectFlexParameters.SetSize(n);
        rectFlexScales.SetSize(n);
        for(int i=0; i<n; i++) {
            rectFlexParameters[i] = rectAllParameters[i];
            rectFlexScales[i] = rectAllScales[i];
        }

        //---------- Ramp Initial TF parameters ------------//
        // ramp scale parameters (all and fixed))
        rampAllScales = convertRectToRampScales(rectAllScales);         // add the Endosteal Transition width
        rampFlexScales = convertRectToRampScales(rectAllScales);

        // parameter mapping
        transformFlexMap.SetSize(registrationMethod->GetNumberOfCombinedParameters(kEndostealRamp)); // use TF with greatest number of parameters
        transformFlexMap.Fill(1.0); // all are flex. 0.0 indicates a parameter is fixed

    }

    void CorticalBone::createRegistrationMethods() {

        Utilities::measureTime(true);

        //---------- Registration (includes TFs and Optimisers) setup ----------------//
        createFittingMethod();

        //---------- Threshold Method Setup ------------//
        createClassifierMethod();

        //---------- Calibration Method Setup ----------//
        createCalibrationMethod();

        //----------- Display / Storage Setup ----------//
        clearPreviousResults();

        //---------- Update model profile settings -------//
        modelProfile->SetEdgeRatio(defaultPeriostealEdgeRatio); // call after the registration/threolding methods have been called to ensure later modification time
        modelProfile->TurnOnCurvedProfiles(mesh, meshTree, meshNormals); // call after mesh is created  todo decide if always on

        if(verbose){cout<<Utilities::getTabString()<<Utilities::getTabString()<<Utilities::measureTime(false)<<" to create the registration methods"<<endl;}
    }

    void CorticalBone::createModelProfile() { // call after image opened
        //-------------- Profile Setup ------------------------//
        modelProfile = ProfileType::New();
        modelProfile->SetInterpolator(interpolator);
        modelProfile->SetProfileLength(profileLength);
        modelProfile->SetNumberOfSamples(numberOfSamples);
        modelProfile->SetEdgeRatio(defaultPeriostealEdgeRatio);
        modelProfile->SetMaximumPeriostealOffset(maximumOffset); // as long as not nan/inf will be valid
        if(calibrationSet) {
            modelProfile->SetCalibration(calP0, calP1, calP2);
        } else {
            modelProfile->GetCalibration(calP0, calP1, calP2);
        }

        ImageType::SpacingType image_spacing = image->GetSpacing();
        if(verbose){cout<<Utilities::getTabString()<<Utilities::getTabString()<<"image spacing "<<image_spacing<<endl;}
        double imageResolution[Dimension] = {image_spacing[0], image_spacing[1], image_spacing[2]};
        modelProfile->SetResolution(imageResolution);
    }

    void CorticalBone::createFittingMethod() {
        registrationMethod = ModelRegistrationType::New();
        registrationMethod->SetProfile(modelProfile);
        registrationMethod->SetMaxBMDDensity(maxBMDValue); // all values intermal to the registration method assumed to be in BMD
        registrationMethod->SetOptimiserSelection(kLMOptimiser);


        //--------- Parameter Creation ---------------//
        createParameters();
    }

    void CorticalBone::createClassifierMethod() {
        classifierMethod = ThresholdRegistrationType::New();

        classifierMethod->SetProfile(modelProfile);
        if(!isnan(maximumOffset)) {
            classifierMethod->TurnOnProfileResampling(maximumOffset);
        } else {
            classifierMethod->TurnOffProfileResampling();
        }
    }

    void CorticalBone::createCalibrationMethod() {
        calibrationMethod = CalibrationType::New();

        calibrationMethod->SetProfile(modelProfile);
    }

    /*------------- Private Setters -------------------*/
    void CorticalBone::setImageBuffer(){
        return; // still unable to make it work.

        if(imageSet && meshSet) {// && bufferSet == false) {

            // Get the extents of the mesh (increase by the amount of the profile line outside the mesh)
            double profileLengthOutsideMesh = profileLength * defaultPeriostealEdgeRatio;
            double meshBounds[6]; mesh->GetBounds(meshBounds);

            ImageType::PointType startPoint;
            startPoint[0]=meshBounds[0] - profileLengthOutsideMesh;
            startPoint[1]=meshBounds[2] - profileLengthOutsideMesh;
            startPoint[2]=meshBounds[4] - profileLengthOutsideMesh;

            ImageType::PointType endPoint;
            endPoint[0]=meshBounds[1] + profileLengthOutsideMesh;
            endPoint[1]=meshBounds[3] + profileLengthOutsideMesh;
            endPoint[2]=meshBounds[5] + profileLengthOutsideMesh;

            // Get the image index of the mesh boundaries
            ImageType::IndexType startIndex, endIndex;
            bool startInside = image->TransformPhysicalPointToIndex(startPoint, startIndex);
            bool endInside = image->TransformPhysicalPointToIndex(endPoint, endIndex);


            ImageType::RegionType region = image->GetBufferedRegion();
            std::cerr<<"Old Region. Start Index: "<<region.GetIndex()<<", Size: "<<region.GetSize()<<std::endl;
            if(!startInside) {
                //cerr<<"start isn't inside"<<endl;
                ImageType::RegionType::IndexType startImageIndex = region.GetIndex();
                startIndex[0] = startImageIndex[0];
                startIndex[1] = startImageIndex[1];
                startIndex[2] = startImageIndex[2];
                //image->TransformIndexToPhysicalPoint(m_startIndex, startPoint);
            }
            if(!endInside) {
                //cerr<<"end isn't inside"<<endl;
                ImageType::RegionType::IndexType startImageIndex = region.GetIndex();
                ImageType::RegionType::SizeType imageSize = region.GetSize();
                endIndex[0] = startImageIndex[0] + imageSize[0];
                endIndex[1] = startImageIndex[1] + imageSize[1];
                endIndex[2] = startImageIndex[2] + imageSize[2];
            }

            // Set the image extents based upon the image size
            ImageType::RegionType::SizeType regionSize;
            regionSize[0] = endIndex[0] - startIndex[0];
            regionSize[1] = endIndex[1] - startIndex[1];
            regionSize[2] = endIndex[2] - startIndex[2];

            //cerr<<"Image regions reset to size of mesh"<<endl;
            region.SetSize(regionSize);
            region.SetIndex(startIndex);
            // Test case only add Zind=10
            //ImageType::RegionType::SizeType OldRegionSize = region.GetSize(); OldRegionSize[2] = OldRegionSize[2]-10;
            //ImageType::RegionType::IndexType OldRegionIndex = region.GetIndex(); OldRegionIndex[2] = OldRegionIndex[2]+10;
            //region.SetSize(OldRegionSize);
            //region.SetIndex(OldRegionIndex);

            //ImageType::SpacingType image_spacing = image->GetSpacing();
            //bufferOffset[0] = image_spacing[0] * m_startIndex[0];
            //bufferOffset[1] = image_spacing[1] * m_startIndex[1];
            //bufferOffset[2] = image_spacing[2] * m_startIndex[2];

            //return;
            image->SetRequestedRegion(region);
            image->SetBufferedRegion(region);
            image->Update();

        }

    }

    void CorticalBone::setSampling() {

        // set length based on mesh - actually use fixed values set in initialise
        //profileLength = mesh->GetLength()/10;
        //defaultPeriostealEdgeRatio = 1.0/3;

        // set sample number based on profile length and image sampling rate
        ImageType::SpacingType image_spacing = image->GetSpacing();

        // sample at nyquest of dimension with greatest nyquest requirement - sw compatability
        numberOfSamples = ceil(profileLength / (std::min(image_spacing[0], image_spacing[1])/2))+1;
        if (numberOfSamples < 75) { numberOfSamples = 75; } // sw compatability
        if(verbose){cout<<Utilities::getTabString()<<"Number of Samples: "<<numberOfSamples<<";"<<Utilities::getTabString()<<"Image spacing: ["<<image_spacing[0]<<","<<image_spacing[1]<<","<<image_spacing[2]<<"]"<<endl;}
    }

    ParametersType CorticalBone::convertRectToRampParameters(ParametersType rectParameters) {

        int rampDimension = rectParameters.GetSize()+1;
        ParametersType rampParameters(rampDimension);
        for(int i=0;i<rampDimension-1;i++) {
            rampParameters[i] = rectParameters[i];
        }
        rampParameters[rampDimension-1] = 0.5 * (rectParameters[1] - rectParameters[0]); // X_RW (ramp width)

        return rampParameters;
    }

    ParametersType CorticalBone::convertRectToRampScales(ParametersType rectScales) {

        int rampDimension = rectScales.GetSize()+1;
        ParametersType rampScales(rampDimension);
        for(int i=0;i<rampDimension-1;i++) {
            rampScales[i] = rectScales[i];
        }
        rampScales[rampDimension-1] = 1.0/profileLength; // X_RW (ramp width)

        return rampScales;
    }

    ParametersType CorticalBone::getTableRow(vtkSmartPointer<vtkTable> table, vtkIdType rowID) {
        int n = table->GetNumberOfColumns();
        ParametersType row(n);

        for(int i=0; i<n; i++) {
            row[i] = table->GetValue(rowID, i).ToDouble();
        }
        return row;

    }

    /* Processing Methods */
    void CorticalBone::runModelFittingPipeline(vtkIdType meshId, int fittingSchemeIndex, int selectedModel) {

        //----------- Check if the point is in image ---------------//
        if(!registrationMethod->IsPtInsideImage()){//isInsideImageBounds) {
        } else if(fittingSchemeIndex == kStdAFitting) {
            runV1aFittingPipeline();
        } else if(fittingSchemeIndex == kStdBFitting) {
            runV1bFittingPipeline();
        } else if(fittingSchemeIndex == kStdCFitting) {
            runV1cFittingPipeline();
        } else if(fittingSchemeIndex == kStdDFitting) {
            runV1dFittingPipeline();
        } else if(fittingSchemeIndex == kCBMV2AFitting) {
            runV2aFittingPipeline(); // todo - consider error check to ensure that the sigma parameter exists
        } else if(fittingSchemeIndex == kCBMV2BFitting) {
            runV2bFittingPipeline(); // todo - consider error check to ensure that the sigma parameter exists
        } else if(fittingSchemeIndex == kCBMV2CFitting) {
            runV2cFittingPipeline(); // todo - consider error check to ensure that the sigma parameter exists
        } else if(fittingSchemeIndex == kCBMV2DFitting) {
            runV2dFittingPipeline(); // todo - consider error check to ensure that the sigma parameter exists
        } else if(fittingSchemeIndex == kUnconstrainedAFitting) {
            runV3aFittingPipeline(selectedModel);
        } else if(fittingSchemeIndex == kUnconstrainedBFitting) {
            runV3bFittingPipeline(selectedModel);
        } else if(fittingSchemeIndex == kUnconstrainedCFitting) {
            runV3cFittingPipeline(selectedModel);
        } else {
            cerr<<"Error: fittingSchemeIndex out of bounds in CorticalBone::runModelFittingPipeline()"<<endl;
            throw std::logic_error( "Error: fittingSchemeIndex out of bounds in CorticalBone::runModelFittingPipeline()" );
        }

        // store values
        storeParameterValues(meshId);
        storeSampledProfiles(meshId);
        storeDisplayValues(meshId);

        // store validity of model
        if (modelFitStatusArray[meshId]==kMasked) {
            cerr<<"Error in: CorticalBone::runModelFIttingPipeline() should not have entered for masked point"<<endl;
        } else if(!registrationMethod->IsPtInsideImage()) {
            modelFitStatusArray[meshId] = kOutOfBounds;
            outOfBoundsCount++;
        } else if(registrationMethod->IsValid()) {
            modelFitStatusArray[meshId] = kValid;
        } else {
            modelFitStatusArray[meshId] = kInvalid;
            failureCount++;
        }
        if(registrationMethod->IsPtInsideImage() && isnan(registrationMethod->GetPeriostealEdgePosistion())) {
            nanCount++;
        }
    }

    void CorticalBone::updateModelProfile(vtkIdType meshId, bool runningOverMesh) { // returns true if location is inside the image
        //----------- Get measurement point and normal -------------//
        double normal[Dimension], location[Dimension];
        meshTree->GetDataSet()->GetPoint(meshId, location);
        meshNormals->GetTuple(meshId, normal);

        //location[0] = location[0] - bufferOffset[0]; // seemed to work for  an offset of Zind=10. Doesn't work in general
        //location[1] = location[1] - bufferOffset[1];
        //location[2] = location[2] - bufferOffset[2];

        if(!runningOverMesh) {
            if(verbose){std::cout<<Utilities::getTabString()<<"Selected Measurement Point"<<Utilities::getTabString()<<" x: "<<location[0]<<", y: "<<location[1]<<", z: "<<location[2]<<std::endl;}
        }

        // decide weather or not to turn on double peak detection
        bool doublePkDetection = false;
        if(meshMask[meshId]==kDblPeak && modelFitStatusArray[meshId]==kUndefined) {
            doublePkDetection = true;
        }

        //----------- Update Profile model -------------//
        modelProfile->SetEdgeAndOrientation(location, normal, doublePkDetection); // offset reset

    }

    void CorticalBone::fitModel(int localModelIndex, ParametersType parameters, int weightingFunction, double scale, bool dynamicWeighting) {

        registrationMethod->SetWeightingMode(weightingFunction, scale);
        registrationMethod->SetDynamicWeighting(dynamicWeighting);
        ParametersType scales = (localModelIndex==kThreeTierRect) ? rectFlexScales : rampFlexScales;
        registrationMethod->Initialize(parameters, scales, localModelIndex);

        try{
            registrationMethod->Update();
        } catch( itk::ExceptionObject & e ) {
            std::cerr << "An error occurred during Optimisation. Location    = " << e.GetLocation()<< "Description = " << e.GetDescription() << std::endl;
        }
    }

    void CorticalBone::runV1aFittingPipeline() {
        //---------- Update Fixed Constraints as necessary ----------//
        if(modelProfile->IsSigmaSet()) {
            setFixedSigma();
        }
        double scale = 1.0/1500.0;
        fitModel(kThreeTierRect, rectFlexParameters, TransformBaseType::kHueristicWeighting, scale);

        
        if(registrationMethod->IsValid()) {
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();
            if(modelIndex==kEndostealRamp) {
                optimisedParameters = convertRectToRampParameters(optimisedParameters);
            }
            scale = (modelIndex==kThreeTierRect) ? 1.0/50.0 : 1.0/100.0;
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kHueristicWeighting, scale); // todo replace -- , true);
        }
    }

    void CorticalBone::runV1bFittingPipeline() {
        //---------- Update Fixed Constraints as necessary ----------//
        if(modelProfile->IsSigmaSet()) {
            setFixedSigma();
        }
        double scale = 1.0/1500.0;
        fitModel(kThreeTierRect, rectFlexParameters, TransformBaseType::kWeightingFunction1, scale);


        if(registrationMethod->IsValid()) {
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();
            if(modelIndex==kEndostealRamp) {
                optimisedParameters = convertRectToRampParameters(optimisedParameters);
            }
            scale = (modelIndex==kThreeTierRect) ? 1.0/50.0 : 1.0/100.0;
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kWeightingFunction1, scale, true);
        }
    }

    void CorticalBone::runV1cFittingPipeline() {
        //---------- Update Fixed Constraints as necessary ----------//
        if(modelProfile->IsSigmaSet()) {
            setFixedSigma();
        }
        double scale = 1.0/1500.0;
        fitModel(kThreeTierRect, rectFlexParameters, TransformBaseType::kWeightingFunction1, scale);


        if(registrationMethod->IsValid()) {
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();
            if(modelIndex==kEndostealRamp) {
                optimisedParameters = convertRectToRampParameters(optimisedParameters);
            }
            scale = (modelIndex==kThreeTierRect) ? 1.0/50.0 : 1.0/100.0;
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kWeightingFunction2, scale, true); // was 25
        }
    }

    void CorticalBone::runV1dFittingPipeline() {
        //---------- Update Fixed Constraints as necessary ----------//
        if(modelProfile->IsSigmaSet()) {
            setFixedSigma();
        }
        double scale = 1.0/1500.0;
        fitModel(kThreeTierRect, rectFlexParameters, TransformBaseType::kWeightingFunction1, scale);


        if(registrationMethod->IsValid()) {
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();
            if(modelIndex==kEndostealRamp) {
                optimisedParameters = convertRectToRampParameters(optimisedParameters);
            }
            scale = (modelIndex==kThreeTierRect) ? 1.0/50.0 : 1.0/100.0;
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kDensityWeighting, scale, true); // kWeightingFunction2
        }
    }

    void CorticalBone::runV2aFittingPipeline() { // grahams cbm v2 - 3 pass CBM V2 with sigma correction
        // reset constraints
        resetCBConstraint(); // reset to global, or removed if no global set
        removeFixedSigma(false);

        // run two pass fit
        double scale = 1.0/1500.0;
        fitModel(kThreeTierRect, rectFlexParameters, TransformBaseType::kHueristicWeighting, scale);
        if(registrationMethod->IsValid()) {
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();
            if(modelIndex==kEndostealRamp) {
                optimisedParameters = convertRectToRampParameters(optimisedParameters);
            }
            scale = (modelIndex==kThreeTierRect) ? 1.0/100.0 : 1.0/100.0;
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kHueristicWeighting, scale, true);
        }

        // perform sigma based CB estimate correction
        if(registrationMethod->IsValid()) {

            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters(); // get optimiser parameters

            // calculate new cb estimate
            double sigmaAdjustedDensity = registrationMethod->CalculateCBMv2SigmaAdjustedCBDensity(modelIndex); // sigma correction applied to selected model


            if(!isCBFixed()) { // is cb not fixed. fix cb and remove from optimised params
                int cbFlexIndex=getFlexIndex(ModelRegistrationType::kYcbParam);
                removeParametersValue(optimisedParameters, cbFlexIndex); // remove the cb value from the optimised values
            }

            scale = (modelIndex==kThreeTierRect) ? 1.0/25.0 : 1.0/50.0;
            setFixedCBDensity(sigmaAdjustedDensity, false); // set the sigma adjusted estaimte of cortical density as fixed
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kHueristicWeighting, scale, true);
        }
    }

    void CorticalBone::runV2bFittingPipeline() { // grahams cbm v2 - rect to ramp, ramp CBM V2 correction and selected - SELECTED METHOD
        // reset constraints
        resetCBConstraint(); // reset to global, or removed if no global set
        removeFixedSigma(false);

        // run two pass fit
        double scale = 1.0/1500.0;
        fitModel(kThreeTierRect, rectFlexParameters, TransformBaseType::kWeightingFunction1, scale);
        if(registrationMethod->IsValid()) {
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();
            if(modelIndex==kEndostealRamp) {
                optimisedParameters = convertRectToRampParameters(optimisedParameters);
            }
            scale = (modelIndex==kThreeTierRect) ? 1.0/100.0 : 1.0/100.0;
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kWeightingFunction1, scale, true);
        }

        // perform sigma based CB estimate correction
        if(registrationMethod->IsValid()) {

            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters(); // get optimiser parameters

            // calculate new cb estimate
            double sigmaAdjustedDensity = registrationMethod->CalculateCBMv2SigmaAdjustedCBDensity(modelIndex); // sigma correction applied to selected model


            if(!isCBFixed()) { // is cb not fixed. fix cb and remove from optimised params
                int cbFlexIndex=getFlexIndex(ModelRegistrationType::kYcbParam);
                removeParametersValue(optimisedParameters, cbFlexIndex); // remove the cb value from the optimised values
            }
            scale = (modelIndex==kThreeTierRect) ? 1.0/25.0 : 1.0/50.0;
            setFixedCBDensity(sigmaAdjustedDensity, false); // set the sigma adjusted estaimte of cortical density as fixed
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kWeightingFunction1, scale, true);
        }
    }

    void CorticalBone::runV2cFittingPipeline() {
        // reset constraints
        resetCBConstraint(); // reset to global, or removed if no global set
        removeFixedSigma(false);

        // run two pass fit
        double scale = 1.0/1500.0;
        fitModel(kThreeTierRect, rectFlexParameters, TransformBaseType::kWeightingFunction1, scale);
        if(registrationMethod->IsValid()) {
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();
            if(modelIndex==kEndostealRamp) {
                optimisedParameters = convertRectToRampParameters(optimisedParameters);
            }
            scale = (modelIndex==kThreeTierRect) ? 1.0/100.0 : 1.0/100.0;
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kWeightingFunction2, scale, true);
        }

        // perform sigma based CB estimate correction
        if(registrationMethod->IsValid()) {

            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters(); // get optimiser parameters

            // calculate new cb estimate
            double sigmaAdjustedDensity = registrationMethod->CalculateCBMv2SigmaAdjustedCBDensity(modelIndex); // sigma correction applied to selected model


            if(!isCBFixed()) { // is cb not fixed. fix cb and remove from optimised params
                int cbFlexIndex=getFlexIndex(ModelRegistrationType::kYcbParam);
                removeParametersValue(optimisedParameters, cbFlexIndex); // remove the cb value from the optimised values
            }

            scale = (modelIndex==kThreeTierRect) ? 1.0/25.0 : 1.0/50.0;
            setFixedCBDensity(sigmaAdjustedDensity, false); // set the sigma adjusted estaimte of cortical density as fixed
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kWeightingFunction2, scale, true);
        }
    }

    void CorticalBone::runV2dFittingPipeline() {
        // reset constraints
        resetCBConstraint(); // reset to global, or removed if no global set
        removeFixedSigma(false);

        // run two pass fit
        double scale = 1.0/1500.0;
        fitModel(kThreeTierRect, rectFlexParameters, TransformBaseType::kWeightingFunction1, scale);
        if(registrationMethod->IsValid()) {
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();
            if(modelIndex==kEndostealRamp) {
                optimisedParameters = convertRectToRampParameters(optimisedParameters);
            }
            scale = (modelIndex==kThreeTierRect) ? 1.0/100.0 : 1.0/100.0;
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kWeightingFunction2, scale, true);
        }

        // perform sigma based CB estimate correction
        if(registrationMethod->IsValid()) {

            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters(); // get optimiser parameters

            // calculate new cb estimate
            double sigmaAdjustedDensity = registrationMethod->CalculateCBMv2SigmaAdjustedCBDensity(modelIndex); // sigma correction applied to selected model


            if(!isCBFixed()) { // is cb not fixed. fix cb and remove from optimised params
                int cbFlexIndex=getFlexIndex(ModelRegistrationType::kYcbParam);
                removeParametersValue(optimisedParameters, cbFlexIndex); // remove the cb value from the optimised values
            }

            scale = (modelIndex==kThreeTierRect) ? 1.0/25.0 : 1.0/50.0;
            setFixedCBDensity(sigmaAdjustedDensity, false); // set the sigma adjusted estaimte of cortical density as fixed
            fitModel(modelIndex, optimisedParameters, TransformBaseType::kDensityWeighting, scale, true);
        }
    }

    void CorticalBone::runV3aFittingPipeline(int selectedModel) { // run constrained rect then unconstrained (selected model)
        //---------- Update Fixed Constraints as necessary ----------//
        if(modelProfile->IsSigmaSet()) {
            setFixedSigma();
        }
        resetCBConstraint();

        double scale = 1.0/20000.0; // todo replace -- 1.0/1500.0;
        fitModel(kThreeTierRect, rectFlexParameters, TransformBaseType::kWeightingFunction1, scale);

        if(registrationMethod->IsValid()){
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();
            if(selectedModel==kEndostealRamp) {
                optimisedParameters = convertRectToRampParameters(optimisedParameters);
            }
            //scale = (modelIndex==kThreeTierRect) ? 1.0/100.0 : 1.0/100.0; // see if smaller scales effects the density estimations
            fitModel(selectedModel, optimisedParameters, TransformBaseType::kWeightingFunction1, scale); // todo replace -- true);
        }

        if(registrationMethod->IsValid()){ // remove constraints
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();

            // temporarily remove constraints
            if(registrationMethod->isSigmaFixed()){
                double sigma = registrationMethod->GetSigma();
                removeFixedSigma(false);

                int sigmaFlexIndex = getFlexIndex(ModelRegistrationType::kSigmaParam);
                insertParametersValue(optimisedParameters, sigma, sigmaFlexIndex); // remove the cb value from the optimised values

            }
            if(registrationMethod->isCBDensityFixed()){
                double cbDensity = registrationMethod->GetCorticalDensity();
                removeFixedCBDensity(false);

                int cbFlexIndex = getFlexIndex(ModelRegistrationType::kYcbParam);
                insertParametersValue(optimisedParameters, cbDensity, cbFlexIndex); // remove the cb value from the optimised values

            }

            //scale = (modelIndex==kThreeTierRect) ? 1.0/25.0 : 1.0/50.0; // see if smaller scales effects the density estimations
            fitModel(selectedModel, optimisedParameters, TransformBaseType::kWeightingFunction1, scale); // todo replace -- true);
        }

    }

    void CorticalBone::runV3bFittingPipeline(int selectedModel) { // run constained (selected model) then unconstrained (selected model)
        //---------- Update Fixed Constraints as necessary ----------//
        if(modelProfile->IsSigmaSet()) {
            setFixedSigma();
        }
        resetCBConstraint();

        double scale = 1.0/1500.0;
        fitModel(kThreeTierRect, rectFlexParameters, TransformBaseType::kWeightingFunction1, scale);

        if(registrationMethod->IsValid()){
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();
            if(selectedModel==kEndostealRamp) {
                optimisedParameters = convertRectToRampParameters(optimisedParameters);
            }
            scale = (modelIndex==kThreeTierRect) ? 1.0/100.0 : 1.0/100.0;
            fitModel(selectedModel, optimisedParameters, TransformBaseType::kWeightingFunction2, scale, true);
        }

        if(registrationMethod->IsValid()){ // remove constraints
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();

            // temporarily remove constraints
            if(registrationMethod->isSigmaFixed()){
                double sigma = registrationMethod->GetSigma();
                removeFixedSigma(false);

                int sigmaFlexIndex = getFlexIndex(ModelRegistrationType::kSigmaParam);
                insertParametersValue(optimisedParameters, sigma, sigmaFlexIndex); // remove the cb value from the optimised values

            }
            if(registrationMethod->isCBDensityFixed()){
                double cbDensity = registrationMethod->GetCorticalDensity();
                removeFixedCBDensity(false);

                int cbFlexIndex = getFlexIndex(ModelRegistrationType::kYcbParam);
                insertParametersValue(optimisedParameters, cbDensity, cbFlexIndex); // remove the cb value from the optimised values

            }

            scale = (modelIndex==kThreeTierRect) ? 1.0/25.0 : 1.0/50.0;
            fitModel(selectedModel, optimisedParameters, TransformBaseType::kWeightingFunction2, scale, true);
        }
    }

    void CorticalBone::runV3cFittingPipeline(int selectedModel) { // run constained (selected model) then unconstrained (selected model)
        //---------- Update Fixed Constraints as necessary ----------//
        if(modelProfile->IsSigmaSet()) {
            setFixedSigma();
        }
        resetCBConstraint();

        double scale = 1.0/1500.0;
        fitModel(kThreeTierRect, rectFlexParameters, TransformBaseType::kWeightingFunction1, scale);

        if(registrationMethod->IsValid()){
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();
            if(selectedModel==kEndostealRamp) {
                optimisedParameters = convertRectToRampParameters(optimisedParameters);
            }
            scale = (modelIndex==kThreeTierRect) ? 1.0/100.0 : 1.0/100.0;
            fitModel(selectedModel, optimisedParameters, TransformBaseType::kWeightingFunction2, scale, true);
        }

        if(registrationMethod->IsValid()){ // remove constraints
            ParametersType optimisedParameters = registrationMethod->GetLastTransformParameters();

            // temporarily remove constraints
            if(registrationMethod->isSigmaFixed()){
                double sigma = registrationMethod->GetSigma();
                removeFixedSigma(false);

                int sigmaFlexIndex = getFlexIndex(ModelRegistrationType::kSigmaParam);
                insertParametersValue(optimisedParameters, sigma, sigmaFlexIndex); // remove the cb value from the optimised values

            }
            if(registrationMethod->isCBDensityFixed()){
                double cbDensity = registrationMethod->GetCorticalDensity();
                removeFixedCBDensity(false);

                int cbFlexIndex = getFlexIndex(ModelRegistrationType::kYcbParam);
                insertParametersValue(optimisedParameters, cbDensity, cbFlexIndex); // remove the cb value from the optimised values

            }

            scale = (modelIndex==kThreeTierRect) ? 1.0/25.0 : 1.0/50.0;
            fitModel(selectedModel, optimisedParameters, TransformBaseType::kDensityWeighting, scale, true);
        }
    }

    void CorticalBone::resetCBConstraint() {

        // if fixed CB - reset
        if(!isnan(globalFixedCB)) {
            setFixedCBDensity(globalFixedCB, false);
        } else {
            removeFixedCBDensity(false);
        }

    }

    void CorticalBone::turnOnRegistrationSigmaCorrection() {
        double sigma = modelProfile->GetProfileSigma();

        registrationMethod->TurnOnSigmaCorrection(sigma);
    }

    void CorticalBone::turnOffRegistrationSigmaCorrection() {

        registrationMethod->TurnOffSigmaCorrection();
    }

    void CorticalBone::runThresholdingPipeline(vtkIdType meshId) {

        classifierMethod->Process();

        storeParameterValues(meshId);
        storeSampledProfiles(meshId);
        storeDisplayValues(meshId);

        if(modelFitStatusArray[meshId]==kMasked) {
            cerr<<"Error in CorticalBone::runThresholdingPipeline() should not have entered with a masked point"<<endl;
        } else if(!classifierMethod->AreProfilesSufficientlyInsideImage()) {
            modelFitStatusArray[meshId] = kOutOfBounds;
            outOfBoundsCount++;
        } else if(classifierMethod->IsValid()) {
            modelFitStatusArray[meshId] = kValid;
        } else {
            modelFitStatusArray[meshId] = kInvalid;
            failureCount++;
        }
        if(classifierMethod->AreProfilesSufficientlyInsideImage() && isnan(classifierMethod->IsNan())) {
            nanCount++;
        }
    }

    void CorticalBone::runGenerateCalibrationValues(vtkIdType meshId) {
        calibrationMethod->Process();

        storeParameterValues(meshId);
        storeSampledProfiles(meshId);
        storeDisplayValues(meshId);

        if(modelFitStatusArray[meshId]==kMasked) {
            cerr<<"Error in CorticalBone::runGenerateCalibrationValues() should not be entered for a masked point"<<endl;
        } else if(!calibrationMethod->AreProfilesSufficientlyInsideImage()) {
            modelFitStatusArray[meshId] = kOutOfBounds;
            outOfBoundsCount++;
        } else if(calibrationMethod->IsValid()) {
            modelFitStatusArray[meshId] = kValid;
        } else {
            modelFitStatusArray[meshId] = kInvalid;
            failureCount++;
        }
    }

    // fixed sigma constraint
    void CorticalBone::setFixedSigma() {

        if(!modelProfile->IsSigmaSet()) {
            return;
        }

        double sigma = modelProfile->GetProfileSigma();

        registrationMethod->FixParameter(sigma, ModelRegistrationType::kSigmaParam);

        // add value to 'flexParameterArray' - unless already added
        if(transformFlexMap[ModelRegistrationType::kSigmaParam]==1.0) {

            // Note: rect and ramp can be the same as ModelRegistrationType::kCorticalDensityParam index is less than divergance

            int flexIndex=getFlexIndex(ModelRegistrationType::kSigmaParam); // get index to remove

            transformFlexMap[ModelRegistrationType::kSigmaParam]=0.0;

            // remove value
            removeParametersValue(rectFlexParameters, flexIndex);
            removeParametersValue(rectFlexScales, flexIndex);
            removeParametersValue(rampFlexScales, flexIndex);

        }



    }

    void CorticalBone::removeFixedSigma(bool clearDisplayBool) {

        // set map value 0, remove fixed value, insert flex value
        registrationMethod->FreeParameter(ModelRegistrationType::kSigmaParam);

        // reinsert value into 'flexParameterArray' - unless not removed
        if(transformFlexMap[ModelRegistrationType::kSigmaParam]==0.0) {

            // Note: rect and ramp can be the same as ModelRegistrationType::kCorticalDensityParam index is less than divergance

            int flexIndex=getFlexIndex(ModelRegistrationType::kSigmaParam); // get index to remove

            transformFlexMap[ModelRegistrationType::kSigmaParam]=1.0;

            // insert value
            insertParametersValue(rectFlexParameters, rectAllParameters[ModelRegistrationType::kSigmaParam], flexIndex);
            insertParametersValue(rectFlexScales, rectAllScales[ModelRegistrationType::kSigmaParam], flexIndex);
            insertParametersValue(rampFlexScales, rampAllScales[ModelRegistrationType::kSigmaParam], flexIndex);

        }

        if(clearDisplayBool) {
            clearPreviousResults();
        }


    }

    bool CorticalBone::noPhantomCalibration() {
        calP1 = 1.0; calP0 = 0.0; calP2 = 0.0;
        phantomChanged = false; cntrlPts = NULL; calHUVals = calBMDVals = NULL;
        return calibrationSet;
    }

    bool CorticalBone::mindwaysCalibration() {

        //----------- Instructions: --------------//
        // Four options. Either
        // pt0 = bright top. pt1 = dim top, pt2 = bright bottom, pt3 = dim bottom
        // pt0 = bright bottom. pt1 = dim bottom, pt2 = bright top, pt3 = dim top
        // pt0 = dim top. pt1 = bright top, pt2 = dim bottom, pt3 = bright bottom
        // pt0 = dim bottom. pt1 = bright bottom, pt2 = dim top, pt3 = bright top
        //
        // calibration conversion:
        // HU = calSlope x BMD + calIntercept
        // BMD = (HU - calIntercept) / calSlope
        //
        // Stradwin Code: config_dialogs.cpp,
        //                line 3443
        //                void MainFrame::OnDICOMCalibrate(wxCommandEvent &event)


        // check right phantom index
        if(phantomIndex != kMindwaySolidCal) {
            cout<<"Error in CorticalBone::mindwaysCalibration. Incorrect phantom index selected "<<phantomIndex<<endl;
            return calibrationSet;
        }

        // work out bright / dim / top / bottom points
        double tb[Dimension], td[Dimension], bb[Dimension], bd[Dimension];
        if(!checkCalibrationControlPoints(tb,td,bb, bd)) {
            if(verbose){std::cout<<"Invalid Calibration (Mindways Solid K2HPO4): Landmarks are not correctly aligned."<<std::endl;}
            return false;
        }

        // slice info
        double zSpacing = image->GetSpacing()[kZ];
        int numberOfSlices = (int) (((tb[kZ] - bb[kZ]) / zSpacing ) + 1);


        //---- calculate other (in-between) core locations ---//
        // x & y offset vectors
        double xShift[Dimension], yShift[Dimension], dCoreA[Dimension], dCoreB[Dimension];
        vtkMath::Subtract(bd, bb, xShift); // same slice to z=0
        yShift[kX]=-1.0/148.0*xShift[kY];
        yShift[kY]= 1.0/148.0*xShift[kX];
        yShift[kZ]= 1.0/148.0*xShift[kZ];
        // bottom to top direction
        vtkMath::Subtract(tb, bb, dCoreA);  vtkMath::Normalize(dCoreA);
        vtkMath::Subtract(td, bd, dCoreB);  vtkMath::Normalize(dCoreB);

        // Locations of each core bottom
        double coreXYZ[calNumber][Dimension], coredXYZ[calNumber][Dimension];
        coreXYZ[0][kX]=bb[kX]; coreXYZ[0][kY]=bb[kY]; coreXYZ[0][kZ]=bb[kZ];
        coreXYZ[4][kX]=bd[kX]; coreXYZ[4][kY]=bd[kY]; coreXYZ[4][kZ]=bd[kZ];

        coreXYZ[1][kX]=bb[kX] + 0.25*xShift[kX]+3*yShift[kX];  coreXYZ[1][kY]=bb[kY] + 0.25*xShift[kY]+3*yShift[kY];  coreXYZ[1][kZ]=bb[kZ] + 0.25*xShift[kZ]+3*yShift[kZ];
        coreXYZ[2][kX]=bb[kX] + 0.5*xShift[kX]+5*yShift[kX];   coreXYZ[2][kY]=bb[kY] + 0.5*xShift[kY]+5*yShift[kY];   coreXYZ[2][kZ]=bb[kZ] + 0.5*xShift[kZ]+5*yShift[kZ];
        coreXYZ[3][kX]=bb[kX] + 0.75*xShift[kX]+3*yShift[kX];  coreXYZ[3][kY]=bb[kY] + 0.75*xShift[kY]+3*yShift[kY];  coreXYZ[3][kZ]=bb[kZ] + 0.75*xShift[kZ]+3*yShift[kZ];

        // Calculate the direction to travel along each core
        vtkMath::MultiplyScalar(dCoreA, zSpacing / dCoreA[kZ]); vtkMath::MultiplyScalar(dCoreB, zSpacing / dCoreB[kZ]);

        coredXYZ[0][kX]= dCoreA[kX]; coredXYZ[0][kY]= dCoreA[kY]; coredXYZ[0][kZ]= dCoreA[kZ];
        coredXYZ[4][kX]= dCoreB[kX]; coredXYZ[4][kY]= dCoreB[kY]; coredXYZ[4][kZ]= dCoreB[kZ];

        coredXYZ[1][kX]= 0.75 * dCoreA[kX] + 0.25 * dCoreB[kX]; coredXYZ[1][kY]= 0.75 * dCoreA[kY] + 0.25 * dCoreB[kY];
        coredXYZ[1][kZ]= 0.75 * dCoreA[kZ] + 0.25 * dCoreB[kZ];
        coredXYZ[2][kX]= 0.5 * dCoreA[kX] + 0.5 * dCoreB[kX]; coredXYZ[2][kY]= 0.5 * dCoreA[kY] + 0.5 * dCoreB[kY];
        coredXYZ[2][kZ]= 0.5 * dCoreA[kZ] + 0.5 * dCoreB[kZ];
        coredXYZ[3][kX]= 0.25 * dCoreA[kX] + 0.75 * dCoreB[kX]; coredXYZ[3][kY]= 0.25 * dCoreA[kY] + 0.75 * dCoreB[kY];
        coredXYZ[3][kZ]= 0.25 * dCoreA[kZ] + 0.75 * dCoreB[kZ];


        //---- Sample cores to generate hu / bmd pairs

        // work out number of samples for each pt
        ImageType::SpacingType imageSpacing = image->GetSpacing();
        double sampleSquareWidth= cntrlRadius;
        unsigned long nii = (unsigned long) (2 * floor((sampleSquareWidth / imageSpacing[0]) / 2.0) + 1); // number of x samples taken at each cal pt
        unsigned long nij = (unsigned long) (2 * floor((sampleSquareWidth / imageSpacing[1]) / 2.0) + 1); // number of y samples taken at each cal pt
        unsigned long nati = (nii * nij * numberOfSlices); // ensure is odd - n samples at each cal pt

        // set up HU and BMD arrays
        itk::Array<double> huArray = itk::Array<double>(calNumber * nati);
        itk::Array<double> bmdArray = itk::Array<double>(calNumber * nati);

        // calibration corrections
        double bmdCal[calNumber] = {375.8, 157.0, 58.9, -53.4, -51.8};
        double huCal[calNumber] = {923.2, 1119.5, 1103.6, 1057.0, 1012.2};

        // look up HU values
        for(int i=0; i< calNumber; i++) {


            for(int ik=0; ik< numberOfSlices; ik++) {

                // get point: ofset to courner and offset along profile
                double pt[Dimension];
                pt[kX] = coreXYZ[i][kX]+ik*coredXYZ[i][kX]-imageSpacing[0]*(int)(nii /2.0);
                pt[kY] = coreXYZ[i][kY]+ik*coredXYZ[i][kY]-imageSpacing[1]*(int)(nij /2.0);
                pt[kZ] = coreXYZ[i][kZ]+ik*coredXYZ[i][kZ];

                for (int ii = 0; ii < nii; ii++) {
                    for (int ij = 0; ij < nij; ij++) {


                        double huValue = interpolator->Evaluate(pt);
                        huArray[i * nati + ik * nij * nii + ii * nij + ij] = huValue-huCal[i];
                        bmdArray[i * nati + ik * nij * nii + ii * nij + ij] =bmdCal[i];

                        // increment point
                        pt[1] += imageSpacing[1];
                    }
                    pt[0] += imageSpacing[0];
                    pt[1] -= nij * imageSpacing[1]; // reset to starting y values
                }

            }
        }

        // calculate Least Squares Fit- set BMD as X as water corrections to apply in this direction then switch to BMD in Y
        LinearRegressionCalculator::RegressionLine regressionLine = LinearRegressionCalculator::CalculateLinearRegression(bmdArray, huArray);
        double toHUSlope = regressionLine.slope - 0.2174, toHUIntercept = regressionLine.intercept + 999.6;

        // convert to to BMD conversion
        calP0 = -toHUIntercept/toHUSlope; calP1 = 1/toHUSlope; calP2=0.0;

        if(!isnan(calP1) && !isnan(calP0) && !isinf(calP1) && !isinf(calP2)) {
            if(verbose){std::cout<<"Completed Calibration (Mindways Solid K2HPO4) [HU = m x BMD + b]: m = "<< calP1 <<", b = "<<
                        calP0 <<std::endl;}
            calBMDVals=Utilities::convertItkToVtkArray(bmdArray, "BMD values");
            calHUVals=Utilities::convertItkToVtkArray(huArray, "HU Values");
            calibrationSet=true;
        } else {
            if(verbose){std::cout<<"Unable to Complete Calibration (Mindways Solid K2HPO4) [HU = m x BMD + b]: m = "<<
                        calP1 <<", b = "<< calP0 <<std::endl;}
            calP1 = 1.0; calP0 = calP2 = 0.0; cntrlPts = NULL; calHUVals = calBMDVals = cntrlVals= NULL;
        }
        phantomChanged = false;
        return calibrationSet;
    }

    bool CorticalBone::BoneDensityCalibration() {

        // check right phantom index
        if(phantomIndex != kBoneDensityCal) {
            cout<<"Error in CorticalBone::mindwaysCalibration. Incorrect phantom index selected "<<phantomIndex<<endl;
            return calibrationSet;
        }

        // work out bright / dim / top / bottom points
        double tb[Dimension], td[Dimension], bb[Dimension], bd[Dimension];
        if(!checkCalibrationControlPoints(tb,td,bb, bd)) {
            if(verbose){std::cout<<"Invalid Calibration (QRM-BDM): Landmarks are not correctly aligned."<<std::endl;}
            return false;
        }



        // calibration conversion
        // HU = calSlope x BMD + calIntercept
        // BMD = (HU - calIntercept) / calSlope
        calP1 = 1.0; calP0 = 0.0; calP2=0.0; // reset values to no calibration values

        // Stradwin Code: config_dialogs.cpp,
        //                line 3443
        //                void MainFrame::OnDICOMCalibrate(wxCommandEvent &event)


        // slice info
        double zSpacing = image->GetSpacing()[kZ];
        int numberOfSlices = (int) (((tb[kZ] - bb[kZ]) / zSpacing ) + 1);


        //---- calculate other (in-between) core locations ---//
        // x & y offset vectors
        double xShift[Dimension], yShift[Dimension], dCoreA[Dimension], dCoreB[Dimension];
        vtkMath::Subtract(bd, bb, xShift); // same slice to z=0
        yShift[kX]=-1.0/148.0*xShift[kY];
        yShift[kY]= 1.0/148.0*xShift[kX];
        yShift[kZ]= 1.0/148.0*xShift[kZ];
        // bottom to top direction
        vtkMath::Subtract(tb, bb, dCoreA);  vtkMath::Normalize(dCoreA);
        vtkMath::Subtract(td, bd, dCoreB);  vtkMath::Normalize(dCoreB);

        // Locations of each core bottom
        double coreXYZ[calNumber][Dimension], coredXYZ[calNumber][Dimension];
        coreXYZ[0][kX]=bb[kX]; coreXYZ[0][kY]=bb[kY]; coreXYZ[0][kZ]=bb[kZ];
        coreXYZ[2][kX]=bd[kX]; coreXYZ[2][kY]=bd[kY]; coreXYZ[2][kZ]=bd[kZ];

        coreXYZ[1][kX]=bb[kX] + 0.5*xShift[kX]+yShift[kX];
        coreXYZ[1][kY]=bb[kY] + 0.5*xShift[kY]+yShift[kY];
        coreXYZ[1][kZ]=bb[kZ];

        // Calculate the direction to travel along each core
        vtkMath::MultiplyScalar(dCoreA, zSpacing / dCoreA[kZ]); vtkMath::MultiplyScalar(dCoreB, zSpacing / dCoreB[kZ]);

        coredXYZ[0][kX]= dCoreA[kX]; coredXYZ[0][kY]= dCoreA[kY]; coredXYZ[0][kZ]= dCoreA[kZ];
        coredXYZ[2][kX]= dCoreB[kX]; coredXYZ[2][kY]= dCoreB[kY]; coredXYZ[2][kZ]= dCoreB[kZ];

        coredXYZ[1][kX]= 0.5*dCoreA[kX]+0.5*dCoreB[kX];
        coredXYZ[1][kY]= 0.5*dCoreA[kY]+0.5*dCoreB[kY];
        coredXYZ[1][kZ]= 0.5*dCoreA[kZ]+0.5*dCoreB[kZ];


        //---- Sample cores to generate hu / bmd pairs

        // work out number of samples for each pt
        ImageType::SpacingType imageSpacing = image->GetSpacing();
        double sampleSquareWidth= cntrlRadius;
        unsigned long nii = (unsigned long) (2 * floor((sampleSquareWidth / imageSpacing[0]) / 2.0) + 1); // number of x samples taken at each cal pt
        unsigned long nij = (unsigned long) (2 * floor((sampleSquareWidth / imageSpacing[1]) / 2.0) + 1); // number of y samples taken at each cal pt
        unsigned long nati = (nii * nij * numberOfSlices); // ensure is odd - n samples at each cal pt

        // set up HU and BMD arrays
        itk::Array<double> huArray = itk::Array<double>(calNumber * nati);
        itk::Array<double> bmdArray = itk::Array<double>(calNumber * nati);

        double bmdCal[calNumber] = {200.0, 0.0, 100.0};

        // look up HU values
        for(int i=0; i< calNumber; i++) {


            for(int ik=0; ik< numberOfSlices; ik++) {

                // get point: ofset to courner and offset along profile
                double pt[Dimension];
                pt[kX] = coreXYZ[i][kX]+ik*coredXYZ[i][kX]-imageSpacing[0]*(int)(nii /2.0);
                pt[kY] = coreXYZ[i][kY]+ik*coredXYZ[i][kY]-imageSpacing[1]*(int)(nij /2.0);
                pt[kZ] = coreXYZ[i][kZ]+ik*coredXYZ[i][kZ];

                for (int ii = 0; ii < nii; ii++) {
                    for (int ij = 0; ij < nij; ij++) {


                        double huValue = interpolator->Evaluate(pt);
                        huArray[i * nati + ik * nij * nii + ii * nij + ij] = huValue;
                        bmdArray[i * nati + ik * nij * nii + ii * nij + ij] =bmdCal[i];

                        // increment point
                        pt[1] += imageSpacing[1];
                    }
                    pt[0] += imageSpacing[0];
                    pt[1] -= nij * imageSpacing[1]; // reset to starting y values
                }

            }
        }

        // calculate Least Squares Fit
        LinearRegressionCalculator::RegressionLine regressionLine = LinearRegressionCalculator::CalculateLinearRegression(huArray, bmdArray);
        calP0 = regressionLine.intercept; calP1 = regressionLine.slope; calP2=0.0;

        if(!isnan(calP1) && !isnan(calP0) && !isinf(calP1) && !isinf(calP2)) {
            if(verbose){std::cout<<"Completed Calibration (QRM-BDM) [HU = m x BMD + b]: m = "<< calP1 <<", b = "<<
                        calP0 <<std::endl;}
            calBMDVals=Utilities::convertItkToVtkArray(bmdArray, "BMD values");
            calHUVals=Utilities::convertItkToVtkArray(huArray, "HU Values");
            calibrationSet=true;
        } else {
            if(verbose){std::cout<<"Unable to Complete Calibration (QRM-BDM) [HU = m x BMD + b]: m = "<<
                        calP1 <<", b = "<< calP0 <<std::endl;}
            calP1 = 1.0; calP0 = calP2 = 0.0; cntrlPts = NULL; calHUVals = calBMDVals = cntrlVals= NULL;
        }
        phantomChanged = false;
        return calibrationSet;
    }

    bool CorticalBone::ESPCalibration() {


        //----------- Instructions: --------------//
        // L1 = low density. pt1 = cortex, pt2 = processes, pt3 = sponge
        // L2 = mid density. pt4 = cortex, pt5 = processes, pt6 = sponge
        // L3 = high density. pt7 = cortex, pt8 = processes, pt9 = sponge


        calibrationSet = false; // set QRM-EPS

        if(cntrlPts->GetNumberOfTuples() != cntrlNumber || cntrlNumber!=calNumber) {
            cout<<"Error in CorticalBone::ESPCalibration. Incorrect number of caibration points provided or contorl numbers don't match teh calibration numbers"<<endl;
            calHUVals = calBMDVals = cntrlVals= NULL; return calibrationSet;
        }

        // work out number of samples for each pt
        ImageType::SpacingType imageSpacing = image->GetSpacing();
        double sampleSquareWidth= cntrlRadius;
        unsigned long nii = (unsigned long)2*floor((sampleSquareWidth/ imageSpacing[0])/2.0)+1; // number of x samples taken at each cal pt
        unsigned long nij = (unsigned long)2*floor((sampleSquareWidth/ imageSpacing[1])/2.0)+1; // number of y samples taken at each cal pt
        unsigned long nati = nii * nij; // ensure is odd - n samples at each cal pt

        // set up HU and BMD arrays
        itk::Array<double> huArray = itk::Array<double>(cntrlNumber * nati);
        itk::Array<double> bmdArray = itk::Array<double>(cntrlNumber * nati);

        // look up HU values
        for(int i=0; i< cntrlNumber; i++) {
            double pt[Dimension]; cntrlPts->GetTuple(i, pt);

            //offset pt to corner of square
            pt[0]-= imageSpacing[0]*(int)(nii /2.0); pt[1]-= imageSpacing[1]*(int)(nij /2.0);
            double huSum=0;

            for(int ii=0; ii< nii; ii++) {
                for(int ij=0; ij< nij; ij++) {

                    if (!interpolator->IsInsideBuffer(pt)) {
                        if (verbose) {
                            std::cout << "Invalid Calibration (ESP): Landmark [" << i <<
                            "] is outside the image extents." << std::endl;
                        }
                        calP1 = 0.0; calP0 = 1.0; calP2 = 0.0;
                        cntrlPts = NULL;
                        return false;
                    }
                    double huValue=interpolator->Evaluate(pt);
                    huArray[i * nati + ii * nij + ij] = huValue;
                    huSum+=huValue;

                    // increment point
                    pt[1]+= imageSpacing[1];
                }
                pt[0]+= imageSpacing[0];
                pt[1]-= nij *imageSpacing[1]; // reset to starting y values
            }
        }

        // set bmd values
        for(int i=0; i< nati; i++) {
            bmdArray[0*  nati +i] = 400;
            bmdArray[1*  nati +i] = 400;
            bmdArray[2*  nati +i] = 50;
            bmdArray[3*  nati +i] = 0;
            bmdArray[4*  nati +i] = 800;
            bmdArray[5*  nati +i] = 400;
            bmdArray[6*  nati +i] = 100;
            bmdArray[7*  nati +i] = 0;
            bmdArray[8*  nati +i] = 800;
            bmdArray[9*  nati +i] = 400;
            bmdArray[10* nati +i] = 200;
            bmdArray[11* nati +i] = 0;
        }

        // calculate Least Squares Fit
        LinearRegressionCalculator::RegressionLine regressionLine = LinearRegressionCalculator::CalculateLinearRegression(huArray, bmdArray);
        calP1 =regressionLine.slope; calP0 =regressionLine.intercept; calP2 = 0.0;

        if(!isnan(calP1) && !isnan(calP0) && !isinf(calP1) && !isinf(calP2)) {
            if(verbose){std::cout<<"Completed Calibration (European Spine Phantom) [HU = m x BMD + b]: m = "<< calP1 <<", b = "<<
                        calP0 <<std::endl;}
            calBMDVals=Utilities::convertItkToVtkArray(bmdArray, "BMD values");
            calHUVals=Utilities::convertItkToVtkArray(huArray, "HU Values");
            calibrationSet=true;
        } else {
            if(verbose){std::cout<<"Unable to Complete Calibration (European Spine Phantom) [HU = m x BMD + b]: m = "<<
                        calP1 <<", b = "<< calP0 <<std::endl;}
            calP1 = 1.0; calP0 = calP2 = 0.0; cntrlPts = NULL; calHUVals = calBMDVals = cntrlVals= NULL;
        }
        phantomChanged = false;
        return calibrationSet;

    }

    bool CorticalBone::ControlPointCalibration() {

        // 1st step: work out the number of samples for each point
        int ni = cntrlNumber; // number of user specified points

        calibrationSet = false;
        cntrlVals = vtkSmartPointer<vtkDoubleArray>::New(); cntrlVals->SetNumberOfValues(ni);

        if(cntrlPts->GetNumberOfTuples() != ni) {
            cout<<"Error in CorticalBone::ManualControlPoints: Incorrect number of caibration points provided"<<endl;
            return calibrationSet;
        }

        if (verbose) {
            std::cout << "----------------- Calibration (ManualControlPts) --------------------" << std::endl;
        }

        // work out number of samples for each pt
        ImageType::SpacingType imageSpacing = image->GetSpacing();
        double sampleSquareWidth= 2*cntrlRadius;
        int nii = 2*floor((sampleSquareWidth/ imageSpacing[0])/2.0)+1; // number of x samples taken at each cal pt
        int nij = 2*floor((sampleSquareWidth/ imageSpacing[1])/2.0)+1; // number of y samples taken at each cal pt

        // 2nd step: sample with in each control point and average values at each
        for(int i=0; i< ni; i++) {

            // get control point and offset pt to corner of square
            double pt[Dimension]; cntrlPts->GetTuple(i, pt);
            if (verbose) {
                std::cout << "Cntrl["<<i<<"], x: "<<pt[0]<<", y: "<<pt[1]<<", z: "<<pt[2];
            }
            pt[0]-= imageSpacing[0]*(int)(nii /2.0); pt[1]-= imageSpacing[1]*(int)(nij /2.0);

            double valueSum=0;
            for(int ii=0; ii< nii; ii++) {
                for(int ij=0; ij< nij; ij++) {

                    if (!interpolator->IsInsideBuffer(pt)) {
                        if (verbose) {
                            std::cout <<std::endl<<Utilities::getTabString()<< "Invalid Cntrl is outside the image extents." << std::endl;
                        }
                        calP1 = 0.0; calP0 = 1.0; calP2 = 0.0; cntrlPts = NULL; cntrlVals = NULL; calHUVals = calBMDVals = cntrlVals= NULL;
                        return false;
                    }
                    valueSum += interpolator->Evaluate(pt);

                    // increment point
                    pt[1]+= imageSpacing[1];
                }
                pt[0]+= imageSpacing[0];
                pt[1]-= nij * imageSpacing[1]; // reset to starting y values
            }
            cntrlVals->SetValue(i, valueSum/(nii*nij));
            if (verbose) {
                std::cout << ", value = "<<valueSum/(nii*nij)<<endl;
            }
        }

        if (verbose) {
            std::cout << "------------------------------------------------------------------------------" << std::endl;
        }

        calibrationSet=true; phantomChanged = false;

        return calibrationSet;
    }

    bool CorticalBone::checkCalibrationControlPoints(double topBrightPt[Dimension], double topDimPt[Dimension], double bottomBrightPt[Dimension], double bottomDimPt[Dimension]) {

        //----------- Instructions: --------------//
        // Four options. Either
        // pt0 = bright top. pt1 = dim top, pt2 = bright bottom, pt3 = dim bottom
        // pt0 = bright bottom. pt1 = dim bottom, pt2 = bright top, pt3 = dim top
        // pt0 = dim top. pt1 = bright top, pt2 = dim bottom, pt3 = bright bottom
        // pt0 = dim bottom. pt1 = bright bottom, pt2 = dim top, pt3 = bright top
        //

        calibrationSet = false; calHUVals = calBMDVals = NULL; calP1 = 1.0; calP0 = 0.0; calP2 = 0.0;

        // check cal pts size
        if(cntrlPts->GetNumberOfTuples() != cntrlNumber) {
            cout<<"Error in CorticalBone::checkCalibrationControlPoints. Incorrect number of caibration points provided"<<endl;
            return calibrationSet;
        }

        // get control points
        double pt0[Dimension], pt1[Dimension], pt2[Dimension], pt3[Dimension];
        cntrlPts->GetTuple(0, pt0); cntrlPts->GetTuple(1, pt1);
        cntrlPts->GetTuple(2, pt2); cntrlPts->GetTuple(3, pt3);

        // ensure valid slice selections
        if(pt0[kZ]!=pt1[kZ] || pt2[kZ]!=pt3[kZ]) { // tops & bottoms in the same respective slices
            if(verbose){std::cout<<"Invalid Calibration (checkCalibrationControlPoints): Bottom, top or both sets of landmarks are not in the same slice."<<std::endl;}
            calP0 = 0.0; calP1 = 1.0; calP2 = 0.0; calHUVals = cntrlVals= NULL; return false;
        } else if(pt0[kZ]==pt2[kZ]) { // check tops and bottoms aren't in the same slice
            if(verbose){std::cout<<"Invalid Calibration (checkCalibrationControlPoints): Bottom and top landmarks are in the same slice."<<std::endl;}
            return false;
        }

        // check for parallel alignment
        double dCoreA[Dimension], dCoreB[Dimension];
        vtkMath::Subtract(pt2, pt0, dCoreA);  vtkMath::Normalize(dCoreA); // don't care if top to bottom, or bottom to top
        vtkMath::Subtract(pt3, pt1, dCoreB);  vtkMath::Normalize(dCoreB);
        double lcross[Dimension]; vtkMath::Cross(dCoreA, dCoreB, lcross); // cross should be close to zero
        if(vtkMath::Norm(lcross, Dimension) > 0.1) {
            if(verbose){std::cout<<"Invalid Calibration (checkCalibrationControlPoints): Bright and dim landmark sets are not parallel."<<std::endl;}
            return false;
        }


        // check landmark values are in image bounds
        InterpolatorType::Pointer interpolator = InterpolatorType::New();
        interpolator->SetInputImage(image);
        if(!interpolator->IsInsideBuffer(pt0) || !interpolator->IsInsideBuffer(pt1) || !interpolator->IsInsideBuffer(pt2) || !interpolator->IsInsideBuffer(pt3)) {
            if(verbose){std::cout<<"Invalid Calibration (checkCalibrationControlPoints): One of the landmarks is outside the image extents."<<std::endl;}
            return false;
        }

        // check the sets of bright / dim values do not vary by more than 10%
        double pt0Val, pt1Val, pt2Val, pt3Val;
        pt0Val = interpolator->Evaluate(pt0); pt1Val = interpolator->Evaluate(pt1);
        pt2Val = interpolator->Evaluate(pt2); pt3Val = interpolator->Evaluate(pt3);
        if( 2*std::abs(pt0Val-pt2Val)/(pt0Val+pt2Val) > 0.1 || 2*std::abs(pt1Val-pt3Val)/(pt1Val+pt3Val)  > 0.1 ) {
            if(verbose){std::cout<<"Warming: Calibration (checkCalibrationControlPoints): The landmark values in each set vary by more than 10%."<<std::endl;}
            if(verbose){std::cout<<"CoreA[0]="<<pt0Val<<", CoreA[1]="<<pt2Val<<", CoreB[0]="<<pt1Val<<", CoreB[1]="<<pt3Val<<std::endl;}
            //return false;
        }


        // work out bright / dim / top / bottom points
        double tb[Dimension], td[Dimension], bb[Dimension], bd[Dimension];
        if(pt0Val > pt1Val && pt2Val > pt3Val && pt0[2] > pt2[2]) { // pt0/pt2 are bright and pt0/pt1 are on top
            tb[0]=pt0[0]; tb[1]=pt0[1]; tb[2]=pt0[2];
            td[0]=pt1[0]; td[1]=pt1[1]; td[2]=pt1[2];
            bb[0]=pt2[0]; bb[1]=pt2[1]; bb[2]=pt2[2];
            bd[0]=pt3[0]; bd[1]=pt3[1]; bd[2]=pt3[2];
        } else if(pt0Val > pt1Val && pt2Val > pt3Val && pt0[2] < pt2[2]) { // pt0/pt2 are bright and pt2/pt3 are on top
            tb[0]=pt2[0]; tb[1]=pt2[1]; tb[2]=pt2[2];
            td[0]=pt3[0]; td[1]=pt3[1]; td[2]=pt3[2];
            bb[0]=pt0[0]; bb[1]=pt0[1]; bb[2]=pt0[2];
            bd[0]=pt1[0]; bd[1]=pt1[1]; bd[2]=pt1[2];
        } else if(pt0Val < pt1Val && pt2Val < pt3Val && pt0[2] > pt2[2]) { // pt1/pt3 are bright and pt0/pt1 are on top
            tb[0]=pt1[0]; tb[1]=pt1[1]; tb[2]=pt1[2];
            td[0]=pt0[0]; td[1]=pt0[1]; td[2]=pt0[2];
            bb[0]=pt3[0]; bb[1]=pt3[1]; bb[2]=pt3[2];
            bd[0]=pt2[0]; bd[1]=pt2[1]; bd[2]=pt2[2];
        } else if(pt0Val < pt1Val && pt2Val < pt3Val && pt0[2] < pt2[2]) { // pt1/pt3 are bright and pt2/pt3 are on top
            tb[0]=pt3[0]; tb[1]=pt3[1]; tb[2]=pt3[2];
            td[0]=pt2[0]; td[1]=pt2[1]; td[2]=pt2[2];
            bb[0]=pt1[0]; bb[1]=pt1[1]; bb[2]=pt1[2];
            bd[0]=pt0[0]; bd[1]=pt0[1]; bd[2]=pt0[2];
        } else {
            if(verbose){std::cout<<"Invalid Calibration (checkCalibrationControlPoints): One set of landmarks not consistantly brighter than the other."<<std::endl;}
            return false;
        }


        topBrightPt[kX]=tb[kX]; topBrightPt[kY]=tb[kY]; topBrightPt[kZ]=tb[kZ];
        topDimPt[kX]=td[kX]; topDimPt[kY]=td[kY]; topDimPt[kZ]=td[kZ];
        bottomBrightPt[kX]=bb[kX]; bottomBrightPt[kY]=bb[kY]; bottomBrightPt[kZ]=bb[kZ];
        bottomDimPt[kX]=bd[kX]; bottomDimPt[kY]=bd[kY]; bottomDimPt[kZ]=bd[kZ];

        return true;
    }

    /* Utility Methods */
    void CorticalBone::writeModelParameterToOutput(vtkIdType pointId, vtkSmartPointer<vtkTable> parameterTable, int registrationIndex) {

        vtkIdType n = parameterTable->GetNumberOfColumns();
        if(verbose){cout<<Utilities::getTabString()<<"Parameters: ";}
        for(int i=0; i<n; i++) {

            double value = parameterTable->GetValue(pointId, i).ToDouble();
            if(registrationIndex<=kEndostealRamp) {
                if(verbose){cout<<registrationMethod->GetParameterNameShort(i)<<"="<<value<<",    ";}
            } else if(registrationIndex==kHighResClassifier) {
                if(verbose){cout<< classifierMethod->GetParameterNameShort(i)<<"="<<value<<",    ";}
            } else if(registrationIndex==kCalibration) {
                if(verbose){cout<< calibrationMethod->GetParameterNameShort(i)<<"="<<value<<",    ";}
            } else {
                cerr<<"Invalid modelIndex in CorticalBone::writeModelParameterToOutput()"<<endl;
            }
        }
        if(verbose){cout<<endl;}
    }

    void CorticalBone::insertParametersValue(ParametersType & array, CoordinateType value, unsigned int index) {

        ParametersType newArray(array.GetSize()+1);

        for(int i=0; i<array.GetSize(); i++) {

            if(i<index) {
                newArray[i]=array[i];
            } else {
                newArray[i+1]=array[i];
            }
        }
        newArray[index] = value;

        array = newArray;

    }

    void CorticalBone::removeParametersValue(ParametersType & array, unsigned int index) {

        ParametersType newArray(array.GetSize()-1);

        for(int i=0; i<newArray.GetSize(); i++) {

            if(i<index) {
                newArray[i]=array[i];
            } else {
                newArray[i]=array[i+1];
            }
        }

        array = newArray;

    }

    unsigned int  CorticalBone::getFlexIndex(unsigned int absoluteIndex) {
        int flexIndex=0; // get index of the fixed cb value
        for(int i=0; i<absoluteIndex; i++) {
            flexIndex += transformFlexMap[i];
        }

        return flexIndex;
    }

    void CorticalBone::writeModelResults(vtkIdType pointId) {

        int n;

        if(modelIndex<=kEndostealRamp) {

            if(verbose){cout<<Utilities::getTabString()<<"Display Values:   ";}
            n = getNumberOfDisplays();
            for(int i=0; i<n; i++) {
                if(verbose){cout<<registrationMethod->GetDisplayNameShort(i, parametersImported)<<"="<<
                                                                                                  modelDisplayTable->GetValue(pointId, i).ToDouble()<<"    ";}
            } if(verbose){cout<<endl;}

        } else if(modelIndex==kHighResClassifier) {

            if(verbose){cout<<Utilities::getTabString()<<"Display Values:   ";}
            n = getNumberOfDisplays();
            for(int i=0; i<n; i++) {
                if(verbose){cout<< classifierMethod->GetDisplayNameShort(i)<<"="<< modelDisplayTable->GetValue(pointId, i).ToDouble()<<"    ";}
            } if(verbose){cout<<endl;}

            if(verbose){cout<<Utilities::getTabString()<<"Curve Values:   ";}
            n = getNumberOfParameters();
            for(int i=0; i<n; i++) {
                if(verbose){cout<< classifierMethod->GetParameterNameShort(i)<<"="<<modelParametersTable->GetValue(pointId, i).ToDouble()<<"    ";}
            } if(verbose){cout<<endl;}
        } else if(modelIndex==kCalibration) {
            if(verbose){cout<<Utilities::getTabString()<<"Display Values:   ";}
            n = getNumberOfDisplays();
            for(int i=0; i<n; i++) {
                if(verbose){cout<< calibrationMethod->GetDisplayNameShort(i)<<"="<< modelDisplayTable->GetValue(pointId, i).ToDouble()<<"    ";}
            } if(verbose){cout<<endl;}

        } else {
            cout<<"Error in CorticalBone::writeModelResults; invalid parameter"<<endl;
        }

    }

    /* Getters */

    std::string CorticalBone::getPhantomName(int index) {

        if(index == kMindwaySolidCal) {

            std::ostringstream p0Stream; std::ostringstream p1Stream;
            p0Stream << calP0; std::string p0String = p0Stream.str();
            p1Stream << calP1; std::string p1String = p1Stream.str();

            return "Mindways Solid K2HPO4: p0="+p0String+", p1="+p1String;
        } else if(index == kBoneDensityCal) {

            std::ostringstream p0Stream; std::ostringstream p1Stream;
            p0Stream << calP0; std::string p0String = p0Stream.str();
            p1Stream << calP1; std::string p1String = p1Stream.str();

            return "Bone Density: p0="+p0String+", p1="+p1String;
        } else if(index == kManualQuadraticCal) {
            std::ostringstream p0Stream, p1Stream, p2Stream;
            p0Stream << calP0; std::string p0String = p0Stream.str();
            p1Stream << calP1; std::string p1String = p1Stream.str();
            p2Stream << calP2; std::string p2String = p2Stream.str();
            return "Manual Quadratic Calibration: p0="+p0String+", p1="+p1String+", p2="+p2String;

        } else if(index == kEuropeanSpineCal) {

            std::ostringstream p0Stream, p1Stream;
            p0Stream << calP0; std::string p0String = p0Stream.str();
            p1Stream << calP1; std::string p1String = p1Stream.str();
            return "European Spine Phantom: p0="+p0String+", p1="+p1String;

        } else if(index == kManualLinearCal) {
            std::ostringstream p0Stream, p1Stream;
            p0Stream << calP0; std::string p0String = p0Stream.str();
            p1Stream << calP1; std::string p1String = p1Stream.str();
            return "Manual Linear Calibration: p0="+p0String+", p1="+p1String;

        } else if (index == kManualControlPtsCal) {
            return "Manual Control Point Calibration";
        } else if(index == kNoCal) {
            return "No Calibration Phantom";
        } else {
            std::cerr<<"Error: CorticalBone::getPhantomName() invalid phantom index = "<<index<<std::endl;
            return "Invalid Selection";
        }
    }

    std::string CorticalBone::getModelName(int index) {
        if(!meshSet || !imageSet) {
            return "N/A";
        } else if(index == kThreeTierRect) {
            return "Three Tier Rectangular";
        } else if(index == kEndostealRamp) {
            return "Endosteal Ramp Model";
        } else if(index == kHighResClassifier) {
            return "High Res Classifier";
        } else if(index == kCalibration) {
            return "Calibration";
        } else if(index == kFabricated) {
            return "Fabricated Values";
        } else {
            std::cerr<<"Error: CorticalBone::getModelName() invalid model index: "<<index<<std::endl;
            return "Invalid Selection";
        }
    }

    std::string CorticalBone::getModelName() {
        return getModelName(modelIndex);
    }

    std::string CorticalBone::getSchemeName() { // TODO - combine with 'GetModelSchemeSelection()'

        if(!meshSet || !imageSet) {
            return "N/A";
        } else if(modelIndex<=kEndostealRamp) {
            if(fittingSchemeIndex== kStdAFitting) {
                return "Std A Fitting";
            } else if(fittingSchemeIndex== kStdBFitting) {
                return "Std B Fitting";
            } else if(fittingSchemeIndex== kStdCFitting) {
                return "Std C Fitting";
            } else if(fittingSchemeIndex== kStdDFitting) {
                return "Std D Fitting";
            } else if(fittingSchemeIndex==kCBMV2AFitting) {
                return "CBM Version A";
            } else if(fittingSchemeIndex==kCBMV2BFitting) {
                return "CBM Version B";
            } else if(fittingSchemeIndex==kCBMV2CFitting) {
                return "CBM Version C";
            } else if(fittingSchemeIndex==kCBMV2DFitting) {
                return "CBM Version D";
            } else if(fittingSchemeIndex==kUnconstrainedAFitting) {
                return "Unconstrained Version A";
            } else if(fittingSchemeIndex==kUnconstrainedBFitting) {
                return "Unconstrained Version B";
            } else if(fittingSchemeIndex==kUnconstrainedCFitting) {
                return "Unconstrained Version C";
            } else if(fittingSchemeIndex==kCBSmoothingAFitting) {
                return "Smoothing Version A";
            } else if(fittingSchemeIndex==kCBSmoothingBFitting) {
                return "Smoothing Version B";
            } else if(fittingSchemeIndex==kCBSmoothingCFitting) {
                return "Smoothing Version C";
            } else if(fittingSchemeIndex==kCBSmoothingDFitting) {
                return "Smoothing Version D";
            } else {
                std::cerr<<"Error: CorticalBone::getSchemeName() invalid model index: "<<fittingSchemeIndex<<std::endl;
                return "Invalid Selection";
            }
        } else {
            return classifierMethod->GetModelSchemeName();
        }
    }

    std::string CorticalBone::getOptimisationName() {

        if(!meshSet || !imageSet) {
            return "N/A";
        }

        int index = registrationMethod->GetOptimiserSelection();
        if(index == kLMOptimiser) {
            return "LM Optimiser";
        } else if(index == kPowellOptimiser) {
            return "Powell Optimiser";
        } else if(index == kEvolutionaryOptimiser) {
            return "Evolutionary Optimiser";
        } else {
            std::cerr<<"Error: CorticalBone::getOptimisationName() invalid optimsation index = "<<index<<std::endl;
            return "Invalid Selection";
        }
    }

    std::string CorticalBone::getCBDensityState() {

        if(!meshSet || !imageSet) {
            return "N/A";
        } else if(registrationMethod->GetFWHMMode()) {
            return "FWHM CB Density";
        } else if(!isnan(globalFixedCB)) { // fixed CB

            std::ostringstream fixedCBStream; std::string fixedCBString;
            fixedCBStream << globalFixedCB; fixedCBString = fixedCBStream.str();
            return "Fixed CB Density = "+fixedCBString;
        } else {
            return "Free CB Density";
        }

    }

    std::string CorticalBone::getParameterName(int index) {
        std::string parameterName;
        if(modelIndex<= kEndostealRamp) {
            parameterName=registrationMethod->GetParameterName(index);
        } else if (modelIndex==kHighResClassifier) {
            parameterName=classifierMethod->GetParameterName(index);
        } else if (modelIndex==kCalibration) {
            parameterName=calibrationMethod->GetParameterName(index);
        } else {
            cerr<<"Invalid index in CorticalBone::getParameterName "<<index<<endl;
            parameterName=std::string("");
        }
        return parameterName;
    }

    int CorticalBone::getNumberOfParameters() {

        int n;
        if(modelIndex<= kEndostealRamp) {
            n=registrationMethod->GetNumberOfCombinedParameters(modelIndex);
        } else if (modelIndex==kHighResClassifier) {
            n=classifierMethod->GetNumberOfParameterValues();
        } else if (modelIndex==kCalibration) {
            n=calibrationMethod->GetNumberOfParameterValues();
        } else {
            cerr<<"Invalid index in CorticalBone::getNumberOfParameters "<<modelIndex<<endl;
            n=kInvalid;
        }
        return n;
    }

    std::string CorticalBone::getDisplayName(int index) {
        std::string parameterName;
        if(modelIndex<= kEndostealRamp) {
            parameterName=registrationMethod->GetDisplayName(index, parametersImported);
        } else if (modelIndex==kHighResClassifier) {
            parameterName=classifierMethod->GetDisplayName(index);
        } else if (modelIndex==kCalibration) {
            parameterName=calibrationMethod->GetDisplayName(index);
        } else {
            cerr<<"Invalid index in CorticalBone::getParameterName "<<index<<endl;
            parameterName=std::string("");
        }
        return parameterName;
    }

    int CorticalBone::getNumberOfDisplays() {

        int n;
        if(modelIndex<= kEndostealRamp) {
            n=registrationMethod->GetNumberOfDisplayValues(parametersImported);
        } else if (modelIndex==kHighResClassifier) {
            n=classifierMethod->GetNumberOfDisplayValues(parametersImported);
        } else if (modelIndex==kCalibration) {
            n=calibrationMethod->GetNumberOfDisplayValues(parametersImported);
        } else {
            cerr<<"Invalid index in CorticalBone::getNumberOfParameters "<<modelIndex<<endl;
            n=kInvalid;
        }
        return n;
    }

}

