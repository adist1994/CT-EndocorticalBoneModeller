/* 
 * File:   corticalbone.h
 * Author: rap58
 *
 * Created on 02 October 2014, 13:11
 */

#ifndef CORTICALBONE_H
#define	CORTICALBONE_H

/* General notes
 * 1. contains both itk and vtk objects
 * 2. pass in only file name to open and save files
 * 3. requires an interface for viewing openGL windows */

// itk + namespace
#include <itkImage.h>

// mesh type headers
#include <vtkVRMLImporter.h>
#include <vtkSmartPointer.h>

#include <itkImageToVTKImageFilter.h>
#include <itkImageSeriesReader.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkQCTImageIO.h>
#include <itkLinearInterpolateImageFunction.h> 
#include <itkMinimumMaximumImageCalculator.h>

#include <itkModelRegistrationMethodLocal.h>

#include <itkArray2D.h>


#include "vtkPolyData.h"
#include "profileClassifier.h"
#include "ProfileCalibration.h"

//------------------ class predefine ------------------------//
class vtkKdTreePointLocator;

class vtkSphereSource;
class vtkLineSource;
class vtkTable;
class vtkColorTransferFunction;

//------------------ type defines -------------------//

// Image
typedef signed short PixelType;
typedef double CoordinateType;
const unsigned int Dimension = 3;
typedef itk::Image< PixelType, Dimension > ImageType; // CT Image pixel details
typedef itk::Point< CoordinateType, Dimension > PointType;
typedef itk::ImageSeriesReader< ImageType > ImageReaderType;
typedef itk::ImageToVTKImageFilter<ImageType>       ConnectorType;
typedef itk::LinearInterpolateImageFunction< ImageType, float > InterpolatorType; //TODO - BSplineInterpolateImageFunction
typedef itk::MinimumMaximumImageCalculator< ImageType > ImageMinMaxCalculatorType;
//DICOM
typedef itk::GDCMImageIO  ImageDICOMIOType; // DICOM v3 class = GDCMImageIO
typedef itk::GDCMSeriesFileNames DICOMNamesGeneratorType; // Slice names given directory
typedef std::vector<std::string>    FileNamesContainer; // File names
//QCT
//#include <itkRawImageIO.h>
//typedef itk::RawImageIO< PixelType, Dimension > ImageRawIOType;
typedef itk::QCTImageIO< PixelType, Dimension > ImageQCTIOType;
// TIF
//typedef itk::NumericSeriesFileNames NumericNameGeneratorType;
// MODEL
// other
typedef itk::ProfileSpatialObject< Dimension >   ProfileType;
typedef itk::ModelTransformBase<CoordinateType, 1> TransformBaseType;
typedef TransformBaseType::ParametersType ParametersType; // itk::Array<double>
// registration methods
typedef itk::ModelRegistrationMethod<CoordinateType, Dimension> ModelRegistrationType;
typedef itk::ProfileClassifier ThresholdRegistrationType;
typedef itk::ProfileCalibration CalibrationType;

namespace itk {

    class CorticalBone  : public ProcessObject {

        // no silent errors. Always raise a warning or error.
        /* TODO - Priority
         * 1. make HRTC threshold selection robust to italy plastic wrapping
         * 2. include Ycb (and later Xp) variances in the calculation of smoothing weights
         * 3. smooth and then constrain Xp (and Ycb)
         * 4. support 2-peak model over femoral head joint space (defined by a mask of the mesh)
         *  */
        /* TODO - Long Term
         * 1. consider alignment of imported profiles
         * 2. support specification of profile length
         * 3. decided on a consistent approach to calibration (return calibrated results, require calibrated fixed cb)??
         * 4. Consider image value for profiles outside of image -> Single/MultipleModelMetric::ResampleImageProfile()
         * 5. Make it so that all of my model processing methods have a common inheritant
         */

        /* Define API
         * 1. Getters
         *  a. Single points written to double* argument
         *  b. Multiple points returned using vtkPoints
         *  c. throw exception if return variable unset
         * 2. Setters
         *  a. set class variable - state change
         *  b. throw exception if write error - void return
         * 3. Run
         *  a. use set class variables to run an operation - state change
         *  b. throw exception if run error (unset class variables) - void return */


    public:

        typedef CorticalBone                            Self;
        typedef ProcessObject                           Superclass;
        typedef SmartPointer< Self >                    Pointer;
        typedef SmartPointer< const Self >              ConstPointer;

        // standard class macros
        itkNewMacro(Self);
        itkTypeMacro(CorticalBone, ProcessObject);

        CorticalBone(); // constructor
        void setVerbose(bool verboseIn);
        bool reset();

        //---- I/O ---//
        // loaders
        bool loadMesh(std::string fileName);
        bool loadImage(std::string fileName);
        bool loadImage(std::vector<std::string> nameList, double zSpacing = nan("1"), double zOffset = nan("1"));
        bool loadValueArrays(std::string baseFileName);

        // savers
        bool saveMesh(std::string fileName);
        bool saveValueArrays(std::string baseFileName);
        bool saveCalibration(std::string baseFileName);


        //---- setters ---//
        // set values
        void setFixedCBDensity(double fixedCBDensity, bool setGlobalFixedDensity=true);
        void setClassifierThresholdInfo(double softTissue, double corticalBone, double threshold, int classifierIndex);
        void setClassifierThresholdInfo(double softTissue, double corticalBone, double threshold, double weight, int classifierIndex);
        bool setFixedSampleNumber(int sampleNumber);
        void setSmoothingRadius(double);
        bool setCalibration(double p0, double p1, double p2);
        bool setCalibrationPtGeometry(double radius, int number);
        bool setImportProfiles(std::string fileName);

        // remove values
        void removeFixedCBDensity(bool clearGlobalFixedDensity=true);
        void removeFixedSampleMode();
        void removeThresholds();
        void removeImportedProfile();

        // selections
        void setDisplayIndex(int index);
        void setPhantomIndex(int index);
        void setModelIndex(int index);
        void setOptimiserIndex(int index);
        void setFittingSchemeIndex(int index);

        //--- getters ---//
        int getState(); // 0=none open, 1=mesh open, 2=image open, ? cal or if phantom index changed, 3=both open, ? cal, 4=processed

        double getProfileLength();
        double getPeriostealOffset();

        // values set by user
        double getFixedCB();
        int getSampleNumber();
        bool getFixedSigma(double &x, double &y, double &z);
        double getSmoothingValue();
        bool getProfileAveragingValues(double &x, double &y, double &z);
        bool getCalibrationValues(double& p0, double& p1, double& p2);
        bool getCalibrationPoints(int &phantomType, vtkSmartPointer<vtkDoubleArray> &pts);
        bool getCalibrationPhantom(std::string& phantomType, int& phantomIndex);
        bool getCalibrationPtGeometry(double &radius, int &number);
        vtkSmartPointer<vtkDoubleArray> getControlPointValues();
        bool getClassifierInfo(double &softTissue, double &corticalBone, double &threshold, int &classifierIndex, std::string &name);
        bool getClassifierInfo(double &softTissue, double &corticalBone, double &threshold, double &weight, int &classifierIndex, std::string &name);
        bool getclassifierName(std::string &name);
        std::string getDisplayName();

        // seclections set by user
        void getModelSelection(int& modelIndex, std::string& modelName);
        void getOptimiserSelection(int& optimiserIndex, std::string& optimiserName);
        void getSchemeSelection(int& schemeIndex, std::string& schemeName);
        int getDisplayIndex();
        bool getClassifierIndex(double &index);

        // 3D display API
        vtkSmartPointer<vtkPolyData> getMesh();
        vtkSmartPointer<vtkPolyData> getPeriostealMesh();
        ImageType::Pointer getDICOM();

        // mesh measurement API
        void getDisplayRange(double &min, double &max);
        vtkSmartPointer<vtkDoubleArray> getDisplayArray();

        // pt measurement API
        vtkIdType getClosestPoint(double point[Dimension]);
        void getMeasurementLocation(double location[Dimension]);
        void getProfileStart(double start[Dimension]);
        void getProfileEnd(double end[Dimension]);
        vtkSmartPointer<vtkPoints> getProfilePoints();

        // profile display API
        vtkSmartPointer<vtkTable> getImageProfileTable(vtkIdType pointId);
        vtkSmartPointer<vtkTable> getImportedProfileTable(vtkIdType pointId);
        vtkSmartPointer<vtkTable> getDisplayModelTable(vtkIdType pointId);
        vtkSmartPointer<vtkTable> getWeightsTable(vtkIdType pointId);
        vtkSmartPointer<vtkTable> getImportedDisplayModelTable(vtkIdType pointId);
        vtkSmartPointer<vtkTable> getProcessingModelTable(vtkIdType pointId);
        vtkSmartPointer<vtkTable> getImportedProcessingModelTable(vtkIdType pointId);
        vtkSmartPointer<vtkTable> getErrorAreaTable(vtkIdType pointId);


        //--- get state ---//
        bool isSigmaFixed();
        bool isCBFixed();
        bool isSetToFWHMMode();
        bool isMeshSet();
        bool isImageSet();
        bool isMeshMeasured();
        bool isParametersImported();
        bool isClassifierMode();
        bool areResultsLoaded();
        bool areThresholdsSet();
        bool isSampleNumberFixed();
        bool isProfileAveragingOn();
        bool isSmoothingOn();
        bool isCalibrated();

        //------ change state -----//
        void turnOnFWHMMode();
        void turnOffFWHMMode();
        void turnOffFixedSigma();
        void turnOnFixedSigma(double x, double y, double z);
        void turnOffProfileAveraging();
        void turnOnProfileAveraging(double x, double y, double z);
        void turnOffSmoothingMode();
        void turnOnSmoothingMode();

        //--- Processing ---//
        bool calculateThresholds(int thresholdIndex);
        bool calculateThresholds(int thresholdIndex, double weight);
        bool runModellingOverMesh(); // fit CBM model to each point on a mesh
        bool runModellingAtPoint(vtkIdType pointId);  // fit CBM model to a single point - return true = valid model, false - invalid model
        bool runCalibration(vtkSmartPointer<vtkDoubleArray> calPts = vtkSmartPointer<vtkDoubleArray>::New());

        //------- Enums -------//
        typedef enum { // model type
            kThreeTierRect=0,
            kEndostealRamp,
            kHighResClassifier,
            kCalibration,
            kFabricated,
        } modelOptions;
        typedef enum { // fitting schemes
            kStdAFitting =0, // dynamic hueristic weighting
            kStdBFitting, // dynamic narrow weighting
            kStdCFitting, // dynamic wide weighting
            kStdDFitting, // dynamic narrow->wide weighting
            kCBMV2AFitting, // published CBMv2 (both fittings to selected model) - Static Hueristic weighting function
            kCBMV2BFitting, // published CBMv2 (both fittings to selected model) - Static Narrow weighting function
            kCBMV2CFitting, // published CBMv2 (both fittings to selected model) - Static Narrow then wide weighting
            kCBMV2DFitting, // published CBMv2 (both fittings to selected model) - Dynamic Narrow then wide weighting
            kUnconstrainedAFitting, // constrained to unconstrained - 1st constrained rect, the unconstrained (selected model)
            kUnconstrainedBFitting, // constrained to unconstrained - 1st constrained (seleced model), the unconstrained (selected model)
            kUnconstrainedCFitting, // constrained to unconstrained - 1st constrained (seleced model), the unconstrained (selected model)
            kCBSmoothingAFitting, // CB smoothing of kUnconstrainedAFitting (rect to selected based unconstrained) to precision/distance/error weighted smoothing using selected model
            kCBSmoothingBFitting, // CB smoothing of kUnconstrainedBFitting (selected to selected based unconstrained) to precision/distance/error weighted smoothing using selected model
            kCBSmoothingCFitting, // CB smoothing of kUnconstrainedAFitting (rect to rect based unconstrained) to precision/distance/error weighted smoothing using selected model
            kCBSmoothingDFitting, // kCBSmoothingBFitting followed by a second cb constraint smoothing
        } fittingSchemeOptions;
        typedef enum { // run status
            kMasked=-3,
            kOutOfBounds, // -2,
            kInvalid, // -1
            kUndefined,
            kValid,
            kDblPeak,
        } GeneralState;
        typedef enum { // calibration
            kNoCal=0,
            kMindwaySolidCal,
            kBoneDensityCal,
            kEuropeanSpineCal,
            kManualControlPtsCal,
            kManualLinearCal,
            kManualQuadraticCal,
            kFileSpecifiedCal,
        } calibrationOptions;

    private:

        //--- Private Enumeration's ---//
        typedef enum { // model type
            kLMOptimiser=0,  // TODO - wrap into the registraion class
            kPowellOptimiser,
            kEvolutionaryOptimiser,
        } optimiserOptions;
        typedef enum { // calibration
            kX=0,
            kY,
            kZ,
        } DimensionIndicies;


        //--- Setup ---//
        // general
        void initialiseSettings();

        // object creation / setup
        bool setupImage();
        void createParameters();
        void createRegistrationMethods(); // call once image and mesh added
        void createModelProfile();
        void createFittingMethod();
        void createClassifierMethod();
        void createCalibrationMethod();




        //--- loading ---//
        // Meshes
        bool openVRML(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn);
        bool openVRML_SW(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn);
        bool openSTL(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn);
        bool openPLY(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn);
        bool openOBJ(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn);

        // Images
        bool openDICOM(std::string filePath);
        bool openRAW(std::string fileName);
        bool openQCT(std::string fileName);
        void openTIF(std::vector<std::string> fileNames, double zSpacing, double zOffset);

        // masks
        bool openMeshMask(std::string fileName);
        bool openMaskTxtFile(std::string fileName);
        bool openMaskSelFile(std::string fileName);

        // results
        bool openParameters(std::string fileName);
        bool openImportParameters(std::string fileName);
        bool openDisplays(std::string fileName);
        bool openProfiles(std::string fileName,vtkSmartPointer<vtkDoubleArray> &array,std::string name,int modelIndex);


        //--- saving ---//
        // Meshes
        bool saveVRML(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn);
        bool saveSTL(std::string fileName, vtkSmartPointer<vtkPolyData> meshIn);
        bool saveOBJ(std::string filename, vtkSmartPointer<vtkPolyData> meshIn);

        // results
        bool saveProfiles(std::string fileName, vtkSmartPointer<vtkDoubleArray> array, std::string description);


        //--- processing ---//

        // CBM
        void fitModel(int localModelIndex, ParametersType parameters, int index=TransformBaseType::kHueristicWeighting, double scale=1.0/25000.0, bool dynamicWeighting=false);

        // pipelines
        void runModelFittingPipeline(vtkIdType meshId, int schemeIndex, int modelIndex);
        void runThresholdingPipeline(vtkIdType meshId);
        void runGenerateCalibrationValues(vtkIdType meshId);
        void updateModelProfile(vtkIdType meshId, bool runningOverMesh=true);

        // over mesh
        bool runCBMOverMeshWithCBSmoothing(); // fit CBM model to each point on a mesh
        bool runCBMOverMeshWithNoSmoothing(); // fit CBM model to each point on a mesh
        bool runCalibrationOverMesh();
        bool runClassifierOverMesh(); // fit HRTC model to each point on a mesh
        bool estimateMinBasedThresholdsOverMesh(int thresholdIndex);
        bool estimateMedianBasedThresholdsOverMesh(int thresholdIndex);

        // CBM schemes
        void runV1aFittingPipeline(); // rect then selected model - narrow weighting
        void runV1bFittingPipeline(); // rect then selected model - wide weighting
        void runV1cFittingPipeline(); // rect then selected model - narrow->wide weighting
        void runV1dFittingPipeline(); // rect then selected model - wide weighting
        void runV2aFittingPipeline(); // selected model - narrow weighting
        void runV2bFittingPipeline(); // selected model - wide weighting
        void runV2cFittingPipeline(); // selected model - narrow->wide weighting
        void runV2dFittingPipeline(); // selected model - dynamic hueristic weighting
        void runV3aFittingPipeline(int selectedModel); // rect then selected model - narrow weighting
        void runV3bFittingPipeline(int selectedModel); // rect then selected model - wide weighting
        void runV3cFittingPipeline(int selectedModel); // rect then selected model - narrow->wide weighting

        // post-processing results
        void calculateImportErrors(vtkIdType pointId);
        double calculateImportOffset(vtkIdType pointId);

        // mesh processing
        bool calculatePeriostealMesh();
        void resetPeriostealMesh();
        bool updatePeriostealMeshNode(vtkIdType meshNode);

        // calibration curves
        bool mindwaysCalibration();
        bool BoneDensityCalibration();
        bool noPhantomCalibration();
        bool ESPCalibration();
        bool ControlPointCalibration();

        bool checkCalibrationControlPoints(double topBrightPt[Dimension], double topDimPt[Dimension], double bottomBrightPt[Dimension], double bottomDimPt[Dimension]);

        //--- Results ---//
        // text out
        void writeModelResults(vtkIdType pointId);
        void writeModelParameterToOutput(vtkIdType pointId, vtkSmartPointer<vtkTable> parameterTable, int registrationIndex);

        // clear results
        void clearPreviousResults();
        void clearCalibration();
        void clearParametersTable();
        void clearMeshAttributeTable();
        void clearImageProfilesArray();
        void maskResultArrays();

        // store results
        void storeDisplayValues(vtkIdType pointId);
        void storeParameterValues(vtkIdType pointId);
        void storeSampledProfiles(vtkIdType pointId);


        //--- Getters ---//
        // names given indicies
        std::string getPhantomName(int index);
        std::string getModelName(int index);
        std::string getModelName();
        std::string getSchemeName();
        std::string getOptimisationName();
        std::string getCBDensityState();
        std::string getParameterName(int index);
        std::string getDisplayName(int index);

        // result info
        int getNumberOfParameters();
        int getNumberOfDisplays();


        //--- Setters ---//
        void setImageBuffer(); // base upon the size of the mesh
        void setSampling();
        void setModelToImport(vtkIdType pointId); // todo - high res only generalise or remove
        void setFixedSigma();


        //-- Get State ---//
        bool isImportedModelValid(vtkIdType pointId); // TODO - generalise to classifier/threshold models or remove

        //--- Change State ---//
        void turnOnRegistrationSigmaCorrection();
        void turnOffRegistrationSigmaCorrection();
        void resetCBConstraint();
        void removeFixedSigma(bool clearDisplayBool = true);


        //--- Utilities ---//
        // opening
        bool readTableValues(ifstream& file, vtkSmartPointer<vtkTable>& table);
        bool readArrayValues(ifstream& file, vtkSmartPointer<vtkDoubleArray>& array);

        // saving
        bool saveTable(std::string fileName, vtkSmartPointer<vtkTable> table, std::string description);
        bool saveArray(std::string fileName, vtkSmartPointer<vtkDoubleArray> array, std::string header_line1, std::string header_line2);

        // parameter conversion
        ParametersType convertRectToRampParameters(ParametersType rectParameters);
        ParametersType convertRectToRampScales(ParametersType rectScales);
        ParametersType getTableRow(vtkSmartPointer<vtkTable> table, vtkIdType rowId);

        // parameter manipulation
        void insertParametersValue(ParametersType & array, CoordinateType value, unsigned int index);
        void removeParametersValue(ParametersType & array, unsigned int index);
        unsigned int getFlexIndex(unsigned int allIndex);

        //---- Private Objects ---/
        bool verbose;
        // image components
        ImageType::Pointer image;
        InterpolatorType::Pointer interpolator;
        PixelType maxImageValue;
        const PixelType maxBMDValue = 1600; // units are mg/cm^3
        vtkSmartPointer<vtkKdTreePointLocator> meshTree;
        vtkSmartPointer<vtkDataArray> meshNormals;
        //  ImageType::PointType bufferOffset; // used in an attempt to implement the buffered region
        // mesh Components
        vtkSmartPointer<vtkPolyData> mesh, periostealMesh;
        ParametersType meshMask;

        // model components
        ParametersType rectFlexParameters;
        ParametersType rectFlexScales;
        ParametersType rectAllScales;
        ParametersType rectAllParameters;
        ParametersType transformFlexMap;

        ParametersType rampFlexScales;
        ParametersType rampAllScales;

        ProfileType::Pointer modelProfile;

        // Registration classes
        ThresholdRegistrationType::Pointer classifierMethod;
        ModelRegistrationType::Pointer registrationMethod;
        CalibrationType::Pointer calibrationMethod;

        // storing parameter/profile/display values
        ParametersType modelFitStatusArray;
        vtkSmartPointer<vtkTable> modelParametersTable;
        vtkSmartPointer<vtkTable> modelDisplayTable;
        vtkSmartPointer<vtkDoubleArray> imagePositionsArray;
        vtkSmartPointer<vtkDoubleArray> imageProfilesArray;
        vtkSmartPointer<vtkDoubleArray> imageClassificationArray; // classifier - classified as not bone, tb bone or bone
        vtkSmartPointer<vtkDoubleArray> imagePercentageArray; // classifier - % above or below a threshold

        // importing parameter/profile values
        vtkSmartPointer<vtkTable> importedParametersTable;
        vtkSmartPointer<vtkDoubleArray> importedPositionsArray;
        vtkSmartPointer<vtkDoubleArray> importedProfilesArray;
        vtkSmartPointer<vtkDoubleArray> importedClassificationArray; // classifier - classified as not bone, tb bone or bone
        vtkSmartPointer<vtkDoubleArray> importedPercentageArray; // classifier - % above or below a threshold

        // calibration components - only updated when 'runCalibration' called so all in sync
        double calP1, calP0, calP2; double cntrlRadius; unsigned long cntrlNumber, calNumber;
        vtkSmartPointer<vtkDoubleArray> cntrlPts;
        vtkSmartPointer<vtkDoubleArray> cntrlVals, calHUVals, calBMDVals;

        // imported parameters density scaling parameters
        double p0, p1, p2;

        // constraints
        double globalFixedCB;

        // state
        bool imageSet, meshSet, meshMeasured, calibrationSet, phantomChanged, parametersImported, sampleNumberFixed, resampleAboutPeriostealEdge;
        int displayIndex, phantomIndex, modelIndex, fittingSchemeIndex, importedModelIndex;

        // profile info
        int numberOfSamples;
        double profileLength, defaultPeriostealEdgeRatio, maximumOffset;

        // classifier thresholds
        double thresholdWeight;

        // smoothing radius
        double smoothingRadius; bool smoothingSet;

        // previous run
        int nanCount, failureCount, outOfBoundsCount, maskedCount;

        // distinguish between loaded and generated results
        bool loadedResults;

    };

}

#endif	/* CORTICALBONE_H */

