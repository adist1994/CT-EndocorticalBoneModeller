/* 
 * File:   visualisationframe.h
 * Author: rap58
 *
 * Created on 01 October 2014, 16:56
 */

#ifndef VISUALISATIONFRAME_H
#define	VISUALISATIONFRAME_H

#include <QFrame>
#include <vtkVRMLImporter.h>
#include <vtkSmartPointer.h>

#include "corticalbone.h"

class QVTKWidget;
class CorticalBone;
class vtkResliceImageViewer;
class vtkImagePlaneWidget;
class vtkDICOMImageReader;
class vtkPolyData;
class vtkImageData;



//---------------- set type defs and constants --------------------// todo - come up with a better way of specifyting these values
static const int sliceAxis = 2;
static const int dimension = 3;

class QVisualisationFrame : public QFrame
{
Q_OBJECT

public:
    QVisualisationFrame(QWidget *parent = 0);
    bool reset();
    bool setImage(QString baseFileName);
    bool setMesh(QString *fileName);
    bool saveMesh(QString *fileName);
    bool runCalibration();
    bool setCalibration(double p0, double p1, double p2);
    bool runModellingOverMesh();
    bool saveValueArrays(QString *fileName);
    bool loadValueArrays(QString fileName);
    bool saveCalibration(QString *fileName);
    bool saveDisplays(QString *sliceName, QString *profileName, QString *threeDimName);
    bool openParameters(QString *fileName);


    void disablePtMeasures();
    void enablePtMeasures();
    void enableCalibrationMode();
    void disableCalibrationMode();
    void setDisplayParameter(int index);
    void setDisplayMesh(int index);
    void setCalibrationPhantom(int index);
    void setModelFunction(int index);
    void setOptimiser(int index);
    void setfittingScheme(int index);
    void setFixedCBDensity(double CBDensity);
    void turnOnFWHMMode();
    void turnOffFWHMMode();
    bool isSetToFWHMMode();
    void removeFixedCBDensity();
    int getState();

    // fixed sigma value
    void turnOffFixedSigma();
    void turnOnFixedSigma(double x, double y, double z);
    bool isSigmaFixed();
    void getFixedSigma(double &x, double &y, double &z);

    // thresholds
    bool calculateThresholds(int classifierThresholdIndex);
    bool calculateThresholds(int classifierThresholdIndex, double weight);
    bool setClassifierThresholdInfo(double stDensity, double cbDensity, double thresholdDensity, int thresholdIndex);
    bool setClassifierThresholdInfo(double stDensity, double cbDensity, double thresholdDensity, double weight, int thresholdIndex);
    bool getClassifierThresholdInfo(double &stDensity, double &cbDensity, double &thresholdDensity, int &classifierThresholdIndex, QString& classifierThresholdIndexName);
    bool getClassifierThresholdInfo(double &stDensity, double &cbDensity, double &thresholdDensity, double &weight, int &classifierThresholdIndex, QString& classifierThresholdIndexName);
    void removeThresholds();
    bool areThresholdsSet();

    // profile averaging
    void turnOnProfileAveraging(double x, double y, double z);
    void turnOffProfileAvergaing();
    bool isProfileAveragingOn();
    bool getProfileAveragingValues(double &x, double &y, double &z);

    // fixed sample number
    void removeFixedSampledMode();
    bool isSampleNumberFixed();
    int getSampleNumber();
    bool setFixedSampleNumber(int sampleNumber);

    // model info
    void getModelSelection(int& modelIndex, QString& modelName);
    void getOptimiserSelection(int& optimiserIndex, QString& optimiserName);
    void getSchemeSelection(int& schemeIndex, QString& schemeName);

    // smoothing
    void setSmoothingRadius(double);
    void turnOffSmoothingMode();
    void turnOnSmoothingMode();
    bool isSmoothingOn();
    double getSmoothingValue();

    // get / set calibration values
    bool isCalibrated();
    vtkSmartPointer<vtkDoubleArray> getCalibrationPoints();
    vtkSmartPointer<vtkDoubleArray> getCalibrationValues();
    bool getCalibrationValues(double& p0, double& p1, double& p2);
    bool getCalibrationPhantomType(std::string& phantomType, int& phantomIndex);
    bool setCalibrationPoints(int phantomIndex, vtkSmartPointer<vtkDoubleArray> pts);
    bool getCalibrationPtGeometry(double &radius, int &number);
    bool setCalibrationPtGeometry(double radius, int number);

    // import profiles
    void removeImportedProfile();

    bool isFixedCB();
    double getFixedCB();

    bool isImageSet();
    bool isMeshSet();
    bool isParametersImported();
    bool areResultsLoaded();
    bool isMeshMeasured();

private slots:
    //void insertOperation(const QString &operation);

private:

    QVTKWidget *sliceView;
    QVTKWidget *threeDimView;
    QVTKWidget *profileView;

    //vtkSmartPointer< SliceInteractorStyle > sliceStyle;
    //vtkSmartPointer< ThreeDimInteractorStyle > threeDimStyle;

    vtkSmartPointer<vtkRenderer> sliceRenderer;
    vtkSmartPointer<vtkRenderer> threeDimRenderer;
    vtkSmartPointer<vtkRenderer> profileRenderer;

    void createObjects();
    QDockWidget * createSliceView();
    QDockWidget * createThreeDimView();
    QDockWidget * createProfileView();
    void createStyles();

    void updateProfileView();
    void setupProfileView();



    //QVTKWidget *testView;
    //QDockWidget * createTestView();

    itk::CorticalBone::Pointer corticalBone;


};

#endif	/* VISUALISATIONFRAME_H */

