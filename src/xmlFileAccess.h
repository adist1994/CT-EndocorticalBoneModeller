/* 
 * File:   xmlFileAccess.h
 * Author: rap58
 *
 * Created on 04 April 2015, 14:01
 */

#ifndef XMLFILEACCESS_H
#define	XMLFILEACCESS_H

/* General notes
 * 1. read and write xml files
 * 2. accessed by MainWindow and ScriptFrame
 * 3. contains a snapshot of project information either:  
 *    a. read from an xml file (to be accessed by the MainWindow or ScriptFrame)
 *    b. set by the MainWindow or ScriptFile (to be written to an xml file)
 * 4. it reset values before reading an xml file or after writing an xml file
 * 
 * 
 *  */

// itk + namespace

#include "corticalbone.h"
#include <QObject>
#include <QtCore/qxmlstream.h>
//#include <QXmlStreamReader>
//#include <QXmlStreamWriter>

//------------------ class predefs------------------------//
class QXmlStreamWriter;
class QXmlStreamReader;
//class QString;

//------------------ type defs -------------------//

class xmlFileAccess : QObject {

    // no silent errors. Always raise a warning or error.
    /* TODO  
     * 1. 
     */

    /* Define API 
     * 1. Getters
     *  a. QT objects and primitives
     * 2. Setters
     *  a. Set file paths
     * 3. Open / Save
     *  a. 
     * */


public:


    xmlFileAccess(bool verboseIn=false); // constructor

    void reset();

    //-- setters

    // project details
    void setImageName(QString imagePath);
    void setMeshName(QString meshFile);
    void setImportBaseName(QString importBaseName);

    void setCBConstraintInfo(int cbConstraintTypeIn, double cbFixedValueIn = nan("1"));
    void setCBConstraintInfo(int cbConstraintTypeIn, QString cbFixedFilterIn = QString());
    void setCalibrationName(int index, QString name);
    void setCalibrationPoints(vtkSmartPointer<vtkDoubleArray> calPtsIn);
    void setCalibrationPoints(vtkSmartPointer<vtkDoubleArray> calPtsIn, vtkSmartPointer<vtkDoubleArray> calValsIn);
    void setCalibrationPtGeometry(double radius, int number);
    void setCalibrationParameters(double p0, double p1, double p2);
    void setCalibrationParameterFile(QString fileName);
    void setProfileAveraging(double average[Dimension]);
    void setSampleNumber(int sampleNumberIn);
    void setSmoothing(double smoothingValueIn);
    void setSigmaConstraint(double sigma[Dimension]);
    void setClassifierThresholdInfo(double st, double cb, double threshold, int percentIndex, QString classifierNameIn);
    void setClassifierThresholdInfo(double st, double cb, double threshold, int percentIndex, QString classifierNameIn, double weight);
    void setModelInfo(int modelIndexIn, QString modelNameIn, int schemeIndexIn, QString schemeNameIn, int optimiserIndexIn, QString optimiserNameIn);


    // script details
    void setScriptProjectInfo(QString projectPath, QString projectFilter);
    void setScriptImageInfo(QString imageFilter, QString imageExtnName, int imageExtnIndex);
    void setScriptMeshInfo(QString meshFilter, QString meshExtnName, int meshExtnIndex);
    void setClassifierThresholdInfo(int optimiserIndexIn, QString optimiserNameIn);
    void setClassifierThresholdInfo(int optimiserIndexIn, QString optimiserNameIn, double weight);
    void setScriptActiveTab(int tabIndex, QString tabName);

    //-- getters
    QString getProjectName();
    QString getImagePath();
    QString getMeshName();
    QString getImportBaseName();

    int getSampleNumber();
    double getSmoothing();
    void getCBConstraintInfo(int& cbConstraintTypeIn, double& cbFixedValueIn); // project
    void getCBConstraintInfo(int& cbConstraintTypeIn, QString& cbFixedFilterIn); // script
    bool getCalibrationName(int& index, QString& name);
    vtkSmartPointer<vtkDoubleArray> getCalibrationPoints();
    vtkSmartPointer<vtkDoubleArray> getCalibrationValues();
    bool getCalibrationPtGeometry(double &radius, int &number);
    bool getCalibrationParameters(double& p0, double& p1, double& p2);
    bool getCalibrationParameterFile(QString& name);
    bool getProfileAveraging(double average[Dimension]);
    bool getSigmaConstraint(double sigma[Dimension]);
    bool getClassifierThresholds(double &st, double &cb, double &threshold, int &classifierIndexIn);
    bool getClassifierThresholdWeight(double &weight);
    bool getModelInfo(int &modelIndexIn, QString &modelNameIn, int &schemeIndexIn, QString &schemeNameIn, int &optimiserIndexIn, QString &optimiserNameIn);

    bool isImageSet();
    bool isMeshSet();
    bool isImportSet();
    bool isSampleNumberSet();
    bool isSmoothingSet();
    bool isAveragingSet();
    bool areCalibrationParametersSet();
    bool isCalibrationFileSet();
    bool isSigmaSet();
    bool isClassifierThresholdIndexSet();
    bool areClassifierThresholdsSet();
    bool isClassifierThresholdWeightSet();
    bool isModelSet();

    // script details
    void getScriptProjectInfo(QString &scriptProjectPathIn, QString &scriptProjectFilterIn);
    void getScriptImageInfo(QString &imageFilter, QString &imageExtnName, int &imageExtnIndex);
    void getScriptMeshInfo(QString &meshFilter, QString &meshExtnName, int &imageExtnIndex);
    bool getClassifierThresholdIndex(int &thresholdIndexIn, QString &thresholdNameIn);
    void getScriptActiveTab(int &tabIndex, QString &tabName);

    // open
    bool openProjectFile(QString projectFileName);
    bool openScriptFile(QString filename);

    // save
    bool saveProjectFile(QString projectFileName, bool isScript);
    bool saveScriptFile(QString scriptFileName);

    enum{
        kCBConstraintNone=0,
        kCBConstraintFixed,
        kCBConstraintFWHM,
    } FA_CBConstraintType;

    enum{
        kCreationScript=0,
        kViewScript,
        kProcessScript,
    } FA_ScriptType;

private:

    void resetState();

    QString getCBConstraintName();

    // save
    bool saveProjectFileUser(QXmlStreamWriter &xmlWriter);
    bool saveProjectFileScript(QXmlStreamWriter &xmlWriter);

    // open
    bool openCreationScript(QXmlStreamReader &xmlReader);
    bool openViewScript(QXmlStreamReader &xmlReader);
    bool openProcessScript(QXmlStreamReader &xmlReader);

    bool readImagePath(QXmlStreamReader &xmlReader);
    bool readMeshPath(QXmlStreamReader &xmlReader);
    bool readImportPath(QXmlStreamReader &xmlReader);
    bool readSampleNumber(QXmlStreamReader &xmlReader);
    bool readSmoothing(QXmlStreamReader &xmlReader);
    bool readCalibration(QXmlStreamReader &xmlReader);
    bool readCBConstraint(QXmlStreamReader &xmlReader);
    bool readCBScriptConstraint(QXmlStreamReader &xmlReader);
    bool readProfileAveraging(QXmlStreamReader &xmlReader);
    bool readGlobalSigma(QXmlStreamReader &xmlReader);
    bool readClassifierInfo(QXmlStreamReader &xmlReader);
    bool readModelInfo(QXmlStreamReader &xmlReader);

    bool readScriptProjectInfo(QXmlStreamReader &xmlReader);
    bool readScriptImageInfo(QXmlStreamReader &xmlReader);
    bool readScriptMeshInfo(QXmlStreamReader &xmlReader);
    bool readScriptCalibrationScales(QXmlStreamReader &xmlReader);
    bool readScriptActiveTab(QXmlStreamReader &xmlReader);

    bool writeImagePath(QXmlStreamWriter &xmlWriter);
    bool writeMeshPath(QXmlStreamWriter &xmlWriter);
    bool writeImportPath(QXmlStreamWriter &xmlWriter);
    bool writeSampleNumber(QXmlStreamWriter &xmlWriter);
    bool writeSmoothing(QXmlStreamWriter &xmlWriter);
    bool writeCalibration(QXmlStreamWriter &xmlWriter);
    bool writeCBConstraint(QXmlStreamWriter &xmlWriter);
    bool writeCBScriptConstraint(QXmlStreamWriter &xmlWriter);
    bool writeProfileAveraging(QXmlStreamWriter &xmlWriter);
    bool writeGlobalSigma(QXmlStreamWriter &xmlWriter);
    bool writeClassifierInfo(QXmlStreamWriter &xmlWriter);
    bool writeModelInfo(QXmlStreamWriter &xmlWriter);

    bool writeScriptProjectInfo(QXmlStreamWriter &xmlWriter);
    bool writeScriptImageInfo(QXmlStreamWriter &xmlWriter);
    bool writeScriptMeshInfo(QXmlStreamWriter &xmlWriter);
    bool writeScriptCalibrationScales(QXmlStreamWriter &xmlWriter);
    bool writeScriptActiveTab(QXmlStreamWriter &xmlWriter);

    bool verbose;

    // project variables
    bool imageSet, meshSet;
    bool importSet;
    bool sampleNumberSet, smoothingSet, calibrationParametersSet, calibrationFileSet, averagingSet, cbConstraintSet, sigmaSet, classifierThresholdsSet, classifierThresholdIndexSet, classifierThresholdWeightSet, modelSet;

    QString imageFileName, meshFilePath, importBasePath;
    int sampleNumber;
    int cbConstraintIndex; double fixedCBValue; QString fixedCBFilter;
    double profileAverages[Dimension];
    QString calPhantomName; int calPhantomIndex;
    vtkSmartPointer<vtkDoubleArray> calPts;
    vtkSmartPointer<vtkDoubleArray> calVals;
    double calP0, calP1, calP2; double calRadius; int calNumber;
    QString calFile;
    double sigmaValues[Dimension];
    double stDensity, cbDensity, thresholdDensity;
    double classifierWeight;
    double smoothingValue;

    // script variables
    QString scriptTabName;
    QString scriptProjectPath, scriptProjectFilter;
    QString scriptImageFilter, scriptImageExtn;
    QString scriptMeshFilter, scriptMeshExtn;
    int scriptImageExtnIndex, scriptMeshExtnIndex, scriptTabIndex;

    // settings details
    int modelIndex, optimiserIndex, schemeIndex, classifierIndex;
    QString modelName, optimiserName, schemeName, classifierName;

};



#endif	/* XMLFILEACCESS_H */

