/* includes */
#include "xmlFileAccess.h"

#include <QtWidgets>
#include <iostream>
#include <corticalbone.h>


/* Public Methods */
xmlFileAccess::xmlFileAccess(bool verboseIn) {

    verbose = verboseIn;
    resetState();
}

void xmlFileAccess::reset() {
    resetState();
}

//-- setters [Project info]
void xmlFileAccess::setImageName(QString imagePath) {
    imageSet=true;
    imageFileName = imagePath;
}

void xmlFileAccess::setMeshName(QString meshPath) {
    meshSet=true;
    meshFilePath = meshPath;
}

void xmlFileAccess::setImportBaseName(QString importPath) {
    importSet=true;
    importBasePath = importPath;
}

void xmlFileAccess::setSigmaConstraint(double sigma[Dimension]) {
    sigmaSet=true;
    sigmaValues[0] = sigma[0];
    sigmaValues[1] = sigma[1];
    sigmaValues[2] = sigma[2];
}

void xmlFileAccess::setClassifierThresholdInfo(double st, double cb, double threshold, int classifierIndexIn,
                                               QString classifierNameIn) { // for project

    if(classifierIndexIn==itk::ClassifierTransform::kMedianManual) {
        cerr<<"Error in xmlFileAccess::setClassifierThresholdInfo a threshold mode kMedianManual used without a threshold weight value"<<endl;
        classifierThresholdIndexSet = classifierThresholdsSet = classifierThresholdWeightSet = false;
    } else {
        classifierThresholdsSet = classifierThresholdIndexSet = true; classifierThresholdWeightSet = false;
        stDensity = st;
        cbDensity = cb;
        thresholdDensity = threshold;
        classifierIndex = classifierIndexIn;
        classifierName = classifierNameIn;
    }
}

void xmlFileAccess::setClassifierThresholdInfo(double st, double cb, double threshold, int classifierIndexIn,
                                               QString classifierNameIn, double weight) { // for project
    if(classifierIndexIn==itk::ClassifierTransform::kMedianManual) {
        classifierThresholdIndexSet = classifierThresholdsSet = classifierThresholdWeightSet = true;
        stDensity = st; cbDensity = cb; thresholdDensity = threshold;
        classifierIndex = classifierIndexIn; classifierName = classifierNameIn;
        classifierWeight = weight;
    } else {
        cerr<<"Error in xmlFileAccess::setClassifierThresholdInfo a threshold mode other than kMedianManual associated with a threshold weight value"<<endl;
        classifierThresholdIndexSet = classifierThresholdsSet = classifierThresholdWeightSet = false;
    }
}

void xmlFileAccess::setModelInfo(int modelIndexIn, QString modelNameIn, int schemeIndexIn, QString schemeNameIn, int optimiserIndexIn, QString optimiserNameIn) {
    modelIndex = modelIndexIn; modelName = modelNameIn;
    schemeIndex = schemeIndexIn; schemeName = schemeNameIn;
    optimiserIndex = optimiserIndexIn; optimiserName=optimiserNameIn;
    modelSet = true;

}

void xmlFileAccess::setCBConstraintInfo(int cbConstraintTypeIn, double cbFixedValueIn) {
    cbConstraintSet = true;
    cbConstraintIndex = cbConstraintTypeIn;
    fixedCBValue = cbFixedValueIn;
}

void xmlFileAccess::setCBConstraintInfo(int cbConstraintTypeIn, QString cbFixedFilterIn) {
    cbConstraintSet = true;
    cbConstraintIndex = cbConstraintTypeIn;
    fixedCBFilter = cbFixedFilterIn;
}

void xmlFileAccess::setCalibrationName(int index, QString name) {
    calPhantomIndex = index;
    calPhantomName = name;
    if(calPhantomIndex==itk::CorticalBone::kFileSpecifiedCal) {
        calibrationFileSet = true;
    } else {
        calibrationParametersSet = true;
    }
}

void xmlFileAccess::setCalibrationPoints(vtkSmartPointer<vtkDoubleArray> calPtsIn) {
    calibrationParametersSet = true;
    calPts=calPtsIn;
}

void xmlFileAccess::setCalibrationPoints(vtkSmartPointer<vtkDoubleArray> calPtsIn, vtkSmartPointer<vtkDoubleArray> calValsIn) {
    calibrationParametersSet = true;
    calPts=calPtsIn; calVals=calValsIn;
}

void xmlFileAccess::setCalibrationPtGeometry(double radius, int number) {
    calibrationParametersSet = true;
    calRadius = radius; calNumber = number;
}

void xmlFileAccess::setCalibrationParameters(double p0, double p1, double p2) {
    calibrationParametersSet = true;
    calP0 = p0;
    calP1 = p1;
    calP2 = p2;
}

void xmlFileAccess::setCalibrationParameterFile(QString fileName) {
    calibrationFileSet = true;
    calFile = fileName;
}

void xmlFileAccess::setProfileAveraging(double average[Dimension]) {
    averagingSet=true;
    profileAverages[0] = average[0]; profileAverages[1] = average[1]; profileAverages[2] = average[2];

}

void xmlFileAccess::setSampleNumber(int sampleNumberIn) {
    sampleNumberSet = true;
    sampleNumber = sampleNumberIn;
}

void xmlFileAccess::setSmoothing(double smoothingValueIn) {
    smoothingSet = true;
    smoothingValue = smoothingValueIn;
}

//-- setters [Scripting Info]
void xmlFileAccess::setScriptProjectInfo(QString projectPath, QString projectFilter) {
    scriptProjectPath = projectPath;
    scriptProjectFilter = projectFilter;

}

void xmlFileAccess::setScriptImageInfo(QString imageFilter, QString imageExtnName, int imageExtnIndex) {
    scriptImageFilter = imageFilter;
    scriptImageExtn = imageExtnName;
    scriptImageExtnIndex = imageExtnIndex;

}

void xmlFileAccess::setScriptMeshInfo(QString meshFilter, QString meshExtnName, int meshExtnIndex) {
    scriptMeshFilter = meshFilter;
    scriptMeshExtn = meshExtnName;
    scriptMeshExtnIndex = meshExtnIndex;

}

void xmlFileAccess::setClassifierThresholdInfo(int classifierIndexIn, QString classifierNameIn) {
    if(classifierIndexIn==itk::ClassifierTransform::kMedianManual) {
        cerr<<"Error in xmlFileAccess::setClassifierThresholdInfo a threshold mode kMedianManual used without a threshold weight value"<<endl;
        classifierThresholdIndexSet = classifierThresholdsSet = classifierThresholdWeightSet = false;
    } else {
        classifierIndex = classifierIndexIn;
        classifierName = classifierNameIn;
        classifierThresholdIndexSet = true; classifierThresholdsSet = classifierThresholdWeightSet = false;
    }
}

void xmlFileAccess::setClassifierThresholdInfo(int classifierIndexIn, QString classifierNameIn, double weight) {
    if(classifierIndexIn==itk::ClassifierTransform::kMedianManual) {
        classifierIndex = classifierIndexIn;
        classifierName = classifierNameIn;
        classifierWeight = weight;
        classifierThresholdIndexSet = classifierThresholdWeightSet = true; classifierThresholdsSet = false;
    } else {
        cerr<<"Error in xmlFileAccess::setClassifierThresholdInfo a threshold mode other than kMedianManual used with a threshold weight value"<<endl;
        classifierThresholdIndexSet = classifierThresholdsSet = classifierThresholdWeightSet = false;
    }
}

void xmlFileAccess::setScriptActiveTab(int tabIndex, QString tabName) {
    scriptTabName = tabName;
    scriptTabIndex = tabIndex;
}

//-- getters [Project Info]
QString xmlFileAccess::getImagePath() {
    return imageFileName;
}

QString xmlFileAccess::getMeshName() {
    return meshFilePath;
}

QString xmlFileAccess::getImportBaseName() {
    return importBasePath;
}

int xmlFileAccess::getSampleNumber() {
    return sampleNumber;
}

double xmlFileAccess::getSmoothing() {
    return smoothingValue;
}

void xmlFileAccess::getCBConstraintInfo(int& cbConstraintTypeIn, double& cbFixedValueIn) {
    cbConstraintTypeIn = cbConstraintIndex;
    if(cbConstraintIndex == kCBConstraintFixed) {
        cbFixedValueIn = fixedCBValue;
    } else {
        cbFixedValueIn = nan("1");
    }
}

void xmlFileAccess::getCBConstraintInfo(int& cbConstraintTypeIn, QString& cbFixedFilterIn) {
    cbConstraintTypeIn = cbConstraintIndex;
    if(cbConstraintIndex == kCBConstraintFixed) {
        cbFixedFilterIn = fixedCBFilter;
    } else {
        cbFixedFilterIn = QString();
    }
}

bool xmlFileAccess::getCalibrationName(int& index, QString& name) {
    index = calPhantomIndex;
    name = calPhantomName;
    return calibrationFileSet || calibrationParametersSet;
}

vtkSmartPointer<vtkDoubleArray> xmlFileAccess::getCalibrationPoints() {
    if(calibrationParametersSet) {
        return calPts;
    } else {
        vtkSmartPointer<vtkDoubleArray> pts = vtkSmartPointer<vtkDoubleArray>::New();
        return pts;
    }
}

vtkSmartPointer<vtkDoubleArray> xmlFileAccess::getCalibrationValues() {
    if(calibrationParametersSet && calPhantomIndex == itk::CorticalBone::kManualControlPtsCal) {
        return calVals;
    } else {
        vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();
        return vals;
    }
}

bool xmlFileAccess::getCalibrationParameters(double& p0, double& p1, double& p2) {
    p0 = calP0;
    p1 = calP1;
    p2 = calP2;
    return calibrationParametersSet;
}

bool xmlFileAccess::getCalibrationParameterFile(QString& name) {
    name = calFile;
    return calibrationFileSet;
}

bool xmlFileAccess::getCalibrationPtGeometry(double &radius, int &number) {
    radius = calRadius; number = calNumber;
    return calibrationParametersSet && calPhantomIndex == itk::CorticalBone::kManualControlPtsCal;
}

bool xmlFileAccess::getProfileAveraging(double average[Dimension]) {
    average[0] = profileAverages[0]; average[1] = profileAverages[1]; average[2]=profileAverages[2];
    return averagingSet;
}

bool xmlFileAccess::getSigmaConstraint(double sigma[Dimension]) {

    sigma[0] = sigmaValues[0]; sigma[1] = sigmaValues[1]; sigma[2] = sigmaValues[2];

    return sigmaSet;
}

bool xmlFileAccess::getClassifierThresholds(double &st, double &cb, double &threshold, int &classifierIndexIn) {
    st = stDensity; cb =cbDensity;
    threshold = thresholdDensity; classifierIndexIn = classifierIndex;
    return classifierThresholdsSet;
}

bool xmlFileAccess::getClassifierThresholdWeight(double &weight) {
    weight = classifierWeight;
    return classifierThresholdWeightSet;
}

bool xmlFileAccess::getModelInfo(int &modelIndexIn, QString &modelNameIn, int &schemeIndexIn, QString &schemeNameIn, int &optimiserIndexIn, QString &optimiserNameIn) {
    modelIndexIn = modelIndex; modelNameIn = modelName;
    schemeIndexIn = schemeIndex; schemeNameIn = schemeName;
    optimiserIndexIn = optimiserIndex; optimiserNameIn = optimiserName;
    return modelSet;
}


// Getters - states
bool xmlFileAccess::isImageSet() {
    return imageSet;
}

bool xmlFileAccess::isMeshSet() {
    return meshSet;
}

bool xmlFileAccess::isImportSet() {
    return importSet;
}

bool xmlFileAccess::isSampleNumberSet() {
    return sampleNumberSet;
}

bool xmlFileAccess::isSmoothingSet() {
    return smoothingSet;
}

bool xmlFileAccess::isAveragingSet() {
    return averagingSet;
}

bool xmlFileAccess::areCalibrationParametersSet() {
    return calibrationParametersSet;
}

bool xmlFileAccess::isCalibrationFileSet() {
    return calibrationFileSet;
}

bool xmlFileAccess::isSigmaSet() {
    return sigmaSet;
}

bool xmlFileAccess::isClassifierThresholdIndexSet() {
    return classifierThresholdIndexSet;
}

bool xmlFileAccess::areClassifierThresholdsSet() {
    return classifierThresholdsSet;
}

bool xmlFileAccess::isClassifierThresholdWeightSet() {
    return classifierThresholdWeightSet;
}

bool xmlFileAccess::isModelSet() {
    return modelSet;
}

//-- getters [Scripting Info]
void xmlFileAccess::getScriptProjectInfo(QString &scriptProjectPathIn, QString &scriptProjectFilterIn) {
    scriptProjectPathIn = scriptProjectPath;
    scriptProjectFilterIn = scriptProjectFilter;
}

void xmlFileAccess::getScriptImageInfo(QString &imageFilter, QString &imageExtnName, int &imageExtnIndex) {
    imageFilter = scriptImageFilter;
    imageExtnName = scriptImageExtn;
    imageExtnIndex = scriptImageExtnIndex;
}

void xmlFileAccess::getScriptMeshInfo(QString &meshFilter, QString &meshExtnName, int &imageExtnIndex) {
    meshFilter = scriptMeshFilter;
    meshExtnName = scriptMeshExtn;
    imageExtnIndex = scriptMeshExtnIndex;
}

bool xmlFileAccess::getClassifierThresholdIndex(int &classifierthresholdIndexIn, QString &classifierNameIn) {
    classifierthresholdIndexIn = classifierIndex;
    classifierNameIn = classifierName;
    return classifierThresholdIndexSet;
}

void xmlFileAccess::getScriptActiveTab(int &tabIndex, QString &tabName) {
    tabIndex = scriptTabIndex; tabName = scriptTabName;
}

// open files
bool xmlFileAccess::openProjectFile(QString projectFileName) {

    resetState();
    if(verbose){std::cout<<"<b>----- Begin Parsing Project File -----</b>"<<std::endl;}
    QFile* file = new QFile(projectFileName);

    //---- check file can be written -----//
    if (!file->open(QFile::ReadOnly)) {
        QString msg = QString("Failed to open %1\n%2").arg(projectFileName).arg(file->errorString());
        //QMessageBox::warning(this, QString("Error"), msg);
        if(verbose){std::cout<<msg.toStdString()<<std::endl;}
        return false;
    }

    QXmlStreamReader xmlReader;
    xmlReader.setDevice(file);
    xmlReader.readNextStartElement();

    // TODO - make it so parsing is completed then values are set
    bool status = true;

    if(verbose){std::cout<<Utilities::getTabString()<<"file path: "<<projectFileName.toStdString()<<std::endl;}

    while(!xmlReader.atEnd() && !xmlReader.hasError() && status == true) {
        /* Read next element.*/
        xmlReader.readNextStartElement();

        if(xmlReader.name() == "CBM_ProjectFile" && xmlReader.isEndElement()) { // TODO perform start element + "CBMFileType" checks before the while loop
            // do nothing
            break;
        } else if(xmlReader.name() == "ImageFilePath" && xmlReader.isStartElement()) { // read in image file
            status = readImagePath(xmlReader) && status;
        } else if(xmlReader.name() == "MeshFilePath" && xmlReader.isStartElement()) { // read in image file
            status = readMeshPath(xmlReader) && status;
        } else if(xmlReader.name() == "ModelInformation" && xmlReader.isStartElement()) {
            status = readModelInfo(xmlReader) && status;
        } else if(xmlReader.name() == "ImportParameters" && xmlReader.isStartElement()) { // read in image file
            status = readImportPath(xmlReader) && status;
        } else if(xmlReader.name() == "FixedSampleNumber" && xmlReader.isStartElement()) {
            status = readSampleNumber(xmlReader) && status;
        } else if(xmlReader.name() == "Calibration" && xmlReader.isStartElement()) { // dig further and get cal info
            status = readCalibration(xmlReader) && status;
        } else if(xmlReader.name() == "CBConstraint" && xmlReader.isStartElement()) {
            status = readCBConstraint(xmlReader) && status;
        } else if(xmlReader.name() == "ProfileAveraging" && xmlReader.isStartElement()) {
            status = readProfileAveraging(xmlReader) && status;
        } else if(xmlReader.name() == "GlobalSigma" && xmlReader.isStartElement()) {
            readGlobalSigma(xmlReader);
        } else if(xmlReader.name() == "SmoothingValue" && xmlReader.isStartElement()) {
            status = readSmoothing(xmlReader) && status;
        } else if(xmlReader.name() == "ClassifierInformation" && xmlReader.isStartElement()) {
            readClassifierInfo(xmlReader);
        } else if ( xmlReader.isStartElement()) { // set display settings, specify model settings
            if(verbose){std::cout<<"\tUnsupported command ignored: "<<xmlReader.name().toString().toStdString()<<std::endl;}
        }
    }

    // handle any xml errors
    if(xmlReader.hasError()) {
        //QString msg = QString("File: %1\n%2").arg(projectFileName).arg(xmlReader.errorString());
        //QMessageBox::critical(this, "Error reading in project file: ", msg, QMessageBox::Ok);
        std::cerr<<xmlReader.lineNumber()<<Utilities::getTabString()<<"\tWarning: "<<xmlReader.errorString().toStdString()<<std::endl;
        status = false;
    }
    if(verbose){std::cout<<"<b>----- Finish Parsing Project File -----</b>"<<std::endl;}
    return status;
}

bool xmlFileAccess::openScriptFile(QString filename) {
    resetState();

    QFile* file = new QFile(filename);

    //---- check file can be written -----//
    if (!file->open(QFile::ReadOnly)) {
        QString msg = QString("Failed to open %1\n%2").arg(filename).arg(file->errorString());
        //QMessageBox::warning(this, QString("Error"), msg);
        if(verbose){std::cout<<msg.toStdString()<<std::endl;}
        return false;
    }

    QXmlStreamReader xmlReader;
    xmlReader.setDevice(file);
    xmlReader.readNext();

    // read start document
    while(!xmlReader.isStartElement()) {
        xmlReader.readNext();
    }

    if(xmlReader.name() != "CBM_Script") {
        std::clog<<"Invalid script file. Possibly of an old now obsolete formate. Invalid header name = "<<xmlReader.name().toString().toStdString()<<std::endl;
        return -1;
    }

    std::clog<<Utilities::getTabString()<<"script file: "<<filename.toStdString()<<std::endl;

    while(!xmlReader.atEnd()) {
        xmlReader.readNextStartElement();

        if(xmlReader.name() == "ProjectInformation" && xmlReader.isStartElement()) {
            readScriptProjectInfo(xmlReader);
        } else if(xmlReader.name() == "ActiveTab" && xmlReader.isStartElement()) {
            readScriptActiveTab(xmlReader);
        } else if(xmlReader.name() == "Creation_Tab") {
            openCreationScript(xmlReader);
        } else if(xmlReader.name() == "View_Tab") {
            openViewScript(xmlReader);
        } else if(xmlReader.name() == "Process_Tab") {
            openProcessScript(xmlReader);
        } else if(xmlReader.isStartElement()) {
            std::clog<<Utilities::getTabString()<<Utilities::getTabString()<<"Warning unrecognised entry '"<<xmlReader.name().toString().toStdString()<<"' in Script"<<std::endl;
        }
    }

    return true;


}

bool xmlFileAccess::openCreationScript(QXmlStreamReader &xmlReader) {

    while(!(xmlReader.name() == "Creation_Tab" && xmlReader.isEndElement())) {
        xmlReader.readNextStartElement();
        QString tokenName = xmlReader.name().toString();
        if(tokenName == "ImageInformation" && xmlReader.isStartElement()) {
            readScriptImageInfo(xmlReader);
        } else if(tokenName == "MeshInformation" && xmlReader.isStartElement()) {
            readScriptMeshInfo(xmlReader);
        } else if(tokenName == "FixedSampleNumber" && xmlReader.isStartElement()) {
            readSampleNumber(xmlReader);
        } else if(tokenName == "CBConstraint" && xmlReader.isStartElement()) {
            readCBScriptConstraint(xmlReader);
        } else if(tokenName == "GlobalSigma" && xmlReader.isStartElement()) {
            readGlobalSigma(xmlReader);
        } else if(tokenName == "ProfileAveraging" && xmlReader.isStartElement()) {
            readProfileAveraging(xmlReader);
        } else if(tokenName == "GlobalSigma" && xmlReader.isStartElement()) {
            readGlobalSigma(xmlReader);
        } else if(tokenName == "ImportParameters" && xmlReader.isStartElement()) {
            readImportPath(xmlReader);
        } else if(xmlReader.isStartElement()) {
            std::clog<<Utilities::getTabString()<<Utilities::getTabString()<<"Warning unrecognised entry '"<<tokenName.toStdString()<<"' in Creation Section of Script"<<std::endl;
        }
    }
    return true;
}

bool xmlFileAccess::openViewScript(QXmlStreamReader &xmlReader) {

    while(!(xmlReader.name() == "View_Tab" && xmlReader.isEndElement())) {
        xmlReader.readNextStartElement();
        QString tokenName = xmlReader.name().toString();
        if(tokenName == "Calibration" && xmlReader.isStartElement()) {
            readScriptCalibrationScales(xmlReader);
        } else if(xmlReader.isStartElement()) {
            std::clog<<Utilities::getTabString()<<Utilities::getTabString()<<"Warning unrecognised entry '"<<tokenName.toStdString()<<"' in View Section of Script"<<std::endl;
        }
    }
    return true;
}

bool xmlFileAccess::openProcessScript(QXmlStreamReader &xmlReader) {

    while(!(xmlReader.name() == "Process_Tab" && xmlReader.isEndElement())) {
        xmlReader.readNextStartElement();
        QString tokenName = xmlReader.name().toString();
        if(tokenName == "ModelInformation" && xmlReader.isStartElement()) {
            readModelInfo(xmlReader);
        } if(tokenName == "ClassifierInformation" && xmlReader.isStartElement()) {
            readClassifierInfo(xmlReader);
        } else if(tokenName == "SmoothingValue" && xmlReader.isStartElement()) {
            readSmoothing(xmlReader);
        } else if(xmlReader.isStartElement()) {
            std::clog<<Utilities::getTabString()<<Utilities::getTabString()<<"Warning unrecognised entry '"<<tokenName.toStdString()<<"' in Process Section of Script"<<std::endl;
        }
    }
    return true;
}


// save files
bool xmlFileAccess::saveProjectFile(QString projectFileName, bool isScript) {

    QFile* file = new QFile(projectFileName);

    //---- check file can be written -----//
    if (!file->open(QFile::WriteOnly)) {
        QString msg = QString("Failed to open %1\n%2").arg(projectFileName).arg(file->errorString());
        //QMessageBox::warning(this, QString("Error"), msg);
        if(verbose){cout<<msg.toStdString()<<endl;}
        cerr<<msg.toStdString()<<endl;
        return false;
    }

    //---- open xml streamer ------//
    QXmlStreamWriter xmlWriter(file);
    xmlWriter.setAutoFormatting(true);
    xmlWriter.writeStartDocument();

    // write file body
    xmlWriter.writeStartElement("CBM_ProjectFile");
    if(isScript) {
        saveProjectFileScript(xmlWriter);
    } else {
        saveProjectFileUser(xmlWriter);

    }

    /***********************
     * TODO
     * 1. specify model properties: optimiser type, model type
     * 2. display options: mesh/volume, point measurement mode
     *  */

    xmlWriter.writeEndElement();
    xmlWriter.writeEndDocument();

    //---- close streamer ----//
    file->close();

    resetState();


    return true;

}

bool xmlFileAccess::saveScriptFile(QString scriptFileName) {
    std::clog<<Utilities::getTabString()<<"Saving Script file:"<<scriptFileName.toStdString().c_str()<<endl;

    // write stored values to file
    QFile* file = new QFile(scriptFileName);

    //---- check file can be written -----//
    if (!file->open(QFile::WriteOnly)) {
        QString msg = QString("Failed to open %1\n%2").arg(scriptFileName).arg(file->errorString());
        //QMessageBox::warning(this, QString("Error"), msg);
        if(verbose){cout<<msg.toStdString()<<endl;}
        return false;
    }

    //---- open xml streamer ------//
    QXmlStreamWriter xmlWriter(file);
    xmlWriter.setAutoFormatting(true);
    xmlWriter.writeStartDocument();

    // write file body
    xmlWriter.writeStartElement("CBM_Script");

    // general info
    writeScriptProjectInfo(xmlWriter);
    writeScriptActiveTab(xmlWriter);

    // creation info
    xmlWriter.writeStartElement("Creation_Tab");
    writeScriptImageInfo(xmlWriter);
    writeScriptMeshInfo(xmlWriter);
    if(importSet) {
        writeImportPath(xmlWriter);
    }
    if(sampleNumberSet) {
        writeSampleNumber(xmlWriter);
    }
    // constraints
    writeCBScriptConstraint(xmlWriter); // regardless of type write info
    if(sigmaSet) {
        writeGlobalSigma(xmlWriter);
    }
    if(averagingSet) {
        writeProfileAveraging(xmlWriter);
    }
    xmlWriter.writeEndElement();

    // view info
    xmlWriter.writeStartElement("View_Tab");
    if(calibrationParametersSet || calibrationFileSet) {
        writeScriptCalibrationScales(xmlWriter);
    }
    xmlWriter.writeEndElement();

    // process info
    xmlWriter.writeStartElement("Process_Tab");
    writeModelInfo(xmlWriter);
    if(classifierThresholdIndexSet) {
        writeClassifierInfo(xmlWriter);
    }
    if(smoothingSet) {
        writeSmoothing(xmlWriter);
    }
    xmlWriter.writeEndElement();


    xmlWriter.writeEndDocument();

    // close file stream and reset state
    file->close();
    resetState();
    
    return true;
}

/* Private Methods */
void xmlFileAccess::resetState() {

    // project variables
    imageSet= meshSet= importSet= sampleNumberSet= smoothingSet= calibrationParametersSet = calibrationFileSet = averagingSet=false;
    sigmaSet = classifierThresholdsSet = modelSet = false;
    classifierThresholdIndexSet = classifierThresholdWeightSet = false;
    imageFileName = meshFilePath = importBasePath = QString();

    cbConstraintSet = false;
    cbConstraintIndex = 0; // defaults to none

    sampleNumber = -1;
    smoothingValue = nan("1");
    fixedCBValue = nan("0"); fixedCBFilter = QString();

    profileAverages[0] = profileAverages[1] = profileAverages[2] = nan("0");

    calPhantomName = QString(); calPhantomIndex = -1;
    calPts = vtkSmartPointer<vtkDoubleArray>::New();
    calVals = vtkSmartPointer<vtkDoubleArray>::New();

    calP1 = calP0 = calP2 = nan("0"); calRadius=nan("1"); calNumber=-1;

    // script variables
    scriptProjectPath = QString(); scriptProjectFilter = QString();
    scriptTabIndex = -1; scriptTabName = QString();
    scriptImageFilter = scriptImageExtn = QString(); scriptImageExtnIndex = -1;
    scriptMeshFilter = scriptMeshExtn = QString(); scriptMeshExtnIndex = -1;
    modelIndex = -1; modelName = QString();
    optimiserIndex = -1; optimiserName = QString();
    schemeIndex = -1; schemeName = QString();
    classifierIndex = -1; classifierName = QString();

}

QString xmlFileAccess::getCBConstraintName() {
    QString typeName;
    if(cbConstraintIndex==kCBConstraintNone) {
        typeName=tr("None");
    } else if(cbConstraintIndex==kCBConstraintFixed) {
        typeName=tr("Fixed");
    } else if(cbConstraintIndex==kCBConstraintFWHM) {
        typeName=tr("FWHM");
    } else {
        if(verbose){cout<<"Error in xmlFileAccess::writeScriptCBConstraint. Invalid CB constraint index."<<endl;}
    }

    return typeName;
}

bool xmlFileAccess::saveProjectFileUser(QXmlStreamWriter &xmlWriter) {

    //---- write attributes ----//
    // image file name
    if(imageSet) {
        writeImagePath(xmlWriter);
    }
    // mesh file name
    if(meshSet) {
        writeMeshPath(xmlWriter);
    }
    // model info set
    if(modelSet) {
        writeModelInfo(xmlWriter);
    }
    // fixed sample number
    if(sampleNumberSet) {
        writeSampleNumber(xmlWriter);
    }
    // imported parameters
    if(importSet) {
        writeImportPath(xmlWriter);
    }
    // calibration (positions, type, values)
    if(calibrationParametersSet) {
        writeCalibration(xmlWriter);
    }
    // in case of fixed CB or set to FWHM mode - don't both if set to None as this is the default
    if(cbConstraintSet) {
        writeCBConstraint(xmlWriter);
    }
    // in case of HR profile averaging
    if(averagingSet) {
        writeProfileAveraging(xmlWriter);
    }
    // sigma set
    if(sigmaSet) {
        writeGlobalSigma(xmlWriter);
    }
    // smoothing set
    if(smoothingSet) {
        writeSmoothing(xmlWriter);
    }
    // HR thresholds set
    if(classifierThresholdIndexSet) {
        writeClassifierInfo(xmlWriter);
    }
    return true;

}

bool xmlFileAccess::saveProjectFileScript(QXmlStreamWriter &xmlWriter) {


    //---- write attributes ----//
    writeImagePath(xmlWriter);
    writeMeshPath(xmlWriter);

    if(modelSet) {
        writeModelInfo(xmlWriter);
    }

    if(sampleNumberSet) {
        writeSampleNumber(xmlWriter);
    }

    if(cbConstraintSet) {
        writeCBConstraint(xmlWriter);
    }

    if(averagingSet) {
        writeProfileAveraging(xmlWriter);
    }

    if(sigmaSet) {
        writeGlobalSigma(xmlWriter);
    }

    if(smoothingSet) {
        writeSmoothing(xmlWriter);
    }

    if(importSet) {
        writeImportPath(xmlWriter);
    }

    if(calibrationParametersSet) {
        writeCalibration(xmlWriter);
    }

    if(classifierThresholdIndexSet) {
        writeClassifierInfo(xmlWriter);
    }

    /***********************
     * TODO
     * 1. specify model properties: optimiser type, model type
     * 2. display options: mesh/volume, point measurement mode
     *  */

    return true;

}

bool xmlFileAccess::writeImagePath(QXmlStreamWriter &xmlWriter) {
    xmlWriter.writeTextElement("ImageFilePath", imageFileName );
    return true;
}

bool xmlFileAccess::writeMeshPath(QXmlStreamWriter &xmlWriter) {
    xmlWriter.writeTextElement("MeshFilePath", meshFilePath );
    return true;
}

bool xmlFileAccess::writeImportPath(QXmlStreamWriter &xmlWriter) {
    xmlWriter.writeStartElement("ImportParameters");

    xmlWriter.writeTextElement("Path", importBasePath);
    xmlWriter.writeEndElement();
    return true;
}

bool xmlFileAccess::writeSampleNumber(QXmlStreamWriter &xmlWriter) {
    xmlWriter.writeTextElement("FixedSampleNumber", QString::number(sampleNumber) );
    return true;
}

bool xmlFileAccess::writeSmoothing(QXmlStreamWriter &xmlWriter) {
    xmlWriter.writeTextElement("SmoothingValue", QString::number(smoothingValue) );
    return true;
}

bool xmlFileAccess::writeCalibration(QXmlStreamWriter &xmlWriter) {
    xmlWriter.writeStartElement("Calibration");

    xmlWriter.writeStartElement("Phantom");
    xmlWriter.writeAttribute("Name", calPhantomName);
    xmlWriter.writeCharacters(QString::number(calPhantomIndex));
    xmlWriter.writeEndElement();

    if(calPhantomIndex>=itk::CorticalBone::kMindwaySolidCal && calPhantomIndex<=itk::CorticalBone::kManualControlPtsCal) {
        xmlWriter.writeTextElement("radius", QString::number(calRadius));
        xmlWriter.writeTextElement("number", QString::number(calNumber));
    }

    vtkIdType numberOfPts=calPts->GetNumberOfTuples(); // todo ensure give a valid response for zero tuples

    for(int i=0; i<numberOfPts; i++) {
        double pt[3]; calPts->GetTuple(i, pt);
        xmlWriter.writeStartElement("Point"); xmlWriter.writeAttribute("Number",  QString::number(i));
        xmlWriter.writeTextElement("x", QString::number(pt[0]) );
        xmlWriter.writeTextElement("y", QString::number(pt[1]) );
        xmlWriter.writeTextElement("z", QString::number(pt[2]) );
        if(calPhantomIndex==itk::CorticalBone::kManualControlPtsCal) {
            xmlWriter.writeTextElement("value", QString::number(calVals->GetValue(i)) );
        }
        xmlWriter.writeEndElement();

    }

    if(calPhantomIndex!=itk::CorticalBone::kManualControlPtsCal) {
        xmlWriter.writeTextElement("p0", QString::number(calP0, 'e', 9));
        xmlWriter.writeTextElement("p1", QString::number(calP1, 'e', 9));
    }
    if(calPhantomIndex==itk::CorticalBone::kManualQuadraticCal) {
        xmlWriter.writeTextElement("p2", QString::number(calP2, 'e', 9));
    }

    xmlWriter.writeEndElement();
    return true;
}

bool xmlFileAccess::writeCBConstraint(QXmlStreamWriter &xmlWriter) {
    xmlWriter.writeStartElement("CBConstraint");

    QString typeName = getCBConstraintName();

    xmlWriter.writeStartElement("Type");
    xmlWriter.writeAttribute("Name", typeName);
    xmlWriter.writeCharacters(QString::number(cbConstraintIndex));
    xmlWriter.writeEndElement();

    if(cbConstraintIndex==kCBConstraintFixed) {
        xmlWriter.writeTextElement("Value", QString::number(fixedCBValue, 'g', 10));
    }

    xmlWriter.writeEndElement();
    return true;
}

bool xmlFileAccess::writeCBScriptConstraint(QXmlStreamWriter &xmlWriter) {
    xmlWriter.writeStartElement("CBConstraint");

    QString typeName = getCBConstraintName();

    xmlWriter.writeStartElement("Type");
    xmlWriter.writeAttribute("Name", typeName);
    xmlWriter.writeCharacters(QString::number(cbConstraintIndex));
    xmlWriter.writeEndElement();

    if(cbConstraintIndex==kCBConstraintFixed) {
        xmlWriter.writeTextElement("ValueFilter", fixedCBFilter);
    }

    xmlWriter.writeEndElement();
    return true;
}

bool xmlFileAccess::writeProfileAveraging(QXmlStreamWriter &xmlWriter) {

    xmlWriter.writeStartElement("ProfileAveraging");
    xmlWriter.writeTextElement("x", QString::number(profileAverages[0]) );
    xmlWriter.writeTextElement("y", QString::number(profileAverages[1]) );
    xmlWriter.writeTextElement("z", QString::number(profileAverages[2]) );
    xmlWriter.writeEndElement();

    return true;
}

bool xmlFileAccess::writeGlobalSigma(QXmlStreamWriter &xmlWriter) {

    xmlWriter.writeStartElement("GlobalSigma");
    xmlWriter.writeTextElement("x", QString::number(sigmaValues[0]) );
    xmlWriter.writeTextElement("y", QString::number(sigmaValues[1]) );
    xmlWriter.writeTextElement("z", QString::number(sigmaValues[2]) );
    xmlWriter.writeEndElement();

    return true;
}

bool xmlFileAccess::writeClassifierInfo(QXmlStreamWriter &xmlWriter) {

    xmlWriter.writeStartElement("ClassifierInformation");

    xmlWriter.writeStartElement("Selection");
    xmlWriter.writeAttribute("Name", classifierName);
    xmlWriter.writeCharacters(QString::number(classifierIndex));
    xmlWriter.writeEndElement();

    if(classifierThresholdWeightSet) {
        xmlWriter.writeTextElement("Weight", QString::number(classifierWeight));
    }

    if (classifierThresholdsSet) {
        xmlWriter.writeStartElement("Thresholds");
        xmlWriter.writeTextElement("ST", QString::number(stDensity));
        xmlWriter.writeTextElement("CB", QString::number(cbDensity));
        xmlWriter.writeTextElement("Threshold", QString::number(thresholdDensity));
        xmlWriter.writeEndElement();
    }
    xmlWriter.writeEndElement();

    return true;
}

bool xmlFileAccess::writeScriptProjectInfo(QXmlStreamWriter &xmlWriter) {

    xmlWriter.writeStartElement("ProjectInformation");

    xmlWriter.writeTextElement("Path", scriptProjectPath );
    xmlWriter.writeTextElement("Filter", scriptProjectFilter );

    xmlWriter.writeEndElement();

    return true;
}

bool xmlFileAccess::writeScriptImageInfo(QXmlStreamWriter &xmlWriter) {

    xmlWriter.writeStartElement("ImageInformation");

    xmlWriter.writeTextElement("Filter", scriptImageFilter);

    xmlWriter.writeStartElement("Extension");
    xmlWriter.writeAttribute("Name", scriptImageExtn);
    xmlWriter.writeCharacters(QString::number(scriptImageExtnIndex));
    xmlWriter.writeEndElement();

    xmlWriter.writeEndElement();

    return true;

}

bool xmlFileAccess::writeScriptMeshInfo(QXmlStreamWriter &xmlWriter) {

    xmlWriter.writeStartElement("MeshInformation");

    xmlWriter.writeTextElement("Filter", scriptMeshFilter);

    xmlWriter.writeStartElement("Extension");
    xmlWriter.writeAttribute("Name", scriptMeshExtn);
    xmlWriter.writeCharacters(QString::number(scriptMeshExtnIndex));
    xmlWriter.writeEndElement();

    xmlWriter.writeEndElement();

    return true;
}

bool xmlFileAccess::writeScriptCalibrationScales(QXmlStreamWriter &xmlWriter) {
    xmlWriter.writeStartElement("Calibration");

    // write phantom type
    xmlWriter.writeStartElement("Phantom");
    xmlWriter.writeAttribute("Name", calPhantomName);
    xmlWriter.writeCharacters(QString::number(calPhantomIndex));
    xmlWriter.writeEndElement();

    if(calPhantomIndex==itk::CorticalBone::kManualControlPtsCal) {
        xmlWriter.writeTextElement("radius", QString::number(calRadius));
        xmlWriter.writeTextElement("number", QString::number(calNumber));
    }

    if(calPhantomIndex==itk::CorticalBone::kFileSpecifiedCal) {
        xmlWriter.writeTextElement("ParameterFile", calFile);
    }

    // write calibration values
    if(calPhantomIndex==itk::CorticalBone::kManualLinearCal || calPhantomIndex==itk::CorticalBone::kManualQuadraticCal) {
        xmlWriter.writeTextElement("p0", QString::number(calP0));
        xmlWriter.writeTextElement("p1", QString::number(calP1));
    }
    if(calPhantomIndex==itk::CorticalBone::kManualQuadraticCal) {
        xmlWriter.writeTextElement("p2", QString::number(calP2));
    }

    xmlWriter.writeEndElement();
    return calibrationParametersSet;
}

bool xmlFileAccess::writeScriptActiveTab(QXmlStreamWriter &xmlWriter) {

    xmlWriter.writeStartElement("ActiveTab");

    xmlWriter.writeAttribute("Name", scriptTabName);
    xmlWriter.writeCharacters(QString::number(scriptTabIndex));

    xmlWriter.writeEndElement();

    return true;

}

bool xmlFileAccess::writeModelInfo(QXmlStreamWriter &xmlWriter) {
    xmlWriter.writeStartElement("ModelInformation");

    xmlWriter.writeStartElement("Model");
    xmlWriter.writeAttribute("Name", modelName);
    xmlWriter.writeCharacters(QString::number(modelIndex));
    xmlWriter.writeEndElement();

    xmlWriter.writeStartElement("Optimiser");
    xmlWriter.writeAttribute("Name", optimiserName);
    xmlWriter.writeCharacters(QString::number(optimiserIndex));
    xmlWriter.writeEndElement();

    xmlWriter.writeStartElement("Scheme");
    xmlWriter.writeAttribute("Name", schemeName);
    xmlWriter.writeCharacters(QString::number(schemeIndex));
    xmlWriter.writeEndElement();

    xmlWriter.writeEndElement();

    return true;
}


bool xmlFileAccess::readImagePath(QXmlStreamReader &xmlReader) {
    imageSet = true;
    imageFileName = xmlReader.readElementText();
    if(verbose){cout<<Utilities::getTabString()<<"Image Path: "<<imageFileName.toStdString()<<endl;}
    return true;
}

bool xmlFileAccess::readMeshPath(QXmlStreamReader &xmlReader) {
    meshSet = true;
    meshFilePath = xmlReader.readElementText();
    if(verbose){cout<<Utilities::getTabString()<<"Mesh Path: "<<meshFilePath.toStdString()<<endl;}
    return true;
}

bool xmlFileAccess::readImportPath(QXmlStreamReader &xmlReader) {

    xmlReader.readNextStartElement();
    importBasePath = xmlReader.readElementText();
    if(!importBasePath.isEmpty()) {
        importSet = true;
        if(verbose){cout<<Utilities::getTabString()<<"Import Path: "<<importBasePath.toStdString()<<endl;}
    } else {
        return false;
    }

    return true;

}

bool xmlFileAccess::readSampleNumber(QXmlStreamReader &xmlReader) {
    sampleNumberSet = true;
    sampleNumber = xmlReader.readElementText().toInt();
    if(verbose){cout<<Utilities::getTabString()<<"Sample Number: "<<sampleNumber<<endl;}
    return sampleNumberSet;
}

bool xmlFileAccess::readSmoothing(QXmlStreamReader &xmlReader) {
    smoothingSet = true;
    smoothingValue = xmlReader.readElementText().toDouble();
    if(verbose){cout<<Utilities::getTabString()<<"Smoothing Value: "<<sampleNumber<<endl;}
    return smoothingSet;
}

bool xmlFileAccess::readCalibration(QXmlStreamReader &xmlReader) {
    if(verbose){std::cout<<"//--------- Calibration Mode --------//"<<std::endl;}
    calibrationParametersSet =false;

    //--- phantom type ---//
    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Phantom") {
        if(verbose){std::cout<<"Error: Calibration phantom missing."<<std::endl;}
        return calibrationParametersSet;
    }
    calPhantomName = xmlReader.attributes().value("Name").toString();
    calPhantomIndex = xmlReader.readElementText().toInt();
    if(verbose){std::cout<<Utilities::getTabString()<<"Phantom Name: "<<calPhantomName.toStdString()<<std::endl;}

    //--- Read in Control Point information ---//
    if(calPhantomIndex>=itk::CorticalBone::kMindwaySolidCal && calPhantomIndex<=itk::CorticalBone::kManualControlPtsCal) {

        // control point dimensions
        xmlReader.readNextStartElement();
        if(xmlReader.name() != "radius") {
            if(verbose){std::cout<<"Error: Calibration radius missing."<<std::endl;}
            return calibrationParametersSet;
        }
        calRadius = xmlReader.readElementText().toDouble();
        xmlReader.readNextStartElement();
        if(xmlReader.name() != "number") {
            if(verbose){std::cout<<"Error: Calibration number missing."<<std::endl;}
            return calibrationParametersSet;
        }
        calNumber = xmlReader.readElementText().toInt();

        // work out number of points
        if (calPhantomIndex==itk::CorticalBone::kManualControlPtsCal){
            calVals->SetNumberOfValues(calNumber);
        }
        calPts->SetNumberOfComponents(Dimension); calPts->SetNumberOfTuples(calNumber);

        //--- read in landmark locations ---//
        for(int i=0; i<calNumber; i++) {
            xmlReader.readNextStartElement();
            if(xmlReader.name() != "Point") {
                if(verbose){std::cout<<"Error: Point ["<<i<<"] missing."<<std::endl;}
                return calibrationParametersSet;
            }
            double pt[Dimension];

            xmlReader.readNextStartElement(); pt[0] = xmlReader.readElementText().toDouble();
            xmlReader.readNextStartElement(); pt[1] = xmlReader.readElementText().toDouble();
            xmlReader.readNextStartElement(); pt[2] = xmlReader.readElementText().toDouble();
            calPts->SetTuple(i, pt);
            if(verbose){std::cout<<Utilities::getTabString()<<Utilities::getTabString()<<"Point["<<i<<"]: ["<<pt[0]<<", "<<pt[1]<<", "<<pt[2]<<"]"<<std::endl;}

            // if manual mode read in the average value
            if (calPhantomIndex==itk::CorticalBone::kManualControlPtsCal) {
                xmlReader.readNextStartElement();
                if(xmlReader.name() != "value") {
                    if(verbose){std::cout<<"Error: Point ["<<i<<"] missing value."<<std::endl;}
                    return calibrationParametersSet;
                }
                calVals->SetValue(i, xmlReader.readElementText().toDouble());
            }

            xmlReader.readNextStartElement(); // goes to closing tag
        }
    }



    //--- read in calibration values ---//
    if(calPhantomIndex!=itk::CorticalBone::kManualControlPtsCal) {
        xmlReader.readNextStartElement();
        if(xmlReader.name() != "p0") {
            if(verbose){std::cout<<"Error: Calibration p0 missing."<<std::endl;}
            return calibrationParametersSet;
        }
        calP0 = xmlReader.readElementText().toDouble();
        xmlReader.readNextStartElement();
        if(xmlReader.name() != "p1") {
            if(verbose){std::cout<<"Error: Calibration p1 missing."<<std::endl;}
            return false;
        }
        calP1 = xmlReader.readElementText().toDouble();
    } else {
        calP0 = 0.0;
        calP1 = 0.0;
    }

    if(calPhantomIndex==itk::CorticalBone::kManualQuadraticCal) {
        xmlReader.readNextStartElement();
        if(xmlReader.name() != "p2") {
            if(verbose){std::cout<<"Error: Calibration p2 missing."<<std::endl;}
            return calibrationParametersSet;
        }
        calP2 = xmlReader.readElementText().toDouble();
    } else {
        calP2 = 0.0;
    }

    // TODO check the slope and intercept are valid
    if(verbose){std::cout<<"//--------- Calibration Completed --------//"<<std::endl;}//phantomName.toStdString()<<std::endl;
    calibrationParametersSet = true;
    return calibrationParametersSet;

}

bool xmlFileAccess::readCBConstraint(QXmlStreamReader &xmlReader) {

    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Type") {
        if(verbose){cout<<"Error in xmlFileAccess::readScriptCBConstraint() missing 'Type'"<<endl;}
        return false;
    }

    QString cbConstraintName = xmlReader.attributes().value("Name").toString();
    cbConstraintIndex = xmlReader.readElementText().toInt();
    if(verbose){cout<<Utilities::getTabString()<<"CB Constraint: type="<<cbConstraintName.toStdString();}

    if(cbConstraintIndex==kCBConstraintFixed) {
        xmlReader.readNextStartElement();
        if(xmlReader.name() != "Value") {
            if(verbose){cout<<"Error in xmlFileAccess::readCBConstraint() missing 'ValueFilter'"<<endl;}
            return false;
        }
        fixedCBValue = xmlReader.readElementText().toDouble();
        if(verbose){cout<<",  fixed value="<<fixedCBValue;}
    }
    if(verbose){cout<<endl;}
    return true;
}

bool xmlFileAccess::readCBScriptConstraint(QXmlStreamReader &xmlReader) {

    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Type") {
        if(verbose){cout<<"Error in xmlFileAccess::readScriptCBConstraint() missing 'Type'"<<endl;}
        return false;
    }

    QString cbConstraintName = xmlReader.attributes().value("Name").toString();
    cbConstraintIndex = xmlReader.readElementText().toInt();
    if(verbose){cout<<Utilities::getTabString()<<"CB Constraint: type="<<cbConstraintName.toStdString();}

    if(cbConstraintIndex==kCBConstraintFixed) {
        xmlReader.readNextStartElement();
        if(xmlReader.name() != "ValueFilter") {
            if(verbose){cout<<"Error in xmlFileAccess::readScriptCBConstraint() missing 'ValueFilter'"<<endl;}
            return false;
        }
        fixedCBFilter = xmlReader.readElementText();
        if(verbose){cout<<",  fixed value="<<fixedCBFilter.toStdString().c_str();}
    }
    if(verbose){cout<<endl;}
    return true;
}

bool xmlFileAccess::readProfileAveraging(QXmlStreamReader &xmlReader) {

    averagingSet = true;

    xmlReader.readNextStartElement(); profileAverages[0] = xmlReader.readElementText().toDouble();
    xmlReader.readNextStartElement(); profileAverages[1] = xmlReader.readElementText().toDouble();
    xmlReader.readNextStartElement(); profileAverages[2] = xmlReader.readElementText().toDouble();
    if(verbose){std::cout<<Utilities::getTabString()<<"Profile Averaging: x="<<profileAverages[0]<<", y="<<profileAverages[1]<<", z="<<profileAverages[2]<<std::endl;}

    return true;
}

bool xmlFileAccess::readGlobalSigma(QXmlStreamReader &xmlReader) {
    sigmaSet = true;

    xmlReader.readNextStartElement(); sigmaValues[0] = xmlReader.readElementText().toDouble();
    xmlReader.readNextStartElement(); sigmaValues[1] = xmlReader.readElementText().toDouble();
    xmlReader.readNextStartElement(); sigmaValues[2] = xmlReader.readElementText().toDouble();
    if(verbose){std::cout<<Utilities::getTabString()<<"Global Sigma: x="<<sigmaValues[0]<<", y="<<sigmaValues[1]<<", z="<<sigmaValues[2]<<std::endl;}

    return sigmaSet;

}

bool xmlFileAccess::readClassifierInfo(QXmlStreamReader &xmlReader) {
    classifierThresholdIndexSet = true;

    xmlReader.readNextStartElement();
    classifierName = xmlReader.attributes().value("Name").toString();
    classifierIndex = xmlReader.readElementText().toInt();

    if(verbose){std::cout<<Utilities::getTabString()<<"Classifier Threshold Info: id="<<classifierIndex<<", name="<<classifierName.toStdString().c_str();}

    if(classifierIndex==itk::ClassifierTransform::kMedianManual) {
        xmlReader.readNextStartElement(); classifierWeight = xmlReader.readElementText().toDouble();
        classifierThresholdWeightSet = true;
        if(verbose){std::cout<<", weight="<<classifierWeight;}
    }


    xmlReader.readNextStartElement(); // either read in the classifier info end element or the <Thresholds> start element
    if(xmlReader.name() == "Thresholds") {

        xmlReader.readNextStartElement(); stDensity = xmlReader.readElementText().toDouble();
        xmlReader.readNextStartElement(); cbDensity = xmlReader.readElementText().toDouble();
        xmlReader.readNextStartElement(); thresholdDensity = xmlReader.readElementText().toDouble();
        if(verbose){std::cout<<", ST="<<stDensity<<", CB="<<cbDensity<<", threshold="<<thresholdDensity;}
        classifierThresholdsSet = true;
    }

    if(verbose){std::cout<<std::endl;}

    return classifierThresholdIndexSet;
}



bool xmlFileAccess::readScriptProjectInfo(QXmlStreamReader &xmlReader) {

    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Path") {
        if(verbose){cout<<"Error in xmlFileAccess::readScriptProjectInfo() missing 'Path'"<<endl;}
        return false;
    }
    scriptProjectPath = xmlReader.readElementText();

    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Filter") {
        if(verbose){cout<<"Error in xmlFileAccess::readScriptProjectInfo() missing 'Filter'"<<endl;}
        return false;
    }
    scriptProjectFilter = xmlReader.readElementText();

    return true;
}

bool xmlFileAccess::readScriptImageInfo(QXmlStreamReader &xmlReader) {
    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Filter") {
        if(verbose){cout<<"Error in xmlFileAccess::readScriptImageInfo() missing 'Filter'"<<endl;}
        return false;
    }
    scriptImageFilter = xmlReader.readElementText();

    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Extension") {
        if(verbose){cout<<"Error in xmlFileAccess::readScriptImageInfo() missing 'Extn'"<<endl;}
        return false;
    }
    scriptImageExtn = xmlReader.attributes().value("Name").toString();
    scriptImageExtnIndex = xmlReader.readElementText().toInt();

    return true;
}

bool xmlFileAccess::readScriptMeshInfo(QXmlStreamReader &xmlReader) {
    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Filter") {
        if(verbose){cout<<"Error in xmlFileAccess::readScriptMeshInfo() missing 'Filter'"<<endl;}
        return false;
    }
    scriptMeshFilter = xmlReader.readElementText();

    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Extension") {
        if(verbose){cout<<"Error in xmlFileAccess::readScriptMeshInfo() missing 'Extn'"<<endl;}
        return false;
    }
    scriptMeshExtn = xmlReader.attributes().value("Name").toString();
    scriptMeshExtnIndex = xmlReader.readElementText().toInt();

    return true;
}

bool xmlFileAccess::readScriptCalibrationScales(QXmlStreamReader &xmlReader) {

    calibrationParametersSet = calibrationFileSet = false;

    //--- phantom type ---//
    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Phantom") {
        if(verbose){std::cout<<"Error: Calibration phantom missing."<<std::endl;}
        return calibrationParametersSet;
    }
    calPhantomName = xmlReader.attributes().value("Name").toString();
    calPhantomIndex = xmlReader.readElementText().toInt();
    if(verbose){std::cout<<Utilities::getTabString()<<"Phantom Name: "<<calPhantomName.toStdString()<<std::endl;}


    if(calPhantomIndex==itk::CorticalBone::kFileSpecifiedCal) {
        xmlReader.readNextStartElement();
        if(xmlReader.name() != "ParameterFile") {
            if(verbose){std::cout<<"Error: Calibration parameter file name missing."<<std::endl;}
            return calibrationFileSet;
        }
        calFile = xmlReader.readElementText();
        if(verbose){std::cout<<Utilities::getTabString()<<"Parameter File: "<<calFile.toStdString()<<std::endl;}

        calibrationFileSet = true;

    } else if(calPhantomIndex==itk::CorticalBone::kManualControlPtsCal) {

        xmlReader.readNextStartElement();
        if (xmlReader.name() != "radius") {
            if (verbose) { std::cout << "Error: Calibration radius missing." << std::endl; }
            return calibrationParametersSet;
        }
        calRadius = xmlReader.readElementText().toDouble();
        xmlReader.readNextStartElement();
        if (xmlReader.name() != "number") {
            if (verbose) { std::cout << "Error: Calibration number missing." << std::endl; }
            return calibrationParametersSet;
        }
        calNumber = xmlReader.readElementText().toInt();

        calibrationParametersSet = true;

    } else if (calPhantomIndex != itk::CorticalBone::kManualControlPtsCal) {

        xmlReader.readNextStartElement();
        if (xmlReader.name() != "p0") {
            if (verbose) { std::cout << "Error: Calibration p0 missing." << std::endl; }
            return calibrationParametersSet;
        }
        calP0 = xmlReader.readElementText().toDouble();
        xmlReader.readNextStartElement();
        if (xmlReader.name() != "p1") {
            if (verbose) { std::cout << "Error: Calibration p1 missing." << std::endl; }
            return false;
        }
        calP1 = xmlReader.readElementText().toDouble();

        calP2 = 0.0;

        calibrationParametersSet = true;

    } else if (calPhantomIndex != itk::CorticalBone::kManualControlPtsCal) {
        xmlReader.readNextStartElement();
        if (xmlReader.name() != "p0") {
            if (verbose) { std::cout << "Error: Calibration p0 missing." << std::endl; }
            return calibrationParametersSet;
        }
        calP0 = xmlReader.readElementText().toDouble();
        xmlReader.readNextStartElement();
        if (xmlReader.name() != "p1") {
            if (verbose) { std::cout << "Error: Calibration p1 missing." << std::endl; }
            return false;
        }
        calP1 = xmlReader.readElementText().toDouble();

        xmlReader.readNextStartElement();
        if (xmlReader.name() != "p2") {
            if (verbose) { std::cout << "Error: Calibration p2 missing." << std::endl; }
            return calibrationParametersSet;
        }
        calP2 = xmlReader.readElementText().toDouble();

        calibrationParametersSet = true;

    } else {
        if (verbose) { std::cout << "Error: Calibration INVALID Phantom Index." << std::endl; }
        return false;
    }
    return calibrationParametersSet || calibrationFileSet;
}

bool xmlFileAccess::readScriptActiveTab(QXmlStreamReader &xmlReader) {

    scriptTabName = xmlReader.attributes().value("Name").toString();
    scriptTabIndex = xmlReader.readElementText().toInt();
    return true;
}

bool xmlFileAccess::readModelInfo(QXmlStreamReader &xmlReader) {

    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Model") {
        if(verbose){cout<<"Error in xmlFileAccess::readModelInfo() missing 'Model'"<<endl;}
        return modelSet;
    }
    modelName = xmlReader.attributes().value("Name").toString();
    modelIndex = xmlReader.readElementText().toInt();


    xmlReader.readNextStartElement();
    if (xmlReader.name() != "Optimiser") {
        if (verbose) { cout << "Error in xmlFileAccess::readModelInfo() missing 'Optimiser'" << endl; }
        return modelSet;
    }
    optimiserName = xmlReader.attributes().value("Name").toString();
    optimiserIndex = xmlReader.readElementText().toInt();

    xmlReader.readNextStartElement();
    if(xmlReader.name() != "Scheme") {
        if(verbose){cout<<"Error in xmlFileAccess::readModelInfo() missing 'Scheme'"<<endl;}
        return modelSet;
    }
    schemeName = xmlReader.attributes().value("Name").toString();
    schemeIndex = xmlReader.readElementText().toInt();

    modelSet = true;

    if(verbose){std::cout<<Utilities::getTabString()<<"Model Info: id="<<modelIndex<<" name="<<modelName.toStdString().c_str()
                <<", Scheme Info: id="<<schemeIndex<<" name="<<schemeName.toStdString().c_str()
                <<", Optimiser Info: id="<<optimiserIndex<<" name="<<optimiserName.toStdString().c_str()<<std::endl;}

    return modelSet;

}
