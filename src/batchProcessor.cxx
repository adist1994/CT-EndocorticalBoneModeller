//
// Created by rap58 on 31/08/15.
//

#include "batchProcessor.h"
#include <QtWidgets> // todo remove dependance on Qt


BatchProcessor::BatchProcessor(std::string scriptFileName, int processCode, int startIndexIn, int endIndexIn) {

    bool status;
    std::clog<<"Name="<<scriptFileName.c_str()<<", code="<<processCode<<", start="<<startIndexIn<<", end="<<endIndexIn<<std::endl;

    // create class members
    initialiseGeneralValues(startIndexIn, endIndexIn);

    // read the script
    status = openScript(scriptFileName);

    if(!status) {
        std::cerr<<"Invalid script file or filename = "<<scriptFileName.c_str()<<endl;
        return;
    }

    std::clog<<"Begin Batch Processing Script: "<<scriptFileName.c_str()<<std::endl;

    // begin running the script
    status = runCreationScript(); // always run 'create script' inbetween any specified indicies

    if(status && processCode!=kRunProcess) { // run the process script
        status = runProcessScript();
    }

    if(status) {
        std::clog<<"Successfully Completed Batch Processing of Script"<<std::endl;
    } else {
        std::clog<<"Unsuccessfully Batch Processing of Script"<<std::endl;
    }

}

void BatchProcessor::initialiseGeneralValues(int startIndexIn, int endIndexIn) {

    if(startIndexIn!=-1 && endIndexIn!=-1) {
        startIndex = startIndexIn; endIndex = endIndexIn;
        indiciesSet=true;
    } else {
        startIndex = -1; endIndex = -1;
        indiciesSet=false;
    }
    xmlFileAccessor = new xmlFileAccess();
    corticalBone = itk::CorticalBone::New();
}

void BatchProcessor::initialiseScriptValues() {


    scriptProjectPath = scriptProjectFilter = QString();
    scriptImageFilter = scriptMeshFilter = QString();
    scriptImageExtn = scriptMeshExtn = QString();
    scriptProjectExtn=QString(".xml"); scriptTxtExtn=QString(".txt");
    scriptFixedCBFilter = QString();
    scriptImportFile = QString();

    scriptImageExtnIndex = -1; scriptMeshExtnIndex = -1;
    scriptCBConstraintIndex = -1; scriptSampleNumber = -1;
    scriptModelIndex = -1; scriptOptimiserIndex = -1; scriptClassifierIndex = -1;

    // bools
    scriptFixedCBFilterSet = scriptSampleNumberSet = scriptAveragingSet = scriptImportSet = calibrationSet = false;

    // scalars
    scriptProfileAverages[0] = nan("1"); scriptProfileAverages[1] = nan("1"); scriptProfileAverages[2] = nan("1");
    scriptDensityScales[0] = nan("1"); scriptDensityScales[1] = nan("1"); scriptDensityScales[2] = nan("1");
}

void BatchProcessor::initialiseProjectValues() {
    classifierThresholdsCalculated = classifierThresholdLevelSet =false;
    projectImageName = projectMeshName = "";
    projectImportPath = "";
}

bool BatchProcessor::openScript(std::string fileName) {

    initialiseScriptValues();

    QString file = QString(fileName.c_str());

    bool status = xmlFileAccessor->openScriptFile(file);

    if(status) {

        // general infomation
        xmlFileAccessor->getScriptProjectInfo(scriptProjectPath, scriptProjectFilter);

        // get creation tab information
        xmlFileAccessor->getScriptImageInfo(scriptImageFilter, scriptImageExtn, scriptImageExtnIndex);
        xmlFileAccessor->getScriptMeshInfo(scriptMeshFilter, scriptMeshExtn, scriptMeshExtnIndex);

        xmlFileAccessor->getCBConstraintInfo(scriptCBConstraintIndex, scriptFixedCBFilter);
        scriptFixedCBFilterSet = (scriptCBConstraintIndex==xmlFileAccess::kCBConstraintFixed);

        scriptSampleNumberSet = xmlFileAccessor->isSampleNumberSet();
        if(scriptSampleNumberSet) {
            scriptSampleNumber = xmlFileAccessor->getSampleNumber();
        }

        scriptAveragingSet = xmlFileAccessor->isAveragingSet();
        if(scriptAveragingSet) {
            xmlFileAccessor->getProfileAveraging(scriptProfileAverages);
        }

        scriptImportSet = xmlFileAccessor->isImportSet();
        if(scriptImportSet) {
            scriptImportFile = xmlFileAccessor->getImportBaseName();
        }

        calibrationSet= xmlFileAccessor->areCalibrationParametersSet();
        if(calibrationSet) {
            double p0, p1, p2;
            xmlFileAccessor->getCalibrationParameters(p0, p1, p2);
            scriptDensityScales[0]=p0; scriptDensityScales[1]=p1; scriptDensityScales[2]=p2;

            xmlFileAccessor->getCalibrationName(scriptCalibrationIndex, calibrationPhantomName);
        }

        // get process tab information
        QString modelName;
        xmlFileAccessor->getModelInfo(scriptModelIndex, modelName, scriptSchemeIndex, scriptSchemeName, scriptOptimiserIndex,
                                      scriptOptimiserName);

        if(xmlFileAccessor->areClassifierThresholdsSet()) {
            QString classifierName;
            xmlFileAccessor->getClassifierThresholdIndex(scriptClassifierIndex, classifierName);
        }

    } else {
        std::clog<<"In QScriptFrame::loadScript() script index is invalid"<<endl;
    }

    std::clog<<"<b>---- Script Loaded ----</b>"<<endl;

    return status;

}

bool BatchProcessor::generateBatchFileLists(int code) {

    if(code==kCreate) {

        QStringList directoryList = generateDirectoryList();

        // filter the lists to include only those folders with only one image, mesh and settings file
        scriptProjectNameList = scriptImageNameList = scriptMeshNameList = scriptFixedCBFileList = scriptImportedFileList = QStringList();
        for(int i=0; i<directoryList.size(); i++) {

            QString tempImageString, tempFixedCBString; QStringList tempMeshList;

            QDir directory = QDir(directoryList.at(i));

            // get images, meshes and possibly settings
            directory.setNameFilters(QStringList()<<"*"+scriptImageFilter+"*"+scriptImageExtn);
            if(directory.entryList().size()==1 && scriptImageExtnIndex==kRawExtn) {
                tempImageString = directory.entryList().first();
            } else if(scriptMeshExtnIndex==kDicomExtn) {

                if(directory.entryList().size()>=1) {
                    tempImageString = directory.entryList().first(); // no longer need for just directory //QFileInfo(directory.entryList().first()).absolutePath()+QDir::separator();
                } else { // make one final check for files with no extension
                    directory.setNameFilters(QStringList()<<"*"+scriptImageFilter+"*");
                    QStringList tempImageNoExtnList = directory.entryList();
                    bool isDicomDir = false;
                    for(int j=0; j<tempImageNoExtnList.size(); j++) {
                        if(!tempImageNoExtnList.at(j).contains(QString("."))) {
                            tempImageString = directory.entryList().at(j);
                            isDicomDir = true; break;
                        }
                    }
                    if(!isDicomDir) {
                        continue;
                    }
                }

            } else {

                continue;
            }

            directory.setNameFilters(QStringList()<<"*"+scriptMeshFilter+"*"+scriptMeshExtn);
            tempMeshList = directory.entryList();
            if(tempMeshList.size()<1) {
                continue;
            }

            if(scriptCBConstraintIndex==xmlFileAccess::kCBConstraintFixed) {
                directory.setNameFilters(QStringList()<<"*"+scriptFixedCBFilter+"*"+scriptTxtExtn);
                QStringList tempFixedCBFileList = directory.entryList();

                if(directory.entryList().size() == 1) {
                    tempFixedCBString = tempFixedCBFileList.first();
                } else {
                    continue;
                }
            }

            for(int meshIndex=0; meshIndex<tempMeshList.size(); meshIndex++) {

                QString path = directoryList.at(i)+QDir::separator();

                scriptImageNameList.append(path+tempImageString);
                scriptMeshNameList.append(path+tempMeshList.at(meshIndex));
                if(scriptCBConstraintIndex==xmlFileAccess::kCBConstraintFixed) {
                    scriptFixedCBFileList.append(path+tempFixedCBString);
                }

                QString combinedName;
                QString imageName = QFileInfo(tempImageString).baseName();
                QString meshName = QFileInfo(tempMeshList.at(meshIndex)).baseName();
                if(imageName.contains(meshName)) {
                    combinedName = imageName;
                } else if(meshName.contains(imageName)) {
                    combinedName = meshName;
                } else {
                    combinedName = imageName+QString("_")+meshName;
                }
                QString projectString=path+combinedName+scriptProjectFilter+scriptProjectExtn;
                scriptProjectNameList.append(projectString);
            }
        }

        // create a import mapping list if set
        if(scriptImportSet) {
            QFile importTextFile(scriptImportFile);
            if (importTextFile.open(QIODevice::ReadOnly))
            {
                QTextStream textStream(&importTextFile);
                while (true) {
                    QString line = textStream.readLine();
                    if (line.isNull()) {
                        break;
                    } else {
                        scriptImportedFileList.append(line);
                    }
                }
            }
            if(scriptImportedFileList.size() != scriptProjectNameList.size()) {
                std::cerr<<"Miss-match between import mapping file (n="<<scriptImportedFileList.size()<<") and the create script (n="<<scriptProjectNameList.size()<<")"<<endl;
                scriptProjectNameList = scriptImageNameList = scriptMeshNameList = scriptFixedCBFileList = scriptImportedFileList = QStringList();
            }
        }

    } else if(code==kView) {

        QString projectFilter = QString("*")+scriptProjectFilter+QString("*")+scriptProjectExtn;
        scriptProjectNameList = generateFileList(projectFilter);

        scriptImageNameList = scriptMeshNameList = scriptFixedCBFileList = scriptImportedFileList = QStringList();


    } else if(code==kProcess) {

        QString projectFilter = QString("*") + scriptProjectFilter + QString("*") + scriptProjectExtn;
        scriptProjectNameList = generateFileList(projectFilter);

        scriptImageNameList = scriptMeshNameList = scriptFixedCBFileList = scriptImportedFileList = QStringList();

        // generate results list
        scriptResultsNameList = QStringList();
        for (int i = 0; i < scriptProjectNameList.size(); i++) {

            QFileInfo projectFileInfo = QFileInfo(scriptProjectNameList.at(i));
            QString resultString = projectFileInfo.absolutePath() + QDir::separator() + projectFileInfo.baseName();

            scriptResultsNameList << resultString;
        }

    }
    return true;
}

QStringList BatchProcessor::generateDirectoryList() { // no filter input as filters name and only looking at directories

    QDirIterator directoryIterator(scriptProjectPath, QStringList(), QDir::Dirs, QDirIterator::Subdirectories);

    QStringList dirList = QStringList(scriptProjectPath);

    int i=0;
    while(directoryIterator.hasNext() && i<MAX_NUMBER_OF_FILES_TO_PROCESS) {
        QString directory = directoryIterator.next();
        if(!directory.endsWith(QString(".")) && !directory.endsWith(QString(".."))) {
            dirList.append(directory);
            i++;
        }

    }

    if(i>=MAX_NUMBER_OF_FILES_TO_PROCESS) {
        std::clog<<"Max number of directories reached for path"<<scriptProjectPath.toStdString()<<" others ignored."<<std::endl;
    }

    return dirList;
}

QStringList BatchProcessor::generateFileList(QString filter) {

    QDirIterator directoryIterator(scriptProjectPath, QStringList()<<filter, QDir::Files, QDirIterator::Subdirectories);

    QStringList fileList = QStringList();

    int i=0;
    while(directoryIterator.hasNext() && i<MAX_NUMBER_OF_FILES_TO_PROCESS) {
        fileList<<directoryIterator.next();
        i++;
    }

    if(i>=MAX_NUMBER_OF_FILES_TO_PROCESS) {
        std::clog<<"Max number of files reached for filter"<<filter.toStdString()<<" others ignored."<<std::endl;
    }

    return fileList;
}

bool BatchProcessor::runCreationScript() {
    generateBatchFileLists(kCreate);

    std::clog<<"Start Creation Script Task"<<endl;

    double fixedCBValue = -1;
    startIndex = (indiciesSet) ? startIndex : 0;
    endIndex = (indiciesSet) ? endIndex : scriptProjectNameList.size();
    for(int i=startIndex; i<endIndex; i++) {

        std::clog<<i<<"    "<<scriptProjectNameList.at(i).toStdString().c_str()<<endl;

        xmlFileAccessor->setImageName(scriptImageNameList.at(i));
        xmlFileAccessor->setMeshName(scriptMeshNameList.at(i));
        xmlFileAccessor->setModelInfo(scriptModelIndex, QString(), 0, QString(), scriptOptimiserIndex, QString());

        if(scriptSampleNumberSet) {
            xmlFileAccessor->setSampleNumber(scriptSampleNumber);
        }

        if(scriptCBConstraintIndex==xmlFileAccess::kCBConstraintFixed) {
            QFile cbValueFile(scriptFixedCBFileList.at(i));
            QTextStream cbValueStream(&cbValueFile);
            if (cbValueFile.open(QIODevice::ReadOnly)) {
                fixedCBValue = cbValueStream.readLine().toDouble();
                cbValueFile.close();
            } else {
                fixedCBValue=nan("1");
            }  //cout<<"parsed file: "<<fixedCBFileList.at(i).toStdString()<<". string="<<cbString.toStdString()<<", value="<<fixedCBValue<<endl;
        }
        xmlFileAccessor->setCBConstraintInfo(scriptCBConstraintIndex, fixedCBValue);

        if(scriptAveragingSet) {
            xmlFileAccessor->setProfileAveraging(scriptProfileAverages);
        }

        if(scriptImportSet) {
            xmlFileAccessor->setImportBaseName(scriptImportedFileList.at(i));
        }

        if(calibrationSet) {
            xmlFileAccessor->setCalibrationParameters(scriptDensityScales[0], scriptDensityScales[1], scriptDensityScales[2]);
            xmlFileAccessor->setCalibrationName(scriptCalibrationIndex, calibrationPhantomName);
        }
        xmlFileAccessor->saveProjectFile(scriptProjectNameList.at(i), true);

    }

    std::clog<<"Finished Creation Script Task"<<endl;

    return true;

}

bool BatchProcessor::runProcessScript() {
    generateBatchFileLists(kProcess);

    std::clog<<"Start Processing Script Task"<<endl;

    startIndex = (indiciesSet) ? startIndex : 0;
    endIndex = (indiciesSet) ? endIndex : scriptProjectNameList.size();
    for(int i=startIndex; i<=endIndex; i++) {
        std::clog << i << "    " << scriptProjectNameList.at(i).toStdString().c_str() << endl;
        openProject(scriptProjectNameList.at(i).toStdString());
        corticalBone->setModelIndex(scriptModelIndex);
        if (scriptModelIndex >= itk::CorticalBone::kHighResClassifier && !classifierThresholdsCalculated) { // checked if set when loading file
            corticalBone->calculateThresholds(scriptClassifierIndex);
        }

        bool status = corticalBone->runModellingOverMesh();
        if(!status) {
            std::cerr<<"Warning File "<<i<<" Failed to Correctly Process Over Mesh"<<std::endl;
        } else {
            corticalBone->saveValueArrays(scriptResultsNameList.at(i).toStdString());
            saveProject(scriptProjectNameList.at(i).toStdString()); // save so if any new info gained it is saved
        }
        corticalBone->reset();
    }

    std::clog<<"Finished Processing Script Task"<<endl;

    return true;

}

bool BatchProcessor::openProject(std::string fileName) {

    initialiseProjectValues();

    bool status = xmlFileAccessor->openProjectFile(QString(fileName.c_str()));
    bool meshSet = xmlFileAccessor->isMeshSet(), imageSet = xmlFileAccessor->isImageSet();

    if(imageSet && status) {
        projectImageName = xmlFileAccessor->getImagePath().toStdString();
        status = corticalBone->loadImage(projectImageName);
    } else { // required
        status = false;
    }
    if(meshSet && status) {
        projectMeshName = xmlFileAccessor->getMeshName().toStdString();
        status = corticalBone->loadMesh(projectMeshName);
    } else { // required
        status = false;
    }

    if(xmlFileAccessor->isSampleNumberSet() && status) {
        int localSampleNumber = xmlFileAccessor->getSampleNumber();
        status = corticalBone->setFixedSampleNumber(localSampleNumber);
    }
    if(xmlFileAccessor->areCalibrationParametersSet() && status) {
        int phantomIndex; QString localPhantomName;
        status = xmlFileAccessor->getCalibrationName(phantomIndex, localPhantomName);
        corticalBone->setPhantomIndex(phantomIndex);

        vtkSmartPointer<vtkDoubleArray> pts = xmlFileAccessor->getCalibrationPoints();
        status = corticalBone->runCalibration(pts);

    }
    if(meshSet && imageSet  && status) {
        int modelIndex, schemeIndex, optimiserIndex; QString modelName, schemeName, optimiserName;
        xmlFileAccessor->getModelInfo(modelIndex, modelName, schemeIndex, schemeName, optimiserIndex, optimiserName);
        corticalBone->setModelIndex(modelIndex);
        corticalBone->setOptimiserIndex(optimiserIndex);
        corticalBone->setFittingSchemeIndex(schemeIndex);
    }

    int cbConstraintIndex; double fixedCBValue;
    xmlFileAccessor->getCBConstraintInfo(cbConstraintIndex, fixedCBValue);
    if(cbConstraintIndex==xmlFileAccess::kCBConstraintFWHM && status) {
        corticalBone->turnOnFWHMMode();
    } else if(cbConstraintIndex==xmlFileAccess::kCBConstraintFixed && status) {
        corticalBone->setFixedCBDensity(fixedCBValue);
    }

    if(xmlFileAccessor->isAveragingSet() && status) {
        double averages[3];
        status = xmlFileAccessor->getProfileAveraging(averages);
        corticalBone->turnOnProfileAveraging(averages[0], averages[1], averages[2]);
    }

    if(xmlFileAccessor->isSigmaSet() && status) {
        double sigma[3];
        status = xmlFileAccessor->getSigmaConstraint(sigma);
        corticalBone->turnOnFixedSigma(sigma[0], sigma[1], sigma[2]);
    }

    if(xmlFileAccessor->areClassifierThresholdsSet() && status) { // thresolds already calculated
        double stDensity, cbDensity, thresholdDensity; int classifierIndex;
        xmlFileAccessor->getClassifierThresholds(stDensity, cbDensity, thresholdDensity, classifierIndex);
        classifierThresholdsCalculated = classifierThresholdLevelSet = true;
    } else if(xmlFileAccessor->isClassifierThresholdIndexSet()) { // only threshold level selected, still too calculate
        int classifierIndex; QString classifiedIndexName;
        xmlFileAccessor->getClassifierThresholdIndex(classifierIndex, classifiedIndexName); // todo load just classifier index if set
        classifierThresholdLevelSet = true;
    }

    if(xmlFileAccessor->isImportSet() && status) {
        projectImportPath = xmlFileAccessor->getImportBaseName().toStdString();
        status = corticalBone->setImportProfiles(projectImportPath);
    }

    return status;
}

bool BatchProcessor::saveProject(std::string fileName) {
//---- set attributes ----//
    bool imageSet=corticalBone->isImageSet(), meshSet=corticalBone->isMeshSet();
    if(imageSet) {  // image file name
        xmlFileAccessor->setImageName(QString(projectImageName.c_str()));
    }

    if(meshSet) { // mesh file name
        xmlFileAccessor->setMeshName(QString(projectMeshName.c_str()));
    }

    // model info
    if(imageSet && meshSet) {
        int modelIndex, optimiserIndex, schemeIndex;
        std::string modelName, optimiserName, schemeName;
        corticalBone->getModelSelection(modelIndex, modelName);
        corticalBone->getOptimiserSelection(optimiserIndex, optimiserName);
        corticalBone->getSchemeSelection(schemeIndex, schemeName);
        xmlFileAccessor->setModelInfo(modelIndex, QString(modelName.c_str()), schemeIndex, QString(schemeName.c_str()), optimiserIndex, QString(optimiserName.c_str()));
    }

    // fixed sample number
    if(corticalBone->isSampleNumberFixed()) {
        int fixedSampleNumber = corticalBone->getSampleNumber();
        xmlFileAccessor->setSampleNumber(fixedSampleNumber);
    }

    // imported parameters
    if(corticalBone->isParametersImported()) {
        xmlFileAccessor->setImportBaseName(QString(projectImportPath.c_str()));

    }

    // calibration (positions, type, values)
    if(corticalBone->isCalibrated()) {

        std::string phantomName; int phantomIndex;
        corticalBone->getCalibrationPhantom(phantomName, phantomIndex);

        xmlFileAccessor->setCalibrationName(phantomIndex, QString(phantomName.c_str()));

        vtkSmartPointer<vtkDoubleArray> pts;
        corticalBone->getCalibrationPoints(phantomIndex, pts);
        xmlFileAccessor->setCalibrationPoints(pts);

        double p0, p1, p2;
        corticalBone->getCalibrationValues(p0, p1, p2);
        xmlFileAccessor->setCalibrationParameters(p0, p1, p2);

    }

    // in case of fixed CB or set to FWHM mode - don't both if set to None as this is the default
    int cbConstraintIndex; double fixedCBValue=nan("1");
    if(corticalBone->isSetToFWHMMode()) {
        cbConstraintIndex = xmlFileAccess::kCBConstraintFWHM;
    } else if(corticalBone->isCBFixed()) {
        cbConstraintIndex = xmlFileAccess::kCBConstraintFixed;
        fixedCBValue = corticalBone->getFixedCB();
    } else {
        cbConstraintIndex = xmlFileAccess::kCBConstraintNone;
    }
    xmlFileAccessor->setCBConstraintInfo(cbConstraintIndex, fixedCBValue);


    if(corticalBone->isProfileAveragingOn()) { // in case of HR profile averaging
        double profileAverages[3];
        corticalBone->getProfileAveragingValues(profileAverages[0], profileAverages[1], profileAverages[2]);
        xmlFileAccessor->setProfileAveraging(profileAverages);
    }

    if(corticalBone->isSigmaFixed()) { // fixed average radius
        double sigma[3];
        corticalBone->getFixedSigma(sigma[0], sigma[1], sigma[2]);
        xmlFileAccessor->setSigmaConstraint(sigma);
    }

    if(corticalBone->areThresholdsSet()) { // classifier thresholds calculated
        double stDensity, cbDensity, thresholdDensity; std::string classifierName; int classifierIndex;
        corticalBone->getClassifierInfo(stDensity, cbDensity, thresholdDensity, classifierIndex, classifierName);
        corticalBone->getclassifierName(classifierName);
        xmlFileAccessor->setClassifierThresholdInfo(stDensity, cbDensity, thresholdDensity, classifierIndex,
                                                    QString(classifierName.c_str()));
    }

    // save file
    xmlFileAccessor->saveProjectFile(QString(fileName.c_str()), false); // false as saving an existing project, as opposed to generating one
    
    return true;
}