/****************************************************************************
**
** Copyright (C) 2013 Digia Plc and/or its subsidiary(-ies).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Digia Plc and its Subsidiary(-ies) nor the names
**     of its contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include <QtWidgets>
#include <iostream>
#include <istream>
#include <sstream>

#include "mainwindow.h"
#include "scriptframe.h"
#include "qdebugstream.h"
#include "visualisationframe.h"
#include "controlframe.h"
//#include "projectparser.h"
#include "xmlFileAccess.h"

/* Public Methods */
MainWindow::MainWindow() {

    // custom private functions
    createDockWindows();
    createActions();
    createMenus();
    createToolBars();
    createStatusBar();
    setState(0);

    //imageFilePath = QString();
    //meshFilePath = QString();

    readSettings(); // set user preferences

    curName = QString();
    // establish a signal-slot connection
    //connect(debugText->document(), SIGNAL(contentsChanged()), this, SLOT(documentWasModified()));

    setCurrentFile("");
    setUnifiedTitleAndToolBarOnMac(true);

    upToDate = true;

}

void MainWindow::closeEvent(QCloseEvent *event) {
    std::cout << "closing file" << std::endl;
    writeSettings();
}

void MainWindow::saveResults(QString fileName) {
    std::cout<<"<b>----- Save CBM Values -----</b>"<<std::endl;
    bool status = visualisationFrame->saveValueArrays(&fileName);

    if (status) {
        updateStatusBar(tr("Data Files saved"));
    }
    std::cout<<"<b>----- Finished Save Values -----</b>"<<std::endl;
}

void MainWindow::saveDisplay(QString fileName) {
    std::cout<<"<b>----- Save CBM Displays -----</b>"<<std::endl;

    QString sliceName = QString(fileName); QString profileName = QString(fileName); QString threeDimName = QString(fileName);

    if(fileName.endsWith(".png")) {
        sliceName.replace(QString(".png"), QString("_slice.png"));
        profileName.replace(QString(".png"), QString("_profile.png"));
        threeDimName.replace(QString(".png"), QString("_threeDim.png"));
    } else if(fileName.endsWith(".eps")) {
        sliceName.replace(QString(".eps"), QString("_slice.eps"));
        profileName.replace(QString(".eps"), QString("_profile.eps"));
        threeDimName.replace(QString(".eps"), QString("_threeDim.eps"));
    }



    bool status = visualisationFrame->saveDisplays(&sliceName, &profileName, &threeDimName); // todo is & nessecary?

    if (status) {
        updateStatusBar(tr("Displays saved"));
    }
    std::cout<<"<b>----- Finished Save Displays -----</b>"<<std::endl;
}

QPlainTextEdit* MainWindow::getScriptOutputEditor() {
    return scriptText;
}

//------- Key Press Events -----------//
void MainWindow::keyPressEvent(QKeyEvent *event) {
    if(event->key() == Qt::Key_N && scriptControlFrame->isViewRunning()) {
        scriptControlFrame->nextProjectFile();
    }
}

//------- Tool bar / menu slots -------//
void MainWindow::openScriptSlot() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Script"), scriptFilePath, tr("Text Files (*.xml)"));
    if (!fileName.isEmpty()) {

        // update current path
        //curPath = QFileInfo(fileName).absolutePath();
        //loadScript(fileName);
        bool status = scriptControlFrame->loadScript(fileName);

        if(status) {
            scriptFilePath = fileName;
        }

        // todo - put under 'openCBM' when implemeted
        setCurrentFile(fileName);
        updateStatusBar(tr("File loaded"));

        upToDate = false;

    }
}

void MainWindow::openImageSlot() {

//  // get directory
//  QString filePath = QFileDialog::getExistingDirectory(this, tr("Open DICOM Directory"), imageFilePath);
    // get file [.dcm, .mhd]
    QString defaultFilter = tr("CT Files (*.dcm *.mhd *.tif *.QCT)");
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open CT Image"), imageFilePath, tr("CT Files (*.dcm *.mhd *.tif *.QCT);;All files (*)"), &defaultFilter);
    if (!fileName.isEmpty()) {

        // TODO add support for vtk files
        openImage(fileName);

        upToDate = false;

    } else {
        std::clog << "Error - selected file is empty" << std::endl;
    }

}

void MainWindow::openMeshSlot() {

    // get file [.wrl, .stl, .obj]
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), meshFilePath, tr("Mesh Files (*.wrl *.stl *.obj *.ply)"));
    if (!fileName.isEmpty()) {
        openMesh(fileName);

        upToDate = false;

    } else {
        std::clog << "Error - selected file is empty" << std::endl;
    }

}

void MainWindow::openFileSlot() {

    // check to see if there is possibly unsaved work
    closeFileSlot();

    // get file [.xml]
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Project"), curPath, tr("Project Files (*.xml)"));
    if (!fileName.isEmpty()) {

        bool status = openProject(fileName);

        upToDate = true;

        //setCurrentFile(fileName);
        //updateStatusBar(tr("Project File loaded"));

    } else {
        std::clog << "Error - selected filename is empty" << std::endl;
    }

}

void MainWindow::openParamsSlot() {

    // todo - check state and raise error if this should not be called

    importProfile();

}

void MainWindow::saveFileSlot() {

    // get file [.xml]
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Project"), curPath.append(QDir::separator()).append(curName), tr("Project Files (*.xml)"));

    if (!fileName.isEmpty()) {

        if(!fileName.endsWith(".xml")) {
            fileName.append(".xml");
        }

        if(QFileInfo(fileName).exists()) { // if it already exists save the selected file
            curPath = QFileInfo(fileName).absolutePath();
            curName = QFileInfo(fileName).baseName();
        } else { // if it doesn't exist create a directory and save it their
            curPath = QFileInfo(fileName).absolutePath()+QDir::separator()+QFileInfo(fileName).baseName()+QDir::separator();
            curName = QFileInfo(fileName).baseName();
            fileName = curPath+curName+tr(".xml");
            QDir dir = QDir(); dir.mkdir(curPath);
        }

        //fileParser->writeFile(fileName);
        saveProject(fileName);

        upToDate = true;

        /*QFile file(fileName);

        if (!file.open(QFile::WriteOnly)) {
              QString msg = tr("Failed to open %1\n%2").arg(fileName).arg(file.errorString());
              QMessageBox::warning(this, tr("Error"), msg);
              return;
        }*/

        /*QByteArray geo_data = saveGeometry();
        QByteArray layout_data = saveState();

        bool ok = file.putChar((uchar)geo_data.size());
        if (ok)
            ok = file.write(geo_data) == geo_data.size();
        if (ok)
            ok = file.write(layout_data) == layout_data.size();

        if (!ok) {
            QString msg = tr("Error writing to %1\n%2")
                            .arg(fileName)
                            .arg(file.errorString());
            QMessageBox::warning(this, tr("Error"), msg);
            return;
        }*/
        updateStatusBar(tr("Project File saved"));
    } else {
        std::clog << "Error - selected filename is empty" << std::endl;
    }

}

void MainWindow::saveMeshSlot() {

    // save file [.wrl, .stl, .obj]
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Mesh File"), meshFilePath, tr("Mesh Files (*.obj *.stl *.wrl)"));

    if (!fileName.isEmpty()) {

        if(!fileName.endsWith(".stl") && !fileName.endsWith(".obj") && !fileName.endsWith(".wrl")) {
            fileName.append(".obj");
        }

        bool status = visualisationFrame->saveMesh(&fileName);

        if(status) {updateStatusBar(tr("Mesh File saved"));}
    } else {
        std::clog << "Error - selected filename is empty" << std::endl;
    }
}

void MainWindow::saveValuesSlot() {

    QString valuesPath = curPath.append(QDir::separator()).append(curName);

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Data Arrays"), valuesPath, tr("Data Files (*.txt)"));

    if (!fileName.isEmpty()) {

        // do something
        saveResults(fileName);

    } else {
        std::clog << "Error - selected filename is empty" << std::endl;
    }
}

void MainWindow::saveDisplaySlot() {
    QString valuesPath = curPath.append(QDir::separator()).append(curName);

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Display"), valuesPath, tr("Image Files (*.png *.eps)"));

    if (!fileName.isEmpty()) {

        if(!fileName.endsWith(".png") && !fileName.endsWith(".eps")) {
            fileName.append(".png");
        }

        // do something
        saveDisplay(fileName);

    } else {
        std::clog << "Error - selected filename is empty" << std::endl;
    }
}

void MainWindow::closeFileSlot() {
    // check for changes & ask to save
    if(visualisationFrame->getState() != 0) {
        if(!upToDate) {
            QMessageBox::StandardButton reply;
            reply = QMessageBox::question(this, "Potentially Unsaved Changes", "Do you want to save the current project before closing it?", QMessageBox::Yes|QMessageBox::No);
            if (reply == QMessageBox::Yes) {
                //qDebug() << "Yes was clicked";
                saveFileSlot();
                QApplication::quit();
            } else {
                //qDebug() << "Yes was *not* clicked";
            }
        }

        // close file
        closeProject();
        setCurrentFile("");


        // reset curName
        curName = curPath;

        upToDate = true;

    } else {
        //cerr<<"MainWindow::closeFileSlot Nothing open, so nothing to close."<<endl; // disable if nothing open
    }

}

// Tab Contorl Slots
void MainWindow::controlTabChanged(int index) {
    if(index==0) { // display projects
        scriptEditor->hide(); visualisationFrame->show();
    } else if(index==1) { // display script
        scriptEditor->show(); visualisationFrame->hide();
    } else {
        cerr<<"Error in MainWindow::controlTabChanged() invaid index entered="<<index<<endl;
    }
}

// control display slots
void MainWindow::showDisplaySlot(bool state) {
    projectControlFrame->setDisplayBoxVisibilities(state);
}

void MainWindow::showCalibrationSlot(bool state) {
    projectControlFrame->setCalibrationBoxVisibilities(state);
}

void MainWindow::showConstraintsSlot(bool state) {
    projectControlFrame->setConstraintBoxVisibilities(state);
}

void MainWindow::showThresholdSlot(bool state) {
    projectControlFrame->setThresholdBoxVisibilities(state);
}

void MainWindow::showProfileSettingsSlot(bool state) {
    projectControlFrame->setProfileSettingsBoxVisibilities(state);
}

void MainWindow::showImportSlot(bool state) {
    projectControlFrame->setImportBoxVisibilities(state);
}

// control frame slots
void MainWindow::disablePtMeasures() {

    this->visualisationFrame->disablePtMeasures(); // or could call updateMeasurementEnables();
    updateState();
}

void MainWindow::updatePtMeasures() {
    updateMeasurementEnables();
    updateState();
}

void MainWindow::runModellingOverMesh() {
    std::cout<<"<b>----- Run CBM Over Mesh -----</b>"<<std::endl;
    bool status = this->visualisationFrame->runModellingOverMesh();

    std::cout<<"<b>----- Finish Running CBM -----</b>"<<std::endl;

    updateGUIEntries();
    updateState();
    updateMeasurementEnables();

}

void MainWindow::runCalibration() {

    bool status;
    if(this->projectControlFrame->isManualCalibrationMode()) {
        double p0, p1, p2;
        this->projectControlFrame->getCalibrationValues(p0, p1, p2);
        status = this->visualisationFrame->setCalibration(p0, p1, p2);
    } else {
        status = this->visualisationFrame->runCalibration();
    }

    if(status) {
        double p0, p1, p2;
        this->visualisationFrame->getCalibrationValues(p0, p1, p2);
        this->projectControlFrame->setCalibrationScalingValues(p0, p1, p2);
    } else {
        this->projectControlFrame->removeCalibrationScalingValues();
    }

    updateState();
    updateMeasurementEnables();

}

void MainWindow::documentWasModified(bool modified) {
    setWindowModified(modified);
}

void MainWindow::togglePtMeasures(bool buttonState) {

    if (buttonState) {
        this->visualisationFrame->enablePtMeasures();
    } else {
        this->visualisationFrame->disablePtMeasures();
    }
}

void MainWindow::toggleCalibrationMode(bool buttonState){

    if (buttonState) {
        this->visualisationFrame->enableCalibrationMode();
    } else {
        this->visualisationFrame->disableCalibrationMode();
    }

    updateMeasurementEnables();
    updateState();
}

void MainWindow::displayParameterChanged(int index) {
    this->visualisationFrame->setDisplayParameter(index);
}

void MainWindow::displayMeshChanged(int index) {
    this->visualisationFrame->setDisplayMesh(index);
}

void MainWindow::calibrationPhantomChanged(int index) {

    this->visualisationFrame->setCalibrationPhantom(index);

    if(projectControlFrame->isManualCalibrationPtsMode() && this->projectControlFrame->isManualCalibrationPtsUpToDate()) { // set
        double radius; int number;
        projectControlFrame->getCalibrationRadius(radius);
        projectControlFrame->getCalibrationNumber(number);

        this->visualisationFrame->setCalibrationPtGeometry(radius, number);
    }

    toggleCalibrationMode(projectControlFrame->isCalibrationReady());

    updateMeasurementEnables();
    updateState();
}

void MainWindow::updateCalibrationRadius() {

    if(this->projectControlFrame->isManualCalibrationPtsUpToDate()) { // set
        double radius; projectControlFrame->getCalibrationRadius(radius);
        int number; projectControlFrame->getCalibrationNumber(number);
        this->visualisationFrame->setCalibrationPtGeometry(radius, number);
    } else {
        this->visualisationFrame->disableCalibrationMode();
    }

    toggleCalibrationMode(projectControlFrame->isCalibrationReady());

    updateMeasurementEnables();
    updateState();
}

void MainWindow::updateCalibrationNumber() {

    if(this->projectControlFrame->isManualCalibrationPtsUpToDate()) { // set
        double radius; projectControlFrame->getCalibrationRadius(radius);
        int number; projectControlFrame->getCalibrationNumber(number);
        this->visualisationFrame->setCalibrationPtGeometry(radius, number);
    } else {
        this->visualisationFrame->disableCalibrationMode();
    }

    toggleCalibrationMode(projectControlFrame->isCalibrationReady());

    updateMeasurementEnables();
    updateState();
}

void MainWindow::disableCalibration() {
    this->visualisationFrame->disableCalibrationMode();

    updateMeasurementEnables();
    updateState();
}

void MainWindow::modelFunctionChanged(int index) {
    this->visualisationFrame->setModelFunction(index);

    // update what control frame visibilities
    bool showThresholds = projectControlFrame->isClassifierMode();
    showThresholdAct->setChecked(showThresholds);

    bool showConstraints = projectControlFrame->isConstraintMode();
    showConstraintsAct->setChecked(showConstraints);

    bool showImports = !projectControlFrame->isCalibrationMode();
    showImportsAct->setChecked(showImports);


    updateGUIEntries();
    updateMeasurementEnables();
    updateState();
}

void MainWindow::optimiserChanged(int index) {
    this->visualisationFrame->setOptimiser(index);

    updateState();
    // updateMeasurementEnables();
}

void MainWindow::fittingSchemeChanged(int index) {

    if(index!=-1) { // only call scheme change if it's been set to an actual value as apposed to be set to empty
        this->visualisationFrame->setfittingScheme(index);

        updateGUIEntries();
        updateState();
        updateMeasurementEnables();
    }

}

// sampling
void MainWindow::toggleSampleNumber(bool buttonState) {


    if(buttonState) { // disabling button presses taken car of by Control Frame

        if(projectControlFrame->isSampleNumberSet()) {
            sampleNumberChanged();
        }
    } else {
        visualisationFrame->removeFixedSampledMode();
    }

    updateMeasurementEnables();
    updateState();

}

void MainWindow::sampleNumberChanged() {
    // change sample manually

    int sampleNumber = this->projectControlFrame->getSampleNumber();

    bool status = this->visualisationFrame->setFixedSampleNumber(sampleNumber);

    updateMeasurementEnables();
    updateState();
}

// smoothing
void MainWindow::smoothingChanged() {

    // change sample manually

    double smoothingRadius = this->projectControlFrame->getSmoothingRadius();

    this->visualisationFrame->setSmoothingRadius(smoothingRadius);

    updateMeasurementEnables();
    updateState();

}

void MainWindow::toggleSmoothing(bool buttonState) {

    if(buttonState) { // remove fixed sample value
        visualisationFrame->turnOnSmoothingMode();
    } else {
        // do nothing as taken care of in Control frame
        visualisationFrame->turnOffSmoothingMode();
    }

    updateMeasurementEnables();
    updateState();

}


// fix CB Density
void MainWindow::toggleCBDensityConstraint(bool state) {

    if(state) {

        // check if fwhm or cb fixed
        if(this->projectControlFrame->isFWHMMode()) {
            turnOnFWHMCBDensity(); // update state/enables taken care of in this method
        } else {
            turnOnFixCBDensity(); // update state/enables taken care of in this method
        }


    } else {
        this->visualisationFrame->turnOffFWHMMode();
        this->visualisationFrame->removeFixedCBDensity();

        updateMeasurementEnables();
        updateState();
    }

}

void MainWindow::turnOnFixCBDensity() {

    this->visualisationFrame->turnOffFWHMMode();

    fixedCBDensityChanged();

    updateMeasurementEnables();
    updateState();

}

void MainWindow::turnOnFWHMCBDensity() {

    this->visualisationFrame->turnOnFWHMMode();

    updateMeasurementEnables();
    updateState();

}

void MainWindow::fixedCBDensityChanged() {

    double CBDensity;
    bool status = this->projectControlFrame->getCBDensity(CBDensity);

    if(status) {
        this->visualisationFrame->setFixedCBDensity(CBDensity);
    }

    updateMeasurementEnables();
    updateState();
}


// fix sigma
void MainWindow::toggleFixedSigma(bool sigmaState) {

    if(sigmaState) {
        // turn on fixed mode
        sigmaChanged(); // update state/enables taken care of in this method
    } else {
        // remove fixed value
        this->visualisationFrame->turnOffFixedSigma();

        updateMeasurementEnables();
        updateState();
    }
}

void MainWindow::sigmaChanged() {

    double x, y, z;
    bool allSet = this->projectControlFrame->getSigmaValues(x, y, z);

    if(allSet) {
        this->visualisationFrame->turnOnFixedSigma(x, y, z);
    }

    updateMeasurementEnables();
    updateState();
}

// thresholds
void MainWindow::thresholdSelectorChanged(int index) {
    // remove thresholds + update gui
    visualisationFrame->removeThresholds();
    updateGUIEntries();
    updateMeasurementEnables();
    updateState();

}

void MainWindow::setThresholdIndex(int index) {
    // called by script
    if(index != projectControlFrame->getThresholdSelection()) {
        visualisationFrame->removeThresholds();
        projectControlFrame->setThresholdSelection(index);
        updateGUIEntries();
        updateMeasurementEnables();
        updateState();
    }

}

void MainWindow::runThresholdCalculation() {

    bool status;
    if(this->projectControlFrame->isManualThresholdsMode()) {
        double stThreshold, cbThreshold, threshold; int classifierThresholdIndex;
        this->projectControlFrame->getThresholdValues(stThreshold, cbThreshold, threshold, classifierThresholdIndex);
        status = this->visualisationFrame->setClassifierThresholdInfo(stThreshold, cbThreshold, threshold, classifierThresholdIndex);
    } else if(this->projectControlFrame->isMedianManualThresholdsMode()) {
        int classifierThresholdIndex = projectControlFrame->getThresholdSelection();
        double weight;
        status = projectControlFrame->getThresholdWeight(weight);
        status = visualisationFrame->calculateThresholds(classifierThresholdIndex, weight);
    } else {
        int classifierThresholdIndex = projectControlFrame->getThresholdSelection();
        status = visualisationFrame->calculateThresholds(classifierThresholdIndex);
    }

    updateGUIEntries();
    updateMeasurementEnables();
    updateState();
}

// profile averaging
void MainWindow::toggleProfileAverging(bool state) {


    if(state) { // if profile averaging enabled check if specified

        ProfileAveragingValuesChanged();

    } else { // if unchecked turn off profile averaging
        this->visualisationFrame->turnOffProfileAvergaing();

        updateMeasurementEnables();
        updateState();
    }
}

void MainWindow::ProfileAveragingValuesChanged() {

    double x, y, z;
    bool allSet = this->projectControlFrame->getAveragingValues(x, y, z);

    if(allSet) {
        this->visualisationFrame->turnOnProfileAveraging(x, y, z);
    }

    updateMeasurementEnables();
    updateState();
}


// import profiles / parameters
void MainWindow::importProfile() {
    // get file [.txt]

    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Model Parameters File"), importParamsPath, tr("Model Parameters Files (*.txt)"));
    if (!fileName.isEmpty()) {
        bool status = openImportParams(fileName);

    } else {
        projectControlFrame->setImportPath(QString());

        std::clog << "Error - selected file is empty" << std::endl;
    }

    updateMeasurementEnables();
    updateState();
}

void MainWindow::removeImportedProfile() {
    projectControlFrame->setImportPath(QString());
    visualisationFrame->removeImportedProfile();

    updateMeasurementEnables();
    updateState();
}


//--------------open/close/save files ---------------//
bool MainWindow::openMesh(QString fileName) {

    // load file and update the display

    bool status = visualisationFrame->setMesh(&fileName);

    // updateMeasurementEnables();
    updateState();

    if (status) {
        meshFilePath = fileName;
        updateStatusBar(tr("Mesh loaded")); // todo - put under 'openCBM' when implemeted
    } else {
        cout<<"ERROR: Failed to load mesh"<<endl;
    }
    return status;

}

bool MainWindow::openImage(QString fileName) {

    bool status = visualisationFrame->setImage(fileName);

    // updateMeasurementEnables();
    updateState();

    if (status) {
        imageFilePath = fileName; // TODO check successful loading before setting this

        updateStatusBar(tr("CT Image loaded"));
    } else {
        cout<<"ERROR: Failed to load image"<<endl;
    }
    return status;
}

bool MainWindow::openImportParams(QString fileName) { // load importparams


    bool status = visualisationFrame->openParameters(&fileName);

    if(status) {
        importParamsPath = fileName;

        QString baseName = QFileInfo(importParamsPath).baseName();
        QString baseNameStub = baseName.left(baseName.lastIndexOf(tr("_")));

        projectControlFrame->setImportPath(baseNameStub);
        updateStatusBar(tr("Parameters File loaded"));

    } else {
        projectControlFrame->setImportPath(QString());
        updateStatusBar(tr("Error: Unable to Load Import File"));
        std::cerr << "Error - Unable to Load Import File" << std::endl;
    }

    // updateMeasurementEnables();
    updateState();

    return status;
}

bool MainWindow::saveProject(QString fileName) {

    xmlFileAccessor->reset(); // ensure state reset

    //---- set attributes ----//
    bool meshSet = visualisationFrame->isMeshSet();
    bool imageSet = visualisationFrame->isImageSet();

    if(imageSet) {  // image file name
        xmlFileAccessor->setImageName(imageFilePath);
    }


    if(meshSet) { // mesh file name
        xmlFileAccessor->setMeshName(meshFilePath);
    }

    // fixed sample number
    if(visualisationFrame->isSampleNumberFixed()) {
        int fixedSampleNumber = visualisationFrame->getSampleNumber();
        xmlFileAccessor->setSampleNumber(fixedSampleNumber);
    }

    // model info
    if(imageSet && meshSet) {
        int modelIndex, optimiserIndex, schemeIndex;
        QString modelName, optimiserName, schemeName;
        visualisationFrame->getModelSelection(modelIndex, modelName);
        visualisationFrame->getOptimiserSelection(optimiserIndex, optimiserName);
        visualisationFrame->getSchemeSelection(schemeIndex, schemeName);
        xmlFileAccessor->setModelInfo(modelIndex, modelName, schemeIndex, schemeName, optimiserIndex, optimiserName);
    }

    // imported parameters
    if(visualisationFrame->isParametersImported()) {
        xmlFileAccessor->setImportBaseName(importParamsPath);
    }

    // calibration (positions, type, values)
    if(imageSet && this->visualisationFrame->isCalibrated()) {
        std::string phantomType; int phantomIndex;
        this->visualisationFrame->getCalibrationPhantomType(phantomType, phantomIndex);

        xmlFileAccessor->setCalibrationName(phantomIndex, QString(phantomType.c_str()));

        if(phantomIndex>=itk::CorticalBone::kMindwaySolidCal && phantomIndex<=itk::CorticalBone::kEuropeanSpineCal) {
            double radius; int number;
            this->visualisationFrame->getCalibrationPtGeometry(radius, number);
            xmlFileAccessor->setCalibrationPtGeometry(radius, number);

            vtkSmartPointer<vtkDoubleArray> pts = this->visualisationFrame->getCalibrationPoints();
            xmlFileAccessor->setCalibrationPoints(pts);

        } else if(phantomIndex==itk::CorticalBone::kManualControlPtsCal) {
            double radius; int number;
            this->visualisationFrame->getCalibrationPtGeometry(radius, number);
            xmlFileAccessor->setCalibrationPtGeometry(radius, number);

            vtkSmartPointer<vtkDoubleArray> pts = this->visualisationFrame->getCalibrationPoints();
            vtkSmartPointer<vtkDoubleArray> vals = this->visualisationFrame->getCalibrationValues();
            xmlFileAccessor->setCalibrationPoints(pts, vals);
        }

        if(phantomIndex!=itk::CorticalBone::kManualControlPtsCal) {
            double p0, p1, p2;
            this->visualisationFrame->getCalibrationValues(p0, p1, p2);
            xmlFileAccessor->setCalibrationParameters(p0, p1, p2);
        }
    }

    if(imageSet && meshSet) { // values that only have meaning if the image and mesh have bothbeen loaded

        // in case of fixed CB or set to FWHM mode - don't both if set to None as this is the default
        int cbConstraintIndex; double fixedCBValue=nan("1");
        if(this->visualisationFrame->isSetToFWHMMode()) {
            cbConstraintIndex = xmlFileAccess::kCBConstraintFWHM;
        } else if(visualisationFrame->isFixedCB()) {
            cbConstraintIndex = xmlFileAccess::kCBConstraintFixed;
            fixedCBValue = this->visualisationFrame->getFixedCB();
        } else {
            cbConstraintIndex = xmlFileAccess::kCBConstraintNone;
        }
        xmlFileAccessor->setCBConstraintInfo(cbConstraintIndex, fixedCBValue);


        if(this->visualisationFrame->isProfileAveragingOn()) { // in case of HR profile averaging

            double profileAverages[3];
            this->visualisationFrame->getProfileAveragingValues(profileAverages[0], profileAverages[1], profileAverages[2]);
            xmlFileAccessor->setProfileAveraging(profileAverages);

        }
        // smoothing radius for smoothed mode fitting
        if(visualisationFrame->isSmoothingOn()) {
            double smoothingRadius = visualisationFrame->getSmoothingValue();
            xmlFileAccessor->setSmoothing(smoothingRadius);
        }

        // fixed sample number
        if(visualisationFrame->isSigmaFixed()) {
            double sigma[3];
            this->visualisationFrame->getFixedSigma(sigma[0], sigma[1], sigma[2]);
            xmlFileAccessor->setSigmaConstraint(sigma);
        }

        // thresholds
        if(projectControlFrame->isClassifierMode()) {
            int thresholdMode; QString thresholdModeName;
            if(projectControlFrame->isMedianManualThresholdsMode() && visualisationFrame->areThresholdsSet()) {
                double notBoneDensity, boneDensity, thresholdDensity, weight;
                this->visualisationFrame->getClassifierThresholdInfo(notBoneDensity, boneDensity, thresholdDensity, weight, thresholdMode, thresholdModeName);
                xmlFileAccessor->setClassifierThresholdInfo(notBoneDensity, boneDensity, thresholdDensity, thresholdMode, thresholdModeName, weight);
            } else if(projectControlFrame->isMedianManualThresholdsMode() && !visualisationFrame->areThresholdsSet()) {
                thresholdMode = projectControlFrame->getThresholdSelection();
                thresholdModeName = projectControlFrame->getClassifierThresholdName();
                double weight; projectControlFrame->getThresholdWeight(weight);
                xmlFileAccessor->setClassifierThresholdInfo(thresholdMode, thresholdModeName, weight);
            } else if(visualisationFrame->areThresholdsSet()) {
                double notBoneDensity, boneDensity, thresholdDensity;
                this->visualisationFrame->getClassifierThresholdInfo(notBoneDensity, boneDensity, thresholdDensity, thresholdMode, thresholdModeName);
                xmlFileAccessor->setClassifierThresholdInfo(notBoneDensity, boneDensity, thresholdDensity, thresholdMode, thresholdModeName);
            } else {
                thresholdMode = projectControlFrame->getThresholdSelection();
                thresholdModeName = projectControlFrame->getClassifierThresholdName();
                xmlFileAccessor->setClassifierThresholdInfo(thresholdMode, thresholdModeName);
            }

        }
    }


    // save file
    xmlFileAccessor->saveProjectFile(fileName, false);

    // if locally generated results save them
    if(!visualisationFrame->areResultsLoaded() && visualisationFrame->isMeshMeasured()) {
        visualisationFrame->saveValueArrays(&fileName);
    }
    // save results
    if(imageSet && this->visualisationFrame->isCalibrated()) {
        this->visualisationFrame->saveCalibration(&fileName);
    }

    /***********************
     * TODO
     * 1. specify model properties: optimiser type, model type
     * 2. display options: mesh/volume, point measurement mode
     *  */


    return true;

}

bool MainWindow::openProject(QString fileName){

    xmlFileAccessor->openProjectFile(fileName); // reset included

    // TODO - make it so parsing is completed then values are set
    bool status = true;

    bool imageSet = xmlFileAccessor->isImageSet(), meshSet = xmlFileAccessor->isMeshSet();

    if(status && imageSet) {
        status = openImage(xmlFileAccessor->getImagePath());
    }
    if(status && meshSet) {
        status = openMesh(xmlFileAccessor->getMeshName());
    }
    if(status && imageSet && meshSet) {
        int modelIndex, schemeIndex, optimiserIndex;
        QString modelName, schemeName, optmiserName;
        xmlFileAccessor->getModelInfo(modelIndex, modelName, schemeIndex, schemeName, optimiserIndex, optmiserName);
        projectControlFrame->setModelIndex(modelIndex); // automatically triggers: modelFunctionChanged(modelIndex);
        projectControlFrame->setOptimiserIndex(optimiserIndex); // automatically triggers: optimiserChanged(modelIndex);
        projectControlFrame->setSchemeIndex(schemeIndex); // automatically triggers: fittingSchemeChanged(modelIndex);

    }

    if(status && xmlFileAccessor->isSampleNumberSet()) {
        int sampleNumber = xmlFileAccessor->getSampleNumber();
        visualisationFrame->setFixedSampleNumber(sampleNumber);
        projectControlFrame->setFixedSampleNumber(QString::number(sampleNumber));
    }
    if(status && xmlFileAccessor->areCalibrationParametersSet()) {
        int phantomIndex; QString phantomName;
        xmlFileAccessor->getCalibrationName(phantomIndex, phantomName);
        this->visualisationFrame->setCalibrationPhantom(phantomIndex);
        this->projectControlFrame->setCalibrationPhantom(phantomIndex);

        if( phantomIndex==itk::CorticalBone::kManualControlPtsCal ) {
            vtkSmartPointer<vtkDoubleArray> pts = xmlFileAccessor->getCalibrationPoints();
            double radius; int number;
            xmlFileAccessor->getCalibrationPtGeometry(radius, number);

            this->projectControlFrame->setCalibrationRadius(radius);
            this->projectControlFrame->setCalibrationNumber(number);

            this->visualisationFrame->setCalibrationPtGeometry(radius,number);
            this->visualisationFrame->setCalibrationPoints(phantomIndex, pts);

        } else if(phantomIndex==itk::CorticalBone::kManualLinearCal || phantomIndex==itk::CorticalBone::kManualQuadraticCal ) {
            double p0, p1, p2;
            xmlFileAccessor->getCalibrationParameters(p0, p1, p2);
            this->projectControlFrame->setCalibrationScalingValues(p0, p1, p2);
        } else if(phantomIndex>=itk::CorticalBone::kMindwaySolidCal && phantomIndex<=itk::CorticalBone::kEuropeanSpineCal) {
            vtkSmartPointer<vtkDoubleArray> pts = xmlFileAccessor->getCalibrationPoints();

            if (!this->visualisationFrame->setCalibrationPoints(phantomIndex, pts)) {
                //status = false;
                this->projectControlFrame->removeCalibrationScalingValues();
                std::cout << "error: setting calibration points" << std::endl;
            } else { // update the cal scale displays
                double p0, p1, p2;
                this->visualisationFrame->getCalibrationValues(p0, p1, p2);
                this->projectControlFrame->setCalibrationScalingValues(p0, p1, p2);
            }
        }

        runCalibration(); this->visualisationFrame->disableCalibrationMode();

    }

    int cbConstraintIndex; double fixedCBValue;
    xmlFileAccessor->getCBConstraintInfo(cbConstraintIndex, fixedCBValue);
    if(status && cbConstraintIndex==xmlFileAccess::kCBConstraintFWHM) {

        this->projectControlFrame->setFWHMCBMode();
        this->turnOnFWHMCBDensity();
    } else if(status && cbConstraintIndex==xmlFileAccess::kCBConstraintFixed) {

        this->visualisationFrame->setFixedCBDensity(fixedCBValue);
        this->projectControlFrame->setFixedCBValue(QString::number(fixedCBValue));
    }

    if(status && xmlFileAccessor->isAveragingSet()) {
        double averages[3];
        xmlFileAccessor->getProfileAveraging(averages);
        this->visualisationFrame->turnOnProfileAveraging(averages[0], averages[1], averages[2]);
        this->projectControlFrame->setAveragingValues(QString::number(averages[0]), QString::number(averages[1]), QString::number(averages[2]));
    }

    if(status && xmlFileAccessor->isSmoothingSet()) {
        double smoothingValue = xmlFileAccessor->getSmoothing();
        visualisationFrame->setSmoothingRadius(smoothingValue);
        projectControlFrame->setSmoothingRadius(QString::number(smoothingValue));
        visualisationFrame->turnOnSmoothingMode();
    }

    if(status && xmlFileAccessor->isSigmaSet()) {
        double sigma[3];
        xmlFileAccessor->getSigmaConstraint(sigma);
        this->visualisationFrame->turnOnFixedSigma(sigma[0], sigma[1], sigma[2]);
        this->projectControlFrame->setSigmaValues(QString::number(sigma[0]), QString::number(sigma[1]), QString::number(sigma[2]));
    }

    if(status && xmlFileAccessor->areClassifierThresholdsSet()) { // thresolds already calculated
        double stDensity, cbDensity, thresholdDensity; int classifierIndex;
        xmlFileAccessor->getClassifierThresholds(stDensity, cbDensity, thresholdDensity, classifierIndex);
        if(xmlFileAccessor->isClassifierThresholdWeightSet()) {
            double weight; xmlFileAccessor->getClassifierThresholdWeight(weight);
            this->visualisationFrame->setClassifierThresholdInfo(stDensity, cbDensity, thresholdDensity, weight,
                                                                 classifierIndex);
            projectControlFrame->setThresholdValues(QString::number(stDensity), QString::number(cbDensity), QString::number(thresholdDensity), QString::number(weight), classifierIndex);
        } else {
            this->visualisationFrame->setClassifierThresholdInfo(stDensity, cbDensity, thresholdDensity,
                                                                 classifierIndex);
            this->projectControlFrame->setThresholdValues(QString::number(stDensity), QString::number(cbDensity), QString::number(thresholdDensity), classifierIndex);
        }

    } else if(status && xmlFileAccessor->isClassifierThresholdIndexSet()) { // only threshold level selected, still too calculate
        int classifierIndex; QString classifiedIndexName;
        xmlFileAccessor->getClassifierThresholdIndex(classifierIndex, classifiedIndexName);
        if(xmlFileAccessor->isClassifierThresholdWeightSet()) {
            double weight; xmlFileAccessor->getClassifierThresholdWeight(weight);
            projectControlFrame->setThresholdWeight(QString::number(weight), classifierIndex);
        } else {
            projectControlFrame->setThresholdSelection(classifierIndex);
        }
    }

    if(status && xmlFileAccessor->isImportSet()) {
        openImportParams(xmlFileAccessor->getImportBaseName());
    }

    if(status && meshSet && imageSet) {

        // check for previously save results
        QString baseName = QFileInfo(fileName).absolutePath() + QDir::separator() + QFileInfo(fileName).baseName(); // fileNameStub = fileName.substr(0, fileName.rfind('.'));
        QString extn = QString(".txt");

        if (QFileInfo(baseName+"_Parameters"+extn).exists() && QFileInfo(baseName+"_Displays"+extn).exists() && QFileInfo(baseName+"_ImageProfiles"+extn).exists()) {
            // try load any previously saved results
            if (visualisationFrame->loadValueArrays(fileName)) { // TODO - allow this to be called by the user
                cout << Utilities::getTabString() << "Prevously saved results loaded" << endl;
            }
        }
    }

    // update the state
    updateGUIEntries();
    updateState();
    updateMeasurementEnables();

    if(status) {
        curPath = QFileInfo(fileName).absolutePath();
        curName = QFileInfo(fileName).baseName();
        setCurrentFile(fileName);
        updateStatusBar(tr("Project File loaded"));
    } else {
        updateStatusBar(tr("Error: Unable to Load Project File"));
        std::clog << "Error - Unable to Load Project File" << std::endl;
    }

    return status;

}

bool MainWindow::closeProject() {
    cout<<"closeProject"<<endl;
    // close the file
    projectControlFrame->initialiseSettings();
    visualisationFrame->reset();

    if(!scriptControlFrame->isProcessRunning() && !scriptControlFrame->isViewRunning() ) {
        outputText->clear();
    }

    //updateMeasurementEnables(); //taken care of in the 'initialiseSettings()' and reset() calls
    updateState();
    return true;
}

void MainWindow::setScriptFileName(QString fileName) {
    scriptFilePath = fileName;
}

QString MainWindow::getScriptFileName() {
    return scriptFilePath;
}


//---------------- update states --------//
void MainWindow::setCurrentFile(const QString &fileName) {

    // set name of script-file and indicates when it has been modified
    curFile = fileName;
    scriptText->document()->setModified(false);
    setWindowModified(false);

    QString shownName = curFile;
    if (curFile.isEmpty())
        shownName = "untitled.xml";
    setWindowFilePath(shownName);
}

void MainWindow::updateState() {
    int state = visualisationFrame->getState();
    projectControlFrame->setState(state);
    this->setState(state);
}

void MainWindow:: updateMeasurementEnables() {

    if(this->projectControlFrame->isPtMeasuresEnabled()) {
        this->visualisationFrame->enablePtMeasures();
    } else {
        this->visualisationFrame->disablePtMeasures();
    }

}

void MainWindow::updateGUIEntries() {

    if(this->visualisationFrame->areThresholdsSet()) { // note thresholds cannot be set by the user; thus update only from the CB
        double stDensity, cbDensity, thresholdDensity; int classifierThresholdIndex; QString classifierThresholdName;
        this->visualisationFrame->getClassifierThresholdInfo(stDensity, cbDensity, thresholdDensity, classifierThresholdIndex, classifierThresholdName);
        projectControlFrame->setThresholdValues(QString::number(stDensity), QString::number(cbDensity), QString::number(thresholdDensity), classifierThresholdIndex);
    } else {
        projectControlFrame->removeThresholdValues();
    }
}

void MainWindow::setState(int state) {

    if(scriptControlFrame->isProcessRunning()) {
        scriptControlFrame->setDisabled(true); scriptEditor->setDisabled(true);
        projectControlFrame->setDisabled(true);
        return;
    } else {
        scriptControlFrame->setEnabled(true); scriptEditor->setEnabled(true);
        projectControlFrame->setEnabled(true);
    }

    if (state==0 ) { // image & mesh not open

        openScriptAct->setEnabled(true);
        openImageAct->setEnabled(true);
        openMeshAct->setEnabled(true);
        openProjAct->setEnabled(true);
        openParamsAct->setDisabled(true);
        saveProjAct->setDisabled(true);
        saveMeshAct->setDisabled(true);
        saveDisplayAct->setDisabled(true);
        saveValuesAct->setDisabled(true);
        closeProjAct->setDisabled(true);

        controlTab->setCurrentIndex(1); outputTab->setCurrentIndex(1); //centralTab->setCurrentIndex(1);
        scriptControlFrame->setEnabled(true); scriptEditor->setEnabled(true);

    } else if (state==1) { // image not open; mesh  open

        openScriptAct->setDisabled(true);
        openImageAct->setEnabled(true);
        openMeshAct->setDisabled(true);
        openProjAct->setDisabled(true);
        openParamsAct->setDisabled(true);
        saveProjAct->setEnabled(true);
        saveMeshAct->setEnabled(true);
        saveDisplayAct->setEnabled(true);
        saveValuesAct->setDisabled(true);
        closeProjAct->setEnabled(true);

        controlTab->setCurrentIndex(0); outputTab->setCurrentIndex(0); // centralTab->setCurrentIndex(0);


        if(scriptControlFrame->isViewRunning()) {
            scriptControlFrame->setEnabled(true); scriptEditor->setEnabled(true);
        } else {
            scriptControlFrame->setDisabled(true); scriptEditor->setDisabled(true);
        }

    } else if (state==2) { // if the image open; mesh not open

        openScriptAct->setDisabled(true);
        openImageAct->setDisabled(true);
        openMeshAct->setEnabled(true);
        openProjAct->setDisabled(true);
        openParamsAct->setDisabled(true);
        saveProjAct->setEnabled(true);
        saveMeshAct->setDisabled(true);
        saveDisplayAct->setEnabled(true);
        saveValuesAct->setDisabled(true);
        closeProjAct->setEnabled(true);

        controlTab->setCurrentIndex(0); outputTab->setCurrentIndex(0); // centralTab->setCurrentIndex(0);
        scriptControlFrame->setDisabled(true); scriptEditor->setDisabled(true);

        if(scriptControlFrame->isViewRunning()) {
            scriptControlFrame->setEnabled(true); scriptEditor->setEnabled(true);
        } else {
            scriptControlFrame->setDisabled(true); scriptEditor->setDisabled(true);
        }

    } else if (state==3) { // if the image; mesh open; calibration don't care; measurements not made

        openScriptAct->setDisabled(true);
        openImageAct->setDisabled(true);
        openMeshAct->setDisabled(true);
        openProjAct->setDisabled(true);
        openParamsAct->setEnabled(true);
        saveProjAct->setEnabled(true);
        saveMeshAct->setEnabled(true);
        saveDisplayAct->setEnabled(true); // todo only allow when pt measure has been made
        saveValuesAct->setDisabled(true);
        closeProjAct->setEnabled(true);

        controlTab->setCurrentIndex(0); outputTab->setCurrentIndex(0); // centralTab->setCurrentIndex(0);
        scriptControlFrame->setDisabled(true); scriptEditor->setDisabled(true);

        if(scriptControlFrame->isViewRunning()) {
            scriptControlFrame->setEnabled(true); scriptEditor->setEnabled(true);
        } else {
            scriptControlFrame->setDisabled(true); scriptEditor->setDisabled(true);
        }

    } else if (state==4) { // if the image; mesh open; calibration don't care; measurements made

        openScriptAct->setDisabled(true);
        openImageAct->setDisabled(true);
        openMeshAct->setDisabled(true);
        openProjAct->setDisabled(true);
        openParamsAct->setEnabled(true);
        saveProjAct->setEnabled(true);
        saveMeshAct->setEnabled(true);
        saveDisplayAct->setEnabled(true);
        saveValuesAct->setEnabled(true);
        closeProjAct->setEnabled(true);

        controlTab->setCurrentIndex(0); outputTab->setCurrentIndex(0); // centralTab->setCurrentIndex(0);
        scriptControlFrame->setDisabled(true); scriptEditor->setDisabled(true);

        if(scriptControlFrame->isViewRunning()) {
            scriptControlFrame->setEnabled(true); scriptEditor->setEnabled(true);
        } else {
            scriptControlFrame->setDisabled(true); scriptEditor->setDisabled(true);
        }

    } else {
        std::cout<<"incorrect state selected"<<std::endl;
    }

}

bool MainWindow::updateStatusBar(QString message) {
    statusBar()->showMessage(message, 2000);
    return true;
}


//----------- set up methods --------------//
void MainWindow::createControlFrames() {

    // set up control frames
    projectControlFrame = new QControlFrame(this);
    scriptControlFrame = new QScriptFrame(this, scriptEditor, xmlFileAccessor);

    controlTab = new QTabWidget(this);
    controlTab->addTab(projectControlFrame, tr("Project Controls"));
    controlTab->addTab(scriptControlFrame, tr("Script Controls"));
    controlTab->setCurrentIndex(1); controlTabChanged(1);

    QDockWidget *controlDock = new QDockWidget(tr("Controls"), this);
    controlDock->setStyleSheet("QDockWidget { font: bold }"); // bold title 'Controls'
    controlDock->setAllowedAreas(Qt::LeftDockWidgetArea);
    //QScrollArea *controlScrollFrame; controlScrollFrame->setWidget(projectControlFrame);

    controlDock->setWidget(controlTab);
    addDockWidget(Qt::LeftDockWidgetArea, controlDock);


}

void MainWindow::createXmlFileAccessor() {
    xmlFileAccessor = new xmlFileAccess(true);
}

void MainWindow::createDisplayWorkspace() {

    // set up displays
    visualisationFrame = new QVisualisationFrame(this);
    scriptEditor = new QPlainTextEdit(this);

    QGroupBox* displayBox = new QGroupBox(tr("Display"));
    QHBoxLayout *displayLayout = new QHBoxLayout;
    displayLayout->addWidget(visualisationFrame);
    displayLayout->addWidget(scriptEditor);
    displayBox->setLayout(displayLayout);

    setCentralWidget(displayBox);
}

void MainWindow::createOutputFrame() {

    // create dock - todo consider removing docking capabilities
    QDockWidget *dock = new QDockWidget(this);
    dock->setAllowedAreas(Qt::BottomDockWidgetArea);
    dock->setContentsMargins(0,0,0,0);

    /* create output frame */
    outputTab = new QTabWidget(this);

    outputText = new QPlainTextEdit(outputTab);
    scriptText = new QPlainTextEdit(outputTab);

    scriptText->setReadOnly(true);
    outputText->setReadOnly(true);

    outputTab->addTab(outputText, tr("Project Output"));
    outputTab->addTab(scriptText, tr("Script Output"));

    // add dock to MainWindow
    dock->setWidget(outputTab);
    addDockWidget(Qt::BottomDockWidgetArea, dock);

    // redirect IO - debug
    scriptStream = new QDebugStream (std::clog, scriptText);

    // redirect IO - output
    outputStream = new QDebugStream (std::cout, outputText);

}

void MainWindow::createDockWindows() {

    createXmlFileAccessor();
    createOutputFrame();
    createDisplayWorkspace();
    createControlFrames();

    // change the display when the control tab is changed
    connect(controlTab, SIGNAL(currentChanged(int)), this, SLOT(controlTabChanged(int)));

}

void MainWindow::createActions() {
    openScriptAct = new QAction(QIcon(":/images/openScript.png"), tr("&Open Script"), this);
    //openAct->setShortcuts(QKeySequence::Open);
    openScriptAct->setStatusTip(tr("Open an existing script file"));
    connect(openScriptAct, SIGNAL(triggered()), this, SLOT(openScriptSlot()));

    openImageAct = new QAction(QIcon(":/images/openImage.png"), tr("Open &Image"), this);
    openImageAct->setStatusTip(tr("Open an existing Image directory"));
    connect(openImageAct, SIGNAL(triggered()), this, SLOT(openImageSlot()));

    openMeshAct = new QAction(QIcon(":/images/openMesh.png"), tr("Open &Mesh"), this);
    openMeshAct->setStatusTip(tr("Open an existing Mesh directory"));
    connect(openMeshAct, SIGNAL(triggered()), this, SLOT(openMeshSlot()));

    openProjAct = new QAction(QIcon(":/images/open.png"), tr("&Open"), this);
    openProjAct->setStatusTip(tr("Open an existing Project"));
    connect(openProjAct, SIGNAL(triggered()), this, SLOT(openFileSlot()));

    openParamsAct = new QAction(QIcon(":/images/openParameters.png"), tr("Open &Parameters"), this);
    openParamsAct->setStatusTip(tr("Open model parameters"));
    connect(openParamsAct, SIGNAL(triggered()), this, SLOT(openParamsSlot()));

    closeProjAct = new QAction(QIcon(":/images/close.png"), tr("&Close"), this);
    closeProjAct->setStatusTip(tr("Close the existing Project"));
    connect(closeProjAct, SIGNAL(triggered()), this, SLOT(closeFileSlot()));

    saveProjAct = new QAction(QIcon(":/images/save.png"), tr("&Save"), this);
    saveProjAct->setStatusTip(tr("Save the current Project"));
    connect(saveProjAct, SIGNAL(triggered()), this, SLOT(saveFileSlot()));

    saveMeshAct = new QAction(QIcon(":/images/saveMesh.png"), tr("Save &Mesh"), this);
    saveMeshAct->setStatusTip(tr("Save an existing Mesh directory"));
    connect(saveMeshAct, SIGNAL(triggered()), this, SLOT(saveMeshSlot()));

    saveDisplayAct = new QAction(QIcon(":/images/saveDisplay.png"), tr("Save &Display"), this);
    saveDisplayAct->setStatusTip(tr("Save the current Displays"));
    connect(saveDisplayAct, SIGNAL(triggered()), this, SLOT(saveDisplaySlot()));

    saveValuesAct = new QAction(QIcon(":/images/saveValues.png"), tr("Save &Values"), this);
    saveValuesAct->setStatusTip(tr("Save the current mesh Values"));
    connect(saveValuesAct, SIGNAL(triggered()), this, SLOT(saveValuesSlot()));


    // show displays
    showDisplayAct = new QAction( tr("Show &Display Options"), this);
    showDisplayAct->setCheckable(true); showDisplayAct->setChecked(true);
    showDisplayAct->setStatusTip(tr("Shows the 'Display' options in the Control Frame"));
    connect(showDisplayAct, SIGNAL(toggled(bool)), this, SLOT(showDisplaySlot(bool)));

    showCalibrationAct = new QAction( tr("Show &Calibration Options"), this);
    showCalibrationAct->setCheckable(true); showCalibrationAct->setChecked(true);
    showCalibrationAct->setStatusTip(tr("Shows the 'Calibration' options in the Control Frame"));
    connect(showCalibrationAct, SIGNAL(toggled(bool)), this, SLOT(showCalibrationSlot(bool)));

    showConstraintsAct = new QAction( tr("Show &Constraint Options"), this);
    showConstraintsAct->setCheckable(true); showConstraintsAct->setChecked(true);
    showConstraintsAct->setStatusTip(tr("Shows the 'Constraint' options in the Control Frame"));
    connect(showConstraintsAct, SIGNAL(toggled(bool)), this, SLOT(showConstraintsSlot(bool)));

    showThresholdAct = new QAction( tr("Show &Threshold Options"), this);
    showThresholdAct->setCheckable(true); showThresholdAct->setChecked(false);
    showThresholdAct->setStatusTip(tr("Shows the 'Threshold' options in the Control Frame"));
    connect(showThresholdAct, SIGNAL(toggled(bool)), this, SLOT(showThresholdSlot(bool)));
    showThresholdSlot(false);

    showProfileSettingsAct = new QAction( tr("Show &Profile Settings Options"), this);
    showProfileSettingsAct->setCheckable(true); showProfileSettingsAct->setChecked(true);
    showProfileSettingsAct->setStatusTip(tr("Shows the 'Profile Settings' options in the Control Frame"));
    connect(showProfileSettingsAct, SIGNAL(toggled(bool)), this, SLOT(showProfileSettingsSlot(bool)));

    showImportsAct = new QAction( tr("Show &Import Options"), this);
    showImportsAct->setCheckable(true); showImportsAct->setChecked(true);
    showImportsAct->setStatusTip(tr("Shows the 'Import' options in the Control Frame"));
    connect(showImportsAct, SIGNAL(toggled(bool)), this, SLOT(showImportSlot(bool)));

}

void MainWindow::createMenus() {
    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(openScriptAct);
    fileMenu->addAction(openImageAct);
    fileMenu->addAction(openMeshAct);
    fileMenu->addAction(openProjAct);
    fileMenu->addAction(openParamsAct);
    fileMenu->addAction(closeProjAct);

    fileMenu->addAction(saveProjAct);
    fileMenu->addAction(saveMeshAct);
    fileMenu->addAction(saveDisplayAct);
    fileMenu->addAction(saveValuesAct);

    controlMenu = menuBar()->addMenu(tr("&Controls"));
    controlMenu->addAction(showDisplayAct);
    controlMenu->addAction(showCalibrationAct);
    controlMenu->addAction(showConstraintsAct);
    controlMenu->addAction(showThresholdAct);
    controlMenu->addAction(showProfileSettingsAct);
    controlMenu->addAction(showImportsAct);
}

void MainWindow::createToolBars() {
    fileToolBar = addToolBar(tr("File"));
    fileToolBar->addAction(openScriptAct);
    fileToolBar->addAction(openImageAct);
    fileToolBar->addAction(openMeshAct);
    fileToolBar->addAction(openProjAct);
    fileToolBar->addAction(openParamsAct);
    fileToolBar->addAction(closeProjAct);

    fileToolBar->addAction(saveProjAct);
    fileToolBar->addAction(saveMeshAct);
    fileToolBar->addAction(saveDisplayAct);
    fileToolBar->addAction(saveValuesAct);
}

void MainWindow::createStatusBar() {
    statusBar()->showMessage(tr("Ready"));
}

void MainWindow::readSettings() {

    // todo - set location of the .ini file
    QSettings settings("QtProject", "Application Example");
    QPoint pos = settings.value("pos", QPoint(200, 200)).toPoint();
    QSize size = settings.value("size", QSize(400, 400)).toSize();

    QString defaultPath =  QDir::currentPath();
    curPath = settings.value("curPath", defaultPath).toString();
    imageFilePath = settings.value("imagePath", defaultPath).toString();
    meshFilePath = settings.value("meshPath", defaultPath).toString();
    scriptFilePath = settings.value("scriptFilePath", defaultPath).toString();

    importParamsPath = settings.value("importParamsPath", defaultPath).toString();

    QString scriptPath = settings.value("scriptFilePath", defaultPath).toString();
    scriptControlFrame->setProjectPath(scriptPath);

    resize(size);
    move(pos);
}

void MainWindow::writeSettings() {
    QSettings settings("QtProject", "Application Example");
    settings.setValue("pos", pos());
    settings.setValue("size", size());
    settings.setValue("curPath", curPath);
    settings.setValue("imagePath", imageFilePath);
    settings.setValue("meshPath", meshFilePath);
    settings.setValue("importParamsPath", importParamsPath);
    settings.setValue("scriptFilePath", scriptFilePath);

}
