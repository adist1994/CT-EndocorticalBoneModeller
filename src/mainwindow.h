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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class QAction;
class QMenu;
class QPlainTextEdit;
class QFrame;
class QPushButton;
class QTabWidget;
class QListWidget;
class QDebugStream;
class QScriptFrame;
class QVisualisationFrame;
class QControlFrame;
class ProjectParser;
class xmlFileAccess;

class MainWindow : public QMainWindow
{
Q_OBJECT

public:
    MainWindow();
    void saveResults(QString fileName);
    void saveDisplay(QString fileName);
    bool openProject(QString fileName);
    bool saveProject(QString fileName);
    void setScriptFileName(QString fileName);
    QString getScriptFileName();
    bool closeProject();
    bool updateStatusBar(QString message);
    QPlainTextEdit* getScriptOutputEditor();
    void setThresholdIndex(int percentIndex); // allows script to set the value in control frame

protected:
    void closeEvent(QCloseEvent *event);

private slots:
    void openScriptSlot();
    void openImageSlot();
    void openMeshSlot();
    void openFileSlot();
    void openParamsSlot();
    void saveFileSlot();
    void saveDisplaySlot();
    void saveValuesSlot(); // values
    void saveMeshSlot();
    void closeFileSlot();

    void controlTabChanged(int index);

    void showDisplaySlot(bool state);
    void showCalibrationSlot(bool state);
    void showConstraintsSlot(bool state);
    void showThresholdSlot(bool state);
    void showProfileSettingsSlot(bool state);
    void showImportSlot(bool state);


public slots:
    void documentWasModified(bool modified);  // todo - use
    // measurements
    void runModellingOverMesh();
    void togglePtMeasures(bool state);
    void disablePtMeasures();
    void updatePtMeasures();
    void modelFunctionChanged(int index);
    void optimiserChanged(int index);
    void fittingSchemeChanged(int index);
    // display
    void displayParameterChanged(int index);
    void displayMeshChanged(int index);
    // calibration
    void runCalibration();
    void toggleCalibrationMode(bool state);
    void calibrationPhantomChanged(int index);
    void updateCalibrationRadius();
    void updateCalibrationNumber();
    void disableCalibration();
    // constrained CB Density
    void toggleCBDensityConstraint(bool cbState);
    void turnOnFixCBDensity();
    void turnOnFWHMCBDensity();
    void fixedCBDensityChanged();
    // constrained sigma
    void toggleFixedSigma(bool state);
    void sigmaChanged();
    // thresholds
    void thresholdSelectorChanged(int index);
    void runThresholdCalculation();
    // profile settings - Sample number
    void toggleSampleNumber(bool state);
    void sampleNumberChanged();
    // profile settings - smoothing
    void smoothingChanged();
    void toggleSmoothing(bool state);
    // profile settings - Profile Averaging
    void ProfileAveragingValuesChanged();
    void toggleProfileAverging(bool state);
    // Import Profiles / Parameters
    void importProfile();
    void removeImportedProfile();

protected:
    void keyPressEvent(QKeyEvent *);

private:

    // high level setup
    void createDockWindows();
    void createActions();
    void createMenus();
    void createToolBars();
    void createStatusBar();

    void readSettings();
    void writeSettings();

    // updates
    void updateState();
    void updateMeasurementEnables();
    void updateGUIEntries();

    // tasks in response to user inputs   
    bool openMesh(QString fileName);
    bool openImage(QString fileName);
    bool openImportParams(QString fileName);

    // todo - change to indicate if general state has been changed 
    void setCurrentFile(const QString &fileName);
    void setState(int state);

    // control panel 
    void createControlFrames();
    QControlFrame *projectControlFrame;
    QScriptFrame *scriptControlFrame;
    QTabWidget *controlTab;

    // visualisation and script panel
    void createDisplayWorkspace();
    QVisualisationFrame *visualisationFrame;
    QPlainTextEdit *scriptEditor;

    void createXmlFileAccessor();
    xmlFileAccess *xmlFileAccessor;

    // output panel
    void createOutputFrame();
    QTabWidget *outputTab;
    QPlainTextEdit *scriptText;
    QPlainTextEdit *outputText;
    QDebugStream *scriptStream;
    QDebugStream *outputStream;

    // file paths and other setting values - todo wrap up in some data structure
    QString curFile;
    QString curPath;

    QString imageFilePath;
    QString meshFilePath;
    QString importParamsPath;
    QString scriptFilePath;
    QString curName;

    QMenu *fileMenu;
    QToolBar *fileToolBar;
    QAction *openScriptAct;
    QAction *openImageAct;
    QAction *openMeshAct;
    QAction *openProjAct;
    QAction *openParamsAct;
    QAction *saveProjAct;
    QAction *saveMeshAct;
    QAction *saveValuesAct;
    QAction *saveDisplayAct;
    QAction *closeProjAct;

    QMenu *controlMenu;
    QAction *showDisplayAct;
    QAction *showCalibrationAct;
    QAction *showConstraintsAct;
    QAction *showThresholdAct;
    QAction *showProfileSettingsAct;
    QAction *showImportsAct;

    bool upToDate;
};

#endif
