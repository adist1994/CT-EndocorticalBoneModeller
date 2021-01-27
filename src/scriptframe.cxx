#include <QtWidgets>
#include <iostream>
#include <istream>
#include <sstream>

#include "scriptframe.h"
#include "xmlFileAccess.h"
#include "mainwindow.h"


QScriptFrame::QScriptFrame(MainWindow *mainWindow, QPlainTextEdit* scriptEditorIn, xmlFileAccess *xmlFileWriterIn) : QFrame(mainWindow) {

    parentObject = mainWindow;
    scriptEditor = scriptEditorIn;
    xmlFileAccessor = xmlFileWriterIn;

    createObjects();
    createConnections();

    projectPath = QString();
    reset();
}

QString QScriptFrame::getProjectPath() {
    return projectPath;
}

void QScriptFrame::setProjectPath(QString filePath) {
    projectPath = filePath;
}

bool QScriptFrame::isViewRunning() {
    return runningViewScript;
}

bool QScriptFrame::isProcessRunning() {
    return runningProcessScript;
}

// private functions
bool QScriptFrame::reset() {

    //- clear output display
    parentObject->getScriptOutputEditor()->clear();

    //--- state
    runningProcessScript = runningViewScript= false;
    projectDirectorySet = projectFilterSet = false;
    imageFilterSet = meshFilterSet = false;
    importFileSet = false;
    tabsUpToDate = false;
    sampleNumberSet = smoothingSet = fixedCBFilterSet = xAvgSet = yAvgSet = zAvgSet = false;
    includeSampleNumber=includeSmoothing=includeAveraging=false;
    includeCalibration = p0ScaleSet = p1ScaleSet = p2ScaleSet = false;
    calibrationIsInFile = pFileSet = false;
    stThresholdSet = cbThresholdSet = thresholdSet = includeThresholds = includeThresholdWeight = false;
    cbConstraintIndex=xmlFileAccess::kCBConstraintNone; fixedCBSet=fixedCBIsNumber=false; fixedCBValue=nan("1");
    viewIteratorIndex = -1;
    sigmaRequired = sigmaXSet = sigmaYSet = sigmaZSet = false; sigmaUpToDate = true;

    projectExtnString=tr(".xml"); txtExtnString=tr(".txt");
    meshExtnString=tr(".obj"); imageExtnString=tr(".dcm");

    projectNameList = imageNameList = meshNameList = fixedCBFileList = importedFileList = calibrationFileList = QStringList();

    projectFilterString=QString();
    imageFilterString=meshFilterString=fixedCBFilterString=QString();
    sampleNumberInt=-1;
    profileAverages[0]=profileAverages[1]=profileAverages[2]=-1;
    calibrationScales[0]= calibrationScales[0]= calibrationScales[0]=-1;
    thresholds[0]= thresholds[0]= thresholds[0]=nan("1");
    thresholdWeight=nan("1");

    // set object states
    scriptDirectoryEdit->setText(tr(""));
    cbDensityButton->setChecked(false); profileAveragingButton->setChecked(false); sigmaFixedButton->setChecked(false);
    imageExtnSelector->setCurrentIndex(0); meshExtnSelector->setCurrentIndex(0);
    projectFilterEdit->setText(projectFilterString);
    imageFilterEdit->setText(imageFilterString); meshFilterEdit->setText(meshFilterString);

    cbDensityFixedEdit->setText(QString()); importFileEdit->setText(QString());
    xAvgEdit->setText(QString()); yAvgEdit->setText(QString()); zAvgEdit->setText(QString());
    p0ScaleEdit->setText(QString()); p1ScaleEdit->setText(QString()); p2ScaleEdit->setText(QString());
    pFileEdit->setText(QString());
    stThreshEdit->setText(QString()); cbThreshEdit->setText(QString()); thresholdEdit->setText(QString());
    thresholdWeightEdit->setText(QString());
    sigmaXEdit->setText(QString()); sigmaYEdit->setText(QString()); sigmaYEdit->setText(QString());

    cbConstraintUpToDate=true;

    modelTypeSelector->setCurrentIndex(0); schemeTypeSelector->setCurrentIndex(0);
    optimiserTypeLabel->show(); optimiserTypeSelector->show();
    classifierLevelLabel->hide(); thresholdLevelSelector->hide();

    xmlFileAccessor->reset();

    updateObjectEnables();
    return true;
}

void QScriptFrame::createObjects() {

    // script editor
    //scriptEditor = new QPlainTextEdit(this);
    scriptEditor->setReadOnly(true);
    QFont font; // font.setFamily("Courier"); font.setStyleHint(QFont::Monospace); font.setFixedPitch(true); font.setPointSize(10);
    QFontMetrics metrics(font); const int tabWidth = 4;  // 4 characters
    scriptEditor->setTabStopWidth(tabWidth * metrics.width(' '));

    // script directory
    scriptDirectoryButton = new QPushButton(tr("Script Directory"));
    scriptDirectoryEdit = new QLineEdit;
    scriptDirectoryEdit->setEnabled(false); scriptDirectoryEdit->setReadOnly(true);

    projectFilterLabel = new QLabel(tr("Project Name Filter"));
    projectFilterEdit = new QLineEdit;

    QGroupBox *projectNameBox = new QGroupBox();
    QGridLayout *projectHLayout = new QGridLayout;
    projectHLayout->addWidget(scriptDirectoryButton, 0, 0);
    projectHLayout->addWidget(scriptDirectoryEdit, 0, 1);
    projectHLayout->addWidget(projectFilterLabel, 1, 0);
    projectHLayout->addWidget(projectFilterEdit, 1, 1);
    projectNameBox->setLayout(projectHLayout);
    projectNameBox->setContentsMargins(0,0,0,0);

    // create tab
    createProjectsBox = new QGroupBox();

    // create - string tags
    imageFilterLabel = new QLabel(tr("Image Filter"));
    meshFilterLabel = new QLabel(tr("Mesh Filter"));
    imageFilterEdit = new QLineEdit;
    meshFilterEdit = new QLineEdit;

    imageExtnLabel = new QLabel(tr("Image Extn"));
    imageExtnSelector = new QComboBox(this);
    imageExtnSelector->addItem("DICOM");
    imageExtnSelector->addItem("RAW");
    imageExtnSelector->addItem("TIFF");
    imageExtnSelector->addItem("QCT");

    meshExtnLabel = new QLabel(tr("Mesh Extn"));
    meshExtnSelector = new QComboBox(this);
    meshExtnSelector->addItem("OBJ");
    meshExtnSelector->addItem("VRML");
    meshExtnSelector->addItem("STL");
    meshExtnSelector->addItem("PLY");

    sampleNumberButton = new QCheckBox(tr("Sample Number"));
    sampleNumberEdit = new QLineEdit;
    sampleNumberEdit->setValidator( new QIntValidator(0, 1000, this) );

    // averaging
    profileAveragingButton = new QCheckBox(tr("Averaging"));

    xAvgEdit = new QLineEdit; yAvgEdit = new QLineEdit; zAvgEdit = new QLineEdit;
    xAvgEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );
    yAvgEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );
    zAvgEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );

    xAvgLabel = new QLabel(tr("X")); yAvgLabel = new QLabel(tr("Y")); zAvgLabel = new QLabel(tr("Z"));

    QGroupBox *createTopBox = new QGroupBox(); // no text
    QGridLayout *createTopLayout = new QGridLayout;
    createTopLayout->addWidget(imageFilterLabel, 0, 0);
    createTopLayout->addWidget(imageFilterEdit, 0, 1, 1, 6); // row, col, row span, col span
    createTopLayout->addWidget(meshFilterLabel, 1, 0);
    createTopLayout->addWidget(meshFilterEdit, 1, 1, 1, 6);
    createTopLayout->addWidget(imageExtnLabel, 2, 0);
    createTopLayout->addWidget(imageExtnSelector, 2, 1, 1, 6);
    createTopLayout->addWidget(meshExtnLabel, 3, 0);
    createTopLayout->addWidget(meshExtnSelector, 3, 1, 1, 6);
    createTopLayout->addWidget(sampleNumberButton, 4, 0);
    createTopLayout->addWidget(sampleNumberEdit, 4, 2, 1, 5);

    createTopLayout->addWidget(profileAveragingButton, 5, 0);
    createTopLayout->addWidget(xAvgLabel, 5, 1); createTopLayout->addWidget(xAvgEdit, 5, 2);
    createTopLayout->addWidget(yAvgLabel, 5, 3); createTopLayout->addWidget(yAvgEdit, 5, 4);
    createTopLayout->addWidget(zAvgLabel, 5, 5); createTopLayout->addWidget(zAvgEdit, 5, 6);
    createTopBox->setLayout(createTopLayout);
    createTopBox->setContentsMargins(0,0,0,0);

    // constraints

    // cb
    cbDensityButton = new QCheckBox(); cbDensityLabel = new QLabel(tr("Y<sub>CB</sub>"));
    cbDensityFixedButton = new QRadioButton(tr("Fixed"));
    cbDensityFWHMButton = new QRadioButton(tr("FWHM"));
    cbDensityFixedEdit = new QLineEdit; // tag of file to parse of max cb value

    // sigma
    sigmaFixedButton = new QCheckBox(tr("σ"));
    sigmaXLabel = new QLabel(tr("X")); sigmaYLabel = new QLabel(tr("Y")); sigmaZLabel = new QLabel(tr("Z"));
    sigmaXEdit = new QLineEdit; sigmaYEdit = new QLineEdit; sigmaZEdit = new QLineEdit;
    sigmaXEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );
    sigmaYEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );
    sigmaZEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );

    QGroupBox *constraintsBox = new QGroupBox(tr("Model Constraints"));
    QGridLayout *constraintsLayout = new QGridLayout;
    constraintsLayout->addWidget(cbDensityButton, 0,0);
    constraintsLayout->addWidget(cbDensityLabel, 0,1);
    constraintsLayout->addWidget(cbDensityFWHMButton, 0,2, 1,2);
    constraintsLayout->addWidget(cbDensityFixedButton, 0,4, 1,2);
    constraintsLayout->addWidget(cbDensityFixedEdit, 0,6, 1,2);
    constraintsLayout->addWidget(sigmaFixedButton, 1,0, 1,2);
    constraintsLayout->addWidget(sigmaXLabel, 1,2); constraintsLayout->addWidget(sigmaXEdit, 1,3);
    constraintsLayout->addWidget(sigmaYLabel, 1,4); constraintsLayout->addWidget(sigmaYEdit, 1,5);
    constraintsLayout->addWidget(sigmaZLabel, 1,6); constraintsLayout->addWidget(sigmaZEdit, 1,7);
    constraintsBox->setLayout(constraintsLayout);


    // import files
    ImportFileButton = new QPushButton(tr("Import Map"));
    importFileEdit = new QLineEdit;
    importFileEdit->setEnabled(false); importFileEdit->setReadOnly(true);


    QGroupBox *importBox = new QGroupBox(tr("Importing"));
    QGridLayout *importLayout = new QGridLayout;
    importLayout->addWidget(ImportFileButton, 0, 0, 1, 1); // row, column
    importLayout->addWidget(importFileEdit, 0, 1, 1, 6);
    importBox->setLayout(importLayout);

    // overall layout
    QVBoxLayout *createLayout = new QVBoxLayout;
    createLayout->addWidget(createTopBox);
    createLayout->addWidget(constraintsBox);
    createLayout->addWidget(importBox);
    createProjectsBox->setLayout(createLayout);

    //- view
    viewProjectsBox = new QGroupBox();

    // calibrate files
    calibrationTypeBox = new QComboBox(this);
    calibrationTypeBox->addItem("No Phantom");
    calibrationTypeBox->addItem("Mindways Solid Phantom");
    calibrationTypeBox->addItem("Bone Density Phantom");
    calibrationTypeBox->addItem("European Spine Phantom");
    calibrationTypeBox->addItem("Manual Control Points");
    calibrationTypeBox->addItem("Manual Linear");
    calibrationTypeBox->addItem("Manual Quadratic");
    calibrationTypeBox->addItem("Parameter File");
    calibrationTypeBox->setItemData(1, 0, Qt::UserRole - 1); // dirty hack - to disable todo remove
    calibrationTypeBox->setItemData(2, 0, Qt::UserRole - 1); // dirty hack - to disable todo remove
    calibrationTypeBox->setItemData(3, 0, Qt::UserRole - 1); // dirty hack - to disable todo remove


    p0ScaleEdit = new QLineEdit; p0ScaleEdit->setValidator( new QDoubleValidator(-10000, 10000, 9, this) );
    p1ScaleEdit = new QLineEdit; p1ScaleEdit->setValidator( new QDoubleValidator(-10000, 10000, 9, this) );
    p2ScaleEdit = new QLineEdit; p2ScaleEdit->setValidator( new QDoubleValidator(-10000, 10000, 9, this) );

    p0ScaleLabel = new QLabel(tr("P0")); p1ScaleLabel = new QLabel(tr("P1")); p2ScaleLabel = new QLabel(tr("P2"));

    pFileEdit = new QLineEdit;

    pFileLabel = new QLabel(tr("File Name"));

    QGroupBox *calibrationBox = new QGroupBox(tr("Calibration"));
    QGridLayout *calibrationLayout = new QGridLayout;
    calibrationLayout->addWidget(calibrationTypeBox, 0, 0, 1, 3); // row, column
    calibrationLayout->addWidget(p0ScaleLabel, 1, 0); calibrationLayout->addWidget(p0ScaleEdit, 1, 1);
    calibrationLayout->addWidget(p1ScaleLabel, 1, 2); calibrationLayout->addWidget(p1ScaleEdit, 1, 3);
    calibrationLayout->addWidget(p2ScaleLabel, 1, 4); calibrationLayout->addWidget(p2ScaleEdit, 1, 5);
    calibrationLayout->addWidget(pFileLabel, 2, 0, 1, 3); calibrationLayout->addWidget(pFileEdit, 2, 3, 1, 3);
    calibrationBox->setLayout(calibrationLayout);

    QVBoxLayout *viewLayout = new QVBoxLayout;
    viewLayout->addWidget(calibrationBox);
    viewProjectsBox->setLayout(viewLayout);

    //- process
    processProjectsBox = new QGroupBox();

    modelTypeLabel = new QLabel(tr("Model Type"));
    modelTypeSelector = new QComboBox(this);
    modelTypeSelector->addItem("Three Tier Rect");
    modelTypeSelector->addItem("Endosteal Ramp");
    modelTypeSelector->addItem("HR Classifier");
    modelTypeSelector->addItem("Calibration");

    optimiserTypeLabel = new QLabel(tr("Optimiser Type"));
    optimiserTypeSelector = new QComboBox(this);
    optimiserTypeSelector->addItem("LM Optimiser");
    optimiserTypeSelector->addItem("Powell Optimiser");
    optimiserTypeSelector->addItem("Evolutionary Optimiser");

    classifierLevelLabel = new QLabel(tr("Classifier Threshold Level"));
    thresholdLevelSelector = new QComboBox(this);
    thresholdLevelSelector->addItem("Very High Threshold [99.99%]");
    thresholdLevelSelector->addItem("High Threshold [99.90%]");
    thresholdLevelSelector->addItem("Medium Threshold [99.0%]");
    thresholdLevelSelector->addItem("Low Threshold [95.0%]");
    thresholdLevelSelector->addItem("Manually Set");
    thresholdLevelSelector->addItem("Median midpoint");
    thresholdLevelSelector->addItem("Median manual");

    thresholdWeightingLabel = new QLabel(tr("ρ<sub>st</sub>"));
    thresholdWeightEdit = new QLineEdit; thresholdWeightEdit->setValidator( new QDoubleValidator(0, 1, 6, this) );

    thresholdWeightBox = new QGroupBox();
    QHBoxLayout *thresholdWeightLayout = new QHBoxLayout;
    thresholdWeightLayout->addWidget(thresholdWeightingLabel); thresholdWeightLayout->addWidget(thresholdWeightEdit);
    thresholdWeightBox->setLayout(thresholdWeightLayout);

    stThreshLabel = new QLabel(tr("ρ<sub>st</sub>")); cbThreshLabel = new QLabel(tr("ρ<sub>cb</sub>")); thresholdLabel = new QLabel(tr("ρ<sub>thresh</sub>"));
    stThreshEdit = new QLineEdit; stThreshEdit->setValidator( new QDoubleValidator(-20000, 20000, 6, this) );
    cbThreshEdit = new QLineEdit; cbThreshEdit->setValidator( new QDoubleValidator(20000, 20000, 6, this) );
    thresholdEdit = new QLineEdit; thresholdEdit->setValidator( new QDoubleValidator(20000, 20000, 6, this) );

    thresholdBox = new QGroupBox();
    QHBoxLayout *thresholdLayout = new QHBoxLayout;
    thresholdLayout->addWidget(stThreshLabel); thresholdLayout->addWidget(stThreshEdit);
    thresholdLayout->addWidget(cbThreshLabel); thresholdLayout->addWidget(cbThreshEdit);
    thresholdLayout->addWidget(thresholdLabel); thresholdLayout->addWidget(thresholdEdit);
    thresholdBox->setLayout(thresholdLayout);

    schemeTypeLabel = new QLabel(tr("Fitting Scheme"));
    schemeTypeSelector = new QComboBox(this);
    schemeTypeSelector->addItem("std A"); // rect 1st then ramp if ramp selected - huerstic weighting
    schemeTypeSelector->addItem("std B"); // rect 1st then ramp if ramp selected - narrow weighting
    schemeTypeSelector->addItem("std C"); // rect 1st then ramp if ramp selected - wide weighting
    schemeTypeSelector->addItem("std D"); // rect 1st then ramp if ramp selected - narrow to wide weighting
    schemeTypeSelector->addItem("CBMV2a"); // (selected model fit 2x + sigma correction) - static huristic weighting
    schemeTypeSelector->addItem("CBMV2b"); // (selected model fit 2x + sigma correction) - static narrow weighting
    schemeTypeSelector->addItem("CBMV2c"); // (selected model fit 2x + sigma correction) - static narrow then wide weighting
    schemeTypeSelector->addItem("CBMV2d"); // (selected model fit 2x + sigma correction) - dynamic narrow then wide weighting
    schemeTypeSelector->addItem("Unconstrained a"); // rect followed by ramp/rect followed by ramp/rect no constraints
    schemeTypeSelector->addItem("Unconstrained b");
    schemeTypeSelector->addItem("Unconstrained c");
    schemeTypeSelector->addItem("Smoothing a"); // V4 then smooth the CB estimates then V1
    schemeTypeSelector->addItem("Smoothing b");
    schemeTypeSelector->addItem("Smoothing c");
    schemeTypeSelector->addItem("Smoothing d");

    smoothingButton = new QCheckBox(tr("Smoothing Value"));
    smoothingEdit = new QLineEdit;
    smoothingEdit->setValidator( new QDoubleValidator(0, 10, 3, this) );

    rangeLabel = new QLabel(tr("Processing Range"));
    startSpinBox = new QSpinBox(this);
    startSpinBox->setRange(0, 0);
    startSpinBox->setSingleStep(1);
    startSpinBox->setValue(0);
    stopSpinBox = new QSpinBox(this);
    stopSpinBox->setRange(0, 0);
    stopSpinBox->setSingleStep(1);
    stopSpinBox->setValue(0);


    // formatting
    QGridLayout *processingLayout = new QGridLayout;
    processingLayout->addWidget(modelTypeLabel, 0, 0); // row, column
    processingLayout->addWidget(modelTypeSelector, 0, 1, 1,2);
    processingLayout->addWidget(optimiserTypeLabel, 1, 0);
    processingLayout->addWidget(optimiserTypeSelector, 1, 1, 1,2);
    processingLayout->addWidget(classifierLevelLabel, 2, 0);
    processingLayout->addWidget(thresholdLevelSelector, 2, 1, 1, 2);
    processingLayout->addWidget(thresholdWeightBox, 3,0,1,3);
    processingLayout->addWidget(thresholdBox, 4,0,1,3);
    processingLayout->addWidget(schemeTypeLabel, 5, 0);
    processingLayout->addWidget(schemeTypeSelector, 5, 1, 1,2);
    processingLayout->addWidget(smoothingButton, 6, 0);
    processingLayout->addWidget(smoothingEdit, 6, 1, 1,2);
    processingLayout->addWidget(rangeLabel, 7, 0);
    processingLayout->addWidget(startSpinBox, 7, 1);
    processingLayout->addWidget(stopSpinBox, 7, 2);
    processProjectsBox->setLayout(processingLayout);

    // tab controls
    scriptControlsTab = new QTabWidget(this);
    scriptControlsTab->addTab(createProjectsBox, tr("Create"));
    scriptControlsTab->addTab(viewProjectsBox, tr("View"));
    scriptControlsTab->addTab(processProjectsBox, tr("Process"));

    //--- overall controls
    loadButton = new QPushButton(tr("Load"));
    saveButton = new QPushButton(tr("Save"));
    runButton = new QPushButton(tr("Run"));
    nextButton = new QPushButton(tr("Next"));


    // formatting
    QGroupBox *buttonBox = new QGroupBox();
    QHBoxLayout *buttonLayout = new QHBoxLayout();
    buttonLayout->addWidget(loadButton);
    buttonLayout->addWidget(saveButton);
    buttonLayout->addWidget(runButton);
    buttonLayout->addWidget(nextButton);
    buttonBox->setLayout(buttonLayout);

    //--- control formatting
    QGroupBox *scriptControlBox = new QGroupBox(tr("Script Controls"));
    QVBoxLayout *scriptControlLayout = new QVBoxLayout();
    //scriptControlLayout->addWidget(projectDirectoryButton);
    scriptControlLayout->addWidget(projectNameBox);
    scriptControlLayout->addWidget(scriptControlsTab);
    scriptControlLayout->addWidget(buttonBox);
    scriptControlBox->setLayout(scriptControlLayout);

    //----------- Overall Formatting
    QGridLayout *layout = new QGridLayout(this);
    layout->addWidget(scriptControlBox, 0,0); //layout->setColumnStretch(1, 10);
    //layout->addWidget(scriptEditor, 0,1);  //layout->setColumnStretch(1, 15);
    setLayout(layout);
}

void QScriptFrame::createConnections() {

    // general directory and filtering
    connect(scriptDirectoryButton, SIGNAL(clicked()), this, SLOT(setDirectory()));
    connect(projectFilterEdit, SIGNAL(returnPressed()), this, SLOT(projectFilterChanged()));
    connect(projectFilterEdit, SIGNAL(textChanged(const QString &)), this, SLOT(projectFilterChanging()));

    // general buttons
    connect(runButton, SIGNAL(clicked()), this, SLOT(runScript()));
    connect(saveButton, SIGNAL(clicked()), this, SLOT(saveScript()));
    connect(loadButton, SIGNAL(clicked()), this, SLOT(loadScript()));
    connect(nextButton, SIGNAL(clicked()), this, SLOT(nextProjectFile()));

    // general tans
    connect(scriptControlsTab, SIGNAL(currentChanged(int)), this, SLOT(tabChanged()));

    // create tab
    connect(imageFilterEdit, SIGNAL(returnPressed()), this, SLOT(imageFilterChanged()));
    connect(imageFilterEdit, SIGNAL(textChanged(const QString &)), this, SLOT(imageFilterChanging()));
    connect(meshFilterEdit, SIGNAL(returnPressed()), this, SLOT(meshFilterChanged()));
    connect(meshFilterEdit, SIGNAL(textChanged(const QString &)), this, SLOT(meshFilterChanging()));
    connect(imageExtnSelector, SIGNAL(currentIndexChanged(int)), this, SLOT(imageExtnChanged(int)));
    connect(meshExtnSelector, SIGNAL(currentIndexChanged(int)), this, SLOT(meshExtnChanged(int)));

    connect(sampleNumberEdit, SIGNAL(returnPressed()), this, SLOT(sampleNumberChanged()));
    connect(sampleNumberEdit, SIGNAL(textChanged(const QString &)), this, SLOT(sampleNumberChanging()));
    connect(sampleNumberButton, SIGNAL(clicked(bool)), this, SLOT(toggleSampleNumber(bool)));

    // constraints
    // cb
    connect(cbDensityButton, SIGNAL(clicked(bool)), this, SLOT(toggleCBConstraintState(bool)));
    connect(cbDensityFixedEdit, SIGNAL(returnPressed()), this, SLOT(cbValueFilterChanged()));
    connect(cbDensityFixedEdit, SIGNAL(textChanged(const QString &)), this, SLOT(cbValueFilterChanging()));
    connect(cbDensityFixedButton, SIGNAL(clicked(bool)), this, SLOT(toggleFixedCB(bool)));
    connect(cbDensityFWHMButton, SIGNAL(clicked(bool)), this, SLOT(toggleFWHMCB(bool)));
    // sigma
    connect(sigmaFixedButton, SIGNAL(clicked(bool)), this, SLOT(toggleFixedSigma(bool)));
    connect(sigmaXEdit , SIGNAL(returnPressed()),this, SLOT(fixedSigmaXChanged()));
    connect(sigmaXEdit, SIGNAL(textChanged(const QString &)), this, SLOT(fixedSigmaXChanging()));
    connect(sigmaYEdit , SIGNAL(returnPressed()),this, SLOT(fixedSigmaYChanged()));
    connect(sigmaYEdit, SIGNAL(textChanged(const QString &)), this, SLOT(fixedSigmaYChanging()));
    connect(sigmaZEdit , SIGNAL(returnPressed()),this, SLOT(fixedSigmaZChanged()));
    connect(sigmaZEdit, SIGNAL(textChanged(const QString &)), this, SLOT(fixedSigmaZChanging()));


    // averaging
    connect(xAvgEdit, SIGNAL(returnPressed()), this, SLOT(setXAveragingValue()));
    connect(xAvgEdit, SIGNAL(textChanged(const QString &)), this, SLOT(xAveragingChanging()));
    connect(yAvgEdit, SIGNAL(returnPressed()), this, SLOT(setYAveragingValue()));
    connect(yAvgEdit, SIGNAL(textChanged(const QString &)), this, SLOT(yAveragingChanging()));
    connect(zAvgEdit, SIGNAL(returnPressed()), this, SLOT(setZAveragingValue()));
    connect(zAvgEdit, SIGNAL(textChanged(const QString &)), this, SLOT(zAveragingChanging()));
    connect(profileAveragingButton, SIGNAL(clicked(bool)), this, SLOT(toggleAveraging(bool)));

    connect(ImportFileButton, SIGNAL(clicked()), this, SLOT(setImportFile()));


    // view
    connect(p0ScaleEdit, SIGNAL(returnPressed()), this, SLOT(setp0ScaleValue()));
    connect(p0ScaleEdit, SIGNAL(textChanged(const QString &)), this, SLOT(p0ScaleChanging()));
    connect(p1ScaleEdit, SIGNAL(returnPressed()), this, SLOT(setp1ScaleValue()));
    connect(p1ScaleEdit, SIGNAL(textChanged(const QString &)), this, SLOT(p1ScaleChanging()));
    connect(p2ScaleEdit, SIGNAL(returnPressed()), this, SLOT(setp2ScaleValue()));
    connect(p2ScaleEdit, SIGNAL(textChanged(const QString &)), this, SLOT(p2ScaleChanging()));
    connect(pFileEdit, SIGNAL(returnPressed()), this, SLOT(setpFileTextChanged()));
    connect(pFileEdit, SIGNAL(textChanged(const QString &)), this, SLOT(pFileTextChanging()));
    connect(calibrationTypeBox, SIGNAL(currentIndexChanged(int)), this, SLOT(calibrationTypeChanged(int)));




    // process tab
    connect(modelTypeSelector, SIGNAL(currentIndexChanged(int)),this,SLOT(modelSelectionChanged(int)));
    connect(startSpinBox, SIGNAL(valueChanged(int)),this,SLOT(processMinRangeChanged(int)));
    connect(stopSpinBox, SIGNAL(valueChanged(int)),this,SLOT(processMaxRangeChanged(int)));

    connect(schemeTypeSelector, SIGNAL(currentIndexChanged(int)),this,SLOT(schemeChanged(int)));

    connect(smoothingEdit, SIGNAL(returnPressed()), this, SLOT(smoothingChanged()));
    connect(smoothingEdit, SIGNAL(textChanged(const QString &)), this, SLOT(smoothingChanging()));
    connect(smoothingButton, SIGNAL(clicked(bool)), this, SLOT(toggleSmoothing(bool)));

    connect(stThreshEdit, SIGNAL(returnPressed()), this, SLOT(stThresholdChanged()));
    connect(stThreshEdit, SIGNAL(textChanged(const QString &)), this, SLOT(stThresholdChanging()));
    connect(cbThreshEdit, SIGNAL(returnPressed()), this, SLOT(cbThresholdChanged()));
    connect(cbThreshEdit, SIGNAL(textChanged(const QString &)), this, SLOT(cbThresholdChanging()));
    connect(thresholdEdit, SIGNAL(returnPressed()), this, SLOT(thresholdChanged()));
    connect(thresholdEdit, SIGNAL(textChanged(const QString &)), this, SLOT(thresholdChanging()));
    connect(thresholdWeightEdit, SIGNAL(returnPressed()), this, SLOT(thresholdWeightChanged()));
    connect(thresholdWeightEdit, SIGNAL(textChanged(const QString &)), this, SLOT(thresholdWeightChanging()));
    connect(thresholdLevelSelector, SIGNAL(currentIndexChanged(int)), this, SLOT(thresholdLevelChanged(int)));

}

bool QScriptFrame::createProject(int index) { // create a single project file - assumes path exists

    std::clog << Utilities::getTabString() <<index<<Utilities::getTabString()<<"Create"<<Utilities::getTabString()<<projectNameList.at(index).toStdString().c_str()<<endl;

    xmlFileAccessor->reset();
    xmlFileAccessor->setImageName(imageNameList.at(index));
    xmlFileAccessor->setMeshName(meshNameList.at(index));

    // set classifier threshold if specified
    if (modelTypeSelector->currentIndex() == itk::CorticalBone::kHighResClassifier) {
        int classifierThresholdIndex = thresholdLevelSelector->currentIndex();
        QString classifierThresholdName = thresholdLevelSelector->currentText();
        if(includeThresholds && includeThresholdWeight) {
            xmlFileAccessor->setClassifierThresholdInfo(thresholds[0], thresholds[1], thresholds[2],
                                                        classifierThresholdIndex, classifierThresholdName, thresholdWeight);
        } else if(includeThresholds && !includeThresholdWeight) {
            xmlFileAccessor->setClassifierThresholdInfo(thresholds[0], thresholds[1], thresholds[2],
                                                        classifierThresholdIndex, classifierThresholdName);
        } else if (!includeThresholds && includeThresholdWeight){
            xmlFileAccessor->setClassifierThresholdInfo(classifierThresholdIndex, classifierThresholdName, thresholdWeight);
        } else {
            xmlFileAccessor->setClassifierThresholdInfo(classifierThresholdIndex, classifierThresholdName);
        }
    }

    if (includeSampleNumber) {
        xmlFileAccessor->setSampleNumber(sampleNumberInt);
    }

    double projectFixedCBValue=nan("1");
    if (cbConstraintIndex == xmlFileAccess::kCBConstraintFixed && !fixedCBIsNumber) {
        QFile cbValueFile(fixedCBFileList.at(index));
        QTextStream cbValueStream(&cbValueFile);
        if (cbValueFile.open(QIODevice::ReadOnly)) {
            projectFixedCBValue = cbValueStream.readLine().toDouble();
            cbValueFile.close();
        } else {
            projectFixedCBValue = nan("1");
        }  //cout<<"parsed file: "<<fixedCBFileList.at(i).toStdString()<<". string="<<cbString.toStdString()<<", value="<<fixedCBValue<<endl;
    } else if (cbConstraintIndex == xmlFileAccess::kCBConstraintFixed && fixedCBIsNumber) {
        projectFixedCBValue=fixedCBValue;
    }
    xmlFileAccessor->setCBConstraintInfo(cbConstraintIndex, projectFixedCBValue);

    if(sigmaUpToDate&&sigmaFixedButton->isChecked()) {
        double sigma[3]; sigma[0] = sigmaXEdit->text().toDouble(); sigma[1] = sigmaYEdit->text().toDouble();
        sigma[2] = sigmaZEdit->text().toDouble();
        xmlFileAccessor->setSigmaConstraint(sigma);
    }

    if (includeAveraging) {
        xmlFileAccessor->setProfileAveraging(profileAverages);
    }

    if (importFileSet) {
        cerr<<"importFileList(i)="<<importedFileList.at(index).toStdString().c_str()<<", n="<<importedFileList.size()<<endl;
        xmlFileAccessor->setImportBaseName(importedFileList.at(index));
    }

    if (includeCalibration & !calibrationIsInFile) {
        xmlFileAccessor->setCalibrationParameters(calibrationScales[0], calibrationScales[1], calibrationScales[2]);
        xmlFileAccessor->setCalibrationName(calibrationTypeBox->currentIndex(), calibrationTypeBox->currentText());
    } else if(includeCalibration & calibrationIsInFile) {
        QFile pCalFile(calibrationFileList.at(index));
        QTextStream pCalStream(&pCalFile);
        if (pCalFile.open(QIODevice::ReadOnly)) {
            double calibrationScalesLocal[3];
            calibrationScalesLocal[0] = pCalStream.readLine().remove(0,3).toDouble();
            calibrationScalesLocal[1] = pCalStream.readLine().remove(0,3).toDouble();
            calibrationScalesLocal[2] = pCalStream.readLine().remove(0,3).toDouble();
            pCalFile.close();
            xmlFileAccessor->setCalibrationParameters(calibrationScalesLocal[0], calibrationScalesLocal[1], calibrationScalesLocal[2]);
            xmlFileAccessor->setCalibrationName(itk::CorticalBone::kManualQuadraticCal, calibrationTypeBox->itemText(itk::CorticalBone::kManualQuadraticCal));
        } else {
            projectFixedCBValue = nan("1");
        }  //cout<<"parsed file: "<<fixedCBFileList.at(i).toStdString()<<". string="<<cbString.toStdString()<<", value="<<fixedCBValue<<endl;
    }

    // set processor values
    // model info
    int modelIndex = modelTypeSelector->currentIndex();
    QString modelName = modelTypeSelector->currentText();
    int schemeIndex = schemeTypeSelector->currentIndex();
    QString schemeName = schemeTypeSelector->currentText();
    int optimiserIndex = optimiserTypeSelector->currentIndex();
    QString optimiserName = optimiserTypeSelector->currentText();
    xmlFileAccessor->setModelInfo(modelIndex, modelName, schemeIndex, schemeName, optimiserIndex,
                                  optimiserName);

    // smoothing value
    if (includeSmoothing) {
        xmlFileAccessor->setSmoothing(smoothingValue);
    }

    xmlFileAccessor->saveProjectFile(projectNameList.at(index), true);

    return true;

}

bool QScriptFrame::viewProject(int index) {
    cerr<<"View Script: File "<<viewIteratorIndex<<endl;
    std::clog<<Utilities::getTabString()<<viewIteratorIndex<<Utilities::getTabString()<<"View"<<Utilities::getTabString()<<projectNameList.at(viewIteratorIndex).toStdString().c_str()<<endl;
    if(!QFileInfo(projectNameList.at(viewIteratorIndex)).exists()) {
        QDir dir = QDir(); dir.mkdir(QFileInfo(projectNameList.at(index)).absolutePath());
        createProject(viewIteratorIndex);
    }
    parentObject->openProject(projectNameList.at(viewIteratorIndex));

    return true;
}

bool QScriptFrame::processProject(int index) {
    std::clog<<Utilities::getTabString()<<index<<Utilities::getTabString()<<"Process"<<Utilities::getTabString()<<projectNameList.at(index).toStdString().c_str()<<endl;
    cerr<<"Process Script: File "<<index<<endl;

    std::time_t currentTime = std::time(nullptr);
    cerr<<"    "<<std::asctime(std::localtime(&currentTime));

    if(!QFileInfo(projectNameList.at(index)).exists()) {
        QDir dir = QDir(); dir.mkdir(QFileInfo(projectNameList.at(index)).absolutePath());
        createProject(index);
    }

    parentObject->openProject(projectNameList.at(index));

    if(modelTypeSelector->currentIndex()==itk::CorticalBone::kHighResClassifier && !includeThresholds) { // don't try regardless as will just fail if bogus values
        parentObject->runThresholdCalculation();
    }
    // set disable all
    parentObject->runModellingOverMesh();
    // set disable all
    parentObject->saveProject(projectNameList.at(index)); // save so if any new info gained it is saved
    parentObject->closeProject();

    cerr<<endl;

    return true;
}

void QScriptFrame::updateObjectEnables() {

    this->setEnabled(true);

    updateEditorEnables(); // call first to update 'tabsUpToDate'


    if(runningViewScript) {

        scriptDirectoryButton->setDisabled(true);
        projectFilterLabel->setDisabled(true);
        projectFilterEdit->setDisabled(true);

        loadButton->setDisabled(true);
        saveButton->setDisabled(true);
        runButton->setDisabled(true);
        nextButton->setEnabled(true);

        disableTabs();
    } else if(runningProcessScript) {

        this->setDisabled(true);

    } else if(!projectDirectorySet) { // Nothing set
        // nothing selected - enable proj location & open
        scriptDirectoryButton->setEnabled(true);
        projectFilterLabel->setDisabled(true);
        projectFilterEdit->setDisabled(true);

        loadButton->setEnabled(true);
        saveButton->setDisabled(true);
        runButton->setDisabled(true);
        nextButton->setDisabled(true);

        disableTabs();

    } else if(projectDirectorySet && !tabsUpToDate) { // directory set, other options not set

        scriptDirectoryButton->setEnabled(true);
        projectFilterLabel->setEnabled(true);
        projectFilterEdit->setEnabled(true);

        loadButton->setEnabled(true);
        saveButton->setDisabled(true);
        runButton->setDisabled(true);
        nextButton->setDisabled(true);

        enableTabs();

    } else if (projectDirectorySet && tabsUpToDate) { // directory and options set, ready to run

        scriptDirectoryButton->setEnabled(true);
        projectFilterLabel->setEnabled(true);
        projectFilterEdit->setEnabled(true);

        loadButton->setEnabled(true);
        saveButton->setEnabled(true);
        runButton->setEnabled(true);
        nextButton->setDisabled(true);

        enableTabs();

    } else {
        std::clog<<"Invalid script state entered"<<endl; std::cerr<<"Invalid script state entered"<<endl;
    }
}

void QScriptFrame::enableTabs() {

    scriptControlsTab->setEnabled(true);

    if(includeSampleNumber) {
        sampleNumberEdit->setEnabled(true);
    } else {
        sampleNumberEdit->setDisabled(true);
    }

    if(cbDensityButton->isChecked()) {
        cbDensityFixedButton->setEnabled(true); cbDensityFWHMButton->setEnabled(true);
        if(cbDensityFixedButton->isChecked()) {
            cbDensityFixedEdit->setEnabled(true);
        } else {
            cbDensityFixedEdit->setDisabled(true);
        }
    } else {
        cbDensityFixedButton->setDisabled(true); cbDensityFWHMButton->setDisabled(true);
    }

    if(includeAveraging) {
        xAvgEdit->setEnabled(true); yAvgEdit->setEnabled(true); zAvgEdit->setEnabled(true);
    } else {
        xAvgEdit->setDisabled(true); yAvgEdit->setDisabled(true); zAvgEdit->setDisabled(true);
    }

    if(importFileSet) {
    } else {
    }

    if(includeCalibration && !calibrationIsInFile) {

        p0ScaleEdit->setEnabled(true); p1ScaleEdit->setEnabled(true);
        if(calibrationTypeBox->currentIndex()==itk::CorticalBone::kManualLinearCal) {
            p2ScaleEdit->setDisabled(true);
        } else {
            p2ScaleEdit->setEnabled(true);
        }
        p0ScaleLabel->setVisible(true); p1ScaleLabel->setVisible(true); p2ScaleLabel->setVisible(true);
        p0ScaleEdit->setVisible(true); p1ScaleEdit->setVisible(true); p2ScaleEdit->setVisible(true);


        pFileEdit->setVisible(false); pFileLabel->setVisible(false);
        pFileEdit->setDisabled(true); pFileLabel->setDisabled(true);

    } else if(includeCalibration && calibrationIsInFile) {

        pFileEdit->setVisible(true); pFileLabel->setVisible(true);
        pFileEdit->setEnabled(true); pFileLabel->setEnabled(true);

        p0ScaleLabel->setVisible(false); p1ScaleLabel->setVisible(false); p2ScaleLabel->setVisible(false);
        p0ScaleEdit->setVisible(false); p1ScaleEdit->setVisible(false); p2ScaleEdit->setVisible(false);
        p0ScaleEdit->setDisabled(true); p1ScaleEdit->setDisabled(true); p2ScaleEdit->setDisabled(true);

    } else {

        pFileEdit->setVisible(false); pFileLabel->setVisible(false);
        pFileEdit->setDisabled(true); pFileLabel->setDisabled(true);

        p0ScaleLabel->setVisible(false); p1ScaleLabel->setVisible(false); p2ScaleLabel->setVisible(false);
        p0ScaleEdit->setVisible(false); p1ScaleEdit->setVisible(false); p2ScaleEdit->setVisible(false);
        p0ScaleEdit->setDisabled(true); p1ScaleEdit->setDisabled(true); p2ScaleEdit->setDisabled(true);
    }

    if(includeThresholdWeight) {
        thresholdWeightBox->setVisible(true); thresholdWeightBox->setEnabled(true);
    } else {
        thresholdWeightBox->setVisible(false); thresholdWeightBox->setDisabled(true);
    }

    if(includeThresholds) {
        thresholdBox->setVisible(true); thresholdBox->setEnabled(true);
    } else {
        thresholdBox->setVisible(false); thresholdBox->setDisabled(true);
    }

    if(includeSmoothing) {
        smoothingEdit->setEnabled(true);
    } else {
        smoothingEdit->setDisabled(true);
    }

}

void QScriptFrame::disableTabs() {

    scriptControlsTab->setDisabled(true);

}

void QScriptFrame::updateEditorEnables() {


    if (!projectDirectorySet || !projectFilterSet) {
        scriptEditor->setDisabled(true); tabsUpToDate = false;
    } else if(isCreateTabUpToDate() && isViewTabUpToDate() && isProcessTabUpToDate()) {
        scriptEditor->setEnabled(true); tabsUpToDate = true;
    } else {
        scriptEditor->setDisabled(true); tabsUpToDate = false;
    }

}

bool QScriptFrame::isCreateTabUpToDate() {

    if(!imageFilterSet || !meshFilterSet) {
        return false;
    }
    if(includeSampleNumber && !sampleNumberSet) {
        return false;
    }
    if(!cbConstraintUpToDate || !sigmaUpToDate) {
        return false;
    }
    if(includeAveraging) {
        if(!xAvgSet || !yAvgSet || !zAvgSet) {
            return false;
        }
    }
    return true;
}

bool QScriptFrame::isViewTabUpToDate() {
    if(includeCalibration && p0ScaleSet && p1ScaleSet && p2ScaleSet && !calibrationIsInFile) {
        return true;
    } else if(includeCalibration && calibrationIsInFile && pFileSet) {
        return true;
    } else if(!includeCalibration) {
        return true;
    } else {
        return false;
    }
}

bool QScriptFrame::isProcessTabUpToDate() {

    if(includeSmoothing && !smoothingSet) {
        return false;
    } else if(!includeThresholds && !includeThresholdWeight) {
        return true;
    } else if(includeThresholds && stThresholdSet && cbThresholdSet && thresholdSet && !includeThresholdWeight) {
        return true;
    } else if(includeThresholds && stThresholdSet && cbThresholdSet && thresholdSet && includeThresholdWeight && thresholdWeightSet) {
        return true; // Should never be a valid combination. Either manual thresholds or else set threshold weight.
    } else if(!includeThresholds && includeThresholdWeight && thresholdWeightSet) {
        return true;
    } else {
        return false;
    }
}

void QScriptFrame::updateScript(bool resetRanges) {

    updateFileLists();
    if(resetRanges) {
        resetRangeLimits();
    }

    updateDisplay();

}

void QScriptFrame::updateDisplay() {

    scriptEditor->clear();

    scriptEditor->appendPlainText(tr("//----- Script to Create New Project Files -----//"));
    scriptEditor->appendPlainText(tr("Script Directory=")+projectPath);

    if(includeSampleNumber) {
        scriptEditor->appendPlainText(tr("Fixed Sample Number=")+QString::number(sampleNumberInt));
    }
    if(includeAveraging) {
        scriptEditor->appendPlainText(tr("Profile Averaging. x=")+QString::number(profileAverages[0])
                                      +tr(", y=")+QString::number(profileAverages[1])+tr(", x=")+QString::number(profileAverages[2]));
    }
    if(includeCalibration && !calibrationIsInFile) {
        scriptEditor->appendPlainText(tr("Density Calibration. P0=")+QString::number(calibrationScales[0])
                                      +tr(", P1=")+QString::number(calibrationScales[1])+tr(", P2=")+QString::number(
                calibrationScales[2]));
    } else if(includeCalibration && calibrationIsInFile) {
        scriptEditor->appendPlainText(tr("Density Calibration. Parameter File=")+calibrationTypeBox->currentText());
    }
    if(includeThresholds) {
        scriptEditor->appendPlainText(tr("Classifier Thresholds (Manually set). P<sub>st</sub>=")+QString::number(thresholds[0])
                                      +tr(", P<sub>cb</sub>=")+QString::number(thresholds[1])+tr(", P<sub>thresh</sub>=")+QString::number(
                thresholds[2]));
    }
    if(includeThresholdWeight) {
        scriptEditor->appendPlainText(tr("Classifier Thresholds (Median Manual Mode). Weight=")+QString::number(thresholdWeight));
    }
    if(includeSmoothing) {
        scriptEditor->appendPlainText(tr("Smoothing Gaussian Sigma=")+QString::number(smoothingValue));
    }
    scriptEditor->appendPlainText(tr("Model Info: Type=")+modelTypeSelector->currentText()
                                  +tr(",Optimiser=")+optimiserTypeSelector->currentText()+tr(", Fitting Scheme=")+schemeTypeSelector->currentText()); // todo - only include optimiser and scheme info if relevent
    if(cbConstraintIndex==xmlFileAccess::kCBConstraintFixed && fixedCBIsNumber) {
        scriptEditor->appendPlainText(tr("Fixed CB Value=") + fixedCBFilterString);
    }


    QString indent = QString::fromStdString(Utilities::getTabString());
    int startIndex = startSpinBox->value(), stopIndex = stopSpinBox->value();
    for(int i=0; i<projectNameList.size(); i++) {

        QString style, spaces;
        if(i<startIndex || i>stopIndex) { // gray out if filtered by display ranges
            style = tr("<font color=\"Grey\"><i>");
        } else {
            style = tr("<font color=\"Black\">");
        }

        // number
        scriptEditor->appendHtml(style+QString::number(i)+tr(" - "));
        QFileInfo projectFileInfo = QFileInfo(projectNameList.at(i));

        // path
        spaces=QString::fromStdString(Utilities::getSpacesString(21));
        scriptEditor->appendHtml(style+indent +tr("Path -")+spaces+projectFileInfo.absolutePath());
        // project file
        spaces=QString::fromStdString(Utilities::getSpacesString(4));
        scriptEditor->appendHtml(style+indent +tr("Proj File -")+spaces+projectFileInfo.fileName());
        // image files
        spaces=QString::fromStdString(Utilities::getSpacesString(10));
        if(imageExtnSelector->currentIndex()==kRawExtn) {
            scriptEditor->appendHtml(style+indent +tr("Image File -")+spaces+QFileInfo(imageNameList.at(i)).fileName());
        } else if(imageExtnSelector->currentIndex()==kQCTExtn) {
            scriptEditor->appendHtml(style+indent +tr("Image File -")+spaces+QFileInfo(imageNameList.at(i)).fileName());
        } else if (imageExtnSelector->currentIndex()==kDicomExtn) {
            scriptEditor->appendHtml(style+indent +tr("Image File -")+spaces+tr("/")+imageFilterString+tr("/*")+imageExtnString);
        } else if (imageExtnSelector->currentIndex()==kTiffExtn) {
            scriptEditor->appendHtml(style+indent +tr("Image File -")+spaces+tr("/")+imageFilterString+tr("/*")+imageExtnString);
        }
        // mesh files
        spaces=QString::fromStdString(Utilities::getSpacesString(12));
        scriptEditor->appendHtml(style+indent +tr("Mesh File -")+spaces+QFileInfo(meshNameList.at(i)).fileName());
        // cb constraint
        spaces=QString::fromStdString(Utilities::getSpacesString(2));
        if(cbConstraintIndex==xmlFileAccess::kCBConstraintFixed && !fixedCBIsNumber) {
            scriptEditor->appendHtml(style+indent+tr("Fixed CB File -") + spaces + QFileInfo(fixedCBFileList.at(i)).fileName());
        }
        // calibration parameter file
        else if(includeCalibration && calibrationIsInFile) {
            scriptEditor->appendHtml(style+indent+tr("Calibration File -") + spaces + QFileInfo(calibrationFileList.at(i)).fileName());
        }
        // import path
        spaces=QString::fromStdString(Utilities::getSpacesString(7));
        if(importFileSet) {
            scriptEditor->appendHtml(style+indent +tr("Import Path -")+ QString::fromStdString(Utilities::getSpacesString(4)) +importedFileList.at(i));
        }
        // results
        spaces=QString::fromStdString(Utilities::getSpacesString(8));
        scriptEditor->appendHtml(style+indent +tr("Results File -")+spaces+QFileInfo(resultsNameList.at(i)).baseName()+"_Parameters|Displays|ImageProfiles"+txtExtnString);
    }

}

void QScriptFrame::updateFileLists() {

    QStringList directoryList = generateDirectoryList();

    // filter the lists to include only those folders with only one image, mesh and settings file
    projectNameList = imageNameList = meshNameList = fixedCBFileList = importedFileList = resultsNameList = calibrationFileList = QStringList();
    for(int i=0; i<directoryList.size(); i++) {

        QString tempImageString, tempFixedCBString, tempCalibrationString; QStringList tempMeshList;

        QDir directory = QDir(directoryList.at(i));


        // check first for meshes in the directories
        directory.setNameFilters(QStringList()<<"*"+meshFilterString+"*"+meshExtnString);
        tempMeshList = directory.entryList();
        if(tempMeshList.size()<1) {
            continue;
        }

        // get images, meshes and possibly settings
        tempImageString = getImageFile(directoryList.at(i));

        if(tempImageString.isEmpty()) {
            continue;
        }

        if(cbConstraintIndex==xmlFileAccess::kCBConstraintFixed && !fixedCBIsNumber) {
            directory.setNameFilters(QStringList()<<"*"+fixedCBFilterString+"*"+txtExtnString);
            QStringList tempFixedCBFileList = directory.entryList();

            if(directory.entryList().size() == 1) {
                tempFixedCBString = tempFixedCBFileList.first();
            } else {
                continue;
            }
        }

        if(includeCalibration && calibrationIsInFile) {
            directory.setNameFilters(QStringList()<<"*"+pFileEdit->text()+"*"+txtExtnString);
            QStringList tempCalibrationFileList = directory.entryList();

            if(directory.entryList().size() == 1) {
                tempCalibrationString = tempCalibrationFileList.first();
            } else {
                continue;
            }
        }

        for(int meshIndex=0; meshIndex<tempMeshList.size(); meshIndex++) {

            QString path = directoryList.at(i)+QDir::separator();

            imageNameList.append(path+tempImageString);
            meshNameList.append(path+tempMeshList.at(meshIndex));
            if(cbConstraintIndex==xmlFileAccess::kCBConstraintFixed && !fixedCBIsNumber) {
                fixedCBFileList.append(path+tempFixedCBString);
            }

            if(includeCalibration && calibrationIsInFile) {
                calibrationFileList.append(path + tempCalibrationString);
            }

            QString combinedName;
            QString imageName = QFileInfo(tempImageString).baseName();
            QString meshName = QFileInfo(tempMeshList.at(meshIndex)).baseName();
            if(imageName.contains(meshName)) {
                combinedName = imageName;
            } else if(meshName.contains(imageName)) {
                combinedName = meshName;
            } else {
                combinedName = imageName+tr("_")+meshName;
            }
            QString projectString=path+combinedName+projectFilterString+QDir::separator()+combinedName+projectFilterString+projectExtnString; // implicitly included in project folder if creating or opening a newfile
            projectNameList.append(projectString);
        }
    }

    // create a import mapping list if set
    if(importFileSet) { // lr path, hr path
        QFile importTextFile(importFileString);
        importedFileList = QStringList(projectNameList); //QStringList tempImportList = QStringList(projectNameList);
        int numberOfImportMatches=0;
        if (importTextFile.open(QIODevice::ReadOnly)) {
            QTextStream textStream(&importTextFile);
            while (true) {
                QString line = textStream.readLine();
                if (line.isNull()) {
                    break;
                } else {

                    QStringList importPair = line.split( " " );

                    if(importPair.size()!=2) {
                        cerr<<"incorrect import pair in the import mapping file: "<<line.toStdString().c_str()<<endl;
                        break; // error
                    }

                    int index = projectNameList.indexOf(importPair.at(0));
                    if(index!=-1) { // valid match, otherwise ignore
                        importedFileList.replace(index, importPair.at(1));
                        numberOfImportMatches++;
                    }
                }

                if(numberOfImportMatches==projectNameList.size()) {
                    break; //all matches made
                }
            }
        }
        if(numberOfImportMatches != projectNameList.size()) { //importedFileList.size() != projectNameList.size()) {
            cerr<<"Miss-match between import mapping file (n="<<numberOfImportMatches<<") and the create script (n="<<projectNameList.size()<<")"<<endl;
            projectNameList = imageNameList = meshNameList = fixedCBFileList = importedFileList = resultsNameList = calibrationFileList = QStringList();
        }
    }


    // generate results list
    resultsNameList = QStringList();
    for(int i=0; i<projectNameList.size(); i++) {

        QFileInfo projectFileInfo = QFileInfo(projectNameList.at(i));
        QString resultString=projectFileInfo.absolutePath()+QDir::separator()+projectFileInfo.baseName();

        resultsNameList<<resultString;
    }

}

QString QScriptFrame::getImageFile(QString dirName) {

    QDir directory = QDir(dirName);

    // get images, meshes and possibly settings
    directory.setNameFilters(QStringList()<<"*"+imageFilterString+"*"+imageExtnString);
    if(directory.entryList().size()==1 && imageExtnSelector->currentIndex()==kRawExtn) {
        return directory.entryList().first();
    } else if(directory.entryList().size()==1 && imageExtnSelector->currentIndex()==kQCTExtn) {
        return directory.entryList().first();
    } else if(imageExtnSelector->currentIndex()==kDicomExtn && directory.entryList().size()>=1) {
        return directory.entryList().first(); // no longer need for just directory //QFileInfo(directory.entryList().first()).absolutePath()+QDir::separator();
    } else if(imageExtnSelector->currentIndex()==kDicomExtn) { // make one final check for files with no extension or folders of the appropiate name

        // check for a folder with files with no extension (no '.')
        directory.setNameFilters(QStringList()<<"*"+imageFilterString+"*");
        QStringList tempImageNoExtnList = directory.entryList();

        for(int j=0; j<tempImageNoExtnList.size(); j++) {

            // try to locate .dcm file that are missing the .dcm extension
            if(!tempImageNoExtnList.at(j).contains(tr(".")) && !tempImageNoExtnList.at(j).contains(tr(".tif")) && QFileInfo(tempImageNoExtnList.at(j)).isFile()) { // both no '.' in name and a file not a directory

                return directory.entryList().at(j);
            }
        }

        // one final check - check if contains a sub dir

            // otherwise check for subdirectories with the filter - match in its intirety
        QDir subdirectoryFilter(dirName); subdirectoryFilter.setNameFilters(QStringList()<<imageFilterString);
        subdirectoryFilter.setFilter(QDir::Dirs | QDir::NoSymLinks | QDir::NoDotAndDotDot);

        if (subdirectoryFilter.entryList().size() == 1) {

            return subdirectoryFilter.entryList().first()+QDir::separator();
        }

    } else if(imageExtnSelector->currentIndex()==kTiffExtn &&directory.entryList().size()>=1) {  // tiffs in same folder as mesh
        return directory.entryList().first();
    } else if(imageExtnSelector->currentIndex()==kTiffExtn) { // check for a subdirectory with the filter name - that contains tiffs


        QDirIterator tiffDirectorySearch(dirName, QDir::Dirs | QDir::NoSymLinks | QDir::NoDotAndDotDot,
                                         QDirIterator::Subdirectories);

        int numberTiffDirs = 0;
        QString tempImageString;
        while (tiffDirectorySearch.hasNext()) {
            tiffDirectorySearch.next();

            if (!tiffDirectorySearch.filePath().contains(imageFilterString)) {
                continue; // ignore folders not matching the specified image filter
            }

            QDir imageFolderDirectory = QDir(tiffDirectorySearch.filePath());
            imageFolderDirectory.setNameFilters(QStringList() << "*.tif"); // todo support other styles of tiff files

            if (imageFolderDirectory.entryList().size() >= 1) {
                QStringList tiffList = imageFolderDirectory.entryList();
                for (int k = 0; k < tiffList.size(); k++) {
                    QRegularExpression pattern("(\\d+)(.tif)$"); // must have form xxx#.tif
                    if (tiffList.at(k).contains(pattern)) {
                        tempImageString = tiffDirectorySearch.fileName()+QDir::separator()+tiffList.at(k);
                        numberTiffDirs++;
                        break; // break having found file of the correct form
                    }
                }
            }
        }

        if(numberTiffDirs>1) {
            return QString();// too many folders containing tif files. break out of while loop
        }
        return tempImageString;

    }

    return QString();

}

void QScriptFrame::resetRangeLimits() { // call to reset the range limits

    // only adjust if in process tab, otherwise the possible number of files to process is unknown
    int maxValue = std::max(0, projectNameList.size() - 1); // not results but project file lists

    // the min of the little one is always 0
    if(stopSpinBox->maximum() != maxValue) { // if max val have changed

        startSpinBox->blockSignals(true);
        stopSpinBox->blockSignals(true); // turn off signals

        int startValue = 0, endValue = maxValue;
        startSpinBox->setRange(0, endValue);
        stopSpinBox->setRange(startValue, maxValue);
        startSpinBox->setValue(startValue);
        stopSpinBox->setValue(endValue);

        startSpinBox->blockSignals(false);
        stopSpinBox->blockSignals(false); // turn on signals

    }
}

void QScriptFrame::updateSchemeOptions(int index) {

    int numberOfOptions = schemeTypeSelector->count();
    for(int i=numberOfOptions-1; i>=0; i--) {
        schemeTypeSelector->removeItem(i);
    }

    if(index<=itk::CorticalBone::kEndostealRamp) { // model fitting options

        schemeTypeSelector->addItem("std A"); // rect 1st then ramp if ramp selected - huerstic weighting
        schemeTypeSelector->addItem("std B"); // rect 1st then ramp if ramp selected - narrow weighting
        schemeTypeSelector->addItem("std C"); // rect 1st then ramp if ramp selected - wide weighting
        schemeTypeSelector->addItem("std D"); // rect 1st then ramp if ramp selected - narrow to wide weighting
        schemeTypeSelector->addItem("CBMV2a"); // (selected model fit 2x + sigma correction) - static huristic weighting
        schemeTypeSelector->addItem("CBMV2b"); // (selected model fit 2x + sigma correction) - static narrow weighting
        schemeTypeSelector->addItem("CBMV2c"); // (selected model fit 2x + sigma correction) - static narrow then wide weighting
        schemeTypeSelector->addItem("CBMV2d"); // (selected model fit 2x + sigma correction) - dynamic narrow then wide weighting
        schemeTypeSelector->addItem("Unconstrained a"); // rect followed by ramp/rect followed by ramp/rect no constraints
        schemeTypeSelector->addItem("Unconstrained b");
        schemeTypeSelector->addItem("Unconstrained c");
        schemeTypeSelector->addItem("Smoothing a"); // V4 then smooth the CB estimates then V1
        schemeTypeSelector->addItem("Smoothing b");
        schemeTypeSelector->addItem("Smoothing c");
        schemeTypeSelector->addItem("Smoothing d");
        schemeTypeSelector->setCurrentIndex(0);

    } else if(index==itk::CorticalBone::kHighResClassifier) { // classifier options
        schemeTypeSelector->addItem("Global Thresholds"); // todo pull out schemeTypeSelection update options into a separate method
        schemeTypeSelector->addItem("Median Thresholds");
        schemeTypeSelector->setCurrentIndex(0);

    } else if(index==itk::CorticalBone::kCalibration) { // classifier options
        schemeTypeSelector->addItem("Median Upper"); // todo pull out schemeTypeSelection update options into a separate method
        schemeTypeSelector->addItem("Maxium Upper");
        schemeTypeSelector->setCurrentIndex(0);

    }
//    else if(index==itk::CorticalBone::kHighResThresholder) { // thresholder options
//        schemeTypeSelector->addItem("Global, Edge Detection");
//        schemeTypeSelector->addItem("Median, Edge Detection");
//        schemeTypeSelector->addItem("Linear Regression");
//        schemeTypeSelector->addItem("Average Regression");
//        schemeTypeSelector->addItem("Deming Regression");
//        schemeTypeSelector->setCurrentIndex(0);
//
//    } else if(index==itk::CorticalBone::kHighResOptimiser) { // optimiser options
//        schemeTypeSelector->addItem("Endosteal Edges");
//        schemeTypeSelector->setCurrentIndex(0);
//    }
}

void QScriptFrame::processMinRangeChanged(int value) {

    // the end spinner minimum limit is changed
    stopSpinBox->setMinimum(value);
    updateDisplay();
}

void QScriptFrame::processMaxRangeChanged(int value) {
    // the end spinner minimum limit is changed
    startSpinBox->setMaximum(value);
    updateDisplay();
}

QStringList QScriptFrame::generateFileList(QString filter) {

    QDirIterator directoryIterator(projectPath, QStringList()<<filter, QDir::Files, QDirIterator::Subdirectories);

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

QStringList QScriptFrame::generateDirectoryList() { // no filter input as filters name and only looking at directories

    QDirIterator directoryIterator(projectPath, QStringList(), QDir::Dirs, QDirIterator::Subdirectories);

    QStringList dirList = QStringList(projectPath);

    int i=0;
    while(directoryIterator.hasNext() && i<MAX_NUMBER_OF_FILES_TO_PROCESS) {
        QString directory = directoryIterator.next();
        if(!directory.endsWith(tr(".")) && !directory.endsWith(tr(".."))) {
            dirList.append(directory);
            i++;
        }
    }

    if(i>=MAX_NUMBER_OF_FILES_TO_PROCESS) {
        std::clog<<"Max number of directories reached for path"<<projectPath.toStdString()<<" others ignored."<<std::endl;
    }

    return dirList;
}


//------------- Private Slots ------------//

// project name filter
void QScriptFrame::projectFilterChanged() {

    projectFilterSet = true;
    projectFilterString = projectFilterEdit->text();

    updateScript();
    updateObjectEnables();

}

void QScriptFrame::projectFilterChanging() {

    projectFilterSet = false;
    updateObjectEnables();

}

void QScriptFrame::sampleNumberChanging() {
    sampleNumberSet = false;
    updateObjectEnables();
}

void QScriptFrame::sampleNumberChanged() {

    sampleNumberInt = sampleNumberEdit->text().toInt();
    sampleNumberSet = true;

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::toggleSampleNumber(bool state) {
    includeSampleNumber = state;
    sampleNumberSet = false;

    updateScript();
    updateObjectEnables();
}

// cb constraints
void QScriptFrame::cbValueFilterChanged() {
    fixedCBFilterString = cbDensityFixedEdit->text();
    cbConstraintUpToDate = fixedCBSet = true;

    bool status; double value=fixedCBFilterString.toDouble(&status);
    if (status) {
        fixedCBIsNumber=true; fixedCBValue=value;
    } else {
        fixedCBIsNumber=false; fixedCBValue=nan("1");
    }

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::cbValueFilterChanging() {
    cbConstraintUpToDate = fixedCBSet = false;
    updateObjectEnables();
}

void QScriptFrame::toggleFixedCB(bool state) {
    if(state) {
        cbConstraintUpToDate = fixedCBSet;
        cbConstraintIndex=xmlFileAccess::kCBConstraintFixed;
        updateScript();
        updateObjectEnables();
    }
}

void QScriptFrame::toggleCBConstraintState(bool state) {

    if(state) { // update 'uptodate' bool
        if(cbDensityFWHMButton->isChecked()) {
            cbConstraintUpToDate = true;
            cbConstraintIndex=xmlFileAccess::kCBConstraintFWHM;
        } else {
            cbConstraintUpToDate = fixedCBSet;
            cbConstraintIndex=xmlFileAccess::kCBConstraintFixed;
        }
    } else { // disable FWHM and Fixed
        cbConstraintUpToDate = true;
        cbConstraintIndex=xmlFileAccess::kCBConstraintNone;
    }

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::toggleFWHMCB(bool state) {
    if(state) {
        cbConstraintIndex=xmlFileAccess::kCBConstraintFWHM;

        updateScript();
        updateObjectEnables();
    }
}

// sigma
void QScriptFrame::toggleFixedSigma(bool state) {
    if(state) {
        if(sigmaXSet && sigmaYSet && sigmaZSet) {
            sigmaUpToDate = true;
        } else {
            sigmaUpToDate = false;
        }
    } else if(!state && sigmaRequired) {
        sigmaUpToDate = false;
    } else {
        sigmaUpToDate = true;
    }
    updateObjectEnables();
}

void QScriptFrame::fixedSigmaXChanged() {
    sigmaXSet = true;
    if(sigmaXSet && sigmaYSet && sigmaZSet) {
        sigmaUpToDate = true;
    }
}

void QScriptFrame::fixedSigmaXChanging() {
    sigmaXSet = sigmaUpToDate = false;
}

void QScriptFrame::fixedSigmaYChanged() {
    sigmaYSet = true;
    if(sigmaXSet && sigmaYSet && sigmaZSet) {
        sigmaUpToDate = true;
    }
}

void QScriptFrame::fixedSigmaYChanging() {
    sigmaYSet = sigmaUpToDate = false;
}

void QScriptFrame::fixedSigmaZChanged() {
    sigmaZSet = true;
    if(sigmaXSet && sigmaYSet && sigmaZSet) {
        sigmaUpToDate = true;
    }
}

void QScriptFrame::fixedSigmaZChanging() {
    sigmaZSet = sigmaUpToDate = false;
}

// averaging
void QScriptFrame::xAveragingChanging() {
    xAvgSet = false;
    updateObjectEnables();
}

void QScriptFrame::yAveragingChanging() {
    yAvgSet = false;
    updateObjectEnables();
}

void QScriptFrame::zAveragingChanging() {
    zAvgSet = false;
    updateObjectEnables();
}

void QScriptFrame::setXAveragingValue() {
    xAvgSet = true;
    profileAverages[0] = xAvgEdit->text().toDouble();

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::setYAveragingValue() {
    yAvgSet = true;
    profileAverages[1] = yAvgEdit->text().toDouble();

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::setZAveragingValue() {
    zAvgSet = true;
    profileAverages[2] = zAvgEdit->text().toDouble();

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::toggleAveraging(bool state) {
    includeAveraging = state;
    profileAverages[0] = profileAverages[1] = profileAverages[2] = nan("1");

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::setImportFile() {

    // open file dialogue
    QString fileName = QFileDialog::getOpenFileName(this, tr("Select File containing Import Mapping"), projectPath, tr("Mapping Files (*.txt)"));
    if (!fileName.isEmpty()) {

        importFileString = fileName; importFileSet = true;

    } else {
        std::clog << "Message - No import file selected." << std::endl;
        importFileSet = false;
        importFileString = QString();
    }
    importFileEdit->setText(importFileString);

    updateScript();
    updateObjectEnables();

}

// calibration
void QScriptFrame::p0ScaleChanging() {
    p0ScaleSet = false;
    updateObjectEnables();
}

void QScriptFrame::p1ScaleChanging() {
    p1ScaleSet = false;
    updateObjectEnables();
}

void QScriptFrame::p2ScaleChanging() {
    p2ScaleSet = false;
    updateObjectEnables();
}

void QScriptFrame::pFileTextChanging() {
    pFileSet = false;
    updateObjectEnables();
}

void QScriptFrame::setp0ScaleValue() {
    p0ScaleSet = true;
    calibrationScales[0] = p0ScaleEdit->text().toDouble();

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::setp1ScaleValue() {
    p1ScaleSet = true;
    calibrationScales[1] = p1ScaleEdit->text().toDouble();

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::setp2ScaleValue() {
    p2ScaleSet = true;
    calibrationScales[2] = p2ScaleEdit->text().toDouble();

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::setpFileTextChanged() {
    pFileSet = true;

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::calibrationTypeChanged(int index) {

    if(index==itk::CorticalBone::kManualLinearCal) {
        includeCalibration = true; calibrationIsInFile = false;
        p2ScaleEdit->setText(QString::number(0.0)); p2ScaleSet = true; calibrationScales[2] = 0.0;
    } else if(index==itk::CorticalBone::kManualQuadraticCal) {
        includeCalibration = true; calibrationIsInFile = false;
    } else if(index==itk::CorticalBone::kFileSpecifiedCal) {
        includeCalibration = calibrationIsInFile = true;
    } else {
        includeCalibration = calibrationIsInFile = false;
    }

    //densityScales[0] = densityScales[1] = densityScales[2] = nan("1"); - leave unchanged. only store if includeDensityScaling is true

    updateScript();
    updateObjectEnables();
}

// threholds
void QScriptFrame::thresholdLevelChanged(int index) {

    if(index==itk::LinearTransform::kManualThreshold) {
        includeThresholds = true; includeThresholdWeight = false;
    } else if(index==itk::LinearTransform::kMedianManual) {
        includeThresholdWeight = true; includeThresholds = false;
    }else {
        includeThresholds = includeThresholdWeight = false;
    }

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::stThresholdChanged() {
    stThresholdSet = true;
    thresholds[0] = stThreshEdit->text().toDouble();

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::cbThresholdChanged() {
    cbThresholdSet = true;
    thresholds[1] = cbThreshEdit->text().toDouble();

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::thresholdChanged() {
    thresholdSet = true;
    thresholds[2] = thresholdEdit->text().toDouble();

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::stThresholdChanging() {

    stThresholdSet = false;
    updateObjectEnables();
}

void QScriptFrame::cbThresholdChanging() {
    cbThresholdSet = false;
    updateObjectEnables();
}

void QScriptFrame::thresholdChanging() {
    thresholdSet = false;
    updateObjectEnables();
}

void QScriptFrame::thresholdWeightChanged() {
    thresholdWeightSet = true;
    thresholdWeight = thresholdWeightEdit->text().toDouble();

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::thresholdWeightChanging() {
    thresholdWeightSet = false;
    updateObjectEnables();
}

// image name filter
void QScriptFrame::imageFilterChanged() {

    imageFilterSet = true;
    imageFilterString = imageFilterEdit->text();

    updateScript();
    updateObjectEnables();

}

void QScriptFrame::imageFilterChanging() {

    imageFilterSet = false;
    updateObjectEnables();
}

// mesh name filter
void QScriptFrame::meshFilterChanged() {

    meshFilterSet = true;
    meshFilterString = meshFilterEdit->text();

    updateScript();
    updateObjectEnables();

}

void QScriptFrame::meshFilterChanging() {

    meshFilterSet = false;
    updateObjectEnables();
}

// process model selector - and other processing options
void QScriptFrame::modelSelectionChanged(int index) {
    if(index<=itk::CorticalBone::kEndostealRamp) { // hide classifier
        includeThresholds = includeThresholdWeight = false;
        optimiserTypeLabel->show(); optimiserTypeSelector->show();
        classifierLevelLabel->hide(); thresholdLevelSelector->hide();
    } else if(index==itk::CorticalBone::kHighResClassifier) { // hide fitting
        optimiserTypeLabel->hide(); optimiserTypeSelector->hide();
        classifierLevelLabel->show(); thresholdLevelSelector->show();
    } else if(index<=itk::CorticalBone::kCalibration) {
        includeThresholds = includeThresholdWeight = false;
        optimiserTypeLabel->hide(); optimiserTypeSelector->hide();
        classifierLevelLabel->hide(); thresholdLevelSelector->hide();
    } else {
        cerr<<"Error invalid model index selected in QScriptFrame::modelSelectionChanged"<<endl;
    }
    updateSchemeOptions(index);

}

void QScriptFrame::schemeChanged(int index) {
    if(modelTypeSelector->currentIndex()<=itk::CorticalBone::kEndostealRamp && (schemeTypeSelector->currentIndex()>=itk::CorticalBone::kCBMV2AFitting&&schemeTypeSelector->currentIndex()<=itk::CorticalBone::kCBMV2DFitting)) {
        sigmaRequired=true;

    } else {
        sigmaRequired=false;
    }

    toggleFixedSigma(sigmaFixedButton->isChecked());
}

void QScriptFrame::smoothingChanging() {
    smoothingSet = false;
    updateObjectEnables();
}

void QScriptFrame::smoothingChanged() {

    smoothingValue = smoothingEdit->text().toDouble();
    smoothingSet = true;

    updateScript();
    updateObjectEnables();
}

void QScriptFrame::toggleSmoothing(bool state) {
    includeSmoothing = state;
    smoothingSet = false;

    updateScript();
    updateObjectEnables();
}

//--- dropdown menus ---//
void QScriptFrame::imageExtnChanged(int index) {

    if(index==kDicomExtn) {
        imageExtnString=tr(".dcm"); // todo other cases
    } else if (index==kRawExtn) {
        imageExtnString=tr(".mhd");
    } else if (index==kQCTExtn) {
        imageExtnString=tr(".QCT");
    } else if (index==kTiffExtn) {
        imageExtnString=tr(".tif");
    } else {
        imageExtnString=tr("");
        std::cerr<<"Invalid index in QScriptFrame::imageExtnChanged()"<<std::endl;
    }
    updateScript();
    updateObjectEnables();

}

void QScriptFrame::meshExtnChanged(int index) {

    if(index==kObjExtn) {
        meshExtnString=tr(".obj"); // todo other cases
    } else if (index==kVrmlExtn) {
        meshExtnString=tr(".wrl");
    } else if (index==kStlExtn) {
        meshExtnString=tr(".stl");
    } else if (index==kPlyExtn) {
        meshExtnString=tr(".ply");
    } else {
        meshExtnString=tr("");
        std::cerr<<"Invalid index in QScriptFrame::imageExtnChanged()"<<std::endl;
    }
    updateScript();
    updateObjectEnables();

}

//--- buttons ---//
void QScriptFrame::setDirectory() {

    // open file dialogue
    QString filePath = QFileDialog::getExistingDirectory(this, tr("Select Directory in which to Run Script"), projectPath);
    if (!filePath.isEmpty()) {

        projectPath = filePath; projectDirectorySet = true;
        scriptDirectoryEdit->setText(projectPath);

    } else {
        std::clog << "Error - selected file is empty" << std::endl;
        projectDirectorySet = false;
        scriptDirectoryEdit->setText(tr(""));
    }
    updateScript();
    updateObjectEnables();

}

void QScriptFrame::saveScript() {

    QString fileName, scriptPath = parentObject->getScriptFileName();

    // get file [.xml]
    fileName = QFileDialog::getSaveFileName(this, tr("Save Script File"), scriptPath, tr("Script Files (*.xmls)"));

    if (fileName.isEmpty()) {
        std::clog<<"In QScriptFrame::saveScript() no file name selected so nothing opened."<<std::endl;
        return;
    }

    if(!fileName.endsWith(".xmls")) {
        fileName.append(".xmls");
    }

    std::clog<<"<b>---- Save Script ----</b>"<<endl;

    // reset script
    xmlFileAccessor->reset();

    // general script info
    xmlFileAccessor->setScriptProjectInfo(projectPath, projectFilterString);
    xmlFileAccessor->setScriptActiveTab(scriptControlsTab->currentIndex(), getActiveTabName());

    // Creation Tab Specific Information
    xmlFileAccessor->setScriptImageInfo(imageFilterString, imageExtnString, imageExtnSelector->currentIndex());
    xmlFileAccessor->setScriptMeshInfo(meshFilterString, meshExtnString, meshExtnSelector->currentIndex());
    if(includeSampleNumber) {
        xmlFileAccessor->setSampleNumber(sampleNumberInt);
    }

    // constraints
    xmlFileAccessor->setCBConstraintInfo(cbConstraintIndex, fixedCBFilterString);
    if(sigmaUpToDate&&sigmaFixedButton->isChecked()) {
            double sigma[3]; sigma[0] = sigmaXEdit->text().toDouble(); sigma[1] = sigmaYEdit->text().toDouble();
            sigma[2] = sigmaZEdit->text().toDouble();
            xmlFileAccessor->setSigmaConstraint(sigma);
    }
    if(includeAveraging) {
        xmlFileAccessor->setProfileAveraging(profileAverages);
    }
    if(importFileSet) {
        xmlFileAccessor->setImportBaseName(importFileString);
    }


    // View tab specific information
    if(includeCalibration && p0ScaleSet && p1ScaleSet && p2ScaleSet  && !calibrationIsInFile) {
        xmlFileAccessor->setCalibrationParameters(calibrationScales[0], calibrationScales[1], calibrationScales[2]);
        xmlFileAccessor->setCalibrationName(calibrationTypeBox->currentIndex(), calibrationTypeBox->currentText());
    } else if(includeCalibration && calibrationIsInFile && pFileSet) {
        xmlFileAccessor->setCalibrationParameterFile(pFileEdit->text());
        xmlFileAccessor->setCalibrationName(calibrationTypeBox->currentIndex(), calibrationTypeBox->currentText());
    }

    // Processing Tab Specific Information
    int modelIndex = modelTypeSelector->currentIndex(); QString modelName = modelTypeSelector->currentText();
    int schemeIndex = schemeTypeSelector->currentIndex(); QString schemeName = schemeTypeSelector->currentText();
    int optimiserIndex = optimiserTypeSelector->currentIndex(); QString optimiserName = optimiserTypeSelector->currentText();
    xmlFileAccessor->setModelInfo(modelIndex, modelName, schemeIndex, schemeName, optimiserIndex, optimiserName);

    if(modelTypeSelector->currentIndex()==itk::CorticalBone::kHighResClassifier) {

        if(includeThresholds && stThresholdSet && cbThresholdSet && thresholdSet && !includeThresholdWeight) {
            xmlFileAccessor->setClassifierThresholdInfo(thresholds[0], thresholds[1], thresholds[2],
                                                        thresholdLevelSelector->currentIndex(),
                                                        thresholdLevelSelector->currentText());
        } else if(includeThresholds && stThresholdSet && cbThresholdSet && thresholdSet && includeThresholdWeight && thresholdWeightSet) {
            xmlFileAccessor->setClassifierThresholdInfo(thresholds[0], thresholds[1], thresholds[2],
                                                        thresholdLevelSelector->currentIndex(),
                                                        thresholdLevelSelector->currentText(), thresholdWeight);
        } else if(!includeThresholds && includeThresholdWeight && thresholdWeightSet) {
            xmlFileAccessor->setClassifierThresholdInfo(thresholdLevelSelector->currentIndex(),
                                                        thresholdLevelSelector->currentText(), thresholdWeight);
        } else {
            xmlFileAccessor->setClassifierThresholdInfo(thresholdLevelSelector->currentIndex(),
                                                        thresholdLevelSelector->currentText());
        }
    }

    if(includeSmoothing && smoothingSet) {
        xmlFileAccessor->setSmoothing(smoothingValue);
    }

    xmlFileAccessor->saveScriptFile(fileName);

    updateScript();
    updateObjectEnables();

    parentObject->updateStatusBar(tr("Script File Saved"));
    std::clog<<"<b>---- Script Saved ----</b>"<<endl;
}

void QScriptFrame::loadScript() {

    QString nameTag = parentObject->getScriptFileName();//QString(projectPath) + QDir::separator();

    // get file [.xml]
    QString fileName = QFileDialog::getOpenFileName(this, tr("Load Script File"), nameTag, tr("Script Files (*.xmls)"));

    if (fileName.isEmpty()) {
        std::clog<<"In QScriptFrame::loadScript() no file name selected so nothing opened."<<std::endl;
        return;
    }

    if(!fileName.endsWith(".xmls")) {
        fileName.append(".xmls");
    }

    bool status = loadScript(fileName);

    if(status) { parentObject->setScriptFileName(fileName); }
}

bool QScriptFrame::loadScript(QString &fileName) {

    reset();

    std::clog<<"<b>---- Load Script ----</b>"<<endl;
    QFile *file = new QFile(fileName);
    if (!file->open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"), tr("Cannot read file %1:\n%2.").arg(fileName).arg(file->errorString()));
        return false;
    }

    bool status = xmlFileAccessor->openScriptFile(fileName); // return the index saved in the file = the active pane when saved

    if(status) {
        // get general information
        QString activeTabName; int activeTabIndex;
        xmlFileAccessor->getScriptActiveTab(activeTabIndex, activeTabName);
        scriptControlsTab->setCurrentIndex(activeTabIndex);

        xmlFileAccessor->getScriptProjectInfo(projectPath, projectFilterString);
        projectFilterEdit->setText(projectFilterString); scriptDirectoryEdit->setText(projectPath);
        projectDirectorySet = true; projectFilterSet = true;

        // get creation tab information
        int imageExtnIndex, meshExtnIndex;
        xmlFileAccessor->getScriptImageInfo(imageFilterString, imageExtnString, imageExtnIndex);
        xmlFileAccessor->getScriptMeshInfo(meshFilterString, meshExtnString, meshExtnIndex);

        imageFilterEdit->setText(imageFilterString); imageExtnSelector->setCurrentIndex(imageExtnIndex);
        meshFilterEdit->setText(meshFilterString); meshExtnSelector->setCurrentIndex(meshExtnIndex);

        imageFilterSet = true; meshFilterSet = true;

        // constraints
        xmlFileAccessor->getCBConstraintInfo(cbConstraintIndex, fixedCBFilterString);
        if(cbConstraintIndex==xmlFileAccess::kCBConstraintNone) {
            cbDensityButton->setChecked(false);
        } else if(cbConstraintIndex==xmlFileAccess::kCBConstraintFixed) {
            cbDensityButton->setChecked(true); cbDensityFixedButton->setChecked(true);

            cbDensityFixedEdit->setText(fixedCBFilterString);
            cbValueFilterChanged();

        } else {
            cbDensityButton->setChecked(true); cbDensityFWHMButton->setChecked(true);
        }

        if(xmlFileAccessor->isSigmaSet()) {
            sigmaFixedButton->setChecked(true);

            double sigma[3];
            xmlFileAccessor->getSigmaConstraint(sigma);

            sigmaXEdit->setText(QString::number(sigma[0]));
            sigmaYEdit->setText(QString::number(sigma[1]));
            sigmaZEdit->setText(QString::number(sigma[2]));

            sigmaUpToDate = sigmaXSet = sigmaYSet = sigmaZSet = true;
        }

        if(xmlFileAccessor->isSampleNumberSet()) {
            sampleNumberButton->setChecked(true);  // todo remove set enabled
            sampleNumberInt = xmlFileAccessor->getSampleNumber();
            sampleNumberEdit->setText(QString::number(sampleNumberInt));
            sampleNumberSet = includeSampleNumber = true;
        }

        if(xmlFileAccessor->isSmoothingSet()) {
            smoothingButton->setChecked(true);
            smoothingValue = xmlFileAccessor->getSmoothing();
            smoothingEdit->setText(QString::number(smoothingValue));
            includeSmoothing = true; smoothingSet = true;
        }

        if(xmlFileAccessor->isAveragingSet()) {
            profileAveragingButton->setChecked(true);
            xmlFileAccessor->getProfileAveraging(profileAverages);
            xAvgEdit->setText(QString::number(profileAverages[0]));
            yAvgEdit->setText(QString::number(profileAverages[1]));
            zAvgEdit->setText(QString::number(profileAverages[2]));
            includeAveraging = xAvgSet = yAvgSet = zAvgSet = true;
        }

        if(xmlFileAccessor->isImportSet()) {
            importFileString = xmlFileAccessor->getImportBaseName();
            importFileEdit->setText(importFileString);
            importFileSet = true;
        }

        if(xmlFileAccessor->areCalibrationParametersSet()) {
            double p0, p1, p2;
            xmlFileAccessor->getCalibrationParameters(p0, p1, p2);
            calibrationScales[0]=p0; calibrationScales[1]=p1; calibrationScales[2]=p2;
            p0ScaleEdit->setText(QString::number(calibrationScales[0]));
            p1ScaleEdit->setText(QString::number(calibrationScales[1]));
            p2ScaleEdit->setText(QString::number(calibrationScales[2]));

            includeCalibration = p0ScaleSet = p1ScaleSet = p2ScaleSet = true;
            int index; QString name;
            xmlFileAccessor->getCalibrationName(index, name);
            calibrationTypeBox->setCurrentIndex(index);
        }
        if(xmlFileAccessor->isCalibrationFileSet()) {
            QString fileName;
            xmlFileAccessor->getCalibrationParameterFile(fileName);
            pFileEdit->setText(fileName);

            includeCalibration = calibrationIsInFile = pFileSet = true;
            int index; QString name;
            xmlFileAccessor->getCalibrationName(index, name);
            calibrationTypeBox->setCurrentIndex(index);
        }

        // get process tab information
        QString modelName, schemeName, optimiserName;
        int modelIndex, schemeIndex, optimiserIndex;
        xmlFileAccessor->getModelInfo(modelIndex, modelName, schemeIndex, schemeName, optimiserIndex, optimiserName);
        if(xmlFileAccessor->isClassifierThresholdIndexSet()){
            int thresholdIndex;
            if(xmlFileAccessor->areClassifierThresholdsSet()) {
                double st, cb, threshold;
                xmlFileAccessor->getClassifierThresholds(st, cb, threshold, thresholdIndex);
                thresholds[0]=st; thresholds[1]=cb; thresholds[2]=threshold;
                stThreshEdit->setText(QString::number(thresholds[0]));
                cbThreshEdit->setText(QString::number(thresholds[1]));
                thresholdEdit->setText(QString::number(thresholds[2]));
                stThresholdSet = cbThresholdSet = thresholdSet = true;
            } else if(xmlFileAccessor->isClassifierThresholdWeightSet()) {
                QString thresholdName;
                xmlFileAccessor->getClassifierThresholdIndex(thresholdIndex, thresholdName);
                xmlFileAccessor->getClassifierThresholdWeight(thresholdWeight);
                thresholdWeightEdit->setText(QString::number(thresholdWeight));
                thresholdWeightSet = true;
            } else {
                QString thresholdName;
                xmlFileAccessor->getClassifierThresholdIndex(thresholdIndex, thresholdName);
            }
            thresholdLevelSelector->setCurrentIndex(thresholdIndex); thresholdLevelChanged(thresholdIndex);
        }
        modelTypeSelector->setCurrentIndex(modelIndex); modelSelectionChanged(modelIndex);
        optimiserTypeSelector->setCurrentIndex(optimiserIndex);
        schemeTypeSelector->setCurrentIndex(schemeIndex);

    } else {
        std::clog<<"In QScriptFrame::loadScript() script index is invalid"<<endl;
    }

    updateEditorEnables();
    updateScript();
    updateObjectEnables();

    parentObject->updateStatusBar(tr("Script File Loaded"));
    std::clog<<"<b>---- Script Loaded ----</b>"<<endl;

    return status;
}

void QScriptFrame::runScript() {

    int startIndex = startSpinBox->value(), stopIndex = stopSpinBox->value();

    if(scriptControlsTab->currentIndex()==kCreateTab) {
        std::clog << "<b>---- Running Create Script ----</b>" << endl;


        for (int i = startIndex; i <= stopIndex; i++) {

            if(QFileInfo(projectNameList.at(i)).exists()) {
                std::clog << Utilities::getTabString() <<i<< "    Ignored " <<projectNameList.at(i).toStdString().c_str()<<" as already exists"<< endl;
            } else {
                QDir dir = QDir(); dir.mkdir(QFileInfo(projectNameList.at(i)).absolutePath());
                createProject(i);
            }
        }

        parentObject->updateStatusBar(tr("Finished Running Script"));
        std::clog << "<b>---- Finished Script ----</b>" << endl;

    } else if(scriptControlsTab->currentIndex()==kViewTab) {

        std::clog<<"<b>---- Running View Script ----</b>"<<endl; std::cout<<"<b>---- Running View Script ----</b>"<<endl;

        runningViewScript=true;
        updateObjectEnables();


        viewIteratorIndex = startIndex;
        viewProject(viewIteratorIndex);

    } else if(scriptControlsTab->currentIndex()==kProcessTab) {

        std::clog<<"<b>---- Running Process Script ----</b>"<<endl; std::cout<<"<b>---- Running Process Script ----</b>"<<endl;

        runningProcessScript = true;
        updateObjectEnables();

        for(int i=startIndex; i<=stopIndex; i++) {

            QString parameterResults = resultsNameList.at(i)+"_Parameters"+txtExtnString;
            QString displayResults = resultsNameList.at(i)+"_Displays"+txtExtnString;
            QString imageResults = resultsNameList.at(i)+"_ImageProfiles"+txtExtnString;
            if(QFileInfo(parameterResults).exists() && QFileInfo(displayResults).exists() && QFileInfo(imageResults).exists()) {
                std::clog << Utilities::getTabString() <<i<< "    Ignored " <<projectNameList.at(i).toStdString().c_str()<<" as already processed"<< endl;
            } else {
                processProject(i);
            }
        }

        runningProcessScript = false;

        parentObject->updateStatusBar(tr("Finished Running Script"));
        std::cout<<"<b>---- Finished Script ----</b>"<<endl;  std::clog<<"<b>---- Finished Script ----</b>"<<endl;

    }

    updateObjectEnables();

}

void QScriptFrame::nextProjectFile() {

    parentObject->closeProject(); // don't save

    int stopIndex = stopSpinBox->value();

    viewIteratorIndex++;
    if(viewIteratorIndex<=stopIndex) {
        viewProject(viewIteratorIndex);
    } else {
        runningViewScript=false;
        parentObject->updateStatusBar(tr("Finished Running Script"));
        std::cout<<"<b>---- Finished Script ----</b>"<<endl;  std::clog<<"<b>---- Finished Script ----</b>"<<endl;
    }
    updateObjectEnables();
}

//--- tabs ---//
void QScriptFrame::tabChanged() {
    updateScript(false); // false so ranges not changed
    updateObjectEnables();
}

QString QScriptFrame::getActiveTabName() {
    if(scriptControlsTab->currentIndex()==kCreateTab) {
        return QString("Creation Tab");
    } else if(scriptControlsTab->currentIndex()==kViewTab) {
        return QString("View Tab");
    } else if(scriptControlsTab->currentIndex()==kProcessTab) {
        return QString("Processing Tab");
    } else {
        return QString("Invalid Tab");
    }
}