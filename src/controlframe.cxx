#include <QtWidgets>
#include <iostream>
#include <istream>
#include <sstream>

#include "controlframe.h"
#include "corticalbone.h"




QControlFrame::QControlFrame(MainWindow *parent) : QScrollArea(parent) {

    createObjects();
    createConnections();
    initialiseSettings();

    // initilise button enable/disabling
    setState(0);

    parentObject = parent;

}

void QControlFrame::initialiseSettings() {

    // measurement settings
    calibrationPointsButton->setChecked(false);
    modelFunctionSelector->setCurrentIndex(0);
    optimiserSelector->setCurrentIndex(0);
    fittingSelector->setCurrentIndex(0); measurementsUpToDate = false;
    sigmaRequired = smoothingRequired = measurementsRequired = averagingRequired = false;

    // calibration settings
    ptMeasureButton->setChecked(false);
    calibrationPhantomSelector->setCurrentIndex(0);
    p0Edit->setText(QString()); p1Edit->setText(QString()); p2Edit->setText(QString());
    calScalesUpToDate = true; calUpToDate = true; p0CalSet = p1CalSet = p2CalSet = false; calRadiusSet = false;
    calRadiusEdit->clear();
    calRadiusEdit->setHidden(true); calRadiusLabel->setHidden(true);
    calNumberEdit->clear();
    calNumberEdit->setHidden(true); calNumberLabel->setHidden(true);

    // display settings
    displayParameterSelector->setCurrentIndex(-1);
    displayMeshSelector->setCurrentIndex(0);

    // sample number
    sampleNumberButton->setChecked(false); sampleNumberEdit->setText(QString());
    sampleNumberUpToDate = true; sampleNumberSet = false; // set state after reseting GUI values

    // smoothing
    smoothingButton->setChecked(false); smoothingEdit->setText(QString());
    smoothingUpToDate = true; smoothingSet = false;

    // cb density constraint
    cbDensityButton->setChecked(false); cbDensityFWHMButton->setChecked(true);
    cbDensityFixedEdit->setText(QString());
    cBMDUpToDate = true;  cbDensitySet = false;

    // sigma constraint
    sigmaFixedButton->setChecked(false);
    sigmaXEdit->setText(QString()); sigmaYEdit->setText(QString()); sigmaZEdit->setText(QString());
    sigmaUpToDate = true; sigmaXSet = sigmaYSet = sigmaZSet = false;
    fittingMode = true; calibrationMode = classifierMode = false;

    // thresholds constraints
    stThresholdEdit->setText(QString()); cbThresholdEdit->setText(QString()); thresholdEdit->setText(QString());
    modelMode = itk::CorticalBone::kThreeTierRect;

    // profile averaging
    profileAveragingButton->setChecked(false);
    xAveragingEdit->setText(QString()); yAveragingEdit->setText(QString()); zAveragingEdit->setText(QString());
    averagingUpToDate = true; xAveragingSet = yAveragingSet = zAveragingSet = false;

    // import
    importPathEdit->setText(QString());

    // adjust settings
    modelFunctionSelector->setCurrentIndex(0); // todo remove double up

}

void QControlFrame::createObjects() {

    // measure CBM
    ptMeasureButton = new  QCheckBox(tr("&Point Measurements"));
    runOverMeshButton = new QPushButton(tr("&Run Over Mesh"));

    modelFunctionSelector = new QComboBox(this); // Model function used to fit the model
    modelFunctionSelector->addItem("Three Tier Rect");
    modelFunctionSelector->addItem("Endosteal Ramp");
    modelFunctionSelector->addItem("HR Classifier");
    modelFunctionSelector->addItem("Calibration");
    modelFunctionSelector->setCurrentIndex(0);

    optimiserSelector = new QComboBox(this); // Optimiser used to fit the model
    optimiserSelector->addItem("LM Optimiser");
    optimiserSelector->addItem("Powell Optimiser");
    optimiserSelector->addItem("Evolutionary Optimiser");
    optimiserSelector->setCurrentIndex(0);

    fittingSelector = new QComboBox(this); // fitting scheme used to fit the model
    fittingSelector->addItem("std A"); // rect 1st then ramp if ramp selected - hueristic weighting
    fittingSelector->addItem("std B"); // rect 1st then ramp if ramp selected - narrow weighting
    fittingSelector->addItem("std C"); // rect 1st then ramp if ramp selected - wide weighting
    fittingSelector->addItem("std D"); // rect 1st then ramp if ramp selected - narrow to wide weighting
    fittingSelector->addItem("CBMV2a"); // (selected model fit 2x + sigma correction) - static huristic weighting
    fittingSelector->addItem("CBMV2b"); // (selected model fit 2x + sigma correction) - static narrow weighting
    fittingSelector->addItem("CBMV2c"); // (selected model fit 2x + sigma correction) - static narrow then wide weighting
    fittingSelector->addItem("CBMV2d"); // (selected model fit 2x + sigma correction) - dynamic narrow then wide weighting
    fittingSelector->addItem("Unconstrained a"); // rect followed by ramp/rect followed by ramp/rect no constraints
    fittingSelector->addItem("Unconstrained b");
    fittingSelector->addItem("Unconstrained c");
    fittingSelector->addItem("Smoothing a"); // V4 then smooth the CB estimates then V1
    fittingSelector->addItem("Smoothing b");
    fittingSelector->addItem("Smoothing c");
    fittingSelector->addItem("Smoothing d");
    fittingSelector->setCurrentIndex(0);

    QGroupBox *cmbBox = new QGroupBox(tr("CBM Measures"));
    QGridLayout *cbmLayout = new QGridLayout();
    cbmLayout->addWidget(ptMeasureButton, 0,0);
    cbmLayout->addWidget(runOverMeshButton, 0,1, 1,2);
    cbmLayout->addWidget(modelFunctionSelector, 1,0);
    cbmLayout->addWidget(optimiserSelector, 1,1);
    cbmLayout->addWidget(fittingSelector, 1,2);
    cmbBox->setLayout(cbmLayout); //cmbBox->setContentsMargins(0,0,0,0); // left, top, right, bottom


    // calibration // todo consider text = QLabel(tr("Quadratically Scale Imported Densities [y=p<sub>0</sub>+p<sub>1</sub>x+p<sub>2</sub>x<sup>2</sup>]"));
    calibrationPointsButton = new  QCheckBox(tr("&Set Calibration Points"));
    calibrateButton = new QPushButton(tr("&Calibrate Phantom"));

    calibrationPhantomSelector = new QComboBox(this);
    calibrationPhantomSelector->addItem("No Phantom");
    calibrationPhantomSelector->addItem("Mindways Solid Phantom");
    calibrationPhantomSelector->addItem("Bone Density Phantom");
    calibrationPhantomSelector->addItem("European Spine Phantom");
    calibrationPhantomSelector->addItem("Manual Control Points");
    calibrationPhantomSelector->addItem("Manual Linear");
    calibrationPhantomSelector->addItem("Manual Quadratic");
    calibrationPhantomSelector->setCurrentIndex(0);

    p0Label = new QLabel(tr("P<sub>0</sub>")); p1Label = new QLabel(tr("P<sub>1</sub>")); p2Label = new QLabel(tr("P<sub>2</sub>"));
    p0Edit = new QLineEdit; p1Edit = new QLineEdit; p2Edit = new QLineEdit;
    p0Edit->setValidator( new QDoubleValidator(-5000, 5000, 9, this) );
    p1Edit->setValidator( new QDoubleValidator(-5000, 5000, 9, this) );
    p2Edit->setValidator( new QDoubleValidator(-5000, 5000, 9, this) );

    calRadiusLabel = new QLabel(tr("Radius"));
    calRadiusEdit = new QLineEdit; calRadiusEdit->setValidator( new QDoubleValidator(0, 5, 5, this) );

    calNumberLabel = new QLabel(tr("#"));
    calNumberEdit = new QLineEdit; calNumberEdit->setValidator( new QIntValidator(1, 100, this) );

    calibrationBox = new QGroupBox(tr("BMD Calibration"));
    QGridLayout *calibrationLayout = new QGridLayout();
    calibrationLayout->addWidget(calibrationPointsButton, 0,0, 1,3);
    calibrationLayout->addWidget(calibrateButton, 0,3, 1,3);
    calibrationLayout->addWidget(calibrationPhantomSelector, 1,0, 1,6);
    calibrationLayout->addWidget(p0Label, 2,0); calibrationLayout->addWidget(p0Edit, 2,1);
    calibrationLayout->addWidget(p1Label, 2,2); calibrationLayout->addWidget(p1Edit, 2,3);
    calibrationLayout->addWidget(p2Label, 2,4); calibrationLayout->addWidget(p2Edit, 2,5);
    calibrationLayout->addWidget(calRadiusLabel, 3,0, 1, 3); calibrationLayout->addWidget(calRadiusEdit, 3,3);
    calibrationLayout->addWidget(calNumberLabel, 3,4); calibrationLayout->addWidget(calNumberEdit, 3,5);
    calibrationBox->setLayout(calibrationLayout); //calibrationBox->setContentsMargins(0,0,0,0); // left, top, right, bottom


    // display
    displayMeshSelector = new QComboBox(this); // select mesh or volume
    displayMeshSelector->addItem("Mesh");
    displayMeshSelector->addItem("Periosteal");
    displayMeshSelector->addItem("Volume");

    displayParameterSelector = new QComboBox(this); // select parameters
    displayParameterSelector->addItem("Mean Cortical Thickness");
    displayParameterSelector->addItem("Dense Cortical Thickness");
    displayParameterSelector->addItem("Endosteal Thickness");
    displayParameterSelector->addItem("Periosteal Position");
    displayParameterSelector->addItem("Endosteal Position");
    displayParameterSelector->addItem("Endosteal CB Position");
    displayParameterSelector->addItem("Endosteal TB Position");
    displayParameterSelector->addItem("Cortical Density");
    displayParameterSelector->addItem("Trabecular Density");
    displayParameterSelector->addItem("Soft Tissue Density");
    displayParameterSelector->addItem("Mass SA");
    displayParameterSelector->addItem("Sigma");
    displayParameterSelector->addItem("Mean Error (ABS)");
    displayParameterSelector->setCurrentIndex(-1);

    displayBox = new QGroupBox(tr("Display"));
    QHBoxLayout *displayLayout = new QHBoxLayout();
    displayLayout->addWidget(displayMeshSelector);
    displayLayout->addWidget(displayParameterSelector);
    displayBox->setLayout(displayLayout); //displayBox->setContentsMargins(0,0,0,0); // left, top, right, bottom


    // model Constraints
    cbDensityButton = new QCheckBox(); cbDensityLabel = new QLabel(tr("Y<sub>CB</sub>"));
    cbDensityFixedButton = new QRadioButton(tr("Fixed"));
    cbDensityFWHMButton = new QRadioButton(tr("FWHM"));

    cbDensityFixedEdit = new QLineEdit;
    cbDensityFixedEdit->setValidator( new QDoubleValidator(0, 10000, 6, this) ); // TODO - ModelTransformBase::MAX_BMD

    sigmaFixedButton = new QCheckBox(tr("σ"));
    sigmaXLabel = new QLabel(tr("X")); sigmaYLabel = new QLabel(tr("Y")); sigmaZLabel = new QLabel(tr("Z"));
    sigmaXEdit = new QLineEdit; sigmaYEdit = new QLineEdit; sigmaZEdit = new QLineEdit;
    sigmaXEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );
    sigmaYEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );
    sigmaZEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );

    constraintsBox = new QGroupBox(tr("Model Constraints"));
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
    constraintsBox->setLayout(constraintsLayout); //cbDensityBox->setContentsMargins(0,0,0,0);

    // model thresholds
    thresholdsSelector = new QComboBox(this);
    thresholdsSelector->addItem("Very High Threshold [99.99%]");
    thresholdsSelector->addItem("High Threshold [99.9%]");
    thresholdsSelector->addItem("Medium Threshold [99.0%]");
    thresholdsSelector->addItem("Low Threshold [95.0%]");
    thresholdsSelector->addItem("Manually Set");
    thresholdsSelector->addItem("Median midpoint");
    thresholdsSelector->addItem("Median manual");
    calculateThresholdsButton = new  QPushButton(tr("CalculateThresholds"));
    thresholdSTLabel = new QLabel(tr("ρ<sub>st</sub>")); thresholdCBLabel = new QLabel(tr("ρ<sub>cb</sub>")); thresholdLabel = new QLabel(tr("ρ<sub>thresh</sub>"));
    stThresholdEdit = new QLineEdit; cbThresholdEdit = new QLineEdit; thresholdEdit = new QLineEdit;
    stThresholdEdit->setValidator( new QDoubleValidator(-20000, 20000, 6, this) );
    cbThresholdEdit->setValidator(new QDoubleValidator(-20000, 20000, 6, this) );
    thresholdEdit->setValidator( new QDoubleValidator(-20000, 20000, 6, this) );

    thresholdWeightingLabel = new QLabel(tr("Median weighting"));
    thresholdWeightingEdit = new QLineEdit;
    thresholdWeightingEdit->setValidator( new QDoubleValidator(0, 1, 3, this) );

    thresholdsBox = new QGroupBox(tr("Model Thresholds"));
    QGridLayout *thresholdsLayout = new QGridLayout;
    thresholdsLayout->addWidget(thresholdsSelector, 0,0, 1,4);
    thresholdsLayout->addWidget(calculateThresholdsButton, 0,4, 1,2);
    thresholdsLayout->addWidget(thresholdWeightingLabel, 1,0, 1,2);
    thresholdsLayout->addWidget(thresholdWeightingEdit, 1,2, 1,4);
    thresholdsLayout->addWidget(thresholdSTLabel, 2, 0); thresholdsLayout->addWidget(stThresholdEdit, 2, 1);
    thresholdsLayout->addWidget(thresholdCBLabel, 2, 2); thresholdsLayout->addWidget(cbThresholdEdit, 2, 3);
    thresholdsLayout->addWidget(thresholdLabel, 2,4); thresholdsLayout->addWidget(thresholdEdit, 2,5);
    thresholdsBox->setLayout(thresholdsLayout); //cbDensityBox->setContentsMargins(0,0,0,0);


    // profile settings - sample #
    sampleNumberButton = new  QCheckBox(tr("#"));
    sampleNumberEdit = new QLineEdit;
    sampleNumberEdit->setValidator( new QIntValidator(0, 1000, this) );

    // profile averaging
    profileAveragingButton = new QCheckBox(); profileAveragingLabel = new QLabel(tr("<b>x\u0305</b>"));
    profileAveragingButton->setChecked(false);

    xAveragingEdit = new QLineEdit; yAveragingEdit = new QLineEdit; zAveragingEdit = new QLineEdit;
    xAveragingEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );
    yAveragingEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );
    zAveragingEdit->setValidator( new QDoubleValidator(0, 10, 6, this) );

    xAveragingLabel = new QLabel(tr("X")); yAveragingLabel = new QLabel(tr("Y")); zAveragingLabel = new QLabel(tr("Z"));

    // smoothing
    smoothingButton = new  QCheckBox(tr("Smoothing σ"));
    smoothingEdit = new QLineEdit;
    smoothingEdit->setValidator( new QDoubleValidator(0, 10, 2, this) );

    // layout
    profileSettingsBox = new QGroupBox(tr("Profile Settings"));
    QGridLayout *profileSettingsLayout = new QGridLayout;
    profileSettingsLayout->addWidget(sampleNumberButton, 0, 0, 1,2);
    profileSettingsLayout->addWidget(sampleNumberEdit, 0, 2, 1,8);
    profileSettingsLayout->addWidget(profileAveragingButton, 1,0); profileSettingsLayout->addWidget(profileAveragingLabel, 1,1);
    profileSettingsLayout->addWidget(xAveragingLabel, 1,2); profileSettingsLayout->addWidget(xAveragingEdit, 1,3, 1,2);
    profileSettingsLayout->addWidget(yAveragingLabel, 1,5); profileSettingsLayout->addWidget(yAveragingEdit, 1,6, 1,2);
    profileSettingsLayout->addWidget(zAveragingLabel, 1,8); profileSettingsLayout->addWidget(zAveragingEdit, 1,9);
    profileSettingsLayout->addWidget(smoothingButton, 2, 0, 1,4);
    profileSettingsLayout->addWidget(smoothingEdit, 2, 4, 1,6);
    profileSettingsBox->setLayout(profileSettingsLayout);  //hrAveragingBox->setContentsMargins(0,0,0,0); // left, top, right, bottom


    // import profiles / parameters
    setImportButton = new  QPushButton(tr("Set Import"));
    removeImportButton = new  QToolButton; removeImportButton->setIcon(QIcon(":/images/clearText"));
    removeImportButton->setToolButtonStyle(Qt::ToolButtonIconOnly);
    importPathEdit = new QLineEdit;
    importPathEdit->setEnabled(false); importPathEdit->setReadOnly(true);

    importBox = new QGroupBox(tr("Import External Profiles and Parameters"));
    QGridLayout *importLayout = new QGridLayout;
    importLayout->addWidget(setImportButton, 0,0, 1,2);
    importLayout->addWidget(importPathEdit, 0,2, 1,3);
    importLayout->addWidget(removeImportButton, 0,5);
    importBox->setLayout(importLayout);  //importBox->setContentsMargins(0,0,0,0); // left, top, right, bottom


    // -------------- overall layout -----------------//
    QVBoxLayout *overall_layout = new QVBoxLayout;
    overall_layout->addWidget(cmbBox);
    overall_layout->addWidget(displayBox);
    overall_layout->addWidget(calibrationBox);
    overall_layout->addWidget(thresholdsBox);
    overall_layout->addWidget(constraintsBox);
    overall_layout->addWidget(profileSettingsBox);
    overall_layout->addWidget(importBox);
    overall_layout->setSizeConstraint(QLayout::SetFixedSize);

    QGroupBox *overallBox = new QGroupBox();
    overallBox->setLayout(overall_layout);
    //overallBox->setContentsMargins(0,0,0,0);

    setWidget(overallBox);

    //setLayout(overall_layout);

}

void QControlFrame::updateFittingSelectorOptions() {

    int numberOfOptions = fittingSelector->count();
    for(int i=numberOfOptions-1; i>=0; i--) {
        fittingSelector->removeItem(i);
    }

    if(modelMode<=itk::CorticalBone::kEndostealRamp) {
        fittingSelector->addItem("std A"); // rect 1st then ramp if ramp selected - hueristic weighting
        fittingSelector->addItem("std B"); // rect 1st then ramp if ramp selected - narrow weighting
        fittingSelector->addItem("std C"); // rect 1st then ramp if ramp selected - wide weighting
        fittingSelector->addItem("std D"); // rect 1st then ramp if ramp selected - narrow to wide weighting
        fittingSelector->addItem("CBMV2a"); // (selected model fit 2x + sigma correction) - static huristic weighting
        fittingSelector->addItem("CBMV2b"); // (selected model fit 2x + sigma correction) - static narrow weighting
        fittingSelector->addItem("CBMV2c"); // (selected model fit 2x + sigma correction) - static narrow then wide weighting
        fittingSelector->addItem("CBMV2d"); // (selected model fit 2x + sigma correction) - dynamic narrow then wide weighting
        fittingSelector->addItem("Unconstrained a"); // rect followed by ramp/rect followed by ramp/rect no constraints
        fittingSelector->addItem("Unconstrained b");
        fittingSelector->addItem("Unconstrained c");
        fittingSelector->addItem("Smoothing a"); // V4 then smooth the CB estimates then V1
        fittingSelector->addItem("Smoothing b");
        fittingSelector->addItem("Smoothing c");
        fittingSelector->addItem("Smoothing d");
        fittingSelector->setCurrentIndex(0);

    } else if(modelMode==itk::CorticalBone::kHighResClassifier) {
        fittingSelector->addItem("Global Thresholds");
        fittingSelector->addItem("Median Thresholds");
        fittingSelector->setCurrentIndex(0);

    }else if(modelMode==itk::CorticalBone::kCalibration) {
        fittingSelector->addItem("Median Upper");
        fittingSelector->addItem("Maximum Upper");
        fittingSelector->setCurrentIndex(0);

    }
//    else if(modelMode==itk::CorticalBone::kHighResOptimiser) {
//        fittingSelector->addItem("Endosteal Edges");
//        fittingSelector->setCurrentIndex(0);
//    }
}

void QControlFrame::updateDisplaySelectorOptions() {

    int numberOfOptions = displayParameterSelector->count();
    for(int i=numberOfOptions-1; i>=0; i--) {
        displayParameterSelector->removeItem(i);
    }

    if(modelMode<=itk::CorticalBone::kEndostealRamp) {
        displayParameterSelector->addItem("Total Thickness");
        displayParameterSelector->addItem("Cortical Thickness");
        displayParameterSelector->addItem("Endosteal Thickness");
        displayParameterSelector->addItem("Periosteal Position");
        displayParameterSelector->addItem("Endosteal Position");
        displayParameterSelector->addItem("Endosteal CB Position");
        displayParameterSelector->addItem("Endosteal TB Position");
        displayParameterSelector->addItem("Cortical Density");
        displayParameterSelector->addItem("Trabecular Density");
        displayParameterSelector->addItem("Soft Tissue Density");
        displayParameterSelector->addItem("Mass SA");
        displayParameterSelector->addItem("Sigma");
        displayParameterSelector->addItem("Mean Error (ABS)");

        if(!importPathEdit->text().isEmpty()) {
            displayParameterSelector->addItem("Imported vs Local Image Profile Diff"); // actually the Mean Error between the imported image profiles and the local image profiles
            displayParameterSelector->addItem("Imported vs Local Image Profile STD"); // actually the Mean Error between the imported image profiles and the local image profiles
            displayParameterSelector->addItem("Imported Image vs Local Model Error"); // actually the Mean Error between the imported image profiles and the local blurred model profiles
        }

    } else if(modelMode==itk::CorticalBone::kHighResClassifier) {

        displayParameterSelector->addItem("Not Bone Threshold");
        displayParameterSelector->addItem("Cortical Bone Threshold");
        displayParameterSelector->addItem("Threshold Density");
        displayParameterSelector->addItem("Macro Classification");
        displayParameterSelector->addItem("Not Bone Density");
        displayParameterSelector->addItem("Cortical Bone Density");
        displayParameterSelector->addItem("Dense Bone Density");
        displayParameterSelector->addItem("Cortical + Endo Density");
        displayParameterSelector->addItem("Trabecular Bone Density");
        displayParameterSelector->addItem("Cortical Thickness");
        displayParameterSelector->addItem("Endocortical Thickness");
        displayParameterSelector->addItem("Endocortical Slope");
        displayParameterSelector->addItem("Porosity");
        displayParameterSelector->addItem("# Profiles");

    } else if(modelMode==itk::CorticalBone::kCalibration) {

        displayParameterSelector->addItem("Lower Value");
        displayParameterSelector->addItem("Upper Value");
    }
//    else if(modelMode==itk::CorticalBone::kHighResOptimiser) {
//
//        displayParameterSelector->addItem("Periosteal Edge");
//        displayParameterSelector->addItem("Endosteal Edge Start");
//        displayParameterSelector->addItem("Endosteal Edge End");
//        displayParameterSelector->addItem("Soft Tissue Linearity");
//        displayParameterSelector->addItem("Cortical Bone Linearity");
//        displayParameterSelector->addItem("Endosteal Linearity");
//        displayParameterSelector->addItem("Trabecular Linearity");
//        displayParameterSelector->addItem("Percentage Soft Tissue");
//        displayParameterSelector->addItem("Percentage Cortical Bone");
//        displayParameterSelector->addItem("Percentage Trabecular Bone");
//        displayParameterSelector->addItem("Percentage Medullary");
//    }

    displayParameterSelector->setCurrentIndex(-1);
}

void QControlFrame::createConnections() {  // NOTE - slots called in the order in which the connections were made

    // measurements
    connect(runOverMeshButton, SIGNAL(clicked()), this->parent(), SLOT(runModellingOverMesh()));

    connect(ptMeasureButton, SIGNAL(clicked(bool)), this->parent(), SLOT(togglePtMeasures(bool)));

    connect(modelFunctionSelector , SIGNAL(currentIndexChanged(int)),this,SLOT(modelFunctionChanged(int)));
    connect(modelFunctionSelector , SIGNAL(currentIndexChanged(int)),this->parent(),SLOT(modelFunctionChanged(int)));

    connect(optimiserSelector , SIGNAL(currentIndexChanged(int)),this->parent(),SLOT(optimiserChanged(int)));

    connect(fittingSelector, SIGNAL(currentIndexChanged(int)), this, SLOT(fittingSchemeChanged(int)));
    connect(fittingSelector , SIGNAL(currentIndexChanged(int)),this->parent(),SLOT(fittingSchemeChanged(int)));


    // displays
    connect(displayParameterSelector , SIGNAL(currentIndexChanged(int)),this->parent(),SLOT(displayParameterChanged(int)));
    connect(displayMeshSelector , SIGNAL(currentIndexChanged(int)),this->parent(),SLOT(displayMeshChanged(int)));


    // calibration
    connect(calibrateButton, SIGNAL(clicked()), this->parent(), SLOT(runCalibration()));

    connect(calibrationPointsButton, SIGNAL(clicked(bool)), this->parent(), SLOT(toggleCalibrationMode(bool)));
    connect(calibrationPhantomSelector , SIGNAL(currentIndexChanged(int)),this->parent(),SLOT(calibrationPhantomChanged(int)));
    connect(calibrationPhantomSelector , SIGNAL(currentIndexChanged(int)),this,SLOT(calibrationPhantomChanged(int)));

    connect(calRadiusEdit, SIGNAL(returnPressed()), this, SLOT(calRadiusChanged()));
    connect(calRadiusEdit, SIGNAL(returnPressed()), this->parent(), SLOT(updateCalibrationRadius()));
    connect(calRadiusEdit, SIGNAL(textChanged(const QString &)), this, SLOT(calRadiusChanging()));
    connect(calRadiusEdit, SIGNAL(textChanged(const QString &)), this->parent(), SLOT(disableCalibration()));

    connect(calNumberEdit, SIGNAL(returnPressed()), this, SLOT(calNumberChanged()));
    connect(calNumberEdit, SIGNAL(returnPressed()), this->parent(), SLOT(updateCalibrationNumber()));
    connect(calNumberEdit, SIGNAL(textChanged(const QString &)), this, SLOT(calNumberChanging()));
    connect(calNumberEdit, SIGNAL(textChanged(const QString &)), this->parent(), SLOT(disableCalibration()));

    connect(p0Edit , SIGNAL(returnPressed()),this,SLOT(p0ValueChanged()));
    connect(p1Edit , SIGNAL(returnPressed()),this,SLOT(p1ValueChanged()));
    connect(p2Edit , SIGNAL(returnPressed()),this,SLOT(p2ValueChanged()));

    connect(p0Edit , SIGNAL(returnPressed()),this->parent(),SLOT(updatePtMeasures())); //editing of one finished
    connect(p1Edit , SIGNAL(returnPressed()),this->parent(),SLOT(updatePtMeasures())); //editing of one finished
    connect(p2Edit , SIGNAL(returnPressed()),this->parent(),SLOT(updatePtMeasures())); //editing of one finished

    connect(p0Edit, SIGNAL(textChanged(const QString &)), this, SLOT(p0ValueChanging()));
    connect(p1Edit, SIGNAL(textChanged(const QString &)), this, SLOT(p1ValueChanging()));
    connect(p2Edit, SIGNAL(textChanged(const QString &)), this, SLOT(p2ValueChanging()));

    connect(p0Edit , SIGNAL(textChanged(const QString &)),this->parent(),SLOT(disablePtMeasures())); //editing of one finished
    connect(p1Edit , SIGNAL(textChanged(const QString &)),this->parent(),SLOT(disablePtMeasures())); //editing of one finished
    connect(p2Edit , SIGNAL(textChanged(const QString &)),this->parent(),SLOT(disablePtMeasures())); //editing of one finished


    // sample number setting
    connect(sampleNumberButton, SIGNAL(clicked(bool)), this, SLOT(toggleSampleNumber(bool)));
    connect(sampleNumberButton, SIGNAL(clicked(bool)), this->parent(), SLOT(toggleSampleNumber(bool)));

    connect(sampleNumberEdit, SIGNAL(returnPressed()), this, SLOT(sampleNumberChanged()));
    connect(sampleNumberEdit, SIGNAL(returnPressed()), this->parent(), SLOT(sampleNumberChanged()));

    connect(sampleNumberEdit, SIGNAL(textChanged(const QString &)), this, SLOT(sampleNumberChanging()));
    connect(sampleNumberEdit, SIGNAL(textChanged(const QString &)), this->parent(), SLOT(disablePtMeasures()));

    // smoothing setting
    connect(smoothingButton, SIGNAL(clicked(bool)), this, SLOT(toggleSmoothing(bool)));
    connect(smoothingButton, SIGNAL(clicked(bool)), this->parent(), SLOT(toggleSmoothing(bool)));

    connect(smoothingEdit, SIGNAL(returnPressed()), this, SLOT(smoothingChanged()));
    connect(smoothingEdit, SIGNAL(returnPressed()), this->parent(), SLOT(smoothingChanged()));

    connect(smoothingEdit, SIGNAL(textChanged(const QString &)), this, SLOT(smoothingChanging()));
    connect(smoothingEdit, SIGNAL(textChanged(const QString &)), this->parent(), SLOT(disablePtMeasures()));


    // set fixed CB Density connections
    connect(cbDensityFixedEdit , SIGNAL(returnPressed()),this,SLOT(fixedCBDensityChanged())); //editingFinished
    connect(cbDensityFixedEdit , SIGNAL(returnPressed()),this->parent(),SLOT(fixedCBDensityChanged())); //editingFinished

    connect(cbDensityFixedEdit , SIGNAL(textChanged(const QString &)),this,SLOT(fixedCBDensityChanging()));
    connect(cbDensityFixedEdit, SIGNAL(textChanged(const QString &)), this->parent(), SLOT(disablePtMeasures()));

    connect(cbDensityButton, SIGNAL(clicked(bool)), this, SLOT(toggleCBConstraintDensity(bool)));
    connect(cbDensityButton, SIGNAL(clicked(bool)), this->parent(), SLOT(toggleCBDensityConstraint(bool)));

    connect(cbDensityFixedButton, SIGNAL(clicked(bool)), this, SLOT(turnOnFixedCBDensity()));
    connect(cbDensityFixedButton, SIGNAL(clicked(bool)), this->parent(), SLOT(turnOnFixCBDensity()));

    connect(cbDensityFWHMButton, SIGNAL(clicked(bool)), this, SLOT(turnOnFWHMCBDensity()));
    connect(cbDensityFWHMButton, SIGNAL(clicked(bool)), this->parent(), SLOT(turnOnFWHMCBDensity()));


    // set fixed CB Density connections
    connect(sigmaXEdit , SIGNAL(returnPressed()),this, SLOT(fixedSigmaXChanged()));
    connect(sigmaXEdit , SIGNAL(returnPressed()),this->parent(),SLOT(sigmaChanged())); //editingFinished

    connect(sigmaXEdit, SIGNAL(textChanged(const QString &)), this, SLOT(fixedSigmaXChanging()));
    connect(sigmaXEdit, SIGNAL(textChanged(const QString &)), this->parent(), SLOT(disablePtMeasures()));

    connect(sigmaYEdit , SIGNAL(returnPressed()),this, SLOT(fixedSigmaYChanged()));
    connect(sigmaYEdit , SIGNAL(returnPressed()),this->parent(),SLOT(sigmaChanged())); //editingFinished

    connect(sigmaYEdit, SIGNAL(textChanged(const QString &)), this, SLOT(fixedSigmaYChanging()));
    connect(sigmaYEdit, SIGNAL(textChanged(const QString &)), this->parent(), SLOT(disablePtMeasures()));

    connect(sigmaZEdit , SIGNAL(returnPressed()),this, SLOT(fixedSigmaZChanged()));
    connect(sigmaZEdit , SIGNAL(returnPressed()),this->parent(),SLOT(sigmaChanged())); //editingFinished

    connect(sigmaZEdit, SIGNAL(textChanged(const QString &)), this, SLOT(fixedSigmaZChanging()));
    connect(sigmaZEdit, SIGNAL(textChanged(const QString &)), this->parent(), SLOT(disablePtMeasures()));

    connect(sigmaFixedButton, SIGNAL(clicked(bool)), this, SLOT(toggleFixedSigma(bool)));
    connect(sigmaFixedButton, SIGNAL(clicked(bool)), this->parent(), SLOT(toggleFixedSigma(bool)));

    // set thresholds connections
    connect(thresholdsSelector , SIGNAL(currentIndexChanged(int)),this->parent(),SLOT(thresholdSelectorChanged(int)));
    connect(calculateThresholdsButton , SIGNAL(clicked()),this->parent(),SLOT(runThresholdCalculation()));

    connect(stThresholdEdit, SIGNAL(returnPressed()), this, SLOT(stThresholdChanged()));
    connect(cbThresholdEdit, SIGNAL(returnPressed()), this, SLOT(cbThresholdChanged()));
    connect(thresholdEdit , SIGNAL(returnPressed()),this,SLOT(thresholdChanged()));

    connect(stThresholdEdit, SIGNAL(textChanged(const QString &)), this, SLOT(stThresholdChanging()));
    connect(cbThresholdEdit, SIGNAL(textChanged(const QString &)), this, SLOT(cbThresholdChanging()));
    connect(thresholdEdit, SIGNAL(textChanged(const QString &)), this, SLOT(thresholdChanging()));

    connect(stThresholdEdit , SIGNAL(textChanged(const QString &)),this->parent(),SLOT(disablePtMeasures())); //editing of one finished
    connect(cbThresholdEdit , SIGNAL(textChanged(const QString &)),this->parent(),SLOT(disablePtMeasures())); //editing of one finished
    connect(thresholdEdit , SIGNAL(textChanged(const QString &)),this->parent(),SLOT(disablePtMeasures())); //editing of one finished

    connect(thresholdWeightingEdit, SIGNAL(returnPressed()), this, SLOT(thresholdWeightingChanged()));
    connect(thresholdWeightingEdit, SIGNAL(textChanged(const QString &)), this, SLOT(thresholdWeightingChanging()));
    connect(thresholdWeightingEdit , SIGNAL(textChanged(const QString &)),this->parent(),SLOT(disablePtMeasures()));

    // set HR profile averaging connections
    connect(xAveragingEdit , SIGNAL(returnPressed()),this,SLOT(setXAveragingValue()));
    connect(yAveragingEdit , SIGNAL(returnPressed()),this,SLOT(setYAveragingValue()));
    connect(zAveragingEdit , SIGNAL(returnPressed()),this,SLOT(setZAveragingValue()));

    connect(xAveragingEdit , SIGNAL(returnPressed()),this->parent(),SLOT(ProfileAveragingValuesChanged())); //editing of one finished
    connect(yAveragingEdit , SIGNAL(returnPressed()),this->parent(),SLOT(ProfileAveragingValuesChanged())); //editing of one finished
    connect(zAveragingEdit , SIGNAL(returnPressed()),this->parent(),SLOT(ProfileAveragingValuesChanged())); //editing of one finished

    connect(xAveragingEdit, SIGNAL(textChanged(const QString &)), this, SLOT(xAveragingChanging()));
    connect(yAveragingEdit, SIGNAL(textChanged(const QString &)), this, SLOT(yAveragingChanging()));
    connect(zAveragingEdit, SIGNAL(textChanged(const QString &)), this, SLOT(zAveragingChanging()));

    connect(xAveragingEdit , SIGNAL(textChanged(const QString &)),this->parent(),SLOT(disablePtMeasures()));
    connect(yAveragingEdit , SIGNAL(textChanged(const QString &)),this->parent(),SLOT(disablePtMeasures()));
    connect(zAveragingEdit , SIGNAL(textChanged(const QString &)),this->parent(),SLOT(disablePtMeasures()));

    connect(profileAveragingButton, SIGNAL(clicked(bool)), this, SLOT(toggleProfileAveraging(bool)));
    connect(profileAveragingButton, SIGNAL(clicked(bool)), this->parent(), SLOT(toggleProfileAverging(bool)));


    // Import Parameters
    connect(setImportButton, SIGNAL(clicked(bool)), this->parent(), SLOT(importProfile()));
    connect(removeImportButton, SIGNAL(clicked(bool)), this->parent(), SLOT(removeImportedProfile()));

}

void QControlFrame::setCalibrationPhantom(int index) {
    calibrationPhantomSelector->setCurrentIndex(index);
}

void QControlFrame::setFixedSampleNumber(QString value) {
    sampleNumberButton->setChecked(true);
    sampleNumberEdit->setText(value);

    sampleNumberSet = sampleNumberUpToDate = true;
}

void QControlFrame::setSmoothingRadius(QString value) {
    smoothingButton->setChecked(true);
    smoothingEdit->setText(value);

    smoothingSet = smoothingUpToDate = true;
}

void QControlFrame::setModelIndex(int index) {
    modelFunctionSelector->setCurrentIndex(index);
}

void QControlFrame::setOptimiserIndex(int index) {
    optimiserSelector->setCurrentIndex(index);
}

void QControlFrame::setSchemeIndex(int index) {
    fittingSelector->setCurrentIndex(index);
}

void QControlFrame::setFixedCBValue(QString value) {

    cbDensityButton->setChecked(true);

    cbDensityFixedButton->setChecked(true);
    cbDensityFixedEdit->setText(value);

    cbDensitySet = cBMDUpToDate = true;
}

void QControlFrame::setSigmaValues(QString x, QString y, QString z) {

    sigmaFixedButton->setChecked(true);
    sigmaXEdit->setText(x); sigmaYEdit->setText(y); sigmaZEdit->setText(z);
    sigmaXSet = sigmaYSet = sigmaZSet = sigmaUpToDate = true;
}

void QControlFrame::setThresholdSelection(int index) {
    removeThresholdValues();
    thresholdsSelector->setCurrentIndex(index);

}

void QControlFrame::setThresholdValues(QString stThresh, QString cbThresh, QString threshold, int thresholdIndex) { // todo read in and display the percentage

    if(thresholdIndex==itk::ClassifierTransform::kMedianManual && threshWeightSet) {
        thresholdsSelector->blockSignals(true);
        thresholdsSelector->setCurrentIndex(thresholdIndex);
        thresholdsSelector->blockSignals(false);
        stThresholdEdit->setText(stThresh); cbThresholdEdit->setText(cbThresh); thresholdEdit->setText(threshold);
        thresholdsUpToDate = threshSTSet = threshCBSet = thresholdSet = true;
    } else if(thresholdIndex==itk::ClassifierTransform::kMedianManual && !threshWeightSet) {
        cerr<<"Error in QControlFrame::setThresholdValues kMedianManual mode, but weight not set."<<endl;
        removeThresholdValues();
    } else if(thresholdIndex >= itk::LinearTransform::kVLowPercent && thresholdIndex <= itk::LinearTransform::kMedianManual) {
        thresholdsSelector->blockSignals(true);
        thresholdsSelector->setCurrentIndex(thresholdIndex);
        thresholdsSelector->blockSignals(false);
        stThresholdEdit->setText(stThresh); cbThresholdEdit->setText(cbThresh); thresholdEdit->setText(threshold);
        thresholdsUpToDate = threshSTSet = threshCBSet = thresholdSet = true;
    } else {
        cerr<<"Error in QControlFrame::setThresholdValues invalid threshold mode: "<<thresholdIndex<<endl;
        removeThresholdValues();
    }
}

void QControlFrame::setThresholdValues(QString stThresh, QString cbThresh, QString threshold, QString weight, int thresholdIndex) {

    if(thresholdIndex==itk::ClassifierTransform::kMedianManual) {
        thresholdsSelector->blockSignals(true);
        thresholdsSelector->setCurrentIndex(thresholdIndex);
        thresholdsSelector->blockSignals(false);
        stThresholdEdit->setText(stThresh); cbThresholdEdit->setText(cbThresh); thresholdEdit->setText(threshold);
        thresholdWeightingEdit->setText(weight);
        thresholdsUpToDate = threshSTSet = threshCBSet = thresholdSet = threshWeightSet = true;
    } else {
        cerr<<"Error in QControlFrame::setThresholdValues kMedianManual mode, but weight not set."<<endl;
        removeThresholdValues();
    }
}

void QControlFrame::setThresholdWeight(QString weight, int thresholdIndex) {

    removeThresholdValues();
    if(thresholdIndex==itk::LinearTransform::kMedianManual) {
        thresholdsSelector->blockSignals(true);
        thresholdsSelector->setCurrentIndex(thresholdIndex);
        thresholdsSelector->blockSignals(false);
        thresholdWeightingEdit->setText(weight);
        threshWeightSet = true;
    } else {
        cerr<<"Error in QControlFrame::setThresholdWeight - not kMedianManual mode"<<endl;
    }
}

void QControlFrame::removeThresholdValues() {
    stThresholdEdit->setText(QString()); cbThresholdEdit->setText(QString()); thresholdEdit->setText(QString());
    thresholdsUpToDate = false;
}

void QControlFrame::setFWHMCBMode() {

    cbDensityButton->setChecked(true);
    cbDensityFWHMButton->setChecked(true);

    cBMDUpToDate = true;
}

void QControlFrame::setAveragingValues(QString x, QString y, QString z) {
    profileAveragingButton->setChecked(true);
    xAveragingEdit->setText(x);
    yAveragingEdit->setText(y);
    zAveragingEdit->setText(z);
    xAveragingSet = yAveragingSet = zAveragingSet = averagingUpToDate = true;
}

void QControlFrame::setImportPath(QString importPath) {
    importPathEdit->setText(importPath);
    updateDisplaySelectorOptions();
}

void QControlFrame::setCalibrationScalingValues(double p0, double p1, double p2) {
    p0Edit->setText(QString::number(p0));
    p1Edit->setText(QString::number(p1));
    if(calibrationPhantomSelector->currentIndex()!=itk::CorticalBone::kManualLinearCal) {
        p2Edit->setText(QString::number(p2));
    } else {
        p2Edit->setText(QString::number(0.0));
    }
    p0CalSet = p1CalSet = p2CalSet = calScalesUpToDate = calUpToDate = true;

}

void QControlFrame::setCalibrationRadius(double radius) {
    calRadiusEdit->setText(QString::number(radius));
    calRadiusSet=true;
}

void QControlFrame::setCalibrationNumber(double number) {
    calNumberEdit->setText(QString::number(number));
    calNumberSet=true;
}

void QControlFrame::removeCalibrationScalingValues() {
    p0Edit->setText(QString(""));
    p1Edit->setText(QString(""));
    p2Edit->setText(QString(""));
    p0CalSet = p1CalSet = p2CalSet = calScalesUpToDate = false;

}

void QControlFrame::setState(int state) {

    measurementsUpToDate = (state == 4); // are measurements currently up to date

    if (state==0 || state==1) { // image not open; mesh don't care

        displayParameterSelector->setDisabled(true);
        displayMeshSelector->setDisabled(true);

        disableMeasurements();
        disableMeasurementSettings();
        disableCalibration();

        disableCBConstraints();
        disableSigmaConstraints();
        disableThresholds();
        disableSampleNumber();
        disableSmoothing();
        disableAveraging();

        disableImports();

    } else if (state==2) { // if the image open; mesh not open

        displayParameterSelector->setDisabled(true);

        displayMeshSelector->setEnabled(true);
        if(displayMeshSelector->currentIndex()==1) {
            displayMeshSelector->setCurrentIndex(0);
        }
        displayMeshSelector->setItemData(1, 0, Qt::UserRole - 1); // dirty hack - disable todo remove

        disableMeasurements();
        disableMeasurementSettings();

        enableCalibration();
        disableCBConstraints();
        disableSigmaConstraints();
        disableThresholds();
        disableSampleNumber();
        disableSmoothing();
        disableAveraging();

        disableImports();



    } else if (state==3) { // if the image; mesh open; calibration don't care; measurements not made

        displayParameterSelector->setDisabled(true);
        displayParameterSelector->setCurrentIndex(-1);

        updateMeasureEnables();
        enableMeasurementSettings();

        enableCalibration();
        enableCBConstraints();
        updateThresholds();
        enableSigmaConstraints();
        enableAveraging();
        enableSampleNumber();

        if(smoothingRequired) {
            enableSmoothing();
        } else {
            disableSmoothing();
        }

        displayMeshSelector->setEnabled(true);
        if(displayMeshSelector->currentIndex()==1) {
            displayMeshSelector->setCurrentIndex(0);
        }
        // disable periosteal mesh
        displayMeshSelector->setItemData(1, 0, Qt::UserRole - 1); // dirty hack - disable todo remove

        enableImports();

    } else if (state==4) { // if the image; mesh open; calibration don't care; measurements made

        displayParameterSelector->setEnabled(true);

        displayMeshSelector->setEnabled(true);
        // enable periosteal mesh
        displayMeshSelector->setItemData(1, 33, Qt::UserRole - 1); // dirty hack - enable todo remove

        updateMeasureEnables();
        enableMeasurementSettings();

        enableCalibration();
        enableCBConstraints();
        updateThresholds();
        enableSigmaConstraints();
        enableAveraging();
        enableSampleNumber();

        if(smoothingRequired) {
            enableSmoothing();
        } else {
            disableSmoothing();
        }


        enableImports();

    } else {
        std::cout<<"incorrect state selected"<<std::endl;
    }

}

//------------- Getters ------------//
int QControlFrame::getSampleNumber() {
    return sampleNumberEdit->text().toInt();
}

double QControlFrame::getSmoothingRadius() {
    return smoothingEdit->text().toDouble();
}

bool QControlFrame::getCBDensity(double &density) {
    if(cbDensityFixedButton->isChecked() && cbDensitySet) {
        density = cbDensityFixedEdit->text().toDouble();
        return true;
    } else {
        return false;
    }
}

bool QControlFrame::getSigmaValues(double &x, double &y, double &z) {

    if(sigmaXSet && sigmaYSet && sigmaZSet && sigmaFixedButton->isChecked()) {
        x = sigmaXEdit->text().toDouble();
        y = sigmaYEdit->text().toDouble();
        z = sigmaZEdit->text().toDouble();
        return true;
    } else {
        return false;
    }
}

bool QControlFrame::getCalibrationValues(double &p0, double &p1, double &p2) { // calibration not nessecarily run - but values selected if in manual mode

    if(calibrationPhantomSelector->currentIndex()==itk::CorticalBone::kManualControlPtsCal) {
        p0=-1; p1=-1; p2=-1;
    } else if(p0CalSet && p1CalSet && p2CalSet && calScalesUpToDate) {
        p0=p0Edit->text().toDouble();
        p1=p1Edit->text().toDouble();
        p2=p2Edit->text().toDouble();
    }
    return calScalesUpToDate;
}

bool QControlFrame::getCalibrationRadius(double& radius) {

    if(calibrationPhantomSelector->currentIndex()==itk::CorticalBone::kManualControlPtsCal && calRadiusSet) {
        radius=calRadiusEdit->text().toDouble();
    } else {
        radius = nan("1");
    }
    return calRadiusSet;
}

bool QControlFrame::getCalibrationNumber(int& number) {
    if(calibrationPhantomSelector->currentIndex()==itk::CorticalBone::kManualControlPtsCal && calNumberSet) {
        number=calNumberEdit->text().toInt();
    } else {
        number = -1;
    }
    return calNumberSet;
}

int QControlFrame::getThresholdSelection() {
    return thresholdsSelector->currentIndex();
}

QString QControlFrame::getClassifierThresholdName() {
    return thresholdsSelector->currentText();
}

bool QControlFrame::getThresholdValues(double &st, double &cb, double &threshold, int &thresholdIndex) {

    thresholdIndex = thresholdsSelector->currentIndex();

    if(threshSTSet && threshCBSet && thresholdSet) {
        st = stThresholdEdit->text().toDouble();
        cb = cbThresholdEdit->text().toDouble();
        threshold = thresholdEdit->text().toDouble();
    } else {
        st = nan("1"); cb = nan("1"); threshold = nan("1");
    }
    return thresholdsUpToDate;
}

bool QControlFrame::getThresholdWeight(double &weight) {
    if(thresholdsSelector->currentIndex()==itk::ClassifierTransform::kMedianManual && threshWeightSet) {
        weight = thresholdWeightingEdit->text().toDouble();
        return true;
    } else {
        weight = nan("1");
        return false;
    }
}

bool QControlFrame::getAveragingValues(double &x, double &y, double &z) {
    if(xAveragingSet && yAveragingSet && zAveragingSet && profileAveragingButton->isChecked()) {
        x = xAveragingEdit->text().toDouble();
        y = yAveragingEdit->text().toDouble();
        z = zAveragingEdit->text().toDouble();
        return true;
    } else {
        return false;
    }
}

bool QControlFrame::getImportDensityScalingValues(double &p0, double &p1, double &p2) {
    if(p0CalSet && p1CalSet && p2CalSet) {
        p0 = p0Edit->text().toDouble();
        p1 = p1Edit->text().toDouble();
        p2 = p2Edit->text().toDouble();
        return true;
    } else {
        return false;
    }
}

//----------- visibility ------------//
void QControlFrame::setDisplayBoxVisibilities(bool visibility) {
    if(visibility) {
        displayBox->show();
    } else {
        displayBox->hide();
    }
}

void QControlFrame::setCalibrationBoxVisibilities(bool visibility) {
    if(visibility) {
        calibrationBox->show();
    } else {
        calibrationBox->hide();
    }
}

void QControlFrame::setConstraintBoxVisibilities(bool visibility) {
    if(visibility) {
        constraintsBox->show();
    } else {
        constraintsBox->hide();
    }
}

void QControlFrame::setThresholdBoxVisibilities(bool visibility) {
    if(visibility) {
        thresholdsBox->show();
    } else {
        thresholdsBox->hide();
    }
}

void QControlFrame::setProfileSettingsBoxVisibilities(bool visibility) {
    if(visibility) {
        profileSettingsBox->show();
    } else {
        profileSettingsBox->hide();
    }
}

void QControlFrame::setImportBoxVisibilities(bool visibility) {
    if(visibility) {
        importBox->show();
    } else {
        importBox->hide();
    }
}


//-------------- state getters ---------------//
bool QControlFrame::isFWHMMode() {
    return cbDensityFWHMButton->isChecked();
}

bool QControlFrame::isFixedCBMode() {
    return cbDensityFixedButton->isChecked();
}

bool QControlFrame::isPtMeasuresEnabled() {
    bool upToDate = sampleNumberUpToDate && averagingUpToDate
                    && calScalesUpToDate && calUpToDate && smoothingUpToDate && (!measurementsRequired || measurementsUpToDate);

    upToDate = (fittingMode || calibrationMode) ? upToDate && cBMDUpToDate && sigmaUpToDate : upToDate &&
                                                                                              thresholdsUpToDate;

    return upToDate && ptMeasureButton->isChecked();
}

bool QControlFrame::isConstraintMode() {
    return modelMode <= itk::CorticalBone::kEndostealRamp;
}

bool QControlFrame::isCalibrationMode() {
    return modelMode == itk::CorticalBone::kCalibration;
}

bool QControlFrame::isManualCalibrationMode() {
    return calibrationPhantomSelector->currentIndex() == itk::CorticalBone::kManualLinearCal || calibrationPhantomSelector->currentIndex() == itk::CorticalBone::kManualQuadraticCal;
}

bool QControlFrame::isManualCalibrationPtsMode() {
    return calibrationPhantomSelector->currentIndex() == itk::CorticalBone::kManualControlPtsCal;
}

bool QControlFrame::isManualCalibrationPtsUpToDate() {
    return calibrationPhantomSelector->currentIndex() == itk::CorticalBone::kManualControlPtsCal && calRadiusSet && calNumberSet;
}

bool QControlFrame::isCalibrationReady() {

    if(calibrationPhantomSelector->currentIndex() == itk::CorticalBone::kManualLinearCal || calibrationPhantomSelector->currentIndex() == itk::CorticalBone::kManualQuadraticCal) {
        return calScalesUpToDate;
    } else if(calibrationPhantomSelector->currentIndex() == itk::CorticalBone::kManualControlPtsCal) {
        return calRadiusSet && calNumberSet && calibrationPointsButton->isChecked();
    } else {
        return calibrationPointsButton->isChecked();
    }
}

bool QControlFrame::isManualThresholdsMode() {
    return thresholdsSelector->currentIndex()==itk::LinearTransform::kManualThreshold;
}

bool QControlFrame::isMedianManualThresholdsMode() {
    return thresholdsSelector->currentIndex()==itk::LinearTransform::kMedianManual;
}

bool QControlFrame::isClassifierMode() {
    return modelMode == itk::CorticalBone::kHighResClassifier;
}

//bool QControlFrame::isThresholdingMode() {
//    return modelMode == itk::CorticalBone::kHighResThresholder;
//}

bool QControlFrame::isSampleNumberSet() {
    return sampleNumberUpToDate;
}

bool QControlFrame::isThresholdingSet() {
    return thresholdsUpToDate;
}

//----------- measurement slot(s) ------------//
void QControlFrame::modelFunctionChanged(int index) {

    if(index >= itk::CorticalBone::kThreeTierRect && index <= itk::CorticalBone::kEndostealRamp) { // if changing

        modelMode = index;
        fittingMode = true; classifierMode = false; calibrationMode = false;

    } else if(index == itk::CorticalBone::kHighResClassifier) {// && index <= itk::CorticalBone::kHighResOptimiser) { // if changing

        modelMode = index;
        fittingMode = false; classifierMode = true; calibrationMode = false;

    } else if(index == itk::CorticalBone::kCalibration) {

        modelMode = index;
        fittingMode = false; classifierMode = false; calibrationMode = true;
        parentObject->removeImportedProfile(); // remove any imported values
    } else {
        cerr<<"Error in QControlFrame::modelFunctionChanged() invalid model index entered "<<index<<" set to kThreeTierRect"<<endl;
        modelMode = itk::CorticalBone::kThreeTierRect;
    }

    updateFittingSelectorOptions();
    updateDisplaySelectorOptions();
}

void QControlFrame::fittingSchemeChanged(int index) {

    if(modelMode<=itk::CorticalBone::kEndostealRamp) {
        if(index >= itk::CorticalBone::kCBMV2AFitting && index <= itk::CorticalBone::kCBMV2DFitting) { // is set to three pass require sigma set to enable measurements
            sigmaRequired = true;
        } else {
            sigmaRequired = false;
        }
        averagingRequired = false;
        if(index >= itk::CorticalBone::kCBSmoothingAFitting && index <= itk::CorticalBone::kCBSmoothingDFitting) {
            measurementsRequired = smoothingRequired = true;
            toggleSmoothing(smoothingButton->isChecked());
        } else {
            measurementsRequired = smoothingRequired = false;
            toggleSmoothing(false); // is disabled therefore force to be enabled
        }
        toggleFixedSigma(sigmaFixedButton->isChecked());

    } else if(modelMode==itk::CorticalBone::kHighResClassifier) {

        measurementsRequired = smoothingRequired = false; averagingRequired = true;

    } else if(modelMode==itk::CorticalBone::kCalibration) {
        // do nothing
    }
//    else if(modelMode==itk::CorticalBone::kHighResThresholder) {
//
//        measurementsRequired = smoothingRequired = false;
//
//    }
}

//----------- Sample Number Slots ------------//
void QControlFrame::sampleNumberChanging() {

    sampleNumberUpToDate = sampleNumberSet = false;

}

void QControlFrame::sampleNumberChanged() {

    sampleNumberUpToDate = sampleNumberSet = true;

}

void QControlFrame::toggleSampleNumber(bool state) {

    if(state) { // is the sample number set?
        sampleNumberUpToDate = sampleNumberSet;
    } else {
        sampleNumberUpToDate = true;
    }
}

//---------- Smoothing Slots -------------//
void QControlFrame::toggleSmoothing(bool state) {

    if(state) { // is the sample number set?
        smoothingUpToDate = smoothingSet;
    } else {
        smoothingUpToDate = !smoothingRequired;
    }

}

void QControlFrame::smoothingChanging() {
    smoothingUpToDate = smoothingSet = false;
}

void QControlFrame::smoothingChanged() {
    smoothingUpToDate = smoothingSet = true;
}

//------------ Fix CB Slots -------------// 
void QControlFrame::fixedCBDensityChanging() {

    // disable make measure until done editing
    cbDensitySet = cBMDUpToDate = false;
}

void QControlFrame::fixedCBDensityChanged() {
    cbDensitySet = cBMDUpToDate = true;
}

void QControlFrame::turnOnFixedCBDensity() {

    cBMDUpToDate = cbDensitySet;

}

void QControlFrame::toggleCBConstraintDensity(bool state) {
    // disable - CB updated from a separate slot in MainWindow

    if(state) { // update 'uptodate' bool
        if(cbDensityFWHMButton->isChecked()) {
            cBMDUpToDate = true;
        } else {
            cBMDUpToDate = cbDensitySet;
        }
    } else { // disable FWHM and Fixed

        cBMDUpToDate = true;
    }
}

void QControlFrame::turnOnFWHMCBDensity() {
    cBMDUpToDate = true;
}

//------------ Free Sigma Slots -------------//
void QControlFrame::toggleFixedSigma(bool state) {
    if(state) {
        if(sigmaXSet && sigmaYSet && sigmaZSet) {
            sigmaUpToDate = true;
        } else {
            sigmaUpToDate = false;
        }
    } else {
        if(sigmaRequired) {
            sigmaUpToDate = false;
        } else {
            sigmaUpToDate = true;
        }
    }

}

void QControlFrame::fixedSigmaXChanged() {
    sigmaXSet = true;
    if(sigmaXSet && sigmaYSet && sigmaZSet) {
        sigmaUpToDate = true;
    }
}

void QControlFrame::fixedSigmaXChanging() {
    sigmaXSet = sigmaUpToDate = false;
}

void QControlFrame::fixedSigmaYChanged() {
    sigmaYSet = true;
    if(sigmaXSet && sigmaYSet && sigmaZSet) {
        sigmaUpToDate = true;
    }
}

void QControlFrame::fixedSigmaYChanging() {
    sigmaYSet = sigmaUpToDate = false;
}

void QControlFrame::fixedSigmaZChanged() {
    sigmaZSet = true;
    if(sigmaXSet && sigmaYSet && sigmaZSet) {
        sigmaUpToDate = true;
    }
}

void QControlFrame::fixedSigmaZChanging() {
    sigmaZSet = sigmaUpToDate = false;
}

//---------- Threshold Slots ----------------------//
void QControlFrame::stThresholdChanged() {
    threshSTSet = true;
    if(threshSTSet && threshCBSet && thresholdSet) {
        updateThresholds();
    }
}

void QControlFrame::cbThresholdChanged() {
    threshCBSet = true;
    if(threshSTSet && threshCBSet && thresholdSet) {
        updateThresholds();
    }
}

void QControlFrame::thresholdChanged() {
    thresholdSet = true;
    if(threshSTSet && threshCBSet && thresholdSet) {
        updateThresholds();
    }
}

void QControlFrame::thresholdWeightingChanged() {
    threshWeightSet = true; updateThresholds();
}

void QControlFrame::stThresholdChanging() {
    threshSTSet = thresholdsUpToDate = false;
}

void QControlFrame::cbThresholdChanging() {
    threshCBSet = thresholdsUpToDate = false;
}

void QControlFrame::thresholdChanging() {
    thresholdSet = thresholdsUpToDate = false;
}

void QControlFrame::thresholdWeightingChanging() {
    threshWeightSet = thresholdsUpToDate = false;
}


//---------- HR Averaging Slots ------------------//
void QControlFrame::xAveragingChanging() {

    xAveragingSet = averagingUpToDate = false;

    // disable make measure until done editing
    disableMeasurements();
}

void QControlFrame::yAveragingChanging() {

    yAveragingSet = averagingUpToDate = false;
}

void QControlFrame::zAveragingChanging() {

    zAveragingSet = averagingUpToDate = false;
}

void QControlFrame::toggleProfileAveraging(bool state) {

    // enable / disable - HR averaging from a separate slot in MainWindow
    if(state) {
        if(xAveragingSet && yAveragingSet && zAveragingSet) {
            averagingUpToDate = true;
        } else {
            averagingUpToDate = false;
        }
    } else {
        averagingUpToDate = !averagingRequired;
    }
}

void QControlFrame::setXAveragingValue() {
    xAveragingSet = true;
    if(xAveragingSet && yAveragingSet && zAveragingSet) {
        averagingUpToDate = true;
    }
}

void QControlFrame::setYAveragingValue() {
    yAveragingSet = true;
    if(xAveragingSet && yAveragingSet && zAveragingSet) {
        averagingUpToDate = true;
    }
}

void QControlFrame::setZAveragingValue() {
    zAveragingSet = true;
    if(xAveragingSet && yAveragingSet && zAveragingSet) {
        averagingUpToDate = true;
    }
}

//----------- Calibration Value Changed -----------------//
void QControlFrame::calibrationPhantomChanged(int) {

    if(calibrationPhantomSelector->currentIndex()==itk::CorticalBone::kNoCal) {
        calUpToDate=true;
    } else {
        calUpToDate=false;
    }
    if(calibrationPhantomSelector->currentIndex()==itk::CorticalBone::kManualLinearCal) {
        p2Edit->setText(QString::number(0.0)); p2CalSet=true;
        calScalesUpToDate = p0CalSet && p1CalSet && p2CalSet;
    } else if(calibrationPhantomSelector->currentIndex()==itk::CorticalBone::kManualControlPtsCal) {
    } else {
        calScalesUpToDate = p0CalSet && p1CalSet && p2CalSet;
    }
}

void QControlFrame::p0ValueChanged() {
    p0CalSet =true; calUpToDate = false;
    if(p0CalSet && p1CalSet && p2CalSet) {
        calScalesUpToDate = true;
    }
}

void QControlFrame::p1ValueChanged() {
    p1CalSet =true; calUpToDate = false;
    if(p0CalSet && p1CalSet && p2CalSet) {
        calScalesUpToDate = true;
    }
}

void QControlFrame::p2ValueChanged() {
    p2CalSet = true; calUpToDate = false;
    if(p0CalSet && p1CalSet && p2CalSet) {
        calScalesUpToDate = true;
    }
}

void QControlFrame::p0ValueChanging() {
    p0CalSet = calScalesUpToDate = false; calUpToDate = false;
}

void QControlFrame::p1ValueChanging() {
    p1CalSet = calScalesUpToDate = false; calUpToDate = false;
}

void QControlFrame::p2ValueChanging() {
    p2CalSet = calScalesUpToDate = false; calUpToDate = false;
}

void QControlFrame::calRadiusChanged() {
    calRadiusSet = true;
}

void QControlFrame::calRadiusChanging() {
    calRadiusSet = false;
}

void QControlFrame::calNumberChanged() {
    calNumberSet=true;
}

void QControlFrame::calNumberChanging() {
    calNumberSet=false;
}


//------------- Enables ---------------//
void QControlFrame::updateMeasureEnables() {

    bool upToDate = sampleNumberUpToDate && averagingUpToDate && calScalesUpToDate && calUpToDate && smoothingUpToDate;

    if(fittingMode) {
        upToDate = upToDate && cBMDUpToDate && sigmaUpToDate;
    } else if(classifierMode) {
        upToDate = upToDate && thresholdsUpToDate;
    } else if(calibrationMode) {
        // do nothing
    }

    if(upToDate) {
        enableMeasurements();
    } else {
        disableMeasurements();
    }

}

void QControlFrame::disableMeasurements() {
    runOverMeshButton->setDisabled(true);
    ptMeasureButton->setDisabled(true);
}

void QControlFrame::enableMeasurements() {
    runOverMeshButton->setEnabled(true);
    if(!measurementsRequired || measurementsUpToDate) {
        ptMeasureButton->setEnabled(true);
    } else {
        ptMeasureButton->setDisabled(true);
    }
}

void QControlFrame::disableMeasurementSettings() {
    modelFunctionSelector->setDisabled(true);
    optimiserSelector->setDisabled(true);
    fittingSelector->setDisabled(true);
}

void QControlFrame::enableMeasurementSettings() {
    modelFunctionSelector->setEnabled(true);
    if(fittingMode) {
        optimiserSelector->setEnabled(true);
        fittingSelector->setEnabled(true);
    } else if(classifierMode) {
        optimiserSelector->setDisabled(true);
        fittingSelector->setEnabled(true);
    } else if(calibrationMode) {
        optimiserSelector->setDisabled(true);
        fittingSelector->setEnabled(true);
    } else { // same as fitting mode
        cerr<<"invalid mode combinations in ControlFrame::enableMeasurementSettings()"<<endl;
        optimiserSelector->setEnabled(true);
        fittingSelector->setEnabled(true);
    }
}

void QControlFrame::disableImports() {
    setImportButton->setDisabled(true);
    removeImportButton->setDisabled(true);
    importPathEdit->setDisabled(true);
}

void QControlFrame::enableImports() {

    setImportButton->setEnabled(true);

    if(importPathEdit->text().isEmpty()) { // not set
        removeImportButton->setDisabled(true);
        importPathEdit->setDisabled(true);

    } else { // import path set
        removeImportButton->setEnabled(true);
        importPathEdit->setDisabled(true);

    }

}

void QControlFrame::disableCBConstraints() {
    cbDensityButton->setDisabled(true); cbDensityLabel->setDisabled(true);
    cbDensityFixedButton->setDisabled(true);
    cbDensityFWHMButton->setDisabled(true);
    cbDensityFixedEdit->setDisabled(true);
}

void QControlFrame::enableCBConstraints() {

    if(!fittingMode) {
        disableCBConstraints();
        return;
    }

    cbDensityButton->setEnabled(true); cbDensityLabel->setEnabled(true);
    if(cbDensityButton->isChecked()) {
        cbDensityFixedButton->setEnabled(true);
        cbDensityFWHMButton->setEnabled(true);
        if(cbDensityFixedButton->isChecked()) {
            cbDensityFixedEdit->setEnabled(true);
        } else {
            cbDensityFixedEdit->setDisabled(true);
        }
    } else {
        cbDensityFixedButton->setDisabled(true);
        cbDensityFWHMButton->setDisabled(true);
        cbDensityFixedEdit->setDisabled(true);
    }
}

void QControlFrame::disableSigmaConstraints() {
    sigmaFixedButton->setDisabled(true);
    sigmaXLabel->setDisabled(true);
    sigmaXEdit->setDisabled(true);
    sigmaYLabel->setDisabled(true);
    sigmaYEdit->setDisabled(true);
    sigmaZLabel->setDisabled(true);
    sigmaZEdit->setDisabled(true);
}

void QControlFrame::enableSigmaConstraints() {

    if(!fittingMode) {
        disableSigmaConstraints();
        return;
    }

    sigmaFixedButton->setEnabled(true);
    if(sigmaFixedButton->isChecked()) {
        sigmaXLabel->setEnabled(true);
        sigmaXEdit->setEnabled(true);
        sigmaYLabel->setEnabled(true);
        sigmaYEdit->setEnabled(true);
        sigmaZLabel->setEnabled(true);
        sigmaZEdit->setEnabled(true);
        sigmaZLabel->setEnabled(true);
        sigmaZEdit->setEnabled(true);
    } else {
        sigmaXLabel->setDisabled(true);
        sigmaXEdit->setDisabled(true);
        sigmaYLabel->setDisabled(true);
        sigmaYEdit->setDisabled(true);
        sigmaZLabel->setDisabled(true);
        sigmaZEdit->setDisabled(true);
        sigmaZLabel->setDisabled(true);
        sigmaZEdit->setDisabled(true);
    }

}

void QControlFrame::updateThresholds() {

    if(fittingMode || calibrationMode) {
        disableThresholds();
    } else {
        enableThresholds();
    }
}

void QControlFrame::disableThresholds() {
    thresholdsSelector->setDisabled(true);
    calculateThresholdsButton->setDisabled(true);
    stThresholdEdit->setDisabled(true); stThresholdEdit->setDisabled(true);
    thresholdCBLabel->setDisabled(true); cbThresholdEdit->setDisabled(true);
    thresholdLabel->setDisabled(true); thresholdEdit->setDisabled(true);
}

void QControlFrame::enableThresholds() {

    thresholdsSelector->setEnabled(true);

    if(thresholdsUpToDate) { // if thresholds already measured or set
        calculateThresholdsButton->setDisabled(true);
    } else if( thresholdsSelector->currentIndex()==itk::ClassifierTransform::kManualThreshold && threshSTSet && threshCBSet && thresholdSet) { // if manual threshold mode and manual values not set
        calculateThresholdsButton->setEnabled(true); // at least one threshold not set
    } else if (thresholdsSelector->currentIndex()==itk::ClassifierTransform::kManualThreshold) {
        calculateThresholdsButton->setEnabled(false);
    } else if (thresholdsSelector->currentIndex()==itk::ClassifierTransform::kMedianManual && threshWeightSet) {
        calculateThresholdsButton->setEnabled(true);
    } else if (thresholdsSelector->currentIndex()==itk::ClassifierTransform::kMedianManual) {
        calculateThresholdsButton->setEnabled(false); // weight not set
    } else {
        calculateThresholdsButton->setEnabled(true);
    }

    if(thresholdsSelector->currentIndex()!=itk::LinearTransform::kManualThreshold) {
        thresholdSTLabel->setDisabled(true); stThresholdEdit->setDisabled(true);
        thresholdCBLabel->setDisabled(true); cbThresholdEdit->setDisabled(true);
        thresholdLabel->setDisabled(true); thresholdEdit->setDisabled(true);
    } else {
        thresholdSTLabel->setEnabled(true); stThresholdEdit->setEnabled(true);
        thresholdCBLabel->setEnabled(true); cbThresholdEdit->setEnabled(true);
        thresholdLabel->setEnabled(true); thresholdEdit->setEnabled(true);
    }
    if(thresholdsSelector->currentIndex()!=itk::LinearTransform::kMedianManual) {
        thresholdWeightingLabel->setDisabled(true); thresholdWeightingEdit->setDisabled(true);
        thresholdWeightingLabel->setVisible(false); thresholdWeightingEdit->setVisible(false);
    } else {
        thresholdWeightingLabel->setEnabled(true); thresholdWeightingEdit->setEnabled(true);
        thresholdWeightingLabel->setVisible(true); thresholdWeightingEdit->setVisible(true);
    }
}

void QControlFrame::disableAveraging() {
    profileAveragingButton->setDisabled(true); profileAveragingLabel->setDisabled(true);
    xAveragingLabel->setDisabled(true); xAveragingEdit->setDisabled(true);
    yAveragingLabel->setDisabled(true); yAveragingEdit->setDisabled(true);
    zAveragingLabel->setDisabled(true); zAveragingEdit->setDisabled(true);
}

void QControlFrame::enableAveraging() {
    profileAveragingButton->setEnabled(true); profileAveragingLabel->setEnabled(true);
    if(profileAveragingButton->isChecked()) {
        xAveragingLabel->setEnabled(true); xAveragingEdit->setEnabled(true);
        yAveragingLabel->setEnabled(true); yAveragingEdit->setEnabled(true);
        zAveragingLabel->setEnabled(true); zAveragingEdit->setEnabled(true);
    } else {
        xAveragingLabel->setDisabled(true); xAveragingEdit->setDisabled(true);
        yAveragingLabel->setDisabled(true); yAveragingEdit->setDisabled(true);
        zAveragingLabel->setDisabled(true); zAveragingEdit->setDisabled(true);
    }
}

void QControlFrame::disableSampleNumber() {
    sampleNumberButton->setDisabled(true);
    sampleNumberEdit->setDisabled(true);
}

void QControlFrame::enableSampleNumber() {

    sampleNumberButton->setEnabled(true);
    if(sampleNumberButton->isChecked()) {
        sampleNumberEdit->setEnabled(true);
    } else {
        sampleNumberEdit->setDisabled(true);
    }

}

void QControlFrame::disableSmoothing() {
    smoothingButton->setDisabled(true);
    smoothingEdit->setDisabled(true);
}

void QControlFrame::enableSmoothing() {
    smoothingButton->setEnabled(true);
    if(smoothingButton->isChecked()) {
        smoothingEdit->setEnabled(true);
    } else {
        smoothingEdit->setDisabled(true);
    }
}

void QControlFrame::disableCalibration() {
    calibrateButton->setDisabled(true);
    calibrationPointsButton->setDisabled(true);
    calibrationPhantomSelector->setDisabled(true);

    p0Edit->setDisabled(true); p1Edit->setDisabled(true); p2Edit->setDisabled(true);
    p0Label->setDisabled(true); p1Label->setDisabled(true); p2Label->setDisabled(true);
}

void QControlFrame::enableCalibration() {
    calibrationPhantomSelector->setEnabled(true);

    if(calibrationPhantomSelector->currentIndex()==itk::CorticalBone::kManualControlPtsCal) {

        if(calRadiusSet && calNumberSet) {
            calibrateButton->setEnabled(true);
            calibrationPointsButton->setEnabled(true);
        } else {
            calibrateButton->setDisabled(true);
            calibrationPointsButton->setDisabled(true);
        }

        calRadiusLabel->setEnabled(true); calRadiusEdit->setEnabled(true);
        calRadiusLabel->setVisible(true); calRadiusEdit->setVisible(true);
        calNumberLabel->setEnabled(true); calNumberEdit->setEnabled(true);
        calNumberLabel->setVisible(true); calNumberEdit->setVisible(true);

        p0Edit->setDisabled(true); p1Edit->setDisabled(true); p2Edit->setDisabled(true);
        p0Edit->setHidden(true); p1Edit->setHidden(true); p2Edit->setHidden(true);
        p0Label->setDisabled(true); p1Label->setDisabled(true); p2Label->setDisabled(true);
        p0Label->setHidden(true); p1Label->setHidden(true); p2Label->setHidden(true);
    } else if(calibrationPhantomSelector->currentIndex()==itk::CorticalBone::kManualLinearCal) {

        if(p0CalSet && p1CalSet && p2CalSet && calScalesUpToDate && !calUpToDate) {
            calibrateButton->setEnabled(true);
        } else {
            calibrateButton->setDisabled(true);
        }

        calibrationPointsButton->setDisabled(true);
        p0Edit->setEnabled(true); p1Edit->setEnabled(true); p2Edit->setDisabled(true);
        p0Label->setEnabled(true); p1Label->setEnabled(true); p2Label->setDisabled(true);
        p0Edit->setVisible(true); p1Edit->setVisible(true); p2Edit->setVisible(true);
        p0Label->setVisible(true); p1Label->setVisible(true); p2Label->setVisible(true);

        calRadiusLabel->setDisabled(true); calRadiusEdit->setDisabled(true);
        calRadiusLabel->setHidden(true); calRadiusEdit->setHidden(true);
        calNumberLabel->setDisabled(true); calNumberEdit->setDisabled(true);
        calNumberLabel->setHidden(true); calNumberEdit->setHidden(true);

    } else if (calibrationPhantomSelector->currentIndex()==itk::CorticalBone::kManualQuadraticCal) {

        if(p0CalSet && p1CalSet && p2CalSet && calScalesUpToDate && !calUpToDate) {
            calibrateButton->setEnabled(true);
        } else {
            calibrateButton->setDisabled(true);
        }

        calibrationPointsButton->setDisabled(true);
        p0Edit->setEnabled(true); p1Edit->setEnabled(true); p2Edit->setEnabled(true);
        p0Label->setEnabled(true); p1Label->setEnabled(true); p2Label->setEnabled(true);
        p0Edit->setVisible(true); p1Edit->setVisible(true); p2Edit->setVisible(true);
        p0Label->setVisible(true); p1Label->setVisible(true); p2Label->setVisible(true);

        calRadiusLabel->setDisabled(true); calRadiusEdit->setDisabled(true);
        calRadiusLabel->setHidden(true); calRadiusEdit->setHidden(true);
        calNumberLabel->setDisabled(true); calNumberEdit->setDisabled(true);
        calNumberLabel->setHidden(true); calNumberEdit->setHidden(true);

    } else if(calibrationPhantomSelector->currentIndex()==itk::CorticalBone::kNoCal) {
        calibrateButton->setDisabled(true);
        calibrationPointsButton->setDisabled(true); calibrationPointsButton->setChecked(false);
        p0Edit->setDisabled(true); p1Edit->setDisabled(true); p2Edit->setDisabled(true);
        p0Label->setDisabled(true); p1Label->setDisabled(true); p2Label->setDisabled(true);
        p0Edit->setVisible(true); p1Edit->setVisible(true); p2Edit->setVisible(true);
        p0Label->setVisible(true); p1Label->setVisible(true); p2Label->setVisible(true);

        calRadiusLabel->setDisabled(true); calRadiusEdit->setDisabled(true);
        calRadiusLabel->setHidden(true); calRadiusEdit->setHidden(true);
        calNumberLabel->setDisabled(true); calNumberEdit->setDisabled(true);
        calNumberLabel->setHidden(true); calNumberEdit->setHidden(true);

    } else  {
        calibrateButton->setEnabled(true);
        calibrationPointsButton->setEnabled(true);
        p0Edit->setDisabled(true); p1Edit->setDisabled(true); p2Edit->setDisabled(true);
        p0Label->setDisabled(true); p1Label->setDisabled(true); p2Label->setDisabled(true);
        p0Edit->setVisible(true); p1Edit->setVisible(true); p2Edit->setVisible(true);
        p0Label->setVisible(true); p1Label->setVisible(true); p2Label->setVisible(true);

        calRadiusLabel->setDisabled(true); calRadiusEdit->setDisabled(true);
        calRadiusLabel->setHidden(true); calRadiusEdit->setHidden(true);
        calNumberLabel->setDisabled(true); calNumberEdit->setDisabled(true);
        calNumberLabel->setHidden(true); calNumberEdit->setHidden(true);

    }
}

  
