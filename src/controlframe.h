/* 
 * File:   controlFrame.h
 * Author: rap58
 *
 * Created on 02 October 2014, 10:11
 */

#ifndef CONTROLFRAME_H
#define	CONTROLFRAME_H

#include <QScrollArea>
#include "mainwindow.h"

class QPushButton;
class QCheckBox;
class QRadioButton;
class QComboBox;
class QLineEdit;
class QLabel;
class QToolButton;
class QGroupBox;

class QControlFrame : public QScrollArea
{
Q_OBJECT

public:
    QControlFrame(MainWindow *parent = 0);
    void initialiseSettings();

    // setters
    void setState(int state);
    void setCalibrationPhantom(int index);
    void setCalibrationScalingValues(double p0, double p1, double p2);
    void setCalibrationRadius(double radius);
    void setCalibrationNumber(double number);
    void removeCalibrationScalingValues();
    void setFixedSampleNumber(QString value);
    void setFixedCBValue(QString value);
    void setSigmaValues(QString x, QString y, QString z);
    void setThresholdValues(QString st, QString cb, QString threshold, QString weight, int thresholdIndex);
    void setThresholdValues(QString st, QString cb, QString threshold, int thresholdIndex);
    void setThresholdWeight(QString weight, int thresholdIndex);
    void setThresholdSelection(int index); // used by script to adjust the percent index
    void removeThresholdValues();
    void setFWHMCBMode();
    void setAveragingValues(QString x, QString y, QString z);
    void setImportPath(QString);
    void setSmoothingRadius(QString value);
    void setModelIndex(int index);
    void setOptimiserIndex(int index);
    void setSchemeIndex(int index);

    // getters
    int getSampleNumber();
    double getSmoothingRadius();
    bool getCBDensity(double &density);
    bool getSigmaValues(double &x, double &y, double &z);
    bool getCalibrationValues(double &p0, double &p1, double &p2);
    bool getCalibrationRadius(double& radius);
    bool getCalibrationNumber(int& number);
    int getThresholdSelection();
    QString getClassifierThresholdName();
    bool getThresholdValues(double &st, double &cb, double &threshold, int &thresholdIndex);
    bool getThresholdWeight(double &weight);
    bool getAveragingValues(double &x, double &y, double &z);
    bool getImportDensityScalingValues(double &p0, double &p1, double &p2);

    // visibility
    void setDisplayBoxVisibilities(bool visibility);
    void setCalibrationBoxVisibilities(bool visibility);
    void setConstraintBoxVisibilities(bool visibility);
    void setThresholdBoxVisibilities(bool visibility);
    void setProfileSettingsBoxVisibilities(bool visibility);
    void setImportBoxVisibilities(bool visibility);


    // get state
    bool isFWHMMode();
    bool isFixedCBMode();
    bool isPtMeasuresEnabled();
    bool isConstraintMode();
    bool isClassifierMode();
    bool isCalibrationMode();
    bool isManualCalibrationMode();
    bool isManualCalibrationPtsMode();
    bool isManualCalibrationPtsUpToDate();
    bool isCalibrationReady();
    bool isManualThresholdsMode();
    bool isMedianManualThresholdsMode();
    bool isSampleNumberSet();
    bool isThresholdingSet();

private slots:
    // measurements
    void modelFunctionChanged(int index);
    void fittingSchemeChanged(int index);
    // sample number
    void sampleNumberChanging();
    void sampleNumberChanged();
    void toggleSampleNumber(bool state);
    // smoothing
    void toggleSmoothing(bool state);
    void smoothingChanging();
    void smoothingChanged();
    // fixed CB
    void fixedCBDensityChanging();
    void fixedCBDensityChanged();
    void turnOnFixedCBDensity();
    void toggleCBConstraintDensity(bool state);
    void turnOnFWHMCBDensity();
    // fixed sigma
    void toggleFixedSigma(bool state);
    void fixedSigmaXChanged();
    void fixedSigmaXChanging();
    void fixedSigmaYChanged();
    void fixedSigmaYChanging();
    void fixedSigmaZChanged();
    void fixedSigmaZChanging();
    // set thresholds
    void stThresholdChanged();
    void cbThresholdChanged();
    void thresholdChanged();
    void stThresholdChanging();
    void cbThresholdChanging();
    void thresholdChanging();
    void thresholdWeightingChanging();
    void thresholdWeightingChanged();
    //void thresholdSelectorChanged(int index);
    // HR averaging slots
    void xAveragingChanging();
    void yAveragingChanging();
    void zAveragingChanging();
    void toggleProfileAveraging(bool state);
    void setXAveragingValue();
    void setYAveragingValue();
    void setZAveragingValue();
    // calibration slots
    void calibrationPhantomChanged(int);
    void p0ValueChanged();
    void p1ValueChanged();
    void p2ValueChanged();
    void p0ValueChanging();
    void p1ValueChanging();
    void p2ValueChanging();
    void calRadiusChanged();
    void calRadiusChanging();
    void calNumberChanged();
    void calNumberChanging();

private:

    MainWindow *parentObject;

    // measurements
    QCheckBox *ptMeasureButton;
    QPushButton *runOverMeshButton;
    QComboBox *modelFunctionSelector; bool fittingMode, classifierMode, calibrationMode; int modelMode;
    QComboBox *optimiserSelector;
    QComboBox *fittingSelector; bool sigmaRequired, measurementsRequired, smoothingRequired, averagingRequired;

    // calibration
    QCheckBox *calibrationPointsButton;
    QPushButton *calibrateButton;
    QComboBox *calibrationPhantomSelector;
    QGroupBox *calibrationBox;
    QLineEdit *p0Edit, *p1Edit, *p2Edit;
    QLabel *p0Label, *p1Label, *p2Label;
    bool p0CalSet, p1CalSet, p2CalSet;
    QLineEdit *calRadiusEdit, *calNumberEdit;
    QLabel *calRadiusLabel, *calNumberLabel;
    bool calRadiusSet, calNumberSet;

    // display
    QComboBox *displayParameterSelector;
    QComboBox *displayMeshSelector;
    QGroupBox *displayBox;

    // constraints
    QCheckBox *cbDensityButton; QLabel *cbDensityLabel;
    QRadioButton *cbDensityFWHMButton;
    QRadioButton *cbDensityFixedButton;
    QLineEdit *cbDensityFixedEdit;
    bool cbDensitySet;
    QCheckBox *sigmaFixedButton;
    QLabel *sigmaXLabel, *sigmaYLabel, *sigmaZLabel;
    QLineEdit *sigmaXEdit, *sigmaYEdit, *sigmaZEdit;
    bool sigmaXSet, sigmaYSet, sigmaZSet;
    QGroupBox *constraintsBox;

    // thresholds
    QComboBox *thresholdsSelector; QPushButton *calculateThresholdsButton;
    QLabel *thresholdSTLabel, *thresholdCBLabel, *thresholdLabel;
    QLineEdit *stThresholdEdit, *cbThresholdEdit, *thresholdEdit;
    QLabel *thresholdWeightingLabel; QLineEdit * thresholdWeightingEdit;
    bool threshSTSet, threshCBSet, thresholdSet, threshWeightSet;
    QGroupBox *thresholdsBox;

    // profile settings
    QCheckBox *sampleNumberButton;
    QLineEdit *sampleNumberEdit; bool sampleNumberSet;
    QCheckBox *profileAveragingButton; QLabel *profileAveragingLabel;
    QLabel *xAveragingLabel, *yAveragingLabel, *zAveragingLabel;
    QLineEdit *xAveragingEdit,*yAveragingEdit, *zAveragingEdit;
    bool xAveragingSet, yAveragingSet, zAveragingSet;
    QCheckBox *smoothingButton;
    QLineEdit *smoothingEdit; bool smoothingSet;
    QGroupBox *profileSettingsBox;


    // import
    QPushButton *setImportButton;
    QToolButton *removeImportButton;
    QLineEdit *importPathEdit;
    QGroupBox *importBox;

    bool sampleNumberUpToDate, cBMDUpToDate, sigmaUpToDate, averagingUpToDate, calUpToDate,
            calScalesUpToDate, smoothingUpToDate, measurementsUpToDate, thresholdsUpToDate;


    void createObjects();
    void createConnections();

    void updateFittingSelectorOptions();
    void updateDisplaySelectorOptions();

    void updateMeasureEnables();
    void disableMeasurements();
    void enableMeasurements();

    void disableMeasurementSettings();
    void enableMeasurementSettings();

    void disableImports();
    void enableImports();

    void disableCBConstraints();
    void enableCBConstraints();

    void disableSigmaConstraints();
    void enableSigmaConstraints();

    void updateThresholds();
    void disableThresholds();
    void enableThresholds();

    void disableAveraging();
    void enableAveraging();

    void disableSampleNumber();
    void enableSampleNumber();

    void disableSmoothing();
    void enableSmoothing();

    void disableCalibration();
    void enableCalibration();



    // display results - model parameters


};


#endif	/* CONTROLFRAME_H */

