/* 
* File:   scriptframe.h
* Author: rap58
*
* Created on 01 October 2014, 10:16
*/

#ifndef SCRIPTFRAME_H
#define	SCRIPTFRAME_H

#include <QFrame>

class QPlainTextEdit;
class QSpinBox;
class QListWidget;
class QFile;
class QString;
class QGroupBox;
class QStringList;
class xmlFileAccess;
class MainWindow;

class QScriptFrame : public QFrame
{
  Q_OBJECT

  public:
    QScriptFrame(MainWindow *parent, QPlainTextEdit *scriptEditorIn, xmlFileAccess *xmlFileWriterIn);

    QString getProjectPath();
    void setProjectPath(QString filePath);

    bool isViewRunning();
    bool isProcessRunning();

  public slots:
    bool loadScript(QString &fileName);
    void nextProjectFile();

  private slots:
    // overall slots
    void setDirectory();
    void saveScript();
    void loadScript();
    void runScript();
    void projectFilterChanged();
    void projectFilterChanging();

    // tab slots
    void tabChanged();

    // create tab slots
    void imageFilterChanged();
    void imageFilterChanging();
    void meshFilterChanged();
    void meshFilterChanging();
    void imageExtnChanged(int index);
    void meshExtnChanged(int index);
    void sampleNumberChanging();
    void sampleNumberChanged();
    // cb constraints
    void cbValueFilterChanged();
    void cbValueFilterChanging();
    void toggleFixedCB(bool state);
    void toggleCBConstraintState(bool state);
    void toggleFWHMCB(bool state);
    void toggleSampleNumber(bool state);
    // sigma constraints
    void toggleFixedSigma(bool state);
    void fixedSigmaXChanged();
    void fixedSigmaXChanging();
    void fixedSigmaYChanged();
    void fixedSigmaYChanging();
    void fixedSigmaZChanged();
    void fixedSigmaZChanging();
    // averaging
    void xAveragingChanging();
    void yAveragingChanging();
    void zAveragingChanging();
    void toggleAveraging(bool state);
    void setXAveragingValue();
    void setYAveragingValue();
    void setZAveragingValue();
    void setImportFile();
    void p0ScaleChanging();
    void p1ScaleChanging();
    void p2ScaleChanging();
    void pFileTextChanging();
    void calibrationTypeChanged(int index);
    void setp0ScaleValue();
    void setp1ScaleValue();
    void setp2ScaleValue();
    void setpFileTextChanged();

    // process tab slots
    void modelSelectionChanged(int index);
    void processMinRangeChanged(int value);
    void processMaxRangeChanged(int value);

    void schemeChanged(int index);

    void thresholdLevelChanged(int index);
    void stThresholdChanged();
    void cbThresholdChanged();
    void thresholdChanged();
    void stThresholdChanging();
    void cbThresholdChanging();
    void thresholdChanging();
    void thresholdWeightChanged();
    void thresholdWeightChanging();

    void toggleSmoothing(bool state);
    void smoothingChanging();
    void smoothingChanged();

  private:

    enum {
        kCreateTab=0,
        kViewTab,
        kProcessTab,
    } tabIndicies;
    enum {
        kDicomExtn=0,
        kRawExtn,
        kTiffExtn,
        kQCTExtn,
    } imageExtns;
    enum {
        kObjExtn=0,
        kVrmlExtn,
        kStlExtn,
        kPlyExtn,
    } meshExtns;

    // methods
    bool reset();

    void createObjects();
    void createConnections();

    bool createProject(int index);
    bool viewProject(int index);
    bool processProject(int index);

    void updateObjectEnables();
    void updateEditorEnables();
    bool isCreateTabUpToDate();
    bool isViewTabUpToDate();
    bool isProcessTabUpToDate();

    void enableTabs();
    void disableTabs();

    void updateScript(bool resetRanges=true);
    void updateDisplay();
    void updateFileLists();
    QString getImageFile(QString dirName);

    void resetRangeLimits();
    void updateSchemeOptions(int index);

    QString getActiveTabName();

    QStringList generateFileList(QString filter);
    QStringList generateDirectoryList();

    MainWindow *parentObject;

    //--- GUI objects

    //- script text display
    QPlainTextEdit *scriptEditor;

    //- script controls
    QTabWidget *scriptControlsTab;

    // overall
    QPushButton *scriptDirectoryButton;
    QLineEdit *scriptDirectoryEdit;
    QLabel *projectFilterLabel;
    QLineEdit *projectFilterEdit;

    QPushButton *loadButton, *saveButton, *runButton, *nextButton;

    // create
    QGroupBox *createProjectsBox;
    QLabel *imageFilterLabel, *meshFilterLabel;
    QLineEdit *imageFilterEdit, *meshFilterEdit;
    QLabel *imageExtnLabel, *meshExtnLabel;
    QComboBox *imageExtnSelector, *meshExtnSelector;

    QCheckBox *sampleNumberButton; QLineEdit *sampleNumberEdit;

    // constraints
    // cb density
    QCheckBox *cbDensityButton; QLabel *cbDensityLabel;
    QRadioButton *cbDensityFixedButton, *cbDensityFWHMButton;
    QLineEdit *cbDensityFixedEdit;
    // sigma
    QCheckBox *sigmaFixedButton;
    QLabel *sigmaXLabel, *sigmaYLabel, *sigmaZLabel;
    QLineEdit *sigmaXEdit, *sigmaYEdit, *sigmaZEdit;
    bool sigmaXSet, sigmaYSet, sigmaZSet;

    // averaging
    QCheckBox *profileAveragingButton;
    QLineEdit *xAvgEdit, *yAvgEdit, *zAvgEdit;
    QLabel *xAvgLabel, *yAvgLabel, *zAvgLabel;

    QPushButton *ImportFileButton; QLineEdit *importFileEdit;

    // view
    QGroupBox *viewProjectsBox;

    QComboBox *calibrationTypeBox;
    QLineEdit *p0ScaleEdit, *p1ScaleEdit, *p2ScaleEdit;
    QLabel *p0ScaleLabel, *p1ScaleLabel, *p2ScaleLabel;
    QLabel *pFileLabel; QLineEdit *pFileEdit;
    bool calibrationIsInFile;

    // process
    QGroupBox *processProjectsBox;
    QLabel *modelTypeLabel, *optimiserTypeLabel, *classifierLevelLabel, *schemeTypeLabel, *rangeLabel;
    QComboBox *modelTypeSelector, *optimiserTypeSelector, *schemeTypeSelector;

    QComboBox *thresholdLevelSelector;
    QLineEdit *stThreshEdit, *cbThreshEdit, *thresholdEdit;
    QLabel *stThreshLabel, *cbThreshLabel, *thresholdLabel;
    QLabel *thresholdWeightingLabel; QLineEdit *thresholdWeightEdit;
    QGroupBox *thresholdBox, *thresholdWeightBox;

    QCheckBox *smoothingButton; QLineEdit *smoothingEdit;
    QSpinBox *startSpinBox, *stopSpinBox;

    //--- script objects
    xmlFileAccess *xmlFileAccessor;

    QString projectFilterString;
    QString imageFilterString, meshFilterString;
    QString imageExtnString, meshExtnString;
    QString projectExtnString, txtExtnString;
    QString fixedCBFilterString;
    QString importFileString;

    int sampleNumberInt; double smoothingValue; double profileAverages[3]; double calibrationScales[3]; double thresholds[3]; double thresholdWeight;
    double fixedCBValue;

    // list of project file names
    QStringList projectNameList;
    QStringList imageNameList, meshNameList, fixedCBFileList, importedFileList, calibrationFileList;
    QStringList resultsNameList;

    QString projectPath;

    int viewIteratorIndex;

    static const int MAX_NUMBER_OF_FILES_TO_PROCESS = 10000;



    //--- state
    bool projectDirectorySet, tabsUpToDate;
    bool runningViewScript, runningProcessScript;
    bool imageFilterSet, meshFilterSet, fixedCBFilterSet;
    bool sampleNumberSet, xAvgSet, yAvgSet, zAvgSet;
    bool importFileSet;
    bool p0ScaleSet, p1ScaleSet, p2ScaleSet, pFileSet;
    bool stThresholdSet, cbThresholdSet, thresholdSet;
    bool thresholdWeightSet;
    bool projectFilterSet;
    bool smoothingSet;
    bool fixedCBIsNumber, fixedCBSet, cbConstraintUpToDate, sigmaUpToDate; // as opposed to the file set
    bool sigmaRequired;

    bool includeSampleNumber, includeAveraging, includeCalibration, includeSmoothing, includeThresholds, includeThresholdWeight; int cbConstraintIndex;


};

#endif	/* SCRIPTFRAME_H */

