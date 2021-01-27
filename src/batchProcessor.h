//
// Created by rap58 on 31/08/15.
//

#ifndef MYPROJECT_BATCHPROCESSOR_H
#define MYPROJECT_BATCHPROCESSOR_H

#include <iostream>
#include <xmlFileAccess.h>
#include <QStringList>

//class QStringList;
class xmlFileAccess;

class BatchProcessor {

public:
    BatchProcessor(std::string scriptFileName, int processCode=kRunProcess, int startIndex=-1, int endIndex=-1);

    enum{ // view not supported in commandline mode
        kRunCreate=0, // just creates files
        kRunProcess, // (re-)creates files and processes them
    };

private:

    enum{
        kCreate=0,
        kView,
        kProcess,
    };

    enum {
        kDicomExtn=0,
        kRawExtn,
    };

    enum {
        kObjExtn=0,
        kVrmlExtn,
        kStlExtn,
    };

    void initialiseGeneralValues(int startIndexIn, int endIndexIn);
    void initialiseScriptValues();
    void initialiseProjectValues();

    bool openScript(std::string fileName);
    bool generateBatchFileLists(int code);
    QStringList generateDirectoryList();
    QStringList generateFileList(QString filter);


    bool runCreationScript();
    bool runProcessScript();
    bool openProject(std::string fileName);
    bool saveProject(std::string fileName);


    //--- general values ---//
    bool indiciesSet;
    int startIndex, endIndex;
    xmlFileAccess* xmlFileAccessor;
    itk::CorticalBone::Pointer corticalBone;


    //--- script specific members ---//

    // lists
    QStringList scriptProjectNameList;
    QStringList scriptImageNameList, scriptMeshNameList, scriptFixedCBFileList, scriptImportedFileList;
    QStringList scriptResultsNameList;

    // strings
    QString scriptProjectPath, scriptProjectFilter;
    QString scriptImageFilter, scriptMeshFilter;
    QString scriptImageExtn, scriptMeshExtn, scriptProjectExtn, scriptTxtExtn;
    QString scriptFixedCBFilter;
    QString scriptImportFile;
    QString scriptModelName, scriptSchemeName, scriptOptimiserName;
    QString calibrationPhantomName;

    // ints
    int scriptImageExtnIndex, scriptMeshExtnIndex;
    int scriptCBConstraintIndex, scriptSampleNumber;
    int scriptModelIndex, scriptOptimiserIndex, scriptClassifierIndex, scriptSchemeIndex;
    int scriptCalibrationIndex;

    // bools
    bool scriptFixedCBFilterSet;
    bool scriptSampleNumberSet;
    bool scriptAveragingSet;
    bool scriptImportSet, calibrationSet;

    // scalars
    double scriptProfileAverages[3], scriptDensityScales[3];

    static const int MAX_NUMBER_OF_FILES_TO_PROCESS = 1000;

    //--- project specific members ---//
    bool classifierThresholdsCalculated, classifierThresholdLevelSet;
    std::string projectImageName, projectMeshName;
    std::string projectImportPath;

};


#endif //MYPROJECT_BATCHPROCESSOR_H
