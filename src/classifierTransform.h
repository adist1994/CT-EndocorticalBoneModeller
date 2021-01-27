/* 
 * File:   thresholderTransform.h
 * Author: rap58
 *
 * Created on 17 June 2015, 10:43
 */

#ifndef CLASSIFIERTRANSFORM_H
#define	CLASSIFIERTRANSFORM_H

// includes
#include <vtkSmartPointer.h>

#include <itkImage.h>
#include <itkProfileSpatialObjectLocal.h>
#include <itkLinearInterpolateImageFunction.h> 
//#include <itkBSplineInterpolateImageFunction.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <vtkTable.h>
#include <itkOptimizerParameters.h>
#include <vtkIntArray.h>

#include <linearTransform.h>

// defines

namespace itk {

    class ClassifierTransform: public LinearTransform {


    public:

        typedef ClassifierTransform Self;
        typedef LinearTransform                                         Superclass;
        typedef SmartPointer< Self >                                    Pointer;
        typedef SmartPointer< const Self >                              ConstPointer;

        // standard class macros
        itkNewMacro(Self);
        itkTypeMacro(ClassifierTransform, LinearTransform);

        typedef signed short                                            PixelType;
        typedef double                                                  ScalarType;
        const static unsigned int                                       Dimension = 3;
        typedef itk::TimeStamp                                          TimeType;
        typedef vtkSmartPointer<vtkDoubleArray>                         ArrayType;
        typedef itk::Array<double>                                      itkArrayType;
        typedef itk::Array<int>                                      itkIntArrayType;
        typedef vtkSmartPointer<vtkIntArray>                            IntArrayType;
        typedef vtkSmartPointer<vtkTable>                               TableType;
        typedef itk::OptimizerParameters< ScalarType >                  ParametersType;
        typedef  itk::Array2D<ScalarType> itkArray2DType;

        typedef enum {
            kInvalid=-1,
            kGlobal=0,
            kMedian,
        } Classifier;

        typedef enum { // classification indecies. Not altering the order of these will effect the median filter
            kST=0,
            kTB,
            kCB,
            kCTB, // intracortical-trabecular bone
            kCMC, // intracortical Medullary Canal
            kMC,  // Medullary
        } TissueTypes;

        typedef enum {
            kNoCB,        // no cortical bone
            kOneCBNoTB,   // one cortical layer, no trabecular bone
            kOneCBSomeTB,     // one cortical layer, trabecular present
            kMultiCBNoTB, // multiple cortical layers, no trabecular bone
            kMultiCBSomeTB,   // multiple cortical layers, trabecular present
        } BoneTypes;

        typedef enum { // display
            kNotBoneThreshold =0,
            kCorticalThresold,
            kThreshold, // surface thickness for profileProperties, MSA density for display
            kClassification,
            kNotBoneDensity,
            kCorticalDensity,
            kDenseDensity,
            kTotalCorticalDensity,
            kTrabecularDensity, // intra-cortical trabecular bone
            kCorticalWidth,
            kEndocorticalWidth,
            kEndocorticalSlope, // total cortical bone ratio
            kPorosity,
            kProfileNumber,
        } displayOptions;
        typedef enum { // param
            kSTMean=0,
            kCBMean,
            kDBMean,
            kECMean,
            kTotMean,
            kTBMean,
            kSTSlope,
            kCBSlope,
            kDBSlope,
            kECSlope,
            kTotSlope,
            kTBSlope,
            kSTIntercept,
            kCBIntercept,
            kDBIntercept,
            kECIntercept,
            kTotIntercept,
            kTBIntercept,
            kSTSTD,
            kCBSTD,
            kDBSTD,
            kECSTD,
            kTotSTD,
            kTBSTD,
            kSTStart,
            kCBStart,
            kDBStart,
            kECStart,
            kTotStart,
            kTBStart,
            kSTEnd,
            kCBEnd,
            kDBEnd,
            kECEnd,
            kTotEnd,
            kTBEnd,
            kPorosityParameter,
        } parameterOptions;

        ClassifierTransform();

        bool Threshold() ITK_OVERRIDE;

        void SetProfileRangeLimits(ScalarType maximum);

        // return mode name
        std::string  GetModeName() ITK_OVERRIDE;

        ProfilePropertiesType* GetProfileProperties();

        itkArrayType GetDisplayValues() const ITK_OVERRIDE;
        std::string GetDisplayName(int index) const ITK_OVERRIDE;
        std::string GetDisplayNameShort(int index) const ITK_OVERRIDE;
        itkArray2DType GetDisplayRanges(bool importSet=false) const ITK_OVERRIDE;

        itkArrayType GetClassifications() const;
        itkArrayType GetPercentages() const;

        itkArrayType GetParameterValues() const ITK_OVERRIDE;
        std::string GetParameterName(int index) const ITK_OVERRIDE;
        std::string GetParameterNameShort(int index) const ITK_OVERRIDE;

        // getters for display values
        TableType GetProfileDisplayTable() const ITK_OVERRIDE;
        TableType GetModelDisplayTable() const ITK_OVERRIDE;
        TableType GetModelDisplayTable(ArrayType parametersIn, double offset) const ITK_OVERRIDE;
        TableType GetClassifierTable() const ITK_OVERRIDE;


    private:
        bool CalculatePercentages();
        bool CalculateTissueSectionsDirect(); // find he 1st pore right of the 1/2 location between xP and the 1st xE estimate
        bool CalculateTissueSectionsStrict(); // find he 1st pore right of the max % - stricter
        bool PercentageBasedClassifing();

        void CalculateTissueRegressionLines() ITK_OVERRIDE;

        void MedianFilter(itkIntArrayType classArray, itkArrayType valueArray, int filterSize);

        void CalculateValidity() ITK_OVERRIDE;
        ScalarType CalculateAbsMeanError() const ITK_OVERRIDE;

        void ClearResults() ITK_OVERRIDE;

        itkArrayType m_tissueValues, m_filteredTissueValues, m_percentArray;
        itkIntArrayType m_tissueTypes, m_filteredTissueTypes;

        int m_classificationType;

        RegressionLineType m_DBLine, m_totLine; // cb = cortical bone, db = dense bone (i.e. without endocortical portion of cb)

        ScalarType m_stIndices[2], m_cbIndices[2], m_tbIndices[2], m_ecIndices[2], m_dbIndices[2], m_totIndices[2];

        ScalarType m_porosity;

        ScalarType m_profileRange;

    };
}
#endif	/* CLASSIFIERTRANSFORM_H */

