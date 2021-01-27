/* 
 * File:   profileThresholder.h
 * Author: rap58
 *
 * Created on 21 May 2015, 10:24
 */

#ifndef PROFILECLASSIFIER_H
#define	PROFILECLASSIFIER_H

// includes
#include <vtkSmartPointer.h>

#include <itkImage.h>
#include <itkProfileSpatialObjectLocal.h>
#include <itkLinearInterpolateImageFunction.h>
#include "ProfileProcessor.h"
//#include <itkBSplineInterpolateImageFunction.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkLevenbergMarquardtOptimizer.h>
#include <vtkTable.h>

#include <classifierTransform.h>

// defines

namespace itk {

    class ProfileClassifier : public ProfileProcessor {

    public:

        typedef ProfileClassifier Self;
        typedef ProfileProcessor                          Superclass;
        typedef SmartPointer< Self >                   Pointer;
        typedef SmartPointer< const Self >             ConstPointer;

        // standard class macros
        itkNewMacro(Self);
        itkTypeMacro(ProfileClassifier, ProfileProcessor);

        typedef double                                                  ScalarType;
        const static unsigned int                                       Dimension = 3;
        typedef itk::ProfileSpatialObject<Dimension>                    ProfileType;
        typedef itk::TimeStamp                                          TimeType;
        typedef vtkSmartPointer<vtkDoubleArray>                         ArrayType;
        typedef itk::Array<double>                                      itkArrayType;
        typedef itk::Array<int>                                         itkIntArrayType;
        typedef  itk::Array2D<ScalarType>                               itkArray2DType;
        typedef vtkSmartPointer<vtkTable>                               TableType;
        typedef itk::LinearTransform::Pointer                           TransformPointer;

        ProfileClassifier();

        void SetModelScheme(int index) ITK_OVERRIDE;
        int  GetModelScheme() ITK_OVERRIDE;
        std::string GetModelSchemeName();
        void ClearThresholds();
        bool Process() ITK_OVERRIDE;
        bool InitialiseDensities();

        void TurnOnProfileResampling(ScalarType maximum);
        void TurnOffProfileResampling();
        bool GetResamplingMode();


        bool IsValid() const ITK_OVERRIDE;
        bool IsUpToDate() const ITK_OVERRIDE;
        bool AreProfilesSufficientlyInsideImage() const ITK_OVERRIDE;
        bool AreThresholdsSet() const;
        bool IsNan() const;
        //bool AreThresholdsRequired() const;

        TableType GetImageDataTable() const ITK_OVERRIDE;
        TableType GetDisplayModelTable() const ITK_OVERRIDE;
        TableType GetDisplayModelTable(vtkSmartPointer<vtkDoubleArray> parametersIn, double offset) const;
        TableType GetProcessingModelTable() const ITK_OVERRIDE;

        itkArrayType GetDisplayValues() const ITK_OVERRIDE;
        unsigned int GetNumberOfDisplayValues(bool includeImportedOptions=false) const ITK_OVERRIDE;
        std::string GetDisplayName(int index) const ITK_OVERRIDE;
        std::string GetDisplayNameShort(int index) const ITK_OVERRIDE;
        itkArray2DType GetDisplayRanges(bool includeImportedOptions=false) const ITK_OVERRIDE;

        itkArrayType GetClassifications() const;
        itkArrayType GetPercentages() const;

        itkArrayType GetParameterValues() const ITK_OVERRIDE;
        unsigned int GetNumberOfParameterValues() const ITK_OVERRIDE;
        std::string GetParameterName(int index) const ITK_OVERRIDE;
        std::string GetParameterNameShort(int index) const ITK_OVERRIDE;

        // getters for display values
        int GetClassifierLevelIndex();
        void SetClassifierLevelIndex(int percentIndex);
        std::string GetClassifierLevelName();
        bool SetClassifierLevel(ScalarType stDensity, ScalarType cbDensity, ScalarType threshold);
        bool GetClassifierThresolds(ScalarType &stDensity, ScalarType &cbDensity, ScalarType &threshold);
        ScalarType GetErrorMean() const;

    private:

        void InitaliseProfileArrays();

        // GenerateData() throw ( ExceptionObject );
        //virtual void  GenerateData() ITK_OVERRIDE;

        bool ReadyToProcess(bool thresholdsRequired=true);

        ClassifierTransform::Pointer m_classifierTF;

        bool m_readyToProcess;
        bool m_resampleAboutPeriostealEdge;
        ScalarType m_MaximumOffset;

    };

}
#endif	/* PROFILECLASSIFIER_H */

