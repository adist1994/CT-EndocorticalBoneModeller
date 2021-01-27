//
// Created by rap58 on 20/10/15.
//

#ifndef MYPROJECT_PROFILECALIBRATION_H
#define MYPROJECT_PROFILECALIBRATION_H


#include "ProfileProcessor.h"
#include <itkSmartPointer.h>

namespace itk {

    class ProfileCalibration : public ProfileProcessor {

    public:

        typedef ProfileCalibration Self;
        typedef ProfileProcessor                          Superclass;
        typedef SmartPointer< Self >                   Pointer;
        typedef SmartPointer< const Self >             ConstPointer;

        // standard class macros
        itkNewMacro(Self);
        itkTypeMacro(ProfileThresholder, ProfileProcessor);

        typedef double                                                  ScalarType;
        const static unsigned int                                       Dimension = 3;
        typedef itk::ProfileSpatialObject<Dimension>                    ProfileType;
        typedef itk::TimeStamp                                          TimeType;
        typedef vtkSmartPointer<vtkDoubleArray>                         ArrayType;
        typedef itk::Array<double>                                      itkArrayType;
        typedef vtkSmartPointer<vtkIntArray>                            IntArrayType;
        typedef  itk::Array2D<ScalarType> itkArray2DType;
        typedef vtkSmartPointer<vtkTable>                               TableType;
        typedef itk::LinearTransform::Pointer                            TransformPointer;
        //typedef itk::BSplineInterpolateImageFunction< ImageType, float, float> InterpolatorType; // also window thresholding. read up in ITKManual


        typedef enum { // display & parameters
            kLowMedian=0,
            kHighMedian,
        } displayOptions;

        typedef enum { // display & parameters
            kMedian=0,
            kMaximum,
        } operationMode;

        ProfileCalibration();

        bool Process() ITK_OVERRIDE;

        bool IsValid() const ITK_OVERRIDE;
        bool IsUpToDate() const ITK_OVERRIDE;
        bool AreProfilesSufficientlyInsideImage() const ITK_OVERRIDE;

        void SetModelScheme(int index) ITK_OVERRIDE;
        int  GetModelScheme() ITK_OVERRIDE;

        TableType GetImageDataTable() const ITK_OVERRIDE;
        TableType GetDisplayModelTable() const ITK_OVERRIDE;
        TableType GetProcessingModelTable() const ITK_OVERRIDE;

        itkArrayType GetDisplayValues() const ITK_OVERRIDE;
        unsigned int GetNumberOfDisplayValues(bool includeImportedOptions=false) const ITK_OVERRIDE;
        std::string GetDisplayName(int index) const ITK_OVERRIDE;
        std::string GetDisplayNameShort(int index) const ITK_OVERRIDE;
        itkArray2DType GetDisplayRanges(bool includeImportedOptions=false) const ITK_OVERRIDE;

        itkArrayType GetParameterValues() const ITK_OVERRIDE;
        unsigned int GetNumberOfParameterValues() const ITK_OVERRIDE;
        std::string GetParameterName(int index) const ITK_OVERRIDE;
        std::string GetParameterNameShort(int index) const ITK_OVERRIDE;

    private:

        bool CalculateCalibrationValues();

        // GenerateData() throw ( ExceptionObject );
        //virtual void  GenerateData() ITK_OVERRIDE;

        const static unsigned int m_numberOfParameters=2;

        int m_modeIndex;

        ScalarType m_highMedian, m_lowMedian;

        LinearRegressionCalculator m_linearRegressionCalculator;

        bool m_readyToProcess;

    };

}


#endif //MYPROJECT_PROFILECALIBRATION_H
