//
// Created by rap58 on 20/10/15.
//

#ifndef MYPROJECT_PROFILEPROCESSOR_H
#define MYPROJECT_PROFILEPROCESSOR_H

#include "itkProcessObject.h"
#include <vtkSmartPointer.h>

#include <itkImage.h>
#include <itkProfileSpatialObjectLocal.h>
#include <itkLinearInterpolateImageFunction.h>
//#include <itkBSplineInterpolateImageFunction.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkLevenbergMarquardtOptimizer.h>
#include <vtkTable.h>

#include <classifierTransform.h>

// defines

    namespace itk {

        class ProfileProcessor : public ProcessObject {

        public:

            typedef ProfileProcessor Self;
            typedef ProcessObject                          Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            typedef double                                                  ScalarType;
            const static unsigned int                                       Dimension = 3;
            typedef itk::ProfileSpatialObject<Dimension>                    ProfileType;
            typedef itk::TimeStamp                                          TimeType;
            typedef vtkSmartPointer<vtkDoubleArray>                         ArrayType;
            typedef itk::Array<double>                                      itkArrayType;
            typedef itk::Array<int>                                         itkIntArrayType;
            typedef  itk::Array2D<ScalarType>                               itkArray2DType;
            typedef vtkSmartPointer<vtkTable>                               TableType;

            ProfileProcessor();

            void SetProfile(ProfileType::Pointer profile);

            virtual void SetModelScheme(int index) = 0;

            virtual bool Process() = 0;

            virtual bool IsValid() const = 0;
            virtual bool IsUpToDate() const = 0;
            virtual bool AreProfilesSufficientlyInsideImage() const = 0;
            virtual bool IsPtInsideImage(bool includingTheStart=false) const;


            itkArrayType GetPositions();
            itkArrayType GetImageSamples(); // return the mean value if multiple

            virtual int  GetModelScheme() = 0;

            virtual TableType GetImageDataTable() const = 0;
            virtual TableType GetDisplayModelTable() const = 0;
            virtual TableType GetProcessingModelTable() const = 0;

            virtual itkArrayType GetDisplayValues() const = 0;
            virtual unsigned int GetNumberOfDisplayValues(bool includeImportedOptions=false) const = 0;
            virtual std::string GetDisplayName(int index) const = 0;
            virtual std::string GetDisplayNameShort(int index) const = 0;
            virtual itkArray2DType GetDisplayRanges(bool includeImportedOptions=false) const = 0;

            virtual itkArrayType GetParameterValues() const = 0;
            virtual unsigned int GetNumberOfParameterValues() const = 0;
            virtual std::string GetParameterName(int index) const = 0;
            virtual std::string GetParameterNameShort(int index) const = 0;

        protected:

            // GenerateData() throw ( ExceptionObject );
            //virtual void  GenerateData() ITK_OVERRIDE;

            virtual void UpdateState();

            ProfileType::Pointer m_profile;
            bool m_readyToProcess;
            bool m_requireAllInsideImage;

            mutable TimeType  m_ProcessedTimeStamp;


        };

    }

#endif //MYPROJECT_PROFILEPROCESSOR_H
