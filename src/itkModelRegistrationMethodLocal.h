/* 
 * File:   itkModelRegistrationMethodLocal.h
 * Author: rap58
 *
 * Created on 25 May 2015, 15:21
 */

#ifndef ITKMODELREGISTRATIONMETHODLOCAL_H
#define	ITKMODELREGISTRATIONMETHODLOCAL_H


#include "itkProcessObject.h"
#include "itkImage.h"
#include "itkImageToSpatialObjectMetric.h"
#include "itkMultipleValuedNonLinearOptimizer.h"
#include "itkDataObjectDecorator.h"
#include <itkOptimisationCallback.h>
#include <itkMultipleModelMetricLocal.h>
#include <itkSingleModelMetricLocal.h>
#include <itkArray2D.h>

#include <itkThreeTierRectangularTransformLocal.h>
#include <itkEndostealRampTransformLocal.h>

#include <itkModelTransformBaseLocal.h>

#include <itkPowellOptimizer.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkLevenbergMarquardtOptimizer.h>
#include <vtkTable.h>
#include "itkArray.h"

//class CorticalBone;

namespace itk
{

    template< typename TScalar, unsigned int TDimension>
    class ModelRegistrationMethod:public ProcessObject {
    public:
        // standard class typedefs
        typedef ModelRegistrationMethod Self;
        typedef ProcessObject                          Superclass;
        typedef SmartPointer< Self >                   Pointer;
        typedef SmartPointer< const Self >             ConstPointer;

        // template
        typedef           TScalar                       ScalarType;
        const static unsigned int                       Dimension = TDimension;

        // other typedefs
        typedef itk::Array<ScalarType> ArrayType;

        // standard class macros
        itkNewMacro(Self);
        itkTypeMacro(ModelRegistrationMethod, ProcessObject);

        // pipeline type defs
        typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;


        // profile
        typedef itk::ProfileSpatialObject<3>                    ProfileType;

        // transform
        typedef ModelTransformBase<ScalarType, 1>                 TransformType;
        typedef  typename TransformType::Pointer                  TransformPointer;

        typedef itk::ThreeTierRectangularTransform< ScalarType >  RectTransformType;
        typedef  typename RectTransformType::Pointer              RectTransformPointer;
        typedef itk::EndostealRampTransform< ScalarType >         RampTransformType;
        typedef  typename RampTransformType::Pointer              RampTransformPointer;

        typedef itk::PowellOptimizer                                  PowellOptimiserType;
        typedef itk::OnePlusOneEvolutionaryOptimizer                  EvoluntionaryOptimiserType;
        typedef itk::LevenbergMarquardtOptimizer                      LMOptimiserType;

        //  Type of the metric.
        typedef MultipleModelMetric<TransformType>     MultipleMetricType;
        typedef typename MultipleMetricType::Pointer MetricMultiplePointer;

        typedef SingleModelMetric<TransformType>     SingleMetricType;
        typedef typename SingleMetricType::Pointer MetricSinglePointer;

        // Type for the output: Using Decorator pattern for enabling  the Transform to be passed in the data pipeline
        typedef  DataObjectDecorator< TransformType >      TransformOutputType;
        typedef typename TransformOutputType::Pointer      TransformOutputPointer;
        typedef  typename TransformType::ParametersType ParametersType;
        typedef  typename TransformType::PointArrayType PointArrayType;

        //  Type of the optimiser.
        typedef   MultipleValuedNonLinearOptimizer MultipleOptimizerType; // Note - generic one used
        typedef   SingleValuedNonLinearOptimizer SingleOptimizerType;

        typedef OptimiserCallback<MultipleOptimizerType>          MultipleOptimiserCallbackType;
        typedef typename MultipleOptimiserCallbackType::Pointer   MultipleCallbackOptimiserPointer;
        typedef OptimiserCallback<SingleOptimizerType>            SingleOptimiserCallbackType;
        typedef typename SingleOptimiserCallbackType::Pointer     SingleCallbackOptimiserPointer;

        typedef vtkSmartPointer<vtkTable>                               TableType;

        /** Smart Pointer type to a DataObject. */
        typedef typename DataObject::Pointer DataObjectPointer;

        typedef enum{
            kLevenbergMarquardtOptimizer =0,
            kPowellOptimizer,
            kOnePlusOneEvolutionaryOptimizer,
        } OptimiserIndex;
        typedef enum{
            kThreeTierRectangularModel =0,
            kEndostealRampModel,
        } ModelIndex;

        typedef enum { // display
            kBlank = -1,
            kTotalThicknessDisp =0, // returned by the TF's
            kCorticalThicknessDisp=1,
            kEndoThicknessDisp=2,
            kPeriPositionDisp=3,
            kEndoPositionDisp=4,
            kEndoCBPositionDisp=5,
            kEndoTBPositionDisp=6,
            kCorticalBoneDensityDisp=7,
            kTrabeculaeDensityDisp=8,
            kSoftTissueDensityDisp=9,
            kSurfaceThicknessDisp=10, // surface thickness for profileProperties, MSA density for display
            kMassSADensityDisp=10,
            kSigmaDisp=11,
            kErrorDisp=12,
            kImportImageBias=13, // calculated by the CB
            kImportImageSTD=14,
            kImportRMSError=15,
        } displayOptions;
        typedef enum { // param
            kXpParam=0,
            kxEParam,
            kYtbParam,
            kYcbParam,
            kYstParam,
            kSigmaParam,
            kWidthParam,
        } parameterOptions;


        //------------ General Methods ---------------//
        void Initialize() throw ( ExceptionObject ); // sets interconnects. call after all setters called
        void Initialize(ParametersType initialTransformParameters, ParametersType transformParameterScales, int modelIndex) throw ( ExceptionObject ); // sets interconnects. call after all setters called

        void TurnOffOptimizerObservation();
        void TurnOnOptimizerObservation();

        // set state
        void TurnOnFWHMMode(void) const;
        void TurnOffFWHMMode(void) const;
        bool GetFWHMMode() { return m_MultipleMetric->GetFWHMMode(); }; // both set so can call either multi or single

        void SetDynamicWeighting(bool status) const;
        bool GetDynamicWeighting() const;

        void SetWeightingMode(int index, double scale) const;
        void GetWeightingMode(int &index, double &scale) const;

        void TurnOnSigmaCorrection(double sigma);
        void TurnOffSigmaCorrection();
        ScalarType CalculateCBMv2SigmaAdjustedCBDensity(int modelIndex);

        virtual void FixParameter(const ScalarType value, const unsigned int index);
        virtual void FreeParameter(const unsigned int index);

        // get state
        bool IsValid() const;
        bool IsPtInsideImage(double point[3]) const;
        bool IsPtInsideImage() const;

        bool isSigmaFixed() const;
        bool isCBDensityFixed() const;

        //------------ Setters / Getters ------------//
        void SetProfile( ProfileType::Pointer spatialMovingObject);
        ProfileType::Pointer GetProfile();

        void SetOptimiserSelection(int index);
        int GetOptimiserSelection() const;

        void SetModelSelection(int index);
        int GetModelSelection() const;

        void SetInitialTransformParameters( ParametersType initialTransformParameters);
        itkGetConstReferenceMacro(InitialTransformParameters, ParametersType);

        void SetTransformParameterScales( ParametersType transformParameterScales );

        // get results
        itk::Array<double> GetDisplayValues() const;
        int GetNumberOfDisplayValues(bool importSet=false) const;
        std::string GetDisplayName(int index, bool importSet=false) const;
        std::string GetDisplayNameShort(int index, bool importSet=false) const;
        PointArrayType GetDisplayRanges(bool importSet=false) const;


        ScalarType GetCorticalDensity() const;
        ScalarType GetMeanCorticalDensity() const;
        ScalarType GetPeriostealEdgePosistion() const;
        ScalarType GetSigma() const;

        ScalarType GetCorticalBonePrecision() const;

        void SetMaxBMDDensity(ScalarType maxDensity);
        ScalarType GetMaxDensity() const;

        ScalarType GetErrorSum() const;
        ScalarType GetErrorMean() const;

        ParametersType GetParameterErrors() const;

        void SetCombinedParameters(ParametersType parameters);
        ParametersType GetCombinedParameters() const;
        std::string GetParameterName(int index) const;
        std::string GetParameterNameShort(int index) const;

        ScalarType GetNumberOfCombinedParameters() const;
        ScalarType GetNumberOfCombinedParameters(int modelIndex) const;

        itkGetConstReferenceMacro(LastTransformParameters, ParametersType);

        itk::Array<double> GetSampledProfilePositions() const;
        itk::Array<double> GetSampledImageValues() const; // calculates BMD values each time from stored HU values
        vtkSmartPointer<vtkDoubleArray> GetSampledModelValues() const; // calculates BMD values each time from stored HU values
        vtkSmartPointer<vtkDoubleArray> GetUnblurredModelValues() const;
        PointArrayType GetModelPoints();

        // for import error calculations
        ScalarType GetModelValue(ScalarType x);

        // display tables
        TableType GetImageDataTable() const;
        TableType GetDisplayModelTable() const;
        TableType GetProcessorModelTable() const;
        TableType GetDisplayModelTable(ParametersType parameters, int transformIndex, double offset) const;
        TableType GetWeightTable() const; // check correct weights

        //--------------- Pipeline Necessities --------------------//
        const TransformOutputType * GetOutput() const;

        using Superclass::MakeOutput;
        virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;
        virtual ModifiedTimeType GetMTime() const ITK_OVERRIDE;

    protected:
        ModelRegistrationMethod();
        virtual ~ModelRegistrationMethod() {}
        virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

        // Method invoked by the pipeline in order to trigger the computation of the registration. //
        virtual void  GenerateData() ITK_OVERRIDE;

        virtual void CreateSelectedOptimiser();

        virtual TransformPointer GetSelectedTransform(int modelIndex) const;

        ParametersType m_InitialTransformParameters;
        ParametersType m_LastTransformParameters;
        ParametersType m_TransformParameterScales;

    private:
        ModelRegistrationMethod(const Self &); //purposely not implemented
        void operator=(const Self &);                         //purposely not implemented

        MetricMultiplePointer                     m_MultipleMetric;
        MultipleOptimizerType::Pointer            m_MutlipleOptimiser;
        MultipleCallbackOptimiserPointer          m_MultipleOptimiserCallback;

        MetricSinglePointer                     m_SingleMetric;
        SingleOptimizerType::Pointer            m_SingleOptimiser;
        SingleCallbackOptimiserPointer          m_SingleOptimiserCallback;

        ProfileType::Pointer                    m_Profile;

        RectTransformPointer                    m_RectTransform;
        RampTransformPointer                    m_RampTransform;

        ScalarType                              m_Sigma;
        bool                                    m_SigmaCorrectionSet;

        int                                     m_OptimiserIndex;
        int                                     m_ModelIndex;

        const static unsigned int m_NumberOfDisplays = 13;

//    bool m_multipleOptimiserMode;
//    
//    int modelIndex, optimiserIndex;
//    
    };
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkModelRegistrationMethodLocal.hxx"
#endif

#endif	/* ITKMODELREGISTRATIONMETHODLOCAL_H */

