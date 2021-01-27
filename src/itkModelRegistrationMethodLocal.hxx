/* 
 * File:   itkModelRegistrationMethod.hxx
 * Author: rap58
 *
 * Created on 25 May 2015, 15:22
 */

#ifndef ITKMODELREGISTRATIONMETHOD_HXX
#define	ITKMODELREGISTRATIONMETHOD_HXX

#include "itkModelRegistrationMethodLocal.h"
#include "itkThreeTierRectangularTransformLocal.h"
#include <utilities.h>

#include <itkProfileSpatialObjectLocal.h>
#include <itkPowellOptimizer.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkLevenbergMarquardtOptimizer.h>
#include <itkNormalVariateGenerator.h>


namespace itk
{
    /** Constructor */
    template<typename TScalar, unsigned int TDimension >
    ModelRegistrationMethod< TScalar, TDimension >
    ::ModelRegistrationMethod() {
        this->SetNumberOfRequiredOutputs(1);    // for the Transform
        m_OptimiserIndex = 0;
        m_ModelIndex = 0;

        m_Profile                     = ITK_NULLPTR;          // has to be provided by the user.

        m_RectTransform               = RectTransformType::New();
        m_RampTransform               = RampTransformType::New();

        m_MutlipleOptimiser           = ITK_NULLPTR;
        m_SingleOptimiser             = ITK_NULLPTR;

        m_MultipleMetric              = MultipleMetricType::New();
        m_MultipleOptimiserCallback   = MultipleOptimiserCallbackType::New();
        m_MultipleOptimiserCallback->TurnOffOptimizerObservation();

        m_SingleMetric                = SingleMetricType::New();
        m_SingleOptimiserCallback     = SingleOptimiserCallbackType::New();
        m_SingleOptimiserCallback->TurnOffOptimizerObservation();

        CreateSelectedOptimiser();

        TransformPointer transform    = GetSelectedTransform(m_ModelIndex);
        m_InitialTransformParameters  = ParametersType(transform->GetNumberOfParameters());
        m_LastTransformParameters     = ParametersType(transform->GetNumberOfParameters());

        m_InitialTransformParameters.Fill(0.0f);
        m_LastTransformParameters.Fill(0.0f);

        m_TransformParameterScales    = ParametersType(transform->GetNumberOfParameters());
        m_TransformParameterScales.Fill(1.0f);

        m_Sigma                     = nan("1");
        m_SigmaCorrectionSet        = false;

        TransformOutputPointer transformDecorator = itkDynamicCastInDebugMode< TransformOutputType * >(this->MakeOutput(0).GetPointer() );
        this->ProcessObject::SetNthOutput( 0, transformDecorator.GetPointer() );
    }

    // Required to Run the filter
    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::Initialize() throw ( ExceptionObject ) {

        TransformPointer transform = GetSelectedTransform(m_ModelIndex);

        if (!m_Profile)  {
            itkExceptionMacro(<< "Profile was not set"); // includes Interpolator
            std::cerr<<"Error FixedImage or FixedProfile was not set"<<std::endl;
        }

        if ( m_TransformParameterScales.Size() != transform->GetNumberOfParameters() ) {
            itkExceptionMacro(<< "TF parameters scales of wrong dimension. It should be the same size as TF.");
            std::cerr<<"TF parameters scales of wrong dimension. It should be the same size as TF."<<std::endl;
        }

        CreateSelectedOptimiser();


        transform->SetParameters(m_InitialTransformParameters);

        // update profile offset - profile makes checks fr allowable distance
        ScalarType offset = transform->GetPeriostealEdgePosistion();
        m_Profile->SetPeriostealOffset(offset);

        if(m_OptimiserIndex == kLevenbergMarquardtOptimizer) {

            m_MultipleMetric->SetMovingTransform(transform);

            m_MutlipleOptimiser->SetInitialPosition(m_InitialTransformParameters);

            // create connections now that all is added
            m_MutlipleOptimiser->SetCostFunction(m_MultipleMetric);
            m_MutlipleOptimiser->SetInitialPosition(m_InitialTransformParameters);
            m_MutlipleOptimiser->SetScales(m_TransformParameterScales);

            m_MultipleMetric->Initialize();
        } else {

            m_SingleMetric->SetMovingTransform(transform);

            m_SingleOptimiser->SetInitialPosition(m_InitialTransformParameters);

            // create connections now that all is added
            m_SingleOptimiser->SetCostFunction(m_SingleMetric);
            m_SingleOptimiser->SetInitialPosition(m_InitialTransformParameters);
            m_SingleOptimiser->SetScales(m_TransformParameterScales);

            m_SingleMetric->Initialize();

        }

        // Connect the transform to the Decorator.
        TransformOutputType *transformOutput = static_cast< TransformOutputType * >( this->ProcessObject::GetOutput(0) );
        transformOutput->Set( transform );
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::Initialize(ParametersType initialTransformParameters, ParametersType transformParameterScales, int modelIndex) throw ( ExceptionObject ) {
        SetModelSelection(modelIndex);
        SetTransformParameterScales(transformParameterScales);
        SetInitialTransformParameters(initialTransformParameters);
        Initialize();
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::PrintSelf(std::ostream & os, Indent indent) const {
        Superclass::PrintSelf(os, indent);
        os << indent << "Metric: " << m_MultipleMetric.GetPointer() << std::endl;
        os << indent << "Optimizer: " << m_MutlipleOptimiser.GetPointer() << std::endl;
        os << indent << "Rect Transform: " << m_RectTransform.GetPointer() << std::endl;
        os << indent << "Ramp Transform: " << m_RampTransform.GetPointer() << std::endl;
        os << indent << "Moving SpatialObject: " << m_Profile.GetPointer() << std::endl;
        os << indent << "Initial Transform Parameters: " << m_InitialTransformParameters << std::endl;
        os << indent << "Last    Transform Parameters: " << m_LastTransformParameters << std::endl;
        os << indent << "Transform Parameter Scales: " << m_TransformParameterScales << std::endl;
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::GenerateData() {

        if(m_OptimiserIndex == kLevenbergMarquardtOptimizer) {

            // metric reset iteration count
            m_MultipleMetric->RestartIterationCount();

            try {  // do the optimization
                m_MutlipleOptimiser->StartOptimization();
            }  catch ( ExceptionObject & err ) { // An error has occurred in the optimisation.

                std::cerr<<"m_Optimizer->StartOptimization() exception caught"<<std::endl;
                m_LastTransformParameters = m_MutlipleOptimiser->GetCurrentPosition(); // Update the parameters
                throw err; // Pass exception to caller
            }

            // get the results
            m_LastTransformParameters = m_MutlipleOptimiser->GetCurrentPosition();
            TransformPointer transform = GetSelectedTransform(m_ModelIndex); transform->SetParameters(m_LastTransformParameters);
            m_Profile->SetPeriostealOffset(transform->GetPeriostealEdgePosistion()); m_MultipleMetric->Initialize();
            m_MultipleMetric->GetValue(m_LastTransformParameters); // update transform to optimal solution
        } else {

            // metric reset iteration count
            m_SingleMetric->RestartIterationCount();

            try {  // do the optimization
                m_SingleOptimiser->StartOptimization();
            }  catch ( ExceptionObject & err ) { // An error has occurred in the optimisation.

                std::cerr<<"m_Optimizer->StartOptimization() exception caught"<<std::endl;
                m_LastTransformParameters = m_SingleOptimiser->GetCurrentPosition(); // Update the parameters
                throw err; // Pass exception to caller
            }

            // get the results
            m_LastTransformParameters = m_SingleOptimiser->GetCurrentPosition();
            TransformPointer transform = GetSelectedTransform(m_ModelIndex); transform->SetParameters(m_LastTransformParameters);
            m_Profile->SetPeriostealOffset(transform->GetPeriostealEdgePosistion()); m_MultipleMetric->Initialize();
            m_SingleMetric->GetValue(m_LastTransformParameters); // update transform to optimal solution

        }


    }

    template<typename TScalar, unsigned int TDimension >
    const typename ModelRegistrationMethod< TScalar, TDimension >::TransformOutputType *
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetOutput() const {
        return static_cast< const TransformOutputType * >( this->ProcessObject::GetOutput(0) );
    }

    template<typename TScalar, unsigned int TDimension >
    DataObject::Pointer ModelRegistrationMethod< TScalar, TDimension >
    ::MakeOutput(DataObjectPointerArraySizeType output) {
        switch ( output ) {
            case 0:
                return TransformOutputType::New().GetPointer();
                break;
            default:
            itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs");
                return ITK_NULLPTR;
        }
    }

    template<typename TScalar, unsigned int TDimension >
    ModifiedTimeType ModelRegistrationMethod< TScalar, TDimension >
    ::GetMTime() const {
        ModifiedTimeType mtime = Superclass::GetMTime();
        ModifiedTimeType m;

        TransformPointer transform = GetSelectedTransform(m_ModelIndex);

        // get time last changed - either self or component object. Remember bigger = newer
        if ( transform ) {
            m = transform->GetMTime();
            mtime = ( m > mtime ? m : mtime );
        }

        if ( m_Profile ) {
            m = m_Profile->GetMTime();
            mtime = ( m > mtime ? m : mtime );
        }

        if ( m_MultipleMetric && m_OptimiserIndex == kLevenbergMarquardtOptimizer) {
            m = m_MultipleMetric->GetMTime();
            mtime = ( m > mtime ? m : mtime );
        } else if(m_SingleMetric) {
            m = m_SingleMetric->GetMTime();
            mtime = ( m > mtime ? m : mtime );
        }

        if ( m_MutlipleOptimiser && m_OptimiserIndex == kLevenbergMarquardtOptimizer)  {
            m = m_MutlipleOptimiser->GetMTime();
            mtime = ( m > mtime ? m : mtime );
        } else if ( m_SingleOptimiser )  {
            m = m_SingleOptimiser->GetMTime();
            mtime = ( m > mtime ? m : mtime );
        }

        return mtime;
    }

    // Getters and Setters
    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::SetProfile( ProfileType::Pointer spatialMovingObject) {
        this->m_Profile = spatialMovingObject;
        m_MultipleMetric->SetProfile(m_Profile);
        m_SingleMetric->SetProfile(m_Profile);

        m_RectTransform->SetProfileLength(m_Profile->GetProfileLength());
        m_RampTransform->SetProfileLength(m_Profile->GetProfileLength());

        m_RectTransform->SetProfileEdgeProperties(m_Profile->GetPeriostealEdgeRatio(), m_Profile->GetMaximumPeriostealOffset());
        m_RampTransform->SetProfileEdgeProperties(m_Profile->GetPeriostealEdgeRatio(), m_Profile->GetMaximumPeriostealOffset());
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ProfileType::Pointer
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetProfile() {
        this->m_Profile;
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::SetOptimiserSelection(int index) {
        m_OptimiserIndex = index;
    }

    template<typename TScalar, unsigned int TDimension >
    int ModelRegistrationMethod< TScalar, TDimension >
    ::GetOptimiserSelection() const {
        return m_OptimiserIndex;
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::SetModelSelection(int index) {
        m_ModelIndex = index;
    }

    template<typename TScalar, unsigned int TDimension >
    int ModelRegistrationMethod< TScalar, TDimension >
    ::GetModelSelection() const {
        return m_ModelIndex;
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::SetInitialTransformParameters( ParametersType initialTransformParameters) {

        TransformPointer transform = GetSelectedTransform(m_ModelIndex);

        if(initialTransformParameters != m_InitialTransformParameters || initialTransformParameters != transform->GetParameters()) {

            int parameterDimension = transform->GetNumberOfParameters();

            if ( initialTransformParameters.Size() != transform->GetNumberOfParameters() ) {
                itkExceptionMacro(<< "Initial TF parameters of wrong dimension. It should be the same size as TF.");
                std::cerr<<"Initial TF parameters of wrong dimension. It should be the same size as TF."<<std::endl;
            }

            this->m_InitialTransformParameters = initialTransformParameters;

            m_LastTransformParameters.SetSize(parameterDimension);

            for(int i=0; i< parameterDimension; i++) {
                m_LastTransformParameters[i] = m_InitialTransformParameters[i];
            }

        }

    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::SetTransformParameterScales( ParametersType transformParameterScales) {
        if(m_TransformParameterScales!= transformParameterScales) {
            this->m_TransformParameterScales = transformParameterScales;
        } else {
            //std::cerr<<"Warning parameter scales are unchanged"<<endl;
        }
    }

    // create objects
    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::CreateSelectedOptimiser() {

        if(m_OptimiserIndex == kLevenbergMarquardtOptimizer) {

            LMOptimiserType::Pointer optimiser = LMOptimiserType::New();
            optimiser->UseCostFunctionGradientOff();
            //LMOptimiser->GetOptimizer()->get_covariance();

            double gradientTolerance        =  1e-4;  // Gradient magnitude tolerance
            double positionTolerance        =  1e-4;  // Search space tolerance
            double epsilonFunctionValue     =  1e-5;  // Step 1e-2 gives lower duration and failure rate
            int    maxNumberOfIterations    =   2000;  // Maximum number of iterations

            optimiser->SetNumberOfIterations( maxNumberOfIterations );
            optimiser->SetValueTolerance( positionTolerance );
            optimiser->SetGradientTolerance( gradientTolerance );
            optimiser->SetEpsilonFunction( epsilonFunctionValue );
            //optimiser->GetOptimizer()->set_verbose(false); // TODO - turn off 'get JtJ warnings'

            m_MutlipleOptimiser = (MultipleOptimizerType::Pointer)optimiser.GetPointer();

            m_MultipleOptimiserCallback->SetOptimizer( m_MutlipleOptimiser );

        } else if(m_OptimiserIndex == kPowellOptimizer) {

            PowellOptimiserType::Pointer optimiser = PowellOptimiserType::New();
            optimiser->MaximizeOff();
            optimiser->SetStepLength( 1e-4 );         // increase towards 1 to reduce the processing time
            optimiser->SetStepTolerance( 1e-4 );
            optimiser->SetValueTolerance( 1e-5 );
            optimiser->SetMaximumIteration( 2000 );

            m_SingleOptimiser = (SingleOptimizerType::Pointer)optimiser.GetPointer();

            m_SingleOptimiserCallback->SetOptimizer( m_SingleOptimiser );

        } else if(m_OptimiserIndex == kOnePlusOneEvolutionaryOptimizer) {

            itk::Statistics::NormalVariateGenerator::Pointer generator = itk::Statistics::NormalVariateGenerator::New();
            generator->Initialize(12345); // TODO - decide on seedpoint and improve parameters
            EvoluntionaryOptimiserType::Pointer optimiser = EvoluntionaryOptimiserType::New();
            optimiser->MaximizeOff();
            optimiser->SetNormalVariateGenerator( generator );
            optimiser->SetInitialRadius(0.1);
            optimiser->SetGrowthFactor(1.0); optimiser->SetShrinkFactor(0.95);
            optimiser->SetEpsilon(1e-4);                               // 0.001
            optimiser->SetMaximumIteration( 500 );

            m_SingleOptimiser = (SingleOptimizerType::Pointer)optimiser.GetPointer();

            m_SingleOptimiserCallback->SetOptimizer( m_SingleOptimiser );

        } else {
            std::cerr<<"Error invalid optimiser index"<<std::endl;
        }

    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::TransformPointer
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetSelectedTransform(int modelIndex) const {
        TransformPointer transform;
        if(modelIndex == kThreeTierRectangularModel) {
            transform = m_RectTransform;
        } else {
            transform = m_RampTransform;
        }

        return transform;
    }


    // get results
    template<typename TScalar, unsigned int TDimension >
    itk::Array<double>  ModelRegistrationMethod< TScalar, TDimension >
    ::GetDisplayValues() const {

        // if invalid expect values to be set to nan except in the case that they are fixed
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        ParametersType values =  transform->GetDisplayValues();

        values[kErrorDisp] = GetErrorMean();

        return values;
    }

    template<typename TScalar, unsigned int TDimension >
    int ModelRegistrationMethod< TScalar, TDimension >
    ::GetNumberOfDisplayValues(bool includeImportedOptions) const {
        if(includeImportedOptions) {
            return m_NumberOfDisplays+3;
        } else {
            return m_NumberOfDisplays;
        }
    }

    template<typename TScalar, unsigned int TDimension >
    std::string ModelRegistrationMethod< TScalar, TDimension >
    ::GetDisplayName(int index, bool importSet) const {

        if(index ==kBlank) {
            return "Blank";
        } else if(index == kTotalThicknessDisp) {
            return "Average Cortical Thickness";
        } else if(index == kCorticalThicknessDisp) {
            return "Dense Cortical Thickness";
        } else if(index == kEndoThicknessDisp) {
            return "Endosteal Thickness";
        } else if(index == kPeriPositionDisp) {
            return "Periosteal Position";
        } else if(index == kEndoPositionDisp) {
            return "Endosteal Position";
        } else if(index == kEndoCBPositionDisp) {
            return "Inner Endosteal Position";
        } else if(index == kEndoTBPositionDisp) {
            return "Outer Endosteal Position";
        } else if(index == kCorticalBoneDensityDisp) {
            return "Cortical Density";
        } else if(index == kTrabeculaeDensityDisp) {
            return "Trabecular Density";
        } else if(index == kSoftTissueDensityDisp) {
            return "Soft Tissue Density";
        } else if(index == kMassSADensityDisp) {
            return "Mass Density SA";
        } else if(index == kSigmaDisp) {
            return "Sigma";
        } else if(index == kErrorDisp) {
            return "ABS Mean Error";
        } else if(index == kImportImageBias && importSet) {
            return "Import Bias";
        } else if(index == kImportImageSTD && importSet) {
            return "Import STD";
        } else if(index == kImportRMSError && importSet) {
            return "Import RMS";
        } else {
            std::cerr<<"Error: ModelRegistrationMethod::getDisplayName() invalid display index: "<<index<<endl;
            return "Invalid Selection";
        }
    }

    template<typename TScalar, unsigned int TDimension >
    std::string ModelRegistrationMethod< TScalar, TDimension >
    ::GetDisplayNameShort(int index, bool importSet) const {

        if(index == kBlank) {
            return "";
        } else if(index == kTotalThicknessDisp) {
            return "T<sub>Avg</sub>";
        } else if(index == kCorticalThicknessDisp) {
            return "T<sub>CB</sub>";
        } else if(index == kEndoThicknessDisp) {
            return "T<sub>Endo</sub>"; // Y<sub>CB</sub>
        } else if(index == kPeriPositionDisp) {
            return "X<sub>Peri</sub>";
        } else if(index == kEndoPositionDisp) {
            return "X<sub>Endo</sub>";
        } else if(index == kEndoCBPositionDisp) {
            return "X<sub>InnerEndo</sub>";
        } else if(index == kEndoTBPositionDisp) {
            return "X<sub>OuterEndo</sub>";
        } else if(index == kCorticalBoneDensityDisp) {
            return "Y<sub>CB</sub>";
        } else if(index == kTrabeculaeDensityDisp) {
            return "Y<sub>TB</sub>";
        } else if(index == kSoftTissueDensityDisp) {
            return "Y<sub>ST</sub>";
        } else if(index == kMassSADensityDisp) {
            return "MSA";
        } else if(index == kSigmaDisp) {
            return "\u03C3";
        } else if(index == kErrorDisp) {
            return "ABS ME";
        }  else if(index == kImportImageBias && importSet) {
            return "Import Bias";
        } else if(index == kImportImageSTD && importSet) {
            return "Import STD";
        } else if(index == kImportRMSError && importSet) {
            return "Import RMS";
        } else {
            std::cerr<<"Error: CorticalBone::getDisplayNameShort() invalid display index: "<<index<<endl;
            return "Invalid Selection";
        }
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::PointArrayType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetDisplayRanges(bool includeImportedOptions) const {

        ScalarType profileLength = m_Profile->GetProfileLength();
        int n = (includeImportedOptions) ? m_NumberOfDisplays+3 : m_NumberOfDisplays;
        PointArrayType displayRanges(n, 2);// +3 for imported data

        displayRanges.put(kTotalThicknessDisp,0,0.0);        displayRanges.put(kTotalThicknessDisp,1,6.0);
        displayRanges.put(kCorticalThicknessDisp,0,0.0);       displayRanges.put(kCorticalThicknessDisp,1,6.0);
        displayRanges.put(kEndoThicknessDisp,0,0.0);           displayRanges.put(kEndoThicknessDisp,1,6.0);
        displayRanges.put(kPeriPositionDisp,0,0-profileLength/8.0);      displayRanges.put(kPeriPositionDisp,1,profileLength/8.0); // todo get from profileModel
        displayRanges.put(kEndoPositionDisp,0,0-profileLength/8.0);      displayRanges.put(kEndoPositionDisp,1,profileLength/4.0);
        displayRanges.put(kEndoCBPositionDisp,0,0-profileLength/8.0);      displayRanges.put(kEndoCBPositionDisp,1,profileLength/4.0);
        displayRanges.put(kEndoTBPositionDisp,0,0-profileLength/8.0);      displayRanges.put(kEndoTBPositionDisp,1,profileLength/4.0);

        // in BMD values - mg/cm^3
        displayRanges.put(kCorticalBoneDensityDisp,0,400);      displayRanges.put(kCorticalBoneDensityDisp,1,1600.0);
        displayRanges.put(kTrabeculaeDensityDisp,0,-200.0);     displayRanges.put(kTrabeculaeDensityDisp,1,600.0);
        displayRanges.put(kSoftTissueDensityDisp,0,-200.0);     displayRanges.put(kSoftTissueDensityDisp,1,400.0);
        displayRanges.put(kSurfaceThicknessDisp,0,0.0);         displayRanges.put(kSurfaceThicknessDisp,1,8000.0);
        displayRanges.put(kSigmaDisp,0,0.0);                    displayRanges.put(kSigmaDisp,1,3.0);
        displayRanges.put(kErrorDisp,0,0.0);                    displayRanges.put(kErrorDisp,1,200.0); // 10% of CB

        if(includeImportedOptions) {
            displayRanges.put(kErrorDisp + 1, 0, -200.0);
            displayRanges.put(kErrorDisp + 1, 1, 200.0); // +/-10% CB Max
            displayRanges.put(kErrorDisp + 2, 0, 0.0);
            displayRanges.put(kErrorDisp + 2, 1, 600.0); // +/-30% CB Max
            displayRanges.put(kErrorDisp + 3, 0, 0.0);
            displayRanges.put(kErrorDisp + 3, 1, 600.0); // 30% CB Max
        }


        return displayRanges;
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ScalarType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetCorticalDensity() const {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->GetCorticalDensity();
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ScalarType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetMeanCorticalDensity() const {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->GetMeanCorticalDensity();
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ScalarType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetPeriostealEdgePosistion() const {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->GetPeriostealEdgePosistion();
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ScalarType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetSigma() const {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->GetSigma();
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ScalarType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetCorticalBonePrecision() const {
        if(m_OptimiserIndex==kLevenbergMarquardtOptimizer) {
            itk::LevenbergMarquardtOptimizer::Pointer optimiser = static_cast<itk::LevenbergMarquardtOptimizer*>(m_MutlipleOptimiser.GetPointer());

            //const vnl_matrix<double> covarianceMatrix = optimiser->GetOptimizer()->get_covariance(); // not implemented!

            const vnl_matrix<double> inverseCovarianceApproximation = optimiser->GetOptimizer()->get_JtJ();
            const vnl_vector<double> diagonal = inverseCovarianceApproximation.get_diagonal();
            TransformPointer transform = GetSelectedTransform(m_ModelIndex);

            // Todo - select cortical density or periosteal edge as the parameter to use for precision weighting

            int cbDensityIndex = transform->GetCorticalDensityFlexIndex(), periostealEdgeIndex = transform->GetPeriostealEdgeFlexIndex();
            ScalarType precision = (periostealEdgeIndex!=-1) ? diagonal[periostealEdgeIndex] : nan("1"); // inverse if using the cortical density precision
            //ScalarType precision = (cbDensityIndex!=-1) ? diagonal[cbDensityIndex] : nan("1"); // inverse if using the cortical density precision

            //cerr<<"CB density index="<<cbDensityIndex<<", precision="<<cbDensityPrecision<<", diag range=["<<diagonal.min_value()<<","<<diagonal.max_value()<<"]"<<endl;
            return precision;
        } else {
            return 1.0; //nan("1"); //todo is there an estimate of precision that can be provided?
        }
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::SetMaxBMDDensity(ScalarType maxDensity) {

        // note should be the image value equvalent to 1600mg/cm^3 given the current calibration. Calibration performed by cortical bone though
        m_RampTransform->SetMaxDensity(maxDensity);
        m_RectTransform->SetMaxDensity(maxDensity);
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ScalarType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetMaxDensity() const {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->GetMaxDensity();
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ScalarType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetErrorSum() const {

        if(m_OptimiserIndex == kLevenbergMarquardtOptimizer) {
            return this->m_MultipleMetric->GetAbsErrorSum(this->m_LastTransformParameters);
        } else {
            return this->m_SingleMetric->GetAbsErrorSum(this->m_LastTransformParameters);
        }
        //return this->m_Metric->GetMeanAbsError(this->m_LastTransformParameters);
        //return this->m_Metric->GetSqrErrorSum(this->m_LastTransformParameters);
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ScalarType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetErrorMean() const {
        //return this->m_Metric->GetSqrErrorMean(this->m_LastTransformParameters);
        if(m_OptimiserIndex == kLevenbergMarquardtOptimizer) {
            return this->m_MultipleMetric->GetAbsErrorMean(this->m_LastTransformParameters);
        } else {
            return this->m_SingleMetric->GetAbsErrorMean(this->m_LastTransformParameters);
        }
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ParametersType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetParameterErrors() const {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->GetParameterErrors();
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::SetCombinedParameters(ParametersType parameters) {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        transform->SetCombinedParameters(parameters);
        SetInitialTransformParameters(transform->GetParameters());
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ParametersType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetCombinedParameters() const {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->GetCombinedParameters();
    }

    template<typename TScalar, unsigned int TDimension >
    std::string ModelRegistrationMethod< TScalar, TDimension >
    ::GetParameterNameShort(int index) const {

        if(index == kXpParam) {
            return "X\u209A";
        } if(index == kxEParam) {
            return "X\u2091";
        } if(index == kYtbParam) {
            return "Y\u209Cb";
        } if(index == kYcbParam) {
            return "Ycb"; // possibly specific to QT
        } if(index == kYstParam) {
            return "Y\u209B\u209C";
        } if(index == kSigmaParam) {
            return "\u03C3";
        } if(index == kWidthParam) {
            return "X\u2091w";
        } else {
            std::cerr<<"Error: ModelRegistrationMethod::GetParameterNameShort() invalid index = "<<index<<endl;
            return "Invalid Selection";
        }
    }

    template<typename TScalar, unsigned int TDimension >
    std::string ModelRegistrationMethod< TScalar, TDimension >
    ::GetParameterName(int index) const {

        if(index == kXpParam) {
            return "Periosteal Edge Parameter";
        } if(index == kxEParam) {
            return "Endosteal Edge Parameter";
        } if(index == kYtbParam) {
            return "Trabecular Bone Density Parameter";
        } if(index == kYcbParam) {
            return "Cortical Bone Density Parameter";
        } if(index == kYstParam) {
            return "Soft Tissue Density Parameter";
        } if(index == kSigmaParam) {
            return "Sigma Parameter";
        } if(index == kWidthParam) {
            return "Endosteal Width Parameter";
        } else {
            std::cerr<<"Error: ModelRegistrationMethod::GetParameterName() invalid index = "<<index<<endl;
            return "Invalid Selection";
        }
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ScalarType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetNumberOfCombinedParameters() const {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->GetNumberOfCombinedParameters();
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::ScalarType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetNumberOfCombinedParameters(int modelIndex) const {
        TransformPointer transform = GetSelectedTransform(modelIndex);
        return transform->GetNumberOfCombinedParameters();
    }

    // the GetSampledXXProfileValues reutrns empty NaN arrays if iamge is out of bounds
    template<typename TScalar, unsigned int TDimension >
    itk::Array<double> ModelRegistrationMethod< TScalar, TDimension >
    ::GetSampledProfilePositions() const {

        // sample if even out of bounds
        ArrayType array = m_Profile->GetPositions();
        return array;

    }

    template<typename TScalar, unsigned int TDimension >
    itk::Array<double> ModelRegistrationMethod< TScalar, TDimension >
    ::GetSampledImageValues() const {

        if(!IsPtInsideImage()) {

            itk::Array<double> array = itk::Array<double>(1); array.Fill(nan("1")); // name - "Image Profile"
            return array;

        } else {
            ArrayType array = m_Profile->GetValues();
            return array;
        }
    }

    template<typename TScalar, unsigned int TDimension >
    vtkSmartPointer<vtkDoubleArray> ModelRegistrationMethod< TScalar, TDimension >
    ::GetSampledModelValues() const { // calibration corrected = units are mg/cm^3

        ArrayType array;
        if(!IsPtInsideImage()) {

            vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
            array->SetName("Blurred Model"); array->SetNumberOfValues(1); array->SetValue(0, nan("1"));
            return array;

        } else if(m_OptimiserIndex == kLevenbergMarquardtOptimizer) {
            array = this->m_MultipleMetric->GetModelProfileSamples();
        } else if(m_OptimiserIndex != kLevenbergMarquardtOptimizer) {
            array = this->m_SingleMetric->GetModelProfileSamples();
        }

        return Utilities::convertItkToVtkArray(array, "Blurred Model");
    }

    template<typename TScalar, unsigned int TDimension >
    vtkSmartPointer<vtkDoubleArray> ModelRegistrationMethod< TScalar, TDimension >
    ::GetUnblurredModelValues() const {

        ArrayType array;
        if(m_OptimiserIndex == kLevenbergMarquardtOptimizer) {
            array = this->m_MultipleMetric->GetUnblurredModelSamples();
        } else if(m_OptimiserIndex != kLevenbergMarquardtOptimizer) {
            array = this->m_SingleMetric->GetUnblurredModelSamples();
        }
        return Utilities::convertItkToVtkArray(array, "Model");
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::PointArrayType ModelRegistrationMethod< TScalar, TDimension >
    ::GetModelPoints() {

        if(!IsPtInsideImage()) {
            PointArrayType array(1, 2); array.Fill(nan("1"));
            return array;
        } else {
            TransformPointer transform = GetSelectedTransform(m_ModelIndex);
            return transform->GetModelPoints(this->m_Profile->GetStartXValue());
        }
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod < TScalar, TDimension >::ScalarType  ModelRegistrationMethod< TScalar, TDimension >
    ::GetModelValue(ScalarType x) {

        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->GetUntransformedPoint(x);
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::TableType  ModelRegistrationMethod< TScalar, TDimension >
    ::GetImageDataTable() const {

        TableType table = TableType::New();
        table->SetNumberOfRows(m_Profile->GetNumberOfSamples());
        if(IsPtInsideImage()) {

            vtkSmartPointer<vtkDoubleArray> positions = Utilities::convertItkToVtkArray(GetSampledProfilePositions(), "Positions");
            vtkSmartPointer<vtkDoubleArray> values = Utilities::convertItkToVtkArray(GetSampledImageValues(), "Image Values");

            table->AddColumn(positions);
            table->AddColumn(values);
        }
        return table;
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::TableType  ModelRegistrationMethod< TScalar, TDimension >
    ::GetDisplayModelTable() const {

        TableType table = TableType::New();
        table->SetNumberOfRows(m_Profile->GetNumberOfSamples());
        if(IsPtInsideImage()) {

            TransformPointer transform = GetSelectedTransform(m_ModelIndex);

            itk::Array2D<ScalarType> modelPoints = transform->GetModelPoints(this->m_Profile->GetStartXValue());

            int n = modelPoints.rows();

            table->SetNumberOfRows(n);

            vtkSmartPointer<vtkDoubleArray> positions = vtkSmartPointer<vtkDoubleArray>::New();
            vtkSmartPointer<vtkDoubleArray> values = vtkSmartPointer<vtkDoubleArray>::New();
            positions->SetNumberOfValues(n); positions->SetName("Positions");
            values->SetNumberOfValues(n); values->SetName("Model Values");

            for(int i=0; i<n; i++) {

                positions->SetValue(i,modelPoints.get(i,0));
                double density = modelPoints.get(i,1);

                values->SetValue(i,density);
            }

            table->AddColumn(positions); table->AddColumn(values);
        }
        return table;
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::TableType  ModelRegistrationMethod< TScalar, TDimension >
    ::GetDisplayModelTable(ParametersType parameters, int transformIndex, double offset) const {

        TransformPointer transform = GetSelectedTransform(transformIndex);
        itk::Array2D<ScalarType> modelPoints = transform->GetModelPoints(this->m_Profile->GetStartXValue(), parameters);

        int n = modelPoints.rows();

        vtkSmartPointer<vtkDoubleArray> positions = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> values = vtkSmartPointer<vtkDoubleArray>::New();
        positions->SetNumberOfValues(n); positions->SetName("Positions");
        values->SetNumberOfValues(n); values->SetName("Imported Model Values");

        for(int i=0; i<n; i++) {

            positions->SetValue(i,modelPoints.get(i,0)-offset);

            double density = modelPoints.get(i,1);
            values->SetValue(i,density);
        }

        // create table
        TableType table = TableType::New();
        table->SetNumberOfRows(n);
        table->AddColumn(positions); table->AddColumn(values);
        return table;
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar, TDimension >::TableType  ModelRegistrationMethod< TScalar, TDimension >
    ::GetWeightTable() const {

        TableType table = TableType::New();
        table->SetNumberOfRows(m_Profile->GetNumberOfSamples());

        if(!IsPtInsideImage()) { return table; }


        vtkSmartPointer<vtkDoubleArray> positions;
        vtkSmartPointer<vtkDoubleArray> weights;
        if(m_OptimiserIndex==kLevenbergMarquardtOptimizer) {
            positions = Utilities::convertItkToVtkArray(GetSampledProfilePositions(), "Positions");
            weights = Utilities::convertItkToVtkArray(m_MultipleMetric->GetWeightingArray(), "Weights");


        } else {
            positions = Utilities::convertItkToVtkArray(GetSampledProfilePositions(), "Positions");
            weights = Utilities::convertItkToVtkArray(m_SingleMetric->GetWeightingArray(), "Weights");

        }

        double maxValue=0; int length=weights->GetSize();
        for(int i=0; i<length; i++) {
            double weightValue=weights->GetValue(i);
            maxValue=(maxValue>weightValue)?maxValue:weightValue;
        }
        double scale = 1000/maxValue;
        for(int i=0; i<length; i++) {
            weights->SetValue(i, weights->GetValue(i)*scale);
        }

        table->AddColumn(positions);
        table->AddColumn(weights);

        return table;
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod< TScalar,  TDimension>::TableType
    ModelRegistrationMethod< TScalar, TDimension >
    ::GetProcessorModelTable() const {
        TableType table = TableType::New();
        table->SetNumberOfRows(m_Profile->GetNumberOfSamples());
        if(IsPtInsideImage()) {
            vtkSmartPointer<vtkDoubleArray> positions = Utilities::convertItkToVtkArray(GetSampledProfilePositions(), "Positions");
            vtkSmartPointer<vtkDoubleArray> values = GetSampledModelValues();

            table->AddColumn(positions);
            table->AddColumn(values);
        }

        return table;
    }

    // get state
    template<typename TScalar, unsigned int TDimension >
    bool ModelRegistrationMethod< TScalar, TDimension >
    ::IsValid() const {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->IsValid();
    }

    template<typename TScalar, unsigned int TDimension >
    bool ModelRegistrationMethod< TScalar, TDimension >
    ::IsPtInsideImage(double point[3]) const {
        return m_Profile->IsPtInsideImage(point);
    }

    template<typename TScalar, unsigned int TDimension >
    bool ModelRegistrationMethod< TScalar, TDimension >
    ::IsPtInsideImage() const { // is mesh edge inside image bounds
        return m_Profile->IsPtInsideImage();
    }

    template<typename TScalar, unsigned int TDimension >
    bool ModelRegistrationMethod< TScalar, TDimension >
    ::isSigmaFixed() const {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->isSigmaFixed();
    }

    template<typename TScalar, unsigned int TDimension >
    bool ModelRegistrationMethod< TScalar, TDimension >
    ::isCBDensityFixed() const {
        TransformPointer transform = GetSelectedTransform(m_ModelIndex);
        return transform->isCorticalDensityFixed();
    }

    //------------ General Methods ------------------//
    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::TurnOffOptimizerObservation(){
        m_MultipleOptimiserCallback->TurnOffOptimizerObservation();
        m_SingleOptimiserCallback->TurnOffOptimizerObservation();
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::TurnOnOptimizerObservation() {
        m_MultipleOptimiserCallback->TurnOnOptimizerObservation();
        m_SingleOptimiserCallback->TurnOnOptimizerObservation();
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::TurnOnFWHMMode(void) const {
        m_MultipleMetric->TurnOnFWHMMode();
        m_SingleMetric->TurnOnFWHMMode();
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::TurnOffFWHMMode(void) const {
        m_MultipleMetric->TurnOffFWHMMode();
        m_SingleMetric->TurnOffFWHMMode();
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::SetDynamicWeighting(bool status) const {
        m_MultipleMetric->SetDynamicWeighting(status);
        m_SingleMetric->SetDynamicWeighting(status);

    }

    template<typename TScalar, unsigned int TDimension >
    bool ModelRegistrationMethod< TScalar, TDimension >
    ::GetDynamicWeighting() const {
        if(m_OptimiserIndex==kLevenbergMarquardtOptimizer) {
            return m_MultipleMetric->GetDynamicWeighting();
        } else {
            return m_SingleMetric->GetDynamicWeighting();
        }
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::SetWeightingMode(int index, double scale) const {
        m_MultipleMetric->SetWeightingMode(index, scale);
        m_SingleMetric->SetWeightingMode(index, scale);
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::GetWeightingMode(int& index, double& scale) const {
        if(m_OptimiserIndex==kLevenbergMarquardtOptimizer) {
            index = m_MultipleMetric->GetWeightingMode();
            scale = m_MultipleMetric->GetWeightingScale();
        } else {
            index = m_SingleMetric->GetWeightingScale();
            scale = m_SingleMetric->GetWeightingScale();
        }
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::TurnOnSigmaCorrection(double sigma) {

        m_SigmaCorrectionSet = true;
        m_Sigma = sigma;

    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::TurnOffSigmaCorrection() {
        m_SigmaCorrectionSet = false;
        m_Sigma = nan("1");
    }

    template<typename TScalar, unsigned int TDimension >
    typename ModelRegistrationMethod < TScalar, TDimension >::ScalarType  ModelRegistrationMethod< TScalar, TDimension >
    ::CalculateCBMv2SigmaAdjustedCBDensity(int modelIndex) { // calculate the sigma adjusted cortical density

        TransformPointer transform = GetSelectedTransform(modelIndex);
        ScalarType peakProfileValue;
        if(m_OptimiserIndex == kLevenbergMarquardtOptimizer) {
            peakProfileValue = m_MultipleMetric->GetFWHMCBDensity();
        } else {
            peakProfileValue = m_SingleMetric->GetFWHMCBDensity();
        }

        ScalarType sigma = m_Profile->GetProfileSigma();
        return transform->CalculateSigmaAdjustedCorticalDensity(peakProfileValue, sigma);

    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::FixParameter(const ScalarType value, const unsigned int index) {

        m_RectTransform->FixParameter(value, index);
        m_RampTransform->FixParameter(value, index);
    }

    template<typename TScalar, unsigned int TDimension >
    void ModelRegistrationMethod< TScalar, TDimension >
    ::FreeParameter(const unsigned int index) {

        m_RectTransform->FreeParameter(index);
        m_RampTransform->FreeParameter(index);
    }
} // end namespace itk

#endif	/* ITKMODELREGISTRATIONMETHOD_HXX */

