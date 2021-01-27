/* 
 * File:   itkMultipleModelMetricLocal.hxx
 * Author: rap58
 *
 * Created on 19 November 2014, 12:51
 */

#ifndef ITKMULTIPLEMODELMETRICLOCAL_HXX
#define	ITKMULTIPLEMODELMETRICLOCAL_HXX

#include "itkMultipleModelMetricLocal.h"
#include "corticalbone.h"

namespace itk
{
    // Constructor
    template< typename TTransform >
    MultipleModelMetric<TTransform>
    ::MultipleModelMetric() :  Superclass()   {
        m_sampledImageValues = ArrayType(1); m_sampledImageValues.Fill(nan("1"));
        m_sampledModelValues = ArrayType(1); m_sampledModelValues.Fill(nan("1"));
        m_sampledPositions = ArrayType(1); m_sampledPositions.Fill(nan("1"));

        m_sampledWeights = ArrayType(1); m_sampledWeights.Fill(nan("1"));

        m_fwhmModeSet = false;
        m_weightingMode = TransformType::kWeightingFunction1; m_dynamicWeightingMode = false;

        m_iterationCount = 0;
    }

    template< typename TTransform >
    typename MultipleModelMetric<TTransform>::MeasureType
    MultipleModelMetric<TTransform>
    ::GetValue( const ParametersType & parameters ) const {

        //std::cerr<<"Profile Size m: "<<this->m_ProfileObject->GetSizeModificationTime()<<", This Sample m: "<<this->m_ProfileSizeTimeStamp<<std::endl;
        if(!this->m_Profile|| !this->m_Transform ) {
            std::cerr<<"ERROR: m_FixedImage or m_ProfileObject or m_Transform has not been assigned."<<std::endl; // m_Interpolator included in m_Profile
        }

        // update TF with new parameters
        if(this->m_TransformTimeStamp < this->m_Transform->GetMTime() || this->m_Transform->GetParameters() != parameters || (m_sampledModelValues.GetSize() == 1 && m_sampledModelValues[0] == nan("1")) ) {
            this->m_Transform->SetParameters(parameters); // reset TF parameters, resample model

            if(m_fwhmModeSet) {
                this->SetFWHMCBDensity(parameters);
            }

            ResampleModelProfile();
            if(m_dynamicWeightingMode) { ResampleModelWeights(); }
            this->m_TransformTimeStamp.Modified();
        }

        // view variables
        //cout<<"MMM: xP="<<parameters[0]<<",xE="<<parameters[1]<<",yTB="<<parameters[2]<<",yST="<<parameters[3]<<",sig="<<parameters[4];//<<endl; ",RW="<<parameters[6]<< ",yCB="<<parameters[3]<<

        ParametersType parameterErrors = this->m_Transform->GetParameterErrors();

        unsigned int n = this->m_Profile->GetNumberOfSamples();
        unsigned int nError = parameterErrors.GetSize();

        MeasureType linearDifference;
        linearDifference.SetSize(n+nError);

        for(int i=0; i<parameters.GetSize(); i++) {
            if(std::isinf(parameters[i]) || std::isnan(parameters[i])) { // check for NaN or infinity
                linearDifference.Fill(std::numeric_limits<double>::max());
                //cout<<"\t,abs="<<std::numeric_limits<double>::max()<<"\t";//<<endl;
                return linearDifference;
            }
        }

        // calculate linear difference at each point
        //double absError=0;
        for(int i = 0; i < n; i++) {
            ScalarType elementError;
            if(isnan(this->m_sampledImageValues[i])) {
                elementError = 0; // if out of image bounds to not include in error considerations
            } else {
                //elementError = (this->m_sampledImageValues->GetValue(i)-m_sampledModelValues->GetValue(i)) * this->m_sampledWeights->GetValue(i);
                elementError = (this->m_sampledImageValues[i]-m_sampledModelValues[i]) * this->m_sampledWeights[i];
            }
            //ScalarType elementError = (this->m_sampledImageValues->GetValue(i)-m_sampledModelValues->GetValue(i)) * this->m_sampledWeights->GetValue(i);

            if(std::isinf(elementError) || std::isnan(elementError)) { // ensure a valid value
                elementError=std::numeric_limits<double>::max();
            }

            linearDifference.SetElement(i, elementError);
            //absError+=fabs(elementError);
        }

        //cerr<<"Parameters:"<<parameters<<", Parameter Errors:"<<parameterErrors<<", Error Sum:"<<absError<<endl;//", Linear Difference:"<<linearDifference<<endl;
        for(int i=0; i<nError; i++) {

            double parameterError = parameterErrors[i];
            if(std::isinf(parameterError) || std::isnan(parameterError)) { // ensure a valid value
                parameterError=std::numeric_limits<double>::max();
            }

            linearDifference.SetElement(n+i, parameterError);
            //absError+=fabs(parameterError);
        }
//        cerr<<"Parameters:"<<parameters<<", Parameter Errors"<<parameterErrors<<", Error Sum:"<<absError<<", Linear Difference:"<<linearDifference<<endl;
        m_iterationCount++;
        //this->m_LastTransformParameters = parameters;
        //cout<<"\t,abs="<<absError<<"\t";//this->GetAbsErrorSum(parameters)<<"\t";//<<endl;
        //cout<<"MER: xP="<<parameterErrors[0]<<",xE="<<parameterErrors[1]<<",Th="<<parameterErrors[2]<<",CB="<<parameterErrors[3]<<",TB="<<parameterErrors[4]<<",ST="<<parameterErrors[5]<<",sig="<<parameterErrors[6]<<endl;
        return linearDifference;
    }

    template< typename TTransform >
    void MultipleModelMetric<TTransform>
    ::Initialize() throw ( ExceptionObject ) {
        Superclass::Initialize();

        // ensure profile positions are up-to-date
        if(this->m_ProfileSizeTimeStamp < this->m_Profile->GetSizeModificationTime() || (m_sampledPositions.GetSize() == 1 && m_sampledPositions[0] == nan("1")) ) {
            m_sampledPositions = this->m_Profile->GetPositions(); // get resampeld profile
            this->m_ProfileSizeTimeStamp.Modified();
        }

        // ensure profile image samples are up-to-date
        if( this->m_ProfileTimeStamp < this->m_Profile->GetPositionModificationTime() || (m_sampledImageValues.GetSize() == 1 && m_sampledImageValues[0] == nan("1")) ) {
            m_sampledImageValues = this->m_Profile->GetValues();
            this->m_ProfileTimeStamp.Modified();
        }

        // ensure transform weights are up to date (requires transform parameters to have been set)
        ResampleModelWeights(); // always call as the tf may have changed

    }

    template< typename TTransform >
    void  MultipleModelMetric<TTransform>
    ::GetDerivative( const ParametersType &, DerivativeType & ) const {
        std::cout<<"Error - not implemented: 'LinearModelMetric::GetDerivative(...)"<<std::endl;
        std::cerr<<"Error - not implemented: 'LinearModelMetric::GetDerivative(...)"<<std::endl;
    }

    template< typename TTransform >
    void  MultipleModelMetric<TTransform>
    ::RestartIterationCount() const {
        m_iterationCount=0;
    }

    template< typename TTransform >
    typename MultipleModelMetric<TTransform>::ScalarType
    MultipleModelMetric<TTransform>
    ::GetIterationCount() const {
        return m_iterationCount;
    }

    template< typename TTransform >
    void  MultipleModelMetric<TTransform>
    ::ResampleModelProfile(void) const {

        int numberOfSamples = this->m_Profile->GetNumberOfSamples();

        // reset sample array
        m_sampledModelValues = ArrayType(numberOfSamples);

        // get sample parameters
        ScalarType increment = this->m_Profile->GetIncrement();
        ScalarType startPosition = this->m_Profile->GetStartXValue();

        ScalarType modelValue;

        // resample the image along the profile
        for (int i = 0; i < numberOfSamples; ++i)  {
            // get value
            modelValue = this->m_Transform->TransformPoint(startPosition + i * increment);
            m_sampledModelValues[i] = modelValue;
        }

    }

    template< typename TTransform >
    void  MultipleModelMetric<TTransform>
    ::ResampleModelWeights(void) const { // todo - consider efficiency

        double start, end; this->m_Profile->GetSampledExtents(start, end);
        this->m_Transform->setStart(start); this->m_Transform->setStop(end);

        // requires upto date profile and up to date transform
        int numberOfSamples = this->m_Profile->GetNumberOfSamples();

        // reset sample array
        m_sampledWeights = ArrayType(numberOfSamples);

        ScalarType startPosition = this->m_Profile->GetStartXValue();
        this->m_Transform->GetWeightingFunction(m_sampledWeights, startPosition, m_weightingMode, m_weightingScale);

    }

    template< typename TTransform >
    typename MultipleModelMetric<TTransform>::ArrayType
    MultipleModelMetric<TTransform>
    ::GetModelProfileSamples(void) const {

        // explictly check to see if the model values are up to date (i.e. sampled with the current TF parameters)
        if(this->m_TransformTimeStamp < this->m_Transform->GetMTime()) {

            ResampleModelProfile();
            if(m_dynamicWeightingMode) { ResampleModelWeights(); }
            this->m_TransformTimeStamp.Modified();
        }

        return this->m_sampledModelValues;
    }

    template< typename TTransform >
    typename MultipleModelMetric<TTransform>::ArrayType
    MultipleModelMetric<TTransform>
    ::GetUnblurredModelSamples(void) const {

        int numberOfSamples = this->m_Profile->GetNumberOfSamples();

        // create array to store values in
        ArrayType unblurredModel = ArrayType(numberOfSamples);

        // get sample parameters
        ScalarType increment = this->m_Profile->GetIncrement();

        // resample the image along the profile
        for (int i = 0; i < numberOfSamples; ++i)  {
            // get value
            ScalarType modelValue = this->m_Transform->GetUntransformedPoint(i * increment);
            unblurredModel[i] = modelValue;
        }
        return unblurredModel;
    }

    template< typename TTransform >
    typename MultipleModelMetric< TTransform >::ArrayType
    MultipleModelMetric< TTransform >
    ::GetWeightingArray(void) const {
        // explictly check to see if the model values are up to date (i.e. sampled with the current TF parameters)
        if(this->m_ProfileTimeStamp < this->m_Profile->GetPositionModificationTime()) {

            ResampleModelWeights();
        } else if(m_dynamicWeightingMode && this->m_TransformTimeStamp < this->m_Transform->GetMTime()) {
            ResampleModelProfile();
            ResampleModelWeights();
            this->m_TransformTimeStamp.Modified();
        }


        return this->m_sampledWeights;
    }

    template< typename TTransform >
    typename MultipleModelMetric<TTransform>::ScalarType
    MultipleModelMetric<TTransform>
    ::GetAbsErrorSum(ParametersType parameters) const {

        // explictly check to see if the model values are up to date (i.e. sampled with the current TF parameters)
        if(this->m_TransformTimeStamp < this->m_Transform->GetMTime()) {

            ResampleModelProfile();
            if(m_dynamicWeightingMode) { ResampleModelWeights(); }
            this->m_TransformTimeStamp.Modified();
        }

        ScalarType absErrorSum = 0;

        int startIndex, endIndex;
        this->m_Profile->GetSampledExtentIndices(startIndex, endIndex);
        for(int i=startIndex; i<=endIndex; ++i) {
            absErrorSum += fabs(m_sampledImageValues[i] -  m_sampledModelValues[i]);
        }

        return absErrorSum;
    }

    template< typename TTransform >
    typename MultipleModelMetric<TTransform>::ScalarType
    MultipleModelMetric<TTransform>
    ::GetAbsErrorMean(ParametersType parameters) const {
        int numberOfInboundsSamples = this->m_Profile->GetNumberOfInboundsSamples();
        return GetAbsErrorSum(parameters) / numberOfInboundsSamples;

    }

    template< typename TTransform >
    typename MultipleModelMetric<TTransform>::ScalarType
    MultipleModelMetric<TTransform>
    ::GetAbsErrorSumWeighted(ParametersType parameters) const {

        // explictly check to see if the model values are up to date (i.e. sampled with the current TF parameters)
        if(this->m_TransformTimeStamp < this->m_Transform->GetMTime()) {

            ResampleModelProfile();
            if(m_dynamicWeightingMode) { ResampleModelWeights(); }
            this->m_TransformTimeStamp.Modified();
        }

        ScalarType absErrorSum = 0;

        int startIndex, endIndex;
        this->m_Profile->GetSampledExtentIndices(startIndex, endIndex);
        for(int i=startIndex; i<=endIndex; ++i) {
            //absErrorSum += fabs(m_sampledImageValues->GetValue(i) -  m_sampledModelValues->GetValue(i)) * m_sampledWeights->GetValue(i);
            absErrorSum += fabs(m_sampledImageValues[i] -  m_sampledModelValues[i]) * m_sampledWeights[i];
        }

        return absErrorSum;
    }

    template< typename TTransform >
    typename MultipleModelMetric<TTransform>::ScalarType
    MultipleModelMetric<TTransform>
    ::GetSqrErrorMean(ParametersType parameters) const {
        int numberOfSamples = this->m_Profile->GetNumberOfSamples();
        return GetSqrErrorSum(parameters) / numberOfSamples;
    }

    template< typename TTransform >
    typename MultipleModelMetric<TTransform>::ScalarType
    MultipleModelMetric<TTransform>
    ::GetSqrErrorSum(ParametersType parameters) const {
        int numberOfSamples = this->m_Profile->GetNumberOfSamples();
        // explictly check to see if the model values are up to date (i.e. sampled with the current TF parameters)
        if(this->m_TransformTimeStamp < this->m_Transform->GetMTime()) {

            ResampleModelProfile();
            if(m_dynamicWeightingMode) { ResampleModelWeights(); }
            this->m_TransformTimeStamp.Modified();
        }

        ScalarType absErrorSum = 0;

        int startIndex, endIndex;
        this->m_Profile->GetSampledExtentIndices(startIndex, endIndex);
        for(int i=startIndex; i<=absErrorSum; ++i) {
            absErrorSum += pow(m_sampledImageValues[i] -  m_sampledModelValues[i], 2);
        }

        return absErrorSum;
    }

    template< typename TTransform >
    unsigned int  MultipleModelMetric<TTransform>
    ::GetNumberOfValues(void) const {
        if ( !this->m_Profile ) {
            itkExceptionMacro(<< "Moving point set has not been assigned");
            std::cerr<<"Error in MultipleModelMetric->GetNumberOfValues() m_ProfileObject not allocated"<<std::endl;
        }
        return this->m_Profile->GetNumberOfSamples()+this->m_Transform->GetNumberOfErrors();
    }

    template< typename TTransform >
    void  MultipleModelMetric<TTransform>
    ::TurnOnFWHMMode(void) const {
        m_fwhmModeSet = true;
    }

    template< typename TTransform >
    void  MultipleModelMetric<TTransform>
    ::TurnOffFWHMMode(void) const {
        m_fwhmModeSet = false;
    }

    template< typename TTransform >
    void  MultipleModelMetric< TTransform >
    ::SetDynamicWeighting(bool status) const {
        m_dynamicWeightingMode = status;
    }

    template< typename TTransform >
    bool  MultipleModelMetric< TTransform >
    ::GetDynamicWeighting(void) const {
        return m_dynamicWeightingMode;
    }

    template< typename TTransform >
    void  MultipleModelMetric<TTransform>
    ::SetFWHMCBDensity(ParametersType parameters) const { // TODO tidy this up

        //cout<<"FWHM: xP="<<parameters[0]<<", xE="<<parameters[1]<<", yTB="<<parameters[2]<<", yCB="<<parameters[3]<<" ,yST="<<parameters[4]<<", sigma="<<parameters[5]<<endl;

        ScalarType maxDensity = GetFWHMCBDensity();

        this->m_Transform->FixCorticalDensity(maxDensity); // assumes cbhas already been fixed - otherwise SetParams in the next line will be of the wrong size
        this->m_Transform->SetParameters(parameters); // re-run to ensure all parameters consistently defined

        //cout<<"FWHM: xP="<<parameters[0]<<", xE="<<parameters[1]<<", yTB="<<parameters[2]<<" ,yST="<<parameters[3]<<", sigma="<<parameters[4]<<", yCB="<<maxDensity;

        return; // validity bool
    }

    template< typename TTransform >
    typename MultipleModelMetric<TTransform>::ScalarType
    MultipleModelMetric<TTransform>
    ::GetFWHMCBDensity() const {

        ScalarType maxDensity = -DBL_MAX;
        if (this->m_Transform->IsValid()) {

            // get start and stop index
            ScalarType increment = this->m_Profile->GetIncrement();
            int indexStart = (int) floor(
                    (this->m_Transform->GetPeriostealEdgePosistion() - this->m_Profile->GetStartXValue()) / increment);
            int indexEnd = (int) ceil(
                    (this->m_Transform->GetEndostealEdgePosistion() - this->m_Profile->GetStartXValue()) / increment);

            // enforce bounds on the start and stop indicies
            int maxSampleIndex = this->m_Profile->GetNumberOfSamples() - 1;

            //cout << "Increment="<<increment<<"; Pedge="<<this->m_Transform->GetPeriostealEdgePosistion()<<", Eedge="<<this->m_Transform->GetEndostealEdgePosistion()<<endl;
            //cout << "Initial Index Start="<<indexStart<<", Index End="<<indexEnd<<endl;

            indexStart = (indexStart < 0) ? 0 : ((indexStart > maxSampleIndex) ? maxSampleIndex : indexStart);
            indexEnd = (indexEnd < 0) ? 0 : ((indexEnd > maxSampleIndex) ? maxSampleIndex : indexEnd);

            //cout << "Final   Index Start="<<indexStart<<", Index End="<<indexEnd<<endl;

            // loop through values in-between and select max updating the max sample value
            int maxIndex = -1;
            for (int i = indexStart; i <= indexEnd; i++) {

                if (maxDensity < m_sampledImageValues[i]) {
                    maxDensity = m_sampledImageValues[i];
                    maxIndex = i;
                }

            }

            //cout<<"start Index="<<indexStart<<", max Index"<<maxIndex<<", end Index="<<indexEnd;

            // if the max is on one side sample until a local maxima is reached or the end of the profile is reached
            if (maxIndex == indexStart) { // corner case indexStart is max and = 0. leave as max - the only reasonable/valid thing to do
                for (int i = indexStart - 1; i >= 0; i--) {
                    if (maxDensity < m_sampledImageValues[i]) {
                        maxDensity = m_sampledImageValues[i];
                        //maxIndex=i;
                    } else {
                        break;
                    }
                }
            } else if (maxIndex == indexEnd) {
                for (int i = indexEnd + 1; i <= maxSampleIndex; i++) {
                    if (maxDensity < m_sampledImageValues[i]) {
                        maxDensity = m_sampledImageValues[i];
                        //maxIndex=i;
                    } else {
                        break;
                    }
                }
            }
            //cout<<", final max Index="<<maxIndex<<endl;

        } else {
            int n = m_sampledImageValues.GetSize();
            for (int i = 0; i < n; i++) {
                if (maxDensity < m_sampledImageValues[i]) {
                    maxDensity = m_sampledImageValues[i];
                }
            }
        }

        return maxDensity;
    }

    template< typename TTransform >
    void  MultipleModelMetric<TTransform>
    ::SetWeightingMode(int selectedIndex, ScalarType scale) const {

        if(selectedIndex>=TransformType::kHueristicWeighting && selectedIndex<=TransformType::kDensityWeighting) {
            m_weightingMode=selectedIndex; m_weightingScale=scale;
        } else {
            m_weightingMode=-1; m_weightingScale=nan("1");
        }
    }

    template< typename TTransform >
    int MultipleModelMetric<TTransform>
    ::GetWeightingMode() const {

        return m_weightingMode;
    }

    template< typename TTransform >
    typename MultipleModelMetric< TTransform >::ScalarType
    MultipleModelMetric<TTransform>
    ::GetWeightingScale() const {
        return m_weightingScale;
    }

}

#endif	/* ITKMULTIPLEMODELMETRICLOCAL_HXX */

