/* 
 * File:   itkModelLinearMetricLocal.hxx
 * Author: rap58
 *
 * Created on 11 November 2014, 15:21
 */

#ifndef ITKMODELLINEARMETRICLOCAL_HXX
#define	ITKMODELLINEARMETRICLOCAL_HXX

#include "itkSingleModelMetricLocal.h"

namespace itk
{
    // Constructor
    template< typename TTransform >
    SingleModelMetric< TTransform >
    ::SingleModelMetric() :  Superclass()   {
        m_sampledImageValues = ArrayType(1); m_sampledImageValues.Fill(nan("1"));
        m_sampledModelValues = ArrayType(1); m_sampledModelValues.Fill(nan("1"));
        m_sampledPositions = ArrayType(1); m_sampledPositions.Fill(nan("1"));

        m_sampledWeights = ArrayType(1); m_sampledWeights.Fill(nan("1"));

        m_fwhmModeSet = false;
        m_dynamicWeightingMode = false; m_weightingMode = TransformType::kWeightingFunction1;

        m_iterationCount = 0;
    }

    template< typename TTransform >
    typename SingleModelMetric< TTransform >::MeasureType
    SingleModelMetric< TTransform >
    ::GetValue( const ParametersType & parameters ) const {

        //std::cerr<<"Profile Size m: "<<this->m_ProfileObject->GetSizeModificationTime()<<", This Sample m: "<<this->m_ProfileSizeTimeStamp<<std::endl;
        if(!this->m_Profile|| !this->m_Transform) { // interpolator in Profile
            std::cerr<<"ERROR: m_FixedImage or m_ProfileObject or m_Transform or m_Interpolator has not been assigned."<<std::endl;
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

        ScalarType totalSquaredDifference = 0; // todo squared error vs absolute error

        // a constraint on the parameter space
        ScalarType parameterError = this->m_Transform->GetParameterErrorSum();

        unsigned int n = this->m_Profile->GetNumberOfSamples();

        for(int i = 0; i < n; i++) {

            if(!isnan(this->m_sampledImageValues[i])) { // ignore out of bounds values
                ScalarType elementError = (std::pow(this->m_sampledImageValues[i]-m_sampledModelValues[i], 2)) * this->m_sampledWeights[i];
                totalSquaredDifference += elementError + parameterError; // parameter error is independent of weighting
            }

        }

        //this->m_LastTransformParameters = parameters;

        //cerr<<"Parameters: xP="<<parameters[0]<<", xE="<<parameters[1]<<", yTB="<<parameters[2]<<", yCB="<<parameters[3]<<", yST="<<parameters[4]<<", Sigma="<<parameters[5];
        //if(parameters.GetSize()==7) { cerr<<"ramp width="<<parameters[6]; }
        //cerr<<", TotalSqrErr="<<totalSquaredDifference<<endl;

        if(std::isnan(totalSquaredDifference) || std::isinf(totalSquaredDifference)) {
            totalSquaredDifference = std::numeric_limits<double>::max();
        }
        m_iterationCount++;

        //this->m_LastTransformParameters = parameters;

        return totalSquaredDifference;
    }

    template< typename TTransform >
    void  SingleModelMetric< TTransform >
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
    void  SingleModelMetric< TTransform >
    ::GetDerivative( const ParametersType &, DerivativeType & ) const {
        std::cout<<"Error - not implemented: 'LinearModelMetric::GetDerivative(...)"<<std::endl;
        std::cerr<<"Error - not implemented: 'LinearModelMetric::GetDerivative(...)"<<std::endl;
    }

    template< typename TTransform >
    void  SingleModelMetric<TTransform>
    ::RestartIterationCount() const {
        m_iterationCount=0;
    }

    template< typename TTransform >
    typename SingleModelMetric<TTransform>::ScalarType
    SingleModelMetric<TTransform>
    ::GetIterationCount() const {
        return m_iterationCount;
    }

    template< typename TTransform >
    void  SingleModelMetric< TTransform >
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
    void  SingleModelMetric< TTransform >
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
    typename SingleModelMetric< TTransform >::ArrayType
    SingleModelMetric< TTransform >
    ::GetModelProfileSamples(void) const {

        // explictly check to see if the model values are up to date (i.e. sampled with the current TF parameters)
        if(this->m_TransformTimeStamp < this->m_Transform->GetMTime()) {

            ResampleModelProfile();
            this->m_TransformTimeStamp.Modified();
        }

        return this->m_sampledModelValues;
    }

    template< typename TTransform >
    typename SingleModelMetric< TTransform >::ArrayType
    SingleModelMetric< TTransform >
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
    typename SingleModelMetric< TTransform >::ArrayType
    SingleModelMetric< TTransform >
    ::GetWeightingArray(void) const {
        // explictly check to see if the model values are up to date (i.e. sampled with the current TF parameters)
        if(this->m_ProfileTimeStamp < this->m_Profile->GetPositionModificationTime()) {

            ResampleModelWeights();
            this->m_ProfileTimeStamp.Modified();
        }

        return this->m_sampledWeights;
    }

    template< typename TTransform >
    typename SingleModelMetric< TTransform >::ScalarType
    SingleModelMetric< TTransform >
    ::GetAbsErrorSum(ParametersType parameters) const {

        // explictly check to see if the model values are up to date (i.e. sampled with the current TF parameters)
        if(this->m_TransformTimeStamp < this->m_Transform->GetMTime()) {

            ResampleModelProfile();
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
    typename SingleModelMetric< TTransform >::ScalarType
    SingleModelMetric< TTransform >
    ::GetAbsErrorMean(ParametersType parameters) const {
        int numberOfInboundsSamples = this->m_Profile->GetNumberOfInboundsSamples();
        return GetAbsErrorSum(parameters) / numberOfInboundsSamples;

    }

    template< typename TTransform >
    typename SingleModelMetric< TTransform >::ScalarType
    SingleModelMetric< TTransform >
    ::GetAbsErrorSumWeighted(ParametersType parameters) const {

        // explictly check to see if the model values are up to date (i.e. sampled with the current TF parameters)
        if(this->m_TransformTimeStamp < this->m_Transform->GetMTime()) {

            ResampleModelProfile();
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
    typename SingleModelMetric< TTransform >::ScalarType
    SingleModelMetric< TTransform >
    ::GetSqrErrorMean(ParametersType parameters) const {
        int numberOfInboundsSamples = this->m_Profile->GetNumberOfInboundsSamples();
        return GetSqrErrorSum(parameters) / numberOfInboundsSamples;
    }

    template< typename TTransform >
    typename SingleModelMetric< TTransform >::ScalarType
    SingleModelMetric< TTransform >
    ::GetSqrErrorSum(ParametersType parameters) const {
        int numberOfSamples = this->m_Profile->GetNumberOfSamples();
        // explictly check to see if the model values are up to date (i.e. sampled with the current TF parameters)
        if(this->m_TransformTimeStamp < this->m_Transform->GetMTime()) {

            ResampleModelProfile();
            this->m_TransformTimeStamp.Modified();
        }

        ScalarType absErrorSum = 0;

        int startIndex, endIndex;
        this->m_Profile->GetSampledExtentIndices(startIndex, endIndex);
        for(int i=startIndex; i<=endIndex; ++i) {
            absErrorSum += pow(m_sampledImageValues[i] -  m_sampledModelValues[i], 2); // * m_sampledWeights[i];
        }

        return absErrorSum;
    }

    template< typename TTransform >
    void  SingleModelMetric< TTransform >
    ::TurnOnFWHMMode(void) const {
        m_fwhmModeSet = true;
    }

    template< typename TTransform >
    void  SingleModelMetric< TTransform >
    ::TurnOffFWHMMode(void) const {
        m_fwhmModeSet = false;
    }

    template< typename TTransform >
    void  SingleModelMetric< TTransform >
    ::SetDynamicWeighting(bool status) const {
        m_dynamicWeightingMode = status;
    }

    template< typename TTransform >
    bool  SingleModelMetric< TTransform >
    ::GetDynamicWeighting(void) const {
        return m_dynamicWeightingMode;
    }

    template< typename TTransform >
    void  SingleModelMetric< TTransform >
    ::SetFWHMCBDensity(ParametersType parameters) const { // TODO tidy this up

        //cout<<"FWHM: xP="<<parameters[0]<<", xE="<<parameters[1]<<", yTB="<<parameters[2]<<", yCB="<<parameters[3]<<" ,yST="<<parameters[4]<<", sigma="<<parameters[5]<<endl;

        ScalarType maxDensity = GetFWHMCBDensity();

        this->m_Transform->FixCorticalDensity(maxDensity); // assumes cbhas already been fixed - otherwise SetParams in the next line will be of the wrong size
        this->m_Transform->SetParameters(parameters); // re-run to ensure all parameters consistently defined

        //cout<<"FWHM: xP="<<parameters[0]<<", xE="<<parameters[1]<<", yTB="<<parameters[2]<<" ,yST="<<parameters[3]<<", sigma="<<parameters[4]<<", yCB="<<maxDensity;

        return; // validity bool
    }

    template< typename TTransform >
    typename SingleModelMetric< TTransform >::ScalarType
    SingleModelMetric< TTransform >
    ::GetFWHMCBDensity() const {

        ScalarType maxDensity = -DBL_MAX;
        if(this->m_Transform->IsValid()) {

            // get start and stop index
            ScalarType increment = this->m_Profile->GetIncrement();
            int indexStart = (int)floor((this->m_Transform->GetPeriostealEdgePosistion()-this->m_Profile->GetStartXValue())/increment);
            int indexEnd = (int)ceil((this->m_Transform->GetEndostealEdgePosistion()-this->m_Profile->GetStartXValue())/increment);

            // enforce bounds on the start and stop indicies
            int maxSampleIndex = this->m_Profile->GetNumberOfSamples()-1;

            //cout << "Increment="<<increment<<"; Pedge="<<this->m_Transform->GetPeriostealEdgePosistion()<<", Eedge="<<this->m_Transform->GetEndostealEdgePosistion()<<endl;
            //cout << "Initial Index Start="<<indexStart<<", Index End="<<indexEnd<<endl;

            indexStart = (indexStart<0) ? 0 : ( (indexStart>maxSampleIndex) ? maxSampleIndex : indexStart);
            indexEnd = (indexEnd<0) ? 0 : ( (indexEnd>maxSampleIndex) ? maxSampleIndex : indexEnd);

            //cout << "Final   Index Start="<<indexStart<<", Index End="<<indexEnd<<endl;

            // loop through values in-between and select max updating the max sample value
            int maxIndex = -1;
            for(int i=indexStart; i<=indexEnd; i++) {

                if(maxDensity<m_sampledImageValues[i]) {
                    maxDensity = m_sampledImageValues[i];
                    maxIndex = i;
                }

            }

            //cout<<"start Index="<<indexStart<<", max Index"<<maxIndex<<", end Index="<<indexEnd;

            // if the max is on one side sample until a local maxima is reached or the end of the profile is reached
            if(maxIndex==indexStart) { // corner case indexStart is max and = 0. leave as max - the only reasonable/valid thing to do
                for(int i=indexStart-1; i>=0; i--) {
                    if(maxDensity<m_sampledImageValues[i]) {
                        maxDensity = m_sampledImageValues[i];
                        //maxIndex=i;
                    } else {
                        break;
                    }
                }
            } else if(maxIndex==indexEnd) {
                for(int i=indexEnd+1; i<=maxSampleIndex; i++) {
                    if(maxDensity<m_sampledImageValues[i]) {
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
                if(maxDensity < m_sampledImageValues[i]) {
                    maxDensity = m_sampledImageValues[i];
                }
            }
        }

        return maxDensity;
    }

    template< typename TTransform >
    void  SingleModelMetric<TTransform>
    ::SetWeightingMode(int selectedIndex, ScalarType scale) const {

        if(selectedIndex>=TransformType::kHueristicWeighting && selectedIndex<=TransformType::kDensityWeighting) {
            m_weightingMode=selectedIndex; m_weightingScale=scale;
        } else {
            m_weightingMode=-1; m_weightingScale=scale;
        }
    }

    template< typename TTransform >
    int SingleModelMetric<TTransform>
    ::GetWeightingMode() const {
        return m_weightingMode;
    }

    template< typename TTransform >
    typename SingleModelMetric< TTransform >::ScalarType
    SingleModelMetric<TTransform>
    ::GetWeightingScale() const {
        return m_weightingScale;
    }

}

#endif	/* ITKMODELLINEARMETRICLOCAL_HXX */

