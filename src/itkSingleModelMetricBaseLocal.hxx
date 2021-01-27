/* 
 * File:   itkModelMetricBaseLocal.hxx
 * Author: rap58
 *
 * Created on 11 November 2014, 15:18
 */

#ifndef ITKMODELMETRICBASELOCAL_HXX
#define	ITKMODELMETRICBASELOCAL_HXX

#include "itkSingleModelMetricBaseLocal.h"

namespace itk
{
    // Constructor
    template<typename TTransform>
    SingleModelMetricBase< TTransform >
    ::SingleModelMetricBase(): m_MatchMeasure(0) {
        m_Profile             = ITK_NULLPTR; // has to be provided by the user.
        m_Transform           = ITK_NULLPTR; // has to be provided by the user.

        // set time stamps zero
        m_ProfileTimeStamp.Modified();
        m_ProfileSizeTimeStamp.Modified();
        m_TransformTimeStamp.Modified();
    }

    // Initialize
    template<typename TTransform>
    void SingleModelMetricBase< TTransform >
    ::Initialize(void) throw ( ExceptionObject ) {
        if ( !m_Transform )    {
            itkExceptionMacro(<< "Transform is not present");
        }
        if ( !m_Profile )    { // include interpolator
            itkExceptionMacro(<< "MovingSpatialObject is not present");
        }

        // If there are any observers on the metric, call them to give the
        // user code a chance to set parameters on the metric
        this->InvokeEvent( InitializeEvent() );
    }

    // PrintSelf
    template<typename TTransform>
    void SingleModelMetricBase< TTransform >
    ::PrintSelf(std::ostream & os, Indent indent) const {
        Superclass::PrintSelf(os, indent);
        os << indent << "Moving Spatial Object: " << m_Profile.GetPointer()  << std::endl;
        os << indent << "Transform:    " << m_Transform.GetPointer()    << std::endl;
        //os << indent << "Last Transform parameters = " << m_LastTransformParameters << std::endl;
    }

    //------------------- Getters and Setters ------------------//
    template<typename TTransform>
    void SingleModelMetricBase< TTransform >::
    SetProfile( const MovingSpatialObjectType * profileObject) {
        this->m_Profile = profileObject;
    }

    template<typename TTransform>
    const typename SingleModelMetricBase< TTransform >::MovingSpatialObjectType *
    SingleModelMetricBase< TTransform >::
    GetProfile() {
        return this->m_Profile;
    }

    template<typename TTransform>
    void SingleModelMetricBase< TTransform >::
    SetMovingTransform( TransformType * transform) {
        this->m_Transform = transform;
    }

    template<typename TTransform>
    typename SingleModelMetricBase< TTransform >::TransformType *
    SingleModelMetricBase< TTransform >::
    GetMovingTransform() {
        return this->m_Transform;
    }

    template<typename TTransform>
    unsigned int  SingleModelMetricBase< TTransform >::
    GetNumberOfParameters() const {
        return this->m_Transform->GetNumberOfParameters();
    }

    template<typename TTransform>
    ModifiedTimeType  SingleModelMetricBase< TTransform >::
    GetMTime() const {
        // return biggest time = most recently changed
        TimeStamp mtime = this->m_ProfileSizeTimeStamp;
        if(this->m_ProfileTimeStamp.GetMTime() > mtime.GetMTime()) {
            mtime = this->m_ProfileTimeStamp;
        }
        if(this->m_TransformTimeStamp.GetMTime() > mtime.GetMTime()) {
            mtime =  this->m_ProfileTimeStamp;
        }

        return mtime.GetMTime();
    }

    template<typename TTransform>
    void SingleModelMetricBase< TTransform >::
    GetValueAndDerivative( const ParametersType & parameters, MeasureType & value, DerivativeType  & derivative ) const {
        value = this->GetValue(parameters);
        this->GetDerivative(parameters, derivative);
    }


} // end namespace itk

#endif	/* ITKMODELMETRICBASELOCAL_HXX */
