/* 
 * File:   itkMultipleModelMetricLocal.h
 * Author: rap58
 *
 * Created on 19 November 2014, 12:16
 */

#ifndef ITKMULTIPLEMODELMETRICBASELOCAL_H
#define	ITKMULTIPLEMODELMETRICBASELOCAL_H

#include "itkMultipleValuedCostFunction.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkLinearInterpolateImageFunction.h"
#include "vnl/vnl_vector_fixed.h"
#include "itkTransform.h"
#include "itkModelTransformBaseLocal.h"
#include "itkProfileSpatialObjectLocal.h"

namespace itk
{
    /** \class ImageToSpatialObjectMetric
     * \brief Computes similarity between a moving spatial obejct
     *        and an Image to be registered
     *
     *  The ImageToSpatialObjectMetric is different from the rest of the
     *  registration framework in ITK regarding the interpretation of the Transform
     *  with respect to the Fixed and Moving objects. In most of the ITK
     *  registration framework, the Transform computed by the optimizer is the one
     *  that maps points from the space of the Fixed object into the space of the
     *  Moving object. This direction of the transform is the one that makes easier
     *  to resample the Moving object into the space of the Fixed object.
     *
     *  In the particular case of the ImageToSpatialObject registration, the
     *  Transform to be computed is the one mapping points from the SpatialObject
     *  into the Image, despite the fact that the SpatialObject is called the
     *  "Moving" object and the image is called the "Fixed" object. This change of
     *  reference system is the consequence of using this type of registration in
     *  applications that are based on Visualization. In the context of such
     *  visualizations it is simpler to think in terms of the Transform that can be
     *  used for displaying the SpatialObject in the appropriate position with
     *  respect to the image. Since this process does not involve resampling, but
     *  providing a Transform to a visualization routine, it is usually more
     *  natural to use the Transform that maps points from the SpatialObject space
     *  the image space.
     *
     *  A full discussion of the Transform directions in the ITK registration
     *  framework can be found in the ITK Software Guide.
     *
     * \ingroup ITKRegistrationCommon
     */

    template<typename TTransform>
    class MultipleModelMetricBase:    public MultipleValuedCostFunction
    {
    public:
        typedef MultipleModelMetricBase Self;
        typedef MultipleValuedCostFunction   Superclass;
        typedef SmartPointer< Self >       Pointer;
        typedef SmartPointer< const Self > ConstPointer;
        static constexpr unsigned int Dimension = 3;

        //----- Set Types ------//

        typedef ProfileSpatialObject<Dimension> MovingSpatialObjectType; //typedef TMovingSpatialObject MovingSpatialObjectType;
        itkStaticConstMacro(ObjectDimension, unsigned int, MovingSpatialObjectType::ObjectDimension);

        /**  Type of the Transform Base class */
        typedef TTransform TransformType;
        typedef typename TransformType::ScalarType      ScalarType;
        typedef typename TransformType::Pointer         TransformPointer;
        typedef typename TransformType::InputPointType  InputPointType;
        typedef typename TransformType::OutputPointType OutputPointType;
        typedef typename TransformType::ParametersType  TransformParametersType;
        typedef typename TransformType::JacobianType    TransformJacobianType;

        /** Typede of the vector type to return derivatives */
        typedef vnl_vector_fixed< double, itkGetStaticConstMacro(ObjectDimension) > VectorType;
        typedef std::list<double>     ListType;

        //------------ Superclass typedefs -----------//
        typedef Superclass::MeasureType MeasureType; // double array in SingleValuedCostFunction
        typedef Superclass::DerivativeType DerivativeType;

        // moving spatial image type defs
        typedef typename MovingSpatialObjectType::Pointer MovingSpatialObjectPointer;
        typedef typename MovingSpatialObjectType::ConstPointer MovingSpatialObjectConstPointer;

        //  ParametersType typedef - position in optimisation search space //
        //typedef Superclass::ParametersType ParametersType;
        typedef typename TransformType::ParametersType ParametersType;

        /** Run-time type information (and related methods). */
        itkTypeMacro(MultipleModelMetricBase, Object);


        //--------- Getters and Setters -----------//
        virtual void SetProfile( const MovingSpatialObjectType * spatialMovingObject);
        const MovingSpatialObjectType * GetProfile();

        virtual void SetMovingTransform( TransformType * spatialMovingObject);
        TransformType * GetMovingTransform();

        virtual ScalarType GetAbsErrorMean(ParametersType parameters) const = 0;
        virtual ScalarType GetAbsErrorSum(ParametersType parameters) const = 0;
        virtual ScalarType GetSqrErrorMean(ParametersType parameters) const = 0;
        virtual ScalarType GetSqrErrorSum(ParametersType parameters) const = 0;
        virtual ScalarType GetAbsErrorSumWeighted(ParametersType parameters) const = 0;

        unsigned int GetNumberOfParameters(void) const ITK_OVERRIDE;

        void GetValueAndDerivative(const ParametersType & parameters,  MeasureType & Value, DerivativeType  & Derivative) const;

        ModifiedTimeType GetMTime() const ITK_OVERRIDE;

        //----------- Misc -------------//
        virtual void Initialize(void) throw ( ExceptionObject ); // Initialise the metric

        //----------- Pure Virtual methods to implement --------------//
        virtual void GetDerivative( const ParametersType &, DerivativeType & ) const ITK_OVERRIDE = 0;
        virtual MeasureType GetValue( const ParametersType & parameters ) const ITK_OVERRIDE = 0;

        virtual void ResampleModelWeights(void) const = 0;

        // Get the last transformation parameters visited by the optimiser. Overload superclass method
        //itkGetConstReferenceMacro(LastTransformParameters, ParametersType);
        virtual unsigned int GetNumberOfValues(void) const ITK_OVERRIDE = 0;

    protected:

        virtual void ResampleModelProfile(void) const = 0;

        MultipleModelMetricBase();
        virtual ~MultipleModelMetricBase() {}
        MultipleModelMetricBase(const Self &) {}
        void operator=(const Self &) {}
        virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

        MeasureType              m_MatchMeasure;
        DerivativeType           m_MatchMeasureDerivatives;

        mutable TransformPointer m_Transform;
        mutable typename MovingSpatialObjectType::ConstPointer m_Profile;

        // Time stamp - use to track updating of profile sampling
        mutable TimeStamp         m_ProfileTimeStamp;
        mutable TimeStamp         m_ProfileSizeTimeStamp;
        mutable TimeStamp         m_TransformTimeStamp;
    };
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultipleModelMetricBaseLocal.hxx"
#endif

#endif	/* ITKMULTIPLEMODELMETRICLOCAL_H */
