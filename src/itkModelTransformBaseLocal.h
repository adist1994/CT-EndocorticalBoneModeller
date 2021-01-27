/* 
 * File:   itkModelTransformBaseLocal2.h
 * Author: rap58
 *
 * Created on 10 November 2014, 17:28
 */

#ifndef ITKMODELTRANSFORMBASELOCAL2_H
#define	ITKMODELTRANSFORMBASELOCAL2_H

#include "itkMacro.h"
#include "itkMatrix.h"
#include "itkTransform.h"
#include "itkArray2D.h"

#include <iostream>

namespace itk
{
/** \class ModelTransformBase
 * \brief Matrix and Offset transformation of a vector space (e.g. space coordinates)
 *
 *
 * \tparam ScalarT            The type to be used for scalar numeric values.  Either
 *    float or double.
 *
 * \tparam NInputDimensions   The number of dimensions of the input vector space.
 *
 * \tparam NOutputDimensions  The number of dimensions of the output vector space.
 *
 * \ingroup ITKTransform
 */

    template < typename TScalar = double,  unsigned int NSpaceDimensions = 1 >
// Number of dimensions in the output space
    class ModelTransformBase :  public Transform<TScalar, NSpaceDimensions, NSpaceDimensions>
    {
    public:
        /** Standard typedefs   */
        typedef ModelTransformBase Self;
        typedef Transform<TScalar,NSpaceDimensions,NSpaceDimensions> Superclass;

        typedef SmartPointer<Self>       Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        /** Dimension of the domain space. */
        itkStaticConstMacro(InputSpaceDimension, unsigned int, NSpaceDimensions);
        itkStaticConstMacro(OutputSpaceDimension, unsigned int, NSpaceDimensions);
        //itkStaticConstMacro( ParametersDimension, unsigned int, NParameterDimensions );

        // Superclass definitions
        typedef typename Superclass::ParametersType      ParametersType;
        typedef typename Superclass::NumberOfParametersType    NumberOfParametersType;
        //typedef typename Superclass::ParametersValueType ParametersValueType;

        typedef typename Superclass::JacobianType JacobianType;
        typedef typename Superclass::ScalarType ScalarType;



        // Other definitions
        typedef Array2D <ScalarType> PointArrayType;

        typedef ScalarType ProfilePosistionType;
        typedef Array<ScalarType> ProfilePosistionArrayType;

        /** Standard vector type for this class   */
        typedef Vector<ScalarType, itkGetStaticConstMacro(InputSpaceDimension)>  InputVectorType;
        typedef Vector<ScalarType, itkGetStaticConstMacro(OutputSpaceDimension)> OutputVectorType;
        //typedef typename OutputVectorType::ValueType OutputVectorValueType;

        /** Standard vnl_vector type for this class   */
        //typedef vnl_vector_fixed<TScalar, itkGetStaticConstMacro(InputSpaceDimension)> InputVnlVectorType;
        //typedef vnl_vector_fixed<TScalar, itkGetStaticConstMacro(OutputSpaceDimension)> OutputVnlVectorType;

        /** Standard coordinate point type for this class   */
        typedef Point<ScalarType, itkGetStaticConstMacro(InputSpaceDimension)> InputPointType;
        //typedef typename InputPointType::ValueType InputPointValueType;
        typedef Point<ScalarType, itkGetStaticConstMacro(OutputSpaceDimension)> OutputPointType;
        //typedef typename OutputPointType::ValueType OutputPointValueType;

        typedef enum { // display
            kBlank = -1,
            kTotalThicknessDisp =0, // returned by the TF's
            kCorticalThicknessDisp=1,
            kEndoThicknessDisp=2,
            kPeriPositionDisp=3,
            kEndoPositionDisp=4,
            kEndoCBPositionDisp=5,
            kEndoTBPositionDisp=6,
            kDenseCorticalBoneDensityDisp =7,
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

        typedef enum { // display
            kHueristicWeighting=0,
            kWeightingFunction1=1, // returned by the TF's
            kWeightingFunction2=2,
            kWeightingFunction3=3,
            kWeightingFunction4=4,
            kDensityWeighting=5,
        } WeightingFunctionType;


        /** Run-time type information (and related methods).   */
        itkTypeMacro(ModelTransformBase, Transform);

        //--------------- Methods ---------------//
        virtual void SetIdentity(void);

        void SetParameters(const ParametersType & parameters) ITK_OVERRIDE;
        const ParametersType & GetParameters(void) const ITK_OVERRIDE;

        virtual void SetFixedParameters(const ParametersType &) ITK_OVERRIDE;
        virtual const ParametersType & GetFixedParameters(void) const ITK_OVERRIDE;

        virtual void FixParameter(const ScalarType value, const unsigned int index);
        virtual void FreeParameter(const unsigned int index);

        virtual bool isSigmaFixed() const = 0;
        virtual bool isCorticalDensityFixed() const = 0;

        virtual void FixCorticalDensity(const ScalarType value) = 0;


        const ParametersType & GetCombinedParameters(void) const;
        void SetCombinedParameters(const ParametersType & parameters);


        virtual NumberOfParametersType GetNumberOfParameters(void) const ITK_OVERRIDE {
            return this->m_numberOfFlexibleParameters;
        }

        virtual NumberOfParametersType GetNumberOfCombinedParameters(void) const {
            return m_combinedParameters.GetSize();
        }

        virtual bool setStart(ScalarType start);
        virtual bool setStop(ScalarType stop);
        virtual void SetProfileEdgeProperties(ScalarType ratio, ScalarType maxOffset);
        virtual ScalarType getStart() {return m_start;};
        virtual ScalarType getStop() {return m_stop;};
        virtual TimeStamp GetProfileModificationTime() const;

        virtual void GetWeightingFunction(ProfilePosistionArrayType & array, ScalarType startPosition, int selectionIndex, ScalarType scale) const;
        virtual void CalculateWeightingFunction1(ProfilePosistionArrayType &array, ScalarType startPosition, ScalarType scale) const = 0;
        virtual void CalculateWeightingFunction2(ProfilePosistionArrayType & array, ScalarType startPosition, ScalarType scale) const = 0;
        virtual void CalculateWeightingFunction3(ProfilePosistionArrayType & array, ScalarType startPosition, ScalarType scale) const = 0;
        virtual void CalculateWeightingFunction4(ProfilePosistionArrayType & array, ScalarType startPosition, ScalarType scale) const = 0;
        virtual void CalculateWeightingDensity(ProfilePosistionArrayType & array, ScalarType startPosition, ScalarType scale) const = 0;
        virtual void CalculateWeightingHeuristic(ProfilePosistionArrayType & array, ScalarType startPosition, ScalarType scale) const = 0;

        virtual ScalarType CalculateSigmaAdjustedCorticalDensity(ScalarType peakSampledValue, ScalarType sigma) = 0;

        // display values
        //virtual ScalarType GetDisplayValues() const = 0;
        virtual ScalarType GetTotalThickness() const = 0;
        virtual ScalarType GetCorticalThickness() const = 0;
        virtual ScalarType GetEndostealThickness() const = 0;
        virtual ScalarType GetCorticalDensity() const = 0;
        virtual ScalarType GetTrabeculaeDensity() const = 0;
        virtual ScalarType GetSoftTissueDensity() const = 0;
        virtual ScalarType GetSurfaceThickness() const = 0;
        virtual ScalarType GetPeriostealEdgePosistion() const = 0;
        virtual ScalarType GetEndostealEdgePosistion() const = 0;
        virtual ScalarType GetSigma() const = 0;

        virtual ScalarType GetMeanCorticalDensity() const = 0;

        virtual int GetCorticalDensityFlexIndex() const = 0;
        virtual int GetPeriostealEdgeFlexIndex() const = 0;

        virtual ParametersType GetDisplayValues() const = 0;

        virtual ScalarType GetParameterErrorSum() const = 0;
        virtual ParametersType GetParameterErrors() const = 0;
        virtual NumberOfParametersType GetNumberOfErrors(void) const = 0;

        virtual PointArrayType GetModelPoints(ScalarType xStartPosition) const = 0;
        virtual PointArrayType GetModelPoints(ScalarType xStartPosition, ParametersType parameters) const = 0;

        virtual void SetMaxDensity(ScalarType maxDensity);
        virtual ScalarType GetMaxDensity();

        virtual bool IsValid() const;

        virtual void SetProfileLength(ScalarType length);

        OutputPointType       TransformPoint(const InputPointType & point) const ITK_OVERRIDE;
        ScalarType TransformPoint(const ScalarType position) const;

        ScalarType GetUntransformedPoint(const ScalarType position) const;

        ModifiedTimeType GetMTime() const ITK_OVERRIDE;

        using Superclass::TransformVector;
        OutputVectorType      TransformVector(const InputVectorType & vector) const ITK_OVERRIDE;
        ProfilePosistionArrayType TransformVector(const ProfilePosistionArrayType & positions) const;

        /** Compute the Jacobian of the transformation
         *
         * This method computes the Jacobian matrix of the transformation.
         * given point or vector, returning the transformed point or
         * vector. The rank of the Jacobian will also indicate if the transform
         * is invertible at this point.
         * Get local Jacobian for the given point
         * \c j will sized properly as needed.
         */
        virtual void ComputeJacobianWithRespectToParameters(const InputPointType  & x, JacobianType & j) const ITK_OVERRIDE;

        /** Get the jacobian with respect to position. This simply returns
         * the current Matrix. jac will be resized as needed, but it's
         * more efficient if it's already properly sized. */
        virtual void ComputeJacobianWithRespectToPosition(const InputPointType  & x, JacobianType & jac) const ITK_OVERRIDE;

        /** Get the jacobian with respect to position. This simply returns
         * the inverse of the current Matrix. jac will be resized as needed, but it's
         * more efficient if it's already properly sized. */
        virtual void ComputeInverseJacobianWithRespectToPosition(const InputPointType  & x, JacobianType & jac) const ITK_OVERRIDE;

        /** Indicates that this transform is linear. That is, given two
         * points P and Q, and scalar coefficients a and b, then
         *
         *           T( a*P + b*Q ) = a * T(P) + b * T(Q)
         */
        virtual bool IsLinear() const ITK_OVERRIDE {
            return true;
        }



    protected:

        virtual ScalarType PeakTransformValue() const = 0;

        virtual ScalarType TransformValue(ScalarType position) const = 0;
        virtual ScalarType UntransformedValue(ScalarType position) const = 0;

        virtual void UpdateInternalParameters() = 0;

        ModelTransformBase(unsigned int combinedNumberOfParameters);
        ModelTransformBase();

        /** Destroy an ModelTransformBase object */
        virtual ~ModelTransformBase();

        /** Print contents of an ModelTransformBase */
        void PrintSelf(std::ostream & s, Indent indent) const ITK_OVERRIDE;

        TimeStamp         m_TimeStamp;

        ScalarType m_MaxDensity; // set from image - or based on image
        static constexpr ScalarType m_MaxErrorMultiplier = 700; // as exp(709) gives just below inf

        bool m_valid;

        ScalarType m_length; // indicates the max value of the parameter x values
        ScalarType m_start, m_stop;
        ScalarType m_maxProfileOffset, m_profileEdgeRatio;

        ScalarType m_numberOfFlexibleParameters;

        ParametersType m_combinedParameterMapping;
        ParametersType m_combinedParameters;

        unsigned int m_NumberOfDisplays;

        ParametersType m_ParameterErrorScales;

        TimeStamp  m_ProfileTimeStamp;

    private:

        ModelTransformBase(const Self & other);
        const Self & operator=(const Self &);



    }; // class ModelTransformBase
}  // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkModelTransformBaseLocal.hxx"
#endif

#endif	/* ITKMODELTRANSFORMBASELOCAL2_H */

