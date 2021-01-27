/* 
 * File:   itkThreeTierRectangularTransformLocal.h
 * Author: rap58
 *
 * Created on 10 November 2014, 17:31
 */

#ifndef ITKTHREETIERRECTANGULARTRANSFORMLOCAL_H
#define	ITKTHREETIERRECTANGULARTRANSFORMLOCAL_H

#include <iostream>
#include "itkModelTransformBaseLocal.h"

namespace itk
{

    template< typename TScalar = double > // Data type for scalars (float or double)
    class ThreeTierRectangularTransform : public itk::ModelTransformBase< TScalar, 1 >  {
    public:

        //------ dimension macros -----------//
        itkStaticConstMacro(SpaceDimension, unsigned int, 1);
        itkStaticConstMacro(ParametersDimension, unsigned int, 6);
        itkStaticConstMacro(OutputDimension, unsigned int, 7);

        /** Standard class typedefs. */
        typedef ThreeTierRectangularTransform                           Self;
        typedef itk::ModelTransformBase< TScalar, SpaceDimension > Superclass;
        typedef itk::SmartPointer< Self >                       Pointer;
        typedef itk::SmartPointer< const Self >                 ConstPointer;

        // Run-time / Smart pointer Macros
        itkNewMacro(Self);
        itkTypeMacro(ThreeTierRectangularTransform, ModelTransformBase);

        typedef enum { // param
            kPeriostealEdgeParam = 0,
            kEnodstealEdgeParam,
            kTrabeculaeDensityParam,
            kCorticalDensityParam,
            kSoftTissueDensityParam,
            kSigmaParam,
        } ParameterEnum;

        /** Parameters type. */
        typedef typename Superclass::ScalarType                 ScalarType;
        typedef typename Superclass::ParametersType             ParametersType;
        typedef typename Superclass::NumberOfParametersType     NumberOfParametersType;
        typedef typename Superclass::PointArrayType             PointArrayType;
        typedef typename Superclass::ProfilePosistionArrayType  ArrayType;

        /** Jacobian type. */
        typedef typename Superclass::JacobianType JacobianType;

        typedef typename Superclass::WeightingFunctionType WeightingFunctionType;

        // / Standard vector type for this class
        typedef typename Superclass::InputVectorType       InputVectorType;
        typedef typename Superclass::OutputVectorType      OutputVectorType;
        //typedef typename Superclass::OutputVectorValueType OutputVectorValueType;


        // / Standard coordinate point type for this class
        typedef typename Superclass::InputPointType  InputPointType;
        typedef typename Superclass::OutputPointType OutputPointType;

        void FixCorticalDensity(const ScalarType value) ITK_OVERRIDE;

        bool isSigmaFixed() const ITK_OVERRIDE;
        bool isCorticalDensityFixed() const ITK_OVERRIDE;

        void SetMaxDensity(ScalarType maxDensity) ITK_OVERRIDE;
        void SetProfileLength(ScalarType length) ITK_OVERRIDE;

        ScalarType GetTotalThickness() const ITK_OVERRIDE;
        ScalarType GetCorticalThickness() const ITK_OVERRIDE;
        ScalarType GetEndostealThickness() const ITK_OVERRIDE;
        ScalarType GetCorticalDensity() const ITK_OVERRIDE;
        ScalarType GetTrabeculaeDensity() const ITK_OVERRIDE;
        ScalarType GetSoftTissueDensity() const ITK_OVERRIDE;
        ScalarType GetSurfaceThickness() const ITK_OVERRIDE;
        ScalarType GetPeriostealEdgePosistion() const ITK_OVERRIDE;
        ScalarType GetEndostealEdgePosistion() const ITK_OVERRIDE;
        virtual ScalarType GetSigma() const ITK_OVERRIDE;

        ScalarType GetMeanCorticalDensity() const ITK_OVERRIDE;

        int GetCorticalDensityFlexIndex() const ITK_OVERRIDE;
        int GetPeriostealEdgeFlexIndex() const ITK_OVERRIDE;

        ParametersType GetDisplayValues() const ITK_OVERRIDE;
        //ParametersType GetDisplayValues() const ITK_OVERRIDE;


        ParametersType GetParameterErrors() const ITK_OVERRIDE;
        NumberOfParametersType GetNumberOfErrors(void) const ITK_OVERRIDE;
        ScalarType GetParameterErrorSum() const ITK_OVERRIDE;//ScalarType multiplicationFactor) const;

        PointArrayType GetModelPoints(ScalarType xStartPosition) const ITK_OVERRIDE;
        PointArrayType GetModelPoints(ScalarType xStartPosition, ParametersType parameters) const ITK_OVERRIDE;

        void CalculateWeightingFunction1(ArrayType &array, ScalarType startPosition, ScalarType scale) const ITK_OVERRIDE;
        void CalculateWeightingFunction2(ArrayType & array, ScalarType startPosition, ScalarType scale) const ITK_OVERRIDE;
        void CalculateWeightingFunction3(ArrayType &array, ScalarType startPosition, ScalarType scale) const ITK_OVERRIDE;
        void CalculateWeightingFunction4(ArrayType & array, ScalarType startPosition, ScalarType scale) const ITK_OVERRIDE;
        void CalculateWeightingDensity(ArrayType & array, ScalarType startPosition, ScalarType scale) const ITK_OVERRIDE;
        void CalculateWeightingHeuristic(ArrayType & array, ScalarType startPosition, ScalarType scale) const ITK_OVERRIDE;

        ScalarType CalculateSigmaAdjustedCorticalDensity(ScalarType peakSampledValue, ScalarType globalSigma) ITK_OVERRIDE;

        /// Compute the Jacobian Matrix of the transformation at one point, allowing for thread-safety. //
        virtual void ComputeJacobianWithRespectToParameters( const InputPointType  & p, JacobianType & jacobian) const ITK_OVERRIDE;


        // Reset the parameters to create and identity transform.
        virtual void SetIdentity(void) ITK_OVERRIDE;

        //ModelTransformCategoryType GetTransformCategory() const;

    protected:
        //ThreeTierRectangularTransform(unsigned int outputSpaceDimension, unsigned int parametersDimension);
        //ThreeTierRectangularTransform(unsigned int parametersDimension);
        ThreeTierRectangularTransform();

        ~ThreeTierRectangularTransform();

        /**
          * Print contents of an ThreeTierRectangularTransform
          */
        void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

    private:

        ScalarType TransformValue(TScalar position) const ITK_OVERRIDE;
        ScalarType UntransformedValue(ScalarType position) const ITK_OVERRIDE;

        ScalarType PeakTransformValue() const ITK_OVERRIDE;

        ScalarType PeriostealError(ScalarType multiplicationFactor) const;
        ScalarType EndostealError(ScalarType multiplicationFactor) const;
        ScalarType ThicknessError(ScalarType multiplicationFactor) const;
        ScalarType CBMaxError(ScalarType multiplicationFactor) const;
        ScalarType CBTBDiffError(ScalarType multiplicationFactor) const;
        ScalarType CBSTDiffError(ScalarType multiplicationFactor) const;
        ScalarType SigmaError(ScalarType multiplicationFactor) const;

        void UpdateInternalParameters() ITK_OVERRIDE;


        ThreeTierRectangularTransform(const Self &); // purposely not implemented
        void operator=(const Self &);   // purposely not implemented

        // profile goes in [trabecular] to out[soft tissue]. y0 = tb, y1 = cb, y2 = st, x0 = endosteal, x1 = perosteal
        ScalarType m_ST, m_CB, m_TB, m_thickness, m_centre, m_sigma; // actual parameter values
        ScalarType m_STdiff, m_TBdiff, m_P, m_E; // values used to calculate blur curve


    }; // class ThreeTierRectangularTransform

}  // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkThreeTierRectangularTransformLocal.hxx"
#endif

#endif	/* ITKTHREETIERRECTANGULARTRANSFORMLOCAL_H */

