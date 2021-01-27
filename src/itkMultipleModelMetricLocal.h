/* 
 * File:   itkMultipleModelMetricLocal.h
 * Author: rap58
 *
 * Created on 19 November 2014, 12:51
 */

#ifndef ITKMULTIPLEMODELMETRICLOCAL_H
#define	ITKMULTIPLEMODELMETRICLOCAL_H

#include <itkMultipleModelMetricBaseLocal.h>
#include <vtkDoubleArray.h>

namespace itk
{

    template <typename TTransform>
    class MultipleModelMetric : public itk::MultipleModelMetricBase<TTransform>
    {
        //  Software Guide : EndCodeSnippet
    public:
        typedef MultipleModelMetric  Self;
        typedef itk::MultipleModelMetricBase< TTransform > Superclass;
        typedef itk::SmartPointer<Self>           Pointer;
        typedef itk::SmartPointer<const Self>     ConstPointer;
        typedef typename Superclass::ParametersType ParametersType;
        typedef typename Superclass::DerivativeType DerivativeType;
        typedef typename Superclass::MeasureType    MeasureType; // measure type = Array<double>
        typedef typename Superclass::MovingSpatialObjectType    MovingSpatialObjectType;
        typedef typename Superclass::ScalarType      ScalarType;
        typedef itk::Array<ScalarType> ArrayType;
        typedef TTransform TransformType;

        itkNewMacro(Self);
        itkTypeMacro(MultipleModelMetricBase, ImageToSpatialObjectMetric);
        itkStaticConstMacro( ParametricSpaceDimension, unsigned int, 3 );


        //---------- Methods-----------//
        MultipleModelMetric();
        virtual ~MultipleModelMetric() {}
        MultipleModelMetric(const Self &) {}

        void GetDerivative( const ParametersType &, DerivativeType & ) const ITK_OVERRIDE;

        // calculates the fitness value - difference between model and image profile
        MeasureType    GetValue( const ParametersType & parameters ) const ITK_OVERRIDE;
        virtual void Initialize(void) throw ( ExceptionObject ) ITK_OVERRIDE; // Initialise the metric

        ArrayType GetModelProfileSamples(void) const;
        ArrayType GetUnblurredModelSamples(void) const;
        ArrayType GetWeightingArray(void) const;

        virtual ScalarType GetAbsErrorMean(ParametersType parameters) const ITK_OVERRIDE;
        virtual ScalarType GetAbsErrorSum(ParametersType parameters) const ITK_OVERRIDE;
        virtual ScalarType GetSqrErrorMean(ParametersType parameters) const ITK_OVERRIDE;
        virtual ScalarType GetSqrErrorSum(ParametersType parameters) const ITK_OVERRIDE;
        virtual ScalarType GetAbsErrorSumWeighted(ParametersType parameters) const ITK_OVERRIDE;

        void TurnOnFWHMMode(void) const;
        void TurnOffFWHMMode(void) const;
        bool GetFWHMMode() { return m_fwhmModeSet; };
        ScalarType GetFWHMCBDensity() const;

        void SetDynamicWeighting(bool status) const;
        bool GetDynamicWeighting(void) const;

        void SetWeightingMode(int selectionIndex, ScalarType scale) const;
        int GetWeightingMode() const;
        ScalarType GetWeightingScale() const;

        unsigned int GetNumberOfValues(void) const ITK_OVERRIDE;

        void ResampleModelWeights(void) const ITK_OVERRIDE;
        ScalarType GetImageValue(ScalarType positionAlongProfile) const;

        void RestartIterationCount(void) const;
        ScalarType GetIterationCount(void) const;

        //------------ overrides ------------//


    private:
        void ResampleModelProfile(void) const ITK_OVERRIDE;

        void SetFWHMCBDensity(ParametersType parameters) const;

        mutable ArrayType m_sampledImageValues;
        mutable ArrayType m_sampledModelValues;
        mutable ArrayType m_sampledPositions;
        mutable ArrayType m_sampledWeights;

        mutable bool m_fwhmModeSet;
        mutable bool m_dynamicWeightingMode;
        mutable int m_weightingMode;
        mutable ScalarType m_weightingScale;

        mutable ScalarType m_iterationCount;

    };

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultipleModelMetricLocal.hxx"
#endif

#endif	/* ITKMULTIPLEMODELMETRICLOCAL_H */
