/* 
 * File:   thresholderTransform.h
 * Author: rap58
 *
 * Created on 17 June 2015, 10:43
 */

#ifndef LINEARTRANSFORM_H
#define	LINEARTRANSFORM_H

// includes
#include "itkProcessObject.h"
#include <vtkSmartPointer.h>

#include <itkImage.h>
#include <itkProfileSpatialObjectLocal.h>
#include <itkLinearInterpolateImageFunction.h>
#include <LinearRegressionCalculator.h>
//#include <itkBSplineInterpolateImageFunction.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <vtkTable.h>
#include <itkOptimizerParameters.h>
#include <vtkIntArray.h>
#include <iostream>

// defines

namespace itk {

    class LinearTransform : public ProcessObject {

    public:

        typedef LinearTransform Self;
        typedef ProcessObject                                           Superclass;
        typedef SmartPointer< Self >                                    Pointer;
        typedef SmartPointer< const Self >                              ConstPointer;

        // standard class macros
        itkTypeMacro(LinearTransform, ProcessObject);

        typedef signed short                                            PixelType;
        typedef double                                                  ScalarType;
        const static unsigned int                                       Dimension = 3;
        typedef itk::TimeStamp                                          TimeType;
        typedef vtkSmartPointer<vtkDoubleArray>                         ArrayType;
        typedef itk::Array<double>                                      itkArrayType;
        typedef itk::Array<int>                                         itkIntArrayType;
        typedef vtkSmartPointer<vtkIntArray>                            IntArrayType;
        typedef vtkSmartPointer<vtkTable>                               TableType;
        typedef itk::OptimizerParameters< ScalarType >                  ParametersType;
        typedef itk::Array2D<ScalarType> itkArray2DType;
        typedef LinearRegressionCalculator::RegressionLine              RegressionLineType;
        typedef LinearRegressionCalculator::ProfileProperties           ProfilePropertiesType;


        typedef enum {
            kVLowPercent=0, // 0.01%
            kLowPercent,    // 0.1%
            kMidPercent,    //   1%
            kHighPercent,   //   5%
            kManualThreshold,
            kMedianMidpoint,
            kMedianManual,
        } ThresholdMode; // threshold generation percentages/modes

        LinearTransform();

        //void SetProfileProperties(ProfileProperties &properties);

        void SetProfileValues(itkArrayType values);
        void SetProfilePositions(itkArrayType positions);
        void SetMultipleProfileValues(itkArray2DType multipleValues);

        static void SetProfilePositionInfo(int startIndex, int endIndex, int numberOfProfiles, double startPosition);
        static void SetProfileMaxPositionInfo(int maxIndex, ScalarType maxValue, itkArrayType maxValues, itkIntArrayType maxIndices);
        static void SetProfileSampleInfo(ScalarType increment, int numberOfSamples, ScalarType length, int edgeIndex);

        void SetClassifierLevels(ScalarType stDensity, ScalarType cbDensity, ScalarType threshold);
        void SetMode(int index);
        int  GetMode();
        virtual std::string GetModeName() = 0;
        void ClearThresholds(); // clear the CB, ST, threshold values and the percentage selection

        virtual bool Threshold() = 0;
        bool EstimateThresholds();

        bool IsValid() const;
        virtual bool IsUpToDate() const;
        bool AreThresholdsSet() const;
        bool ReadyToProcess(bool thresholdsRequired = true) const;

        itkArrayType GetProfilePositions() const;
        itkArrayType GetProfileValues() const;
        itkArray2DType GetMultipleProfileValues() const;
        static ProfilePropertiesType GetProfileProperties();
        bool GetClassifierThresholds(ScalarType &stDensity, ScalarType &cbDensity, ScalarType &Threshold);
        int GetClassifierLevelIndex();
        void SetClassifierLevelIndex(int percentIndex);
        std::string GetClassifierLevelName();

        virtual itkArrayType GetDisplayValues() const = 0;
        virtual unsigned int GetNumberOfDisplayValues(bool includeImportedOptions=false) const;
        virtual std::string GetDisplayName(int index) const = 0;
        virtual std::string GetDisplayNameShort(int index) const = 0;
        virtual itkArray2DType GetDisplayRanges(bool includeImportedOptions=false) const = 0;

        virtual itkArrayType GetParameterValues() const = 0;
        virtual unsigned int GetNumberOfParameters() const;
        virtual std::string GetParameterName(int index) const = 0;
        virtual std::string GetParameterNameShort(int index) const = 0;

        unsigned int GetNumberOfSamples() const;

        // getters for display values
        virtual TableType GetProfileDisplayTable() const = 0;
        virtual TableType GetModelDisplayTable() const = 0;
        virtual TableType GetModelDisplayTable(ArrayType parametersIn, double offset) const = 0;
        virtual TableType GetClassifierTable() const { return NULL; };
        ScalarType GetErrorMean() const;


    protected:

        static LinearRegressionCalculator m_linearRegressionCalculator;

        ScalarType LookupPercentDensity(int index);
        ScalarType LookupPValue(int index);
        ScalarType LookupPercentThreshold(int index);

        virtual void CalculateValidity() = 0;
        virtual ScalarType CalculateAbsMeanError() const = 0;

        virtual void CalculateTissueRegressionLines() = 0;

        virtual void ClearResults() = 0;

        void ClearThresholdValues(); // only clear the CB, ST and threshold values

        unsigned int m_numberOfParameters;
        unsigned int m_numberOfDisplayValues;

        itkArrayType m_positions;
        itkArrayType m_meanProfiles;
        itkArray2DType m_multipleProfiles;
        int m_modeIndex;

        ScalarType m_threshold, m_thresholdPercentage, m_thresholdPValue;
        int m_thresholdIndex;

        ScalarType m_ST, m_CB;
        RegressionLineType m_STLine, m_CBLine, m_ECLine, m_TBLine;

        bool m_transformUpToDate, m_thresholdsSet, m_valid; // todo state to indicate that thresholds have been calculated

        mutable TimeType  m_ProcessedTimeStamp;

    private:

    };
}

#endif	/* LINEARTRANSFORM_H */

