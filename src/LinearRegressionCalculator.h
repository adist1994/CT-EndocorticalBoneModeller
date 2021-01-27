//
// Created by rap58 on 19/10/15.
//

#ifndef MYPROJECT_LINEARREGRESSIONCALCULATOR_H
#define MYPROJECT_LINEARREGRESSIONCALCULATOR_H

// includes
#include <vtkDoubleArray.h>
#include <itkArray.h>
#include <itkArray2D.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>

class LinearRegressionCalculator {

public:

    typedef LinearRegressionCalculator Self;


    typedef double                                                  ScalarType;
    const static unsigned int                                       Dimension = 3;
    typedef vtkSmartPointer<vtkDoubleArray>                         ArrayType;
    typedef itk::Array2D<double>                                    itkArray2DType;
    typedef itk::Array<double>                                      itkArrayType;
    typedef itk::Array<int>                                         itkIntArrayType;
    typedef vtkSmartPointer<vtkIntArray>                            IntArrayType;

    struct RegressionLine {
        ScalarType meanX, meanY;
        ScalarType stdX, stdY;
        ScalarType varX, varY;
        ScalarType covXY, corrXY;
        ScalarType slope, intercept;
        ScalarType rSquared;
        int n;

        RegressionLine() {
            meanX = meanY = nan("1");
            stdX = stdY = nan("1");
            varX = varY = nan("1");
            covXY = corrXY = nan("1");
            slope = intercept = nan("1");
            rSquared = nan("1");
            n=-1;
        }
    } ;

    struct ProfileProperties {
        ScalarType increment, startPosition, length, maxValue;
        int numberOfSamples, numberOfProfiles;
        int startIndex, endIndex, maxIndex, edgeIndex;
        itkArrayType maxValues; itkIntArrayType maxIndices;
        ProfileProperties() {
            increment = startPosition = length = maxValue = nan("1");
            numberOfSamples = numberOfProfiles = -1;
            startIndex = endIndex = maxIndex = edgeIndex = -1;
            maxValues = itk::Array<double>(1); maxValues.Fill(nan("1"));
            maxIndices = itk::Array<int>(1); maxIndices.Fill(-1);
        }
    };

    LinearRegressionCalculator();

    void SetProfilePositionInfo(int startIndex, int endIndex, int numberOfProfiles, double startPosition);
    void SetProfileMaxPositionInfo(int maxIndex, ScalarType maxValue, itkArrayType maxValues, itkIntArrayType maxIndices);
    void SetProfileSampleInfo(ScalarType increment, int numberOfSamples, ScalarType length, int edgeIndex);

    int GetProfileStartIndex();
    int GetProfileEndIndex();
    int GetProfileNumber();
    double GetProfileStartPosition();
    int GetProfileMaxIndex();
    ScalarType GetProfileMaxValue();
    itkArrayType GetProfileMaxValues();
    itkIntArrayType GetProfileMaxIndices();
    ScalarType GetProfileIncrement();
    int GetProfileSampleNumber();
    ScalarType GetProfileLength();
    int GetProfileEdgeIndex();
    ProfileProperties GetProfileProperties();

    // static processing methods
    static RegressionLine CalculateLinearRegression(itkArrayType array1, itkArrayType array2);

    // processing methods
    ScalarType CalculateArrayMedian(itkArrayType array, int startIndex, int endIndex);
    ScalarType CalculateArrayMean(itkArrayType array, int startIndex, int endIndex);
    ScalarType CalculateArraySTD(itkArrayType array, int startIndex, int endIndex, ScalarType mean=nan("1"));
    ScalarType CalculateArrayCovariance(itkArrayType array1, itkArrayType array2, int startIndex, int endIndex, ScalarType mean1=nan("1"), ScalarType mean2=nan("1"));
    ScalarType CalculateArrayCorrelation(itkArrayType array1, itkArrayType array2, int startIndex, int endIndex, ScalarType mean1=nan("1"), ScalarType mean2=nan("1"), ScalarType std1=nan("1"), ScalarType std2=nan("1"));
    ScalarType CalculateRSquared(itkArrayType arrayX, itkArrayType arrayY, int startIndex, int endIndex, ScalarType meanY, ScalarType slope, ScalarType intercept);

    RegressionLine CalculateLinearRegression(itkArrayType array1, itkArrayType array2, int startIndex, int endIndex);
    RegressionLine CalculateAverageRegression(itkArrayType array1, itkArrayType array2, int startIndex, int endIndex);
    RegressionLine CalculateDemingRegression(itkArrayType array1, itkArrayType array2, int startIndex, int endIndex);

    RegressionLine CalculateLinearRegressionMultiple(itkArrayType arrayX, itkArray2DType arraysY, int startIndex, int endIndex);

    ScalarType LinearlyWeightIndex(itkArrayType array, int lowerIndex, ScalarType threshold);

    int CalculateThresholdCrossingIndex(itkArrayType array, int startIndex, int endIndex, ScalarType threshold);
    ScalarType CalculateLinearValue(RegressionLine line, ScalarType position) const;
    ScalarType CalculateWelchsTTest(RegressionLine line1, RegressionLine line2);


    private:

    ProfileProperties m_profileDetails;

};

#endif //MYPROJECT_LINEARREGRESSIONCALCULATOR_H
