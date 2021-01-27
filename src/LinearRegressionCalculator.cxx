//
// Created by rap58 on 19/10/15.
//

#include <vtkSortDataArray.h>
#include "LinearRegressionCalculator.h"

//----------- public ----------------//
LinearRegressionCalculator::LinearRegressionCalculator() {

    m_profileDetails = ProfileProperties();

}

void LinearRegressionCalculator::SetProfilePositionInfo(int startIndex, int endIndex, int numberOfProfiles, double startPosition) {
    m_profileDetails.startIndex = startIndex; m_profileDetails.endIndex = endIndex; m_profileDetails.numberOfProfiles = numberOfProfiles;
    m_profileDetails.startPosition = startPosition;
}

void LinearRegressionCalculator::SetProfileMaxPositionInfo(int maxIndex, ScalarType maxValue, itkArrayType maxValues, itkIntArrayType maxIndices) {
    m_profileDetails.maxIndex = maxIndex; m_profileDetails.maxValue = maxValue; m_profileDetails.maxValues = maxValues; m_profileDetails.maxIndices = maxIndices;
}

void LinearRegressionCalculator::SetProfileSampleInfo(ScalarType increment, int numberOfSamples, ScalarType length, int edgeIndex) {
    m_profileDetails.increment = increment; m_profileDetails.numberOfSamples =  numberOfSamples;
    m_profileDetails.length = length; m_profileDetails.edgeIndex = edgeIndex;
}

int LinearRegressionCalculator::GetProfileStartIndex() {
    return m_profileDetails.startIndex;
}

int LinearRegressionCalculator::GetProfileEndIndex() {
    return m_profileDetails.endIndex;
}

int LinearRegressionCalculator::GetProfileNumber() {
    return m_profileDetails.numberOfProfiles;
}

double LinearRegressionCalculator::GetProfileStartPosition() {
    return m_profileDetails.startPosition;
}

int LinearRegressionCalculator::GetProfileMaxIndex() {
    return m_profileDetails.maxIndex;
}

LinearRegressionCalculator::ScalarType LinearRegressionCalculator::GetProfileMaxValue() {
    return m_profileDetails.maxValue;
}

LinearRegressionCalculator::itkArrayType LinearRegressionCalculator::GetProfileMaxValues() {
    return m_profileDetails.maxValues;
}

LinearRegressionCalculator::itkIntArrayType LinearRegressionCalculator::GetProfileMaxIndices() {
    return m_profileDetails.maxIndices;
}

LinearRegressionCalculator::ScalarType LinearRegressionCalculator::GetProfileIncrement() {
    return m_profileDetails.increment;
}

int LinearRegressionCalculator::GetProfileSampleNumber() {
    return m_profileDetails.numberOfSamples;
}

LinearRegressionCalculator::ScalarType LinearRegressionCalculator::GetProfileLength() {
    return m_profileDetails.length;
}

int LinearRegressionCalculator::GetProfileEdgeIndex() {
    return m_profileDetails.edgeIndex;
}

LinearRegressionCalculator::ProfileProperties LinearRegressionCalculator::GetProfileProperties() {
    return m_profileDetails;
}


//----------- private ---------------//
LinearRegressionCalculator::ScalarType LinearRegressionCalculator::CalculateArrayMedian(itkArrayType array, int startIndex, int endIndex) {
    int n = endIndex - startIndex + 1;

    vtkSmartPointer<vtkDoubleArray> sortArray = vtkSmartPointer<vtkDoubleArray>::New();
    sortArray->SetNumberOfValues(n); // todo - consider a itk only solution

    for(int i=0; i<n; i++) {
        sortArray->SetValue(i,array[i+startIndex]);
    }

    vtkSortDataArray::Sort(sortArray);

    return 0.5 *(sortArray->GetValue(floor((n-1)/2.0)) + sortArray->GetValue(floor(n/2.0)));
}

LinearRegressionCalculator::ScalarType LinearRegressionCalculator::CalculateArrayMean(itkArrayType array, int startIndex, int endIndex) {

    int lowIndex = (endIndex > startIndex) ? startIndex : endIndex;
    int highIndex = (endIndex > startIndex) ? endIndex : startIndex;

    lowIndex = (lowIndex > m_profileDetails.startIndex) ? lowIndex : m_profileDetails.startIndex;
    highIndex = (highIndex < m_profileDetails.endIndex) ? highIndex : m_profileDetails.endIndex;

    int n = highIndex - lowIndex + 1;
    //cerr<<"n="<<n<<endl;
    ScalarType meanEstimate;
    for(int i=0; i<n; i++) {
        meanEstimate += array[i+lowIndex];
        //cerr<<"i="<<i<<"; value="<<array->GetValue(i)<<endl;
    }
    //cerr<<endl<<"mean="<<meanEstimate / ScalarType(n)<<endl;
    return meanEstimate / ScalarType(n);
}

LinearRegressionCalculator::ScalarType LinearRegressionCalculator::CalculateArraySTD(itkArrayType array, int startIndex, int endIndex, ScalarType mean) {

    ScalarType variance = CalculateArrayCovariance(array, array, startIndex, endIndex, mean, mean);

    return sqrt(variance);

}

LinearRegressionCalculator::ScalarType LinearRegressionCalculator::CalculateArrayCovariance(itkArrayType array1, itkArrayType array2, int startIndex, int endIndex, ScalarType mean1, ScalarType mean2) {
    if(isnan(mean1)) {
        mean1 = CalculateArrayMean(array1, startIndex, endIndex);
    }
    if(isnan(mean2)) {
        mean2 = CalculateArrayMean(array2, startIndex, endIndex);
    }

    int n = endIndex - startIndex + 1;

    ScalarType productSum=0;
    for(int i=0; i<n; i++) {
        productSum += array1[i+startIndex] * array2[i+startIndex];
    }
    ScalarType covariance = productSum / n - mean1 * mean2;

    return covariance;
}

LinearRegressionCalculator::ScalarType LinearRegressionCalculator::CalculateArrayCorrelation(itkArrayType array1, itkArrayType array2, int startIndex, int endIndex, ScalarType mean1, ScalarType mean2, ScalarType std1, ScalarType std2) {

    if(isnan(mean1)) {
        mean1 = CalculateArrayMean(array1, startIndex, endIndex);
    }
    if(isnan(mean2)) {
        mean2 = CalculateArrayMean(array2, startIndex, endIndex);
    }

    if(isnan(std1)) {
        std1 = CalculateArraySTD(array1, startIndex, endIndex);
    }
    if(isnan(std2)) {
        std2 = CalculateArraySTD(array2, startIndex, endIndex);
    }

    ScalarType covariance = CalculateArrayCovariance(array1, array2, startIndex, endIndex, mean1, mean2);
    ScalarType correlation = covariance / (std1 * std2);

    return correlation;
}

LinearRegressionCalculator::ScalarType LinearRegressionCalculator::CalculateRSquared(itkArrayType arrayX, itkArrayType arrayY, int startIndex, int endIndex, ScalarType meanY, ScalarType slope, ScalarType intercept) {
    int n = endIndex - startIndex + 1;

    ScalarType numeratorSum=0, denominatorSum=0;
    for(int i=0; i<n; i++) {
        ScalarType x = arrayX[i+startIndex];
        ScalarType y = arrayY[i+startIndex];
        ScalarType f = arrayX[i+startIndex] * slope + intercept;

        numeratorSum += (y - f) * (y - f);
        denominatorSum += (y - meanY) * (y - meanY);
    }

    ScalarType rSquared = 1 - numeratorSum / denominatorSum;

    return rSquared;
}

LinearRegressionCalculator::RegressionLine LinearRegressionCalculator::CalculateLinearRegression(itkArrayType arrayX, itkArrayType arrayY, int startIndex, int endIndex) {
    RegressionLine regressionLine = RegressionLine();

    if(startIndex >=0 && startIndex <= endIndex && endIndex < m_profileDetails.numberOfSamples) { // only attempt to calculate for valid indices

        int lowIndex = (endIndex > startIndex) ? startIndex : endIndex; lowIndex = (lowIndex > m_profileDetails.startIndex) ? lowIndex : m_profileDetails.startIndex;
        int highIndex = (endIndex > startIndex) ? endIndex : startIndex; highIndex = (highIndex < m_profileDetails.endIndex) ? highIndex : m_profileDetails.endIndex;
        regressionLine.n = highIndex - lowIndex + 1;


        // if coincident start & stop indices - nan slope and intercept, known meanX and meany
        regressionLine.meanX = CalculateArrayMean(arrayX, startIndex, endIndex);
        regressionLine.meanY = CalculateArrayMean(arrayY, startIndex, endIndex);

        regressionLine.stdX = CalculateArraySTD(arrayX, startIndex, endIndex, regressionLine.meanX);
        regressionLine.stdY = CalculateArraySTD(arrayY, startIndex, endIndex, regressionLine.meanY);

        regressionLine.corrXY = CalculateArrayCorrelation(arrayX, arrayY, startIndex, endIndex, regressionLine.meanX, regressionLine.meanY, regressionLine.stdX, regressionLine.stdY);

        regressionLine.slope = regressionLine.corrXY * regressionLine.stdY / regressionLine.stdX;
        regressionLine.intercept = regressionLine.meanY - regressionLine.slope * regressionLine.meanX;

        regressionLine.rSquared = CalculateRSquared(arrayX, arrayY, startIndex, endIndex, regressionLine.meanY, regressionLine.slope, regressionLine.intercept);
    }

    return regressionLine;
}

LinearRegressionCalculator::RegressionLine LinearRegressionCalculator::CalculateAverageRegression(itkArrayType arrayX, itkArrayType arrayY, int startIndex, int endIndex) {
    RegressionLine regressionLine = RegressionLine();

    if(startIndex >=0 && startIndex <= endIndex && endIndex < m_profileDetails.numberOfSamples) { // only attempt to calculate for valid indices

        int lowIndex = (endIndex > startIndex) ? startIndex : endIndex; lowIndex = (lowIndex > m_profileDetails.startIndex) ? lowIndex : m_profileDetails.startIndex;
        int highIndex = (endIndex > startIndex) ? endIndex : startIndex; highIndex = (highIndex < m_profileDetails.endIndex) ? highIndex : m_profileDetails.endIndex;
        regressionLine.n = highIndex - lowIndex + 1;

        regressionLine.meanX = CalculateArrayMean(arrayX, startIndex, endIndex);
        regressionLine.meanY = CalculateArrayMean(arrayY, startIndex, endIndex);

        regressionLine.stdX = CalculateArraySTD(arrayX, startIndex, endIndex, regressionLine.meanX);
        regressionLine.stdY = CalculateArraySTD(arrayY, startIndex, endIndex, regressionLine.meanY);

        regressionLine.covXY = CalculateArrayCovariance(arrayX, arrayY, startIndex, endIndex, regressionLine.meanX, regressionLine.meanY);

        regressionLine.corrXY = CalculateArrayCorrelation(arrayX, arrayY, startIndex, endIndex, regressionLine.meanX, regressionLine.meanY, regressionLine.stdX, regressionLine.stdY);

        ScalarType slopeX = regressionLine.corrXY * regressionLine.stdX / regressionLine.stdY;
        ScalarType slopeY = regressionLine.corrXY * regressionLine.stdY / regressionLine.stdX;

        ScalarType interceptX = regressionLine.meanX - slopeX * regressionLine.meanY;
        ScalarType interceptY = regressionLine.meanY - slopeY * regressionLine.meanX;


        regressionLine.slope = (slopeX + slopeY) / 2;
        regressionLine.intercept = (interceptX + interceptY) / 2;
    }
    return regressionLine;
}

LinearRegressionCalculator::RegressionLine LinearRegressionCalculator::CalculateDemingRegression(itkArrayType arrayX, itkArrayType arrayY, int startIndex, int endIndex) {
    RegressionLine regressionLine = RegressionLine();

    if(startIndex >=0 && startIndex <= endIndex && endIndex < m_profileDetails.numberOfSamples) { // only attempt to calculate for valid indices

        int lowIndex = (endIndex > startIndex) ? startIndex : endIndex; lowIndex = (lowIndex > m_profileDetails.startIndex) ? lowIndex : m_profileDetails.startIndex;
        int highIndex = (endIndex > startIndex) ? endIndex : startIndex; highIndex = (highIndex < m_profileDetails.endIndex) ? highIndex : m_profileDetails.endIndex;
        regressionLine.n = highIndex - lowIndex + 1;

        ScalarType sigma = 1000.0;

        regressionLine.meanX = CalculateArrayMean(arrayX, startIndex, endIndex);
        regressionLine.meanY = CalculateArrayMean(arrayY, startIndex, endIndex);

        regressionLine.stdX = CalculateArraySTD(arrayX, startIndex, endIndex, regressionLine.meanX);
        regressionLine.stdY = CalculateArraySTD(arrayY, startIndex, endIndex, regressionLine.meanY);

        regressionLine.covXY = CalculateArrayCovariance(arrayX, arrayY, startIndex, endIndex, regressionLine.meanX, regressionLine.meanY);
        regressionLine.varX = CalculateArrayCovariance(arrayX, arrayX, startIndex, endIndex, regressionLine.meanX, regressionLine.meanX);
        regressionLine.varY = CalculateArrayCovariance(arrayY, arrayY, startIndex, endIndex, regressionLine.meanY, regressionLine.meanY);


        regressionLine.corrXY = CalculateArrayCorrelation(arrayX, arrayY, startIndex, endIndex, regressionLine.meanX, regressionLine.meanY, regressionLine.stdX, regressionLine.stdY);

        regressionLine.slope = (regressionLine.varX - sigma * regressionLine.varX + sqrt( (regressionLine.varY - sigma * regressionLine.varX) * (regressionLine.varY - sigma * regressionLine.varX) + 4 * sigma * regressionLine.covXY * regressionLine.covXY) )/ (2 * regressionLine.covXY);
        regressionLine.intercept = regressionLine.meanY - regressionLine.slope * regressionLine.meanX;
    }
    return regressionLine;
}

LinearRegressionCalculator::RegressionLine LinearRegressionCalculator::CalculateLinearRegressionMultiple(itkArrayType arrayX, itkArray2DType arraysY, int startIndex, int endIndex) {
    RegressionLine regressionLine = RegressionLine();

    if(startIndex>=0 && endIndex>=0 && startIndex<m_profileDetails.numberOfSamples && endIndex<m_profileDetails.numberOfSamples) { // only attempt to calculate for valid indices
        int lowIndex = (endIndex > startIndex) ? startIndex : endIndex;
        int highIndex = (endIndex > startIndex) ? endIndex : startIndex;

        lowIndex = (lowIndex >= m_profileDetails.startIndex) ? lowIndex : m_profileDetails.startIndex;
        highIndex = (highIndex <= m_profileDetails.endIndex) ? highIndex : m_profileDetails.endIndex;

        int n = highIndex - lowIndex + 1, m = m_profileDetails.numberOfProfiles;

        ScalarType sumX = 0, sumY = 0, sumXY = 0, sumXX = 0, sumYY = 0;

        for(int i=0; i<n; i++) {
            for(int j=0; j<m; j++) {
                ScalarType x = arrayX[i+lowIndex];
                ScalarType y = arraysY[i+lowIndex][j];

                sumX += x; sumY += y; sumXY += x*y; sumXX += x*x; sumYY += y*y;
            }
        }

        ScalarType meanX = sumX / (n*m), meanY = sumY / (n*m);
        ScalarType stdX = sqrt(sumXX / (n*m) - meanX * meanX);
        ScalarType stdY = sqrt(sumYY / (n*m) - meanY * meanY);

        regressionLine.n = n*m;
        regressionLine.meanX = meanX; regressionLine.stdX = stdX;
        regressionLine.meanY = meanY; regressionLine.stdY = stdY;

        if(lowIndex==highIndex) { // if coincident start & stop indices - nan slope and intercept, known meanX and meany
            return regressionLine;
        }

        ScalarType covXY = sumXY / (n*m) - meanX * meanY, corrXY = covXY / (stdX * stdY);
        ScalarType slope = corrXY * stdY / stdX, intercept = meanY - slope * meanX;


        regressionLine.corrXY = corrXY; regressionLine.covXY = covXY;
        regressionLine.slope = slope; regressionLine.intercept = intercept;

        // calculate 'goodness of fit'
        ScalarType rSum = 0;//, rDenominator = 0;
        for(int i=0; i<n; i++) {
            for(int j=0; j<m; j++) {
                ScalarType y = arraysY[i+lowIndex][j];
                ScalarType f = arrayX[i+lowIndex] * slope + intercept;

                rSum += (y - f) * (y - f);
            }
        }
        regressionLine.rSquared = 1 - rSum / (stdY*stdY*n*m);
    }
    return regressionLine;
}

LinearRegressionCalculator::ScalarType LinearRegressionCalculator::LinearlyWeightIndex(itkArrayType array, int lowerIndex, ScalarType threshold) {

    // lower = smaller index, upper = one greater index

    int upperIndex = lowerIndex+1;
    //cerr<<"linearly weighted index options =["<<lowerIndex<<","<<upperIndex<<"];";


    // check index limits
    if(upperIndex < m_profileDetails.startIndex || lowerIndex > m_profileDetails.endIndex) { // out of bounds
        //cerr<<"(out of bounds)= i"<<nan("1")<<endl;
        return nan("1");
    } else if(upperIndex == m_profileDetails.startIndex) { // lower bound
        //cerr<<"(lower bound)= i"<<upperIndex<<endl;
        return upperIndex;
    } else if(lowerIndex == m_profileDetails.endIndex) { // upper bound
        //cerr<<"(upper bound)= i"<<lowerIndex<<endl;
        return lowerIndex;
    }

    ScalarType lowerValue = array[lowerIndex], upperValue = array[upperIndex];

    if(threshold < lowerValue && threshold < upperValue) { // the threshold is less than both values therefore take the smaller
        double ind = (lowerValue < upperValue) ? lowerIndex : upperIndex;
        //cerr<<"(thresh less than values)= i"<<ind<<endl;
        return (lowerValue < upperValue) ? lowerIndex : upperIndex;
    } else if(threshold > lowerValue && threshold > upperValue) { // the threshold is greater than both values therefore take the larger
        double ind = (lowerValue > upperValue) ? lowerIndex : upperIndex;
        //cerr<<"(thresh greater than values)= i"<<ind<<endl;
        return (lowerValue > upperValue) ? lowerIndex : upperIndex;
    } else { // in bounds and threshold between two values.

        ScalarType weight1 = 1 - fabs(threshold - array[lowerIndex])/fabs(array[lowerIndex] - array[upperIndex]);
        ScalarType weight2 = 1 - fabs(threshold - array[upperIndex])/fabs(array[lowerIndex] - array[upperIndex]);

        //cerr<<"(weighter average)= i"<<weight1 * lowerIndex + weight2 * upperIndex<<endl;
        return weight1 * lowerIndex + weight2 * upperIndex;
    }
}

int LinearRegressionCalculator::CalculateThresholdCrossingIndex(itkArrayType array, int startIndex, int endIndex, ScalarType threshold) {

    int index = endIndex;
    bool positiveCrossing = (array[startIndex] < threshold) ? true : false;

    int lowIndex = (endIndex > startIndex) ? startIndex : endIndex;
    int highIndex = (endIndex > startIndex) ? endIndex : startIndex;
    //cerr<<"threshold="<<threshold<<"; startIndex="<<startIndex<<"; endIndex="<<endIndex<<"; lowIndex="<<lowIndex<<"; highIndex="<<highIndex<<"; positiveCrossing="<<positiveCrossing<<endl;
    for (int i=lowIndex; i<highIndex; i++) {

        int correctedIndex = (endIndex > startIndex) ? i : highIndex - (i - lowIndex);
        bool crossing = (positiveCrossing) ? (array[correctedIndex] > threshold) : (array[correctedIndex] < threshold);
        //cerr<<"i="<<i<<"; correctedIndex="<<correctedIndex<<"; value="<<array->GetValue(correctedIndex)<<"; crossing="<<crossing<<endl;
        if(crossing) {
            index = (endIndex > startIndex) ? correctedIndex-1 : correctedIndex; // always return the index just prior to the profile meeting the threshold
            break;
        }
    }

    return index;
}

LinearRegressionCalculator::ScalarType LinearRegressionCalculator::CalculateLinearValue(RegressionLine line, ScalarType position) const {

    if(isnan(line.slope) || isnan(line.intercept)) {
        return line.meanY;
    } else {
        return line.slope * position + line.intercept;
    }
}

LinearRegressionCalculator::ScalarType LinearRegressionCalculator::CalculateWelchsTTest(RegressionLine line1, RegressionLine line2) {

    ScalarType var1 = line1.stdY * line1.stdY, var2 = line2.stdY * line2.stdY;
    int n1 = line1.n, n2 = line2.n;
    ScalarType t = (line1.meanY - line2.meanY) / sqrt(var1/n1 + var2/n2);
    int v = (var1/n1 + var2/n2) * (var1/n1 + var2/n2)
            / (var1*var1/(n1*n1*(n1-1)) + var2*var2/(n2*n2*(n2-1)));

    cout<<"Student T-test: t="<<t<<", v="<<v<<"given n1="<<n1<<" and n2="<<n2<<endl;

    return t;
}

// static processing methods
LinearRegressionCalculator::RegressionLine LinearRegressionCalculator::CalculateLinearRegression(itkArrayType array1, itkArrayType array2) {
    RegressionLine regressionLine = RegressionLine();
    int n = array1.Size();

    if(n==array2.Size()) { // only attempt to calculate for valid indices

        ScalarType sumX = 0, sumY = 0, sumXY = 0, sumXX = 0, sumYY = 0;

        for(int i=0; i<n; i++) {
            ScalarType x = array1[i];
            ScalarType y = array2[i];

            sumX += x; sumY += y; sumXY += x*y; sumXX += x*x; sumYY += y*y;

        }

        ScalarType meanX = sumX / n, meanY = sumY / n;
        ScalarType stdX = sqrt(sumXX / n - meanX * meanX);
        ScalarType stdY = sqrt(sumYY / n - meanY * meanY);

        regressionLine.n = n;
        regressionLine.meanX = meanX; regressionLine.stdX = stdX;
        regressionLine.meanY = meanY; regressionLine.stdY = stdY;

        if(n==1) { // if coincident start & stop indices - nan slope and intercept, known meanX and meany
            return regressionLine;
        }

        ScalarType covXY = sumXY / n - meanX * meanY, corrXY = covXY / (stdX * stdY);
        ScalarType slope = corrXY * stdY / stdX, intercept = meanY - slope * meanX;

        regressionLine.corrXY = corrXY; regressionLine.covXY = covXY;
        regressionLine.slope = slope; regressionLine.intercept = intercept;

        // calculate 'goodness of fit'
        ScalarType rSum = 0;//, rDenominator = 0;
        for(int i=0; i<n; i++) {
            ScalarType y = array2[i];
            ScalarType f = array1[i] * slope + intercept;

            rSum += (y - f) * (y - f);
        }
        regressionLine.rSquared = 1 - rSum / (stdY*stdY*n);
    }

    return regressionLine;
}