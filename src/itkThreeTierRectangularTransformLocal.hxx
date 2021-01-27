/* 
 * File:   itkThreeTierRectangularTransformLocal.hxx
 * Author: rap58
 *
 * Created on 10 November 2014, 17:31
 */

#ifndef ITKTHREETIERRECTANGULARTRANSFORMLOCAL_HXX
#define	ITKTHREETIERRECTANGULARTRANSFORMLOCAL_HXX

#include "itkThreeTierRectangularTransformLocal.h"
#include "vnl/algo/vnl_svd.h"
#include "itkEndostealRampTransformLocal.h"

namespace itk
{

//----------------- Model Methods ----------//

    // Constructor with default arguments
    template <typename TScalar>
    ThreeTierRectangularTransform<TScalar>
    ::ThreeTierRectangularTransform() :  Superclass(ParametersDimension){

        this->m_NumberOfDisplays = 13;

        m_ST = NumericTraits<ScalarType>::Zero;
        m_CB= NumericTraits<ScalarType>::Zero;
        m_TB = NumericTraits<ScalarType>::Zero;

        m_STdiff = NumericTraits<ScalarType>::Zero;
        m_TBdiff = NumericTraits<ScalarType>::Zero;

        m_sigma = NumericTraits<ScalarType>::Zero;

        m_P = NumericTraits<ScalarType>::Zero;
        m_E = NumericTraits<ScalarType>::Zero;

        m_thickness = NumericTraits<ScalarType>::Zero;

        this->m_valid = false;

        //this->m_Parameters->Fill(0);

    }

    // Destructor
    template <typename TScalar>
    ThreeTierRectangularTransform<TScalar>::
    ~ThreeTierRectangularTransform(){
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::TransformValue(TScalar position) const {
        // calculate model value at position
        double xPNormalised, xENormalised;
        xPNormalised = (position - m_P)/(m_sigma*std::sqrt(2));
        xENormalised = (position - m_E)/(m_sigma*std::sqrt(2));

        double value = m_ST+0.5*(m_STdiff)*(1+erf(xPNormalised))+0.5*(m_TBdiff)*(1+erf(xENormalised));

        return value;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::UntransformedValue(ScalarType position) const {

        if(position < m_P) {
            return m_ST;
        } else if(position < m_E) {
            return m_CB;
        } else {
            return m_TB;
        }
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::PeakTransformValue() const {
        ScalarType x_peak = m_sigma*m_sigma/(m_E-m_P) * log((m_CB-m_ST)/(m_CB-m_TB)) + (m_E+m_P)/2;
        return TransformValue(x_peak);
    }

    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::FixCorticalDensity(const ScalarType value) {
        this->FixParameter(value, kCorticalDensityParam);
    }

    template <typename TScalar>
    bool ThreeTierRectangularTransform<TScalar>
    ::isSigmaFixed() const {
        return this->m_combinedParameterMapping[kSigmaParam] == 0;
    }

    template <typename TScalar>
    bool ThreeTierRectangularTransform<TScalar>
    ::isCorticalDensityFixed() const {
        return this->m_combinedParameterMapping[kCorticalDensityParam] == 0;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::PeriostealError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if(m_P < this->m_start) { // positive thickness required
            error += std::exp(multiplicationFactor * (this->m_start-m_P)) - 1; // so exhibits 1st order (derivative) continuity. Increases continuously from 0
        }

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::EndostealError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        ScalarType maxE = std::min(std::min(m_P, this->m_maxProfileOffset) + this->m_length*(1.0-this->m_profileEdgeRatio), this->m_stop); // calculate given re-sampling that would occur if this were selected
        if (m_E > maxE) {
            error += std::exp(multiplicationFactor * (m_E - maxE)) - 1;
        }

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::ThicknessError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if (m_E < m_P) {
            error += std::exp(multiplicationFactor * (m_P - m_E)) - 1;
        }

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::CBMaxError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if (m_CB > this->m_MaxDensity) {
            error += std::exp(multiplicationFactor * (m_CB - this->m_MaxDensity)) - 1;
        }

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::CBTBDiffError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if (m_CB - m_TB < 0.0) { // cortical density > trabeculae
            error += std::exp(multiplicationFactor * (m_TB - m_CB)) - 1;
        }

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::CBSTDiffError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if (m_CB - m_ST < 0.0) { // cortical density > soft tissue
            error += std::exp(multiplicationFactor * (m_ST - m_CB)) - 1;
        }

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::SigmaError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if (m_sigma < 0.0) { // positive blur required
            error += std::exp(multiplicationFactor * -1.0 * m_sigma) - 1; // positive if too small
        } else if(m_sigma > this->m_length/4.0) { // sw compatability - max sigma < 1/4 profile length - negative if too big
            error -= std::exp(multiplicationFactor * (m_sigma - this->m_length/4.0)) - 1;
        } // NOTE - slope is covered by width and density diffs

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }


    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::UpdateInternalParameters() {
        // update internal variables

        // pre-compute values for transform calculation
        m_P = this->m_combinedParameters[0];
        m_E = this->m_combinedParameters[1];

        m_TB = this->m_combinedParameters[2];
        m_CB = this->m_combinedParameters[3];
        m_ST = this->m_combinedParameters[4];

        m_sigma = this->m_combinedParameters[5];

        m_STdiff = m_CB - m_ST;
        m_TBdiff = m_TB - m_CB;

        m_thickness = m_E-m_P;

        // work out if the values are valid
        ScalarType maxE = std::min(std::min(m_P, this->m_maxProfileOffset) + this->m_length*(1.0-this->m_profileEdgeRatio), this->m_stop); // calculate given re-sampling that would occur if this were selected
        if(isnan(m_P) || isnan(m_E) || isnan(m_TB) || isnan(m_CB) || isnan(m_ST) || isnan(m_sigma)) {
            this->m_valid = false;
        } else if ( m_P>=this->m_start && m_E > m_P && m_E <= maxE && m_CB > m_ST && m_CB > m_TB && m_sigma > 0 && m_CB <= this->m_MaxDensity && m_sigma < this->m_length / 4.0) {
            this->m_valid = true;
        } else {
            this->m_valid = false;
        }

    }

    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::SetMaxDensity(ScalarType maxDensity) {
        Superclass::SetMaxDensity(maxDensity);
        this->m_ParameterErrorScales[3] = this->m_MaxErrorMultiplier/this->m_MaxDensity;
        this->m_ParameterErrorScales[4] = this->m_MaxErrorMultiplier/this->m_MaxDensity;
        this->m_ParameterErrorScales[5] = this->m_MaxErrorMultiplier/this->m_MaxDensity;
    }

    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::SetProfileLength(ScalarType length) {
        Superclass::SetProfileLength(length);
        this->m_ParameterErrorScales[0] = this->m_MaxErrorMultiplier/this->m_length;
        this->m_ParameterErrorScales[1] = this->m_MaxErrorMultiplier/this->m_length;
        this->m_ParameterErrorScales[2] = this->m_MaxErrorMultiplier/(2.0*this->m_length);
        this->m_ParameterErrorScales[6] = this->m_MaxErrorMultiplier/(this->m_length/4.0);
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetTotalThickness() const {
        return m_thickness;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetCorticalThickness() const {
        return m_thickness;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetEndostealThickness() const {
        return 0;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetCorticalDensity() const {
        return m_CB;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetMeanCorticalDensity() const {
        return m_CB;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetTrabeculaeDensity() const {
        return m_TB; // trabecular density
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetSoftTissueDensity() const {
        return m_ST; // trabecular density
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetSurfaceThickness() const {
        return m_thickness;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetPeriostealEdgePosistion() const {
        return m_P;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetEndostealEdgePosistion() const {
        return m_E;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetSigma() const {
        return m_sigma;
    }

    template <typename TScalar>
    int ThreeTierRectangularTransform<TScalar>
    ::GetCorticalDensityFlexIndex() const {
        if(this->m_combinedParameterMapping[kCorticalDensityParam]==0) {
            return -1; // their is no valid cortical density index as it is not free to change
        }

        int flexIndex = 0;
        for(int i=0; i < kCorticalDensityParam; i++) {
            flexIndex+=this->m_combinedParameterMapping[i];
        }

        return flexIndex;
    }

    template <typename TScalar>
    int ThreeTierRectangularTransform<TScalar>
    ::GetPeriostealEdgeFlexIndex() const {
        if(this->m_combinedParameterMapping[kPeriostealEdgeParam]==0) {
            return -1; // their is no valid cortical density index as it is not free to change
        }

        int flexIndex = 0;
        for(int i=0; i < kPeriostealEdgeParam; i++) {
            flexIndex+=this->m_combinedParameterMapping[i];
        }

        return flexIndex;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ParametersType
    ThreeTierRectangularTransform<TScalar>
    ::GetDisplayValues() const {

        ParametersType values(this->m_NumberOfDisplays);
        values.Fill(nan("1"));

        if(this->m_valid || this->m_combinedParameterMapping[0]==0.0) { // if valid or fixed
            values[this->kPeriPositionDisp] = m_P;
        }

        if(this->m_valid || this->m_combinedParameterMapping[1]==0.0) { // if valid or fixed
            values[this->kEndoPositionDisp] = m_E; values[this->kEndoCBPositionDisp] = m_E; values[this->kEndoTBPositionDisp] = m_E;
        }

        if(!isnan(values[this->kPeriPositionDisp]) && !isnan(values[this->kEndoPositionDisp])) { // if valid or fixed
            values[this->kTotalThicknessDisp] = m_thickness;
            values[this->kCorticalThicknessDisp] = m_thickness;
            values[this->kEndoThicknessDisp] = 0.0;
        }

        if(this->m_valid || this->m_combinedParameterMapping[3]==0.0) { // if valid or fixed
            values[this->kDenseCorticalBoneDensityDisp] = m_CB;
        }

        if(this->m_valid || this->m_combinedParameterMapping[4]==0.0) { // if valid or fixed
            values[this->kSoftTissueDensityDisp] = m_ST;
        }

        if(this->m_valid || this->m_combinedParameterMapping[2]==0.0) { // if valid or fixed
            values[this->kTrabeculaeDensityDisp] = m_TB;
        }

        if(!isnan(values[this->kDenseCorticalBoneDensityDisp]) && !isnan(values[this->kTotalThicknessDisp])) { // if valid or fixed
            values[this->kMassSADensityDisp] = m_thickness * m_CB;
        }

        if(this->m_valid || this->m_combinedParameterMapping[5]==0.0) { // if valid or fixed
            values[this->kSigmaDisp] = m_sigma;
        }


        // leave error mean blank as filled by metrix's

        return values;
    }

    // Provides a weighting in an attempt to constrain the parameter space
    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::GetParameterErrorSum() const {

        ParametersType parameterErrors = GetParameterErrors();
        int n = parameterErrors.GetSize();

        ScalarType errorSum = 0;
        for(int i=0; i<n; i++) {
            errorSum += parameterErrors[i]/n;
        }

        if(std::isnan(errorSum) || std::isinf(errorSum)) {
            errorSum = std::numeric_limits<double>::max();
        }

        return errorSum;

    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::PointArrayType
    ThreeTierRectangularTransform<TScalar>
    ::GetModelPoints(ScalarType xStartPosition) const {
        PointArrayType modelPoints(6, 2);

        // Soft Tissue region
        modelPoints.put(0,0,  xStartPosition); modelPoints.put(0,1, m_ST);
        modelPoints.put(1,0,m_P); modelPoints.put(1,1, m_ST);
        // Cortical Bone region
        modelPoints.put(2,0,m_P); modelPoints.put(2,1, m_CB);
        modelPoints.put(3,0,m_E); modelPoints.put(3,1, m_CB);
        // Trabecular Bone region
        modelPoints.put(4,0,m_E); modelPoints.put(4,1, m_TB);
        modelPoints.put(5,0,xStartPosition+this->m_length); modelPoints.put(5,1, m_TB);

        return modelPoints;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::PointArrayType
    ThreeTierRectangularTransform<TScalar>
    ::GetModelPoints(ScalarType xStartPosition, ParametersType parameters) const {
        PointArrayType modelPoints(6, 2);

        ScalarType stDensity = parameters[kSoftTissueDensityParam], cbDensity = parameters[kCorticalDensityParam];
        ScalarType tbDensity = parameters[kTrabeculaeDensityParam];
        ScalarType periEdge = parameters[kPeriostealEdgeParam];
        ScalarType endoEdge = parameters[kEnodstealEdgeParam];

        // Soft Tissue region
        modelPoints.put(0,0,  xStartPosition); modelPoints.put(0,1, stDensity);
        modelPoints.put(1,0, periEdge); modelPoints.put(1,1, stDensity);
        // Cortical Bone region
        modelPoints.put(2,0, periEdge); modelPoints.put(2,1, cbDensity);
        modelPoints.put(3,0, endoEdge); modelPoints.put(3,1, cbDensity);
        // Trabecular Bone region
        modelPoints.put(4,0, endoEdge); modelPoints.put(4,1, tbDensity);
        modelPoints.put(5,0,xStartPosition+this->m_length); modelPoints.put(5,1, tbDensity);

        return modelPoints;
    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ParametersType
    ThreeTierRectangularTransform<TScalar>
    ::GetParameterErrors() const {

        ParametersType parameterErrors(ParametersDimension+1);
        parameterErrors[0] = PeriostealError(this->m_ParameterErrorScales[0]);
        parameterErrors[1] = EndostealError(this->m_ParameterErrorScales[1]);
        parameterErrors[2] = ThicknessError(this->m_ParameterErrorScales[2]);
        parameterErrors[3] = CBTBDiffError(this->m_ParameterErrorScales[3]);
        parameterErrors[4] = CBMaxError(this->m_ParameterErrorScales[4]);
        parameterErrors[5] = CBSTDiffError(this->m_ParameterErrorScales[5]);
        parameterErrors[6] = SigmaError(this->m_ParameterErrorScales[6]);

        return parameterErrors;

    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::NumberOfParametersType
    ThreeTierRectangularTransform<TScalar>
    ::GetNumberOfErrors(void) const {
        return ParametersDimension+1;
    }

    // Print self
    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::PrintSelf(std::ostream & os, Indent indent) const{
        Superclass::PrintSelf(os, indent);
    }

    // Reset the transform to an identity transform
    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::SetIdentity(void){
        this->Superclass::SetIdentity();
        m_ST = NumericTraits<TScalar>::Zero;
        m_STdiff = NumericTraits<TScalar>::Zero;
        m_TBdiff = NumericTraits<TScalar>::Zero;
        m_sigma = NumericTraits<TScalar>::Zero;
        m_P = NumericTraits<TScalar>::Zero;
        m_E = NumericTraits<TScalar>::Zero;

        m_thickness = NumericTraits<TScalar>::Zero;

        m_CB = NumericTraits<TScalar>::Zero;
        m_TB = NumericTraits<TScalar>::Zero;

    }

    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::CalculateWeightingFunction1(ArrayType &array, ScalarType startPosition, ScalarType scale) const { // start position = profile inner end

        int n = array.size(); ScalarType increment = this->m_length / n;

        if(m_P<this->m_start || m_P>this->m_stop) { // invalid m_P so set to constant value
            array.Fill(scale/n);
            return;
        }

        // calculate the overall weighting to ensure it remains a constant area

        ScalarType kbkgnd = 0.1, kline=0.15, kgauss = 1.0-kbkgnd-kline;

        // gaussian values
        ScalarType mean=m_P;
        ScalarType sigma1=1.7, sigma2=2.0;

        // line values
        ScalarType m1=kline/(this->m_profileEdgeRatio*this->m_length), b1;
        ScalarType m2=kline/((this->m_profileEdgeRatio-1.0)*this->m_length), b2;

        if((mean-startPosition)<0) {
            b2=-m1*(this->m_stop); b1=(b2 + m2 * mean) - m1 * mean;
        } else {
            b1=-m1*(this->m_start); b2 = (b1 + m1 * mean) - m2 * mean;
        }

        array.Fill(0.0); // initalise all outside values to empty
        ScalarType areaSum = 0;

        // indicies
        int xPIndex = std::max(0, (int)ceil((mean-startPosition) / increment));
        int startIndex = std::max(0, (int)ceil((this->m_start-startPosition) / increment)); // todo - work out why floor / ceil aren't working
        int endIndex = std::min(n, (int)floor((this->m_stop-startPosition) / increment));

        // calc values
        for(int i=startIndex; i<xPIndex; i++) {
            ScalarType x = startPosition+i*increment;
            ScalarType value = ( kbkgnd + m1*x+b1 + kgauss*exp(-(x-mean)*(x-mean)/(2*sigma1*sigma1)) );
            areaSum += value; array[i] = value;
        }

        for(int i=xPIndex; i<endIndex; i++) {
            ScalarType x = startPosition+i*increment;
            ScalarType value = ( kbkgnd + m2*x+b2 + kgauss*exp(-(x-mean)*(x-mean)/(2*sigma2*sigma2)) );
            areaSum += value; array[i] = value;
        }

        for(int i=startIndex; i<endIndex; i++) {
            array[i] = scale * array[i]/areaSum;
        }

        return;
    }

    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::CalculateWeightingFunction2(ArrayType & array, ScalarType startPosition, ScalarType scale) const  {

        int n = array.size(); ScalarType increment = this->m_length / n;

        if(m_P<this->m_start || m_P>this->m_stop) { // invalid m_P so set to constant value
            array.Fill(scale/n);
            return;
        }

        // calculate the overall weighting to ensure it remains a constant area

        ScalarType kbkgnd = 0.25, kline=0.0, kgauss = 1.0-kbkgnd-kline;

        // gaussian values
        ScalarType mean=m_P;
        ScalarType sigma1=2.0, sigma2=2.5;

        // line values
        ScalarType m1=kline/(this->m_profileEdgeRatio*this->m_length), b1;
        ScalarType m2=kline/((this->m_profileEdgeRatio-1.0)*this->m_length), b2;

        if((mean-startPosition)<0) {
            b2=-m1*(this->m_stop); b1=(b2 + m2 * mean) - m1 * mean;
        } else {
            b1=-m1*(this->m_start); b2 = (b1 + m1 * mean) - m2 * mean;
        }

        array.Fill(0.0); // initalise all outside values to empty
        ScalarType areaSum = 0;

        // indicies
        int xPIndex = std::max(0, (int)ceil((mean-startPosition) / increment));
        int startIndex = std::max(0, (int)ceil((this->m_start-startPosition) / increment)); // todo - work out why floor / ceil aren't working
        int endIndex = std::min(n, (int)floor((this->m_stop-startPosition) / increment));

        // calc values
        for(int i=startIndex; i<xPIndex; i++) {
            ScalarType x = startPosition+i*increment;
            ScalarType value = ( kbkgnd + m1*x+b1 + kgauss*exp(-(x-mean)*(x-mean)/(2*sigma1*sigma1)) );
            areaSum += value; array[i] = value;
        }

        for(int i=xPIndex; i<endIndex; i++) {
            ScalarType x = startPosition+i*increment;
            ScalarType value = ( kbkgnd + m2*x+b2 + kgauss*exp(-(x-mean)*(x-mean)/(2*sigma2*sigma2)) );
            areaSum += value; array[i] = value;
        }

        for(int i=startIndex; i<endIndex; i++) {
            array[i] = scale * array[i]/areaSum;
        }

        return;
    }

    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::CalculateWeightingFunction3(ArrayType &array, ScalarType startPosition, ScalarType scale) const { // start position = profile inner end

        int n = array.size(); ScalarType increment = this->m_length / n;

        if(m_P<this->m_start || m_P>this->m_stop) { // invalid m_P so set to constant value
            array.Fill(scale/n);
            return;
        }

        // calculate the overall weighting to ensure it remains a constant area

        ScalarType kbkgnd = 0.05, kline=0.15, kgauss = 1.0-kbkgnd-kline;

        // gaussian values
        ScalarType mean=m_P;
        ScalarType sigma1=2.0, sigma2=2.35; // sigma1=1.4, sigma2=1.5;

        // line values
        ScalarType m1=kline/(this->m_profileEdgeRatio*this->m_length), b1;
        ScalarType m2=kline/((this->m_profileEdgeRatio-1.0)*this->m_length), b2;

        if((mean-startPosition)<0) {
            b2=-m1*(this->m_stop); b1=(b2 + m2 * mean) - m1 * mean;
        } else {
            b1=-m1*(this->m_start); b2 = (b1 + m1 * mean) - m2 * mean;
        }

        array.Fill(0.0); // initalise all outside values to empty
        ScalarType areaSum = 0;

        // indicies
        int xPIndex = std::max(0, (int)ceil((mean-startPosition) / increment));
        int startIndex = std::max(0, (int)ceil((this->m_start-startPosition) / increment)); // todo - work out why floor / ceil aren't working
        int endIndex = std::min(n, (int)floor((this->m_stop-startPosition) / increment));

        // calc values
        for(int i=startIndex; i<xPIndex; i++) {
            ScalarType x = startPosition+i*increment;
            ScalarType value = ( kbkgnd + m1*x+b1 + kgauss*exp(-(x-mean)*(x-mean)/(2*sigma1*sigma1)) );
            areaSum += value; array[i] = value;
        }

        for(int i=xPIndex; i<endIndex; i++) {
            ScalarType x = startPosition+i*increment;
            ScalarType value = ( kbkgnd + m2*x+b2 + kgauss*exp(-(x-mean)*(x-mean)/(2*sigma2*sigma2)) );
            areaSum += value; array[i] = value;
        }

        for(int i=startIndex; i<endIndex; i++) {
            array[i] = scale * array[i]/areaSum;
        }

        return;
    }

    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::CalculateWeightingFunction4(ArrayType & array, ScalarType startPosition, ScalarType scale) const  {

        int n = array.size(); ScalarType increment = this->m_length / n;

        if(m_P<this->m_start || m_P>this->m_stop) { // invalid m_P so set to constant value
            array.Fill(scale/n);
            return;
        }

        // calculate the overall weighting to ensure it remains a constant area

        ScalarType kbkgnd = 0.1, kline=0.25, kgauss = 1.0-kbkgnd-kline;

        // gaussian values
        ScalarType mean=m_P;
        ScalarType sigma1=2.0, sigma2=2.5;

        // line values
        ScalarType m1=kline/(this->m_profileEdgeRatio*this->m_length), b1;
        ScalarType m2=kline/((this->m_profileEdgeRatio-1.0)*this->m_length), b2;

        if((mean-startPosition)<0) {
            b2=-m1*(this->m_stop); b1=(b2 + m2 * mean) - m1 * mean;
        } else {
            b1=-m1*(this->m_start); b2 = (b1 + m1 * mean) - m2 * mean;
        }

        array.Fill(0.0); // initalise all outside values to empty
        ScalarType areaSum = 0;

        // indicies
        int xPIndex = std::max(0, (int)ceil((mean-startPosition) / increment));
        int startIndex = std::max(0, (int)ceil((this->m_start-startPosition) / increment)); // todo - work out why floor / ceil aren't working
        int endIndex = std::min(n, (int)floor((this->m_stop-startPosition) / increment));

        // calc values
        for(int i=startIndex; i<xPIndex; i++) {
            ScalarType x = startPosition+i*increment;
            ScalarType value = ( kbkgnd + m1*x+b1 + kgauss*exp(-(x-mean)*(x-mean)/(2*sigma1*sigma1)) );
            areaSum += value; array[i] = value;
        }

        for(int i=xPIndex; i<endIndex; i++) {
            ScalarType x = startPosition+i*increment;
            ScalarType value = ( kbkgnd + m2*x+b2 + kgauss*exp(-(x-mean)*(x-mean)/(2*sigma2*sigma2)) );
            areaSum += value; array[i] = value;
        }

        for(int i=startIndex; i<endIndex; i++) {
            array[i] = scale * array[i]/areaSum;
        }

        return;
    }

    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::CalculateWeightingDensity(ArrayType & array, ScalarType startPosition, ScalarType scale) const  {

        int n = array.size(); ScalarType increment = this->m_length / n;

        if(m_P<this->m_start || m_P>this->m_stop || m_E>this->m_stop) { // invalid m_P so set to constant value
            array.Fill(scale/n);
            return;
        }

        // calculate the overall weighting to ensure it remains a constant area

        ScalarType kbkgnd = 0.2, kline=0.05, kgauss = 1.0-kbkgnd-kline;

        // gaussian values
        ScalarType mean1=m_P, mean2=(m_E>m_P) ? m_E : m_P;
        ScalarType sigma1=2.0, sigma2=2.5;

        // line values
        ScalarType m1=kline/(this->m_profileEdgeRatio*this->m_length), m2 = -m1, b1, b2;
        ScalarType  kLineTotal;

        if((mean1-this->m_start)>(this->m_stop-mean2)) {
            b1=-m1 *(this->m_start); b2 = b1 + m1*mean1 - m2*mean2;
            kLineTotal=(mean1-this->m_start)* m1 + kbkgnd + kgauss;
        } else {
            b2=-m2 *(this->m_stop); b1= b2 + m2*mean2 - m1*mean1;
            kLineTotal=(mean2-this->m_stop)* m2 + kbkgnd + kgauss;
        }

        array.Fill(0.0); // initalise all outside values to empty
        ScalarType areaSum = 0;

        // indicies
        int xPIndex = std::max(0, (int)ceil((mean1 -startPosition) / increment));
        int xEIndex = std::min(n, (int)floor((mean2-startPosition) / increment));
        int startIndex = std::max(0, (int)ceil((this->m_start-startPosition) / increment)); // todo - work out why floor / ceil aren't working
        int endIndex = std::min(n, (int)floor((this->m_stop-startPosition) / increment));

        // calc values
        for(int i=startIndex; i<xPIndex; i++) {
            ScalarType x = startPosition+i*increment;
            ScalarType value = ( kbkgnd + m1 *x+b1 + kgauss*exp(-(x- mean1)*(x- mean1)/(2*sigma1*sigma1)) );
            areaSum += value; array[i] = value;
        }


        for(int i=xPIndex; i<xEIndex; i++) {
            areaSum += kLineTotal; array[i] = kLineTotal;
        }

        for(int i=xEIndex; i<endIndex; i++) {
            ScalarType x = startPosition+i*increment;
            ScalarType value = ( kbkgnd + m2 *x+b2 + kgauss*exp(-(x-mean2)*(x-mean2)/(2*sigma2*sigma2)) );
            areaSum += value; array[i] = value;
        }

        for(int i=startIndex; i<endIndex; i++) {
            array[i] = scale * array[i]/areaSum;
        }

        return;
    }

    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::CalculateWeightingHeuristic(ArrayType & array, ScalarType startPosition, ScalarType scale) const  {

        // This weighting function has a peak at static_weighting_location and is designed such that
        // if the peak is at 1/3, then the weight at that point will be 1 and the weight at the two
        // ends will be 0.05 (this can be altered in expression below)
        // It deliberately emphasises the outer cortical edge, since this is the best defined transition
        // which by default is position at 1/3 of the length of the line [from SW]

        int n = array.size(); ScalarType increment = this->m_length / n;

        int startIndex = std::max(0, (int)ceil((startPosition-this->m_start) / increment)); // todo - work out why floor / ceil aren't working
        int endIndex = std::min(n, (int)floor((this->m_stop-startPosition) / increment));


        ScalarType p = m_P+this->m_profileEdgeRatio*this->m_length; // should generally be about 1/3
        float f = (this->m_length/p)*(1.0/0.05 - 1.0);

        ScalarType areaSum = 0;
        for(int i=startIndex; i<endIndex; i++) {

            ScalarType x =  startPosition + i*increment-this->m_start;
            float l1 = (x + p) / this->m_length;
            float l2 = (x - p) / this->m_length;

            array[i] = l1 / (l2 * l2 * f + l1); // area of approximately 6.6
            areaSum+=array[i];
        }

        // could add a normalisation step

        for(int i=startIndex; i<endIndex; i++) {
            array[i] = scale * array[i]/areaSum;
        }


    }

    template <typename TScalar>
    typename ThreeTierRectangularTransform<TScalar>::ScalarType
    ThreeTierRectangularTransform<TScalar>
    ::CalculateSigmaAdjustedCorticalDensity(ScalarType peakSampledValue, ScalarType globalSigma) {
        ScalarType backgroundDensity = std::max(m_ST,m_TB);

        // calculate model max for global sigma - temporarily change the sigma parameter
        ScalarType sigma = m_sigma;
        m_sigma = globalSigma;
        ScalarType peakTransformValue = PeakTransformValue();
        m_sigma = sigma;

        // calculate the sigma adjusted density estimate
        ScalarType density = (m_CB-backgroundDensity) * (peakSampledValue - backgroundDensity) / (peakTransformValue - backgroundDensity) + backgroundDensity;
        density=std::max(std::min(density, this->m_MaxDensity), backgroundDensity); // limit to within reasonable bounds
        return density;
    }

    // Compute transformation Jacobian
    template <typename TScalar>
    void ThreeTierRectangularTransform<TScalar>
    ::ComputeJacobianWithRespectToParameters(const InputPointType & p, JacobianType & j ) const {
        std::cout<<"TODO - Implement 'ThreeTierRectangularTransform::ComputeJacobianWithRespectToParameters'"<<std::endl;
    }

} // namespace

#endif	/* ITKTHREETIERRECTANGULARTRANSFORMLOCAL_HXX */

