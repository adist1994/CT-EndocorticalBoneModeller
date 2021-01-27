/* 
 * File:   itkEndostealRampTransformLocal.hxx
 * Author: rap58
 *
 * Created on 10 November 2014, 17:31
 */

#ifndef ITKENDOSTEALRAMPTRANSFORMLOCAL_HXX
#define	ITKENDOSTEALRAMPTRANSFORMLOCAL_HXX

#include "itkEndostealRampTransformLocal.h"
#include "vnl/algo/vnl_svd.h"
#include "erfIntegral.h"

namespace itk
{

//----------------- Model Methods ----------//

// Constructor with default arguments
    template <typename TScalar>
    EndostealRampTransform<TScalar>
    ::EndostealRampTransform() :  Superclass(ParametersDimension){

        this->m_NumberOfDisplays = 13;

        erfIntLookup = new ERFIntLookup();
        erfIntLookup->setupTable(ERFIntLookup::SimpsonsRule);

        m_ST = NumericTraits<ScalarType>::Zero;
        m_DB = NumericTraits<ScalarType>::Zero;
        m_TB = NumericTraits<ScalarType>::Zero;

        m_sigma = NumericTraits<ScalarType>::Zero;

        m_Ecb = NumericTraits<ScalarType>::Zero;
        m_Etb = NumericTraits<ScalarType>::Zero;
        m_P = NumericTraits<ScalarType>::Zero;
        m_E = NumericTraits<ScalarType>::Zero;

        m_DBThickness = NumericTraits<ScalarType>::Zero;
        m_thickness = NumericTraits<ScalarType>::Zero;
        m_rampWidth = NumericTraits<ScalarType>::Zero;

        this->m_valid = false;

        //this->m_Parameters->Fill(0);

    }

// Destructor
    template <typename TScalar>
    EndostealRampTransform<TScalar>::
    ~EndostealRampTransform(){
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::TransformValue(TScalar position) const {
        // calculate model value at position
        double xEcbNormalised, xEtbNormalised, xPNormalised, constant, scale;
        constant = m_ST + (m_TB - m_DB) / 2;
        scale = m_slope * std::sqrt(2) * m_sigma * 0.5;

        xEcbNormalised = (position - m_Ecb)/(m_sigma*std::sqrt(2));
        xEtbNormalised = (position - m_Etb)/(m_sigma*std::sqrt(2));
        xPNormalised =  (position - m_P) /(m_sigma*std::sqrt(2));

        double value;
        if(m_rampWidth>MINRAMPWIDTH) {
            value = constant
                    + 0.5*(m_DB - m_ST)*(1+erf(xPNormalised))
                    + scale * erfIntLookup->erfIntegralMethod(xEcbNormalised)
                    - scale * erfIntLookup->erfIntegralMethod(xEtbNormalised);
        } else {
            // approximate with rect function in the limit of ramp with approaches zero
            double xENormalised =  (position - m_E) /(m_sigma*std::sqrt(2));
            value = m_ST
                    + 0.5*(m_DB - m_ST)*(1+erf(xPNormalised))
                    + 0.5*(m_TB - m_DB)*(1+erf(xENormalised));
        }

        return value;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::UntransformedValue(ScalarType position) const {

        if(position < m_P) {
            return m_ST;
        } else if(position < m_Ecb) {
            return m_DB;
        } else if(position < m_Etb) {
            double m = (m_TB - m_DB) / (m_Etb - m_Ecb);
            return m_DB + m * (position - m_Ecb);
        } else {
            return m_TB;
        }
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::PeakTransformValue() const {

        // calculate peak value - manually
        ScalarType peakValue = TransformValue(0);
        int numberSamples = 100;
        double increment = this->m_length / numberSamples;
        for(int i=1; i<numberSamples; i++) {
            double x = increment * i;
            double value = TransformValue(x);

            peakValue = (value > peakValue) ? value : peakValue;
        }

        return peakValue;

//  ScalarType x_peak = m_sigma*m_sigma/(m_E-m_P) * log((m_DB-m_ST)/(m_DB-m_TB)) + (m_E+m_P)/2;
//  return TransformValue(x_peak);
    }

    template <typename TScalar>
    void EndostealRampTransform<TScalar>
    ::FixCorticalDensity(const ScalarType value) {
        this->FixParameter(value, kCorticalDensityParam);
    }

    template <typename TScalar>
    bool EndostealRampTransform<TScalar>
    ::isSigmaFixed() const {
        return this->m_combinedParameterMapping[kSigmaParam] == 0;
    }

    template <typename TScalar>
    bool EndostealRampTransform<TScalar>
    ::isCorticalDensityFixed() const {
        return this->m_combinedParameterMapping[kCorticalDensityParam] == 0;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
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
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::EndostealError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        ScalarType maxEtb = std::min(std::min(m_P, this->m_maxProfileOffset) + this->m_length*(1.0-this->m_profileEdgeRatio), this->m_stop);
        if (m_Etb > maxEtb) { // position bounded by profile length required
            error += std::exp(multiplicationFactor * (m_Etb - maxEtb)) - 1;
        }

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::ThicknessError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if (m_Ecb < m_P) { // positive thickness required
            error += std::exp(multiplicationFactor * (m_P-m_Ecb)) - 1;
        }

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::CBMaxError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if (m_DB > this->m_MaxDensity) {
            error += std::exp(multiplicationFactor * (m_DB - this->m_MaxDensity)) - 1;
        }

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::CBTBDiffError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if (m_DB - m_TB < 0.0) { // cortical density > trabeculae
            error += std::exp(multiplicationFactor * (m_TB - m_DB)) - 1;
        }

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::CBSTDiffError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if (m_DB - m_ST < 0.0) { // cortical density > soft tissue
            error += std::exp(multiplicationFactor * (m_ST - m_DB)) - 1;
        }

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::SigmaError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if (m_sigma < 0.0) { // positive blur required
            error += std::exp(multiplicationFactor * -1.0 * m_sigma) - 1; // TODO - make neg or pos depending on if to big or too small
        } else if(m_sigma > this->m_length/4.0) { // sw compatability - max sigma < 1/4 profile length
            error -= std::exp(multiplicationFactor * (m_sigma - this->m_length/4.0)) - 1;
        } // NOTE - slope is covered by width and density diffs

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::RampWidthError(ScalarType multiplicationFactor) const {
        ScalarType error = 0;
        if(m_rampWidth < 0.0) { // ramp width
            error += std::exp(multiplicationFactor * -1.0 * m_rampWidth) - 1;
        } // NOTE - slope is covered by width and density diffs

        /*// todo decide weather or not to keep
        if(fabs(m_DBThickness)/fabs(m_rampWidth) < 0.1) { // ramp width - require to be no more than 10x the cb thickness
            error += std::exp(multiplicationFactor * (fabs(m_rampWidth)/fabs(m_DBThickness))-10) - 1;
        }*/ // NOTE - slope is covered by width and density diffs

        if(std::isnan(error) || std::isinf(error)) { // check for valid error
            error = std::numeric_limits<double>::max();
        }

        return error;
    }

    template <typename TScalar>
    void EndostealRampTransform<TScalar>
    ::UpdateInternalParameters() {
        // update the internal variables
        // pre-compute values for transform calculation
        m_rampWidth = this->m_combinedParameters[6];

        m_P = this->m_combinedParameters[0];
        m_E = this->m_combinedParameters[1];

        m_TB = this->m_combinedParameters[2];
        m_DB = this->m_combinedParameters[3];
        m_ST = this->m_combinedParameters[4];

        m_thickness = m_E-m_P;

        m_Ecb = m_E - 1.0/2.0 * m_rampWidth;
        m_Etb = m_E + 1.0/2.0 * m_rampWidth;

        m_DBThickness = m_Ecb-m_P;

        m_sigma = this->m_combinedParameters[5];

        m_slope = (m_TB - m_DB) / (m_Etb - m_Ecb);

        // view variables
        //cout<<"ERT: xP="<<m_P<<",xE="<<m_E<<",RW="<<m_rampWidth<<",yTB="<<m_TB<<",yCB="<<m_CB<<",yST="<<m_ST<<",sig="<<m_sigma<<endl;


        // work out if the values are valid
        ScalarType maxEtb = std::min(std::min(m_P, this->m_maxProfileOffset) + this->m_length*(1.0-this->m_profileEdgeRatio), this->m_stop); // calculate give re-sampling that would occur if this were selected
        //ScalarType minP = std::max(this->m_start, -this->m_maxProfileOffset); // todo - decide weather of not to limit the offset but m_MaxProfileOffset
        if(isnan(m_P) || isnan(m_E) || isnan(m_TB) || isnan(m_DB) || isnan(m_ST) || isnan(m_sigma) || isnan(m_rampWidth) ) {
            this->m_valid = false;
        } else if ( m_P >= this->m_start && m_Ecb > m_P && m_Etb <= maxEtb && m_DB > m_ST && m_DB > m_TB && m_sigma > 0 && m_DB <= this->m_MaxDensity && m_slope < 0 && m_sigma < this->m_length/4.0) { // sw compatibility
            this->m_valid = true;
        } else {
            this->m_valid = false;
        }

    }

    template <typename TScalar>
    void EndostealRampTransform<TScalar>
    ::SetMaxDensity(ScalarType maxDensity) {
        Superclass::SetMaxDensity(maxDensity);
        this->m_ParameterErrorScales[3] = this->m_MaxErrorMultiplier/this->m_MaxDensity;
        this->m_ParameterErrorScales[4] = this->m_MaxErrorMultiplier/this->m_MaxDensity;
        this->m_ParameterErrorScales[5] = this->m_MaxErrorMultiplier/this->m_MaxDensity;
        //cout<<"Error Scales: ES[3]="<<this->m_ParameterErrorScales[3]<<"ES[4]="<<this->m_ParameterErrorScales[3]<<"ES[5]="<<this->m_ParameterErrorScales[3]<<endl;
    }

    template <typename TScalar>
    void EndostealRampTransform<TScalar>
    ::SetProfileLength(ScalarType length) {
        Superclass::SetProfileLength(length);
        this->m_ParameterErrorScales[0] = this->m_MaxErrorMultiplier/this->m_length;
        this->m_ParameterErrorScales[1] = this->m_MaxErrorMultiplier/this->m_length;
        this->m_ParameterErrorScales[2] = this->m_MaxErrorMultiplier/(2.0*this->m_length);
        this->m_ParameterErrorScales[6] = this->m_MaxErrorMultiplier/(this->m_length/4.0);
        this->m_ParameterErrorScales[7] = this->m_MaxErrorMultiplier/this->m_length;
        //cout<<"Error Scales: ES[0]="<<this->m_ParameterErrorScales[0]<<"ES[1]="<<this->m_ParameterErrorScales[1]<<"ES[2]="<<this->m_ParameterErrorScales[2]<<"ES[6]="<<this->m_ParameterErrorScales[6]<<"ES[7]="<<this->m_ParameterErrorScales[7]<<"exp(max)"<<std::exp(this->m_MaxErrorMultiplier)<<"exp(710)"<<std::exp(710)<<endl;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::GetTotalThickness() const { // todo replace total with cortical
        return m_DBThickness + 0.05 * m_rampWidth;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::GetCorticalThickness() const { // todo replace cortical to dense thickness
        return m_DBThickness;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::GetEndostealThickness() const {
        return m_rampWidth;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::GetCorticalDensity() const {
        return m_DB;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::GetMeanCorticalDensity() const {
        //return m_DB * (m_DBThickness)/(m_DBThickness+m_rampWidth) + 0.5 * (m_DB+m_TB) * (m_rampWidth)/(m_DBThickness+m_rampWidth);

        return m_DB*m_DBThickness/(m_DBThickness+0.5*m_rampWidth) + (0.75*m_DB+0.25*m_TB) * (0.5*m_rampWidth)/(m_DBThickness+0.5*m_rampWidth);

    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::GetTrabeculaeDensity() const {
        return m_TB;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::GetSoftTissueDensity() const {
        return m_ST;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::GetSurfaceThickness() const { // todo remove
        return m_thickness;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::GetPeriostealEdgePosistion() const {
        return m_P;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::GetEndostealEdgePosistion() const {
        return m_E;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::GetSigma() const {
        return m_sigma;
    }

    template <typename TScalar>
    int EndostealRampTransform<TScalar>
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
    int EndostealRampTransform<TScalar>
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
    typename EndostealRampTransform<TScalar>::ParametersType
    EndostealRampTransform<TScalar>
    ::GetDisplayValues() const {

        ParametersType values(this->m_NumberOfDisplays);
        values.Fill(nan("1"));

        if(this->m_valid || this->m_combinedParameterMapping[0]==0.0) { // if valid or fixed
            values[this->kPeriPositionDisp] = m_P;
        }

        if(this->m_valid || this->m_combinedParameterMapping[1]==0.0) { // if valid or fixed
            values[this->kEndoPositionDisp] = m_E; values[this->kEndoCBPositionDisp] = m_Ecb; values[this->kEndoTBPositionDisp] = m_Etb;
        }

        if(!isnan(values[this->kEndoPositionDisp]) && !isnan(values[this->kEndoPositionDisp]) && (this->m_valid || this->m_combinedParameterMapping[6]==0.0)) { // if valid or fixed
            values[this->kTotalThicknessDisp] = m_DBThickness + 0.5 * m_rampWidth;
            values[this->kCorticalThicknessDisp] = m_DBThickness;
            values[this->kEndoThicknessDisp] = m_rampWidth;
        }

        if(this->m_valid || this->m_combinedParameterMapping[3]==0.0) { // if valid or fixed
            values[this->kDenseCorticalBoneDensityDisp] = m_DB;
        }

        if(this->m_valid || this->m_combinedParameterMapping[4]==0.0) { // if valid or fixed
            values[this->kSoftTissueDensityDisp] = m_ST;
        }

        if(this->m_valid || this->m_combinedParameterMapping[2]==0.0) { // if valid or fixed
            values[this->kTrabeculaeDensityDisp] = m_TB;
        }

        if(!isnan(values[this->kDenseCorticalBoneDensityDisp]) && !isnan(values[this->kTotalThicknessDisp])) { // if valid or fixed
            values[this->kMassSADensityDisp] = values[this->kTotalThicknessDisp] * m_DB;
        }

        if(this->m_valid || this->m_combinedParameterMapping[5]==0.0) { // if valid or fixed
            values[this->kSigmaDisp] = m_sigma;
        }


        // leave error mean blank as filled by metrix's

        return values;
    }

// Provides a weighting in an attempt to constrain the parameter space
    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
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
    typename EndostealRampTransform<TScalar>::PointArrayType
    EndostealRampTransform<TScalar>
    ::GetModelPoints(ScalarType xStartPosition) const {
        PointArrayType modelPoints(6, 2);

        // Soft Tissue region
        modelPoints.put(0,0,  xStartPosition); modelPoints.put(0,1, m_ST);
        modelPoints.put(1,0,m_P); modelPoints.put(1,1, m_ST);
        // Cortical Bone region
        modelPoints.put(2,0,m_P); modelPoints.put(2,1, m_DB);
        modelPoints.put(3,0,m_Ecb); modelPoints.put(3,1, m_DB);
        // Trabecular Bone region
        modelPoints.put(4,0,m_Etb); modelPoints.put(4,1, m_TB);
        modelPoints.put(5,0,xStartPosition+this->m_length); modelPoints.put(5,1, m_TB);

        return modelPoints;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::PointArrayType
    EndostealRampTransform<TScalar>
    ::GetModelPoints(ScalarType xStartPosition, ParametersType parameters) const {
        PointArrayType modelPoints(6, 2);

        ScalarType stDensity = parameters[kSoftTissueDensityParam], cbDensity = parameters[kCorticalDensityParam];
        ScalarType tbDensity = parameters[kTrabeculaeDensityParam];
        ScalarType periEdge = parameters[kPeriostealEdgeParam];
        ScalarType endoCBEdge = parameters[EnodstealEdgeParam]-0.5*parameters[kEndostealWidthParam];
        ScalarType endoTBEdge = parameters[EnodstealEdgeParam]+0.5*parameters[kEndostealWidthParam];

        // Soft Tissue region
        modelPoints.put(0,0,  xStartPosition); modelPoints.put(0,1, stDensity);
        modelPoints.put(1,0,periEdge); modelPoints.put(1,1, stDensity);
        // Cortical Bone region
        modelPoints.put(2,0,periEdge); modelPoints.put(2,1, cbDensity);
        modelPoints.put(3,0,endoCBEdge); modelPoints.put(3,1, cbDensity);
        // Trabecular Bone region
        modelPoints.put(4,0,endoTBEdge); modelPoints.put(4,1, tbDensity);
        modelPoints.put(5,0,xStartPosition+this->m_length); modelPoints.put(5,1, tbDensity);

        return modelPoints;
    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::ParametersType
    EndostealRampTransform<TScalar>
    ::GetParameterErrors() const {

        ParametersType parameterErrors(ParametersDimension+1);
        parameterErrors[0] = PeriostealError(this->m_ParameterErrorScales[0]);
        parameterErrors[1] = EndostealError(this->m_ParameterErrorScales[1]);
        parameterErrors[2] = ThicknessError(this->m_ParameterErrorScales[2]);
        parameterErrors[3] = CBTBDiffError(this->m_ParameterErrorScales[3]);
        parameterErrors[4] = CBMaxError(this->m_ParameterErrorScales[4]);
        parameterErrors[5] = CBSTDiffError(this->m_ParameterErrorScales[5]);
        parameterErrors[6] = SigmaError(this->m_ParameterErrorScales[6]);
        parameterErrors[7] = RampWidthError(this->m_ParameterErrorScales[7]);

        //cerr<<"ERRxP="<<parameterErrors[0]<<",xE="<<parameterErrors[1]<<",Th="<<parameterErrors[2]<<",CB="<<parameterErrors[3]<<",TB="<<parameterErrors[4]<<",ST="<<parameterErrors[5]<<",sig="<<parameterErrors[5]<<",RW"<<parameterErrors[6];//<<endl;
        return parameterErrors;

    }

    template <typename TScalar>
    typename EndostealRampTransform<TScalar>::NumberOfParametersType
    EndostealRampTransform<TScalar>
    ::GetNumberOfErrors(void) const {
        return ParametersDimension+1;
    }

// Print self
    template <typename TScalar>
    void EndostealRampTransform<TScalar>
    ::PrintSelf(std::ostream & os, Indent indent) const{
        Superclass::PrintSelf(os, indent);
    }

// Reset the transform to an identity transform
    template <typename TScalar>
    void EndostealRampTransform<TScalar>
    ::SetIdentity(void){
        this->Superclass::SetIdentity();
        m_ST = NumericTraits<ScalarType>::Zero;
        m_DB = NumericTraits<ScalarType>::Zero;
        m_TB = NumericTraits<ScalarType>::Zero;

        m_sigma = NumericTraits<ScalarType>::Zero;

        m_Ecb = NumericTraits<ScalarType>::Zero;
        m_Etb = NumericTraits<ScalarType>::Zero;
        m_P = NumericTraits<ScalarType>::Zero;
        m_E = NumericTraits<ScalarType>::Zero;

        m_DBThickness = NumericTraits<ScalarType>::Zero;
        m_thickness = NumericTraits<ScalarType>::Zero;
        m_rampWidth = NumericTraits<ScalarType>::Zero;


    }

    template <typename TScalar>
    void EndostealRampTransform<TScalar>
    ::CalculateWeightingFunction1(ArrayType &array, ScalarType startPosition, ScalarType scale) const { // start position = profile inner end

        int n = array.size(); ScalarType increment = this->m_length / n;

        if(m_P<this->m_start || m_P>this->m_stop) { // invalid m_P so set to constant value
            array.Fill(scale/n);
            return;
        }

        // calculate the overall weighting to ensure it remains a constant area

        ScalarType kbkgnd = 0.1, kline=0.15, kgauss = 1.0-kbkgnd-kline; // kbkgnd = 0.05, kline=0.15

        // gaussian values
        ScalarType mean=m_P;
        ScalarType sigma1=1.7, sigma2=2.0; //sigma1=2.0, sigma2=2.35

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
    void EndostealRampTransform<TScalar>
    ::CalculateWeightingFunction2(ArrayType & array, ScalarType startPosition, ScalarType scale) const  {

        int n = array.size(); ScalarType increment = this->m_length / n;

        if(m_P<this->m_start || m_P>this->m_stop) { // invalid m_P so set to constant value
            array.Fill(scale/n);
            return;
        }

        // calculate the overall weighting to ensure it remains a constant area

        ScalarType kbkgnd = 0.25, kline=0.0, kgauss = 1.0-kbkgnd-kline; // kbkgnd = 0.10, kline=0.25

        // gaussian values
        ScalarType mean=m_P;
        ScalarType sigma1=2.0, sigma2=2.5; // sigma1=2.0, sigma2=2.5

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
    void EndostealRampTransform<TScalar>
    ::CalculateWeightingFunction3(ArrayType &array, ScalarType startPosition, ScalarType scale) const { // start position = profile inner end

        int n = array.size(); ScalarType increment = this->m_length / n;

        if(m_P<this->m_start || m_P>this->m_stop) { // invalid m_P so set to constant value
            array.Fill(scale/n);
            return;
        }

        // calculate the overall weighting to ensure it remains a constant area

        ScalarType kbkgnd = 0.05, kline=0.15, kgauss = 1.0-kbkgnd-kline; // kbkgnd = 0.1, kline=0.15

        // gaussian values
        ScalarType mean=m_P;
        ScalarType sigma1=2.0, sigma2=2.35; //sigma1=1.4, sigma2=2.0;

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
    void EndostealRampTransform<TScalar>
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
        ScalarType sigma1=2.0, sigma2=2.5; // sigma1=1.4, sigma2=2.0;

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
    void EndostealRampTransform<TScalar>
    ::CalculateWeightingDensity(ArrayType & array, ScalarType startPosition, ScalarType scale) const  {

        int n = array.size(); ScalarType increment = this->m_length / n;

        if(m_P<this->m_start || m_P>this->m_stop || m_E>this->m_stop) { // invalid m_P so set to constant value
            array.Fill(scale/n);
            return;
        }

        // calculate the overall weighting to ensure it remains a constant area

        ScalarType kbkgnd = 0.20, kline=0.05, kgauss = 1.0-kbkgnd-kline;

        // gaussian values
        ScalarType mean1=m_P, mean2=(m_Ecb>m_P) ? m_Ecb : m_P;
        ScalarType sigma1=2.0, sigma2=2.5;//*(m_rampWidth/(m_rampWidth+m_DBThickness)+1); //sigma1=1.5, sigma2=1.6;

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
            ScalarType x = startPosition+i*increment;
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
    void EndostealRampTransform<TScalar>
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
    typename EndostealRampTransform<TScalar>::ScalarType
    EndostealRampTransform<TScalar>
    ::CalculateSigmaAdjustedCorticalDensity(ScalarType peakSampledValue, ScalarType globalSigma) {
        ScalarType backgroundDensity = std::max(m_ST,m_TB);

        // calculate model max for global sigma - temporarily change the sigma parameter
        ScalarType sigma = m_sigma;
        m_sigma = globalSigma;
        ScalarType peakTransformValue = PeakTransformValue();
        m_sigma = sigma;

        // calculate the sigma adjusted density estimate
        ScalarType density = (m_DB -backgroundDensity) * (peakSampledValue - backgroundDensity) / (peakTransformValue - backgroundDensity) + backgroundDensity;
        density = std::max(std::min(density, this->m_MaxDensity), backgroundDensity); // limmit to within reasonable bounds
        return density;
    }

// Compute transformation Jacobian
    template <typename TScalar>
    void EndostealRampTransform<TScalar>
    ::ComputeJacobianWithRespectToParameters(const InputPointType & p, JacobianType & j ) const {
        std::cout<<"TODO - Implement 'EndostealRampTransform::ComputeJacobianWithRespectToParameters'"<<std::endl;
    }

} // namespace

#endif	/* ITKENDOSTEALRAMPTRANSFORMLOCAL_HXX */

