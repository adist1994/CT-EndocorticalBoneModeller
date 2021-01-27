/* 
 * File:   itkModelTransformBaseLocal2.hxx
 * Author: rap58
 *
 * Created on 10 November 2014, 17:29
 */

#ifndef ITKMODELTRANSFORMBASELOCAL2_HXX
#define	ITKMODELTRANSFORMBASELOCAL2_HXX

#include "itkNumericTraits.h"
#include "itkModelTransformBaseLocal.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "itkMath.h"
#include "itkCrossHelper.h"

namespace itk
{

    // define constexpr
    template <typename TScalar, unsigned int NSpaceDimensions>
    constexpr const typename ModelTransformBase<TScalar, NSpaceDimensions>::ScalarType ModelTransformBase<TScalar, NSpaceDimensions>
            ::m_MaxErrorMultiplier;

    //----------------- Model Methods -------------//
// Constructor with default arguments
    // Constructor with default arguments
    template <typename TScalar, unsigned int NSpaceDimensions>
    ModelTransformBase<TScalar, NSpaceDimensions>
    ::ModelTransformBase(unsigned int combinedNumberOfParameters) :
            Superclass(combinedNumberOfParameters) {

        m_TimeStamp.Modified();

        m_length = nan("1");
        m_start = nan("1"); m_stop = nan("1");
        m_maxProfileOffset = nan("1");  m_profileEdgeRatio = nan("1");

        m_numberOfFlexibleParameters = combinedNumberOfParameters;

        m_combinedParameterMapping.SetSize(combinedNumberOfParameters);
        m_combinedParameters.SetSize(combinedNumberOfParameters);

        m_combinedParameterMapping.Fill(1.0);

        m_ParameterErrorScales.SetSize(combinedNumberOfParameters+1);
        m_ParameterErrorScales.Fill(1.0);

        // overwrite fixed param size
        //this->m_FixedParameters.SetSize(0);

        // initialise static variables
        //ScalarType m_MaxDensity;

    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    ModelTransformBase<TScalar, NSpaceDimensions>
    ::ModelTransformBase() : Superclass(){//ParametersDimension) {

        m_TimeStamp.Modified(); m_ProfileTimeStamp.Modified();

        m_length = nan("1");
        m_start = nan("1"); m_stop = nan("1");
        m_maxProfileOffset = nan("1"); m_profileEdgeRatio = nan("1");

        m_numberOfFlexibleParameters = -1;

        m_combinedParameterMapping.SetSize(1);
        m_combinedParameters.SetSize(1);

        m_combinedParameterMapping.Fill(1.0);

        m_ParameterErrorScales.SetSize(1+1);
        m_ParameterErrorScales.Fill(1.0);

    }

    // Destructor
    template <typename TScalar, unsigned int NSpaceDimensions>
    ModelTransformBase<TScalar, NSpaceDimensions>
    ::~ModelTransformBase() { return;  }

    // Print self
    template <typename TScalar, unsigned int NSpaceDimensions>
    void ModelTransformBase<TScalar, NSpaceDimensions>
    ::PrintSelf(std::ostream & os, Indent indent) const {
        Superclass::PrintSelf(os, indent);

        os << indent << "ToDo Implement: 'ModelTransformBase::PrintSelf' " << std::endl;
    }

    // Constructor with explicit arguments
    template <typename TScalar, unsigned int NSpaceDimensions>
    void ModelTransformBase<TScalar, NSpaceDimensions>
    ::SetIdentity(void) {
        this->Modified();
    }

    // Transform a point
    template <typename TScalar, unsigned int NSpaceDimensions>
    typename ModelTransformBase<TScalar, NSpaceDimensions>::ScalarType
    ModelTransformBase<TScalar, NSpaceDimensions>
    ::TransformPoint(const ScalarType position) const {
        //std::cout<<"'TransformPoint' method of 'ModelTransformBase' class"<<std::endl;
        return TransformValue(position);
    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    typename ModelTransformBase<TScalar, NSpaceDimensions>::ScalarType
    ModelTransformBase<TScalar, NSpaceDimensions>
    ::GetUntransformedPoint(const ScalarType position) const {
        return UntransformedValue(position);
    }

    // Transform a point
    template <typename TScalar, unsigned int NSpaceDimensions>
    typename ModelTransformBase<TScalar, NSpaceDimensions>::OutputPointType
    ModelTransformBase<TScalar, NSpaceDimensions>
    ::TransformPoint(const InputPointType & point) const {
        std::cout<<"Error - depreciated method: 'ModelTransformBase::TransformPoint(InputPointType)'"<<std::endl;
        return point;
    }

    // Transform a vector
    template <typename TScalar, unsigned int NSpaceDimensions>
    typename ModelTransformBase<TScalar, NSpaceDimensions>::OutputVectorType
    ModelTransformBase<TScalar, NSpaceDimensions>
    ::TransformVector(const InputVectorType & vect) const {
        std::cout<<"Error - depreciated method: 'ModelTransformBase::TransformVector(InputVectorType),"<<std::endl;
        return vect;
    }

    // Transform a vector
    template <typename TScalar, unsigned int NSpaceDimensions>
    typename ModelTransformBase<TScalar, NSpaceDimensions>::ProfilePosistionArrayType
    ModelTransformBase<TScalar, NSpaceDimensions>
    ::TransformVector(const ProfilePosistionArrayType & positions) const {
        //std::cout<<"'TransformVector' method of 'ModelTransformBase' class"<<std::endl;
        Array<TScalar> outputValues(positions.Size());
        for(int i=0; i<positions.Size(); i++) {
            outputValues[i] = TransformValue(positions[i]);
        }
        return outputValues;
    }

    //----------------- State -----------------------//
    template <typename TScalar, unsigned int NSpaceDimensions>
    bool ModelTransformBase<TScalar, NSpaceDimensions>
    ::IsValid() const {
        return m_valid;
    }


    //------------------ Getters and Setters --------------------//
    template <typename TScalar, unsigned int NSpaceDimensions>
    void  ModelTransformBase< TScalar, NSpaceDimensions >
    ::GetWeightingFunction(ProfilePosistionArrayType & array, ScalarType startPosition, int selectionIndex, ScalarType scale) const {

        if(selectionIndex == kHueristicWeighting) {
            CalculateWeightingHeuristic(array, startPosition, scale);
        } else if(selectionIndex == kWeightingFunction1) {
            CalculateWeightingFunction1(array, startPosition, scale);
        } else if(selectionIndex == kWeightingFunction2) {
            CalculateWeightingFunction2(array, startPosition, scale);
        } else if(selectionIndex == kWeightingFunction3) {
            CalculateWeightingFunction3(array, startPosition, scale);
        } else if(selectionIndex == kWeightingFunction4) {
            CalculateWeightingFunction4(array, startPosition, scale);
        } else if(selectionIndex == kDensityWeighting) {
            CalculateWeightingDensity(array, startPosition, scale);
        } else {
            cerr<<"Error invalid Weighting Function index entered: "<<selectionIndex<<endl;

        }

    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    TimeStamp  ModelTransformBase< TScalar, NSpaceDimensions >
    ::GetProfileModificationTime() const {
        return this->m_ProfileTimeStamp;
    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    void  ModelTransformBase<TScalar, NSpaceDimensions>
    ::SetProfileLength(ScalarType length) {
        m_length = length;
        m_stop = m_length*2.0/3.0; m_start = -m_length*1.0/3.0;
        m_maxProfileOffset = nan("1"); m_profileEdgeRatio = nan("1");
        m_ProfileTimeStamp.Modified();
    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    bool  ModelTransformBase<TScalar, NSpaceDimensions>
    ::setStart(ScalarType start) {

        m_start = start; m_ProfileTimeStamp.Modified();
        if(!isnan(start)) {
            return true;
        } else {
            return false;
        }
    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    bool  ModelTransformBase<TScalar, NSpaceDimensions>
    ::setStop(ScalarType stop) {

        m_stop = stop; m_ProfileTimeStamp.Modified();
        if(!isnan(m_stop)) {
            return true;
        } else {
            return false;
        }
    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    void  ModelTransformBase<TScalar, NSpaceDimensions>
    ::SetProfileEdgeProperties(ScalarType ratio, ScalarType maxOffset) {

        // set ratio
        if(isnan(ratio) || isinf(ratio)) { ratio = 0.5;}
        ratio = ratio<0.0 ? 0.0 : ratio;
        ratio = ratio>1.0 ? 1.0 : ratio;
        m_profileEdgeRatio = ratio;

        // set maximum allowable offset
        if(isnan(maxOffset) || isinf(maxOffset)) { maxOffset = 0.0;}
        m_maxProfileOffset = maxOffset;

    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    void ModelTransformBase<TScalar, NSpaceDimensions>
    ::SetMaxDensity(ScalarType maxDensity) {
        m_MaxDensity = maxDensity;
    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    typename ModelTransformBase<TScalar, NSpaceDimensions>::ScalarType
    ModelTransformBase<TScalar, NSpaceDimensions>
    ::GetMaxDensity() {
        return m_MaxDensity;
    }

    // Set fixed parameters
    template <typename TScalar, unsigned int NSpaceDimensions>
    void  ModelTransformBase<TScalar, NSpaceDimensions>
    ::SetFixedParameters(const ParametersType & fp) {
        if( fp.Size() < m_numberOfFlexibleParameters ) { // todo condiser role of fixed parms
            std::cout<< "Error 'SetFixedParameters' method of 'ModelTransformBase'. Incorrect number of parameters"<<std::endl;
            std::cerr<< "Error 'SetFixedParameters' method of 'ModelTransformBase'. Incorrect number of parameters"<<std::endl;
        }
        this->m_FixedParameters = fp;
    }

    // Set parameters
    template <typename TScalar, unsigned int NSpaceDimensions>
    void ModelTransformBase<TScalar, NSpaceDimensions>
    ::SetParameters(const ParametersType & parameters) {

        if( parameters.Size() < m_numberOfFlexibleParameters ) {
            std::cout<< "Error 'SetParameters' method of 'ModelTransformBase'. Incorrect number of parameters"<<std::endl;
            std::cerr<< "Error 'SetParameters' method of 'ModelTransformBase'. Incorrect number of parameters"<<std::endl;
        }

        // Save parameters. Needed for proper operation of TransformUpdateParameters.
        if( &parameters != &(this->m_Parameters) ) {
            this->m_Parameters = parameters;
            m_TimeStamp.Modified();

            int flexIndex = 0;
            for(int i=0; i<this->m_combinedParameters.GetSize(); i++) {
                if(this->m_combinedParameterMapping[i]==1) {

                    this->m_combinedParameters[i] = this->m_Parameters[flexIndex];
                    flexIndex++;
                }
            }

            UpdateInternalParameters();
            this->Modified();
        }

        itkDebugMacro(<< "After setting parameters ");
    }

    /** Get the Fixed Parameters. */
    template <typename TScalar, unsigned int NSpaceDimensions>
    const typename ModelTransformBase<TScalar, NSpaceDimensions>::ParametersType
    & ModelTransformBase<TScalar, NSpaceDimensions>
    ::GetFixedParameters(void) const {
        return this->m_FixedParameters;
    }

    // Get parameters
    template <typename TScalar, unsigned int NSpaceDimensions>
    const typename ModelTransformBase<TScalar, NSpaceDimensions>::ParametersType
    & ModelTransformBase<TScalar, NSpaceDimensions>
    ::GetParameters(void) const {

        return this->m_Parameters;
    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    const typename ModelTransformBase<TScalar, NSpaceDimensions>::ParametersType
    & ModelTransformBase<TScalar, NSpaceDimensions>
    ::GetCombinedParameters(void) const {
        return this->m_combinedParameters;
    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    void ModelTransformBase<TScalar, NSpaceDimensions>
    ::SetCombinedParameters(const ParametersType & parameters) { // note hard set combined parameters but do not overwrite the norminal fixed parameters info

        if( parameters.Size() < m_combinedParameters.Size() ) {
            std::cout<< "Error 'setCombinedParameters' method of 'ModelTransformBase'. Incorrect number of parameters"<<std::endl;
            std::cerr<< "Error 'setCombinedParameters' method of 'ModelTransformBase'. Incorrect number of parameters"<<std::endl;
        }

        // Save parameters. Needed for proper operation of TransformUpdateParameters.
        if( &parameters != &(this->m_combinedParameters) ) {
            this->m_combinedParameters = parameters;
            // update parameter values too
            int flexIndex = 0;
            for(int i=0; i<this->m_combinedParameters.GetSize(); i++) {
                if(this->m_combinedParameterMapping[i]==1) {

                    this->m_Parameters[flexIndex] = this->m_combinedParameters[i];
                    flexIndex++;
                }
            }

            m_TimeStamp.Modified();
            UpdateInternalParameters();
            this->Modified();
        }
    }


    template <typename TScalar, unsigned int NSpaceDimensions>
    void ModelTransformBase<TScalar, NSpaceDimensions>
    ::FixParameter(const ScalarType value, const unsigned int index) {

        if(index > m_combinedParameters.GetSize()) { // check for valid index
            std::cerr<< "Error 'FixParameter' method of 'ModelTransformBase'. Index is out of bounds."<<std::endl;
        }

        if(m_combinedParameterMapping[index] == 1.0) { // change from flex to fixed
            m_combinedParameterMapping[index] = 0.0;
            m_numberOfFlexibleParameters--;
        } else {
            //std::cerr<< "Warning 'FixParameter' method of 'ModelTransformBase'. Index is already fixed."<<std::endl;
        }
        m_combinedParameters[index] = value;

        // update parameters the number of flexible parameters has changed
        if( this->m_Parameters.size() != m_numberOfFlexibleParameters ) {
            this->m_Parameters.set_size(m_numberOfFlexibleParameters);
            int flexIndex = 0;
            for(int i=0; i<this->m_combinedParameters.GetSize(); i++) {
                if(this->m_combinedParameterMapping[i]==1) {
                    this->m_Parameters[flexIndex] = this->m_combinedParameters[i];
                    flexIndex++;
                }
            }
        }

    }

    template <typename TScalar, unsigned int NSpaceDimensions>
    void ModelTransformBase<TScalar, NSpaceDimensions>
    ::FreeParameter(const unsigned int index) {

        if(index > m_combinedParameters.GetSize()) { // check for valid index
            std::cerr<< "FreeParameter 'FixParameter' method of 'ModelTransformBase'. Index is out of bounds."<<std::endl;
        }

        if(m_combinedParameterMapping[index] == 1.0) {
            // already free
            //std::cerr<< "Error 'FreeParameter' method of 'ModelTransformBase'. Index is already free."<<std::endl;
        } else {
            m_numberOfFlexibleParameters ++;
            m_combinedParameterMapping[index] = 1.0;
        }

        // update parameters the number of flexible parameters has changed
        if( this->m_Parameters.size() != m_numberOfFlexibleParameters ) {
            this->m_Parameters.set_size(m_numberOfFlexibleParameters);
            int flexIndex = 0;
            for(int i=0; i<this->m_combinedParameters.GetSize(); i++) {
                if(this->m_combinedParameterMapping[i]==1) {
                    this->m_Parameters[flexIndex] = this->m_combinedParameters[i];
                    flexIndex++;
                }
            }
        }
    }


    template <typename TScalar, unsigned int NSpaceDimensions>
    ModifiedTimeType ModelTransformBase<TScalar, NSpaceDimensions>
    ::GetMTime() const {
        return this->m_TimeStamp.GetMTime();
    }

    // ---------------- Processing -------------//
    // Compute the Jacobian in one position
    template <typename TScalar, unsigned int NSpaceDimensions>
    void ModelTransformBase<TScalar, NSpaceDimensions>
    ::ComputeJacobianWithRespectToParameters(const InputPointType & p, JacobianType & jacobian) const {
        std::cout<<"TODO Implement: 'ModelTransformBase::ComputeJacobianWithRespectToParameters'"<<std::endl;
    }

    // Return jacobian with respect to position.
    template <typename TScalar, unsigned int NSpaceDimensions>
    void ModelTransformBase<TScalar, NSpaceDimensions>
    ::ComputeJacobianWithRespectToPosition(const InputPointType  &, JacobianType & jac) const{
        std::cout<<"TODO Implement: 'ModelTransformBase::ComputeJacobianWithRespectToPosition'"<<std::endl;
    }

    // Return jacobian with respect to position.
    template <typename TScalar, unsigned int NSpaceDimensions>
    void ModelTransformBase<TScalar, NSpaceDimensions>
    ::ComputeInverseJacobianWithRespectToPosition(const InputPointType  &,  JacobianType & jac) const{
        std::cout<<"TODO Implement: 'ModelTransformBase::ComputeInverseJacobianWithRespectToPosition'"<<std::endl;
    }

} // namespace

#endif	/* ITKMODELTRANSFORMBASELOCAL2_HXX */

