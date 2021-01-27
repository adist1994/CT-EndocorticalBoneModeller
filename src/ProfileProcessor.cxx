//
// Created by rap58 on 20/10/15.
//

#include "ProfileProcessor.h"


// includes
#include <iostream>
#include <math.h> 
#include <istream>
#include <sstream>
#include <vtkTable.h>

#include "corticalbone.h"
#include <vtkSortDataArray.h>

namespace itk
{
    ProfileProcessor::ProfileProcessor() {
        m_readyToProcess = false;

        m_profile = NULL;

        m_ProcessedTimeStamp.Modified();
    }

    void ProfileProcessor::SetProfile(ProfileType::Pointer profile) {

        if(profile != m_profile || m_profile->GetMTime() < m_ProcessedTimeStamp.GetMTime()) {
            m_profile = profile;

            UpdateState();
        }
    }

    bool ProfileProcessor::IsPtInsideImage(bool includingTheStart) const {
        ScalarType insidePoint[Dimension]; m_profile->GetEnd(insidePoint);

        if(includingTheStart) {
            ScalarType outsidePoint[Dimension]; m_profile->GetStart(outsidePoint);
            return m_profile->IsPtInsideImage(outsidePoint) && m_profile->IsPtInsideImage(insidePoint);
        } else {
            ScalarType measurementPoint[Dimension]; m_profile->GetEdge(measurementPoint);
            return m_profile->IsPtInsideImage(measurementPoint) && m_profile->IsPtInsideImage(insidePoint);
        }
    }

    ProfileProcessor::itkArrayType ProfileProcessor::GetPositions() {

        return m_profile->GetPositions();
    }

    ProfileProcessor::itkArrayType ProfileProcessor::GetImageSamples() {

        return m_profile->GetValues();
    }

    void ProfileProcessor::UpdateState() {
        if(m_profile) {// && !isnan(m_stThreshold) && !isnan(m_cbThreshold) && !isnan(m_tbThreshold)) {
            m_readyToProcess = true;
        }
    }

    //void ProfileProcessor::GenerateData() {    Process(); }
}