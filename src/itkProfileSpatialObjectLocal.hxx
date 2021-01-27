/* 
 * File:   ThreeTierRectangularSpatialObject.hxx
 * Author: rap58
 *
 * Created on 06 November 2014, 17:20
 */

#ifndef THREETIERRECTANGULARSPATIALOBJECT_HXX
#define	THREETIERRECTANGULARSPATIALOBJECT_HXX

#include "itkProfileSpatialObjectLocal.h"
#include <linearTransform.h>
#include <cstring>
#include <vtkMath.h>
#include <vtkIdList.h>

namespace itk
{


//-------------------- Base Methods to be ignored---------------//

// Test whether a point is inside or outside the object
//  For computational speed purposes, it is faster if the method does not
//  check the name of the class and the current depth.
    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::IsInside(const PointType & point) const {
        std::cerr<<"Error: 'IsInside' method of 'ThreeTierRectangularSpatialObject' class is not supported"<<std::endl;
        std::cout<<"Error: 'IsInside' method of 'ThreeTierRectangularSpatialObject' class is not supported"<<std::endl;
        return false;
    }

// Test if the given point is inside the ellipse 
    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::IsInside(const PointType & point, unsigned int depth, char *name) const {
        std::cerr<<"Error: 'IsInside' method of 'ThreeTierRectangularSpatialObject' class is not supported"<<std::endl;
        std::cout<<"Error: 'IsInside' method of 'ThreeTierRectangularSpatialObject' class is not supported"<<std::endl;
        return false;
    }

// Compute the bounds of the ellipse
    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::ComputeLocalBoundingBox() const {
        std::cerr<<"Error: 'ComputeLocalBoundingBox' method of 'ThreeTierRectangularSpatialObject' class is not supported"<<std::endl;
        std::cout<<"Error: 'ComputeLocalBoundingBox' method of 'ThreeTierRectangularSpatialObject' class is not supported"<<std::endl;
        return false;
    }

// Returns if the ellipse os evaluable at one point
    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::IsEvaluableAt(const PointType & point, unsigned int depth, char *name) const {
        std::cerr<<"Error: 'IsEvaluableAt' method of 'ThreeTierRectangularSpatialObject' class is not supported"<<std::endl;
        std::cout<<"Error: 'IsEvaluableAt' method of 'ThreeTierRectangularSpatialObject' class is not supported"<<std::endl;
        return false;
    }

// Returns the value at one point
    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::ValueAt(const PointType & point, double & value, unsigned int depth, char *name) const {
        std::cerr<<"Error: 'ValueAt' method of 'ThreeTierRectangularSpatialObject' class is not supported"<<std::endl;
        std::cout<<"Error: 'ValueAt' method of 'ThreeTierRectangularSpatialObject' class is not supported"<<std::endl;
        return false;
    }

//-------------------- Base Methods to be used ------------------//
    // Constructor
    template< unsigned int TDimension >
    ProfileSpatialObject< TDimension >
    ::ProfileSpatialObject() {
        this->SetTypeName("ProfileSpatialObject");
        this->SetDimension(TDimension);
        m_PositionTimeStamp.Modified();
        m_SizeTimeStamp.Modified();

        m_EdgeRatio = 0.5;
        m_ProfileLength = 0;
        m_NumberOfSamples = 0;
        m_PeriostealOffset = 0;
        m_MaxPeriostealOffset = 0.0;

        // sigma values
        m_SigmaSet = false;
        //m_sigma = NULL;
        m_ProfileSigma = nan("1");

        // dbl peak
        m_DblPeakDetectionOn=false;
        m_DblPeakDetected=false;

        // multiple profile values
        m_multipleProfilesSet = m_meanOnlySet = false;
        m_numberProfiles = 1;
        m_Starts = vtkSmartPointer<vtkDoubleArray>::New();
        m_Weights = vtkSmartPointer<vtkDoubleArray>::New();

        m_CurvedProfileMode = false;
        m_ProfilePoints = vtkSmartPointer<vtkPoints>::New();

        // set up sample and position arrays and values
        m_AllSamplesRequired = false; m_RecordMaxValues = m_RecordMinValues = false;

        // default to no calibration
        m_CalibrationType=kNoCal;
        m_P0=0.0; m_P1=1.0; m_P2=0.0;

        // sampling image values
        m_StartIndex = m_EndIndex = m_MaxIndex = m_MinIndex = -1;

    }

    template< unsigned int TDimension >
    ProfileSpatialObject< TDimension >
    ::ProfileSpatialObject(InterpolatorType::Pointer interpolator, ScalarType length, int numberOfSamples) : ProfileSpatialObject() {

        m_ProfileLength = length;
        m_NumberOfSamples = numberOfSamples;

        m_Interpolator = interpolator;

    }

// Destructor 
    template< unsigned int TDimension >
    ProfileSpatialObject< TDimension >
    ::~ProfileSpatialObject() {}

// Print Self function
    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::PrintSelf(std::ostream & os, Indent indent) const {
        Superclass::PrintSelf(os, indent);
        os << "Profile Length: " << m_ProfileLength << std::endl;
        os << "# Of Samples: " << m_NumberOfSamples << std::endl;
        os << "Profile Centre Location: (" << m_MeshEdge[0] <<", " << m_MeshEdge[1] <<", " << m_MeshEdge[2] <<")" << std::endl;
        os << "Profile Normal: " << m_Normal[0] <<", " << m_Normal[1] <<", " << m_Normal[2] <<")" << std::endl;
    }

// Copy the information from another spatial object
    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::CopyInformation(const DataObject *data) {
        // check if we are the same type
        const Self *source = dynamic_cast< const Self * >( data );

        if ( !source )
        {
            std::cout << "CopyInformation: objects are not of the same type" << std::endl;
            return;
        }

        // copy the properties
        Superclass::CopyInformation(data);

        // copy the internal info //TODO - make so compiles get a 'const self' error
        /*ScalarType t1 = source->GetProfileLength();
        unsigned int t2 = source->GetProfileSampleNumber();

        SetProfileLength( t1 );
        SetProfileSampleNumber( source->GetProfileSampleNumber() );
        double location[TDimension], normal[TDimension];
        source->GetCentreAndOrientation(location, normal);
        SetCentreAndOrientation( location, normal );*/
    }

//-------------------- Methods added by Rose ---------------------//
    //--- Setters - Parameters ---//

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::SetInterpolator(InterpolatorType::Pointer interpolator) {
        m_Interpolator =interpolator;
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    :: SetProfileLength(ScalarType profileLength) {

        m_ProfileLength = profileLength;
        m_SizeTimeStamp.Modified();
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    :: SetEdgeRatio(ScalarType edgeRatio) {

        m_EdgeRatio = edgeRatio;
        m_SizeTimeStamp.Modified();
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    :: SetPeriostealOffset(ScalarType offset) {

        offset = (offset>m_MaxPeriostealOffset) ? m_MaxPeriostealOffset : offset;
        offset = (offset<-m_MaxPeriostealOffset) ? -m_MaxPeriostealOffset : offset;

        if(offset == m_PeriostealOffset) {
            return; // already up to date
        }

        m_PeriostealOffset = offset;

        // update the start and end locations
        ScalarType insideOffset[NumberOfDimension];
        ScalarType outsideOffset[NumberOfDimension];

        for(int i=0; i<TDimension; i++) {
            insideOffset[i] = m_Normal[i];
            outsideOffset[i] = -m_Normal[i];
        }

        // calculate the start and end point
        vtkMath::MultiplyScalar(outsideOffset, m_ProfileLength * m_EdgeRatio - offset); // - as how much to offset in outward direct
        vtkMath::MultiplyScalar(insideOffset, m_ProfileLength * (1 - m_EdgeRatio)  + offset);

        vtkMath::Add(m_MeshEdge, outsideOffset, m_Outside); // outside of mesh  - TODO check to ensure this is rigorous
        vtkMath::Add(m_MeshEdge, insideOffset, m_Inside); // inside of mesh

        // note profile positions have been updated
        ResampleProfile();

        m_PositionTimeStamp.Modified();
        m_SizeTimeStamp.Modified();
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::SetNumberOfSamples(unsigned int numberOfSamples) {

        m_NumberOfSamples = numberOfSamples;

        // reset sample array size
        m_Positions=itk::Array<double>(numberOfSamples);

        m_SizeTimeStamp.Modified();
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::SetMaximumPeriostealOffset(ScalarType maximumOffset) {
        if(isnan(maximumOffset) || isinf(maximumOffset)) {
            maximumOffset = 0.0; // ensure a valid value
        }
        m_MaxPeriostealOffset = maximumOffset;
        m_PositionTimeStamp.Modified();
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    :: SetEdgeAndOrientation(ScalarType measurementPoint[NumberOfDimension], ScalarType normal[NumberOfDimension], bool dblPeakDetectionOn) {

        if(m_PeriostealOffset==0 && measurementPoint[0]==m_MeshEdge[0] && measurementPoint[1]==m_MeshEdge[1] && measurementPoint[2]==m_MeshEdge[2]
           && normal[0]==-m_Normal[0] && normal[1]==-m_Normal[1] && normal[2]==-m_Normal[2] && m_SizeTimeStamp>m_PositionTimeStamp) {
            return; // no change needed as nothing changed
        }

        // resets the periosteal offset to the default value of 0
        if(m_PeriostealOffset!=0) {
            m_PeriostealOffset = 0;
            m_SizeTimeStamp.Modified();
        }

        // assumes a normalised normal vector
        ScalarType insideOffset[NumberOfDimension];
        ScalarType outsideOffset[NumberOfDimension];
        for(int i=0; i<TDimension; i++) {
            m_MeshEdge[i] = measurementPoint[i];
            m_Normal[i] = -normal[i]; // assumes input normal to be defined pointing out from the mesh, hence take negative
            m_NormalIncrement[i] = m_Normal[i];
            insideOffset[i] = m_Normal[i];
            outsideOffset[i] = -m_Normal[i];
        }

        // calculate the start and end point
        vtkMath::MultiplyScalar(outsideOffset, m_ProfileLength * m_EdgeRatio - m_PeriostealOffset);
        vtkMath::MultiplyScalar(insideOffset, m_ProfileLength * (1 - m_EdgeRatio) + m_PeriostealOffset);

        vtkMath::Add(m_MeshEdge, outsideOffset, m_Outside); // outside of mesh  - TODO check to ensure this is rigorous
        vtkMath::Add(m_MeshEdge, insideOffset, m_Inside); // inside of mesh

        // calculate the increment vector along the profile
        vtkMath::MultiplyScalar(m_NormalIncrement, m_ProfileLength / (m_NumberOfSamples-1));

        if(this->m_SigmaSet) {
            CalculateProfileSigma();
        }

        // this method only gets called once for each profile - m_DblPeakDetected - is then changed in Resample Profile
        m_DblPeakDetectionOn=m_DblPeakDetected=dblPeakDetectionOn;

        // calculate positions and image samples
        ResampleProfile();

        m_PositionTimeStamp.Modified();

    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::SetResolution(ScalarType resolution[NumberOfDimension]) {
        for(int i=0; i<NumberOfDimension; i++) {
            m_resolution[i] = resolution[i];
        }
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::SetCalibration(ScalarType p0, ScalarType p1, ScalarType p2) {
        m_P0=p0; m_P1=p1; m_P2=p2;
        if(m_P0==0.0 && m_P1==1.0 && m_P2==0.0) {
            m_CalibrationType=kNoCal;
        } else if(m_P2==0.0) {
            m_CalibrationType=kLinearCal;
        } else {
            m_CalibrationType=kQuadraticCal;
        }
    }

    //--- Getters - Parameters ---//
    template< unsigned int TDimension >
    typename ProfileSpatialObject< TDimension >::ScalarType ProfileSpatialObject< TDimension >
    :: GetProfileLength() const {
        return m_ProfileLength;

    }

    template< unsigned int TDimension >
    unsigned int ProfileSpatialObject< TDimension >
    ::GetNumberOfSamples() const {
        return m_NumberOfSamples;

    }

    template< unsigned int TDimension >
    unsigned int ProfileSpatialObject< TDimension >
    ::GetNumberOfInboundsSamples() const {

        if(m_EndIndex==-1 || m_StartIndex==-1) {
            return -1;
        } else {
            return 1 + (m_EndIndex - m_StartIndex);
        }
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    :: GetOrientation(ScalarType normal[NumberOfDimension]) const {
        for(int i=0; i<TDimension; i++) {
            normal[i] = m_Normal[i];
        }
    }

    template< unsigned int TDimension >
    typename ProfileSpatialObject< TDimension >::ScalarType ProfileSpatialObject< TDimension >
    ::GetPeriostealOffset() {
        return m_PeriostealOffset;
    }

    template< unsigned int TDimension >
    typename ProfileSpatialObject< TDimension >::ScalarType ProfileSpatialObject< TDimension >
    :: GetIncrement() const {
        return m_ProfileLength / (m_NumberOfSamples-1);

    }

    template< unsigned int TDimension >
    unsigned int ProfileSpatialObject< TDimension >
    :: GetEdgeIndex() const {
        return m_EdgeRatio * m_NumberOfSamples;
    }

    template< unsigned int TDimension >
    typename ProfileSpatialObject< TDimension >::ScalarType ProfileSpatialObject< TDimension >
    ::GetPeriostealEdgeRatio() const {
        return m_EdgeRatio;
    }

    template< unsigned int TDimension >
    typename ProfileSpatialObject< TDimension >::ScalarType ProfileSpatialObject< TDimension >
    ::GetMaximumPeriostealOffset() const {
        return m_MaxPeriostealOffset;
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::GetEdge(ScalarType measurementPoint[NumberOfDimension]) const {
        for(int i=0; i<TDimension; i++) {
            measurementPoint[i] = m_MeshEdge[i];
        }
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    :: GetNormalIncrement(ScalarType normalIncrement[NumberOfDimension]) const {
        for(int i=0; i<TDimension; i++) {
            normalIncrement[i] = m_NormalIncrement[i];
        }
    }

    template< unsigned int TDimension >
    vtkSmartPointer<vtkPoints> ProfileSpatialObject< TDimension >
    ::GetProfilePoints() const { // get the sampled sample locations

        return m_ProfilePoints;

    }

    template< unsigned int TDimension >
    itk::Array<double> ProfileSpatialObject< TDimension >
    ::GetPositions() const {
        if(m_StartIndex !=-1 && m_EndIndex != -1) {
            return m_Positions;
        } else {
            itk::Array<double> array = itk::Array<double>(m_NumberOfSamples); // todo remove as should be set correctly to nan in sampling method
            array.Fill(nan("1"));
            return array;
        }

    }

    template< unsigned int TDimension >
    itk::Array<double> ProfileSpatialObject< TDimension >
    ::GetPositions(double offset) const {

        itk::Array<double> array = itk::Array<double>(m_NumberOfSamples);
        double start= offset - m_ProfileLength*m_EdgeRatio;
        double increment = GetIncrement();
        for(int i=0; i<m_NumberOfSamples; i++) {
            array[i] = start + i * increment;
        }return array;

    }

    template< unsigned int TDimension >
    itk::Array<double> ProfileSpatialObject< TDimension >
    ::GetValues() const {
        if(m_StartIndex !=-1 && m_EndIndex != -1 && (!m_RecordMaxValues || m_MaxIndex != -1)) {
            return m_ProfileValues;
        } else {
            itk::Array<double> array = itk::Array<double>(m_NumberOfSamples); // todo remove as should be set correctly to nan in sampling method
            array.Fill(nan("1"));
            return array;
        }
    }

    template< unsigned int TDimension >
    itk::Array2D<double> ProfileSpatialObject< TDimension >
    ::GetMultipleValues() const {
        if(m_StartIndex !=-1  && m_EndIndex != -1 && m_multipleProfilesSet) {
            return m_MultipleProfileValues;
        } else {
            itk::Array2D<double> array = itk::Array2D<double>(m_NumberOfSamples,m_numberProfiles); // todo remove as should be set correctly to nan in sampling method
            array.Fill(nan("1"));
            return array;
        }
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::GetSampledExtentIndices(int& startIndex, int& endIndex) const {
        if(m_StartIndex !=-1 || m_EndIndex !=-1) {
            startIndex= m_StartIndex; endIndex= m_EndIndex;
            return true;
        } else{
            return false;
        }
    }

    template< unsigned int TDimension >
    int ProfileSpatialObject< TDimension >
    ::GetStartIndex() {
        return m_StartIndex;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::GetSampledExtents(double& start, double& end) const {
        if(m_StartIndex !=-1 || m_EndIndex !=-1) {
            ScalarType increment = GetIncrement(), startPosition=GetStartXValue();
            start= increment * m_StartIndex + startPosition; end=  increment * m_EndIndex + startPosition;
            return true;
        } else{
            return false;
        }
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::GetSampledMax(int& maxIndex, ScalarType& maxValue) const {
        if(m_MaxIndex!=-1 || isnan(m_MaxValue)) {
            maxIndex= m_MaxIndex; maxValue= m_MaxValue;
            return true;
        } else{
            return false;
        }
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::GetSampledMin(int& maxIndex, ScalarType& maxValue) const {
        if(m_MinIndex!=-1 || isnan(m_MinValue)) {
            maxIndex= m_MinIndex; maxValue= m_MinValue;
            return true;
        } else{
            return false;
        }
    }

    template< unsigned int TDimension >
    itk::Array<int> ProfileSpatialObject< TDimension >
    ::GetMaxIndices() const {
        return m_MaxIndices;
    }

    template< unsigned int TDimension >
    itk::Array<double> ProfileSpatialObject< TDimension >
    ::GetMaxValues() const {
        return m_MaxValues;
    }

    template< unsigned int TDimension >
    itk::Array<int> ProfileSpatialObject< TDimension >
    ::GetMinIndices() const {
        return m_MinIndices;
    }

    template< unsigned int TDimension >
    itk::Array<double> ProfileSpatialObject< TDimension >
    ::GetMinValues() const {
        return m_MinValues;
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::GetStart(ScalarType start[NumberOfDimension]) const {

        for(int i=0; i<TDimension; i++) {
            start[i] = m_Outside[i];
        }
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::GetEnd(ScalarType end[NumberOfDimension]) const {
        for(int i=0; i<TDimension; i++) {
            end[i] = m_Inside[i];
        }
    }
    
    template< unsigned int TDimension >
    typename ProfileSpatialObject< TDimension >::ScalarType ProfileSpatialObject< TDimension >
    ::GetStartXValue() const { // - outside
        return m_PeriostealOffset-m_ProfileLength*m_EdgeRatio;
    }
    
    template< unsigned int TDimension >
    typename ProfileSpatialObject< TDimension >::ScalarType ProfileSpatialObject< TDimension >
    ::GetEndXValue() const { // + inside
        return m_ProfileLength*(1-m_EdgeRatio)+m_PeriostealOffset;
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::GetEdgeAndOrientation(ScalarType measurementPoint[NumberOfDimension], ScalarType normal[NumberOfDimension]) const {
        for(int i=0; i<TDimension; i++) {
            measurementPoint[i] = m_MeshEdge[i];
            normal[i] = m_Normal[i];
        }
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::GetPointOnProfile(ScalarType positionAlongProfile, ScalarType point[NumberOfDimension]) const {
        for(int i=0; i<TDimension; i++) {
            point[i] = m_Outside[i] + m_Normal[i]*positionAlongProfile;
        }
    }

    template< unsigned int TDimension >
    vtkSmartPointer<vtkDoubleArray>  ProfileSpatialObject< TDimension >
    ::GetPointsOnProfile(ScalarType positionAlongProfile) const { // return single point along profile by distance 'positionAlongProfile'
        vtkSmartPointer<vtkDoubleArray> points; points->DeepCopy(m_Starts);


        for(int i=0; i<m_numberProfiles; i++) { // tuple
            for(int j=0; j<TDimension; j++) { // component
                points->SetComponent(i, j, points->GetComponent(i, j) + m_Normal[j]*positionAlongProfile);
            }
        }
        return points;
    }

    template< unsigned int TDimension >
    TimeStamp  ProfileSpatialObject< TDimension >
    ::GetPositionModificationTime() const {
        return this->m_PositionTimeStamp;
    }

    template< unsigned int TDimension >
    TimeStamp  ProfileSpatialObject< TDimension >
    ::GetSizeModificationTime() const {
        return this->m_SizeTimeStamp;
    }

    template< unsigned int TDimension >
    ModifiedTimeType  ProfileSpatialObject< TDimension >
    ::GetMTime() const {
        ModifiedTimeType latestTime = Superclass::GetMTime();

        // gets biggest time = most recently changed
        if(this->m_SizeTimeStamp.GetMTime()>latestTime) {
            latestTime = this->m_SizeTimeStamp.GetMTime();
        } else if(latestTime>this->m_PositionTimeStamp.GetMTime()) {
            latestTime = this->m_PositionTimeStamp.GetMTime();
        }
        return latestTime;
    }

    template< unsigned int TDimension >
    unsigned int ProfileSpatialObject< TDimension >
    ::GetNumberOfProfiles() const {
        return m_numberProfiles;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::GetRadius(ScalarType &x, ScalarType &y, ScalarType &z) {
        if(this->IsProfileAveragingOn()) {
            x=m_radius[0]; y=m_radius[1]; z=m_radius[2];
            return true;
        } else {
            return false;
        }
    }

    template< unsigned int TDimension >
    typename ProfileSpatialObject< TDimension >::ScalarType ProfileSpatialObject< TDimension >
    ::GetProfileSigma() const {
        return m_ProfileSigma;
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::GetCalibration(ScalarType &p0, ScalarType &p1, ScalarType &p2) const {
        p0=m_P0; p1=m_P1; p2=m_P2;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::GetGlobalSigma(ScalarType sigma[NumberOfDimension]) const {
        if(m_SigmaSet) {
            sigma[0] = m_sigma[0]; sigma[1] = m_sigma[1]; sigma[2] = m_sigma[2];
        }
        return m_SigmaSet;
    }

    template< unsigned int TDimension >
    itk::Array<double> ProfileSpatialObject< TDimension >
    ::GetOffsetProfileImageValues(ScalarType overallOffset) { // todo - check works correctly
        // resample image offset to match the imported location


        itk::Array<double> array = itk::Array<double>(m_NumberOfSamples);

        ScalarType increment = GetIncrement();

        for(int i=0; i<m_NumberOfSamples; i++) {

            // calc location along the previoulsy sampled profile. todo - average two normals for curved profiles
            ScalarType offset, normal[Dimension], start[Dimension], samplePt[Dimension];
            int index = (int)(i*increment + overallOffset) / increment;
            if(index<0) { // normal of first two points
                offset = (i*increment + overallOffset);
                ScalarType pt2[Dimension];
                m_ProfilePoints->GetPoint(0, start); m_ProfilePoints->GetPoint(1, pt2);
                vtkMath::Subtract(pt2,start, normal);
                vtkMath::MultiplyScalar(normal, offset/vtkMath::Norm(normal, 3)); // normalise and scale by offset in one go

            } else if(index>=m_NumberOfSamples-1) {
                offset = (i*increment + overallOffset);
                ScalarType pt1[Dimension];
                m_ProfilePoints->GetPoint(m_NumberOfSamples-2, pt1); m_ProfilePoints->GetPoint(m_NumberOfSamples-1, start);
                vtkMath::Subtract(start,pt1, normal);
                vtkMath::MultiplyScalar(normal, offset/vtkMath::Norm(normal, 3)); // normalise and scale by offset in one go
            } else {
                offset = (i*increment + overallOffset);
                ScalarType pt2[Dimension];
                m_ProfilePoints->GetPoint(index, start); m_ProfilePoints->GetPoint(index+1, pt2);
                vtkMath::Subtract(pt2,start, normal);
                vtkMath::MultiplyScalar(normal, offset/vtkMath::Norm(normal, 3)); // normalise and scale by offset in one go
            }

            vtkMath::Add(start, normal, samplePt);

            // sample image
            if(m_multipleProfilesSet) {

                // calculate start points at measurement location
                vtkMath::MultiplyScalar(normal, 1.0/vtkMath::Norm(normal, 3)); // normal now normalised
                CalculateMultipleStarts(samplePt, normal);

                ScalarType value;
                for(int j=0; j<m_numberProfiles; j++) {

                    m_Starts->GetTuple(j, samplePt); // sample pt no longer centre pt

                    // get value
                    if ( this->m_Interpolator->IsInsideBuffer(samplePt)) {
                        value += this->m_Interpolator->Evaluate(samplePt) * m_Weights->GetValue(j);
                    } else {
                        value = nan("1"); break;
                    }
                }
                array[i] = value;

            } else {
                if ( this->m_Interpolator->IsInsideBuffer(samplePt)) { // todo consider including a limit check
                    array[i] = this->m_Interpolator->Evaluate(samplePt);
                } else {
                    array[i] = nan("1");
                }
            }

        }
        return array;
    }

    template< unsigned int TDimension >
    int ProfileSpatialObject< TDimension >
    ::GetMaximumPossibleProfiles() {


        if(!m_multipleProfilesSet) {
            return -1;
        } else if(isnan(m_radius[0])||isnan(m_radius[1])||isnan(m_radius[2])) {
            return -1;
        }

        // calculate the max number of samples in each direction
        double ratio[NumberOfDimension] = {m_radius[0]/m_resolution[0], m_radius[1]/m_resolution[1], m_radius[2]/m_resolution[2]};
        int maxI, maxJ; double resI, resJ;
        if(ratio[0]<ratio[1] && ratio[0]<ratio[2]) {
            maxI=(int)floor(ratio[1]); maxJ=(int)floor(ratio[2]);
            resI = m_resolution[1]; resJ = m_resolution[2];
        } else if(ratio[1]<ratio[0] && ratio[1]<ratio[2]) {
            maxI=(int)floor(ratio[0]); maxJ=(int)floor(ratio[2]);
            resI = m_resolution[0]; resJ = m_resolution[2];
        } else {
            maxI=(int)floor(ratio[0]); maxJ=(int)floor(ratio[1]);
            resI = m_resolution[0]; resJ = m_resolution[1];
        }

        // cycle through positions
        int maxNumberOfProfiles = 0;
        for(int i = -maxI; i<=maxI; i++) {
            for(int j = -maxJ; j<=maxJ; j++) {

                double v[NumberOfDimension];
                double idvi[NumberOfDimension] = {i*resI,0,0};
                double jdvj[NumberOfDimension] = {0,j*resJ,0};

                vtkMath::Add(idvi, jdvj, v); // ensure vector is withing radius for that direction
                double vlength=vtkMath::Norm(v);
                double ijradius[NumberOfDimension]={m_radius[0]*v[0]/vlength,m_radius[1]*v[1]/vlength,0};
                if(vlength <= vtkMath::Norm(ijradius) || vlength == 0) {
                    maxNumberOfProfiles++;
                }
            }
        }
        return maxNumberOfProfiles;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
            ::GetPeriostealMeshNode(double periostealPt[Dimension], double periostealOffset) {
        // todo - add support for calculating the point location along a curved line given an offset

    }


    //-- Setters - Behaviour ---//
    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOnCurvedProfiles(vtkSmartPointer<vtkPolyData> mesh, vtkSmartPointer<vtkKdTreePointLocator> meshTree, vtkSmartPointer<vtkDataArray> meshNormals) {
        m_CurvedProfileMode = true;
        m_Mesh = mesh;
        m_MeshTree=meshTree;
        m_MeshNormals=meshNormals;

        // if current measurement measure the profile points
        if(m_ProfileLength!=0 && m_NumberOfSamples!=0) { // only calculate if properly set up
             CalculateCurvedProfilePoints();
        }
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOffCurvedProfiles() {
        m_CurvedProfileMode = false;
        m_MeshTree=NULL;
        m_MeshNormals=NULL;

        m_ProfilePoints = vtkSmartPointer<vtkPoints>::New();
        m_ProfilePoints->SetNumberOfPoints(2);
        m_ProfilePoints->SetPoint(0, m_Outside[0], m_Outside[1], m_Outside[2]);
        m_ProfilePoints->SetPoint(1, m_Inside[0], m_Inside[1], m_Inside[2]);

        CalculateStraightProfilePoints();
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOnMultipleProfiles(ScalarType fwhmRadius[NumberOfDimension]) {
        m_multipleProfilesSet = true;
        for(int i=0; i<NumberOfDimension; i++) {
            m_radius[i] = fwhmRadius[i];
        }
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOffMultipleProfiles() {
        m_multipleProfilesSet = false;

        m_numberProfiles = 1; m_Starts->Reset(); m_Weights->Reset();

        // default ignore values
        m_MultipleProfileValues=itk::Array2D<double>(1,1); m_MultipleProfileValues.Fill(nan("1"));
        m_MaxIndices=itk::Array<int>(1); m_MaxIndices.Fill(-1);
        m_MaxValues=itk::Array<double>(1); m_MaxValues.Fill(nan("1"));
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOnAllSamplesRequired() {
        m_AllSamplesRequired=true;
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOffAllSamplesRequired() {
        m_AllSamplesRequired=false;
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOnMaxDetection() {
        m_RecordMaxValues = true;

    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOffMaxDetection() {
        m_RecordMaxValues =false;

        // default ignore values
        m_MaxIndex = -1; m_MaxIndices=itk::Array<int>(1); m_MaxIndices.Fill(-1);
        m_MaxValue = nan("1"); m_MaxValues=itk::Array<double>(1); m_MaxValues.Fill(nan("1"));
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOnMinDetection() {
        m_RecordMinValues = true;
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOffMinDetection() {
        m_RecordMinValues =false;

        // default ignore values
        m_MinIndex = -1; m_MinIndices=itk::Array<int>(1); m_MinIndices.Fill(-1);
        m_MinValue = nan("1"); m_MinValues=itk::Array<double>(1); m_MinValues.Fill(nan("1"));
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOnMeanOnly() {
        m_meanOnlySet=true;

        // default ignore values
        m_MultipleProfileValues=itk::Array2D<double>(1,1); m_MultipleProfileValues.Fill(nan("1"));
        m_MaxIndices=itk::Array<int>(1); m_MaxIndices.Fill(-1);
        m_MaxValues=itk::Array<double>(1); m_MaxValues.Fill(nan("1"));
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOffMeanOnly() {
        m_meanOnlySet=false;
    }

    // fixed Sigma
    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOnGlobalSigma(ScalarType sigma[NumberOfDimension]) {
        m_SigmaSet = true;
        m_sigma[0] = sigma[0]; m_sigma[1] = sigma[1]; m_sigma[2] = sigma[2];

    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::TurnOffGlobalSigma() {
        m_SigmaSet = false;
    }

    //--- Getter - Behaviour ---//
    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::IsProfileAveragingOn() const {
        return m_multipleProfilesSet;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::IsSigmaSet() const {
        return m_SigmaSet;
    }

    template< unsigned int TDimension >
    bool  ProfileSpatialObject< TDimension >
    ::IsPtInsideImage(ScalarType point[3]) {
        return this->m_Interpolator->IsInsideBuffer(point);
    }

    template< unsigned int TDimension >
    bool  ProfileSpatialObject< TDimension >
    ::IsPtInsideImage() {
        if(m_StartIndex==-1 || m_EndIndex==-1) {
            return false;
        } else {
            return this->m_Interpolator->IsInsideBuffer(m_MeshEdge);
        }
    }

    //--- Private - Calculations ---//
    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::ResampleProfile() { // todo - dbl pk - add support for edge detection
        CalculatePositions();
        if(m_CurvedProfileMode) { // calculate the profile points
            CalculateCurvedProfilePoints();
            if(m_multipleProfilesSet && !m_meanOnlySet && m_RecordMaxValues && m_RecordMinValues) {
                //cerr<<"Curved Multiple Profiles & Max"<<endl;
                CalculateMultipleCurvedProfileValuesWithMaxAndMin();
            } else if(m_multipleProfilesSet && !m_meanOnlySet && m_RecordMaxValues && !m_RecordMinValues) {
                //cerr<<"Curved Multiple Profiles & Max"<<endl;
                CalculateMultipleCurvedProfileValuesWithMax();
            } else if(m_multipleProfilesSet && m_meanOnlySet && !m_RecordMaxValues && !m_RecordMinValues) {
                //cerr<<"Curved Multiple profiles & No Max & Only Mean"<<endl;
                CalculateCurvedMeanProfileValues();
            } else if(m_multipleProfilesSet && m_meanOnlySet && m_RecordMaxValues && !m_RecordMinValues) {
                //cerr<<"Curved Multiple Profiles & Max & Only Mean"<<endl;
                CalculateCurvedMeanProfileValuesWithMax();
            } else if(!m_multipleProfilesSet && !m_RecordMaxValues && !m_RecordMinValues) {
                //cerr<<"Curved One Profiles & No Max"<<endl;
                CalculateProfileValues();
            } else if(!m_multipleProfilesSet && m_RecordMaxValues && !m_RecordMinValues) {
                //cerr<<"Curved One Profiles & Max"<<endl;
                CalculateProfileValuesWithMax();
            } else {
                cerr<<"Error: Invalid resampling selection in ProfileSpatialObject::ResampleProfile"<<endl;
            }
        } else {
            CalculateStraightProfilePoints();
            if(m_multipleProfilesSet && !m_meanOnlySet && m_RecordMaxValues && m_RecordMinValues) {
                //cerr<<"Straight Multiple Profiles & Max"<<endl;
                CalculateMultipleProfileValuesWithMaxAndMin();
            } else if(m_multipleProfilesSet && !m_meanOnlySet && m_RecordMaxValues && !m_RecordMinValues) {
                //cerr<<"Straight Multiple Profiles & Max"<<endl;
                CalculateMultipleProfileValuesWithMax();
            } else if(m_multipleProfilesSet && m_meanOnlySet && !m_RecordMaxValues && !m_RecordMinValues) {
                //cerr<<"Straight Multiple profiles & No Max & Only Mean"<<endl;
                CalculateMeanProfileValues();
            } else if(m_multipleProfilesSet && m_meanOnlySet && m_RecordMaxValues && !m_RecordMinValues) {
                //cerr<<"Straight Multiple Profiles & Max & Only Mean"<<endl;
                CalculateMeanProfileValuesWithMax();
            } else if(!m_multipleProfilesSet && !m_RecordMaxValues && !m_RecordMinValues) {
                //cerr<<"Straight One Profiles & No Max"<<endl;
                CalculateProfileValues();
            } else if(!m_multipleProfilesSet && m_RecordMaxValues && !m_RecordMinValues) {
                //cerr<<"Straight One Profiles & Max"<<endl;
                CalculateProfileValuesWithMax();
            } else {
                cerr<<"Error: Invalid resampling selection in ProfileSpatialObject::ResampleProfile"<<endl;
            }
        }
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::CalculateMultipleStarts(ScalarType start[Dimension], ScalarType normal[Dimension]) {

        // spacing - set if resolution if not curved, set fixed at max resolution if curved
        if(m_CurvedProfileMode) { // sample at a fixed rate irrespective of orientation
            double spacing=std::min(m_resolution[0],std::min(m_resolution[1],m_resolution[2]));

            // calculate normalised orthogonal vectors in the place perpindicular to the normal
            double dvi[NumberOfDimension]={1.0,0.0,0.0}, dvj[NumberOfDimension];
            if(vtkMath::Dot(dvi, normal) == 1) { // ensure not parallel to begin with
                dvi[0]=0.0; dvi[1]=0.0; dvi[2]=1.0;
            }
            vtkMath::Cross(dvi,normal,dvi); vtkMath::Cross(dvi,normal,dvj);
            vtkMath::MultiplyScalar(dvi, 1.0/vtkMath::Norm(dvi));
            vtkMath::MultiplyScalar(dvj, 1.0/vtkMath::Norm(dvj));


            // calculate the max number of samples in each direction - as could take any oientation later on
            double radius = std::max(m_radius[0],std::max(m_radius[1],m_radius[2]));
            int N = (int)floor(radius/spacing*1.0); // 1.0 = length of dvi or dvj

            // cycle through positions
            vtkSmartPointer<vtkDoubleArray> largeStarts = vtkSmartPointer<vtkDoubleArray>::New();
            largeStarts->SetNumberOfComponents(3);
            largeStarts->SetNumberOfTuples((2*N+1) * (2*N+1));

            m_numberProfiles = 0;
            for(int i = -N; i<=N; i++) {
                for(int j = -N; j<=N; j++) {

                    double v[NumberOfDimension];
                    double idvi[NumberOfDimension] = {i*dvi[0]*spacing,i*dvi[1]*spacing,i*dvi[2]*spacing};
                    double jdvj[NumberOfDimension] = {j*dvj[0]*spacing,j*dvj[1]*spacing,j*dvj[2]*spacing};

                    vtkMath::Add(idvi, jdvj, v); // ensure within radius for that direction
                    double vlength=vtkMath::Norm(v);

                    if(vlength <= radius) { // euclidian distance
                        vtkMath::Add(start, v, v);

                        largeStarts->SetTuple(m_numberProfiles, v);

                        m_numberProfiles++;
                    }
                }
            }

            // copy start point within the radius into the start values
            m_Starts->SetNumberOfComponents(3);
            m_Starts->SetNumberOfTuples(m_numberProfiles);
            m_Weights->SetNumberOfComponents(1);
            m_Weights->SetNumberOfValues(m_numberProfiles);
            for(int i=0; i<m_numberProfiles; i++) {
                m_Starts->SetTuple(i, largeStarts->GetTuple(i));
                m_Weights->SetValue(i, 1.0/m_numberProfiles);
            }

        } else { // sample at the appropiate respolution for the orientation

            // calculate normalised orthogonal vectors in the place perpindicular to the normal
            double dvi[NumberOfDimension]={1.0,0.0,0.0}, dvj[NumberOfDimension];
            if(vtkMath::Dot(dvi, normal) == 1) { // ensure not parallel to begin with
                dvi[0]=0.0; dvi[1]=0.0; dvi[2]=1.0;
            }
            vtkMath::Cross(dvi,normal,dvi); vtkMath::Cross(dvi,normal,dvj);
            vtkMath::MultiplyScalar(dvi, 1.0/vtkMath::Norm(dvi));
            vtkMath::MultiplyScalar(dvj, 1.0/vtkMath::Norm(dvj));


            // calculate the max number of samples in each direction
            double ratio[NumberOfDimension] = {m_radius[0]/m_resolution[0], m_radius[1]/m_resolution[1], m_radius[2]/m_resolution[2]};
            int maxI=(int)floor( std::max(  std::max(  fabs(ratio[0]*dvi[0]), fabs(ratio[1]*dvi[1])), fabs(ratio[2]*dvi[2]))   );
            int maxJ=(int)floor( std::max(  std::max(  fabs(ratio[0]*dvj[0]), fabs(ratio[1]*dvj[1])), fabs(ratio[2]*dvj[2]))   );

            // cycle through positions
            vtkSmartPointer<vtkDoubleArray> largeStarts = vtkSmartPointer<vtkDoubleArray>::New();
            largeStarts->SetNumberOfComponents(3);
            largeStarts->SetNumberOfTuples((2*maxI+1) * (2*maxJ+1));

            m_numberProfiles = 0;
            for(int i = -maxI; i<=maxI; i++) {
                for(int j = -maxJ; j<=maxJ; j++) {

                    double v[NumberOfDimension];
                    double idvi[NumberOfDimension] = {i*dvi[0]*m_resolution[0],i*dvi[1]*m_resolution[1],i*dvi[2]*m_resolution[2]};
                    double jdvj[NumberOfDimension] = {j*dvj[0]*m_resolution[0],j*dvj[1]*m_resolution[1],j*dvj[2]*m_resolution[2]};

                    vtkMath::Add(idvi, jdvj, v); // ensure vector is withing radius for that direction
                    double vlength=vtkMath::Norm(v);
                    double ijradius[NumberOfDimension]={m_radius[0]*v[0]/vlength,m_radius[1]*v[1]/vlength,m_radius[2]*v[2]/vlength};
                    if(vlength <= vtkMath::Norm(ijradius) || vlength == 0) {
                        vtkMath::Add(start, v, v);

                        largeStarts->SetTuple(m_numberProfiles, v);

                        m_numberProfiles++;
                    }
                }
            }

            // copy start point within the radius into the start values
            m_Starts->SetNumberOfComponents(3);
            m_Starts->SetNumberOfTuples(m_numberProfiles);
            m_Weights->SetNumberOfComponents(1);
            m_Weights->SetNumberOfValues(m_numberProfiles);
            for(int i=0; i<m_numberProfiles; i++) {
                m_Starts->SetTuple(i, largeStarts->GetTuple(i));
                m_Weights->SetValue(i, 1.0/m_numberProfiles);
            }
        }

    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::CalculateProfileSigma() {
        double product[NumberOfDimension] = {m_sigma[0]*m_Normal[0], m_sigma[1]*m_Normal[1], m_sigma[2]*m_Normal[2]};
        m_ProfileSigma = vtkMath::Norm(product);
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::CalculateNearestCurvedProfilePoints() {

        m_ProfilePoints = vtkSmartPointer<vtkPoints>::New();
        m_ProfilePoints->SetNumberOfPoints(m_NumberOfSamples);


        if(m_ProfileLength==0 || m_NumberOfSamples==0) {
            return;
        } else if(!m_CurvedProfileMode) {
            cerr<<"Warning ProfileSpatialObject::CalculateCurvedProfilePoints() called when !m_CurvedProfileMode. Instead redirecting call to ProfileSpatialObject::CalculateStraightProfilePoints()"<<endl;
            CalculateStraightProfilePoints();
            return;
        }


        // -------- calculate the outside points: define normal(s) to point outwards

        // point to calc all points
        double pt[3]={m_MeshEdge[0]+m_Normal[0]*m_PeriostealOffset, m_MeshEdge[1]+m_Normal[1]*m_PeriostealOffset, m_MeshEdge[2]+m_Normal[2]*m_PeriostealOffset}; // include any periosteal offset
        double normal[3] = {-m_Normal[0], -m_Normal[1], -m_Normal[2]}; // define pointing outwards
        double inc = GetIncrement(); // set increment to a full iteration distance

        // get initial outside point
        int nOutside = (int)((m_NumberOfSamples-1) * m_EdgeRatio);
        double incOutside = inc * (m_EdgeRatio*(m_NumberOfSamples-1) - ((double)nOutside));
        pt[0] = pt[0]+ incOutside *normal[0]; pt[1] = pt[1]+ incOutside *normal[1]; pt[2] = pt[2]+ incOutside *normal[2];
        m_ProfilePoints->SetPoint(nOutside-1, pt[0], pt[1], pt[2]);

        for(int i=nOutside-2; i>=0; i--) { // -2 as 1st is only a partial inc as sample not nessicarily aligned with the mesh edge
            // iteratively update


            // check closest point - and get normal
            vtkIdType nearestId = m_MeshTree->FindClosestPoint(pt);
            // get nearest normal
            double nearestNormal[3];
            m_MeshNormals->GetTuple(nearestId, nearestNormal); // defined pointing outward
            // calculate the distance weightings
            double nearestPt[3]; m_Mesh->GetPoint(nearestId, nearestPt);
            double distanceVector[3] = {pt[0]-nearestPt[0],pt[1]-nearestPt[1],pt[2]-nearestPt[2]};
            double nearestDistance = vtkMath::Norm(distanceVector);
            double weight =1- inc /(inc+nearestDistance), nearestWeight =1 - weight;




            // calculate new point - and new normal
            double nextNormal[3] = {nearestWeight*nearestNormal[0]+weight*normal[0],nearestWeight*nearestNormal[1]+weight*normal[1],nearestWeight*nearestNormal[2]+weight*normal[2]};
            vtkMath::MultiplyScalar(nextNormal, 1.0/vtkMath::Norm(nextNormal)); // normalise

            double angle = SafeAcos(vtkMath::Dot(normal, nextNormal)); // must be normalised vectors
            double straightDistance;
            if(angle!=0) {
                straightDistance = 2.0 * inc * sin(angle/2.0) / angle;
            } else {
                straightDistance = inc;
            }
            pt[0] = pt[0]+nextNormal[0]*straightDistance; pt[1] = pt[1]+nextNormal[1]*straightDistance; pt[2] = pt[2]+nextNormal[2]*straightDistance;

            // store point and update normal
            m_ProfilePoints->SetPoint(i, pt[0], pt[1], pt[2]);
            normal[0] = nextNormal[0]; normal[1] = nextNormal[1]; normal[2] = nextNormal[2];

        }

        // -------- calculate the inside points: define normal(s) to point outwards

        // point to calc all points
        pt[0]=m_MeshEdge[0]+m_Normal[0]*m_PeriostealOffset; pt[1]=m_MeshEdge[1]+m_Normal[1]*m_PeriostealOffset; pt[2]=m_MeshEdge[2]+m_Normal[2]*m_PeriostealOffset;
        normal[0] = m_Normal[0]; normal[1] = m_Normal[1]; normal[2] = m_Normal[2]; // define pointing inwards

        // loop from mesh edge in
        double incrementInside = GetIncrement() - incOutside;
        pt[0] = pt[0]+incrementInside*normal[0]; pt[1] = pt[1]+incrementInside*normal[1]; pt[2] = pt[2]+incrementInside*normal[2];
        m_ProfilePoints->SetPoint(nOutside, pt[0], pt[1], pt[2]);


        for(int i=nOutside+1; i<m_NumberOfSamples; i++) { // -2 as 1st is only a partial inc as sample not nessicarily aligned with the mesh edge
            // iteratively update

            // check closest point - and get normal
            vtkIdType nearestId = m_MeshTree->FindClosestPoint(pt);
            // get nearest normal
            double nearestNormal[3];
            m_MeshNormals->GetTuple(nearestId, nearestNormal); // defined pointing outward
            vtkMath::MultiplyScalar(nearestNormal,-1); // define pointing inwards
            // calculate the distance weightings
            double nearestPt[3]; m_Mesh->GetPoint(nearestId, nearestPt);
            double distanceVector[3] = {pt[0]-nearestPt[0],pt[1]-nearestPt[1],pt[2]-nearestPt[2]};
            double nearestDistance = vtkMath::Norm(distanceVector);
            double weight =1- inc /(inc+nearestDistance), nearestWeight =1 - weight;

            // create a distance scaled normal
            // calculate new point - and new normal
            double nextNormal[3] = {nearestWeight*nearestNormal[0]+weight*normal[0],nearestWeight*nearestNormal[1]+weight*normal[1],nearestWeight*nearestNormal[2]+weight*normal[2]};
            vtkMath::MultiplyScalar(nextNormal, 1.0/vtkMath::Norm(nextNormal)); // normalise

            double angle = SafeAcos(vtkMath::Dot(normal, nextNormal)); // must be normalised vectors
            double straightDistance;
            if(angle!=0) {
                straightDistance = 2.0 * inc * sin(angle/2.0) / angle;
            } else { // no change in normal
                straightDistance = inc;
            }
            pt[0] = pt[0]+nextNormal[0]*straightDistance; pt[1] = pt[1]+nextNormal[1]*straightDistance; pt[2] = pt[2]+nextNormal[2]*straightDistance;

            // store point and update normal
            m_ProfilePoints->SetPoint(i, pt[0], pt[1], pt[2]);
            normal[0] = nextNormal[0]; normal[1] = nextNormal[1]; normal[2] = nextNormal[2];

        }

        return;
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::CalculateCurvedProfilePoints() {

        // todo - calculate a maximum angle to ensure no possible overlap given multiple profiles

        m_ProfilePoints = vtkSmartPointer<vtkPoints>::New(); m_ProfilePoints->SetNumberOfPoints(m_NumberOfSamples);
        m_ProfileNormals = vtkSmartPointer<vtkPoints>::New(); m_ProfileNormals->SetNumberOfPoints(m_NumberOfSamples);


        if(m_ProfileLength==0 || m_NumberOfSamples==0) {
            return;
        } else if(!m_CurvedProfileMode) {
            cerr<<"Warning ProfileSpatialObject::CalculateCurvedProfilePoints() called when !m_CurvedProfileMode. Instead redirecting call to ProfileSpatialObject::CalculateStraightProfilePoints()"<<endl;
            CalculateStraightProfilePoints();
            return;
        }

        // used by both inside and outside point calculations
        double inc = GetIncrement(); // set increment to a full iteration distance
        double edgePt[3]={m_MeshEdge[0]+m_Normal[0]*m_PeriostealOffset, m_MeshEdge[1]+m_Normal[1]*m_PeriostealOffset, m_MeshEdge[2]+m_Normal[2]*m_PeriostealOffset}; // mesh edge point including any periosteal offset

        // -------- calculate the outside points: define normal(s) to point outwards
        double normal[3] = {-m_Normal[0], -m_Normal[1], -m_Normal[2]}; // define pointing outwards
        int nOutside = (int)((m_NumberOfSamples-1) * m_EdgeRatio);
        double incrementOutside = inc * (m_EdgeRatio*(m_NumberOfSamples-1) - ((double)nOutside));

        // calculate & set the first outside point (not only offset by some portion of an increment)
        double pt[3]={edgePt[0]+ incrementOutside *normal[0], edgePt[1]+ incrementOutside *normal[1], edgePt[2]+
                                                                                                      incrementOutside *normal[2]}; // include any periosteal offset
        m_ProfilePoints->SetPoint(nOutside-1, pt[0], pt[1], pt[2]);
        m_ProfileNormals->SetPoint(nOutside-1, -normal[0], -normal[1], -normal[2]); // define point inwards

        // iteratively calculate the rest of the outside points
        double radius = 3; //incrementOutside;
        for(int i=nOutside-2; i>=0; i--) { // -2 as 1st is only a partial inc as sample not nesessarily aligned with the mesh edge

            // 1st. calculate the normal (weighted average of nearby points & the previous pt)
            // find near surface points
            vtkSmartPointer<vtkIdList> nearPtsList = vtkSmartPointer<vtkIdList>::New();
            m_MeshTree->FindPointsWithinRadius(radius, pt, nearPtsList);
            vtkIdType N = nearPtsList->GetNumberOfIds(); //(vtkIdType)std::min((int)nearPtsList->GetNumberOfIds(),3);
            if(N>0) {
                // calculate the average of the near locations
                double nextNormal[3] = {0, 0, 0};//{normal[0]/inc,normal[1]/inc,normal[2]/inc};
                double distance = 0, distanceMin=DBL_MAX;
                for(vtkIdType j = 0; j < N; j++)
                {
                    vtkIdType index = nearPtsList->GetId(j);
                    double pj[3], nj[3]; // pt j, normal j
                    m_MeshNormals->GetTuple(index, nj);
                    m_Mesh->GetPoint(index, pj);

                    double dj[3] = {pt[0]-pj[0],pt[1]-pj[1],pt[2]-pj[2]};
                    double d = vtkMath::Norm(dj);
                    if(d<distanceMin) {distanceMin=d;}

                    distance += d;
                    nextNormal[0] = nextNormal[0] + nj[0]/d;
                    nextNormal[1] = nextNormal[1] + nj[1]/d;
                    nextNormal[2] = nextNormal[2] + nj[2]/d;
                }
                vtkMath::MultiplyScalar(nextNormal, 1.0/vtkMath::Norm(nextNormal)); // normalise
                distance /= N; // normalise near vectors

                // combine near vectors with previous
                nextNormal[0] = normal[0] / inc + nextNormal[0] / distanceMin;//distance;
                nextNormal[1] = normal[1] / inc + nextNormal[1] / distanceMin;//distance;
                nextNormal[2] = normal[2] / inc + nextNormal[2] / distanceMin;//distance;
                vtkMath::MultiplyScalar(nextNormal, 1.0 / vtkMath::Norm(nextNormal));

                // calculate new point
                double angle = SafeAcos(vtkMath::Dot(normal, nextNormal)); // must be normalised vectors
                double straightDistance;
                if(angle!=0) {
                    straightDistance = 2.0 * inc * sin(angle/2.0) / angle;
                } else {
                    straightDistance = inc;
                }
                double straightNormal[3] = {normal[0]+nextNormal[0], normal[1]+nextNormal[1], normal[2]+nextNormal[2]};
                vtkMath::MultiplyScalar(straightNormal, 1.0/vtkMath::Norm(straightNormal)); // normalise

                pt[0] = pt[0]+straightNormal[0]*straightDistance; pt[1] = pt[1]+straightNormal[1]*straightDistance; pt[2] = pt[2]+straightNormal[2]*straightDistance;

                // update normal
                normal[0] = nextNormal[0]; normal[1] = nextNormal[1]; normal[2] = nextNormal[2];
            } else {
                pt[0] = pt[0] + normal[0] * inc;
                pt[1] = pt[1] + normal[1] * inc;
                pt[2] = pt[2] + normal[2] * inc;
            }
            m_ProfilePoints->SetPoint(i, pt[0], pt[1], pt[2]);
            m_ProfileNormals->SetPoint(i, -normal[0], -normal[1], -normal[2]); // define point inwards

            // update the search distance
            //if(radius<5) {
            //    radius += inc;
            //}


        }

        // -------- calculate the inside points: define normal(s) to point outwards
        normal[0] = m_Normal[0]; normal[1] = m_Normal[1]; normal[2] = m_Normal[2]; // define pointing inwards
        double incrementInside = GetIncrement() - incrementOutside;

        // calculate & set the first inside point (not only offset by some prition of an increment)
        pt[0] = edgePt[0]+incrementInside*normal[0]; pt[1] = edgePt[1]+incrementInside*normal[1]; pt[2] = edgePt[2]+incrementInside*normal[2];
        m_ProfilePoints->SetPoint(nOutside, pt[0], pt[1], pt[2]);
        m_ProfileNormals->SetPoint(nOutside, normal[0], normal[1], normal[2]); // define point inwards

        // iteratively calculate the rest of the outside points
        radius = 6; //incrementInside+0.1;
        for(int i=nOutside+1; i<m_NumberOfSamples; i++) { // -2 as 1st is only a partial inc as sample not nessicarily aligned with the mesh edge

            // 1st. calculate the normal (weighted average of nearby points & the previous pt)
            // find near surface points
            vtkSmartPointer<vtkIdList> nearPtsList = vtkSmartPointer<vtkIdList>::New();
            m_MeshTree->FindPointsWithinRadius(radius, pt, nearPtsList);
            vtkIdType N = nearPtsList->GetNumberOfIds(); //(vtkIdType)std::min((int)nearPtsList->GetNumberOfIds(),3);

            if(N>0) {
                // calculate the average of the near locations
                double nextNormal[3] = {0, 0, 0};//{normal[0]/inc,normal[1]/inc,normal[2]/inc};
                double distance = 0, distanceMin=DBL_MAX;
                for (vtkIdType j = 0; j < N; j++) {
                    vtkIdType index = nearPtsList->GetId(j);
                    double pj[3], nj[3]; // pt j, normal j
                    m_MeshNormals->GetTuple(index, nj);
                    m_Mesh->GetPoint(index, pj);

                    double dj[3] = {pt[0] - pj[0], pt[1] - pj[1], pt[2] - pj[2]};
                    double d = vtkMath::Norm(dj);
                    if(d<distanceMin) {distanceMin=d;}

                    distance += d;
                    nextNormal[0] = nextNormal[0] - nj[0] / d; // subtract to define as pointing inwards
                    nextNormal[1] = nextNormal[1] - nj[1] / d;
                    nextNormal[2] = nextNormal[2] - nj[2] / d;
                }
                vtkMath::MultiplyScalar(nextNormal, 1.0 / vtkMath::Norm(nextNormal));
                distance /= N; // normalise near vectors

                // combine near vectors with previous
                nextNormal[0] = normal[0] / inc + nextNormal[0] / distanceMin;//distance;
                nextNormal[1] = normal[1] / inc + nextNormal[1] / distanceMin;//distance;
                nextNormal[2] = normal[2] / inc + nextNormal[2] / distanceMin;//distance;
                vtkMath::MultiplyScalar(nextNormal,1.0 / vtkMath::Norm(nextNormal));

                // 2nd - calculate new point
                double angle = SafeAcos(vtkMath::Dot(normal, nextNormal)); // must be normalised vectors
                double straightDistance;
                if (angle != 0) {
                    straightDistance = 2.0 * inc * sin(angle / 2.0) / angle;
                } else {
                    straightDistance = inc;
                }
                double straightNormal[3] = {normal[0] + nextNormal[0], normal[1] + nextNormal[1],
                                            normal[2] + nextNormal[2]};
                vtkMath::MultiplyScalar(straightNormal, 1.0 / vtkMath::Norm(straightNormal)); // normalise

                pt[0] = pt[0] + straightNormal[0] * straightDistance;
                pt[1] = pt[1] + straightNormal[1] * straightDistance;
                pt[2] = pt[2] + straightNormal[2] * straightDistance;

                //cout << "i=" << i << ", ni=[" << normal[0] << "," << normal[1] << "," << normal[2] << "], ni+1=[" <<
                //nextNormal[0] << "," << nextNormal[1] << "," << nextNormal[2] << "], N=" << N << ", radius=" <<
                //radius <<endl;

                // store point and update normal
                normal[0] = nextNormal[0];
                normal[1] = nextNormal[1];
                normal[2] = nextNormal[2];
            } else {
                pt[0] = pt[0] + normal[0] * inc;
                pt[1] = pt[1] + normal[1] * inc;
                pt[2] = pt[2] + normal[2] * inc;
            }
            m_ProfilePoints->SetPoint(i, pt[0], pt[1], pt[2]);
            m_ProfileNormals->SetPoint(i, normal[0], normal[1], normal[2]); // define point inwards

            // update the search distance
            //if(radius<5) {
            //    radius += inc;
            //}

        }

        return;
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::CalculateStraightProfilePoints() {

        m_ProfilePoints = vtkSmartPointer<vtkPoints>::New();
        m_ProfilePoints->SetNumberOfPoints(m_NumberOfSamples);

        if(m_ProfileLength==0 || m_NumberOfSamples==0) {
            return;
        } else if(m_CurvedProfileMode) {
            cerr<<"Warning ProfileSpatialObject::CalculateStraightProfilePoints() called when m_CurvedProfileMode. Instead redirecting call to ProfileSpatialObject::CalculateCurvedProfilePoints()"<<endl;
            CalculateCurvedProfilePoints();
            return;
        }

        double pt[3] = {m_Outside[0], m_Outside[1], m_Outside[2]}; // includes any periosteal offset

        for(vtkIdType i=0; i<m_NumberOfSamples; i++) {
            m_ProfilePoints->SetPoint(i, pt);
            vtkMath::Add(pt, m_NormalIncrement, pt);      //pt[0]+=m_NormalIncrement[0]; pt[1]+=m_NormalIncrement[1]; pt[2]+=m_NormalIncrement[2];
        }
    }

    template< unsigned int TDimension >
    void ProfileSpatialObject< TDimension >
    ::CalculatePositions() {

        ScalarType increment = GetIncrement();
        ScalarType startPosition = GetStartXValue();

        // resample the image along the profile
        for (int i = 0; i < m_NumberOfSamples; ++i)  {
            m_Positions[i]= startPosition+i*increment;
        }
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::CalculateProfileValues() { // separate out different calibrations to save on if/else each foor loop iteration

        // reset the size of the arrays storing the profile samples and initialise
        m_ProfileValues=itk::Array<double>(m_NumberOfSamples); m_ProfileValues.Fill(nan("1"));

        // get sample parameters
        ScalarType samplePoint[Dimension];

        // sample the image along the profile
        ScalarType value; m_StartIndex = -1; m_EndIndex = -1;
        if(m_CalibrationType==kNoCal) { // image already in BMD

            for (int i = 0; i < m_NumberOfSamples; ++i)  { // travel along profile

                m_ProfilePoints->GetPoint(i, samplePoint); // get position along profile

                if ( m_Interpolator->IsInsideBuffer(samplePoint)) {
                    value = m_Interpolator->Evaluate(samplePoint); // in BMD

                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; m_EndIndex = i;  // valid profile positions
                } else {
                    value = nan("1");
                }
                m_ProfileValues[i] = value;

            }
        } else if(m_CalibrationType==kLinearCal) { // image values are linearly converted to BMD

            for (int i = 0; i < m_NumberOfSamples; ++i)  { // travel along profile

                m_ProfilePoints->GetPoint(i, samplePoint); // get position along profile

                if ( m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                    value = m_Interpolator->Evaluate(samplePoint);
                    value=m_P0+m_P1*value; // convert to BMD

                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; m_EndIndex = i; // valid profile positions

                } else {
                    value = nan("1");
                }
                m_ProfileValues[i] = value;
            }

        } else if(m_CalibrationType==kQuadraticCal) { // image values are quadratically converted to BMD

            for (int i = 0; i < m_NumberOfSamples; ++i)  { // travel along profile

                m_ProfilePoints->GetPoint(i, samplePoint); // get position along profile

                if ( m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                    value = m_Interpolator->Evaluate(samplePoint);
                    value=m_P0+m_P1*value+m_P2*value*value; // convert to BMD

                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; m_EndIndex = i; // valid profile positions

                } else {
                    value = nan("1");
                }
                m_ProfileValues[i] = value;
            }
        } else {
            cerr<<"Error Invalid calibration value"<<endl;
        }

        // note - that the mutliple profile values only set / accessed unless in multiple profile mode. Ignore otherwise
        return true;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::CalculateProfileValuesWithMax() { // separate out different calibrations to save on if/else each foor loop iteration


        // reset the size of the arrays storing the profile samples and initialise
        m_ProfileValues=itk::Array<double>(m_NumberOfSamples); m_ProfileValues.Fill(nan("1"));

        // get sample parameters
        ScalarType samplePoint[Dimension];

        // sample the image along the profile
        ScalarType value; m_StartIndex = -1; m_EndIndex = -1; m_MaxIndex = 0; m_MaxValue = -DBL_MAX;
        if (m_CalibrationType==kNoCal) { // image already in BMD

            for (int i = 0; i < m_NumberOfSamples; ++i)  { // travel along profile

                m_ProfilePoints->GetPoint(i, samplePoint); // get position along profile

                // get value
                if ( m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds
                    value = m_Interpolator->Evaluate(samplePoint); // in BMD

                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; m_EndIndex = i; // valid profile positions

                    if(value> m_MaxValue) { // update max value and location
                        m_MaxValue = value; m_MaxIndex = i;
                    }

                } else {
                    value = nan("1");
                }
                m_ProfileValues[i] = value;
            }

        } else if(m_CalibrationType==kLinearCal) { // image values are linearly converted to BMD
            for (int i = 0; i < m_NumberOfSamples; ++i)  { // travel along profile

                m_ProfilePoints->GetPoint(i, samplePoint); // get position along profile

                // get value
                if ( m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds
                    value = m_Interpolator->Evaluate(samplePoint);
                    value=m_P0+m_P1*value; // convert to BMD

                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; m_EndIndex = i; // valid profile positions

                    if(value> m_MaxValue) { // update max value and location
                        m_MaxValue = value; m_MaxIndex = i;
                    }

                } else {
                    value = nan("1");
                }
                m_ProfileValues[i] = value;
            }
        } else if(m_CalibrationType==kQuadraticCal) { // image values are quadratically converted to BMD

            for (int i = 0; i < m_NumberOfSamples; ++i)  { // travel along profile

                m_ProfilePoints->GetPoint(i, samplePoint); // get position along profile

                // get value
                if ( m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                    value = m_Interpolator->Evaluate(samplePoint);
                    value=m_P0+m_P1*value+m_P2*value*value; // convert to BMD

                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; m_EndIndex = i; // valid profile positions

                    if(value> m_MaxValue) { // update max value and location
                        m_MaxValue = value; m_MaxIndex = i;
                    }

                } else {
                    value = nan("1");
                }
                m_ProfileValues[i] = value;
            }
        } else {
            cerr<<"Error Invalid calibration value"<<endl;
        }

        // note - that the mutliple profile values only set / accessed unless in multiple profile mode. Ignore otherwise
        return true;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::CalculateMultipleProfileValuesWithMax() {
        // mean values stored in the m_ProfilesArrays

        // get sample parameters
        ScalarType samplePoint[Dimension]; GetStart(samplePoint);

        // arrays to storee values in
        CalculateMultipleStarts(m_Outside, m_Normal);
        vtkSmartPointer<vtkDoubleArray> startPoints = vtkSmartPointer<vtkDoubleArray>::New();
        startPoints->DeepCopy(m_Starts);

        // reset the size of the arrays storing the samples of multiple profiles
        m_ProfileValues=itk::Array<double>(m_NumberOfSamples); m_ProfileValues.Fill(nan("1"));
        m_MultipleProfileValues=itk::Array2D<double>(m_NumberOfSamples, m_numberProfiles);
        m_MaxValues=itk::Array<double>(m_numberProfiles); m_MaxIndices=itk::Array<int>(m_numberProfiles);
        m_MultipleProfileValues.Fill(nan("1")); m_MaxValues.Fill(-DBL_MAX); m_MaxIndices.Fill(-1);

        // sample the image along the profile separated by calibration type
        m_StartIndex = -1; m_EndIndex = m_NumberOfSamples-1; m_MaxIndex = 0; m_MaxValue = -DBL_MAX;
        if(m_CalibrationType==kNoCal) { // image is in BMD
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    startPoints->GetTuple(j, samplePoint); // update to

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint); // image already in BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update overal maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }
        } else if(m_CalibrationType==kLinearCal) {
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    startPoints->GetTuple(j, samplePoint); // update to

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue=m_P0+m_P1*singleValue; // linear conversion to BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update overal maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }

        } else if(m_CalibrationType==kQuadraticCal) {
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    startPoints->GetTuple(j, samplePoint); // update to

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue=m_P0+m_P1*singleValue+m_P2*singleValue*singleValue; // quadratic conversion to BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update overal maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }
        }
        return true;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::CalculateMultipleCurvedProfileValuesWithMax() {
        // mean values stored in the m_ProfilesArrays

        // get sample parameters
        ScalarType profilePoint[Dimension], profileNormal[Dimension];

        // reset the size of the arrays storing the samples of multiple profiles
        CalculateMultipleStarts(m_Outside, m_Normal); // call to calculate m_numberProfiles
        m_ProfileValues=itk::Array<double>(m_NumberOfSamples); m_ProfileValues.Fill(nan("1"));
        m_MultipleProfileValues=itk::Array2D<double>(m_NumberOfSamples, m_numberProfiles);
        m_MaxValues=itk::Array<double>(m_numberProfiles); m_MaxIndices=itk::Array<int>(m_numberProfiles);
        m_MultipleProfileValues.Fill(nan("1")); m_MaxValues.Fill(-DBL_MAX); m_MaxIndices.Fill(-1);

        // sample the image along the profile separated by calibration type
        m_StartIndex = -1; m_EndIndex = m_NumberOfSamples-1; m_MaxIndex = 0; m_MaxValue = -DBL_MAX;
        if(m_CalibrationType==kNoCal) { // image is in BMD
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    double samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // lookup sample point

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint); // image already in BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update overal maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }
        } else if(m_CalibrationType==kLinearCal) {
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    double samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // lookup sample point

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue=m_P0+m_P1*singleValue; // linear conversion to BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update overal maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }

        } else if(m_CalibrationType==kQuadraticCal) {
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    double samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // lookup sample point

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue=m_P0+m_P1*singleValue+m_P2*singleValue*singleValue; // quadratic conversion to BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update overal maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }
        }
        return true;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::CalculateMultipleProfileValuesWithMaxAndMin() {
        // mean values stored in the m_ProfilesArrays

        // get sample parameters
        ScalarType samplePoint[Dimension]; GetStart(samplePoint);

        // arrays to storee values in
        CalculateMultipleStarts(m_Outside, m_Normal);
        vtkSmartPointer<vtkDoubleArray> startPoints = vtkSmartPointer<vtkDoubleArray>::New();
        startPoints->DeepCopy(m_Starts);

        // reset the size of the arrays storing the samples of multiple profiles
        m_ProfileValues=itk::Array<double>(m_NumberOfSamples); m_ProfileValues.Fill(nan("1"));
        m_MultipleProfileValues=itk::Array2D<double>(m_NumberOfSamples, (unsigned int)m_numberProfiles);
        m_MaxValues=itk::Array<double>((unsigned int)m_numberProfiles); m_MaxIndices=itk::Array<int>((unsigned int)m_numberProfiles);
        m_MinValues=itk::Array<double>((unsigned int)m_numberProfiles); m_MinIndices=itk::Array<int>((unsigned int)m_numberProfiles);
        m_MultipleProfileValues.Fill(nan("1"));
        m_MaxValues.Fill(-DBL_MAX); m_MaxIndices.Fill(-1);
        m_MinValues.Fill(DBL_MAX); m_MinIndices.Fill(-1);

        // sample the image along the profile separated by calibration type
        m_StartIndex = -1; m_EndIndex = m_NumberOfSamples-1; m_MaxIndex = 0; m_MaxValue = -DBL_MAX;
        if(m_CalibrationType==kNoCal) { // image is in BMD
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    startPoints->GetTuple(j, samplePoint); // update to

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint); // image already in BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }
                        if(singleValue < m_MinValues[j]) { // update individual profile minimums
                            m_MinValues[j] = singleValue; m_MinIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value > m_MaxValue) { // update overal maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                    if(value < m_MinValue) { // update overall minimum value and index
                        m_MinValue = value; m_MinIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex == m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }
        } else if(m_CalibrationType==kLinearCal) {
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    startPoints->GetTuple(j, samplePoint); // update to

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue=m_P0+m_P1*singleValue; // linear conversion to BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }
                        if(singleValue < m_MinValues[j]) { // update individual profile minimums
                            m_MinValues[j] = singleValue; m_MinIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value > m_MaxValue) { // update overall maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                    if(value < m_MinValue) { // update overall minimum value and index
                        m_MinValue = value; m_MinIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }

        } else if(m_CalibrationType==kQuadraticCal) {
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    startPoints->GetTuple(j, samplePoint); // update to

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue=m_P0+m_P1*singleValue+m_P2*singleValue*singleValue; // quadratic conversion to BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }
                        if(singleValue < m_MinValues[j]) { // update individual profile minimums
                            m_MinValues[j] = singleValue; m_MinIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update overal maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                    if(value < m_MinValue) { // update overall minimum value and index
                        m_MinValue = value; m_MinIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }
        }
        return true;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::CalculateMultipleCurvedProfileValuesWithMaxAndMin() {
        // mean values stored in the m_ProfilesArrays

        // get sample parameters
        ScalarType profilePoint[Dimension], profileNormal[Dimension];

        // reset the size of the arrays storing the samples of multiple profiles
        CalculateMultipleStarts(m_Outside, m_Normal); // call to calculate m_numberProfiles
        m_ProfileValues=itk::Array<double>(m_NumberOfSamples); m_ProfileValues.Fill(nan("1"));
        m_MultipleProfileValues=itk::Array2D<double>(m_NumberOfSamples, m_numberProfiles);
        m_MaxValues=itk::Array<double>(m_numberProfiles); m_MaxIndices=itk::Array<int>(m_numberProfiles);
        m_MinValues=itk::Array<double>(m_numberProfiles); m_MinIndices=itk::Array<int>(m_numberProfiles);
        m_MultipleProfileValues.Fill(nan("1"));
        m_MaxValues.Fill(-DBL_MAX); m_MaxIndices.Fill(-1);
        m_MinValues.Fill(DBL_MAX); m_MinIndices.Fill(-1);

        // sample the image along the profile separated by calibration type
        m_StartIndex = -1; m_EndIndex = m_NumberOfSamples-1;
        m_MaxIndex = 0; m_MaxValue = -DBL_MAX; m_MinIndex = 0; m_MinValue = DBL_MAX;
        if(m_CalibrationType==kNoCal) { // image is in BMD
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    double samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // lookup sample point

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint); // image already in BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }

                        if(singleValue < m_MinValues[j]) { // update individual profile maximums
                            m_MinValues[j] = singleValue; m_MinIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value > m_MaxValue) { // update overall maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                    if(value < m_MinValue) { // update overall maximum value and index
                        m_MinValue = value; m_MinIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex == m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }
        } else if(m_CalibrationType==kLinearCal) {
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    double samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // lookup sample point

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue = m_P0 + m_P1 * singleValue; // linear conversion to BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }
                        if(singleValue < m_MinValues[j]) { // update individual profile maximums
                            m_MinValues[j] = singleValue; m_MinIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value > m_MaxValue) { // update overal maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                    if(value < m_MaxValue) { // update overal maximum value and index
                        m_MinValue = value; m_MinIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }

        } else if(m_CalibrationType==kQuadraticCal) {
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0;
                for(int j=0; j<m_numberProfiles; j++) { // cycle through each profiles

                    double samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // lookup sample point

                    // get value
                    if (!isnan(value) && m_Interpolator->IsInsideBuffer(samplePoint)) {
                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue=m_P0+m_P1*singleValue+m_P2*singleValue*singleValue; // quadratic conversion to BMD

                        value += singleValue * m_Weights->GetValue(j);
                        m_MultipleProfileValues[i][j]= singleValue; // store individual profile values

                        if(singleValue > m_MaxValues[j]) { // update individual profile maximums
                            m_MaxValues[j] = singleValue; m_MaxIndices[j] = i;
                        }
                        if(singleValue < m_MinValues[j]) { // update individual profile maximums
                            m_MinValues[j] = singleValue; m_MinIndices[j] = i;
                        }

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_MultipleProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex == m_NumberOfSamples-1) {
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value > m_MaxValue) { // update overal maximum value and index
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                    if(value < m_MinValue) { // update overal maximum value and index
                        m_MinValue = value; m_MinIndex = i;
                    }

                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1; return true;
                }
            }
        }
        return true;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::CalculateMeanProfileValues() {

        // reset sample array
        m_ProfileValues=itk::Array<double>(m_NumberOfSamples); m_ProfileValues.Fill(nan("1"));

        // get sample parameters
        ScalarType samplePoint[Dimension]; GetStart(samplePoint);

        CalculateMultipleStarts(m_Outside, m_Normal);
        vtkSmartPointer<vtkDoubleArray> startPoints = vtkSmartPointer<vtkDoubleArray>::New();
        startPoints->DeepCopy(m_Starts);

        // separately calculate the mean along the profile depending on calibration type
        m_StartIndex = -1; m_EndIndex = m_NumberOfSamples-1;
        if(m_CalibrationType==kNoCal) { // image in BMD
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    startPoints->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint); // in BMD
                        value += singleValue * m_Weights->GetValue(j);

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = -1; m_EndIndex = -1;
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1;
                    return true;
                }
            }
        } else if(m_CalibrationType==kLinearCal) {

            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    startPoints->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        value += (m_P0+m_P1*singleValue) * m_Weights->GetValue(j); // linear conversion to BMD

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = -1; m_EndIndex = -1;
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // invalid value marking end of sampling
                    m_EndIndex = i-1;
                    return true;
                } else { // keep going looking for the 1st in bounds value
                }
            }
        } else if(m_CalibrationType==kQuadraticCal) {

            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    startPoints->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue = m_P0+m_P1*singleValue+m_P2*singleValue*singleValue;
                        value += (singleValue) * m_Weights->GetValue(j); // quadratic conversion to BMD

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = -1; m_EndIndex = -1;
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // invalid value marking end of sampling
                    m_EndIndex = i-1;
                    return true;
                } else { // keep going looking for the 1st in bounds value
                }
            }
        }
        return true;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::CalculateCurvedMeanProfileValues() {
        // reset sample array
        m_ProfileValues=itk::Array<double>(m_NumberOfSamples); m_ProfileValues.Fill(nan("1"));

        // get sample parameters
        ScalarType profilePoint[Dimension], profileNormal[Dimension];

        // separately calculate the mean along the profile depending on calibration type
        m_StartIndex = -1; m_EndIndex = m_NumberOfSamples-1;
        if(m_CalibrationType==kNoCal) { // image in BMD
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    double samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint); // in BMD
                        value += singleValue * m_Weights->GetValue(j);

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = -1; m_EndIndex = -1;
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // end is out of bounds. stop
                    m_EndIndex = i-1;
                    return true;
                }
            }
        } else if(m_CalibrationType==kLinearCal) {

            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    double samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        value += (m_P0+m_P1*singleValue) * m_Weights->GetValue(j); // linear conversion to BMD

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = -1; m_EndIndex = -1;
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // invalid value marking end of sampling
                    m_EndIndex = i-1;
                    return true;
                } else { // keep going looking for the 1st in bounds value
                }
            }
        } else if(m_CalibrationType==kQuadraticCal) {

            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    double samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue = m_P0+m_P1*singleValue+m_P2*singleValue*singleValue;
                        value += (singleValue) * m_Weights->GetValue(j); // quadratic conversion to BMD

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = -1; m_EndIndex = -1;
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // invalid value marking end of sampling
                    m_EndIndex = i-1;
                    return true;
                } else { // keep going looking for the 1st in bounds value
                }
            }
        }
        return true;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::CalculateMeanProfileValuesWithMax() {
        // currently only take samples if all profiles are inside the image boundary.

        // reset sample array
        m_ProfileValues=itk::Array<double>(m_NumberOfSamples); m_ProfileValues.Fill(nan("1"));

        // get sample parameters
        ScalarType samplePoint[Dimension]; GetStart(samplePoint);

        CalculateMultipleStarts(m_Outside, m_Normal);
        vtkSmartPointer<vtkDoubleArray> startPoints = vtkSmartPointer<vtkDoubleArray>::New();
        startPoints->DeepCopy(m_Starts);

        // separately calculate the mean along the profile depending on calibration type
        m_StartIndex = -1; m_EndIndex = m_NumberOfSamples-1; m_MaxIndex = -1; m_MaxValue = -DBL_MAX;
        if(m_CalibrationType==kNoCal) { // image in BMD
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    startPoints->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint); // in BMD
                        value += singleValue * m_Weights->GetValue(j);

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update max value and location
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // invalid value marking end of sampling
                    m_EndIndex = i-1;
                    return true;
                } else { // keep going looking for the 1st in bounds value
                }
            }
        } else if(m_CalibrationType==kLinearCal) {

            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    startPoints->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        value += (m_P0+m_P1*singleValue) * m_Weights->GetValue(j); // linear conversion to BMD

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update max value and location
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // invalid value marking end of sampling
                    m_EndIndex = i-1;
                    return true;
                } else { // keep going looking for the 1st in bounds value
                }
            }
        } else if(m_CalibrationType==kQuadraticCal) {

            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    startPoints->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue = m_P0+m_P1*singleValue+m_P2*singleValue*singleValue;
                        value += (singleValue) * m_Weights->GetValue(j); // quadratic conversion to BMD

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }

                    // update position along profile AND iterate previous location
                    vtkMath::Add(samplePoint, m_NormalIncrement, samplePoint); startPoints->SetTuple(j, samplePoint);
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update max value and location
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // invalid value marking end of sampling
                    m_EndIndex = i-1;
                    return true;
                } else { // keep going looking for the 1st in bounds value
                }
            }
        }
        return true;
    }

    template< unsigned int TDimension >
    bool ProfileSpatialObject< TDimension >
    ::CalculateCurvedMeanProfileValuesWithMax() {
        // currently only take samples if all profiles are inside the image boundary.

        // reset sample array
        m_ProfileValues=itk::Array<double>(m_NumberOfSamples); m_ProfileValues.Fill(nan("1"));

        // get sample parameters
        ScalarType profilePoint[Dimension], profileNormal[Dimension];

        // separately calculate the mean along the profile depending on calibration type
        m_StartIndex = -1; m_EndIndex = m_NumberOfSamples-1; m_MaxIndex = -1; m_MaxValue = -DBL_MAX;
        if(m_CalibrationType==kNoCal) { // image in BMD
            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    ScalarType samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint); // in BMD
                        value += singleValue * m_Weights->GetValue(j);

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update max value and location
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // invalid value marking end of sampling
                    m_EndIndex = i-1;
                    return true;
                } else { // keep going looking for the 1st in bounds value
                }
            }
        } else if(m_CalibrationType==kLinearCal) {

            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    ScalarType samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        value += (m_P0+m_P1*singleValue) * m_Weights->GetValue(j); // linear conversion to BMD


                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update max value and location
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // invalid value marking end of sampling
                    m_EndIndex = i-1;
                    return true;
                } else { // keep going looking for the 1st in bounds value
                }
            }
        } else if(m_CalibrationType==kQuadraticCal) {

            for (int i = 0; i < m_NumberOfSamples; ++i)  {

                m_ProfilePoints->GetPoint(i, profilePoint); m_ProfileNormals->GetPoint(i, profileNormal);
                CalculateMultipleStarts(profilePoint, profileNormal);

                ScalarType value = 0; // separate by calibration type then cycle through values
                for(int j=0; j<m_numberProfiles; j++) {

                    ScalarType samplePoint[Dimension];
                    m_Starts->GetTuple(j, samplePoint); // update measurement location

                    if (m_Interpolator->IsInsideBuffer(samplePoint)) { // check in bounds

                        ScalarType singleValue = m_Interpolator->Evaluate(samplePoint);
                        singleValue = m_P0+m_P1*singleValue+m_P2*singleValue*singleValue;
                        value += (singleValue) * m_Weights->GetValue(j); // quadratic conversion to BMD

                    } else if(m_AllSamplesRequired) { // reset samples and exit
                        m_ProfileValues.Fill(nan("1")); m_StartIndex = m_EndIndex = m_MaxIndex = -1; m_MaxValue = nan("1");
                        return false;
                    } else { // do nothing as not all are required to be in bounds
                        value=nan("1");
                    }
                }

                if(!isnan(value) && m_EndIndex==m_NumberOfSamples-1) { // value to store
                    m_ProfileValues[i] = value;
                    m_StartIndex = (m_StartIndex >=0) ? m_StartIndex : i; // update valid profile positions

                    if(value> m_MaxValue) { // update max value and location
                        m_MaxValue = value; m_MaxIndex = i;
                    }
                } else if(m_StartIndex >=0 && m_EndIndex ==m_NumberOfSamples-1) { // invalid value marking end of sampling
                    m_EndIndex = i-1;
                    return true;
                } else { // keep going looking for the 1st in bounds value
                }
            }
        }
        return true;
    }

    //--- Private - Utility ---//
    template< unsigned int TDimension >
    typename ProfileSpatialObject< TDimension >::ScalarType ProfileSpatialObject< TDimension >
    ::SafeAcos (double x) {
        if (x < -1.0) x = -1.0 ;
        else if (x > 1.0) x = 1.0 ;
        return acos (x) ;
    }


} // end namespace itk

#endif	/* THREETIERRECTANGULARSPATIALOBJECT_HXX */

