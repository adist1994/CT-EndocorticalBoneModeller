/* 
 * File:   ThreeTierRectangularSpatialObject.h
 * Author: rap58
 *
 * Created on 06 November 2014, 17:20
 */

#ifndef THREETIERRECTANGULARSPATIALOBJECT_H
#define	THREETIERRECTANGULARSPATIALOBJECT_H

#include <vtkKdTreePointLocator.h>
#include <itkSpatialObject.h>
#include <itkAffineTransform.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <itkLinearInterpolateImageFunction.h>

namespace itk
{
/** \class EllipseSpatialObject
 *
 * \brief TODO
 * \ingroup ITKSpatialObjects
 *
 * \wiki
 * \wikiexample{SpatialObjects/EllipseSpatialObject,Ellipse}
 * \endwiki
 */

    template< unsigned int TDimension = 3 >
    class ProfileSpatialObject: public SpatialObject< TDimension >
    {
    public:

        typedef ProfileSpatialObject                         Self;
        typedef double                                       ScalarType;
        typedef SmartPointer< Self >                         Pointer;
        typedef SmartPointer< const Self >                   ConstPointer;
        typedef SpatialObject< TDimension >                  Superclass;
        //typedef SmartPointer< Superclass >                   SuperclassPointer;
        typedef typename Superclass::PointType               PointType;
        //typedef typename Superclass::BoundingBoxType         BoundingBoxType;
        typedef VectorContainer< IdentifierType, PointType > PointContainerType;
        //typedef SmartPointer< PointContainerType >           PointContainerPointer;
        typedef signed short                                            PixelType;
        const static unsigned int                                       Dimension = TDimension;
        typedef itk::Image<PixelType, 3>                        ImageType;
        typedef itk::LinearInterpolateImageFunction< ImageType, float > InterpolatorType;

        typedef FixedArray< double, TDimension > ArrayType;
        itkStaticConstMacro(NumberOfDimension, unsigned int, TDimension);

        itkNewMacro(Self);
        itkTypeMacro(ProfileSpatialObject, SpatialObject);

        /** Set radii via an array of radius values */
        //itkSetMacro(Radius, ArrayType);

        /** Get radii via an array of radius values */
        //itkGetConstReferenceMacro(Radius, ArrayType);

        /** Returns a degree of membership to the object.
         *  That's useful for fuzzy objects. */
        virtual bool ValueAt(const PointType & point, double & value,
                             unsigned int depth = 0,
                             char *name = ITK_NULLPTR) const ITK_OVERRIDE;

        /** Return true if the object provides a method to evaluate the value
         * at the specified point, false otherwise. */
        virtual bool IsEvaluableAt(const PointType & point,
                                   unsigned int depth = 0,
                                   char *name = ITK_NULLPTR) const ITK_OVERRIDE;

        /** Test whether a point is inside or outside the object */
        virtual bool IsInside(const PointType & point,
                              unsigned int depth,
                              char *) const ITK_OVERRIDE;

        /** Test whether a point is inside or outside the object
         *  For computational speed purposes, it is faster if the method does not
         *  check the name of the class and the current depth */
        virtual bool IsInside(const PointType & point) const;

        /** Get the boundaries of a specific object.  This function needs to
         *  be called every time one of the object's components is
         *  changed. */
        virtual bool ComputeLocalBoundingBox() const ITK_OVERRIDE;

        /** Copy the information from another SpatialObject */
        void CopyInformation(const DataObject *data) ITK_OVERRIDE;


        //-------------------- Methods added by Rose ---------------------//
        typedef enum {
            kNoCal=0, // 0.01%
            kLinearCal,    // 0.1%
            kQuadraticCal,   //   5%
        } Calibration;

        ProfileSpatialObject(InterpolatorType::Pointer interpolator, ScalarType length, int numberOfSamples);

        //--- setters - parameters ---//
        void SetInterpolator(InterpolatorType::Pointer interpolator);
        void SetProfileLength(ScalarType profileLength);
        void SetNumberOfSamples(unsigned int numberOfSamples);
        void SetMaximumPeriostealOffset(ScalarType maximumOffset);
        void SetEdgeAndOrientation(ScalarType measurementPoint[NumberOfDimension], ScalarType normal[NumberOfDimension], bool dblPeakDetectionOn=false);
        void SetEdgeRatio(ScalarType ratio);
        void SetPeriostealOffset(ScalarType offset);
        void SetResolution(ScalarType resolution[NumberOfDimension]); // defines spacing of multiple profiles - must be called before SetEdgeAndOrientation
        void SetCalibration(ScalarType p0, ScalarType p1, ScalarType p2);


        //--- setters - behaviour ---//
        void TurnOnCurvedProfiles(vtkSmartPointer<vtkPolyData> mesh, vtkSmartPointer<vtkKdTreePointLocator> meshTree, vtkSmartPointer<vtkDataArray> m_meshNormals);
        void TurnOffCurvedProfiles();
        void TurnOnMultipleProfiles(ScalarType fwhmRadius[NumberOfDimension]); // must be called before SetEdgeAndOrientation
        void TurnOffMultipleProfiles(); // must be called before SetEdgeAndOrientation
        void TurnOnMeanOnly();
        void TurnOffMeanOnly();
        void TurnOnGlobalSigma(ScalarType sigma[NumberOfDimension]);
        void TurnOffGlobalSigma();
        void TurnOnAllSamplesRequired();
        void TurnOffAllSamplesRequired();
        void TurnOnMaxDetection();
        void TurnOffMaxDetection();
        void TurnOnMinDetection();
        void TurnOffMinDetection();

        //---- getters - parameters ---//
        // general
        ScalarType GetProfileLength() const;
        ScalarType GetIncrement() const;
        unsigned int GetNumberOfSamples() const;
        unsigned int GetNumberOfInboundsSamples() const;
        unsigned int GetEdgeIndex() const;
        ScalarType GetPeriostealEdgeRatio() const;
        ScalarType GetMaximumPeriostealOffset() const;
        bool GetRadius(ScalarType &x, ScalarType &y, ScalarType &z); // defines region over which to average the profiles
        bool GetGlobalSigma(ScalarType sigma[NumberOfDimension]) const;
        ScalarType GetProfileSigma() const;
        void GetCalibration(ScalarType &p0, ScalarType &p1, ScalarType &p2) const;
        // specific
        ScalarType GetPeriostealOffset();
        vtkSmartPointer<vtkDoubleArray> GetMultipleProfileStarts() const;
        unsigned int GetNumberOfProfiles() const;
        void GetEdgeAndOrientation(ScalarType measurementPoint[NumberOfDimension], ScalarType normal[NumberOfDimension]) const;
        void GetEdge(ScalarType measurementPoint[NumberOfDimension]) const;
        void GetOrientation(ScalarType normal[NumberOfDimension]) const;
        void GetNormalIncrement(ScalarType normalIncrement[NumberOfDimension]) const;
        void GetStart(ScalarType start[NumberOfDimension]) const;
        void GetEnd(ScalarType end[NumberOfDimension]) const;
        ScalarType GetStartXValue() const;
        ScalarType GetEndXValue() const;
        void GetPointOnProfile(ScalarType positionAlongProfile, ScalarType point[NumberOfDimension]) const;
        // positions
        vtkSmartPointer<vtkDoubleArray> GetPointsOnProfile(ScalarType positionAlongProfile) const;
        vtkSmartPointer<vtkPoints> GetProfilePoints() const;
        itk::Array<double> GetPositions() const;
        itk::Array<double> GetPositions(double offset) const; // return positions with default ratio and specified offset.
        // values
        itk::Array<double> GetValues() const;
        itk::Array2D<double> GetMultipleValues() const;
        bool GetSampledExtentIndices(int& startIndex, int& endIndex) const;
        int GetStartIndex();
        bool GetSampledExtents(double& start, double& end) const;
        bool GetSampledMax(int& maxIndex, ScalarType& maxValue) const;
        itk::Array<int> GetMaxIndices() const;
        itk::Array<double> GetMaxValues() const;
        bool GetSampledMin(int& minIndex, ScalarType& minValue) const;
        itk::Array<int> GetMinIndices() const;
        itk::Array<double> GetMinValues() const;
        itk::Array<double> GetOffsetProfileImageValues(ScalarType positionAlongProfile);
        int GetMaximumPossibleProfiles();
        bool GetPeriostealMeshNode(double periostealPt[Dimension], double periostealOffset);



        //--- getters - behaviour ---//
        bool IsProfileAveragingOn() const;
        bool IsSigmaSet() const;
        bool IsPtInsideImage(ScalarType point[3]);
        bool IsPtInsideImage(); // edge point


        //--- modification time ---//
        TimeStamp GetPositionModificationTime(void) const;
        TimeStamp GetSizeModificationTime(void) const;
        ModifiedTimeType GetMTime() const ITK_OVERRIDE;

    protected:
        ProfileSpatialObject(const Self &); //purposely not implemented
        void operator=(const Self &);       //purposely not implemented

        ProfileSpatialObject(void);
        ~ProfileSpatialObject(void);
    private:
        void ResampleProfile();
        void CalculateMultipleStarts(ScalarType start[Dimension], ScalarType normal[Dimension]);
        void CalculateProfileSigma();
        void CalculateCurvedProfilePoints();
        void CalculateNearestCurvedProfilePoints();
        void CalculateStraightProfilePoints();
        void CalculatePositions();
        bool CalculateProfileValues();
        bool CalculateProfileValuesWithMax();
        bool CalculateMultipleProfileValuesWithMax();
        bool CalculateMultipleCurvedProfileValuesWithMax();
        bool CalculateMultipleProfileValuesWithMaxAndMin();
        bool CalculateMultipleCurvedProfileValuesWithMaxAndMin();
        bool CalculateMeanProfileValues();
        bool CalculateCurvedMeanProfileValues();
        bool CalculateMeanProfileValuesWithMax();
        bool CalculateCurvedMeanProfileValuesWithMax();
        bool CalculateDoublePeakEdgeLocations(); // todo - dbl peak - initialise & use

        ScalarType SafeAcos(double x);

        // general profile parameters
        ScalarType m_EdgeRatio; //ratio of (distance to inside) / (total distance) about the periosteal edge
        ScalarType m_ProfileLength;
        ScalarType m_MaxPeriostealOffset;
        unsigned int m_NumberOfSamples;

        // profile specific parameters
        ScalarType m_PeriostealOffset; // offset between the mesh and periosteal edges - + is an inward direction
        double m_MeshEdge[NumberOfDimension]; // note - mesh is assumed to be located on the periosteal surface.
        double m_Normal[NumberOfDimension]; // note - define pointing into the mesh
        double m_NormalIncrement[NumberOfDimension];
        double m_Inside[NumberOfDimension]; // inside point
        double m_Outside[NumberOfDimension]; // outside point

        // profile sample point locations
        bool m_CurvedProfileMode;
        vtkSmartPointer<vtkPoints> m_ProfilePoints; // points along distance function field
        itk::Array<double> m_Positions;
        vtkSmartPointer<vtkPoints> m_ProfileNormals; // normals along curved profiles
        vtkSmartPointer<vtkKdTreePointLocator> m_MeshTree;
        vtkSmartPointer<vtkDataArray> m_MeshNormals; // note - defined pointing out of the mesh *oposite to m_Normal*
        vtkSmartPointer<vtkPolyData> m_Mesh;

        // profile image values
        bool m_AllSamplesRequired, m_RecordMaxValues, m_RecordMinValues;
        InterpolatorType::Pointer m_Interpolator;
        itk::Array<double> m_ProfileValues;
        itk::Array2D<double> m_MultipleProfileValues;
        int m_StartIndex, m_EndIndex, m_MaxIndex, m_MinIndex;
        ScalarType m_MaxValue, m_MinValue;
        itk::Array<double> m_MaxValues;
        itk::Array<int> m_MaxIndices;
        itk::Array<double> m_MinValues;
        itk::Array<int> m_MinIndices;

        itk::Array<double> m_PkEdgeLocations; // todo - dbl peak - initialise & use

        // calibration
        ScalarType m_P0, m_P1, m_P2;
        int m_CalibrationType;

        // HR - multiple profiles - variables
        double m_radius[NumberOfDimension]; // defines extent over-which profiles are defined
        double m_resolution[NumberOfDimension]; // defines the spacing of the multiple profiles
        vtkSmartPointer<vtkDoubleArray> m_Weights;
        vtkSmartPointer<vtkDoubleArray> m_Starts;
        bool m_multipleProfilesSet, m_meanOnlySet;
        int m_numberProfiles;

        // sigma specific values
        bool m_SigmaSet;
        double m_sigma[NumberOfDimension];
        ScalarType m_ProfileSigma;

        // dbl peak specific values
        bool m_DblPeakDetectionOn;
        bool m_DblPeakDetected;

        static constexpr ScalarType MINDOUBLEPRECISION = 0.0000001; // minimun angle value


        TimeStamp  m_PositionTimeStamp;
        TimeStamp  m_SizeTimeStamp;

        /** Print the object informations in a stream. */
        virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;
    };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkProfileSpatialObjectLocal.hxx"
#endif

#endif	/* THREETIERRECTANGULARSPATIALOBJECT_H */

