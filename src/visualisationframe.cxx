/* 
 * Render the DICOM, Mesh with a 3Dim and slice view and profile view.
 * Notes  - QVTKWidget <- vtkRender <- actor <- vtkObject 
 *        - rest the camera view betweenloading files?  */


#include <QtWidgets>
#include <iostream>
#include <istream>
#include <sstream>

#include <QVTKWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <QVTKInteractor.h>

// render mesh
#include "vtkPolyDataMapper.h"
// load dicom image
#include "vtkImageData.h"
// render image slices
#include <vtkRenderWindow.h>
#include "vtkResliceImageViewer.h"
#include "vtkResliceCursorLineRepresentation.h"
#include "vtkResliceCursorWidget.h"
#include "vtkResliceCursorActor.h"
#include "vtkResliceCursorPolyDataAlgorithm.h"
// render 3d view 
#include <vtkRenderer.h>
#include <vtkRendererCollection.h> 
#include "vtkCellPicker.h"
#include "vtkLookupTable.h"
#include "vtkProperty.h"
#include <vtkImageProperty.h>
// call back classes
#include "vtkPlaneSource.h" 
#include "vtkInteractorStyleImage.h"
#include "vtkImageMapToWindowLevelColors.h"
#include "vtkImagePlaneWidget.h"
// generate contour
#include <vtkPlane.h>
#include <vtkCutter.h>
// display measurement/calibration points
#include <vtkSphereSource.h>
#include <vtkLineSource.h>
// interactor
#include <vtkObjectFactory.h>
#include <vtkRenderViewBase.h>
// chart XY
#include <vtkTable.h>
#include <vtkDoubleArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkChartXY.h>
#include <vtkAxis.h>
#include <vtkPlot.h>
#include <vtkPlotArea.h>
#include <vtkScalarBarActor.h>
#include <vtkPen.h>
#include <vtkBrush.h>
#include <vtkTextProperty.h>
#include <vtkChartLegend.h>

// local headers
#include "visualisationframe.h"
#include "corticalbone.h"

#include "vtkExtractVOI.h"
#include "vtkImageActor.h"
#include "vtkImageMapper3D.h"


// display vtk objects
#include "vtkSphereSource.h"
#include "vtkCellLocator.h"
#include <vtkPointData.h>

// volume render
#include "vtkSmartVolumeMapper.h"
#include "vtkVolumeProperty.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "colourMap.h"
//#include "../../../../vtk/VTK/IO/Export/vtkGL2PSExporter.h"

// save displays
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkDepthSortPolyData.h>
#include <vtkImageSliceMapper.h>
#include <vtkImageStack.h>
#include <vtkSelectionNode.h>
#include <vtkGL2PSExporter.h>
#include <vtkCubeSource.h>

// given a button click get the nearest point
class ThreeDimInteractorStyle : public vtkInteractorStyleTrackballCamera {
public:
    static ThreeDimInteractorStyle* New();
    //vtkTypeRevisionMacro(ThreeDimInteractorStyle,vtkInteractorStyleTrackballCamera);

    // todo - include calibration spheres in display when 'set calibration points' is checked.

    virtual void OnRightButtonDown() {
        if(stateMakePointMeasurements && stateDrawMesh) {

            vtkCellPicker* cellPicker = (vtkCellPicker*)this->Interactor->GetPicker();

            // find point
            double meshOpacity = meshActor->GetProperty()->GetOpacity();
            meshActor->GetProperty()->SetOpacity(1); threeDimRenderer->GetRenderWindow()->Render(); // rerender opaque hack to ensure the cell picker gets the nearest mesh point
            int status = cellPicker->Pick(this->Interactor->GetEventPosition()[0], Interactor->GetEventPosition()[1], 0,  // always zero.
                                          threeDimRenderer);
            meshActor->GetProperty()->SetOpacity(meshOpacity); threeDimRenderer->GetRenderWindow()->Render(); // end of rerender hack


            if(status!=0 && cellPicker->GetActor()!=NULL) { // todo find closest along raycast line

                double picked[Dimension];
                cellPicker->GetPickPosition(picked);

                // display point
                vtkIdType pointId = corticalBone->getClosestPoint(picked);
                bool validity = corticalBone->runModellingAtPoint(pointId);

                // update sphere and line
                double centre[Dimension];// , start[Dimension], end[Dimension];
                corticalBone->getMeasurementLocation(centre);
                //corticalBone->getProfileStart(start);
                //corticalBone->getProfileEnd(end);

                sphereDisplay->SetCenter(centre);
                if(sphereActor->GetVisibility()==0) {
                    sphereActor->VisibilityOn();
                }
                sphereActor->VisibilityOff();
                sphereDisplay->Update();

                lineDisplay->SetPoints(corticalBone->getProfilePoints()); // returns vtkPoints todo - uncomment
                //lineDisplay->SetPoint1(start); lineDisplay->SetPoint2(end); // start = inside [Trabecular bone], end = outside [Soft Tissue]
                if(lineActor->GetVisibility()==0) {
                    lineActor->VisibilityOn();
                }
                lineDisplay->Update();

                updateProfileDisplay(pointId, validity);

            } else {
                if(sphereActor->GetVisibility()==1) {
                    sphereActor->VisibilityOff();
                }
                if(lineActor->GetVisibility()==1) {
                    lineActor->VisibilityOff();
                }
            }
            threeDimRenderer->GetRenderWindow()->Render();
        }
    }

    void updateClippingRange() {

        ((vtkRenderWindow*)threeDimView->GetRenderWindow())->Render();
        threeDimRenderer->ResetCameraClippingRange();
    }

    void setMesh() {

        stateDrawMesh=true;

        mesh = corticalBone->getMesh();

        // create and add the mesh actors
        setupMesh();

        // add the point and sphere actors
        threeDimRenderer->AddActor(sphereActor);
        threeDimRenderer->AddActor(lineActor);

        if(meshDisplaSelection ==kVolumeDisplay) {
            meshActor->VisibilityOff();
            meshScaleBarActor->VisibilityOff();
        }

        threeDimRenderer->ResetCamera();
        ((vtkRenderWindow*)threeDimView->GetRenderWindow())->Render();
        threeDimView->update();
    }

    void setMeshContour(vtkSmartPointer<vtkActor> meshContourActorIn) {
        meshContourActor = meshContourActorIn;
        if(displaySlices) {
            threeDimRenderer->AddActor(meshContourActor);
            threeDimRenderer->ResetCamera();
            ((vtkRenderWindow*)threeDimView->GetRenderWindow())->Render();
            threeDimView->update();
        }
    }

    void setImage(vtkSmartPointer<vtkImageStack> imageStackActorIn, vtkSmartPointer<vtkImageData> imageIn) {

        stateDrawImage=true;

        image = imageIn;

        imageStackActor = imageStackActorIn;

            threeDimRenderer->AddActor(imageStackActor);
            threeDimRenderer->ResetCamera();
        ((vtkRenderWindow*)threeDimView->GetRenderWindow())->Render();
            threeDimView->update();

        this->createImageVolume();

        if(meshDisplaSelection !=kVolumeDisplay) {
            imageVolume->VisibilityOff();
        }
        if(displaySlices) {
            threeDimRenderer->SetBackground(1,1,1); // set backgroung white if in 'figure generation mode'
        }
        //threeDimRenderer->SetBackground(1,1,1);


        // set up cell picker
        vtkCellPicker* cellPicker = (vtkCellPicker*)this->Interactor->GetPicker();
        cellPicker->DeletePickList(imageStackActor);


    }

    void updateMeshColours() {

        // clear display
        meshColourTF->RemoveAllPoints();

        int displayIndex = corticalBone->getDisplayIndex();

        if( corticalBone->isMeshMeasured() && displayIndex != itk::CorticalBone::kInvalid) {
            meshMapper->ScalarVisibilityOn();

            // update mesh scalars and get range
            double min, max; corticalBone->getDisplayRange(min, max);
            vtkSmartPointer<vtkDoubleArray> display = corticalBone->getDisplayArray();


            mesh->GetPointData()->SetScalars(display);

            int n = colourMap.getLength();

            double r,g,b;
            for(int i=0; i<n; i++) {
                double mapValue = min + (max-min) * ((double) i) / ((double) (n-1));

                colourMap.getColour(i,r,g,b);
                meshColourTF->AddRGBPoint(mapValue, r/255.0, g/255.0, b/255.0);
                //cout<<"Colours=["<<r<<","<<g<<","<<b<<"]"<<endl;
            }

            colourMap.getNanColour(r,g,b);
            meshColourTF->SetNanColor(r/255.0, g/255.0, b/255.0);
            colourMap.getNanColour(r,g,b);
            meshColourTF->SetBelowRangeColor(r/255.0, g/255.0, b/255.0);
            colourMap.getNanColour(r,g,b);
            meshColourTF->SetAboveRangeColor(r/255.0, g/255.0, b/255.0);

        } else {
            meshMapper->ScalarVisibilityOff();
        }
        meshColourTF->Build();
    }

    void updateImportedProfiles() {
        stateDisplayImportedProfiles = true;
    }

    void removeImportedProfiles() {
        stateDisplayImportedProfiles = false;
    }

    void initialiseThreeDim(itk::CorticalBone::Pointer corticalBoneIn, QVTKWidget* threeDimViewIn, vtkSmartPointer<vtkRenderer> threeDimRendererIn, QVTKWidget* profileViewIn, bool setSliceDisplay = false) {

        corticalBone = corticalBoneIn;

        // set state
        initialiseStates();

        // setup colourmap
        colourMap = ColourMap(); // TODO define with itk::Array
        meshColourTF = vtkSmartPointer<vtkColorTransferFunction>::New();
        meshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();

        threeDimView = threeDimViewIn;
        profileView = profileViewIn;
        threeDimRenderer = threeDimRendererIn;

        displaySlices = setSliceDisplay;

        // initialise actors
        meshActor = NULL;
        imageStackActor = NULL;
        imageVolume = NULL;

        // set up profile display
        profileContextView = vtkSmartPointer<vtkContextView>::New();
        profileChart = vtkSmartPointer<vtkChartXY>::New();
        profileContextView->GetScene()->AddItem(profileChart);

        
        ((vtkRenderWindow*)profileView->GetRenderWindow())->AddRenderer(profileContextView->GetRenderer());
        profileContextView->SetInteractor(profileView->GetInteractor());
        profileView->SetRenderWindow(profileContextView->GetRenderWindow()); // profileView->SetRenderWindow(vtkGenericOpenGLRenderWindow::SafeDownCast(profileContextView->GetRenderWindow())); //updateProfileDisplay();
        profileChart->GetLegend()->SetHorizontalAlignment(vtkChartLegend::LEFT);
        profileChart->GetLegend()->SetVerticalAlignment(vtkChartLegend::TOP);

        // sphere display
        sphereDisplay = vtkSmartPointer<vtkSphereSource>::New();
        vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        sphereMapper->SetInputConnection(sphereDisplay->GetOutputPort());
        sphereActor = vtkSmartPointer<vtkActor>::New();
        sphereActor->SetMapper(sphereMapper);
        sphereActor->VisibilityOff();
        sphereActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
        sphereActor->VisibilityOff();
        threeDimRenderer->AddActor(sphereActor);

        // line display
        lineDisplay = vtkSmartPointer<vtkLineSource>::New();
        vtkSmartPointer<vtkPolyDataMapper> lineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        lineMapper->SetInputConnection(lineDisplay->GetOutputPort());
        lineActor = vtkSmartPointer<vtkActor>::New();
        lineActor->SetMapper(lineMapper);
        lineActor->VisibilityOff();
        lineActor->GetProperty()->SetLineWidth(2);
        lineActor->GetProperty()->SetColor(0.0, 1.0, 1.0);
        lineActor->VisibilityOff();

    } // todo - override initilize and call from within

    void reset() {

        initialiseStates();

        threeDimRenderer->RemoveActor(imageStackActor);
        imageStackActor = NULL;

        threeDimRenderer->RemoveActor(meshActor);
        meshActor = NULL;

        threeDimRenderer->RemoveActor2D( meshScaleBarActor );
        meshScaleBarActor = NULL;

        meshColourTF = vtkSmartPointer<vtkColorTransferFunction>::New();
        meshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();

        threeDimRenderer->RemoveAllViewProps(); // removes all actors

        disablePointMeasurements();

        mesh = NULL;
        image = NULL;

        // reset view
        ((vtkRenderWindow*)threeDimView->GetRenderWindow())->Render();
        threeDimView->update(); threeDimRenderer->ResetCamera();


    }

    void initialiseStates() {
        // set state
        stateDrawMesh=false;
        stateDrawImage=false;
        stateMakePointMeasurements = false; stateDisplayImportedProfiles = false;
        meshDisplaSelection = kMeshDisplay;
    }

    void enablePointMeasurements() {
        stateMakePointMeasurements=true;
    }

    void disablePointMeasurements() {
        stateMakePointMeasurements=false;
        if(sphereActor->GetVisibility()!=0) {
            sphereActor->VisibilityOff();
            threeDimRenderer->GetRenderWindow()->Render();
        }
        if(lineActor->GetVisibility()!=0) {
            lineActor->VisibilityOff();
            threeDimRenderer->GetRenderWindow()->Render();
        }
        if(profileChart->GetNumberOfPlots()!=0){
            profileChart->ClearPlots();
        }
    }

    void setDisplayParameter(int index) {
        corticalBone->setDisplayIndex(index);

        updateMeshColours();

        meshScaleBarActor->SetTitle(corticalBone->getDisplayName().c_str());
    }

    void setDisplayMesh(int index) { // todo - simplify

        meshDisplaSelection = index;

        // update mesh actor and color maps
        updateMesh();
        updateMeshColours();

        // volume display
        if(stateDrawImage && meshDisplaSelection !=kVolumeDisplay){ // no volume display
            imageVolume->VisibilityOff();
        } else if(stateDrawImage && meshDisplaSelection ==kVolumeDisplay){ // volume display
            imageVolume->VisibilityOff();
        }

        // mesh display
        if(stateDrawMesh && meshDisplaSelection !=kVolumeDisplay) { // mesh display
            meshActor->VisibilityOn();
            meshScaleBarActor->VisibilityOn();
        } else if(stateDrawImage && meshDisplaSelection ==kVolumeDisplay){ // no mesh display
            imageVolume->VisibilityOn();
        }
    }

    void saveDisplays(std::string profileName, std::string threeDimName) {

        ((vtkRenderWindow*)profileView->GetRenderWindow())->Render();
        ((vtkRenderWindow*)threeDimView->GetRenderWindow())->Render();

        // save chart - TODO only if displaying anything
        if(profileName.find(".png") != std::string::npos) {

            vtkSmartPointer<vtkWindowToImageFilter> chartImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
            chartImageFilter->SetInput((vtkRenderWindow*)profileView->GetRenderWindow());
            //chartImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window) // todo uncomment when bug fixed - https://github.com/UV-CDAT/uvcdat/issues/1148
            chartImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
            chartImageFilter->ReadFrontBufferOff(); // read from the back buffer
            chartImageFilter->Update();

            vtkSmartPointer<vtkPNGWriter> chartWriter = vtkSmartPointer<vtkPNGWriter>::New();
            chartWriter->SetFileName(profileName.c_str());
            chartWriter->SetInputConnection(chartImageFilter->GetOutputPort());
            chartWriter->Write();

            // save 3D display with multiple slices
            //return saveTransparentDisplays(threeDimName);

            // save 3D display
            vtkSmartPointer<vtkWindowToImageFilter> threeDimImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
            threeDimImageFilter->SetInput((vtkRenderWindow*)threeDimView->GetRenderWindow());
            threeDimImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
            threeDimImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
            threeDimImageFilter->ReadFrontBufferOff(); // read from the back buffer
            threeDimImageFilter->Update();

            vtkSmartPointer<vtkPNGWriter> threeDimWriter = vtkSmartPointer<vtkPNGWriter>::New();
            threeDimWriter->SetFileName(threeDimName.c_str());
            threeDimWriter->SetInputConnection(threeDimImageFilter->GetOutputPort());
            threeDimWriter->Write();
        } else if(profileName.find(".eps") != std::string::npos) {

            std::string profileNameStub = profileName.substr(0, profileName.find_last_of(".")); // todo - fix profile y label 3 superscript

            vtkSmartPointer<vtkGL2PSExporter> profileWriter = vtkSmartPointer<vtkGL2PSExporter>::New();
            profileWriter->SetRenderWindow((vtkRenderWindow*)profileView->GetRenderWindow());
            profileWriter->SetSortToBSP();
            profileWriter->SetFileFormatToPDF();
            profileWriter->UsePainterSettings();
            profileWriter->CompressOff();
            profileWriter->DrawBackgroundOff();
            profileWriter->TextAsPathOn();
            profileWriter->Write3DPropsAsRasterImageOn();
            profileWriter->SetFilePrefix(profileNameStub.c_str());

            profileWriter->Update();
            profileWriter->Write();

            std::string threeDimNameStub = threeDimName.substr(0, threeDimName.find_last_of(".")); // todo - fix missing line & slice display

            vtkSmartPointer<vtkGL2PSExporter> threeDimWriter = vtkSmartPointer<vtkGL2PSExporter>::New();
            threeDimWriter->SetRenderWindow((vtkRenderWindow*)threeDimView->GetRenderWindow());
            threeDimWriter->SetSortToBSP();
            threeDimWriter->SetFileFormatToPDF();
            threeDimWriter->UsePainterSettings();
            threeDimWriter->CompressOff();
            threeDimWriter->DrawBackgroundOff();
            threeDimWriter->TextAsPathOn();
            threeDimWriter->Write3DPropsAsRasterImageOn();
            threeDimWriter->SetFilePrefix(threeDimNameStub.c_str());

            threeDimWriter->Update();
            threeDimWriter->Write();
        }

    }

    void setSphereRadius(double radius) {
        // update sphere size
        sphereDisplay->SetRadius(radius);
    }

private:

    typedef enum { // model type
        kMeshDisplay=0,
        kPeriostealDisplay,
        kVolumeDisplay,
    } displayOptions;


    void updateProfileDisplay(vtkIdType pointId, bool validity) {

        profileChart->ClearPlots(); // creases previous plots

        bool classifierMode = corticalBone->isClassifierMode();

        if(!classifierMode) {
            profileChart->SetShowLegend(true);
        } else {
            profileChart->SetShowLegend(false);
        }
        profileChart->SetTitle("Cortex Densities in Profile");
        profileChart->GetTitleProperties()->SetFontSize(24);
        profileChart->GetTitleProperties()->SetFontFamilyToArial();

        vtkIdType numberOfRows;
        vtkSmartPointer<vtkTable> table; vtkSmartPointer<vtkPlot> plot;

        // display sampled image profile(s)
        table = corticalBone->getImageProfileTable(pointId);
        numberOfRows = table->GetNumberOfColumns();
        if(numberOfRows>=2) { // the average profile if multiple profiles include
            // individual profiles
            for(int i=2; i<numberOfRows; i++) {
                plot = profileChart->AddPlot(vtkChart::LINE);
                plot->SetInputData(table, 0, i);
                if(i*2==numberOfRows) { // attempt to pull out the central profile
                    plot->SetColor(0.75, 0.0, 0.0);
                    plot->SetWidth(1.0);
                } else { // single sampled profile
                    plot->SetColor(.75, 0.0, 0.0);
                    plot->SetWidth(1.0);
                }
            }
            // averaged overall profile
            plot = profileChart->AddPlot(vtkChart::LINE);
            plot->SetInputData(table, 0, 1);
            if(validity) { plot->SetColor(1.0, 0.0, 0.0);
            } else {plot->SetColor(0.75, 0.0, 0.0);
            }
            plot->SetWidth(2.5);
        }

        // display imported image profile (if relevant)
        table = corticalBone->getImportedProfileTable(pointId);
        numberOfRows = table->GetNumberOfColumns();
        if(numberOfRows==2) {
            plot = profileChart->AddPlot(vtkChart::LINE); // imported image
            plot->SetInputData(table, 0, 1);
            plot->SetColor(1.0, 0.0, 0.0); plot->SetWidth(3.0);
            plot->GetPen()->SetLineType(vtkPen::DASH_LINE);
        } else if(numberOfRows==4) {
            plot = profileChart->AddPlot(vtkChart::LINE); // imported image
            plot->SetInputData(table, 0, 1);
            plot->SetColor(1.0, 0.0, 0.0); plot->SetWidth(2.0);
            plot->GetPen()->SetLineType(vtkPen::DASH_LINE);
            plot = profileChart->AddPlot(vtkChart::LINE); // imported classifications
            plot->SetInputData(table, 0, 2);
            plot->SetColor(0.0, 1.0, 1.0); plot->SetWidth(2.0);
            plot->GetPen()->SetLineType(vtkPen::DASH_LINE);
            plot = profileChart->AddPlot(vtkChart::LINE); // imported percentages
            plot->SetInputData(table, 0, 3);
            plot->SetColor(0.0, 1.0, 0.25); plot->SetWidth(2.0);
            plot->GetPen()->SetLineType(vtkPen::DASH_LINE);
        }

        // display the display model - model (unblurred / best fit lines)
        table = corticalBone->getDisplayModelTable(pointId);
        numberOfRows = table->GetNumberOfColumns();
        for(int i=0; i<numberOfRows; i+=2) {
            plot = profileChart->AddPlot(vtkChart::LINE);
            plot->SetInputData(table, i, i+1);
            if(validity && classifierMode) { plot->SetColor(0.0, 0.0, 1.0);
            } else if(!validity && classifierMode) {plot->SetColor(0.0, 0.0, 0.75);
            } else if(validity && i==0) { plot->SetColor(0.0, 0.0, 1.0);
            } else { plot->SetColor(0.0, 1.0, 1.0);   }
            plot->SetWidth(2.5);
        }

//        // display the weights
//        table = corticalBone->getWeightsTable(pointId);
//        numberOfRows = table->GetNumberOfColumns();
//        if(numberOfRows==2) {
//            plot = profileChart->AddPlot(vtkChart::LINE);
//            plot->SetInputData(table, 0, 1);
//            plot->SetColor(1.0, 0.0, 1.0);
//            plot->SetWidth(2.5); plot->GetPen()->SetLineType(vtkPen::DASH_LINE);
//        }

        // display the imported display model (if relevant)
        table = corticalBone->getImportedDisplayModelTable(pointId);
        numberOfRows = table->GetNumberOfColumns();
        for(int i=0; i<numberOfRows; i+=2) {
            plot = profileChart->AddPlot(vtkChart::LINE);
            plot->SetInputData(table, i, i+1);
            plot->SetColor(0.0, 0.0, 0.75);
            plot->SetWidth(2.5);
            plot->GetPen()->SetLineType(vtkPen::DASH_LINE);
        }

        // display the processing model - classifications and percentages / blurred model
        table = corticalBone->getProcessingModelTable(pointId);
        numberOfRows = table->GetNumberOfColumns();
        if(numberOfRows==2) {
            plot = profileChart->AddPlot(vtkChart::LINE);
            plot->SetInputData(table, 0, 1);
            if(validity) { plot->SetColor(0.0, 1.0, 1.0);
            } else {plot->SetColor(0.0, 0.75, 0.75);
            }
            plot->SetWidth(2.5);
        } else if(numberOfRows==3) {
            plot = profileChart->AddPlot(vtkChart::LINE);
            plot->SetInputData(table, 0, 2); // plot percentages below if present
            if(validity) { plot->SetColor(0.0, 0.75, 0.25);
            } else {plot->SetColor(0.0, 0.5, 0.15);
            }  plot->SetWidth(2.5);
            plot = profileChart->AddPlot(vtkChart::LINE);
            plot->SetInputData(table, 0, 1);
            if(validity) { plot->SetColor(0.0, 1.0, 1.0);
            } else {plot->SetColor(0.0, 0.75, 0.75);
            }  plot->SetWidth(2.5);
        }

        // display the error area
        vtkSmartPointer<vtkTable> areaTable = corticalBone->getErrorAreaTable(pointId);
        numberOfRows = areaTable->GetNumberOfColumns();
        if(numberOfRows==3) {
            vtkSmartPointer<vtkPlotArea> areaPlot = vtkPlotArea::SafeDownCast(profileChart->AddPlot(vtkChart::AREA));
            areaPlot->SetInputData(areaTable);
            areaPlot->SetInputArray(0, "Position");
            areaPlot->SetInputArray(1, "Min Y");
            areaPlot->SetInputArray(2, "Max Y");

            areaPlot->SetColor(1.0, 0.0, 0.0); areaPlot->SetWidth(2.5);
            areaPlot->GetBrush()->SetColor(1.0, 0.0, 0.0);
            areaPlot->GetBrush()->SetOpacityF(.5);

        }

        // set up labels
        profileChart->GetLegend()->GetLabelProperties()->SetFontSize(16);
        profileChart->GetLegend()->GetLabelProperties()->SetFontFamilyToArial();

        vtkSmartPointer<vtkAxis> yAxis = profileChart->GetAxis(vtkAxis::LEFT);
        vtkSmartPointer<vtkAxis> xAxis = profileChart->GetAxis(vtkAxis::BOTTOM);
        xAxis->SetTitle("Position (mm)");
        if(corticalBone->isCalibrated() == true) {
            yAxis->SetTitle("Bone Density (mg/cm\u00B3)");
        } else {
            yAxis->SetTitle("Bone Density (HU)");
        }
        xAxis->GetTitleProperties()->SetFontSize(20);
        xAxis->GetTitleProperties()->BoldOff();
        xAxis->GetTitleProperties()->SetFontFamilyToArial();
        yAxis->GetTitleProperties()->SetFontSize(20);
        yAxis->GetTitleProperties()->SetFontFamilyToArial();
        yAxis->GetTitleProperties()->BoldOff();

        xAxis->GetLabelProperties()->SetFontSize(16);
        xAxis->GetLabelProperties()->SetFontFamilyToArial();
        yAxis->GetLabelProperties()->SetFontSize(16);
        yAxis->GetLabelProperties()->SetFontFamilyToArial();

        /*xAxis->SetBehavior(vtkAxis::FIXED);
        yAxis->SetBehavior(vtkAxis::FIXED);

        double offset=corticalBone->getPeriostealOffset();
        xAxis->SetRange(-6+offset,14+offset);
        yAxis->SetRange(-600,1400);
        yAxis->SetNumberOfTicks(11);*/

//      profileContextView->GetRenderWindow()->SetMultiSamples(0); // looks bad
//      profileContextView->GetInteractor()->Initialize();

    }

    void updateMesh() {

        // set to initial mesh
        if(meshDisplaSelection==kMeshDisplay) { // set to periosteal mesh - must be measured
            mesh = corticalBone->getMesh();
        } else if(meshDisplaSelection==kPeriostealDisplay && corticalBone->isMeshMeasured()) { // if volume - leave unchanged just turn off visibility
            mesh = corticalBone->getPeriostealMesh();
        } else if(meshDisplaSelection==kVolumeDisplay) { // if unmeasure and periosteal mesh type - leave unchanged and print out an error
            // do nothing
            return;
        } else if(meshDisplaSelection==kPeriostealDisplay && !corticalBone->isMeshMeasured()) {
            cerr<<"Error in ThreeDimInteractorStyle::updateMesh() periosteal selected before it has been created."<<endl;
            return;
        } else { // unexpected combination. give a warning
            cerr<<"Error in ThreeDimInteractorStyle::updateMesh() unexpected combination."<<endl;
            return;
        }

        // remove previous mesh actor
        threeDimRenderer->RemoveActor(meshActor);
        threeDimRenderer->RemoveActor2D(meshScaleBarActor);

        // setup the mesh actor
        setupMesh();


    }

    void setupMesh() { // requires mesh to already be set

        // depth sorting the mesh
        vtkSmartPointer<vtkDepthSortPolyData> meshDepthSorter = vtkSmartPointer<vtkDepthSortPolyData>::New();
        meshDepthSorter->SetInputData(mesh); //->SetInputConnection(mesh->GetOutputPort());
        meshDepthSorter->SetDirectionToBackToFront();
        meshDepthSorter->SetVector(1, 1, 1);
        meshDepthSorter->SetCamera(threeDimRenderer->GetActiveCamera());
        meshDepthSorter->SortScalarsOff(); // do not really need this here
        //meshDepthSorter->SetDepthSortModeToBoundsCenter(); //meshDepthSorter->SetDepthSortModeToFirstPoint();
        meshDepthSorter->SetDepthSortModeToParametricCenter();
        meshDepthSorter->Update();

        // mesh depth mapper
        meshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        meshMapper->SetInputConnection(meshDepthSorter->GetOutputPort());
        meshMapper->SetScalarModeToUsePointData();
        meshMapper->SetColorModeToMapScalars();

        // mesh actor
        meshActor = vtkSmartPointer<vtkActor>::New();
        meshActor->SetMapper(meshMapper); //meshActor->GetMapper()->SetInputConnection(meshDepthSorter->GetOutputPort()); // meshActor->SetMapper(meshMapper);
        meshActor->SetVisibility(true); //meshActor->GetProperty()->SetRepresentationToSurface(); meshActor->GetProperty()->LightingOn();
        meshActor->GetProperty()->SetOpacity(0.9);

        // colour bar transform
        meshColourTF = vtkSmartPointer<vtkColorTransferFunction>::New();
        meshColourTF->SetColorSpaceToRGB ();
        meshMapper->SetLookupTable(meshColourTF);

        // scale bar actor
        meshScaleBarActor = vtkSmartPointer<vtkScalarBarActor>::New(); //meshScaleBar->SetLookupTable(meshMapper->GetLookupTable());
        meshScaleBarActor->SetTitle("Blank");
        meshScaleBarActor->SetNumberOfLabels(5);
        meshScaleBarActor->SetLookupTable(meshColourTF);//corticalBone->getMeshMapper()->GetLookupTable());

        // add the actors to the scene
        threeDimRenderer->AddActor ( meshActor );
        threeDimRenderer->AddActor2D( meshScaleBarActor );

        threeDimRenderer->ResetCamera();
        ((vtkRenderWindow*)threeDimView->GetRenderWindow())->Render();
        threeDimView->update();

        // set up cell picker - reset pick selection
        vtkCellPicker* cellPicker = (vtkCellPicker*)this->Interactor->GetPicker();
        cellPicker->PickFromListOn();
        vtkSmartPointer<vtkPropCollection> actors = cellPicker->GetPickList();
        actors->InitTraversal();
        for(int i=0;i<actors->GetNumberOfItems();i++) {
            cellPicker->DeletePickList(actors->GetNextProp());
        }
        cellPicker->AddPickList(meshActor);
    }

    void createImageVolume(){

        // create volume mapper
        vtkSmartPointer<vtkSmartVolumeMapper> volumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
        volumeMapper->SetBlendModeToComposite(); // composite first
        volumeMapper->SetInputData(image);

        // set properties of volume rendering
        vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
        volumeProperty->SetInterpolationTypeToLinear();
        //volumeProperty->SetInterpolationType(VTK_LINEAR_INTERPOLATION);
        volumeMapper->SetBlendModeToComposite();
        volumeProperty->ShadeOn(); //volumeProperty->ShadeOff();
        volumeProperty->SetAmbient(0.2);
        volumeProperty->SetDiffuse(0.9);
        volumeProperty->SetSpecular(0.2);
        volumeProperty->SetSpecularPower(10.0);
        volumeProperty->SetScalarOpacityUnitDistance(0.8919);

        // Create function to define opacity for different tissue types.
        vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = vtkSmartPointer<vtkPiecewiseFunction>::New();
        compositeOpacity->ClampingOn();
        // Create function to map voxel intensities to colours.
        vtkSmartPointer<vtkColorTransferFunction> volumeColour = vtkSmartPointer<vtkColorTransferFunction>::New();
        volumeColour->ClampingOn();
        // midpoint = medium, sharpness- 0=linear, 1=constant

        double air[2] = {std::numeric_limits<double>::min(), 499};
        double soft[2] = {-600, -500}; // The goal is to one colour for flesh (between 500 and 1000)
        double muscle[2] = {50, 60};
        double bone[2]  = {500, std::numeric_limits<double>::max()}; // and another colour for bone (1150 and over).

        const double midpoint=0.5;
        const double sharpness[2]={0,1}; // 0=linear, 1=constant

        // air - no opacity
        /*volumeColour->AddRGBPoint( air[0], 0, 0, 0, midpoint, sharpness[1]);
        volumeColour->AddRGBPoint( air[1], 0, 0, 0, midpoint, sharpness[1]);
        compositeOpacity->AddPoint(air[0], 0, midpoint, sharpness[1] );
        compositeOpacity->AddPoint(air[1], 0, midpoint, sharpness[1] );

        // bone - opaque
        volumeColour->AddRGBPoint( bone[0], 1, 1, 1, midpoint, sharpness[1]);
        volumeColour->AddRGBPoint( bone[1], 1, 1, 1, midpoint, sharpness[1]);
        compositeOpacity->AddPoint(bone[0], 0, midpoint, sharpness[1] );
        compositeOpacity->AddPoint(bone[1], 0, midpoint, sharpness[1] );*/

        // old values
        volumeColour->AddRGBPoint( 0, 0.2, 0.1, 0.1, 0.5, 0.0);
        compositeOpacity->AddPoint(0, 0, 0.5, 0.0 );
        volumeColour->AddRGBPoint( 60, 0.73, 0.6, 0.40, 0.49, 0.61 );
        volumeColour->AddRGBPoint( 100, .72, .36, .68, 0.5, 0.0 );
        compositeOpacity->AddPoint(60, 0, 0.8, .61 );
        compositeOpacity->AddPoint(100, 0, 0.8, .61 );

        volumeColour->AddRGBPoint( 250, .62, 0.36, .28, 0.5, 0.0 );
        volumeColour->AddRGBPoint( 400, 0.85, 0.82, 0.76, .5, 0.0 );
        compositeOpacity->AddPoint(250, 0, 0.5, 0.0 );
        compositeOpacity->AddPoint(400, 0.9, 0.5, 0.0 );

        volumeColour->AddRGBPoint( 800, 0.95, 0.9, 0.8, 0.5, 0.0 );
        compositeOpacity->AddPoint(800, 1, 0.5, 0.0);


        // set colour and opacity
        volumeProperty->SetScalarOpacity(compositeOpacity); // must add first to get colour
        volumeProperty->SetColor(0, compositeOpacity);

        imageVolume = vtkSmartPointer<vtkVolume>::New();
        imageVolume->SetMapper(volumeMapper);
        imageVolume->SetProperty(volumeProperty);
        threeDimRenderer->AddViewProp(imageVolume);

        // 3D texture mode for coverage.
        //volumeMapper->SetRequestedRenderModeToRayCastAndTexture(); // not supported?

        // Software mode, for coverage. Ensures the same regression image on all platforms.
        volumeMapper->SetRequestedRenderModeToRayCast();

    }

    // mesh 3D values - initial mesh
    vtkSmartPointer<vtkPolyData> mesh;
    vtkSmartPointer<vtkActor> meshActor;
    vtkSmartPointer<vtkActor> meshContourActor;

    // periosteal mesh 3D values
    vtkSmartPointer<vtkPolyData> periostealMesh;
    vtkSmartPointer<vtkActor> periostealMeshActor;

    // general mesh values
    vtkSmartPointer<vtkPolyDataMapper> meshMapper;
    vtkSmartPointer<vtkColorTransferFunction> meshColourTF;

    // mesh scale bar
    vtkSmartPointer<vtkScalarBarActor> meshScaleBarActor;

    // image 3D values
    vtkSmartPointer<vtkImageData> image;
    //vtkSmartPointer<vtkImageActor> imageActor;
    vtkSmartPointer<vtkVolume> imageVolume;

    // sphere display values
    vtkSmartPointer<vtkSphereSource> sphereDisplay;
    vtkSmartPointer<vtkActor> sphereActor;

    // line display
    vtkSmartPointer<vtkLineSource> lineDisplay;
    vtkSmartPointer<vtkActor> lineActor;

    // profile displays
    vtkSmartPointer<vtkContextView> profileContextView;
    vtkSmartPointer<vtkChartXY> profileChart;

    // slice displays - for making figures
    vtkSmartPointer<vtkImageStack> imageStackActor;
    vtkSmartPointer<vtkActor> meshSliceActor;

    // state values
    bool stateDrawMesh, stateDrawImage, displaySlices;
    bool stateMakePointMeasurements; bool stateDisplayImportedProfiles;
    int meshDisplaSelection;

    // general values
    itk::CorticalBone::Pointer corticalBone;

    QVTKWidget* threeDimView;
    QVTKWidget* profileView;

    vtkSmartPointer<vtkRenderer> threeDimRenderer;

    ColourMap colourMap;
};
//vtkCxxRevisionMacro(ThreeDimInteractorStyle, "$Revision: 1.1 $");
vtkStandardNewMacro(ThreeDimInteractorStyle);

class SliceInteractorStyle : public vtkInteractorStyleImage {
public:
    static SliceInteractorStyle* New();
    //vtkTypeRevisionMacro(SliceInteractorStyle,vtkInteractorStyleImage);

    virtual void OnMouseWheelForward() { // decrement

        if(drawImage&&!displayFixedSlices) {
            int displaySlice = imageSliceMapper->GetSliceNumber();
            int imageExtents[6]; image->GetExtent(imageExtents);

            if(displaySlice>imageExtents[4]) {
                imageSliceMapper->SetSliceNumber(displaySlice-1);

                if (drawMesh) {
                    updateMeshCut();
                }
                if(calibrationPhantomType !=itk::CorticalBone::kNoCal && calibrationOn) {
                    // update calibration displays
                }

                threeDimInteractor->updateClippingRange();
                ((vtkRenderWindow*)sliceView->GetRenderWindow())->Render();
                sliceRenderer->ResetCameraClippingRange();
            }
        }
    }

    virtual void OnMouseWheelBackward() { // increment

        if(drawImage&&!displayFixedSlices) {
            int imageExtents[6]; image->GetExtent(imageExtents);
            int displaySlice = imageSliceMapper->GetSliceNumber();


            if(displaySlice<imageExtents[5]) {

                imageSliceMapper->SetSliceNumber(displaySlice+1);
                if (drawMesh) {
                    updateMeshCut();
                }

                threeDimInteractor->updateClippingRange();
                ((vtkRenderWindow*)sliceView->GetRenderWindow())->Render();
                sliceRenderer->ResetCameraClippingRange();

            }
        }
    }

    virtual void OnLeftButtonDown() {

        if(calibrationPhantomType <=itk::CorticalBone::kManualControlPtsCal && calibrationPhantomType > itk::CorticalBone::kNoCal && calibrationOn && drawImage) {
            // get point location
            vtkCellPicker* cellPicker = (vtkCellPicker*)this->Interactor->GetPicker();
            int status = cellPicker->Pick(this->Interactor->GetEventPosition()[0], Interactor->GetEventPosition()[1], 0,  // always zero.
                                          this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
            double picked[Dimension];
            cellPicker->GetPickPosition(picked);
            std::cout<<"Selected Cal Point["<<calibrationPointIndex<<"], x: "<<picked[0]<<", y: "<<picked[1]<<", z: "<<picked[2]<<std::endl;

            // create sphere in right place
            cubes[calibrationPointIndex]->SetCenter(picked); // spheres
            cubes[calibrationPointIndex]->Update(); // spheres
            sliceRenderer->GetRenderWindow()->Render();

            calibrationPointIndex++;
            if(calibrationPointIndex >=cubes.size()) { // spheres.size()
                calibrationPointIndex=0; calibrationComplete=true;
            }

        }

        vtkInteractorStyleImage::OnLeftButtonDown();
    }

    virtual void OnMouseMove() { // override to remove delay in application of contrast adjustments to the threeDimView render window

        vtkInteractorStyleImage::OnMouseMove();
        threeDimInteractor->updateClippingRange();
    }

    void setMesh() {

        drawMesh=true;

        // mesh
        mesh = corticalBone->getMesh();

        // create the mesh cutting / display objects
        setupMesh();

        // update renderer
        sliceRenderer->ResetCamera();

        // add contours to three dim view
        threeDimInteractor->setMeshContour(meshActor);

    }

    void setImage(vtkSmartPointer<vtkImageData> imageIn) {

        drawImage=true;
        image = imageIn;

        if(!displayFixedSlices) {

            // Set Up the image actor (Stack)
            imageStackActor = vtkSmartPointer<vtkImageStack>::New();
            imageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
            imageSliceMapper->SetInputData(image);
            imageSliceMapper->SetSliceNumber(imageSliceMapper->GetSliceNumberMinValue()); // could crop by mesh if desired
            vtkSmartPointer<vtkImageSlice> imageSlice = vtkSmartPointer<vtkImageSlice>::New();
            imageSlice->SetMapper(imageSliceMapper);
            imageStackActor->AddImage(imageSlice);
            //imageSlice->GetProperty()->SetLookupTable(imageWindowing); imageSlice->GetProperty()->SetUseLookupTableScalarRange(1); // set to use image windowing

            // set up image windowing
            //double imageRange[2]; image->GetScalarRange(imageRange); // TODO - investiage automatic windowing detection
            //imageSlice->GetProperty()->SetColorLevel((imageRange[1]+imageRange[0])/2);
            //imageSlice->GetProperty()->SetColorWindow((imageRange[1]-imageRange[0])/10);

            imageStackActor->Update();
            sliceRenderer->AddActor(imageStackActor);
            sliceRenderer->ResetCamera();

            if (drawMesh==true) {
                updateMeshCut();
            }
        } else {
            setFixedImageSlices();
        }
        threeDimInteractor->setImage(imageStackActor, image);

        // set calibration spheres - if calibration is one
        if(calibrationPhantomType !=itk::CorticalBone::kNoCal) {
            initiliseSphereLocations();
        }
        if(calibrationOn && calibrationPhantomType <=itk::CorticalBone::kManualControlPtsCal && calibrationPhantomType > itk::CorticalBone::kNoCal) { // should never be turned on before image set
            showSpheres();
        }

    }

    void setDisplayMesh(int index) { // todo - simplify

        meshDisplaSelection = index;

        // update mesh actor and color maps
        updateMesh();


    }

    void initialiseSlice(vtkSmartPointer<ThreeDimInteractorStyle> threeDimInteractorIn, itk::CorticalBone::Pointer corticalBoneIn, QVTKWidget* sliceViewIn, vtkSmartPointer<vtkRenderer> sliceRendererIn, bool setSliceDisplay=false) {

        threeDimInteractor = threeDimInteractorIn;

        corticalBone = corticalBoneIn;
        sliceView = sliceViewIn;
        sliceRenderer = sliceRendererIn;

        displayFixedSlices = setSliceDisplay;

        // set up calibration points
        cubes.clear(); actors.clear(); // spheres.clear();

        // define states
        initialiseStates();

    } // todo - override initilize and call from within

    void startCalibrating() {
        
        calibrationOn=true;

        if(drawImage && calibrationPhantomType <=itk::CorticalBone::kManualControlPtsCal && calibrationPhantomType > itk::CorticalBone::kNoCal && calibrationOn) {
            showSpheres();
            sliceRenderer->GetRenderWindow()->Render(); // todo - remove redundent calls to render window from within nested methods like 'showSpheres()'
        }
    }

    void stopCalibrating() {

        calibrationOn=false;

        // hide any existing spheres
        hideSpheres();
    }

    void setCalibrationMode(int index=0) {

        // set new calibration mode index
        calibrationPhantomType = index; calibrationComplete=false; calibrationPointIndex=0;

        // set up the spheres - only show if calibration mode
        for(int i=0; i<actors.size(); i++) {
            sliceRenderer->RemoveActor(actors[i]);
        }

        // create spheres
        int numberOfCalibrationPoints; double radius;
        corticalBone->getCalibrationPtGeometry(radius, numberOfCalibrationPoints);
        cubes.resize((unsigned long)numberOfCalibrationPoints); actors.resize((unsigned long)numberOfCalibrationPoints);
        // spheres.resize((unsigned long)numberOfCalibrationPoints);
        for(int i=0; i<numberOfCalibrationPoints; i++) {
            cubes[i]=vtkSmartPointer<vtkCubeSource>::New(); // spheres[i]=vtkSmartPointer<vtkSphereSource>::New();

            vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(cubes[i]->GetOutputPort()); // spheres

            actors[i] = vtkSmartPointer<vtkActor>::New();
            actors[i]->SetMapper(mapper);
            actors[i]->VisibilityOff();

            sliceRenderer->AddActor(actors[i]);

            actors[i]->GetProperty()->SetColor(1.0, ((double)i)/(numberOfCalibrationPoints-1), ((double)i)/(numberOfCalibrationPoints-1)); // range from red to white
        }

        if(drawImage && calibrationPhantomType !=itk::CorticalBone::kNoCal) {
            initiliseSphereLocations();
        }

        if(calibrationOn) { // update the calibration values if already set
            startCalibrating();
        }
    }

    vtkSmartPointer<vtkDoubleArray> getCalibrationPts() { // todo get pts as array

        vtkSmartPointer<vtkDoubleArray> points = vtkSmartPointer<vtkDoubleArray>::New();
        unsigned int numberOfSpheres=(unsigned int)cubes.size(); //spheres.size();

        if(calibrationComplete && calibrationPhantomType <=itk::CorticalBone::kManualControlPtsCal && calibrationPhantomType > itk::CorticalBone::kNoCal && numberOfSpheres>0) {
            points->SetNumberOfComponents(Dimension); points->SetNumberOfTuples(numberOfSpheres);
            for(int i=0; i<numberOfSpheres; i++) {
                points->SetTuple(i, cubes[i]->GetCenter()); // spheres
            }
        }
        return points;
    }

    vtkSmartPointer<vtkDoubleArray> getCalibrationValues() { // todo get pts as array
        return corticalBone->getControlPointValues();
    }

    bool setCalibrationPts(int &phantomType, vtkSmartPointer<vtkDoubleArray> points) { // todo set mode and points

        int numberOfCalibrationPoints; double radius;
        corticalBone->getCalibrationPtGeometry(radius, numberOfCalibrationPoints);

        if(numberOfCalibrationPoints!=points->GetNumberOfTuples()) {
            setCalibrationMode(itk::CorticalBone::kNoCal); // todo - ensure the control frame selection is updated accordingly.
            cout<<"Error: In VisualisationFame::setCalibrationPts() numberOfCalibrationPoints does not match the GetNumberOfTuples()"<<endl;
            calibrationComplete=false;
            return false;
        } else if (phantomType<=itk::CorticalBone::kNoCal || phantomType>itk::CorticalBone::kManualControlPtsCal) {
            setCalibrationMode(itk::CorticalBone::kNoCal); // todo - ensure the control frame selection is updated accordingly.
            calibrationComplete=false;
        } else {
            setCalibrationMode(phantomType);
            //startCalibrating();
            for(int i=0; i<points->GetNumberOfTuples(); i++) {
                cubes[i]->SetCenter(points->GetTuple(i)); // spheres
            }
            calibrationComplete=true;
        }
        return calibrationComplete; // return runCalibration();
    }

    bool runCalibration() {
        bool status = false;
        if(calibrationComplete) {
            vtkSmartPointer<vtkDoubleArray> calPts = getCalibrationPts();
            status = corticalBone->runCalibration(calPts);
        } else {
            std::cout<<"Please set all "<<cubes.size()<<" calibration points before calibration can be performed"<<std::endl; // spheres.size()
        }

        return status;
    }

    bool reset() {
        initialiseStates();

        stopCalibrating();

        sliceRenderer->RemoveActor(meshActor);
        sliceRenderer->RemoveActor(imageStackActor);

        sliceRenderer->ResetCamera();
        ((vtkRenderWindow*)sliceView)->Render();
        sliceView->update();
        
        return true;
    }

    void saveDisplays(std::string fileName) {

        sliceRenderer->GetRenderWindow()->Render();

        // save slice view
        if(fileName.find(".png") != std::string::npos) {

            vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
            windowToImageFilter->SetInput(sliceRenderer->GetRenderWindow());
            windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
            windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
            windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
            windowToImageFilter->Update();

            vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
            writer->SetFileName(fileName.c_str());
            writer->SetInputConnection(windowToImageFilter->GetOutputPort());
            writer->Write();
        } else if(fileName.find(".eps") != std::string::npos) {

            std::string rawname = fileName.substr(0, fileName.find_last_of(".")); // todo - fix missing slice display

            vtkSmartPointer<vtkGL2PSExporter> writer = vtkSmartPointer<vtkGL2PSExporter>::New();
            writer->SetRenderWindow(sliceRenderer->GetRenderWindow());
            writer->SetSortToBSP();
            writer->SetFileFormatToPDF();
            writer->UsePainterSettings();
            writer->CompressOff();
            writer->DrawBackgroundOff();
            writer->TextAsPathOn();
            writer->Write3DPropsAsRasterImageOn();
            writer->SetFilePrefix(rawname.c_str());
            writer->Update();
            writer->Write();
        }


    }

private:

    typedef enum { // model type
        kMeshDisplay=0,
        kPeriostealDisplay,
        kVolumeDisplay,
    } displayOptions;

    void updateMeshCut() {

        double zMin, zMax, zIndexRatio;
        int zDisplayIndex, zMinIndex, zMaxIndex;

        zDisplayIndex = imageSliceMapper->GetSliceNumber();//imageActor->GetDisplayExtent()[4];
        zMinIndex = image->GetExtent()[4]; zMaxIndex = image->GetExtent()[5];

        zMin = image->GetBounds()[4]; zMax = image->GetBounds()[5];
        zIndexRatio = ((double) ( zDisplayIndex - zMinIndex )) / ((double) ( zMaxIndex - zMinIndex ));

        double planeOrigin[Dimension]; image->GetCenter(planeOrigin);

        // zDisplay = zMin + (zIndDisplay - zIndMin) / (zIndMax - zIndMin) * (zMin - zMax)
        planeOrigin[2] = zMin + zIndexRatio * (zMax - zMin);

        // cutting plane
        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        plane->SetNormal(0,0,1);
        plane->SetOrigin(planeOrigin);

        meshCutter->SetCutFunction(plane);
        meshCutter->GenerateValues(1, 0, zMax-planeOrigin[2]);

    }

    void updateMesh() {

        // set to initial mesh
        if(meshDisplaSelection==kMeshDisplay) { // set to periosteal mesh - must be measured
            mesh = corticalBone->getMesh();
        } else if(meshDisplaSelection==kPeriostealDisplay && corticalBone->isMeshMeasured()) { // if volume - leave unchanged just turn off visibility
            mesh = corticalBone->getPeriostealMesh();
        } else if(meshDisplaSelection==kVolumeDisplay) { // if unmeasure and periosteal mesh type - leave unchanged and print out an error
            // do nothing
            return;
        } else if(meshDisplaSelection==kPeriostealDisplay && !corticalBone->isMeshMeasured()) {
            cerr<<"Error in ThreeDimInteractorStyle::updateMesh() periosteal selected before it has been created."<<endl;
            return;
        } else { // unexpected combination. give a warning
            cerr<<"Error in ThreeDimInteractorStyle::updateMesh() unexpected combination."<<endl;
            return;
        }

        // remove previous mesh actor
        sliceRenderer->RemoveActor(meshActor);

        // setup the mesh actor
        setupMesh();


    }

    void setupMesh() { // requires mesh to already be set

        // create cutter
        meshCutter = vtkSmartPointer<vtkCutter>::New();
        meshCutter->SetInputData(mesh);

        // cutter mapper
        vtkSmartPointer<vtkPolyDataMapper> cutterMapper = vtkSmartPointer<vtkPolyDataMapper>::New();

        cutterMapper->SetInputConnection( meshCutter->GetOutputPort());
        cutterMapper->ScalarVisibilityOff();

        // actor - slice view
        meshActor = vtkSmartPointer<vtkActor>::New();
        meshActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
        meshActor->GetProperty()->SetInterpolationToFlat();
        meshActor->GetProperty()->SetLineWidth(3);
        meshActor->SetMapper(cutterMapper);

        if(displayFixedSlices) {
            setFixedMeshSlices();
        } else if(drawImage) {
            // create cutting plane and add to cutter
            updateMeshCut();
        } else {
            // generate a temporary cutting plane to avoid error
            double planeOrigin[Dimension]; mesh->GetCenter(planeOrigin);
            vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
            plane->SetNormal(0,0,1);  plane->SetOrigin(planeOrigin);
            meshCutter->SetCutFunction(plane);
        }

        sliceRenderer->AddActor( meshActor );

        //sliceRenderer->ResetCamera();
        ((vtkRenderWindow*)sliceView->GetRenderWindow())->Render();
        sliceView->update();
    }

    void initiliseSphereLocations() { // initialise sphere posistions based upon image bounds

        if(drawImage && calibrationPhantomType <= itk::CorticalBone::kManualControlPtsCal && calibrationPhantomType > itk::CorticalBone::kNoCal) {
            // get image bounds
            double imageBoundary[6];
            image->GetBounds(imageBoundary);
            double centre[] = {imageBoundary[0], imageBoundary[2], imageBoundary[4]};

            double radius; int number; corticalBone->getCalibrationPtGeometry(radius, number);

            for(int i=0; i<cubes.size();i++){ // spheres.size();
                cubes[i]->SetCenter(centre); cubes[i]->SetXLength(2*radius); // spheres, spheres->SetRadius(radius)
                cubes[i]->SetYLength(2*radius); cubes[i]->SetZLength(0.1);

                cubes[i]->Update(); // spheres
                actors[i]->VisibilityOff();

                // update centre posistion for next sphere
                centre[0] = centre[0] + image->GetLength()/50;

            }
        } else if(calibrationPhantomType <= itk::CorticalBone::kNoCal) {
            cerr<<"Error invalid calibrationPhantomIndex provided = "<<calibrationPhantomType<<" in Slice::initiliseSphereLocations"<<endl;
        }
    }

    void initialiseStates() {
        drawMesh=false; drawImage=false;

        calibrationOn=false; calibrationComplete = false;
        calibrationPointIndex=0; calibrationPhantomType =itk::CorticalBone::kNoCal;

        meshDisplaSelection = kMeshDisplay;
    }

    void showSpheres() {

        unsigned int numberOfSpheres=(unsigned int)cubes.size(); // spheres.size();

        for(int i=0; i<numberOfSpheres; i++) {
            actors[i]->VisibilityOn();
        }
        sliceRenderer->GetRenderWindow()->Render();
    }

    void hideSpheres() {

        unsigned int numberOfSpheres=(unsigned int)cubes.size(); //spheres.size();

        for(int i=0; i<numberOfSpheres; i++) {
            actors[i]->VisibilityOff();
        }
        sliceRenderer->GetRenderWindow()->Render();
    }

    void setFixedMeshSlices(){
        double opacity = 0.7; const int numberOfSlices=4;

        // NOTE with hard carsh if incorrect display extents used.
        //int spacing = 20; int displayExtents[6]={150,500,460,690,105,105}; // QCT - bern 22 left
        //int spacing = 200; int displayExtents[6]={200,828,0,1032,250,250}; // HRpQCT - 22L
        int spacing = 200; int displayExtents[6]={0,828,0,1132,250,250}; // HRpQCT - F14R


        // ------------set contoured mesh------------------- //
        meshCutter = vtkSmartPointer<vtkCutter>::New();
        meshCutter->SetInputData(mesh);

        // cutter mapper
        vtkSmartPointer<vtkPolyDataMapper> cutterMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        cutterMapper->SetInputConnection( meshCutter->GetOutputPort());
        //cutterMapper->ScalarVisibilityOff();

        // actor - slice view
        meshActor = vtkSmartPointer<vtkActor>::New();        meshActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
        meshActor->GetProperty()->SetInterpolationToFlat();  meshActor->GetProperty()->SetLineWidth(3);
        meshActor->SetMapper(cutterMapper);
        meshActor->SetVisibility(false);

        // get cutting range
        double zMin, zMax, zIndexRatio;
        int zDisplayIndex, zMinIndex, zMaxIndex;

        zDisplayIndex = displayExtents[4];
        zMinIndex = image->GetExtent()[4]; zMaxIndex = image->GetExtent()[5];

        zMin = image->GetBounds()[4]; zMax = image->GetBounds()[5];
        zIndexRatio = ((double) ( zDisplayIndex - zMinIndex )) / ((double) ( zMaxIndex - zMinIndex ));

        double planeOrigin[Dimension]; image->GetCenter(planeOrigin);

        planeOrigin[2] = zMin + zIndexRatio * (zMax - zMin);

        // cutting plane
        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        plane->SetNormal(0,0,1);
        plane->SetOrigin(planeOrigin);

        double zSliceEnd = ((double)((numberOfSlices-1)*spacing))/((double)(zMaxIndex-zMinIndex))*(zMax-zMin);
        meshCutter->SetCutFunction(plane);
        meshCutter->GenerateValues(numberOfSlices, 0, zSliceEnd);
    }

    void setFixedImageSlices() {
        double opacity = 0.7; int numberOfSlices=4;

        // NOTE with hard crash if incorrect display extents used.
        //int spacing = 20; int displayExtents[6]={150,500,460,690,105,105}; double window=460, level=230; // QCT - bern 22 left
        //int spacing = 200; int displayExtents[6]={200,828,0,1032,250,250}; double window=868, level=3341; // HRpQCT  -  bern 22L
        int spacing = -200; int displayExtents[6]={0,828,0,1132,250,250}; double window=500, level=3641; // HRpQCT  -  bern F14R

        double imageRange[2]; image->GetScalarRange(imageRange); // or could set manually

        // Stack
        imageStackActor = vtkSmartPointer<vtkImageStack>::New();

        // generate a image stack
        for(int i=0; i<numberOfSlices; i++) {

            vtkSmartPointer<vtkImageSliceMapper> imageSliceMapper = vtkSmartPointer<vtkImageSliceMapper>::New();
            imageSliceMapper->SetInputData(image);
            imageSliceMapper->SetSliceNumber(displayExtents[4]); imageSliceMapper->SetCroppingRegion(displayExtents);
            imageSliceMapper->CroppingOn(); displayExtents[4] = displayExtents[5] = displayExtents[5] + spacing;
            vtkSmartPointer<vtkImageSlice> imageSlice = vtkSmartPointer<vtkImageSlice>::New();
            imageSlice->SetMapper(imageSliceMapper); imageSlice->GetProperty()->SetOpacity(.5);
            imageSlice->GetProperty()->SetColorLevel(level); imageSlice->GetProperty()->SetColorWindow(window);
            //imageSlice->SetVisibility(false);

            // Add to stack
            imageStackActor->AddImage(imageSlice);
        }
        sliceRenderer->AddActor(imageStackActor);
        sliceRenderer->ResetCamera();
    }

    // image slice values
    vtkSmartPointer<vtkImageData> image;

    // mesh slicing values
    vtkSmartPointer<vtkPolyData> mesh;
    vtkSmartPointer<vtkActor> meshActor;
    vtkSmartPointer<vtkCutter> meshCutter;

    // state values
    bool drawMesh, drawImage, displayFixedSlices;
    bool calibrationComplete, calibrationOn;
    int calibrationPointIndex, calibrationPhantomType;
    int meshDisplaSelection;

    // slice displays - for making figures
    vtkSmartPointer<vtkImageStack> imageStackActor;
    vtkSmartPointer<ThreeDimInteractorStyle> threeDimInteractor;
    vtkSmartPointer<vtkImageSliceMapper> imageSliceMapper;

    // general values
    itk::CorticalBone::Pointer corticalBone;

    QVTKWidget* sliceView;

    vtkSmartPointer<vtkRenderer> sliceRenderer;

    // sphere display values
    std::vector<vtkSmartPointer<vtkSphereSource> > spheres;
    std::vector<vtkSmartPointer<vtkCubeSource> > cubes;
    std::vector<vtkSmartPointer<vtkActor> > actors;

};
//vtkCxxRevisionMacro(SliceInteractorStyle, "$Revision: 1.1 $");
vtkStandardNewMacro(SliceInteractorStyle);

// 2D display interactor - given a click wheel scroll / and click and scroll

QVisualisationFrame::QVisualisationFrame(QWidget *parent) : QFrame(parent) {

    corticalBone = itk::CorticalBone::New();
    corticalBone->setVerbose(true);
    createObjects();

}

bool QVisualisationFrame::reset() {


    SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
    sliceStyle->reset();
    ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
    threeDimStyle->reset();
    corticalBone->reset();
//  
//  SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
//  sliceStyle->initialiseSlice(corticalBone, threeDimView, sliceView, sliceRenderer, threeDimRenderer);
//  
//  ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
//  threeDimStyle->initialiseThreeDim(corticalBone, threeDimView, threeDimRenderer, profileView);
    return true;
}

bool QVisualisationFrame::setImage(QString fileName) {

    // vtk object for viewing, itk object for processing. initially only vtk.

    // reads a DICOM file - in cortical bone
    // create ITK object - add to corticalBone

    bool status = false;

    // check if .tiff and read in repository if so
    if (fileName.contains(".tif")) {

        // create a fileName copy to stop overwiting fileName
        QString fileNameCopy = QString(fileName.toStdString().c_str());


        // string away everything but the last number
        QRegularExpression pattern("(\\d+)(.tif)$");
        QString searchPattern = fileNameCopy.replace(pattern, "*.tif");
        searchPattern = searchPattern.section("/",-1);
        QString filePath = fileName.section("/",0,-2);

        QDirIterator directoryIterator(filePath, QStringList() << searchPattern, QDir::Files, QDirIterator::NoIteratorFlags);

        //std::cerr<<"filePath="<<filePath.toStdString().c_str()<<std::endl;
        //std::cerr<<"searchPattern="<<searchPattern.toStdString().c_str()<<std::endl<<std::endl;

        QStringList qNameList = QStringList(); fileNameCopy = QString(fileName.toStdString().c_str());
        QRegularExpression namePattern(fileNameCopy.replace(pattern, "(\\d+)(.tif)$"));

        //std::cerr<<"namePattern="<<namePattern.pattern().toStdString().c_str()<<std::endl<<std::endl;

        while(directoryIterator.hasNext()) {
            QString name = directoryIterator.next();
            //std::cerr<<"name="<<name.toStdString().c_str()<<std::endl;
            if(name.contains(namePattern)) { //namePattern.match(name).isValid()
                qNameList<<name;
            } else {
                //std::cerr<<"filtered out="<<name.toStdString().c_str()<<std::endl;
            }

        }
        qNameList.sort(); // sort so added in correct number order

        int startIndex = 0, endIndex = qNameList.size()-1; double zSpacing = nan("1");

        fileNameCopy = QString(fileName.toStdString().c_str());
        //QString seriesDetails = QString(fileNameCopy.toStdString().c_str()); // ensure no shared reference
        //seriesDetails.replace(pattern, "SeriesDetails.txt");
        QString seriesDetails = fileName.section("/",0,-3)+QDir::separator()+QString("SeriesDetails.txt");
        std::cerr<<"seriesDetails="<<seriesDetails.toStdString().c_str()<<std::endl;
        if(QFileInfo(seriesDetails).exists()) {

            QFile fileIn(seriesDetails);
            if(!fileIn.open(QIODevice::ReadOnly)) {
                cout<<"Error reading in .tif SeriesDetails .txt file. Load entire range."<<endl;
            } else { // read in details
                QTextStream streamIn(&fileIn);
                // read in start index
                QString line = streamIn.readLine();
                QStringList fields = line.split(" ");
                startIndex = fields.at(1).toInt();
                // read in end index
                line = streamIn.readLine();
                fields = line.split(" ");
                endIndex = fields.at(1).toInt();
                // read in end
                line = streamIn.readLine();
                fields = line.split(" ");
                zSpacing = fields.at(1).toDouble();

            }
        }

        std::vector<std::string> nameList(endIndex-startIndex+1);
        for (int i=0; i<=endIndex-startIndex; i++) {
            nameList[i]=qNameList.at(i+startIndex).toStdString();
            //cerr<<nameList.at(i).c_str()<<endl;
        }
        //cerr<<endl;

        double zOffset = startIndex*zSpacing;

        status = corticalBone->loadImage(nameList, zSpacing, zOffset);
    } else {
        // set display -> with vtk smartpointer returned from corticalBone
        status = corticalBone->loadImage(fileName.toStdString());
    }

    if(status) {
        Utilities::measureTime(true);
        ImageType::Pointer imageData = corticalBone->getDICOM();
        ConnectorType::Pointer connector = ConnectorType::New();
        connector->SetInput(imageData);  connector->Update();
        vtkSmartPointer<vtkImageData> image = connector->GetOutput();


        SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
        sliceStyle->setImage(image);// calls update image from three dim style
        cout<<Utilities::getTabString()<<Utilities::getTabString()<<Utilities::measureTime(false)<<" to create the image visualisation objects"<<endl;
    }


    return status;


}

bool QVisualisationFrame::setMesh(QString *fileName) {
    // vtk object for viewing, and processing. initially only visualisation


    // reads a Mesh file
    // add to corticalBone
    bool status = corticalBone->loadMesh(fileName->toStdString());

    // set display
    if(status) {
        SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
        sliceStyle->setMesh();
        ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
        threeDimStyle->setMesh();

        vtkSmartPointer<vtkPolyData> mesh = corticalBone->getMesh();
        threeDimStyle->setSphereRadius(corticalBone->getMesh()->GetLength()/100);

    }



    return status;
}

bool QVisualisationFrame::saveMesh(QString *fileName) {
    bool status = corticalBone->saveMesh(fileName->toStdString());
    return status;
}

bool QVisualisationFrame::saveValueArrays(QString *baseFileName) {

    // save data arrays with the base string followed by specific tag
    bool status = corticalBone->saveValueArrays(baseFileName->toStdString());

    return status;
}

bool QVisualisationFrame::saveCalibration(QString *fileName) {

    // save manually selected calibration control points + associated values
    bool status = corticalBone->saveCalibration(fileName->toStdString());

    return status;
}

bool QVisualisationFrame::loadValueArrays(QString fileName) {

    bool status = corticalBone->loadValueArrays(fileName.toStdString());

    return status;

}

bool QVisualisationFrame::saveDisplays(QString *sliceName, QString *profileName, QString *threeDimName) {

    // save the 3D, slice and chart
    SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
    ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
    sliceStyle->saveDisplays(sliceName->toStdString());
    threeDimStyle->saveDisplays(profileName->toStdString(), threeDimName->toStdString());

    return true;

}

bool QVisualisationFrame::openParameters(QString *fileName) {
    bool status = corticalBone->setImportProfiles(fileName->toStdString());

    if(status) {
        ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
        threeDimStyle->updateImportedProfiles();
    }
    return status;
}

void QVisualisationFrame::createObjects() {

    // create the overall layout
    QVBoxLayout *overallLayout = new QVBoxLayout(this);

    // create upper horizontal panel
    QHBoxLayout *topLayout = new QHBoxLayout();
    topLayout->addWidget(createSliceView()); // createSliceView())
    topLayout->addWidget(createThreeDimView());

    // create lower panel
    QVBoxLayout *bottomLayout = new QVBoxLayout();
    bottomLayout->addWidget(createProfileView());

    overallLayout->addLayout(topLayout);
    overallLayout->addLayout(bottomLayout);

    setLayout(overallLayout);

    // create rending styles
    createStyles();

}

void QVisualisationFrame::createStyles() {

    bool sliceDisplay = false; // true if generating figures for papers

    // create actors
    vtkSmartPointer<SliceInteractorStyle> sliceStyle = vtkSmartPointer<SliceInteractorStyle>::New();
    vtkSmartPointer<ThreeDimInteractorStyle> threeDimStyle = vtkSmartPointer<ThreeDimInteractorStyle>::New();

    // Setup slice interactor style
    sliceStyle->initialiseSlice(threeDimStyle, corticalBone, sliceView, sliceRenderer, sliceDisplay);


    vtkSmartPointer<vtkCellPicker> slicePicker = vtkSmartPointer<vtkCellPicker>::New();

    // Setup three dim interactor style
    threeDimStyle->initialiseThreeDim(corticalBone, threeDimView, threeDimRenderer, profileView, sliceDisplay);

    // add styles
    sliceView->GetInteractor()->SetInteractorStyle(sliceStyle);
    sliceView->GetInteractor()->SetPicker(slicePicker);

    threeDimView->GetInteractor()->SetInteractorStyle(threeDimStyle);

}

QDockWidget * QVisualisationFrame::createThreeDimView() {
    // create dock - todo consider removing docking capabilities
    QDockWidget *dock = new QDockWidget(tr("Three Dim"), this);
    dock->setAllowedAreas(Qt::RightDockWidgetArea);

    // create window
    threeDimView = new QVTKWidget(this); // todo see why example created without a pointer

    // Setup rendering
    threeDimRenderer = vtkSmartPointer<vtkRenderer>::New();

    // Add to display
    ((vtkRenderWindow*)threeDimView->GetRenderWindow())->AddRenderer(threeDimRenderer);
    //threeDimView->GetInteractor()->SetInteractorStyle(style);

    // show window
    threeDimView->show();

    // create dock and add to layout
    dock->setWidget(threeDimView);

    return dock;

}

QDockWidget * QVisualisationFrame::createSliceView() {

    // create dock - todo consider removing docking capabilities
    QDockWidget *dock = new QDockWidget(tr("Three Dim"), this);
    dock->setAllowedAreas(Qt::RightDockWidgetArea);
    //dock->setFeatures(QDockWidget::DockWidgetClosable);

    sliceView = new QVTKWidget(this); // todo see why example created without a pointer

    // Setup rendering
    sliceRenderer = vtkSmartPointer<vtkRenderer>::New();

    // Setup interactor style
    //vtkSmartPointer<SliceInteractorStyle> style = vtkSmartPointer<SliceInteractorStyle>::New();
    //style->initialiseSlice(corticalBone, threeDimView, sliceView);

    // Add to display
    ((vtkRenderWindow*)sliceView->GetRenderWindow())->AddRenderer(sliceRenderer);
    //sliceView->GetInteractor()->SetInteractorStyle(style);
    /*//renderWindowInteractor->SetRenderWindow(renderWindow);
    //renderWindowInteractor->Initialize();
    //renderWindowInteractor->Start();*/

    //sliceView->GetInteractor()->SetInteractorStyle(style);

    // show window
    sliceView->show();

    // create dock and add to layout
    dock->setWidget(sliceView);

    return dock;

}

QDockWidget * QVisualisationFrame::createProfileView() {

    // create dock - todo consider removing docking capabilities
    QDockWidget *dock = new QDockWidget(tr("Profile"), this);
    dock->setAllowedAreas(Qt::BottomDockWidgetArea);

    // create window
    profileView = new QVTKWidget(this);

    // Setup rendering
    //profileRenderer = vtkSmartPointer<vtkRenderer>::New();
    //profileView->GetRenderWindow()->AddRenderer(profileRenderer);

    // create data for testing
    // Create a table with some points in it
    /*vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    vtkSmartPointer<vtkDoubleArray> arrX = vtkSmartPointer<vtkDoubleArray>::New();
    arrX->SetName("X Axis");
    table->AddColumn(arrX);

    vtkSmartPointer<vtkDoubleArray> arrC = vtkSmartPointer<vtkDoubleArray>::New();
    arrC->SetName("Cosine");
    table->AddColumn(arrC);

    vtkSmartPointer<vtkDoubleArray> arrS = vtkSmartPointer<vtkDoubleArray>::New();
    arrS->SetName("Sine");
    table->AddColumn(arrS);*/

    //vtkSmartPointer<vtkTable> table = corticalBone->getProfileValues();

    /*// Fill in the table with some example values
    int numPoints = 69;
    double inc = 7.5 / (numPoints-1);
    table->SetNumberOfRows(numPoints);
    for (int i = 0; i < numPoints; ++i)  {
      table->SetValue(i, 0, i * inc);
      table->SetValue(i, 1, cos(i * inc));
      table->SetValue(i, 2, sin(i * inc));
    }*/

    //profileRenderer->SetBackground(1.0, 1.0, 1.0);

    //vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
    //vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
    //view->GetScene()->AddItem(chart);
    /*vtkPlot *line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 1);
    line->SetColor(0, 255, 0, 255);
    line->SetWidth(1.0);
    line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 2);

    line->SetColor(255, 0, 0, 255);
    line->SetWidth(1.0);*/

    //view->SetRenderWindow(profileView->GetRenderWindow());
    //profileView->GetRenderWindow()->AddRenderer(view->GetRenderer());
    //view->SetInteractor(profileView->GetInteractor());
    //profileView->SetRenderWindow(view->GetRenderWindow());
    //view->GetInteractor()->Initialize();
    //view->GetInteractor()->Start();

    //profileView->GetInteractor()->Initialize();
    //profileView->GetInteractor()->Start();

    // create dock and add to layout
    dock->setWidget(profileView);

    return dock;

}

void QVisualisationFrame::setupProfileView() {


}

void QVisualisationFrame::updateProfileView() {



    return ;

}

void QVisualisationFrame::disablePtMeasures() {
    ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
    threeDimStyle->disablePointMeasurements();
}

void QVisualisationFrame::enablePtMeasures() {
    ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
    threeDimStyle->enablePointMeasurements(); // todo return points

}

void QVisualisationFrame::enableCalibrationMode(){
    SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
    sliceStyle->startCalibrating();
}

void QVisualisationFrame::disableCalibrationMode(){
    SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
    sliceStyle->stopCalibrating();
}

bool QVisualisationFrame::runCalibration() {
    SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
    bool status = sliceStyle->runCalibration();

    return status;

}

bool QVisualisationFrame::setCalibration(double p0, double p1, double p2) {
    bool status = corticalBone->setCalibration(p0, p1, p2);
    return status;
}

bool QVisualisationFrame::runModellingOverMesh() {

    bool status = corticalBone->runModellingOverMesh();

    ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
    threeDimStyle->updateMeshColours();

    return status;
}

void QVisualisationFrame::setDisplayParameter(int index) {

    ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
    threeDimStyle->setDisplayParameter(index); // todo return points

}

void QVisualisationFrame::setDisplayMesh(int index) {

    ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
    threeDimStyle->setDisplayMesh(index);
    SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
    sliceStyle->setDisplayMesh(index);

}

void QVisualisationFrame::setCalibrationPhantom(int index) {
    corticalBone->setPhantomIndex(index);
    SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
    sliceStyle->setCalibrationMode(index);
}

void QVisualisationFrame::setModelFunction(int index) {
    corticalBone->setModelIndex(index);
}

void QVisualisationFrame::setOptimiser(int index) {
    corticalBone->setOptimiserIndex(index);
}

void QVisualisationFrame::setfittingScheme(int index) {
    corticalBone->setFittingSchemeIndex(index);
}


//------------- Fix CB Densities -------------//
void QVisualisationFrame::setFixedCBDensity(double CBDensity) {
    corticalBone->setFixedCBDensity(CBDensity);

}

void QVisualisationFrame::turnOnFWHMMode() {
    corticalBone->turnOnFWHMMode();
}

void QVisualisationFrame::turnOffFWHMMode() {
    corticalBone->turnOffFWHMMode();
}

bool QVisualisationFrame::isSetToFWHMMode() {
    return corticalBone->isSetToFWHMMode();
}

void QVisualisationFrame::removeFixedCBDensity() {
    corticalBone->removeFixedCBDensity();

}


//------------- Fix Sigma ------------------//
void QVisualisationFrame::turnOffFixedSigma() {
    corticalBone->turnOffFixedSigma();
}

void QVisualisationFrame::turnOnFixedSigma(double x, double y, double z) {
    corticalBone->turnOnFixedSigma(x, y, z);
}

bool QVisualisationFrame::isSigmaFixed() {
    return corticalBone->isSigmaFixed();
}

void QVisualisationFrame::getFixedSigma(double &x, double &y, double &z) {
    corticalBone->getFixedSigma(x, y, z);
}

//------------ Set Thresholds -----------------//
bool QVisualisationFrame::calculateThresholds(int percentIndex) {
    return corticalBone->calculateThresholds(percentIndex);
}

bool QVisualisationFrame::calculateThresholds(int percentIndex, double weight) {
    return corticalBone->calculateThresholds(percentIndex, weight);
}

bool QVisualisationFrame::setClassifierThresholdInfo(double stDensity, double cbDensity, double thresholdDensity,
                                                     int thresholdIndex) {
    corticalBone->setClassifierThresholdInfo(stDensity, cbDensity, thresholdDensity, thresholdIndex);
    return corticalBone->areThresholdsSet();
}

bool QVisualisationFrame::setClassifierThresholdInfo(double stDensity, double cbDensity, double thresholdDensity, double weight, int thresholdIndex) {
    corticalBone->setClassifierThresholdInfo(stDensity, cbDensity, thresholdDensity, weight, thresholdIndex);
    return corticalBone->areThresholdsSet();
}


bool QVisualisationFrame::getClassifierThresholdInfo(double &stDensity, double &cbDensity, double &thresholdDensity,
                                                     int &classifierIndex, QString& classifierThresholdIndexName) {
    std::string name;
    bool status = corticalBone->getClassifierInfo(stDensity, cbDensity, thresholdDensity, classifierIndex, name);
    classifierThresholdIndexName=QString(name.c_str());
    return status;
}


bool QVisualisationFrame::getClassifierThresholdInfo(double &stDensity, double &cbDensity, double &thresholdDensity, double &weight,
                                                     int &classifierIndex, QString& classifierThresholdIndexName) {
    std::string name;
    bool status = corticalBone->getClassifierInfo(stDensity, cbDensity, thresholdDensity, weight, classifierIndex, name);
    classifierThresholdIndexName=QString(name.c_str());
    return status;
}

void QVisualisationFrame::removeThresholds() {
    corticalBone->removeThresholds();
}

bool QVisualisationFrame::areThresholdsSet() {
    return corticalBone->areThresholdsSet();
}


//------------ Set HR Averaging --------//
void QVisualisationFrame::turnOnProfileAveraging(double x, double y, double z){
    corticalBone->turnOnProfileAveraging(x, y, z);

    ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
    threeDimStyle->setSphereRadius(std::max(x, std::max(y, z)));

}

void QVisualisationFrame::turnOffProfileAvergaing() {
    corticalBone->turnOffProfileAveraging();

    ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
    threeDimStyle->setSphereRadius(corticalBone->getMesh()->GetLength()/100);
}

//--------------- Fixed Sample Number -----------------//
void QVisualisationFrame::removeFixedSampledMode() {
    corticalBone->removeFixedSampleMode();
}

bool QVisualisationFrame::isSampleNumberFixed() {
    return corticalBone->isSampleNumberFixed();
}

int QVisualisationFrame::getSampleNumber() {
    return corticalBone->getSampleNumber();
}

bool QVisualisationFrame::setFixedSampleNumber(int sampleNumber) {
    return corticalBone->setFixedSampleNumber(sampleNumber);
}

//---------------Model Info -----------//
void QVisualisationFrame::getModelSelection(int& modelIndex, QString& modelName) {
    std::string localModelName;
    corticalBone->getModelSelection(modelIndex, localModelName);
    modelName=QString(localModelName.c_str());
}

void QVisualisationFrame::getOptimiserSelection(int& optimiserIndex, QString& optimiserName) {
    std::string localOptimiserName;
    corticalBone->getOptimiserSelection(optimiserIndex, localOptimiserName);
    optimiserName=QString(localOptimiserName.c_str());
}

void QVisualisationFrame::getSchemeSelection(int& schemeIndex, QString& schemeName){
    std::string localSchemeName;
    corticalBone->getSchemeSelection(schemeIndex, localSchemeName);
    schemeName = QString(localSchemeName.c_str());
}

//--------------- Smoothing ---------------------//
void QVisualisationFrame::setSmoothingRadius(double smoothingValue) {
    corticalBone->setSmoothingRadius(smoothingValue);
}

void QVisualisationFrame::turnOffSmoothingMode() {
    corticalBone->turnOffSmoothingMode();
}

void QVisualisationFrame::turnOnSmoothingMode() {
    corticalBone->turnOnSmoothingMode();
}

bool QVisualisationFrame::isSmoothingOn() {
    return corticalBone->isSmoothingOn();
}

double QVisualisationFrame::getSmoothingValue() {
    return corticalBone->getSmoothingValue();
}


//--------------- Set state ---------------//
int QVisualisationFrame::getState() {

    return corticalBone->getState();
}

vtkSmartPointer<vtkDoubleArray> QVisualisationFrame::getCalibrationPoints() {

    SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
    return sliceStyle->getCalibrationPts();
}


vtkSmartPointer<vtkDoubleArray> QVisualisationFrame::getCalibrationValues() {

    SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
    return sliceStyle->getCalibrationValues();
}

bool QVisualisationFrame::getCalibrationValues(double& p0, double& p1, double& p2) {
    bool status = corticalBone->getCalibrationValues(p0, p1, p2);
    return status;
}

bool QVisualisationFrame::getCalibrationPhantomType(std::string& phantomType, int& phantomIndex){
    bool status = corticalBone->getCalibrationPhantom(phantomType, phantomIndex);
    return status;
}

bool QVisualisationFrame::isCalibrated() {

    return corticalBone->isCalibrated();
}

bool QVisualisationFrame::setCalibrationPoints(int phantomIndex, vtkSmartPointer<vtkDoubleArray> pts) {

    SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
    bool status = sliceStyle->setCalibrationPts(phantomIndex, pts);
    return status;
}



bool QVisualisationFrame::getCalibrationPtGeometry(double &radius, int &number) {
    return corticalBone->getCalibrationPtGeometry(radius, number);
}

bool QVisualisationFrame::setCalibrationPtGeometry(double radius, int number) {
    bool status = corticalBone->setCalibrationPtGeometry(radius, number);

    int index; std::string name; corticalBone->getCalibrationPhantom(name, index);

    SliceInteractorStyle* sliceStyle = (SliceInteractorStyle*)sliceView->GetInteractor()->GetInteractorStyle();
    sliceStyle->setCalibrationMode(index);

    return status;
}

void QVisualisationFrame::removeImportedProfile() {
    corticalBone->removeImportedProfile();
    ThreeDimInteractorStyle* threeDimStyle = (ThreeDimInteractorStyle*)threeDimView->GetInteractor()->GetInteractorStyle();
    threeDimStyle->removeImportedProfiles();
}

bool QVisualisationFrame::isFixedCB() {

    return corticalBone->isCBFixed();
}

double QVisualisationFrame::getFixedCB() {
    return corticalBone->getFixedCB();
}

bool QVisualisationFrame::isProfileAveragingOn() {
    return corticalBone->isProfileAveragingOn();
}

bool QVisualisationFrame::getProfileAveragingValues(double &x, double &y, double &z) {
    return corticalBone->getProfileAveragingValues(x, y, z);
}


bool QVisualisationFrame::isImageSet() {
    return corticalBone->isImageSet();
}

bool QVisualisationFrame::isMeshSet() {
    return corticalBone->isMeshSet();
}

bool QVisualisationFrame::isParametersImported() {
    return corticalBone->isParametersImported();
}

bool QVisualisationFrame::areResultsLoaded() {
    return corticalBone->areResultsLoaded(); // todo - check if the project file or results have been renamed or removed
}

bool QVisualisationFrame::isMeshMeasured() {
    return corticalBone->isMeshMeasured();
}
