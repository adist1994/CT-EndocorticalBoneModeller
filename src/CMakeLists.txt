cmake_minimum_required(VERSION 2.8.8) # checked version on uni machine 6/8/15 - was 2.8.12.2
if(POLICY CMP0020)
    cmake_policy(SET CMP0020 NEW)
endif()
if(POLICY CMP0025)
    cmake_policy(SET CMP0025 NEW) # CMake 3.0
endif()
if(POLICY CMP0053)
    cmake_policy(SET CMP0053 NEW) # CMake 3.1
endif()

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

# Find VTK.
find_package(VTK COMPONENTS REQUIRED)
include(${VTK_USE_FILE})

# Messages
#message("ITKPackage:${ITK_LIBRARIES}")
#message("VTKPackage:${VTK_LIBRARIES}")

# include Qt
set(CMAKE_AUTOMOC ON)
find_package(Qt5Widgets REQUIRED)


if("${VTK_QT_VERSION}" STREQUAL "")
    message(FATAL_ERROR "VTK was not built with Qt")
endif()

# Set your files and resources here
set( Srcs
        main.cxx
        mainwindow.cxx
        qdebugstream.cxx
        scriptframe.cxx
        controlframe.cxx
        visualisationframe.cxx
        corticalbone.cxx
        xmlFileAccess.cxx
        erfIntegral.cxx
        itkModelTransformBaseLocal.hxx
        itkThreeTierRectangularTransformLocal.hxx
        itkEndostealRampTransformLocal.hxx
        itkProfileSpatialObjectLocal.hxx
        itkSingleModelMetricBaseLocal.hxx
        itkSingleModelMetricLocal.hxx
        itkMultipleModelMetricBaseLocal.hxx
        itkMultipleModelMetricLocal.hxx
        vtkOBJWriter.cxx # off the Internet. Information in file header
        itkModelRegistrationMethodLocal.hxx
        profileClassifier.cxx
        linearTransform.cxx
        classifierTransform.cxx
        utilities.cxx
        batchProcessor.cxx
        colourMap.cxx
        LinearRegressionCalculator.cxx
        ProfileProcessor.cxx
        ProfileCalibration.cxx  itkQCTImageIO.hxx)

set( Hdrs
        mainwindow.h
        qdebugstream.h
        scriptframe.h
        controlframe.h
        visualisationframe.h
        corticalbone.h
        xmlFileAccess.h
        erfIntegral.h
        itkModelTransformBaseLocal.h
        itkThreeTierRectangularTransformLocal.h
        itkEndostealRampTransformLocal.h
        itkProfileSpatialObjectLocal.h
        itkSingleModelMetricBaseLocal.h
        itkSingleModelMetricLocal.h
        itkOptimisationCallback.h
        itkMultipleModelMetricBaseLocal.h
        itkMultipleModelMetricLocal.h
        vtkOBJWriter.h # off the Internet. Information in file header
        itkModelRegistrationMethodLocal.h
        profileClassifier.h
        linearTransform.h
        classifierTransform.h
        utilities.h
        batchProcessor.h
        colourMap.h
        LinearRegressionCalculator.h
        ProfileProcessor.h
        ProfileCalibration.h itkQCTImageIO.h)

set( Rcs
        application.qrc )

include_directories(
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# CMAKE_AUTOMOC in ON so the MocHdrs will be automatically wrapped.
set(CMAKE_AUTOMOC ON)

qt5_wrap_ui(UI_Srcs ${UIs}) # runs uic on .ui files and saves as US_Srcs
qt5_ADD_RESOURCES( RC_Srcs ${Rcs} ) # run rcc on Rcs saving all paths into RC_Srcs

add_executable(MyProject ${Srcs} ${Hdrs} ${UI_Srcs} ${RC_Srcs} ${MOC_Hdrs})
qt5_use_modules(MyProject Core Gui Widgets)

# perform linking
target_link_libraries(MyProject ${VTK_LIBRARIES} ${ITK_LIBRARIES})

add_definitions(-std=c++11) # added to support clion for mac

