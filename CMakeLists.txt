cmake_minimum_required(VERSION 2.8.8)

project(MyProject)
subdirs(src)

if(TARGET vtkGUISupportQtOpenGL)
  add_subdirectory(GraphicsView)
endif()

#CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.ctest.in"
#  "${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.ctest" @ONLY)
