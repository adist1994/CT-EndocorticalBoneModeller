//
// Created by rap58 on 31/08/15.
//

#ifndef MYPROJECT_UTILITIES_H
#define MYPROJECT_UTILITIES_H

#include <iostream>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <itkArray.h>

class Utilities {


public:
    static std::string getTabString() { static std::string tabString("&nbsp;&nbsp;&nbsp;&nbsp;"); return tabString; }
    static std::string getSpacesString(int length);
    static std::string measureTime(bool startMeasure);
    static vtkSmartPointer<vtkDoubleArray> convertItkToVtkArray(itk::Array<double> arrayIn, std::string name);
    static itk::Array<double> convertVtkToItkArray(vtkSmartPointer<vtkDoubleArray> arrayIn);
    static short swapShort(short input);
    static long swapLong(long input);
    static vtkTypeUInt32 swap32(vtkTypeUInt32 input);
    static vtkTypeFloat32 swap32(vtkTypeFloat32 input);
    static vtkTypeUInt16 swap16(vtkTypeUInt16 input);

    static double CalculateArrayMedian(itk::Array<double> array, int startIndex, int endIndex);
};


#endif //MYPROJECT_UTILITIES_H
