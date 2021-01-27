//
// Created by rap58 on 31/08/15.
//

#include <utilities.h>

#include <chrono>
#include <vtkSortDataArray.h>

std::string Utilities::getSpacesString(int length) {

    std::string spaceString("");
    for(int i=0; i<length; i++) {
        spaceString = spaceString + "&nbsp;";
    }
    return spaceString;
}

using namespace  std::chrono;
std::string Utilities::measureTime(bool startMeasure) {

    static high_resolution_clock::time_point startTime, finishTime;

    if(startMeasure) {
        startTime = high_resolution_clock::now();

        return std::string();

    } else {
        finishTime = high_resolution_clock::now();

        double durationMilliseconds = duration_cast<milliseconds>( finishTime - startTime ).count()%100;
        double durationSeconds = duration_cast<seconds>( finishTime - startTime ).count();

        std::ostringstream durationStream;
        durationStream << durationSeconds << "." << durationMilliseconds << "s";
        return durationStream.str();

    }


}

vtkSmartPointer<vtkDoubleArray> Utilities::convertItkToVtkArray(itk::Array<double> arrayIn, std::string name) {
    int n = arrayIn.Size();
    vtkSmartPointer<vtkDoubleArray> arrayOut = vtkSmartPointer<vtkDoubleArray>::New();
    arrayOut->SetNumberOfValues(n); arrayOut->SetName(name.c_str());
    for(vtkIdType i=0; i<n; i++) {
        arrayOut->SetValue(i, arrayIn[i]);
    }
    return arrayOut;
}

itk::Array<double> Utilities::convertVtkToItkArray(vtkSmartPointer<vtkDoubleArray> arrayIn) {

    vtkIdType n = arrayIn->GetNumberOfTuples();
    itk::Array<double> arrayOut = itk::Array<double>(n);

    for(vtkIdType i=0; i<n; i++) {
        arrayOut[i] = arrayIn->GetValue(i);
    }
    return arrayOut;
}

short Utilities::swapShort(short input) {
    unsigned char b1, b2;

    b1 = (unsigned char)(input & 255);
    b2 = (unsigned char)((input >> 8) & 255);

    return (b1 << 8) + b2;

}

long Utilities::swapLong (long input)
{
    unsigned char b1, b2, b3, b4;

    b1 = (unsigned char)(input & 255);
    b2 = (unsigned char)(( input >> 8 ) & 255);
    b3 = (unsigned char)(( input>>16 ) & 255);
    b4 = (unsigned char)(( input>>24 ) & 255);

    return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
}

vtkTypeUInt32 Utilities::swap32(vtkTypeUInt32 input) {
    unsigned char b1, b2, b3, b4;

    b1 = (unsigned char)(input & 255);
    b2 = (unsigned char)(( input >> 8 ) & 255);
    b3 = (unsigned char)(( input>>16 ) & 255);
    b4 = (unsigned char)(( input>>24 ) & 255);

    return ((vtkTypeUInt32)b1 << 24) + ((vtkTypeUInt32)b2 << 16) + ((vtkTypeUInt32)b3 << 8) + (vtkTypeUInt32)b4;
}

vtkTypeFloat32 Utilities::swap32(vtkTypeFloat32 input) {
    unsigned char b1, b2, b3, b4;

    b1 = (unsigned char)((vtkTypeUInt32)input & 255);
    b2 = (unsigned char)(( (vtkTypeUInt32)input >> 8 ) & 255);
    b3 = (unsigned char)(( (vtkTypeUInt32)input>>16 ) & 255);
    b4 = (unsigned char)(( (vtkTypeUInt32)input>>24 ) & 255);

    return vtkTypeFloat32( ((vtkTypeUInt32)b1 << 24) + ((vtkTypeUInt32)b2 << 16) + ((vtkTypeUInt32)b3 << 8) + b4 );
}

vtkTypeUInt16 Utilities::swap16(vtkTypeUInt16 input) {
    unsigned char b1, b2;

    b1 = (unsigned char)(input & 255);
    b2 = (unsigned char)((input >> 8) & 255);

    return (b1 << 8) + b2;
}

double Utilities::CalculateArrayMedian(itk::Array<double> array, int startIndex, int endIndex) {
    int n = endIndex - startIndex + 1;

    if(n<=0) {
        return nan("1");
    }

    vtkSmartPointer<vtkDoubleArray> sortArray = vtkSmartPointer<vtkDoubleArray>::New();
    sortArray->SetNumberOfValues(n); // todo - consider a itk only solution

    for(int i=0; i<n; i++) {
        sortArray->SetValue(i,array[i+startIndex]);
    }

    vtkSortDataArray::Sort(sortArray);

    return 0.5 *(sortArray->GetValue(floor((n-1)/2.0)) + sortArray->GetValue(floor(n/2.0)));
}