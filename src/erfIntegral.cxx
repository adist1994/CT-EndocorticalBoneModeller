#include <iostream>
#include <math.h>

#include "erfIntegral.h"


ERFIntLookup::ERFIntLookup(void) {
    integralType = -1;
    setupComplete = false;
}

void ERFIntLookup::setupTable(int selection, int lookupSize) {

    integralType = selection;
    intNumPoints = lookupSize;
    erfIntegralArray = itk::Array<double>(intNumPoints);

    intInc = (intEnd-intStart) / (intNumPoints-1);

    double erfSum = 0;

    // --- generate erf integral ---//

    if(integralType==RectangleRule) {

        for (int i = 0; i < intNumPoints; ++i)
        {
            double x = i * intInc + intStart; // position
            erfIntegralArray[i] = erfSum*intInc;

            erfSum += erf(x); // rectangle rule - update sum
        }
        setupComplete = true;

    } else if(integralType==TrapizoidalRule) {

        for (int i = 0; i < intNumPoints; ++i)
        {
            double x = i * intInc + intStart; // position
            erfIntegralArray[i] = erfSum*intInc;//

            erfSum += 0.5*(erf(x) + erf(x+intInc)); // trapezoidal rule - update sum
        }
        setupComplete = true;

    } else if(integralType==SimpsonsRule) {

        for (int i = 0; i < intNumPoints; ++i)
        {
            double x = i * intInc + intStart; // position
            erfIntegralArray[i] = erfSum*intInc;//

            erfSum += 1.0/6.0*(erf(x) + 4.0 * erf(x+0.5*intInc) + erf(x+intInc)); // simpson's rule - update sum
        }
        setupComplete = true;

    } else {
        std::cerr<<"Incorrect selection value for the type of integration to use"<<std::endl;
    }

}

double ERFIntLookup::erfIntegralMethod(double x) {

    if(!setupComplete) {
        std::cerr<<"Lookup table is not set up either invalid integral type selected, or 'setupTable not called."<<std::endl;
        return -1;
    }

    double xIndexDouble = (x -intStart) / (intEnd - intStart) * (intNumPoints-1);
    int xIndex = int(xIndexDouble);
    double remainder = xIndexDouble - xIndex;

    double erfIntegralValue = (xIndex <= 0)
                              ? (-1.0*x + intStart)
                              : ( (xIndex < intNumPoints-1)
                                  ? (1- remainder) * erfIntegralArray[xIndex] + (remainder * erfIntegralArray[xIndex+1])
                                  : (1.0*x + intStart) );

    return erfIntegralValue;

}