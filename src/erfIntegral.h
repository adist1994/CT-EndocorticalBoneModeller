/* 
 * File:   erfIntegral.h
 * Author: rap58
 *
 * Created on 30 December 2014, 12:16
 */

// precomputes and then returns the values of the erf integral between -3.5 and 3.5 at intervals of 10e-3

#ifndef ERFINTEGRAL_H
#define	ERFINTEGRAL_H

#include <itkArray.h>

class vtkDoubleArray;

class ERFIntLookup {

public:
    double erfIntegralMethod(double x);

    ERFIntLookup(void);

    void setupTable(int selection = RectangleRule, int lookupSize = 10001);

    enum { // type to numerical integration
        RectangleRule = 0,
        TrapizoidalRule,
        SimpsonsRule,
    } integrationType;

private:
    static constexpr double intStart = -3.5, intEnd = 3.5;
    int intNumPoints;
    double intInc;

    int integralType;
    bool setupComplete;

    itk::Array<double> erfIntegralArray;


};



#endif	/* ERFINTEGRAL_H */