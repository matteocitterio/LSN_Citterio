#ifndef __Integral__
#define __Integral__

//importing libraries
#include <cmath>
#include <iostream>
#include "functions.h"
#include "random.h"

using namespace std;

class Integral {
    
    public:
        Integral(double a, double b, Functions* f, Random* rnd);
        double arithAverage (unsigned int N);
        double Imp_sampling (unsigned int N, double d_max, Functions* d);

    private:
        double I_min;
        double I_max;
        double I_integral;
        double I_sign;
        Functions* I_f;
        Random* I_rnd;
};


#endif //__Integral__