#include "integral.h"

using namespace std;

Integral :: Integral (double a, double b, Functions* f, Random* rnd){

    /*
    Class Constructor
    */

    I_min = min(a, b);                                                                  // Define the domain limits
    I_max = max(a, b);
    I_f = f;                                                                            // Data member representing the function
    I_rnd = rnd;                                                                        // Data member for the random seed
    if (a>b) I_sign = -1.;      
    else I_sign = 1.;                                                                   // Integral sign data member
    I_integral = 0;                                                                     // Value of the integral
}

double Integral::arithAverage(unsigned int N){

    /*
    Uniform sampling technique.
    Inputs:
    - unigned int N: size of the partition used for the integral evaluation
    */

    I_integral = 0;                                                                     // Reinitialize to zero the value

    for (int i = 0; i < N; i++){

        I_integral += I_f -> Evaluate(I_rnd -> Rannyu(I_min, I_max));                   // Evaluate the function over the partiotion's point

    }

    I_integral = I_sign * (I_max - I_min) * I_integral / (double) N;                    // compute average, add the sign

    return I_integral;
}

double Integral::Imp_sampling(unsigned int N, double d_max, Functions *d){

    /*
    Importance sampling technique.
    Inputs:
    - unigned int N: size of the partition used for the integral evaluation
    - double d_max: maximum of the function within the domain
    - Functions *d: pointer to the function used for importance sampling
    */

    I_integral = 0;                                                                     // Reinitialize to zero the value

    for (int i = 0; i < N; i++){

        double x = I_rnd -> d_prob(I_min, I_max, d_max, d);                             // Sample a point according the d
        I_integral += (I_f->Evaluate(x)) / (d->Evaluate(x));                            // Evaluate the function
    }

    I_integral = I_sign * (I_max - I_min) * I_integral / (double)N;                     // compute the mean and add the sign

    return I_integral;
}