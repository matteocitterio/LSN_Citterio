#include "integral.h"

using namespace std;

//constructor
Integral :: Integral (double a, double b, Functions* f, Random* rnd){
    I_min = min(a, b);
    I_max = max(a, b);
    I_f = f;
    I_rnd = rnd; 
    if (a>b) I_sign = -1.;
    else I_sign = 1.;
    I_integral = 0;
}

//uniform sampling
double Integral::arithAverage(unsigned int N){
    I_integral = 0;
    for (int i = 0; i < N; i++){
        I_integral += I_f -> Evaluate(I_rnd -> Rannyu(I_min, I_max));
    }
    I_integral = I_sign * (I_max - I_min) * I_integral / (double) N; 
    return I_integral;
}

//Importance sampling
double Integral::Imp_sampling(unsigned int N, double d_max, Functions *d){
    I_integral = 0;
    for (int i = 0; i < N; i++){
        double x = I_rnd -> d_prob(I_min, I_max, d_max, d);
        I_integral += (I_f->Evaluate(x)) / (d->Evaluate(x));
    }
    I_integral = I_sign * (I_max - I_min) * I_integral / (double)N;
    return I_integral;
}