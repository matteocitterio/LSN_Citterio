#include "lib.h"

//Function for the blocking average uncertainty
double Error(double averages, double squared, int n){

    if (n ==0){
        return 0;
    }

    else{
        return sqrt((squared - pow(averages, 2)) / (double) n);
    }
}

//Metropolis algorithm for uniform transition probability
void MetropolisUniform(vec &initialPositions, Functions *f, Random *rnd, double c, double &accepted, double &attempted) {

    arma::vec tempPositions(3);

    for (int i = 0; i < 3; i ++ ){
        tempPositions[i] = initialPositions[i] + rnd->Rannyu(-1,1) * c;                   //uniform transition probability
    }

    double alpha = min (1., (f->Evaluate(tempPositions)/f->Evaluate(initialPositions)));

    //accepting the configuration with probability \alpha:
    double p = rnd->Rannyu();
    if (p < alpha) {

        for (int i = 0; i < 3; i++)
        {
            initialPositions[i] = tempPositions[i];
        }
        accepted ++;
    }
    attempted ++;
}

void MetropolisGauss(vec &initialPositions, Functions *f, Random *rnd, double c, double &accepted, double &attempted)
{

    arma::vec tempPositions(3);

    for (int i = 0; i < 3; i++)
    {
        tempPositions[i] = initialPositions[i] + rnd->Gauss(0, c) ;     //Gaussian transition probability
    }

    double alpha = min(1., (f->Evaluate(tempPositions) / f->Evaluate(initialPositions)));
    // cout << "alpha " << alpha << endl;

    // accepting the configuration with probability \alpha:
    double p = rnd->Rannyu();
    // cout << "p " << p << endl;
    if (p < alpha)
    {

        for (int i = 0; i < 3; i++)
        {
            initialPositions[i] = tempPositions[i];
        }
        accepted++;
    }
    attempted++;
}

//run a block of throws without measuring props
void EquilibrateUN(int nblocks,int L, vec &initialPositions, Functions *f, Random *rnd, double c) {
    double accepted = 0.;
    double attempted = 0.;
    for (int j =0; j < nblocks; j++){ 
        for (int i = 0; i < L; i ++){
            MetropolisUniform(initialPositions, f, rnd, c, accepted, attempted);
        }
    }
}

// run a block of throws without measuring props
void EquilibrateGA(int nblocks, int L, vec &initialPositions, Functions *f, Random *rnd, double c)
{
    double accepted = 0.;
    double attempted = 0.;
    for (int j = 0; j < nblocks; j++)
    {
        for (int i = 0; i < L; i++)
        {
            MetropolisGauss(initialPositions, f, rnd, c, accepted, attempted);
        }
    }
}