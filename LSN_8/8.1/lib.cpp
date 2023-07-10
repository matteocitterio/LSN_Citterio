#include "lib.h"

double Error(double averages, double squared, int n) {

    /*
    This computes the blocking average error
    */

    if (n ==0){                                                                             // If it is the first block return 0
        return 0;
    }

    else{
        return sqrt((squared - pow(averages, 2)) / (double) n);
    }
}

void MetropolisUniform(vec &initialPositions, Functions *f, Random *rnd, double c, double &accepted, double &attempted) {

    /*
    Metropolis algorithm for uniform transition probability.
    Inputs:
    - vec & initialPositions: vector of positions passed by reference. If the Metropolis is successful, the vector will be modified
    - Functions* f: pointer to the function used for sampling
    - Random* rnd: pointer to the object of the random generator class
    - double c: tunes the width of the step in order to achieve a desired acceptance rate
    - double &accepted: counter used for computing the acceptance rate in each block
    - double & attempted: counter used for computing the acceptance rate in each block
    */

    arma::vec tempPositions(3);

    for (int i = 0; i < 3; i ++ ) {
        tempPositions[i] = initialPositions[i] + rnd->Rannyu(-1,1) * c;                     // Uniform transition probability
    }

    double alpha = min (1., (f->Evaluate(tempPositions)/f->Evaluate(initialPositions)));    // Compute the metropolis rate

    // Accepting the new configuration with probability \alpha:
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

void MetropolisUniform(double &initialPosition, Functions *f, Random *rnd, double c, double &accepted, double &attempted) {

    /*
    Metropolis algorithm for uniform transition probability 1D.
    Inputs:
    - double & initialPositions: 1D-position passed by reference. If the Metropolis is successful, it will be modified.
    - Functions* f: pointer to the function used for sampling
    - Random* rnd: pointer to the object of the random generator class
    - double c: tunes the width of the step in order to achieve a desired acceptance rate
    - double &accepted: counter used for computing the acceptance rate in each block
    - double & attempted: counter used for computing the acceptance rate in each block
    */

    double tempPosition = initialPosition + rnd->Rannyu(-1, 1) * c;                         // Uniform transition probability

    double alpha = min(1., (f->Evaluate(tempPosition) / f->Evaluate(initialPosition)));     // Compute the metropolis acceptance

    // accepting the configuration with probability \alpha:
    double p = rnd->Rannyu();
    if (p < alpha)
    {

        initialPosition = tempPosition;

        accepted++;
    }
    attempted++;
}

void MetropolisGauss(vec &initialPositions, Functions *f, Random *rnd, double c, double &accepted, double &attempted) {

    /*
    Metropolis algorithm for gaussian transition probability.
    Inputs:
    - vec & initialPositions: vector of positions passed by reference. If the Metropolis is successful, the vector will be modified
    - Functions* f: pointer to the function used for sampling
    - Random* rnd: pointer to the object of the random generator class
    - double c: tunes the width of the step in order to achieve a desired acceptance rate
    - double &accepted: counter used for computing the acceptance rate in each block
    - double & attempted: counter used for computing the acceptance rate in each block
    */

    arma::vec tempPositions(3);

    for (int i = 0; i < 3; i++)
    {
        tempPositions[i] = initialPositions[i] + rnd->Gauss(0, c) ;                         // Gaussian transition probability
    }

    double alpha = min(1., (f->Evaluate(tempPositions) / f->Evaluate(initialPositions)));   // Metropolis rate

    // Accepting the new configuration with probability \alpha:
    double p = rnd->Rannyu();
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

void EquilibrateUN(int nblocks,int L, vec &initialPositions, Functions *f, Random *rnd, double c) {

    /*
    Runs Metropolis wit a uniform transition probability for a number of blocks without measuring properties allowing the system to 'thermalize' towards equilibrium.
    Inputs:
    - int nblocks: number of blocks used for equilibration
    - int L: number of desidered throws in each block
    - vec &intialPositions: vector of positions passed by reference. If the Metropolis is successful, the vector will be modified
    - Functions* f: pointer to the function used for sampling
    - Random* rnd: pointer to the object of the random generator class
    - double c: tunes the width of the step in order to achieve a desired acceptance rate
    */

    double accepted = 0.;                                                                    // Actually i don't really care about computing the acceptance rate, i simply dont want to modify the MetropolisUN function
    double attempted = 0.;
    for (int j =0; j < nblocks; j++) {                                                       // Loop over the blocks
        for (int i = 0; i < L; i ++) {                                                       // Loop over the throws
            MetropolisUniform(initialPositions, f, rnd, c, accepted, attempted);
        }
    }
}

void EquilibrateUN(int nblocks, int L, double &initialPosition, Functions *f, Random *rnd, double c) {

    /*
    1D overload of the previous:: here I accept as input a 1D position instead of a 3D vector
    Inputs:
    - int nblocks: number of blocks used for equilibration
    - int L: number of desidered throws in each block
    - double &intialPositions: 1D.position passed by reference. If the Metropolis is successful, it will be modified
    - Functions* f: pointer to the function used for sampling
    - Random* rnd: pointer to the object of the random generator class
    - double c: tunes the width of the step in order to achieve a desired acceptance rate
    */

    double accepted = 0.;
    double attempted = 0.;
    for (int j = 0; j < nblocks; j++)
    {
        for (int i = 0; i < L; i++)
        {
            MetropolisUniform(initialPosition, f, rnd, c, accepted, attempted);
        }
    }
}

void EquilibrateGA(int nblocks, int L, vec &initialPositions, Functions *f, Random *rnd, double c) {

    /*
    Runs Metropolis wit a gaussian transition probability for a number of blocks without measuring properties allowing the system to 'thermalize' towards equilibrium.
    Inputs:
    - int nblocks: number of blocks used for equilibration
    - int L: number of desidered throws in each block
    - vec &intialPositions: vector of positions passed by reference. If the Metropolis is successful, the vector will be modified
    - Functions* f: pointer to the function used for sampling
    - Random* rnd: pointer to the object of the random generator class
    - double c: tunes the width of the step in order to achieve a desired acceptance rate
    */

    double accepted = 0.;                                                                    // Actually i don't really care about computing the acceptance rate, i simply dont want to modify the MetropolisUN function
    double attempted = 0.;
    for (int j =0; j < nblocks; j++) {                                                       // Loop over the blocks
        for (int i = 0; i < L; i++) {                                                        // Loop over the throws
            MetropolisGauss(initialPositions, f, rnd, c, accepted, attempted);
        }
    }
}