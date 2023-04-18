#ifndef __Library__
#define __Library__

#include <cmath>
#include <armadillo>
#include "functions.h"
#include "random.h"

using namespace std;
using namespace arma;

double Error(double averages, double squared, int n);
void MetropolisUniform(vec &initialPositions, Functions *f, Random *rnd, double c, double &accepted, double &attempted);
void MetropolisGauss(vec &initialPositions, Functions *f, Random *rnd, double c, double &accepted, double &attempted);
void EquilibrateUN(int nblocks, int L, vec &initialPositions, Functions *f, Random *rnd, double c);
void EquilibrateGA(int nblocks, int L, vec &initialPositions, Functions *f, Random *rnd, double c);

#endif //__Library__