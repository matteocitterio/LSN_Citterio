#ifndef __Functions__
#define __Functions__
#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

double Error(vec averages, vec squared, int n);
void SetRandomGen(Random rnd);
std::pair<arma::vec, arma::vec> ProgressiveStatistics(string flag, Random &rnd, int N, int L);

#endif //__Functions__