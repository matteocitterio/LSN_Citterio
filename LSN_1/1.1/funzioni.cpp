#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <armadillo>
#include "funzioni.h"
#include "random.h"

using namespace std;
using namespace arma;

double Error (vec averages, vec squared, int n){

    if (n ==0){
        return 0;
    }

    else{
        return sqrt((squared[n] - pow(averages[n], 2)) / n);
    }

}

//std::pair is needed because i want to return two arma::vec objects

std::pair<arma::vec, arma::vec> ProgressiveStatistics(string flag, Random &rnd, int N, int L){

        if (flag != "averages" && flag != "uncertainty")
    {
        cerr << "PROBLEM: accepted arguments: 'averages' or 'uncertainty' " << endl;
        return std::make_pair(arma::vec(), arma::vec());
    }

    vec averages(N, fill::zeros);
    vec squared(N, fill::zeros);
    vec cumulative_averages(N, fill::zeros);
    vec cumulative_squared(N, fill::zeros);
    vec cumulative_error(N, fill::zeros);

    double temp_average;

    // Double cycle: the outer one cycles on every block, the second cycles on every throw within a block
    for (int j = 0; j < N; j++){

        temp_average = 0;

        if (flag == "averages"){
            for (int i = 0; i < L; i++){
                temp_average += rnd.Rannyu();
            }
        }

        else {
            for (int i = 0; i < L; i++){
                temp_average += pow((rnd.Rannyu() - 0.5), 2);
            }
        }

        temp_average /= L;
        averages(j) = temp_average;
        squared(j) = temp_average * temp_average;
    }

    // Getting the progressive statistics

    double runningSum = 0;
    double runningSquared = 0;

    for (int i = 0; i < int(averages.size()); i++){

        runningSum += averages[i];
        runningSquared += squared[i];
        cumulative_averages(i) = (runningSum / (i + 1));
        cumulative_squared(i) = (runningSquared / (i + 1));
        cumulative_error(i) = (Error(cumulative_averages, cumulative_squared, i));
    }

    return std::make_pair(cumulative_averages, cumulative_error);
}
