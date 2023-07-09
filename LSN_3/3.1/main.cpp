#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "lib.h"
#include "random.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    // Initialize an element of the Random class
    Random *rnd = new Random();
    rnd->RandomRoutine();                                                   // Random settings including setting the seed, ....

    // M thows, N blocks, L throws in each block
    int M = pow(10,6);      
    int N= 1000;
    int L = int(M / N);

    cout << "Working with " << L << " throws in each block" << endl;

    // Define variables of the pricing process
    double r = 0.1;                                                                         // Risk-free interest rate                                       
    double T = 1.;                                                                          // Maturity
    double sigma = 0.25;                                                                    // Volatility
    double K = 100;                                                                         // Strike price
    double S_0 = 100;                                                                       // Initial asset value

    // 3.1 - Direct sampling

    ofstream callFile;
    ofstream putFile;
    callFile.open("3.1.call.txt");
    putFile.open("3.1.put.txt");

    // Define all the auxiliary variables needed for the computation
    double x = 0.;
    double S_T = 0;
    double call = 0;
    double put = 0;
    double runningSumCall = 0.;
    double runningSquaredCall = 0.;
    double runningSumPut = 0.;
    double runningSquaredPut = 0.;
    double errorCall = 0.;
    double errorPut = 0.;


    for (int i = 0; i < N; i++){                                                            // Iterate over the blocks

        call = 0;
        put = 0;

        for (int j = 0; j < L; j++){                                                        // Iterate over the throws in each block
            x= rnd->Gauss(0, T);
            S_T = S_0 * exp( (r - sigma * sigma * 0.5 ) * T + (sigma * x ) );               // Compute the value of the asset at maturity using B&S model
            call += std::max(0., S_T - K) * exp (- r * T);                                  // Compute call price accordingly
            put += std::max(0., K - S_T)* exp (- r * T);                                    // Compute put price accordingly
        }

        call /= L;                                                                          // Average over the number of throws                                  
        put /= L;

        // Call - data blocking

        runningSumCall = ((runningSumCall * i) + call) / (i + 1);                           // at every iteration I have to re-multiplicate the previous division
        runningSquaredCall = ((runningSquaredCall * i) + pow(call, 2)) / (i + 1);
        errorCall = Error(runningSumCall, runningSquaredCall, i);
        callFile << runningSumCall << " " << errorCall << endl;

        // Put - data blocking
        runningSumPut = ((runningSumPut * i) + put) / (i + 1);                              // at every iteration I have to re-multiplicate the previous division
        runningSquaredPut = ((runningSquaredPut * i) + pow(put, 2)) / (i + 1);
        errorPut = Error(runningSumPut, runningSquaredPut, i);
        putFile << runningSumPut << " " << errorPut << endl;

    }

    callFile.close();
    putFile.close();

    // 3.1 - Discrete sampling
    callFile.open("3.1.call_discrete.txt");
    putFile.open("3.1.put_discrete.txt");

    // Redefine variables adding the `step` information 
    int steps = 100;
    double increment = T / double (steps);
    runningSumCall = 0.;
    runningSquaredCall = 0.;
    runningSumPut = 0.;
    runningSquaredPut = 0.;
    errorCall = 0.;
    errorPut = 0.;
    double S_i = 0.;

    for (int i = 0; i < N; i++){                                                               // cycle on N blocks

        call = 0;
        put = 0;

        for (int j = 0; j < L; j++){                                                           // cycle in L throws within the block

            S_i = S_0;

            for (int h = 0; h < steps; h ++){                                                  // cycle on steps       

                x = rnd->Gauss(0, 1);                                                          // compute the single step of the pricing process
                S_i = S_i * exp((r - sigma * sigma * 0.5) * increment + (sigma * x * sqrt(increment)));     

            }

            call += std::max(0., S_i - K) * exp(-r * T);                                       // Compute call price accordingly
            put += std::max(0., K - S_i) * exp(-r * T);                                        // Compute put price accordingly
        }

        call /= L;
        put /= L;

        // Call - data blocking
        runningSumCall = ((runningSumCall * i) + call) / (i + 1);                               // at every iteration I have to re-multiplicate the previous division
        runningSquaredCall = ((runningSquaredCall * i) + pow(call, 2)) / (i + 1);
        errorCall = Error(runningSumCall, runningSquaredCall, i);
        callFile << runningSumCall << " " << errorCall << endl;

        // Put - data blocking
        runningSumPut = ((runningSumPut * i) + put) / (i + 1);                                  // at every iteration I have to re-multiplicate the previous division
        runningSquaredPut = ((runningSquaredPut * i) + pow(put, 2)) / (i + 1);
        errorPut = Error(runningSumPut, runningSquaredPut, i);
        putFile << runningSumPut << " " << errorPut << endl;

    }

    callFile.close();
    putFile.close();

    rnd->SaveSeed();
    return 0;
}    