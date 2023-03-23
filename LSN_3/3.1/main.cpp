#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "lib.h"
#include "random.h"
#include "functions.h"
#include "integral.h"


using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    //Settings for the Random generator class
    Random* rnd = new Random();
    rnd -> RandomRoutine();   //Random settings including seed etc

    // M thows, N blocks, L throws in each block
    int M = pow(10,6);      
    int N= 1000;
    int L = int(M / N);

    cout << "Working with " << L << " throws in each block" << endl;

    //Define variables of the pricing process
    double r = 0.1;
    double T = 1.;
    double sigma = 0.25;
    double K = 100;
    double S_0 = 100;

    // 3.1 - Direct sampling

    ofstream callFile;
    ofstream putFile;
    callFile.open("3.1.call.txt");
    putFile.open("3.1.put.txt");

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


    for (int i = 0; i < N; i++){
        call = 0;
        put = 0;
        for (int j = 0; j < L; j++){
            x= rnd->Gauss(0, T);
            S_T = S_0 * exp( (r - sigma * sigma * 0.5 ) * T + (sigma * x ) );
            call += std::max(0., S_T - K) * exp (- r * T);
            put +=std::max(0., K - S_T)* exp (- r * T);
        }
        call /= L;
        put /= L;
        //call
        runningSumCall = ((runningSumCall * i) + call) / (i + 1); // at every iteration I have to re-multiplicate the previous division
        runningSquaredCall = ((runningSquaredCall * i) + pow(call, 2)) / (i + 1);
        errorCall = Error(runningSumCall, runningSquaredCall, i);
        callFile << runningSumCall << " " << errorCall << endl;
        //put
        runningSumPut = ((runningSumPut * i) + put) / (i + 1); // at every iteration I have to re-multiplicate the previous division
        runningSquaredPut = ((runningSquaredPut * i) + pow(put, 2)) / (i + 1);
        errorPut = Error(runningSumPut, runningSquaredPut, i);
        putFile << runningSumPut << " " << errorPut << endl;
    }

    callFile.close();
    putFile.close();

    // 3.1 - Discrete sampling
    callFile.open("3.1.call_discrete.txt");
    putFile.open("3.1.put_discrete.txt");

    //redefine variables
    int steps = 100;
    double increment = T / double (steps);
    runningSumCall = 0.;
    runningSquaredCall = 0.;
    runningSumPut = 0.;
    runningSquaredPut = 0.;
    errorCall = 0.;
    errorPut = 0.;
    double S_i = 0.;

    for (int i = 0; i < N; i++){                            //cycle on N blocks
        call = 0;
        put = 0;
        for (int j = 0; j < L; j++){                        //cycle in L throws within the block
            S_i = S_0;
            for (int h = 0; h < steps; h ++){               //cycle on steps       
                x = rnd->Gauss(0, 1);
                S_i = S_i * exp((r - sigma * sigma * 0.5) * increment + (sigma * x * sqrt(increment)));
            }
            call += std::max(0., S_i - K) * exp(-r * T);
            put += std::max(0., K - S_i) * exp(-r * T);
        }
        call /= L;
        put /= L;
        // call
        runningSumCall = ((runningSumCall * i) + call) / (i + 1); // at every iteration I have to re-multiplicate the previous division
        runningSquaredCall = ((runningSquaredCall * i) + pow(call, 2)) / (i + 1);
        errorCall = Error(runningSumCall, runningSquaredCall, i);
        callFile << runningSumCall << " " << errorCall << endl;
        // put
        runningSumPut = ((runningSumPut * i) + put) / (i + 1); // at every iteration I have to re-multiplicate the previous division
        runningSquaredPut = ((runningSquaredPut * i) + pow(put, 2)) / (i + 1);
        errorPut = Error(runningSumPut, runningSquaredPut, i);
        putFile << runningSumPut << " " << errorPut << endl;
    }

    callFile.close();
    putFile.close();

    rnd->SaveSeed();
    return 0;
}    