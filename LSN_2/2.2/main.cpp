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
    int M = pow(10,4);      
    int N= 100;
    int Steps = 100;
    int L = int(M / N);

    cout << "Working with " << L << " throws in each block" << endl;

    // 2.2.1

    //cubic lattice random walk
    vec runningSum (Steps);
    vec runningSquared(Steps);
    vec error(Steps);
    vec position(Steps);
    vec r(3);

    for (int i = 0; i < N; i++){                //cycle on blocks

        position.zeros();

        for (int j = 0; j < L; j++){            //cycle on throws in each block
            r.zeros();
            for (int k = 0; k < Steps; k ++){   //cycle on RW steps: adding step wise an entire random walk
                position(k) += (pow(r(0),2) + pow(r(1),2) + pow(r(2),2)); 
                rnd -> Step(r);                 //passed by ref
            }
        }

        for (int k = 0; k < Steps; k++){
            runningSum(k) += sqrt(position(k)/L);  //sqrt(<|r|^2>)
            runningSquared(k)+= (position(k)/L);   //<|r|^2>
        }
    }

    runningSum/=N;
    runningSquared/=N;

    // open output stream
    ofstream outputFile;
    outputFile.open("2.2.1.txt");

    for (int k = 0; k < Steps; k++){
        error(k) = Error(runningSum(k), runningSquared(k), N);
        outputFile << runningSum(k) << " " << error(k) << endl;
    }

    outputFile.close();

    // 2.2.2

    // cubic lattice random walk
    runningSum.zeros();
    runningSquared.zeros();
    error.zeros();
    position.zeros();

    for (int i = 0; i < N; i++){ // cycle on blocks

        position.zeros();

        for (int j = 0; j < L; j++){ // cycle on throws in each block
            r.zeros();
            for (int k = 0; k < Steps; k++){ // cycle on RW steps: adding step wise an entire random walk
                position(k) += (pow(r(0), 2) + pow(r(1), 2) + pow(r(2), 2));
                rnd->CStep(r); // passed by ref
            }
        }

        for (int k = 0; k < Steps; k++)
        {
            runningSum(k) += sqrt(position(k) / L); // sqrt(<|r|^2>)
            runningSquared(k) += (position(k) / L); //<|r|^2>
        }
    }

    runningSum /= N;
    runningSquared /= N;

    // open output stream
    outputFile.open("2.2.2.txt");

    for (int k = 0; k < Steps; k++)
    {
        error(k) = Error(runningSum(k), runningSquared(k), N);
        outputFile << runningSum(k) << " " << error(k) << endl;
    }

    outputFile.close();

    rnd->SaveSeed();
    return 0;
}    