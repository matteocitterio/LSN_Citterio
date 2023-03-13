#include <iostream>
#include <fstream>
#include <string>
#include "lib.h"
#include "random.h"
#include "functions.h"
#include "integral.h"


using namespace std;

int main(int argc, char *argv[]){

    //Settings for the Random generator class

    Random* rnd = new Random();
    rnd -> RandomRoutine();   //Random settings including seed etc

    // M thows, N blocks, L throws in each block
    int M = pow(10,4);
    int N= 300;
    int L = int(M / N);

    cout << "Working with " << L << " throws in each block" << endl;

    // 2.1.1

    //Build up pointers for Cosine class and Integral class
    Cosine *cosine = new Cosine(M_PI * 0.5, M_PI*0.5, 0.);
    Integral *integral = new Integral(0., 1., cosine, rnd);

    //open output stream
    ofstream outputFile;
    outputFile.open("2.1.1.txt");

    double x = 0.;
    double runningSum = 0.;
    double runningSquared = 0.;
    double error = 0.;

    for (int i = 0; i < N; i++){
        x = integral -> arithAverage(L);
        runningSum = ((runningSum * i) + x) / (i +1);       //at every iteration I have to re-multiplicate the previous division
        runningSquared = ((runningSquared * i) + pow(x,2))/ (i+1);
        error = Error(runningSum, runningSquared, i);
        outputFile << runningSum << " " << error << endl;
    }

    outputFile.close();

    //2.1.2

    // Build up pointer for Parabola class
    Parabola *parabola = new Parabola(-1.5,0,1.5);

    // open output stream
    outputFile.open("2.1.2.txt");

    x = 0.;
    runningSum = 0.;
    runningSquared = 0.;
    error = 0.;

    for (int i = 0; i < N; i++){
        x = integral->Imp_sampling(L, 1.5, parabola);
        runningSum = ((runningSum * i) + x) / (i + 1); // at every iteration I have to re-multiplicate the previous division
        runningSquared = ((runningSquared * i) + pow(x, 2)) / (i + 1);
        error = Error(runningSum, runningSquared, i);
        outputFile << runningSum << " " << error << endl;
    }

    outputFile.close();

    rnd->SaveSeed();
    return 0;
}    