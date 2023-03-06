#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "random.h"
#include "funzioni.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    //Settings for the Random generator class

    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open())
    {
        Primes >> p1 >> p2;
    }
    else
        cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while (!input.eof())
        {
            input >> property;
            if (property == "RANDOMSEED")
            {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    }
    else
        cerr << "PROBLEM: Unable to open seed.in" << endl;

    // M thows, N blocks, L throws in each block
    int M = 10000000;
    int N= 200;
    int L = int(M / N);

    cout << L << endl;

    //d spacing between lines, l length of the needle

    double d = 1.5;
    double l = 1.0;

    ofstream outputFile;

    vec cum_av, cum_err;
    string flag = "buffon";

    // tie creates a tuple of references and assigns the corresponding values from the tuple on the right-hand side.
    tie(cum_av, cum_err) = ProgressiveStatistics(flag, rnd, N, L, l, d);

    // Write output file

    outputFile.open("buffon.txt");

    for (int i = 0; i < N; i++){
        outputFile << cum_av[i] << " " << cum_err[i] << endl;
    }

    // Close output file
    outputFile.close();

    rnd.SaveSeed();
    return 0;
}    