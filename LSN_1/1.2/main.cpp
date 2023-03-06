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
    if (input.is_open())
    {
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

    int n = 10000; //number of realizations

    vec N = {1,2,10,100}; //Number of rvs to be added together

    //every file will have 10e4 numbers for each N(i) 
    ofstream standardFile;
    ofstream exponentialFile;
    ofstream lorentianFile;

    standardFile.open("standardDice.txt");
    exponentialFile.open("exponentialDice.txt");
    lorentianFile.open("lorentianDice.txt");

        for (int i = 0; i < int(N.size()); i++){

            for (int j=0; j<n; j++){

            double S_standard = 0;
            double S_exponential = 0;
            double S_lorentian = 0;

            // sum up N(i) rvs
            for (int h = 0; h < N(i); h++){
                S_standard += rnd.Rannyu();
                S_exponential += rnd.Exp(1);
                S_lorentian += rnd.Lorentz(0., 1.);
                }

                S_standard /= double(N(i));
                S_exponential /= double(N(i));
                S_lorentian /= double(N(i));

                // print on file
                standardFile << S_standard << " ";
                exponentialFile << S_exponential << " ";
                lorentianFile << S_lorentian << " ";
            }

        standardFile << endl;
        exponentialFile << endl;
        lorentianFile << endl;
    }
    rnd.SaveSeed();
    return 0;
}