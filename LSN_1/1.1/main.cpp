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

    // M thows, N blocks, L throws in each block

    int M = 1000000;
    int N= 100;
    int L = int(M / N);

    ofstream outputFile;

    //1.1.1

    vec cum_av1, cum_err1;
    string flag = "averages";

    // tie creates a tuple of references and assigns the corresponding values from the tuple on the right-hand side.
    tie(cum_av1, cum_err1) = ProgressiveStatistics(flag , rnd, N, L);

    //Write output file

    outputFile.open("1.1.1.txt");

    for (int i = 0; i < N; i++){
        outputFile << cum_av1[i] << " "<< cum_err1[i] << endl;
    }

    // Close output file
    outputFile.close();

    //1.1.2

    vec cum_av2, cum_err2;
    flag = "uncertainty";
    tie(cum_av2, cum_err2) = ProgressiveStatistics(flag , rnd, N, L);

    //Write output file
    outputFile.open("1.1.2.txt");

    for (int i = 0; i < N; i++){

        outputFile << cum_av2[i] << " " << cum_err2[i] << endl;
    }

    //Close output file
    outputFile.close();

    //1.1.3

    int nbins = 100; //number of bins
    int n = 10000; //throws in each bin
    double expectation = double(n/nbins);
    vec bins = linspace<vec>(0,1,nbins+1); //vector of bins edges
    int repetitions = 100; //number of repetitions of the procedure

    outputFile.open("1.1.3.txt");

    for (int h=0; h<repetitions; h++){

        vec histogram(nbins, fill:: zeros);

        for (int j = 0; j<n; j ++){
            double x = rnd.Rannyu();

            if (x > 0.99){
                histogram(99)+=1;
            }

            else{
                for (int k =0; k<99; k++){
                    if (bins(k) <= x && x < bins(k+1)){
                        histogram(k)+=1;
                    }
                }
            }
        }

        double chi2=0;
        for (int i = 0; i < nbins; i++){
            chi2 += (pow((histogram(i) - expectation), 2) / expectation);
        }

        outputFile<< chi2 << endl;
    }

    outputFile.close();

    //1.1.3_extra
    //I wanted to make sure that my chi2 is compatible with a chi2 pdf with 99 degrees of freedom

    outputFile.open("1.1.3_extra.txt");

    for (int h = 0; h < 100*repetitions; h++)
    {

        vec histogram(nbins, fill::zeros);

        for (int j = 0; j < n; j++)
        {
            double x = rnd.Rannyu();

            if (x > 0.99)
            {
                histogram(99) += 1;
            }

            else
            {
                for (int k = 0; k < 99; k++)
                {
                    if (bins(k) <= x && x < bins(k + 1))
                    {
                        histogram(k) += 1;
                    }
                }
            }
        }

        double chi2 = 0;
        for (int i = 0; i < nbins; i++)
        {
            chi2 += (pow((histogram(i) - expectation), 2) / expectation);
        }

        outputFile << chi2 << endl;
    }

    outputFile.close();

    rnd.SaveSeed();
    return 0;
}