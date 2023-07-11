#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "lib.h"
#include "random.h"
#include "functions.h"
#include "annealing.h"
#include "integral.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    // Initialize an element of the Random class
    Random *rnd = new Random();
    rnd->RandomRoutine();                                                           // Random settings including setting the seed, ....

    // M thows, N blocks, L throws in each block
    int M = pow(10,6);    
    int N= 1000;            
    int L = int(M / N);

    cout << "Working with " << L << " throws in each block" << endl;

    //temperature
    double T = 2.;

    // Build pointers to the objects I need in the SA schema
    TrialWaveFunction *TWF = new TrialWaveFunction(rnd);                            // my trial wave function
    DoubleDwellPotential *DDP = new DoubleDwellPotential();                         // my potential term of the Hamiltonian
    Annealing *ANN = new Annealing(TWF, rnd, DDP);                                  // used for the SA scheme

    //Set the cycle variables:
    ANN->SetN(N);
    ANN->SetBeta(1 / T);
    ANN->SetL(L);

    // File managment
    ofstream Averages;
    ofstream Coordinates;
    ofstream Parameters;

    Parameters.open("Parameters.txt");
    Averages.open("averages.txt");

    Averages << T << " " << ANN->GetEnergy() << " " << ANN->GetError() << " " << ANN-> GetEnergy()<< endl;
    Parameters << T << " " << ANN->GetMu() << " " << ANN->GetSigma() << endl;

    // Auxiliary variables that I will use to output the optimized set of parameters
    double bestEnergy = ANN->GetEnergy();
    arma::vec bestParameters { ANN->GetMu(), ANN->GetSigma()};

    while (T >= 0.01) {                                                             // Loop over the temperature decreasing (annealing)

        cout << "T: " << T << endl;

        //Updates TWF parameters and computes the new energy, if it is accepted by metropolis, it saves the new configuration
        ANN->MetropolisAnnealing();                                                 // it includes equilibration for the TWF with the new set of parameters
        Averages << T << " " << ANN->GetEnergy() << " " << ANN->GetError() << " " << bestEnergy << endl;
        Parameters << T << " " << ANN->GetMu() << " " << ANN->GetSigma() << endl;
        
        // Update temperature 
        T = T * 0.997;
        ANN->SetBeta(1/T);
        if (bestEnergy > ANN->GetEnergy()) {                                        // If a better configuration is found, save it
            bestEnergy = ANN->GetEnergy();
            bestParameters = { ANN->GetMu(),ANN->GetSigma() };
        }
                    
    }

    cout << "Lowest energy found: " << bestEnergy << endl;
    cout << "Optimized parameters: mu = " << bestParameters[0] << " sigma = " << bestParameters[1] << endl;

    Averages.close();
    Parameters.close();

    // The SA scheme is concluded, now that I identified an optimized set of parameters, i compute the corresponding energy with data blocking and a wider statistics
    Averages.open("OptimizedEnergy.txt");
    Coordinates.open("OptimizedCoordinates.txt");

    cout << "Finding the optimized energy" << endl;

    TWF->Set_Mu(bestParameters[0]);                                                 // Setting the optimized parameters
    TWF->Set_Sigma(bestParameters[1]);

    double runningSum = 0.;
    double runningSquared = 0.;
    double error = 0.;
    double integral = 0.;

    TWF->Equilibrate(100, pow(10, 3));                                              // Equilibrate the system before Metropolis

    for (int i = 0; i < N; i++) {                                                   // Loop over the blocks

        integral = 0;

        for (int j = 0; j < L; j++) {                                               // Loop over the throws

            TWF->MetropolisUniform();                                               // Evolve the position through metropolis
            // Now we need to compute the energy
            integral += ((-0.5 * TWF->SecondDerivative(TWF->Get_Position())) / TWF->EvaluateNoModulus(TWF->Get_Position())) + DDP->Evaluate(TWF->Get_Position());
            Coordinates << TWF->Get_Position() << endl;
        }

        // Data blocking - I/O managment
        integral /= L;
        runningSum = ((runningSum * i) + integral) / (i + 1);                       // at every iteration I have to re-multiplicate the previous division
        runningSquared = ((runningSquared * i) + pow(integral, 2)) / (i + 1);
        error = Error(runningSum, runningSquared, i);
        Averages << i << " " << runningSum << " " << error << endl;
    }

    Averages.close();
    Coordinates.close();

    rnd->SaveSeed();
    return 0;
}