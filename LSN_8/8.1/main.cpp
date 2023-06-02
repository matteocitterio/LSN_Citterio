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

    //variables for acceptance rate of the algorithm
    double accepted = 0.;
    double attempted = 0.;

    //Define global variables
    double a0 = 0.0529 * pow(10, -9); // reduced units
    double c = 2 ;       //to tune the acceptance rate
    double integral = 0;             

    //Blocking average variables
    double runningSum = 0.;
    double runningSquared = 0.;
    double error = 0.;

    //File managment
    ofstream Averages;
    ofstream Coordinates;

    double x = 0;       //my position starting point
    TrialWaveFunction *TWF = new TrialWaveFunction();   //my trial wave function
    DoubleDwellPotential *DDP = new DoubleDwellPotential();  //my potential term of the Hamiltonian

    //Open output files
    Averages.open("Averages0.txt");
    Coordinates.open("Coordinates0.txt");

    //Set the parameters
    TWF->Set_Mu(1);
    TWF->Set_Sigma(0.5);

    //equilibrate the system
    cout << "----------------------------" << endl;
    cout <<" Doing equilibration " << endl;
    cout << "----------------------------" << endl;
    EquilibrateUN(100, pow(10,3), x, TWF, rnd, c);
    cout << "----------------------------" << endl;

    for (int i = 0; i < N; i++){
        integral = 0;
        accepted = 0;
        attempted = 0;
        for (int j = 0; j < L; j++){           

                //updates my position according to my trial wave function distribution
                MetropolisUniform(x, TWF, rnd, c, accepted, attempted);
                //Now we need to compute the energy
                integral += (( -0.5 * TWF->SecondDerivative(x) ) / TWF->EvaluateNoModulus(x))  + DDP->Evaluate(x);  //Energy expectation
                Coordinates << x << endl;
        }
        integral /= L; 

        runningSum = ((runningSum * i) + integral) / (i + 1);                  // at every iteration I have to re-multiplicate the previous division
        runningSquared = ((runningSquared * i) + pow(integral, 2)) / (i + 1);
        error = Error(runningSum, runningSquared, i);
        Averages << i <<" "<< runningSum << " " << error << endl;

        if (i%10 ==0){
                cout <<"Block: "<< i << " Uniform Acceptance rate " << accepted / attempted << endl<< endl;
        }
        
    }

    Averages.close();
    Coordinates.close();
    rnd->SaveSeed();
    return 0;
}    