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

    // Initialize an element of the Random class
    Random *rnd = new Random();
    rnd->RandomRoutine();                                                   // Random settings including setting the seed, ....

    // M thows, N blocks, L throws in each block
    int M = pow(10,6);      
    int N= 1000;
    int L = int(M / N);

    cout << "Working with " << L << " throws in each block" << endl;

    // Variables for acceptance rate of the algorithm
    double accepted = 0.;
    double attempted = 0.;

    // Define global variables
    double a0 = 0.0529 * pow(10, -9);                                       // Reduced units
    double c = 1.2 * a0;                                                    // To tune the acceptance rate
    double r = 0;

    //Blocking average variables
    double runningSum = 0.;
    double runningSquared = 0.;
    double error = 0.;

    //File managment
    ofstream Averages;
    ofstream Coordinates;

    // 5.1 - Ground State uniform transition probability

    arma::vec v(3);                                                          // Vector of positions
    
    for (int i = 0; i < 3; i++) {v[i] = a0;}                                 // Initialize the position with initial starting point the bohr radius
    HydrogenGS *HGS = new HydrogenGS();                                      // Initialize a Ground State function object

    Averages.open("1.0.0.averages.Uniform.txt");                             // Open output files
    //Coordinates.open("100UniformCoo.txt");                                 // Uncomment if you want to print out the coordinates

    // Equilibrate the system
    cout << "----------------------------" << endl;
    cout <<" Doing equilibration " << endl;
    cout << "----------------------------" << endl;
    EquilibrateUN(20, pow(10,3), v, HGS, rnd, c);
    cout << "----------------------------" << endl;


    for (int i = 0; i < N; i++){                                              // Loop over the blocks
        r = 0;                                                                // re-initialize position
        accepted = 0;                                                         // re-initialize counters
        attempted = 0;

        for (int j = 0; j < L; j++){                                          // Loop over the throws in each block

                MetropolisUniform(v, HGS, rnd, c, accepted, attempted);       // Metropolis with uniform transition probability
                r += arma::norm(v,2);                                         // Update position
                //Coordinates << v[0] << " " << v[1] << " " << v[2] << endl;  // Uncomment if you want to print out the coordinates
        }

        // Data - blocking and I/O managment
        r /= L;

        runningSum = ((runningSum * i) + r) / (i + 1);                        // at every iteration I have to re-multiplicate the previous division
        runningSquared = ((runningSquared * i) + pow(r, 2)) / (i + 1);
        error = Error(runningSum, runningSquared, i);
        Averages << i <<" "<< runningSum << " " << error << endl;

        if (i%10 ==0){
                cout <<"Block: "<< i << " 100 Uniform Acceptance rate " << accepted / attempted << endl<< endl;
        }
        
    }

    Averages.close();
    //Coordinates.close();                                                    // Uncomment if you want to print out the coordinates

    // 5.1 - Ground State Gaussian Probability transition

    // Reinitialize variables
    for (int i = 0; i < 3; i++){ v[i] = a0;}
    runningSum = 0.;
    runningSquared = 0.;
    error = 0.;

    //Open output files
    Averages.open("1.0.0.averages.Gauss.txt");
    //Coordinates.open("100GaussCoo.txt");

    // Equilibrate the system
    cout << "----------------------------" << endl;
    cout << " Doing equilibration " << endl;
    cout << "----------------------------" << endl;
    EquilibrateGA(20, pow(10, 3), v, HGS, rnd, c);
    cout << "----------------------------" << endl;

    for (int i = 0; i < N; i++) {                                             // Loop over the blocks

        r = 0;                                                                // Re-initialiaze positions and counters
        accepted = 0;
        attempted = 0;
        for (int j = 0; j < L; j++) {                                         // Loop over the throws in each block

                MetropolisGauss(v, HGS, rnd, c, accepted, attempted);         // Metropolis with Gaussian transition
                r += arma::norm(v, 2);
                //Coordinates << v[0] << " " << v[1] << " " << v[2] << endl;
        }

        // Data - blocking
        r /= L;

        runningSum = ((runningSum * i) + r) / (i + 1);                        // At every iteration I have to re-multiplicate the previous division
        runningSquared = ((runningSquared * i) + pow(r, 2)) / (i + 1);
        error = Error(runningSum, runningSquared, i);
        Averages << i << " " << runningSum << " " << error << endl;

        if (i % 10 == 0)
        {
                cout << "Block: " << i << " 100 Gauss Acceptance rate " << accepted / attempted << endl
                     << endl;
        }
    }

    Averages.close();
    //Coordinates.close();

    // 5.1 -  Excited state uniform transition

    c = 2.9 * a0;                                                             // Change step size in order to get 50% acceptance rate

    Hydrogen210 *H210 = new Hydrogen210();                                    // Excited level

    // Set variables to zero:
    for (int i = 0; i < 3; i++){v[i] = a0;}
    runningSum = 0.;
    runningSquared = 0.;
    error = 0.;

    // Open output files
    Averages.open("2.1.0.averages.Uniform.txt");
    //Coordinates.open("210UniformCoo.txt");

    // equilibrate the system
    cout << "----------------------------" << endl;
    cout << " Doing equilibration " << endl;
    cout << "----------------------------" << endl;
    EquilibrateUN(20, pow(10, 3), v, H210, rnd, c);
    cout << "----------------------------" << endl;

    for (int i = 0; i < N; i++) {                                             // Loop over the blocks

        r = 0;                                                                // Reinitialize positions and counters
        accepted = 0;
        attempted = 0;

        for (int j = 0; j < L; j++) {                                         // Loop over throws in each block

                MetropolisUniform(v, H210, rnd, c, accepted, attempted);      // Metropolis with uniform transition
                r += arma::norm(v, 2);                                        // Update position
                //Coordinates << v[0] << " " << v[1] << " " << v[2] << endl;
        }

        // Data - blocking
        r /= L;

        runningSum = ((runningSum * i) + r) / (i + 1);                         // At every iteration I have to re-multiplicate the previous division
        runningSquared = ((runningSquared * i) + pow(r, 2)) / (i + 1);
        error = Error(runningSum, runningSquared, i);
        Averages << i << " " << runningSum << " " << error << endl;

        if (i % 10 == 0)
        {
                cout << "Block: " << i << " 210 Uniform Acceptance rate " << accepted / attempted << endl
                     << endl;
        }
    }

    Averages.close();
    //Coordinates.close();

    // 5.1 - Excited state Gaussian transition

    // Set variables to zero:
    for (int i = 0; i < 3; i++){v[i] = a0;}
    runningSum = 0.;
    runningSquared = 0.;
    error = 0.;

    // Open output files
    Averages.open("2.1.0.averages.Gauss.txt");
    //Coordinates.open("210GaussCoo.txt");

    // equilibrate the system
    cout << "----------------------------" << endl;
    cout << " Doing equilibration " << endl;
    cout << "----------------------------" << endl;
    EquilibrateGA(20, pow(10, 3), v, H210, rnd, c);
    cout << "----------------------------" << endl;

    for (int i = 0; i < N; i++) {                                             // Loop over the blocks

        r = 0;                                                                // Reinitialize positions and counters
        accepted = 0;
        attempted = 0;

        for (int j = 0; j < L; j++) {                                         // Loop over throws in each block

                MetropolisGauss(v, H210, rnd, c, accepted, attempted);        // Metropolis with gaussian transition
                r += arma::norm(v, 2);
                //Coordinates << v[0] << " " << v[1] << " " << v[2] << endl;
        }

        // Data - blocking
        r /= L;

        runningSum = ((runningSum * i) + r) / (i + 1);                        // At every iteration I have to re-multiplicate the previous division
        runningSquared = ((runningSquared * i) + pow(r, 2)) / (i + 1);
        error = Error(runningSum, runningSquared, i);
        Averages << i << " " << runningSum << " " << error << endl;

        if (i % 10 == 0)
        {
                cout << "Block: " << i << " 210 Gauss Acceptance rate " << accepted / attempted << endl
                     << endl;
        }
    }

    Averages.close();
    //Coordinates.close();

    rnd->SaveSeed();
    return 0;
}    