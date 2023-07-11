#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "random.h"
#include "geneticalgorithm.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    //Settings for the Random generator class
    Random* rnd = new Random();
    rnd -> RandomRoutine();                     // Random settings including seed etc

    //Define all the variables and parameters
    ifstream ReadInput;
    int M, NumberOfGenerations, CircleSquare;
    double probSwapMutation, probShiftMutation, probPermutationMutation, probInversionMutation,probCrossover, p;

    //Read input parameters for the algorithm
    ReadInput.open("input.in");
    ReadInput >> CircleSquare;                                                      // Put the cities on a circle or on a square
    ReadInput >> M;                                                                 // Number of chromosomes, i.e. number of initial solutions
    ReadInput >> NumberOfGenerations;                                               // Number of generations
    ReadInput >> probSwapMutation;                                                  // probability of performing a SWAP mutation
    ReadInput >> probShiftMutation;                                                 // probability of performing a SHIFT mutation
    ReadInput >> probPermutationMutation;                                           // probability of performing a PERMUTATION mutation
    ReadInput >> probInversionMutation;                                             // probability of performing a INVERSION mutation
    ReadInput >> probCrossover;                                                     // probability of performing a CROSSOVER
    ReadInput >> p;                                                                 // Exponenent of the loaded die

    // User friendly messages
    cout << "\nProbability of CROSSOVER: " << probCrossover * 100 << "%" << endl;
    cout << "Probability of SWAP mutation: " << probSwapMutation * 100 << "%" << endl;
    cout << "Probability of SHIFT mutation: " << probShiftMutation * 100 << "%" << endl;
    cout << "Probability of PERMUTATION mutation: " << probPermutationMutation * 100 << "%" << endl;
    cout << "Probability of INVERSION mutation: " << probInversionMutation * 100 << "%" << endl;

    cout << "\nNumber of chromosomes: " <<M<<endl;
    cout << "Number of Generations: " << NumberOfGenerations << endl;
    cout << "Loaded die exponent: " << p << endl;

    // GA constructor
    GeneticAlgorithm *GA = new GeneticAlgorithm(rnd,CircleSquare, NumberOfGenerations, M, probSwapMutation, probShiftMutation, probPermutationMutation, probInversionMutation, probCrossover, p);

    // GA initialization
    cout << "\n\nInitializing population and cities..";
    GA->Start();
    cout << " DONE" <<endl;
    
    // GA evolution
    cout << "Evolve"<<endl;
    GA->Evolve();
    cout << "DONE"<<endl;

    rnd->SaveSeed();
    return 0;
}