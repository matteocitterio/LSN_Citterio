#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "random.h"
#include "geneticalgorithm.h"
#include "mpi.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    //MPI initialization
    int size, rank;
    MPI_Init(&argc, &argv);                    // Initialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);      // Get size
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);      // Get current rank
    MPI_Status status;                         // Communications status

    double Tstart = MPI_Wtime();

    //Settings for the Random generator class
    Random* rnd = new Random();
    rnd -> RandomRoutineParallel(rank);                     // Random settings including seed etc

    //Define all the variables and parameters
    ifstream ReadInput;
    int M, NumberOfGenerations, CircleSquare, NMigr;
    double probSwapMutation, probShiftMutation, probPermutationMutation, probInversionMutation,probCrossover, p;

    //Read input parameters for the algorithm
    ReadInput.open("input.in");
    ReadInput >> CircleSquare;                  // Put the cities on a circle or on a square
    ReadInput >> M;                             // Number of chromosomes, i.e. number of initial solutions
    ReadInput >> NumberOfGenerations;           // Number of generations
    ReadInput >> probSwapMutation;              // probability of performing a SWAP mutation
    ReadInput >> probShiftMutation;             // probability of performing a SHIFT mutation
    ReadInput >> probPermutationMutation;       // probability of performing a PERMUTATION mutation
    ReadInput >> probInversionMutation;         // probability of performing a INVERSION mutation
    ReadInput >> probCrossover;                 // probability of performing a CROSSOVER
    ReadInput >> p;                             // Exponenent of the loaded die
    ReadInput >> NMigr;                         // Number of migrations per continent in parallel computing

    if (rank == 0){

        cout << "\nProbability of CROSSOVER: " << probCrossover * 100 << "%" << endl;
        cout << "Probability of SWAP mutation: " << probSwapMutation * 100 << "%" << endl;
        cout << "Probability of SHIFT mutation: " << probShiftMutation * 100 << "%" << endl;
        cout << "Probability of PERMUTATION mutation: " << probPermutationMutation * 100 << "%" << endl;
        cout << "Probability of INVERSION mutation: " << probInversionMutation * 100 << "%" << endl;

        cout << "\nNumber of chromosomes: " << M << endl;
        cout << "Number of Generations: " << NumberOfGenerations << endl;
        cout << "Loaded die exponent: " << p << endl;
        cout << "Number of migration: " << NMigr << endl;
    }

    // GA constructor
    GeneticAlgorithm *GA = new GeneticAlgorithm(rnd,CircleSquare, NumberOfGenerations, M, probSwapMutation, probShiftMutation, probPermutationMutation, probInversionMutation, probCrossover, p);

    cout << "\nRANK: "<<rank<<" Initializing population and cities..";
    GA->Start();
    cout << "RANK: " << rank << " DONE initialize" << endl;

    MPI_Barrier(MPI_COMM_WORLD);                // Wait for all the processes to initialize

    if (size == 1){
        cout << "Serial" << endl;
        GA->Evolve();
        cout << "Done evolving" << endl;
    }
    else{
        cout << "RANK: " << rank << " Evolve" << endl;
        GA->ParallelEvolve(NMigr, MPI_COMM_WORLD, rank, size, status);
        cout << "RANK: " << rank << " DONE evolving" << endl;
    }

    rnd->SaveSeed();

    double Tout = MPI_Wtime();
    double DeltaT = Tout-Tstart;
    cout << "Time: " << DeltaT;

    MPI_Finalize();                             // MPI finilization
    return 0;
}