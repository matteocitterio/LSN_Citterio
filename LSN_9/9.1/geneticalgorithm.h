#ifndef __Genetic__
#define __Genetic__

#include <cmath>
#include <limits>
#include "random.h"
#include <armadillo>

using namespace std;
using namespace arma;

class Chromosome
{

public:
    Chromosome(Random *rnd);                                        // constructor
    ~Chromosome() { ; };                                            // empty destructor

    double Fitness(arma::vec chromosome, arma::mat cities);         // returns the fitness of a given chromosome

    int GetNumberOfGenes() {return m_N+1;};
    void PrintChromosome(arma::vec chromosome);

    arma::vec NewChromosome();
    void CheckChromosome(arma::vec chromosome);                     // checks if the chromosome is an admissabile solution (check the cities are visited once, check the last is the first one)

private:
    Random *m_rnd;
    int m_N;                                                        // Number of genes, i.e. number of cities
    int m_initialcity;

};

class Generation  {

    /*
    DATA STRUCTURES:
        - Chromosomes: a possible solution (e.g. [1,3,4,2,5,1] or program elements or... any data structure)
        - Genes: elements of the chromosome (e.g. the numbers representing our city in the chromosome vector)
        - Alleses: Possible values for the genes, in the city problem a number \in [1,MAX_CITY_NUMBER]

    OPERATORS AND DEFINITIONS:
        - Fitness of a chromosome: measure of goodness of the solution (e.g. L(crhromosome)).
        - Selection: selects chromosomes in the population for reproduction. The fitter the chromosome the more likely it is to be selected
        - Crossover: Once two individuals have been selected, both parents pass their chromosomes onto their offspring.
        - Mutation: conversion of genes from one to another. Helps preventing the population of stagnating (\sim to \beta in a simulation annealing)
        (P.N: crossover and mutation destroy previous solutions)
    */

    public:
        Generation(Random *rnd, int M, double probSwapMutation, double probShiftMutation, double probPermutationMutation, double probInversionMutation, double probCrossover, double p); 
        ~Generation() { delete m_chromo; };                                     // empty destructor

        void CreateInitialPopulation();
        void UpdateFitness();                                                   // Update the fitness of a Generation
        double GetOptimumLoss() {return m_fitness.at(0, 0); };                  // Returns fitness of the fittest solution so far
        void CitiesOnACircle();                                                 // Generates `NumberOfCities` cities randomly distributed over a circle
        void CitiesOnASquare();                                                 // Generates `NumberOfCities` cities randomly distributed INSIDE a square
        void AmericanCities();                                                  // Reads the american capitals file
        void Sort();                                                            // Sorts the m_gene_pool and m_fitness according to the values in m_fitness
        double bestHalfAverage();                                               // Returns the average fitness value of the best half of the population
        arma::mat BestPath();                                                   // Returns the path of the best solution so far

        //Operators
        int Selection() {return int(m_M * pow(m_rnd->Rannyu(), m_p));};         // Loaded die      
        void CrossOver();                                                       // CrossOver operator
        void Mutations(int index);                                              // Mutations operator, takes the index of the chromosome to mutate

    private:
        
        Random *m_rnd;                                                          // seed
        int m_M;                                                                // number of initial chromosomes
        arma::mat m_gene_pool;                                                  // current set of chromosomes
        arma::mat m_new_gen;                                                    // auxiliary mat to store new generation

        Chromosome *m_chromo;                                                   // used to access all the methods of the `Chromosome` class
        arma::mat m_cities;                                                     // contains coordinates of the cities
        arma::mat m_fitness;                                                    // contains fitness of each chromosome
        double m_p;                                                             // Exponent of the loaded die

        double m_prob_swap_mutation;                                            // probability of performing a SWAP mutation
        double m_prob_shift_mutation;                                           // probability of performing a SHIFT mutation
        double m_prob_permutation_mutation;                                     // probability of performing a PERMUTATION mutation
        double m_prob_inversion_mutation;                                       // probability of performing a INVERSION mutation
        double m_prob_crossover;                                                // Probability of CrossOver
};

class GeneticAlgorithm{
    public:
        GeneticAlgorithm(Random *rnd,int CircleSquare, int NumberOfGenerations, int M, double probSwapMutation, double probShiftMutation, double probPermutationMutation, double probInversionMutation, double probCrossover, double p);
        ~GeneticAlgorithm() { delete m_gen; };                                  // empty destructor

        void Evolve();                                                          // Updates the GA
        void Start();                                                           // Draws city at random and creates the first generation of chromosomes and sorts them

    private : 
        
        Random *m_rnd;                                                          // seed
        int m_NumberOfGenerations;                                              // number of initial chromosomes
        int m_circlesquare;                                                     // Flag for cities distribution
        Generation *m_gen;                                                      // used to access all the methods of the `Chromosome` class

};

#endif // __Genetic__