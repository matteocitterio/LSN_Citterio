#include "geneticalgorithm.h"

#include <chrono>
#include <thread>

using namespace std;
using namespace arma;

Chromosome::Chromosome(Random *rnd){

    m_rnd = rnd;                            //set the random generator
    m_N= 34;                                //Fixed number of cities
    m_initialcity = 1;                          //initial city
    //some routine to fill randomly with genes our chromosome
}

void Chromosome::CheckChromosome(arma::vec chromosome){

    if(chromosome.size() != m_N+1){
        throw std::runtime_error("\n\n\n****ERROR: WRONG CHROMOSOME SIZE****\n\n\n");
    }

    //If the sum of the created chromosome is not equal to the gauss sum +1 (of the last city) throw an error
    if (arma::accu(chromosome) != (m_N * (m_N + 1) / 2) + 1)
    {
        throw std::runtime_error("\n\n\n****ERROR: SUM OF THE CHROMOSOME DOES NOT MATCH EXPECTED VALUE****\n\n\n");
    }

    if (chromosome.back() != chromosome(0))
    {
        throw std::runtime_error("\n\n\n****ERROR: ROUNDTRIP CONSTRAINT VIOLATED*****\n\n\n");
    }
}

void Chromosome::PrintChromosome(arma::vec chromosome){
    for (int i=0; i<chromosome.size(); i++){
        cout << chromosome(i) <<" ";
    }
    cout <<endl;
}

double Chromosome::Fitness(arma::vec chromosome, arma::mat cities){

    double len = 0;

    for (int i=1; i < chromosome.n_elem; i++){

        // cout << "i " << i << " /" << chromosome.n_elem << endl;
        // cout << "cities.n_rows " << cities.n_rows <<endl;
        // cout << "chromosome[i] " << chromosome[i] << " \nchromosome[i-1] " << chromosome[i-1] << endl;
        // cout << "cities\n" << cities << endl;
        // cout << "cities.row(chromosome[i]) " << cities.row(chromosome[i] -1) << "\ncities.row(chromosome(i - 1)) " << cities.row(chromosome(i - 1) -1) << endl;

        len += pow((arma::norm(cities.row(chromosome[i] -1) - cities.row(chromosome(i - 1) -1), 2)), 2);
    }

    return len;

}


//Creates a new random chromosome
arma::vec Chromosome::NewChromosome(){
    arma::vec chromosome = arma::linspace(1,m_N+1, m_N+1);                    //Creates a vector of integers of length: number of cities+1, right now the cities are ordered i.e. [1,2,3,..]
    chromosome.back() = m_initialcity;                              //this sets the last city as the first one: round trip constraint

    //Shuffles a subvector of our chromosome:
    chromosome.subvec(1, chromosome.size() - 2) = arma::shuffle(chromosome.subvec(1, chromosome.size() - 2));

    CheckChromosome(chromosome);

    return chromosome;
}


//=============================================================//

Generation::Generation(Random *rnd, int M, double probSwapMutation, double probShiftMutation, double probPermutationMutation, double probInversionMutation, double probCrossover, double p)
{

    m_rnd = rnd;                                                        // set the random generator
    m_M = M;                                                            // size of the gene pool
    m_p = p;                                                            // exponent of the loaded die used in the selection operator
    m_chromo = new Chromosome(rnd);                                     // m_rnd is already a pointer
    m_gene_pool = arma::mat(m_M, m_chromo->GetNumberOfGenes());         // Create the gene pool
    m_cities = arma::mat(m_chromo->GetNumberOfGenes() -1, 2);           // Contains the coordinates of the cities
    m_fitness = arma::mat(m_M,1);                                       // Contains the fitness of each chromosome

    m_prob_swap_mutation = probSwapMutation;                            //Mutation probabilities
    m_prob_shift_mutation = probShiftMutation;
    m_prob_permutation_mutation = probPermutationMutation;
    m_prob_inversion_mutation = probInversionMutation;
    m_prob_crossover = probCrossover;                                   //crossoover probability
}

void Generation::CreateInitialPopulation()
{

    for (int i = 0; i < m_M; i++)
    {
        m_gene_pool.row(i) = m_chromo->NewChromosome().t();
        // m_chromo->PrintChromosome(m_gene_pool.row(i).t());
        m_fitness.row(i) = m_chromo -> Fitness(m_gene_pool.row(i).t(), m_cities);

    }

    // cout << m_fitness << endl;
    // cout << m_cities << endl;
}

void Generation::UpdateFitness(){
    for (int i = 0; i < m_M; i++)
    {
         m_fitness.row(i) = m_chromo->Fitness(m_gene_pool.row(i).t(), m_cities);
    }
}

void Generation::Sort(){

    arma::uvec sorted_indices = arma::sort_index(m_fitness.col(0)); // Sort based on the cost function in m_fitness

    // Apply the sorted indices to m_gene_pool and m_fitness
    m_gene_pool = m_gene_pool.rows(sorted_indices);
    m_fitness = m_fitness.rows(sorted_indices);

    // cout << m_gene_pool << endl;
    // cout << m_fitness << endl;
}

void Generation::CitiesOnACircle()
{
    ofstream Cities;
    Cities.open("CitiesCircle.txt");

    for (int i = 0; i < m_cities.n_rows; i++)
    {

        double theta = m_rnd->Theta() * 4;
        m_cities(i, 0) = cos(theta);
        m_cities(i, 1) = sin(theta);
        Cities << m_cities.row(i) << endl;
    }

    Cities.close();

}

void Generation::CitiesOnASquare(){

    ofstream Cities;
    Cities.open("CitiesSquareSEED1.txt");

    for (int i = 0; i < m_cities.n_rows; i++)
    {
        m_cities(i, 0) = m_rnd->Rannyu();
        m_cities(i, 1) = m_rnd->Rannyu();
        Cities << m_cities.row(i) << endl;
    }
    Cities.close();

}

void Generation::AmericanCities(){
    ifstream Cities;
    cout << "AMERICAN TSP"<<endl;
    double longitude, latitude;
    Cities.open("capitals.dat");
    for (int i = 0; i < m_cities.n_rows; i++)
    {
        Cities >> longitude >> latitude;
        cout << longitude << " "<< latitude <<endl;
        m_cities(i, 0) = longitude;
        m_cities(i, 1) = latitude;
    }
    Cities.close();
}

void Generation::Mutations(int index){

    /*
    There are four types of mutations:
        -)SWAP: randomly selects to genes within a chromosome and swaps them
        -)SHIFT: m cities are shifted of n positions. (n,m) are randomly drawn
        -)PERMUTATION: m neighboring genes are permuted between each other (m) is randomly drawn
        -)INVERSION: the direction of m cities is inverted (m) is drawn randomly

    Each of this mutations has a relative probability which is fixed from the beginning of the algorithm.

    */

    int m;      // number of elements to be permuted
    int start1; // index of the first m contiguous elements to be permuted
    int start2;
    int n;     // number of positions to be shifted
    int start; // first one to be shifted

    int reducedSize = m_chromo->GetNumberOfGenes() - 2;

    //SWAP mutation
    int firstIndex;
    int secondIndex;


    if (m_rnd->Rannyu() <= m_prob_swap_mutation){

        // Sample two random positions: Allowed swaps are between the second element and the second to last element.
        firstIndex = int((m_rnd->Rannyu()) * (m_chromo->GetNumberOfGenes() - 2) + 1);
        secondIndex = int((m_rnd->Rannyu()) * (m_chromo->GetNumberOfGenes() - 2) + 1);

        // cout << "Before" << endl;
        // m_chromo->PrintChromosome(m_gene_pool.row(i).t());
        m_new_gen.row(index).swap_cols(firstIndex, secondIndex);
        // cout << "after" << endl;
        // m_chromo->PrintChromosome(m_gene_pool.row(i).t());
        m_chromo->CheckChromosome(m_new_gen.row(index).t());
        // m_fitness.row(index) = m_chromo->Fitness(m_gene_pool.row(index).t(), m_cities);
    }
    

    //SHIFT mutation
    if (m_rnd->Rannyu() <= m_prob_shift_mutation){

        start = int(m_rnd->Rannyu(0, reducedSize - 2)); // \in [0, L-2] as it is an index and starting from last position makes no sense
        n = int(m_rnd->Rannyu(1, reducedSize - start)); // \in [1,L-start]

        if (reducedSize - (start + n) > 0){

            m = int(m_rnd->Rannyu(1, reducedSize - (start + n)));

            // cout << "Before" << endl; 
            // m_chromo->PrintChromosome(m_gene_pool.row(i).t());
            // cout << "start: " << start << " n: " << n << " m: " << m << endl;

            // auxiliary arrays
            arma::vec temp = m_new_gen.row(index).t();
            arma::vec startingVec =temp.subvec(1, reducedSize);
            arma::vec copyVec = arma::zeros(reducedSize);

            //copy what remains unchanged
            for (int j = 0; j < start; j++){
                copyVec[j] = startingVec[j];
                startingVec[j] = 0;
            }

            //shift the block
            for (int j = start; j < start + m; j++){
                copyVec[j + n] = startingVec[j];
                startingVec[j] = 0;
            }

            //find the elements which have not been moved
            arma::uvec zero_indicesCopy = arma::find(copyVec == 0);
            arma::uvec non_zero_indicesStarting = arma::find(startingVec != 0);


            for (int j = 0; j < zero_indicesCopy.size(); j++){

                copyVec[zero_indicesCopy[j]] = startingVec[non_zero_indicesStarting[j]];
            }

            //save the result
            temp.subvec(1,reducedSize) = copyVec;
            m_new_gen.row(index) = temp.t();
            // cout << "after" << endl;
            // m_chromo->PrintChromosome(m_gene_pool.row(i).t());
            m_chromo->CheckChromosome(m_new_gen.row(index).t());
            // m_fitness.row(index) = m_chromo->Fitness(m_gene_pool.row(index).t(), m_cities);
        }

    }

    // PERMUTATION mutation
    if (m_rnd->Rannyu() <= m_prob_permutation_mutation) {

        start1 = int(m_rnd->Rannyu(1, int( 0.5 * reducedSize +1 ) -1 ));    // \in [0, L/2 -1] as it is an index and starting from last position makes no sense
        m = int(m_rnd->Rannyu(1, int(reducedSize * 0.5) - start1));      // \in [1,L/2-start]
        start2 = int(m_rnd->Rannyu(int( 0.5 * reducedSize +1 ), reducedSize - m -1));       // \in [L/2, L-m-1]  as it is an index

        // cout << "Before" << endl;
        // m_chromo->PrintChromosome(m_gene_pool.row(index).t());
        // cout << "start1: " << start1 << " start2 " << start2 << " m: " << m << endl;

        // auxiliary arrays
        arma::vec temp = m_new_gen.row(index).t();
        arma::vec startingVec1 = temp.subvec(start1, m + start1 -1);
        arma::vec startingVec2 = temp.subvec(start2, m + start2 -1);
        arma::vec concatenated = arma::join_cols(startingVec1, startingVec2);

        //permute
        concatenated = arma::shuffle(concatenated);

        //put the elements in the original vector
        temp.subvec(start1, m + start1 -1) = concatenated.subvec(0,m -1);
        temp.subvec(start2,m + start2 -1) = concatenated.subvec(m,m + m -1);

        m_new_gen.row(index) = temp.t();
        // cout << "after" << endl;
        // m_chromo->PrintChromosome(m_gene_pool.row(index).t());
        // m_chromo->CheckChromosome(m_new_gen.row(index).t());
        // m_fitness.row(index) = m_chromo->Fitness(m_gene_pool.row(index).t(), m_cities);
    }


    // INVERSION mutation
    if (m_rnd->Rannyu() <= m_prob_inversion_mutation){

        start1 = int(m_rnd->Rannyu(1, reducedSize +1 -1)); // \in [1, L -1] as it is an index and starting from last position makes no sense where L is the reduced one
        m = int(m_rnd->Rannyu(1, reducedSize - start1)); //number od cities to invert

        // cout << "Before" << endl;
        // m_chromo->PrintChromosome(m_gene_pool.row(i).t());
        // cout << "start1: " << start1 <<  " m: " << m << endl;

        // auxiliary arrays
        arma::vec temp = m_new_gen.row(index).t();
        arma::vec startingVec = temp.subvec(start1, m + start1 - 1);
        temp.subvec(start1, m + start1 - 1) = arma::reverse(startingVec);

        m_new_gen.row(index) = temp.t();
        // cout << "after" << endl;
        // m_chromo->PrintChromosome(m_gene_pool.row(i).t());
        m_chromo->CheckChromosome(m_new_gen.row(index).t());
        // m_fitness.row(index) = m_chromo->Fitness(m_gene_pool.row(index).t(), m_cities);

    }
}

void Generation::CrossOver(){
    // Crossover mutation

    int momIndex;
    int dadIndex;
    int cutPosition;
    int reducedSize = m_chromo->GetNumberOfGenes() - 2;

    //initialize new generation
    m_new_gen = arma::mat(m_M, m_chromo->GetNumberOfGenes());
    // cout << m_new_gen << endl;

    for (int k = 0; k < int(m_M/2); k+=2) {

        momIndex = Selection(); // Selects two parents with the loaded die (they could be the same, weird :D)
        dadIndex = Selection();

        if (m_rnd->Rannyu() <= m_prob_crossover){

            cutPosition = int(m_rnd->Rannyu(1, reducedSize +1 -1)); //\in [1, L-1]

            //auxiliary vectors
            arma::vec tempMom = m_gene_pool.row(momIndex).t();
            arma::vec tempDad = m_gene_pool.row(dadIndex).t();

            //Split vectors for the crossOver
            arma::vec tempMomSubVec = tempMom.subvec(1, cutPosition);
            arma::vec MissingMom = tempMom.subvec(cutPosition +1, reducedSize);
            // cout <<"Cut index: "<< cutPosition << endl;
            // cout << "mom" << endl;
            // m_chromo->PrintChromosome(tempMom);
            // m_chromo->PrintChromosome(MissingMom);
            // m_chromo->PrintChromosome(tempMomSubVec);
            arma::vec tempDadSubVec = tempDad.subvec(1, cutPosition);
            arma::vec MissingDad = tempDad.subvec(cutPosition + 1, reducedSize);
            // cout << "Dad" << endl;
            // m_chromo->PrintChromosome(tempDad);
            // m_chromo->PrintChromosome(MissingDad);
            // m_chromo->PrintChromosome(tempDadSubVec);

            // auxiliary vectors which will contain the positions of the missing part in the consort
            arma::uvec indexisMum;
            arma::uvec indexisDad;


            for (int j = 0; j < MissingMom.n_elem; j++){
                
                //Look into dad for the missing mum part
                int missingElemMum = MissingMom[j];
                arma::uvec matchedIndexMum = arma::find(tempDad==missingElemMum);
                indexisMum.insert_rows(indexisMum.n_elem, matchedIndexMum);

                //Look into mum for the missing dad part
                int missingElemDad = MissingDad[j];
                arma::uvec matchedIndexDad = arma::find(tempMom == missingElemDad);
                indexisDad.insert_rows(indexisDad.n_elem, matchedIndexDad);
            }

            //Sort them
            indexisMum = arma::sort(indexisMum);
            indexisDad = arma::sort(indexisDad);

            //Put the values in the right order
            for (int j = 0; j < MissingMom.n_elem; j++)
            {
                MissingMom[j] = tempDad[indexisMum[j]];
                MissingDad[j] = tempMom[indexisDad[j]];
            }
            // cout << "New missing mum " <<endl;
            // m_chromo->PrintChromosome(MissingMom);
            // cout << "New missing dad " << endl;
            // m_chromo->PrintChromosome(MissingDad);

            //store the new changes
            tempMom.subvec(cutPosition + 1, reducedSize) = MissingMom;
            tempDad.subvec(cutPosition + 1, reducedSize) = MissingDad;

            // cout << "mom" << endl;
            // m_chromo->PrintChromosome(tempMom);
            // cout << "Dad" << endl;
            // m_chromo->PrintChromosome(tempDad);

            // Routine check
            m_chromo->CheckChromosome(tempMom);
            m_chromo->CheckChromosome(tempDad);

            m_new_gen.row(k) = tempMom.t();
            m_new_gen.row(k+1) = tempDad.t();
            
        }

        else{
            //Just use the same chromosomes
            m_new_gen.row(k) = m_gene_pool.row(momIndex);
            m_new_gen.row(k+1) = m_gene_pool.row(dadIndex);
        }


        //Mutations
        Mutations(k);
        Mutations(k+1);

    }

    //Update the new generation: first fint all the chromosomes which did not undergo a change
    arma::uvec zeroRows = arma::find(arma::all(m_new_gen == 0, 1));
    //Put in this rows the unchanged chromosomes:
    for (int h = 0; h < zeroRows.n_elem; h++){
        m_new_gen.row(zeroRows[h]) = m_gene_pool.row(zeroRows[h]);
    }
    //finally update the current generation
    m_gene_pool = m_new_gen;

    // Compute fitness and sort gene_pool
    UpdateFitness();
    Sort();

}

double Generation::bestHalfAverage(){

    arma::vec firstHalf = m_fitness.rows(0, int(m_fitness.n_rows / 2));
    return arma::mean(firstHalf);
}

arma::mat Generation::BestPath(){
    arma::vec bestPath = m_gene_pool.row(0).t();
    arma::mat bestCities = arma::mat(m_chromo->GetNumberOfGenes() - 1, 2);
    
    for (int i =0; i<bestPath.n_elem-1; i++){
        bestCities(i,0) = m_cities(bestPath[i]-1,0);
        bestCities(i, 1) = m_cities(bestPath[i]-1, 1);
    }
    return bestCities;
}

//=============================================================//

GeneticAlgorithm::GeneticAlgorithm(Random * rnd,int CircleSquare, int NumberOfGenerations, int M, double probSwapMutation, double probShiftMutation, double probPermutationMutation, double probInversionMutation, double probCrossover, double p)
{
    m_circlesquare = CircleSquare;
    m_NumberOfGenerations = NumberOfGenerations;                                                                         // number of evolutions to realize
    m_gen = new Generation(rnd, M, probSwapMutation, probShiftMutation, probPermutationMutation, probInversionMutation, probCrossover, p); // M is the size oof the first generation
    m_rnd = rnd;            
}

void GeneticAlgorithm::Evolve(){
    ofstream Results, Path, Loss;
    Results.open("resultsSquareSEED1.txt");
    Path.open("pathsSquareSEED1.txt");
    Loss.open("optimumLossSquareSEED1.txt");

    for (int i=0; i < m_NumberOfGenerations; i++){
        cout<< "Gen: " << i +1 <<endl;
        m_gen->CrossOver();                 //This makes crossover and mutations all at once
        Results << m_gen->bestHalfAverage() << endl;
        Loss << i << " " << m_gen->GetOptimumLoss()<<endl;
        if (i%25==0){
            arma::mat bestPath = m_gen->BestPath();
            for (int j = 0; j < bestPath.n_rows; j++)
            {
                Path << i << " " << bestPath.row(j) << endl;
            }
        }
        
    }
    Results.close();
    Path.close();
    Loss.close();
        
}

void GeneticAlgorithm::Test()
{
    m_gen->CrossOver();
}

void GeneticAlgorithm::Start(){

    if (m_circlesquare==1){
        m_gen->CitiesOnASquare();       // If you want them distributed over a square
    }
    else if (m_circlesquare==0) {
        m_gen->CitiesOnACircle();
    }
    else{
        m_gen->AmericanCities();
    }

    m_gen->CreateInitialPopulation();   // Creates the initial chromosome population
    m_gen->Sort();                      // Sorts the chromosomes according to path length

}


