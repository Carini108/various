#include "GA.h"

using namespace std;

int main(int argc, char* argv[]) { 

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // initialization
    Input();

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // get cities filename from command line and do a check
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }
    string filename = argv[1];
    vector<City> cities = readCitiesFromFile(filename);
    int number_of_cities = cities.size();

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // display cities
    cout << "The salesman has to visit the following "<< number_of_cities << " cities:" << endl;
    for (const City& city : cities) { // use iterators
        printCity(city);
    }
    cout << "===================================" << endl;
    
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // initialize population
    vector<vector<int>> Population = initializePopulation(population_size, number_of_cities, rnd);
    if (checkPopulation(Population,number_of_cities)) { cout << "Population properly initialized" << endl; } 
    else { cout << "Error during initialization! Exiting... " << endl; exit(0); }
    cout << "Starting population:" << endl;
    printPopulation(Population);
    cout << "===================================" << endl;
    
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // advance generation
    cout << "Commencing evolution process:" << endl;
    for (unsigned int gen = 1; gen<=number_of_generations; gen++ ) {
        //////////////////////////////////////////////////////////////
        // user update: show step
        cout << "Generation " << gen << "/" << number_of_generations << endl;
        //////////////////////////////////////////////////////////////
        // mutations
        for (auto& individual : Population) { // use iterators to mutate
            if (rnd.Rannyu()<P_pair_permutation) // mutate w a certain probability 
            {
                GenMut_PairPerm(individual,rnd,number_of_cities);
            }
            if (rnd.Rannyu()<P_shift)
            {
                GenMut_Shift(individual,rnd,number_of_cities);
            }
            
            if (rnd.Rannyu()<P_inversion) // mutate w a certain probability 
            {
                GenMut_Inversion(individual,rnd,number_of_cities);
            }
            if (rnd.Rannyu()<P_group_permutation) // mutate w a certain probability 
            {
                GenMut_GroupsPermutation(individual, rnd, number_of_cities);
            }
        }
        //////////////////////////////////////////////////////////////
        // population check
        if (!checkPopulation(Population, number_of_cities)){
            cout << "MUTATION WARNING: something went wrong at step "<< gen << endl;
            exit(0);
        }
        //////////////////////////////////////////////////////////////
        // prepare to kill some individuals by assigning unfitness and sorting the population
        vector<double> UP = calculateUnfitnessParameters(Population, cities, which_loss);
        UP = calculateUnfitnessParameters(Population, cities, which_loss);
        sortPopulationByUnfitness(Population, UP);
        UP.clear(); // dispose of the useless variables to free RAM
        //////////////////////////////////////////////////////////////
        // sort probabilities as well
        vector<double> sorted_UP = calculateUnfitnessParameters(Population, cities, which_loss);
        //////////////////////////////////////////////////////////////
        // remove individuals based on random numbers and unfitness parameters
        int number_of_deaths = population_size*0.4; 
        killSomeIndividuals(Population, sorted_UP, number_of_deaths, rnd);
        //////////////////////////////////////////////////////////////
        // replace dead individuals with breeded survivors
        breedOtherIndividuals(Population, number_of_deaths, rnd);
        //////////////////////////////////////////////////////////////
        // sort by increasing fitness / decreasing unfitness
        UP.clear();
        UP = calculateUnfitnessParameters(Population, cities, which_loss);
        sortPopulationByUnfitness(Population, UP);
        UP.clear(); // delete then...
        UP = calculateUnfitnessParameters(Population, cities, which_loss);
        if ( !checkPopulation(Population, number_of_cities) || population_size!=Population.size() ){
            cout << "BREEDING WARNING: something went wrong at step "<< gen << endl;
        }
        //////////////////////////////////////////////////////////////
        // show fittest individual
        cout << "Shortest route found = "<< Loss(Population[population_size-1], cities, which_loss) << endl;
        //////////////////////////////////////////////////////////////
    }

    //////////////////////////////////////////////////////////////
    // print population
    vector<double> UP = calculateUnfitnessParameters(Population, cities, which_loss);
    sortPopulationByUnfitness(Population, UP);
    cout << "Final population" << endl;
    printPopulation(Population);
    // write fittest individual to file
    string output_filename = "bestCIRCLEloss"+to_string(which_loss)+".dat";
    writeFittestToFile(Population.back(), cities, output_filename);
    // show fittest individual
    cout << "Shortest route found = "<< Loss(Population[population_size-1], cities, which_loss) << endl;

    cout << "===================================" << endl;
    cout << "========        END        ========" << endl;
    cout << "===================================" << endl;

    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
// function to read run card
//////////////////////////////////////////////////////////////
void Input(void) {
    cout << "===================================" << endl;
    cout << "======== GENETIC ALGORITHM ========" << endl;
    cout << "========        TSP        ========" << endl;
    cout << "===================================" << endl;
    //////////////////////////////////////////////////////////////
    // read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    if (!Primes.is_open()) {
        cerr << "Error opening Primes! " << endl;
        exit(0);
    }
    Primes >> p1 >> p2;
    Primes.close();
    ifstream input("seed.in");
    if (!input.is_open()) {
        cerr << "Error opening seed.in! " << endl;
        exit(0);
    }
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed, p1, p2);
    input.close();
    //////////////////////////////////////////////////////////////
    // read run card
    ifstream ReadRunCard;
    ReadRunCard.open("run_card.dat");
    if (!ReadRunCard.is_open()) {
        cerr << "Error opening run_card.dat! " << endl;
        exit(0);
    }
    // population size
    ReadRunCard >> population_size;
    cout << "Population size = " << population_size << endl;
    if (population_size<=0) {
        cout << "Select valid population size!" << endl;
        cout << "Exiting..." << endl;
        exit(0);
    }
    // loss function preferred
    ReadRunCard >> which_loss;
    if (which_loss == 1 || which_loss == 2) {
        cout << "Loss function in use = " << which_loss << endl;
    } else {
        cout << "Select a valid loss function! (1 or 2)" << endl;
        cout << "Exiting..." << endl;
        exit(0);
    }
    // probabilities of mutations and crossover
    bool badProbabilities = false;
    // pair perm
    ReadRunCard >> P_pair_permutation;
    cout << "Probability of pair permutation = " << P_pair_permutation << endl;
    if (P_pair_permutation<0||P_pair_permutation>1) { badProbabilities=true; }
    // shift
    ReadRunCard >> P_shift;
    cout << "Probability of shift = " << P_shift << endl;
    if (P_shift<0||P_shift>1) { badProbabilities=true; }
    // group perm
    ReadRunCard >> P_group_permutation;
    cout << "Probability of group permutation = " << P_group_permutation << endl;
    if (P_group_permutation<0||P_group_permutation>1) { badProbabilities=true; }
    // inversion
    ReadRunCard >> P_inversion;
    cout << "Probability of inversion = " << P_inversion << endl;
    if (P_inversion<0||P_inversion>1) { badProbabilities=true; }
    // crossover
    ReadRunCard >> P_crossover;
    cout << "Probability of crossover = " << P_crossover << endl;
    if (P_crossover<0||P_crossover>1) { badProbabilities=true; }
    // check probabilities
    if (badProbabilities) {
        cout << "Select valid probabilities! (between 0 and 1)" << endl;
        cout << "Exiting..." << endl;
        exit(0);
    }
    // number of generations
    ReadRunCard >> number_of_generations;
    cout << "Number of generations = " << number_of_generations << endl;
    if (number_of_generations<=0) {
        cout << "Select valid number of generations!" << endl;
        cout << "Exiting..." << endl;
        exit(0);
    }
    // closing file reading stream
    ReadRunCard.close();
    cout << "===================================" << endl;
}

//////////////////////////////////////////////////////////////
// function to read cities from file
//////////////////////////////////////////////////////////////
vector<City> readCitiesFromFile(const string& filename) {
    vector<City> cities;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening cities file: " << filename << endl;
        exit(0);
    }
    unsigned int id;
    double x, y;
    while (file >> id >> x >> y) {
        cities.emplace_back(id, x, y); // more efficient than move constructor!
    }
    file.close();
    return cities;
}

//////////////////////////////////////////////////////////////
// function to initialize the population
//////////////////////////////////////////////////////////////
vector<vector<int>> initializePopulation(int populationSize, int numberOfCities, Random& rnd) {
    vector<vector<int>> population;
    for (int i = 0; i < populationSize; ++i) {
        vector<int> permutation;
        permutation.push_back(1); // always start with 1 at the first position
        // fill the permutation with unique integers from 2 to numberOfCities
        for (int j = 2; j <= numberOfCities; ++j) {
            bool unique = true;
            int newRandomNumber;
            // Generate a random number until it's unique
            do {
                newRandomNumber = rnd.Rannyu(2, numberOfCities + 1); // Generate a random number from 2 to numberOfCities
                if (find(permutation.begin(), permutation.end(), newRandomNumber) != permutation.end()) { // the result is the last position of the iterator; if the iterator ends up at .end(), the number was not in the list
                    unique = false; // newRandomNumber already exists in the permutation
                } else {
                    unique = true; // newRandomNumber is unique
                }
            } while (!unique);
            permutation.push_back(newRandomNumber);
        }
        population.push_back(permutation);
    }
    return population;
}

//////////////////////////////////////////////////////////////
// function to check whether an individual permutation 
// satisfies the properties of the problem
//////////////////////////////////////////////////////////////
bool checkIndividual(const vector<int>& individual, int numberOfCities) {
    // check whether the size of the individual permutation is correct
    if (static_cast<int>(individual.size()) != numberOfCities) {
        return false;
    }
    // check whether the first element is 1
    if (individual[0] != 1) {
        return false;
    }
    // check whether all integers from 2 to numberOfCities appear exactly once in the permutation
    vector<int> sortedIndividual = individual; // copy the individual permutation to sort
    sort(sortedIndividual.begin(), sortedIndividual.end()); // sort it...
    for (int i = 2; i <= numberOfCities; ++i) { // ...and see if every number is there
        if (sortedIndividual[i - 1] != i) { // easy check if vector is sorted :P
            return false;
        }
    }
    // ok, if all conditions satisfied :)
    return true;
}

//////////////////////////////////////////////////////////////
// function to check the entire population
//////////////////////////////////////////////////////////////
bool checkPopulation(const vector<vector<int>>& population, int numberOfCities) {
    for (const auto& individual : population) {
        if (!checkIndividual(individual, numberOfCities)) {
            return false;
        }
    }
    // if all individuals are okay, so is the population :):):):):)
    return true;
}

//////////////////////////////////////////////////////////////
// function to calculate loss
//////////////////////////////////////////////////////////////
double Loss(const vector<int>& permutation, const vector<City>& cities, const unsigned int which_loss) {
    double totalDistance = 0.0;
    int N = permutation.size();
    // take the sum of distances between consecutive cities
    for (int i = 0; i < N - 1; ++i) {
        unsigned int city1_id = permutation[i];     // once the city is identified 
        unsigned int city2_id = permutation[i + 1]; //by the number in the permutation...
        const City& city1 = cities[city1_id];       //...and the cities 'book' is parsed...
        const City& city2 = cities[city2_id];
        // ... the positions are readily recovered!
        double distance = (which_loss == 1) ? (sqrt(pow(city2.getX() - city1.getX(), 2.) + pow(city2.getY() - city1.getY(), 2.))) : (pow(city2.getX() - city1.getX(), 2.) + pow(city2.getY() - city1.getY(), 2.));
        
        totalDistance += distance;
    }
    // last term: distance from the last city back to the starting city
    unsigned int first_city_id = permutation[0];
    unsigned int last_city_id = permutation[N - 1];
    const City& first_city = cities[first_city_id];
    const City& last_city = cities[last_city_id];
    totalDistance += (which_loss == 1) ? (sqrt(pow(last_city.getX() - first_city.getX(), 2.) + pow(last_city.getY() - first_city.getY(), 2.))) : (pow(last_city.getX() - first_city.getX(), 2.) + pow(last_city.getY() - first_city.getY(), 2.));
    // all terms summed
    return totalDistance;
}

//////////////////////////////////////////////////////////////
// function to print a single individual
//////////////////////////////////////////////////////////////
void printIndividual(const vector<int>& individual) {
    cout << "\t| ";
    for (int city_id : individual) {
        cout << city_id << " ";
    }
    cout << "| " << endl;
}

//////////////////////////////////////////////////////////////
// function to print the entire population
//////////////////////////////////////////////////////////////
void printPopulation(const vector<vector<int>>& population) {
    unsigned int counter = 1;
    for (const auto& individual : population) {
        cout << counter++ << "\t- th path";
        printIndividual(individual);
    }
}

//////////////////////////////////////////////////////////////
// function to print details of a single city
//////////////////////////////////////////////////////////////
void printCity(const City& city) {
    cout << "\tCity ID: " << city.getId() << ",\tX: " << city.getX() << ",\tY: " << city.getY() << endl;
}

//////////////////////////////////////////////////////////////
// function to generate a pair permutation mutation
//////////////////////////////////////////////////////////////
void GenMut_PairPerm(vector<int>& individual, Random& rnd, int numberOfCities) {
    // note: NEVER swap position 0 (there must be always 1)
    int position = rnd.Rannyu(1, numberOfCities - 1); // result: integer between 1 and numberOfCities-2
    swap(individual[position], individual[position + 1]); // swap from STL
    if (!checkIndividual(individual, numberOfCities)) { // check if the mutated individual is valid
        cerr << "Mutation PairPerm resulted in an invalid individual!" << endl;
        // exit(0);
    }
}

//////////////////////////////////////////////////////////////
// function to generate a shift mutation
//////////////////////////////////////////////////////////////
void GenMut_Shift(vector<int>& individual, Random& rnd, int numberOfCities) {
    // ensure m is less than numberOfCities - 1
    int maxM = numberOfCities - 1 - 1; 
    int m = rnd.Rannyu(1, maxM + 1); // extract block size
    int maxN = numberOfCities - m - 1;
    int n = rnd.Rannyu(1, maxN + 1); // extract shift amount
    // define block start position
    int blockStart = rnd.Rannyu(1, numberOfCities - m - n); // start position of the block
    // copy the block to be shifted
    vector<int> block(individual.begin() + blockStart, individual.begin() + blockStart + m);
    // perform the shift by inserting block at the new position and removing the original block
    individual.erase(individual.begin() + blockStart, individual.begin() + blockStart + m);
    for (int i = 0; i < m; i++)
    {
        individual.insert(individual.begin() + blockStart + n + i, block[i] );
    }
    // check if the mutation produced a valid individual
    if (!checkIndividual(individual, numberOfCities)) {
        cout << "Mutation Shift resulted in an invalid individual!" << std::endl;
        exit(0);
    }
}

//////////////////////////////////////////////////////////////
// function to generate a group permutation
//////////////////////////////////////////////////////////////
void GenMut_GroupsPermutation(vector<int>& individual, Random& rnd, int numberOfCities) {
    int maxBlockSize = numberOfCities / 2 - 1; 
    int m = rnd.Rannyu(1, maxBlockSize + 1); // extract block size
    // define two non-overlapping blocks
    int startBlock1 = rnd.Rannyu(1, numberOfCities - 2 * m); // start of first block
    int startBlock2 = rnd.Rannyu(startBlock1 + m, numberOfCities - m); // start of second block
    // ensure blocks are non-overlapping
    if (startBlock1 + m > startBlock2) {
        startBlock2 = startBlock1 + m;
    }
    // copy the blocks
    vector<int> block1(individual.begin() + startBlock1, individual.begin() + startBlock1 + m);
    //cout << "block 1 = " ;
    //printIndividual(block1);
    vector<int> block2(individual.begin() + startBlock2, individual.begin() + startBlock2 + m);
    //cout << "block 2 = " ;
    //printIndividual(block2);
    // perform the swap
    copy(block2.begin(), block2.end(), individual.begin() + startBlock1);
    copy(block1.begin(), block1.end(), individual.begin() + startBlock2);
    // check if the mutation produced a valid individual
    if (!checkIndividual(individual, numberOfCities)) {
        cerr << "Mutation PermutationContiguousCities resulted in an invalid individual!" << std::endl;
        exit(0);
    }
}

//////////////////////////////////////////////////////////////
// function to generate a inversion mutation
//////////////////////////////////////////////////////////////
void GenMut_Inversion(vector<int>& individual, Random& rnd, int numberOfCities) {
    // Define block size
    int maxBlockSize = numberOfCities - 2; // excluding the first city (and the last, otherwise, there is no space for inversion...)
    int m = rnd.Rannyu(2, maxBlockSize + 1); // extract block size
    // Define block start position
    int blockStart = rnd.Rannyu(1, numberOfCities - m); // start position of the block (1 <= blockStart <= numberOfCities-m-1)
    // Perform inversion
    for (int i = 0; i < m / 2; ++i) { // iterate over half of the block size
        swap(individual[blockStart + i], individual[blockStart + m - 1 - i]); // swap elements from the beginning and end of the block
    }
    if (!checkIndividual(individual, numberOfCities)) {
        cout << "Mutation Inversion resulted in an invalid individual!" << endl;
    }
}

//////////////////////////////////////////////////////////////
// swap genetic material between parents
//////////////////////////////////////////////////////////////
void doCrossover(vector<int> &parent1, vector<int> &parent2, Random &rnd) {
  int numberOfCities = parent1.size();
  int numberOfCities_alias = parent2.size();
  if (numberOfCities!=numberOfCities_alias)
  {
    cout << "Uncompatible parents length!" << endl;
    exit(0);
  }  
  // choose a random crossover point
  int crossoverPoint = rnd.Rannyu(1, numberOfCities); // not (0,numberOfCities+1), to ensure both parts are non-empty
  // create temporary vectors to store the new offspring
  vector<int> offspring1(parent1.begin(), parent1.begin() + crossoverPoint);
  vector<int> offspring2(parent2.begin(), parent2.begin() + crossoverPoint);
  // fill in the missing cities from the second part of the opposite parent
  for (int i = 0; i < numberOfCities; ++i) {
    if (find(offspring1.begin(), offspring1.end(), parent2[i]) == offspring1.end()) {
      offspring1.push_back(parent2[i]);
    }
    if (find(offspring2.begin(), offspring2.end(), parent1[i]) == offspring2.end()) {
      offspring2.push_back(parent1[i]);
    }
  }
  if (!checkIndividual(offspring1, numberOfCities) || !checkIndividual(offspring2, numberOfCities)) {
    cout << "Crossover resulted in an invalid individual!" << endl;
    exit(0);
  }
  // overwrite the parents with the new offspring
  parent1 = offspring1;
  parent2 = offspring2;
}

//////////////////////////////////////////////////////////////
// function to calculate unfitness parameters for the population
//////////////////////////////////////////////////////////////
vector<double> calculateUnfitnessParameters(const vector<vector<int>>& population, const vector<City>& cities, unsigned int which_loss) {
    vector<double> unfitnessParameters(population.size());
    double totalLoss = 0.0;
    for (size_t i = 0; i < population.size(); ++i) {
        unfitnessParameters[i] = Loss(population[i], cities, which_loss); // individual loss
        totalLoss += unfitnessParameters[i]; // sum to get total loss (for normalization...)
    }
    for (size_t i = 0; i < population.size(); ++i) { // normalize individual losses to get death probabilities
        unfitnessParameters[i] /= totalLoss;
    }
    return unfitnessParameters;
}

//////////////////////////////////////////////////////////////
// function to print unfitness parameters for the population
//////////////////////////////////////////////////////////////
void printUP(const vector<double>& UnfitnessParameters) {
    cout << "UP:\t( ";
    for (auto x : UnfitnessParameters) {
        cout << x << " ";
    }
    cout << ")" << endl;
}

//////////////////////////////////////////////////////////////
// function to sort the population by decreasing unfitness parameters
//////////////////////////////////////////////////////////////
void sortPopulationByUnfitness(vector<vector<int>>& population, const vector<double>& unfitnessParameters) {
    // Create a vector of indices for the population
    vector<size_t> indices(population.size());
    for (size_t i = 0; i < population.size(); ++i) {
        indices[i] = i;
    }
    // Sort the indices based on unfitness parameters in descending order
    sort(indices.begin(), indices.end(), [&](size_t a, size_t b) { return unfitnessParameters[a] > unfitnessParameters[b]; });
    // Reorder the population based on sorted indices
    vector<vector<int>> sortedPopulation(population.size());
    for (size_t i = 0; i < population.size(); ++i) {
        sortedPopulation[i] = population[indices[i]];
    }
    // Update the original population with the sorted population
    population = sortedPopulation;
}

//////////////////////////////////////////////////////////////
// function to remove individuals from the population based on random numbers and unfitness parameters
//////////////////////////////////////////////////////////////
void killSomeIndividuals(vector<vector<int>>& Population, const vector<double>& sorted_UP, int numToRemove, Random& rnd) {
    //int originalSize = Population.size();
    double cumulativeProbability = 0.0;
    for (int i = 0; i < numToRemove; ++i) {
        double randomNumber = rnd.Rannyu();
        cumulativeProbability = 0.0;
        for (unsigned int j = 0; j < sorted_UP.size(); ++j) {
            cumulativeProbability += sorted_UP[j];
            if (randomNumber < cumulativeProbability) {
                //cout << "pulled the trigger now he's dead"<< endl;
                Population.erase(Population.begin());
                break;
            }
        }
    }
    //cout << "i've killed "<< originalSize-Population.size() << " individuals" << endl;
}

//////////////////////////////////////////////////////////////
// function to breed and add new individuals to replace the dead ones
//////////////////////////////////////////////////////////////
void breedOtherIndividuals(vector<vector<int>>& Population, int numToRemove, Random& rnd) {
    //cout << "initial size " << Population.size() << endl;
    vector<vector<int>> fittestIndividuals(Population.end() - numToRemove, Population.end()); // Extract the fittest individuals
    for (unsigned int i = 0; i < fittestIndividuals.size()-1; i++)
    {
        if (rnd.Rannyu() < P_crossover){
            //cout << "Crossing over..." << endl;
            doCrossover( fittestIndividuals[i],fittestIndividuals[i+1],rnd );
        }
    }
    for (int i = 0; i < numToRemove; ++i) { // Add copies of fittest individuals to the population
        Population.push_back(fittestIndividuals[i]); // Add a random fittest individual
    }
    //cout << "after reproduction" << Population.size() << endl;
}

//////////////////////////////////////////////////////////////
// function to write result to file
//////////////////////////////////////////////////////////////
void writeFittestToFile(const vector<int>& fittest, const vector<City>& cities, const string& filename) {
    ofstream outFile(filename);
    cout << "Writing fittest individual to file "<< filename << endl;
    if (!outFile) {
        cerr << "Error opening file " << filename << " for writing!" << endl;
        return;
    }
    for (int cityIndex : fittest) {
        outFile << cityIndex << "\t";
    }
    outFile.close();
}