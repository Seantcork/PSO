/*

	Authors Sean Cork, Kamaal Palmer, Luca Ostertag-Hill
	Nature Inspired Computation Project 2: PSO
	October 13, 2018

*/

#include <math.h>
#include <fstream>
#include <cstdlib>
#include <stdint.h>
#include <random>
#include <string>
#include <vector>
#include <iostream>
#include <limits>
#include <set>

using namespace std;

#define E 2.71828
const double CONSTRICTION_FACTOR = 0.7298;

//phi values that scales towards pBest and nBest respectively
const double PHI_1 = 2.05;
const double PHI_2 = 2.05;

const string GLOBAL_TOPOLOGY = "gl";
const string RING_TOPOLOGY = "ri";
const string VON_NEUMANN_TOPOLOGY = "vn";
const string RANDOM_TOPOLOGY = "ra";
//neighborhood size for the random topology
const int RANDOM_K = 5;

const string ROSENBROCK_FUNCTION = "rok";
const string ACKLEY_FUNCTION = "ack";
const string RASTRIGIN_FUNCTION = "ras";

//Particle starting range of position by function
const double ROSENBROCK_INIT_MAX_POSITION = 30.0;
const double ROSENBROCK_INIT_MIN_POSITION = 15.0;
const double ACKLEY_INIT_MAX_POSITION = 32.0;
const double ACKLEY_INIT_MIN_POSITION = 16.0;
const double RASTRIGIN_INIT_MAX_POSITION = 5.12
const double RASTRIGIN_INIT_MIN_POSITION = 2.56;

//Particle starting range of velocty by function
const double ROSENBROCK_INIT_MAX_VELOCITY = 2.0;
const double ROSENBROCK_INIT_MIN_VELOCITY = -2.0;
const double ACKLEY_INIT_MAX_VELOCITY = 4.0;
const double ACKLEY_INIT_MIN_VELOCITY = -2.0;
const double RASTRIGIN_INIT_MAX_VELOCITY = 4.0;
const double RASTRIGIN_INIT_MIN_VELOCITY = -2.0;

//Particle min and max velocity by function
const double ROSENBROCK_MAX_VELOCITY = 2.048;
const double ROSENBROCK_MIN_VELOCITY = -2.048;
const double ACKLEY_MAX_VELOCITY = 32.768;
const double ACKLEY_MIN_VELOCITY = -32.768;
const double RASTRIGIN_MAX_VELOCITY = 5.12;
const double RASTRIGIN_MIN_VELOCITY = -5.12;

class Particle {
	/*
		Class to represents a particle in the swarm.

		Attributes:
			- maxVelocity (double): The maximum velocity a particle can move at. Determined
				by the testFunction the user specifies.
			- minVelocity (double): The minimum velocity a particle can move at. Determined
				by the testFunction the user specifies.

			- pbestFitness (double): The best fitness the particle itself has found so far.
			- pBestArray (vector<double>): The position of the particle at its best fitness so far.
			- position (vector<double>): The particle's position in d dimensions.
			- velocity (vector<double>): The particle's velocity in d dimensions.

			- nBestFitness (double): The best fitness of the particle's neighborhood so far.
			- nBestArray (vector<double>): The position of the best fitness of the neighborhood so far.
			- neighborsArray (vector<shared_ptr<Particle> >): An array containing pointers to the
				other particles in this particles neighborhood. 

		Functions:
			- void initParticle(int numDimensions, string testFunction): Function that initializes 
				the particle's min/max velocity, its starting position/velocity, and sets the pBest and 
				nBest to the starting location.

			- void calculateFitness(string testFunction): Determines the fitness of a particle at its
				current position. Takes in which function to use as the evaluation.
			- double evalAckley(): The Ackley evaluation function. Calculates fitness of the particle.
			- double evalRosenbrock(): The Rosenbrock evaluation function. Calculates fitness of the particle.
			- double evalRastrigin(): The Rastrigin evaluation function. Calculates fitness of the particle.

			- void updatePosition(): Uses the particle's current velocity to update its position.
			- void updateVelocity(): Uses the old velocity, pBest, and nBest to update its velocity.
			- void findNeighborhoodBest(): Searches through each particle in the neighborhood to determine
				the maximum neighborhood fitness and the location of that fitness.

	*/

	public:

		double maxVelocity;
		double minVelocity;

		double pBestFitness;
		vector<double> pBestArray;
		vector<double> position;
		vector<double> velocity;
		
		
		double nBestFitness;
		vector<double> nBestArray;
		vector<shared_ptr<Particle> > neighborsArray;

		void initParticle(int numDimensions, string testFunction);	

		void calculateFitness(string testFunction);
		double evalAckley();
		double evalRosenbrock();
		double evalRastrigin();

		void updatePosition();
		void updateVelocity();
		void findNeighborhoodBest();
};


class Swarm {
	/*
		Class that represents a swarm.

		Attributes:
			- swarmSize (int): The user specified size of the swarm.
			- numDimensions (int): The user specified number of dimensions to evaluate 
				the particle in.

			- gBestFitness (double): The best fitness of the swarm so far.
			- gBestArray (vector<double>): The location of the best swarm fitness.

		Functions:
			- void initSwarm(int swarmSize, int numDimensions, string neighborhoodTopology, 
				string testFunction): Function that initializes the swarm of particles and
				calls the user specified topology.
			- void findGlobalBest(): Determines if there is a new global best. If so, it updates
				the gBest fitness and location array.

			- void globalTopology(): Neighborhood topology that includes all particles in the swarm.
			- void ringTopology(): Neighborhood topology that includes the left and right neighbor.
			- void vonNeumanTopology(): Neighborhood topology that includes the left, right, top, 
				and bottom neighbor.
			- void randomTopology(): Neighborhood topology that contains random k-1 other particles. 
				Neighborhood gets updated with probability 0.2 for each particle after each iteration.
	*/

	public:
		vector<shared_ptr<Particle> > swarm;
		
		int swarmSize;
		int numDimensions;

		double gBestFitness;
		vector<double> gBestArray;

		void initSwarm(int swarmSize, int numDimensions, 
			string neighborhoodTopology, string testFunction);
		void findGlobalBest();

		void globalTopology();
		void ringTopology();
		void vonNeumanTopology();
		void randomTopology();
};

/*
	Function that initializes the particle's starting position and velocity for the user
	specified number of dimensions. This and the min/max velocity of the particle is dependent
	on the user specified evaluation functions. The specific ranges for position/velocity and 
	min/max velocity was given (and are constants). The function also initializes the particle's
	pBest and nBest fitness to max (so that all subsequent fitness are better).
*/
void Particle::initParticle(int numDimensions, string testFunction){
	random_device seeder;
	mt19937 engine(seeder());

	/* Checks to see which evaluation function the user specified. */
	if(testFunction.compare(ROSENBROCK_FUNCTION) == 0){
		
		//sets the min and max velocity for that function
		maxVelocity = ROSENBROCK_MAX_VELOCITY;
		minVelocity = ROSENBROCK_MIN_VELOCITY;

		//sets the position and velocity for that function based on the ranges
		uniform_real_distribution<double> genPosition(ROSENBROCK_INIT_MIN_POSITION, ROSENBROCK_INIT_MAX_POSITION);
		uniform_real_distribution<double> genVelocity(ROSENBROCK_INIT_MIN_VELOCITY, ROSENBROCK_INIT_MAX_VELOCITY);
		for(int i = 0; i < numDimensions; i ++){
			position.push_back(genPosition(engine));
			velocity.push_back(genVelocity(engine));
		}
	}

	else if(testFunction.compare(ACKLEY_FUNCTION) == 0){
		maxVelocity = ACKLEY_MAX_VELOCITY;
		minVelocity = ACKLEY_MIN_VELOCITY;

		uniform_real_distribution<double> genPosition(ACKLEY_INIT_MIN_POSITION, ACKLEY_INIT_MAX_POSITION);
		uniform_real_distribution<double> genVelocity(ACKLEY_INIT_MIN_VELOCITY, ACKLEY_INIT_MAX_VELOCITY);
		for(int i = 0; i < numDimensions; i ++){
			position.push_back(genPosition(engine));
			velocity.push_back(genVelocity(engine));

		}
	}
	
	else if(testFunction.compare(RASTRIGIN_FUNCTION) == 0){
		maxVelocity = RASTRIGIN_MAX_VELOCITY;
		minVelocity = RASTRIGIN_MIN_VELOCITY;

		uniform_real_distribution<double> genPosition(RASTRIGIN_INIT_MIN_POSITION, RASTRIGIN_INIT_MAX_POSITION);
		uniform_real_distribution<double> genVelocity(RASTRIGIN_INIT_MIN_VELOCITY, RASTRIGIN_INIT_MAX_VELOCITY);
		for(int i = 0; i < numDimensions; i ++){
			position.push_back(genPosition(engine));
			velocity.push_back(genVelocity(engine));

		}
	}

	/* If the user entered an invalid function name. */
	else {

		cerr << "Optimization Function does not exist" << endl;
	}

	//max out the starting pBest and nBest values
	pBestFitness = numeric_limits<double>::max();
	nBestFitness = numeric_limits<double>::max();
	
	//Particles pBest is set as its initial position
	pBestArray = position;

	//the nBestArray is also at its initial position
	nBestArray = position;
	
}


/*
	Function that calculates and potentially updates the pBest fitness of the particle.
	It checks which evaluation function to apply to determine the fitness. If the new
	fitness is smaller than the pBest fitness, it updates the fitness and the pBest location
	of the fitness.
*/
void Particle::calculateFitness(string testFunction){
	double newFitness;

	/* Determine which evaluation function to run */
	if(testFunction.compare(ROSENBROCK_FUNCTION) == 0) {
		newFitness = evalRosenbrock();

		//PSO wants to find the minimum
		if(newFitness < pBestFitness){
			pBestFitness = newFitness;
			pBestArray = position;
		}

	}

	else if (testFunction.compare(ACKLEY_FUNCTION) == 0) {
		newFitness = evalAckley();

		if(newFitness < pBestFitness){
			pBestFitness = newFitness;
			pBestArray = position;
		}
	}

	else if (testFunction.compare(RASTRIGIN_FUNCTION) == 0){
		newFitness = evalRastrigin();

		if(newFitness < pBestFitness) {
			pBestFitness = newFitness;
			pBestArray = position;
		}
	}

	/* Prints error if the user entered an invalid evaluation function. */
	else {
		cerr << "Optimization Function does not exist" << endl;
	}
}


/*
	Function that calculates a particle's fitness based on the Ackley evaluation function 
	and the user specified number of dimensions.
*/
double Particle::evalAckley () {

    double firstSum = 0.0;
    double secondSum = 0.0;
    double dimensions = double(position.size());

    for(int i = 0; i < dimensions; i++){
    	firstSum+= (position[i] * position[i]);
    }

    for(int i = 0; i < dimensions; i++){
    	secondSum += cos(2 * M_PI * position[i]);
    }
    

    return -20 * exp(-0.2 * sqrt(firstSum/dimensions)) - exp(secondSum/dimensions) + 20.0 + exp(1);
}  

/*
	Function that calculates a particle's fitness based on the Rosenbrock evaluation function
	and the user specified number of dimensions.
*/
double Particle::evalRosenbrock () {
	double sum = 0;
	for(int i = 1; i < position.size() -1; i++){
		sum += (100.0 * pow(position[i+1] -  position[i] * position[i], 2) + pow(position[i] - 1, 2));
	}
	return sum;
}

/*
	Function that calculates a particle's fitness based on the Rastrigin evaluation function
	and the user specified number of dimensions.
*/
double Particle::evalRastrigin () {

	double retVal = 0;

	for(int i = 0; i < position.size(); i++){
		retVal += (pow(position[i], 2) - 10* cos(2* M_PI * position[i]) + 10);
	}
    return retVal;
}


/*
	Function that updates the position of the particle in d dimensions based on the 
	velocity of the particle in each dimension.
*/
void Particle::updatePosition(){
	for(int i = 0; i < position.size(); i++) {
		position.at(i) = position.at(i) + velocity.at(i);
	}

}


/*
	Function that updates the velocity of the particle based on the standard PSO velocity equation
	found in the report. It uses the old velocity and biases from the pBest and nBest locations
	to determine the new velocity.
*/
void Particle::updateVelocity(){
	double newVelocity;
	double pBestBias;
	double nBestBias;


	random_device seeder;
	mt19937 engine(seeder());
	
	for(int i = 0; i < position.size(); i++) { 
		//generate random distribution for phi1 and phi2
		uniform_real_distribution<double> phi1(0,PHI_1);
		uniform_real_distribution<double> phi2(0,PHI_2);

		//calculate the pbestBias by using phi1 and nbestBias with phi2
		pBestBias = phi1(engine) * (pBestArray.at(i) - position.at(i));
		nBestBias = phi2(engine) * (nBestArray.at(i) - position.at(i));

		//calculate the new velocity using the constriction factor and old velocity
		newVelocity = CONSTRICTION_FACTOR * (velocity.at(i) + pBestBias + nBestBias);

		//check to make sure the new velocity isn't greater or less than the min/max values
		if(newVelocity < minVelocity){
			velocity.at(i) = minVelocity;
		}
		else if(newVelocity > maxVelocity){
			velocity.at(i) = maxVelocity;
		}
		else{
			velocity.at(i) = newVelocity;
		}

	}
}


/*
	Function that checks if any particle in the neighborhood has a better fitness than the 
	current nBest fitness. If so, the function updates the nBest fitness and location.
*/
void Particle::findNeighborhoodBest(){
	for(int i = 0; i < neighborsArray.size(); i++){
		if(neighborsArray[i]->pBestFitness < nBestFitness) {
			nBestArray = neighborsArray[i]->position;
			nBestFitness = neighborsArray[i]->pBestFitness;
		}
	}
}


/*
	Function that initializes the swarm object. The function creates the user specified
	number of particles and stores pointers to these particles in a vector. It then 
	calls the user specified topology to create neighborhoods for each of the particles.
*/
void Swarm::initSwarm(int swarmSize, int numDimensions, 
			string neighborhoodTopology, string testFunction){
	
	//set user specified swarm size and max out the global best fitness value
	this->swarmSize = swarmSize;
	gBestFitness = numeric_limits<double>::max();

	//loop that creates swarmSize particles
	for(int i = 0; i < swarmSize; i++){
		shared_ptr<Particle> ptr(new Particle());
		ptr->initParticle(numDimensions, testFunction);
		swarm.push_back(ptr);
	}

	/* Checks which neighborhood topology to apply. */
	if (neighborhoodTopology.compare(GLOBAL_TOPOLOGY) == 0){
		globalTopology();
	}

	else if (neighborhoodTopology.compare(RING_TOPOLOGY) == 0){
		ringTopology();
	}
	else if (neighborhoodTopology.compare(VON_NEUMANN_TOPOLOGY) == 0){
		vonNeumanTopology();
	}
	else if (neighborhoodTopology.compare(RANDOM_TOPOLOGY) == 0){
		randomTopology();
	}
	else {
		cerr << "Topology type not found" << endl;
	}


}

/*
	Function that checks if any particle in the swarm has a better fitness than the 
	current global best fitness. If so, the function updates the global best fitness 
	and location.
*/
void Swarm::findGlobalBest(){
	for (int i = 0; i < swarmSize; i++){
		if (swarm[i]->pBestFitness < gBestFitness){
			gBestArray = swarm[i]->pBestArray;
			gBestFitness = swarm[i]->pBestFitness;
		}
	}

}

/*
	Function that represents the global neighborhood topology. This is called when the swarm
	is initialized to specify the neighborhood of each particle. In the global topology the
	neighborhood of each particle is the entire swarm. 
*/
void Swarm::globalTopology(){
	for(int i = 0; i < swarmSize; i ++){
		for(int j = 0; j < swarmSize; j++){
			swarm[i]->neighborsArray.push_back(swarm.at(j));
		}
	}
}


/*
	Function that represents the ring neighborhood topology. Called by init swarm method.
	In this topology the neighbors of each particle are the particles to the left and right
	of the particle in the swarm vector (when it is initialized).
*/
void Swarm::ringTopology(){
	
	for(int i = 0; i < swarmSize; i++){
		
		//takes care of the first element in the swarm vector
		if (i == 0){
			swarm[i]->neighborsArray.push_back(swarm[swarmSize-1]);
			swarm[i]->neighborsArray.push_back(swarm[i]);
			swarm[i]->neighborsArray.push_back(swarm[i+1]);
		}

		//takes care of the last element in the swarm vector
		else if (i == swarmSize - 1){
			swarm[i]->neighborsArray.push_back(swarm[swarmSize-2]);
			swarm[i]->neighborsArray.push_back(swarm[swarmSize-1]);
			swarm[i]->neighborsArray.push_back(swarm[0]);
		}

		//takes care of all elements in the middle of the swarm vector
		else {
			swarm[i]->neighborsArray.push_back(swarm[i-1]);
			swarm[i]->neighborsArray.push_back(swarm[i]);
			swarm[i]->neighborsArray.push_back(swarm[i+1]);
		}
	}


}


/*
	Function that represents the von Neumann neighborhod topology. Called by the init swarm
	method. In this topology the neighbors of each particle are the particles to the left,
	right, top, and bottom of the particle in the swarm vector (when it is initialized).
	The easiest way to visualize and set the neighbors for this topology is to use
	a 2-d array. However, this is costly and it is possible to use a 1-d array.

*/
void Swarm::vonNeumanTopology(){

	/* 
		Determine possible number of rows and cols that the 1-d swarm array could be
		stored as a 2-d array. This is done to avoid the unnecessary creation of an
		actual 2-d array. 
	*/
	int swarmNumRows = 3;

	//finds a possible number of rows so that row*width = size of the array
	while(swarmSize % swarmNumRows != 0){
		swarmNumRows++;
	}

	int swarmNumCols = swarmSize / swarmNumRows;

	/* 
		Now that the 1-d array can be visualized as a 2-d array, the function finds
		the neighbors using logic.
	*/
	for(int i = 0; i < swarmSize; i++){
		/* Determine each particles left and right neighbors */

		//In case where particle index is in left column of imaginary 2-d array
		if(i % swarmNumCols == 0) {
			swarm[i]->neighborsArray.push_back(swarm[i+swarmNumCols-1]);
			swarm[i]->neighborsArray.push_back(swarm[i+1]);
		}

		//In case where particle index is in right column of imaginary 2-d array
		else if(i % swarmNumCols == swarmNumCols - 1) {
			swarm[i]->neighborsArray.push_back(swarm[i-1]);
			swarm[i]->neighborsArray.push_back(swarm[i-swarmNumCols+1]);
		}

		//In case where particle index is in any middle column of imaginary 2-d array
		else {
			swarm[i]->neighborsArray.push_back(swarm[i-1]);
			swarm[i]->neighborsArray.push_back(swarm[i+1]);
		}

		/* Determine each particles top and bottom neighbors */

		//In case where particle index is in top row of 2-d imaginary array
		if(i / swarmNumCols == 0) {
			swarm[i]->neighborsArray.push_back(swarm[i + ((swarmNumRows - 1) * swarmNumCols)]);
			swarm[i]->neighborsArray.push_back(swarm[i+swarmNumCols]);
		}

		/* In case where particle index is in bottom row of 2-d imaginary array */
		else if(i / swarmNumCols == swarmNumRows - 1) {
			swarm[i]->neighborsArray.push_back(swarm[i-swarmNumCols]);
			swarm[i]->neighborsArray.push_back(swarm[i - ((swarmNumRows - 1) * swarmNumCols)]);
		}

		/* In case where particle index is in any middle row of 2-d imaginary array */
		else {
			swarm[i]->neighborsArray.push_back(swarm[i-swarmNumCols]);
			swarm[i]->neighborsArray.push_back(swarm[i+swarmNumCols]);
		}

		//always push back the particle itself into its neighborhood
		swarm[i]->neighborsArray.push_back(swarm[i]);

	}

}


/*
	Function that represents the random neighborhood topology. Called by the init swarm
	method. The neighborhood of each particle is comprised of k-1 (k is constant set to 5)
	other particles in the swarm without duplicates. The neighborhood of each particle is
	recreated with a probability of 0.2 in each iteration.
*/
void Swarm::randomTopology(){

	random_device seeder;
	mt19937 engine(seeder());
	//gives a random index of a particle in the swarm
	uniform_int_distribution<int> randIndex(0, swarmSize-1);
	//gives a probability between 0 and 1 of neighborhood recreation
	uniform_real_distribution<double> randDouble(0, 1);

	//set is used to determine if particle is already in the neighborhood
	pair< set<int>::iterator, bool> inSet;
	set<int> used; 

	for(int i = 0; i < swarmSize; i++){
		/* 
			Neighborhood is recreated for a particle at the start and in later iterations
			if the probability is met.
		*/
		if(swarm[i]->neighborsArray.size() == 0 || (swarm[i]->neighborsArray.size() != 0 && 
			randDouble(engine) <= 0.2)) {

			//clear out the old neighborhood and add the particle itself to the neighborhood
			swarm[i]->neighborsArray.clear();
			swarm[i]->neighborsArray.push_back(swarm[i]);
			used.insert(i);

			//insert k-1 other random particles
			for(int j = 0; j < RANDOM_K-1; j++){
				
				int index = randIndex(engine);				
				inSet = used.insert(index);
				
				while(inSet.second == false){
					index = randIndex(engine);
					inSet = used.insert(index);
				}
				swarm[i]->neighborsArray.push_back(swarm[index]);
			}

			used.clear();
		}
	}
}

/*
	Function that runs particle swarm optimization on the given user inputted topology, 
	swarm size, iterations, test function, and dimensions. The function creates a swarm
	and then loops for the iteration length. During each loop it updates for each particle
	the neighborhood best, its velocity, position, and calculates its fitness.
*/
void PSO(string neighborhoodTopology, int swarmSize, int numIterations, 
	string testFunction, int numDimensions){
	
	//create and initialize the swarm
	shared_ptr<Swarm> swarmObject(new Swarm());
	swarmObject->initSwarm(swarmSize, numDimensions, neighborhoodTopology, testFunction);

	for(int i = 0; i < numIterations; i++ ){

		//if random topology is selected, recreate it with given probability
		if(neighborhoodTopology.compare(RANDOM_TOPOLOGY) == 0){
			swarmObject->randomTopology();
		}
		for(int j = 0; j < swarmSize; j++){
			swarmObject->swarm.at(j)->findNeighborhoodBest();
			swarmObject->swarm.at(j)->updateVelocity();
			swarmObject->swarm.at(j)->updatePosition();
			swarmObject->swarm.at(j)->calculateFitness(testFunction);
		}
		swarmObject->findGlobalBest();
		if(i % 1000 == 0) {
			cout << "Best Fitness on Iteration " << i << " is " << swarmObject->gBestFitness << endl;
		}
	}
	cout << "Best Overall Fitness found: " << swarmObject->gBestFitness << endl;
}


/*
	Main function for the PSO file. It reads in the user input and calls the PSO function.
*/
int main(int argc, char* argv[]){

	string neighborhoodTopology = string(argv[1]);
	int swarmSize = atoi(argv[2]);
	int numIterations = atoi(argv[3]);
	string testFunction = string(argv[4]);
	int numDimensions = atoi(argv[5]);
	cout << "Topology type: " << neighborhoodTopology << " swarmSize: " << swarmSize << " Number of Iterations: "  << numIterations 
	<< " testFunction: " << testFunction << " numDimensions: " << numDimensions << endl;
	
	PSO(neighborhoodTopology, swarmSize, numIterations, testFunction,numDimensions);


	return 1;

}

