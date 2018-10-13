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
#include <limits>       // std::numeric_limits

using namespace std;


#define E 2.71828
const double CONSTRICTION_FACTOR = 0.7298;

const double C1 = 2.05;

const double C2 = 2.05;

const string GLOBAL_TOPOLOGY = "gl";
const string RING_TOPOLOGY = "ri";
const string VON_NEUMANN_TOPOLOGY = "vn";
const string RANDOM_TOPOLOGY = "ra";

const string ROSENBROCK_FUNCTION = "rok";
const string ACKLEY_FUNCTION = "ack";
const string RASTRIGIN_FUNCTION = "ras";
double evalAckley (vector<double> positions);
double evalRosenbrock (vector<double> position);
double evalRastrigin (vector<double> position);


/*
This class contains all of the necesary components for the Particle Class. This class keeps
track of the particles position, velocity, its personal best position, personal best fitness. It also contains the neighborhood best
for its array and the nbest position.
This class includes the fucntions: Calculate fitness, which calculates the fitness of the position of the particle
with respect to the evaluation function. It also contains the reelvant functions desired for position.


*/
class Particle {
	public:

		double MAX_VELOCITY;
		double MIN_VELOCITY;

		double pBestFitness;
		vector<double> position;
		vector<double> velocity;
		vector<double> pBestArray;
		
		double nBestFitness;
		vector<double> nBestArray;
		vector<shared_ptr<Particle> > neighborsArray;

		void calculateFitness(string testFunction);

		void updatePosition();
		void updateVelocity();

		void findNeighborhoodBest();
		void updateNeighborhoodBest(double currFitness, vector<double> bestFitArray);

		void initParticle(int numDimensions, string testFunction);	
};


class Swarm {

	public:
		vector<shared_ptr<Particle> > swarm;
		
		double swarmSize;
		int numDimensions;

		vector<double> gBestArray;
		double gBestFitness;
		void findGlobalBest();

		void globalTopology();
		void ringTopology();
		void vonNeumanTopology();
		void randomTopology();

		void initSwarm(int swarmSize, int numDimensions, 
			string neighborhoodTopology, string testFunction);
};


/*
	Functions for particle class
*/

//Purpose: initialized the particles position and velocity for the number of dimensions for each given
//Evaluation function
//Parameters: Number of dimensions, the testFunctions.
//return value: None
void Particle::initParticle(int numDimensions, string testFunction){
	random_device seeder;
	mt19937 engine(seeder());

	//initialize particle positions
	//Check for each test function
	if(testFunction.compare(ROSENBROCK_FUNCTION) == 0){
		this->MAX_VELOCITY = 2.048;
		this->MIN_VELOCITY = -2.048;
		//random number generations for each different test fucntions
		uniform_real_distribution<double> genPosition(15.0, 30.0);
		uniform_real_distribution<double> genVelocity(-2.0, 2.0);
		for(int i = 0; i < numDimensions; i ++){
			this->position.push_back(genPosition(engine));
			this->velocity.push_back(genVelocity(engine));
		}
	}

	else if(testFunction.compare(ACKLEY_FUNCTION) == 0){
		uniform_real_distribution<double> genPosition(16.0, 32.0);
		uniform_real_distribution<double> genVelocity(-2.0, 4.0);
		this->MAX_VELOCITY = 32.768;
		this->MIN_VELOCITY = -32.768;
		for(int i = 0; i < numDimensions; i ++){
			this->position.push_back(genPosition(engine));
			this->velocity.push_back(genVelocity(engine));

		}
	}
	
	else if(testFunction.compare(RASTRIGIN_FUNCTION) == 0){
		uniform_real_distribution<double> genPosition(2.56, 5.12);
		uniform_real_distribution<double> genVelocity(-2.0, 4.0);
		this->MAX_VELOCITY = 5.12;
		this->MIN_VELOCITY = -5.12;
		for(int i = 0; i < numDimensions; i ++){
			this->position.push_back(genPosition(engine));
			this->velocity.push_back(genVelocity(engine));

		}
	}

	//if the string is not a function we are evaluating
	else {

		cerr << "Optimization Function does not exist" << endl;
	}


	this->pBestFitness = numeric_limits<double>::max();
	this->nBestFitness = numeric_limits<double>::max();
	//Particles pBest is set as its initial position
	this->pBestArray = this->position;

	//the nBestArray is also at its initial position
	this->nBestArray = this->position;

	
}


/*
Purpose: Calculates the fitness of each individual particle with respect to the testFunction specified.
Parameters: the testFunction
Return value: none

*/
void Particle::calculateFitness(string testFunction){
	double currFitness = numeric_limits<double>::max();
	// Determine which test function to run
	// Evaluate the fitness and update the values if its better than pBest or nBest
	if(testFunction.compare(ROSENBROCK_FUNCTION) == 0) {
		currFitness = evalRosenbrock(this->position);
		if(currFitness < this->pBestFitness){
			this->pBestFitness = currFitness;
			this->pBestArray = this->position;
		}

		//if the current best is better than the current known neighborhood best
		//update it
		if(currFitness < this->nBestFitness) {
			this->nBestFitness = currFitness;
			this->nBestArray = this->position;
			updateNeighborhoodBest(currFitness, this->position);
		}
	}

	else if (testFunction.compare(ACKLEY_FUNCTION) == 0) {
		currFitness = evalAckley(position);
		if(currFitness < this->pBestFitness){
			this->pBestFitness = currFitness;
			this->pBestArray = this->position;
		}
		if(currFitness < this->nBestFitness) {
			this->nBestFitness = currFitness;
			this->nBestArray = this->position;
			updateNeighborhoodBest(currFitness, this->position);
		}
	}

	else if (testFunction.compare(RASTRIGIN_FUNCTION) == 0){
		currFitness = evalRastrigin(position);
		if(currFitness < this->pBestFitness) {
			this->pBestFitness = currFitness;
			this->pBestArray = this->position;
		}
		if(currFitness < this->nBestFitness) {
			this->nBestFitness = currFitness;
			this->nBestArray = this->position;
			updateNeighborhoodBest(currFitness, this->position);
		}
	}
	else {
		cerr << "Optimization Function does not exist" << endl;
	}
	// cout << "currFitness = " << currFitness << endl;

}


//Parameters: None
//Purpose: updates the position with respect to the current velocity
//Return value: 
void Particle::updatePosition(){
	for(int i = 0; i < this->position.size(); i++) {
		// cout << "Position changed to: " << this->position.at(i) + this->velocity.at(i) << endl;
		this->position.at(i) = this->position.at(i) + this->velocity.at(i);
	}

}


//Purpose: update the velcoity of each particle in respect to its initial values
//Parameters: none
//Return value: none
void Particle::updateVelocity(){
	double currVelocity;
	random_device seeder;
	mt19937 engine(seeder());
	for(int i = 0; i < position.size(); i++) { 
		//for E element 
		// This value is really important to getting good values on our runs
		uniform_real_distribution<double> randAcceleration(0,2.0);
		currVelocity = CONSTRICTION_FACTOR * (this->velocity.at(i) +
		(C1*randAcceleration(engine) * (pBestArray.at(i) - this->position.at(i))) + 
		((C2*randAcceleration(engine)) * (nBestArray.at(i) - this->position.at(i)))); 

		if(currVelocity < MIN_VELOCITY){
			this->velocity.at(i) = MIN_VELOCITY;
		}
		else if(currVelocity > MAX_VELOCITY){
			this->velocity.at(i) = MAX_VELOCITY;
		}
		else{
			this->velocity.at(i) = currVelocity;
		}

	}
}


//Purpose: finds the neighboorhood best
//Parameters: None
//Return value none
void Particle::findNeighborhoodBest(){

	vector<double> bestArray;
	double bestFitness = numeric_limits<double>::max();
	for(int i = 0; i < neighborsArray.size(); i++){
		if(this->neighborsArray[i]->pBestFitness < this->nBestFitness) {
			this->nBestArray = this->neighborsArray[i]->position;
			this->nBestFitness = this->neighborsArray[i]->pBestFitness;
		}
	}
	for(int i = 0; i < neighborsArray.size(); i++){
		this->nBestArray = bestArray;
		this->nBestFitness = bestFitness;
	}
}


//Purpose: if a new pbest has been found makes sure its the neighborhood best
//Parameters: the new pbest, and the new array signifiing its position
//Return value: None
void Particle::updateNeighborhoodBest(double bestFitness, vector<double> bestFitArray) { 

	for(int i = 0; i < neighborsArray.size(); i++ ) {
		this->neighborsArray.at(i)->nBestArray = bestFitArray;
		this->neighborsArray.at(i)->nBestFitness = bestFitness;
	}

}


/*
	Functions for swarm class
*/
void Swarm::initSwarm(int swarmSize, int numDimensions, 
			string neighborhoodTopology, string testFunction){
	this->swarmSize = swarmSize;
	this->gBestFitness = numeric_limits<double>::max();
	for(int i = 0; i < swarmSize; i++){
		shared_ptr<Particle> ptr(new Particle());
		ptr->initParticle(numDimensions, testFunction);

		//sets the min and max value for each particle
		this->swarm.push_back(ptr);
	}
}

void Swarm::findGlobalBest(){
	for (int i = 0; i < swarm.size(); i++){
		if (swarm[i]->pBestFitness < gBestFitness){
			this->gBestArray = swarm[i]->pBestArray;
			this->gBestFitness = swarm[i]->pBestFitness;
		}
	}

}

//not sure we need this
void Swarm::globalTopology(){
	for(int i = 0; i < swarm.size(); i ++){
		for(int j = 0; j < swarm.size(); j++){
			this->swarm[i]->neighborsArray.push_back(swarm.at(j));

		}
		swarm[i]->neighborsArray.push_back(swarm.at(i));
	}
}


//This function bases neighborhoods on the 
void Swarm::ringTopology(){
	
	//takes care of all the elemetns between the first and last element
	for(int i = 0; i < swarmSize; i++){
		
		//takes care of first elements
		if (i == 0){
			swarm[i]->neighborsArray.push_back(swarm[swarmSize-1]);
			swarm[i]->neighborsArray.push_back(swarm[i]);
			swarm[i]->neighborsArray.push_back(swarm[i+1]);
		}

		//takes care of last element in the swarm
		else if (i == swarmSize - 1){
			swarm[i]->neighborsArray.push_back(swarm[swarm.size()-2]);
			swarm[i]->neighborsArray.push_back(swarm[swarm.size()-1]);
			swarm[i]->neighborsArray.push_back(swarm[0]);
		}

		else {
			swarm[i]->neighborsArray.push_back(swarm[i-1]);
			swarm[i]->neighborsArray.push_back(swarm[i]);
			swarm[i]->neighborsArray.push_back(swarm[i+1]);
		}
	}


}

void Swarm::vonNeumanTopology(){

}

void Swarm::randomTopology(){
	std::random_device seeder;
	std::mt19937 engine(seeder());


}

//returns distance between particles
double distance(Particle a, Particle b){
	return 0.0;
}


double evalAckley (vector<double> positions) {

    double firstSum = 0.0;
    double secondSum = 0.0;

    for(int i = 0; i < positions.size(); i++){
    	firstSum += pow(positions[i], 2);
    }

    for(int i = 0; i < positions.size(); i++){
    	secondSum += cos(2 * M_PI * positions[i]);
    }

    firstSum = -20.0 * exp(-0.2 * sqrt(firstSum/positions.size()));

    secondSum = -exp(secondSum/positions.size());

    return( firstSum + secondSum + 20 + exp(1));
}  

 //evaluates rosenbrock for the specified number of dimensions
double evalRosenbrock (vector<double> position) {
	double sum = 0;
	for(int i = 1; i < position.size() -1; i++){
		sum += (100.0 * pow(position[i+1] -  position[i] * position[i], 2) + pow(position[i] - 1, 2));
	}
	return sum;
}

 // returns the value of the Rastrigin Function at point (x, y)
 // minimum is 0.0, which occurs at (0.0,...,0.0)
double evalRastrigin (vector<double> position) {

	double retVal = 0;

	for(int i = 0; i < position.size(); i++){
		retVal += (pow(position[i], 2) - 10* cos(2* M_PI * position[i]) + 10);
	}
    return retVal;
}

void PSO(Swarm swarm, int numIterations, string testFunction){
	for(int i = 0; i < numIterations; i++ ){
		for(int j = 0; j < swarm.swarmSize; j++){
			swarm.swarm.at(j)->updateVelocity();
			swarm.swarm.at(j)->updatePosition();
			swarm.swarm.at(j)->calculateFitness(testFunction);
			swarm.swarm.at(j)->updateVelocity();
			swarm.findGlobalBest();
		}
	}
	cout << "Best Fitness found: " << swarm.gBestFitness << endl;

}



int main(int argc, char* argv[]){

	string neighborhoodTopology = string(argv[1]);
	int swarmSize = atoi(argv[2]);
	int numIterations = atoi(argv[3]);
	string testFunction = string(argv[4]);
	int numDimensions = atoi(argv[5]);

	Swarm *swarm = new Swarm;
	swarm->initSwarm(swarmSize, numDimensions, neighborhoodTopology, testFunction);
	PSO(*swarm, numIterations, testFunction);


	return 1;



}

