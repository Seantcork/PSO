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

const double PHI_1 = 2.05;
const double PHI_2 = 2.05;

const string GLOBAL_TOPOLOGY = "gl";
const string RING_TOPOLOGY = "ri";
const string VON_NEUMANN_TOPOLOGY = "vn";
const string RANDOM_TOPOLOGY = "ra";

const string ROSENBROCK_FUNCTION = "rok";
const string ACKLEY_FUNCTION = "ack";
const string RASTRIGIN_FUNCTION = "ras";

/*
This class contains all of the necesary components for the Particle Class. This class keeps
track of the particles position, velocity, its personal best position, personal best fitness. It also contains the neighborhood best
for its array and the nbest position.
This class includes the fucntions: Calculate fitness, which calculates the fitness of the position of the particle
with respect to the evaluation function. It also contains the reelvant functions desired for position.


*/
class Particle {
	public:

		double maxVelocity;
		double minVelocity;

		double pBestFitness;
		vector<double> position;
		vector<double> velocity;
		vector<double> pBestArray;
		
		double nBestFitness;
		vector<double> nBestArray;
		vector<shared_ptr<Particle> > neighborsArray;

		void calculateFitness(string testFunction);
		double evalAckley ();
		double evalRosenbrock ();
		double evalRastrigin ();

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
		maxVelocity = 2.048;
		minVelocity = -2.048;
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
		this->maxVelocity = 32.768;
		this->minVelocity = -32.768;
		for(int i = 0; i < numDimensions; i ++){
			this->position.push_back(genPosition(engine));
			this->velocity.push_back(genVelocity(engine));

		}
	}
	
	else if(testFunction.compare(RASTRIGIN_FUNCTION) == 0){
		uniform_real_distribution<double> genPosition(2.56, 5.12);
		uniform_real_distribution<double> genVelocity(-2.0, 4.0);
		this->maxVelocity = 5.12;
		this->minVelocity = -5.12;
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
	double currFitness = 0;
	// Determine which test function to run
	// Evaluate the fitness and update the values if its better than pBest or nBest
	if(testFunction.compare(ROSENBROCK_FUNCTION) == 0) {
		currFitness = evalRosenbrock(this->position);
		if(currFitness < this->pBestFitness){
			this->pBestFitness = currFitness;
			this->pBestArray = this->position;
		}

	}

	else if (testFunction.compare(ACKLEY_FUNCTION) == 0) {
		currFitness = evalAckley(position);
		if(currFitness < this->pBestFitness){
			this->pBestFitness = currFitness;
			this->pBestArray = this->position;
		}
	}

	else if (testFunction.compare(RASTRIGIN_FUNCTION) == 0){
		currFitness = evalRastrigin(position);
		if(currFitness < this->pBestFitness) {
			this->pBestFitness = currFitness;
			this->pBestArray = this->position;
		}
	}
	else {
		cerr << "Optimization Function does not exist" << endl;
	}
	// cout << "currFitness = " << currFitness << endl;

}

double Particle::evalAckley () {

    double firstSum = 0.0;
    double secondSum = 0.0;
    double dimensions = position.size();

    for(int i = 0; i < position.size(); i++){
    	firstSum+= (position[i] * position[i]);
    }

    for(int i = 0; i < position.size(); i++){
    	secondSum += cos(2 * M_PI * position[i]);
    }
    

    return -20 * exp(-0.2 * sqrt(firstSum/dimensions)) - exp(secondSum/dimensions) + 20.0 + exp(1);
}  

 //evaluates rosenbrock for the specified number of dimensions
double Particle::evalRosenbrock () {
	double sum = 0;
	for(int i = 1; i < position.size() -1; i++){
		sum += (100.0 * pow(position[i+1] -  position[i] * position[i], 2) + pow(position[i] - 1, 2));
	}
	return sum;
}

 // returns the value of the Rastrigin Function at point (x, y)
 // minimum is 0.0, which occurs at (0.0,...,0.0)
double Particle::evalRastrigin () {

	double retVal = 0;

	for(int i = 0; i < position.size(); i++){
		retVal += (pow(position[i], 2) - 10* cos(2* M_PI * position[i]) + 10);
	}
    return retVal;
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

		if(newVelocity < minVelocity){
			this->velocity.at(i) = minVelocity;
		}
		else if(newVelocity > maxVelocity){
			this->velocity.at(i) = maxVelocity;
		}
		else{
			this->velocity.at(i) = newVelocity;
		}

	}
}


//Purpose: finds the neighboorhood best
//Parameters: None
//Return value none
void Particle::findNeighborhoodBest(){
	for(int i = 0; i < neighborsArray.size(); i++){
		if(this->neighborsArray[i]->pBestFitness < this->nBestFitness) {
			this->nBestArray = this->neighborsArray[i]->position;
			this->nBestFitness = this->neighborsArray[i]->pBestFitness;
		}
	}
}


// //Purpose: if a new pbest has been found makes sure its the neighborhood best
// //Parameters: the new pbest, and the new array signifiing its position
// //Return value: None
// void Particle::updateNeighborhoodBest(double bestFitness, vector<double> bestFitArray) { 

// 	for(int i = 0; i < neighborsArray.size(); i++ ) {
// 		this->neighborsArray.at(i)->nBestArray = bestFitArray;
// 		this->neighborsArray.at(i)->nBestFitness = bestFitness;
// 	}

// }


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
	for (int i = 0; i < swarmSize; i++){
		if (swarm[i]->pBestFitness < gBestFitness){
			this->gBestArray = swarm[i]->pBestArray;
			this->gBestFitness = swarm[i]->pBestFitness;
		}
	}

}

//not sure we need this
void Swarm::globalTopology(){
	for(int i = 0; i < swarmSize; i ++){
		for(int j = 0; j < swarmSize; j++){
			this->swarm[i]->neighborsArray.push_back(swarm.at(j));
		}
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

void PSO(string neighborhoodTopology, int swarmSize, int numIterations, string testFunction, int numDimensions){
	
	shared_ptr<Particle> swarmObject(new Swarm());

	swarmObject->initSwarm(swarmSize, numDimensions, neighborhoodTopology, testFunction);

	for(int i = 0; i < numIterations; i++ ){
		for(int j = 0; j < swarm.swarmSize; j++){
			swarmObject->swarm.at(j)->updateVelocity();
			swarmObject->swarm.at(j)->updatePosition();
			swarmObject->swarm.at(j)->calculateFitness(testFunction);
			swarmObject->swarm.at(j)->findNeighborhoodBest();
		}
		swarmObject->findGlobalBest();
	}
	cout << "Best Fitness found: " << swarmObject->gBestFitness << endl;

}



int main(int argc, char* argv[]){

	string neighborhoodTopology = string(argv[1]);
	int swarmSize = atoi(argv[2]);
	int numIterations = atoi(argv[3]);
	string testFunction = string(argv[4]);
	int numDimensions = atoi(argv[5]);

	
	PSO(neighborhoodTopology, swarmSize, numIterations, testFunction,numDimensions);


	return 1;

}

