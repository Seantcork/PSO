/*Authors Sean Cork, Kamaal Palmer, Luca Ostertag-Hill*/
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <stdint.h>
#include <random>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

# define M_PI    3.14159265358979323846

#define E 2.71828
const double CONSTRICTION_FACTOR = 0.7298;

const double P_BEST_ATTRACTION = 0.5;

const double G_BEST_ATTRACTION = 0.5;

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



class Particle {
	public:
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
void Particle::initParticle(int numDimensions, string testFunction){
	random_device seeder;
	mt19937 engine(seeder());

	//initialize particle positions
	if(testFunction.compare(ROSENBROCK_FUNCTION) == 0){
		uniform_real_distribution<double> genPosition(15.0, 30.0);
		uniform_real_distribution<double> genVelocity(-2.0, 2.0);
		for(int i = 0; i < numDimensions; i ++){
			position.push_back(genPosition(engine));
			velocity.push_back(genVelocity(engine));
		}
	}

	else if(testFunction.compare(ACKLEY_FUNCTION) == 0){
		uniform_real_distribution<double> genPosition(16.0, 32.0);
		uniform_real_distribution<double> genVelocity(-2.0, 4.0);
		for(int i = 0; i < numDimensions; i ++){
			position.push_back(genPosition(engine));
			velocity.push_back(genVelocity(engine));

		}
	}
	
	else if(testFunction.compare(RASTRIGIN_FUNCTION) == 0){
		uniform_real_distribution<double> genPosition(2.56, 5.12);
		uniform_real_distribution<double> genVelocity(-2.0, 4.0);
		for(int i = 0; i < numDimensions; i ++){
			position.push_back(genPosition(engine));
			velocity.push_back(genVelocity(engine));

		}
	}

	else {
		cerr << "Optimization Function does not exist" << endl;
	}

	//Particles pBest is set as its initial position
	pBestArray = position;

	
}

void Particle::calculateFitness(string testFunction){
	double currFitness = 0;
	// Determine which test function to run
	// Evaluate the fitness and update the values if its better than pBest or nBest
	if(testFunction.compare(ROSENBROCK_FUNCTION) == 0) {
		currFitness = evalRosenbrock(position);
		if(currFitness > pBestFitness){
			pBestFitness = currFitness;
			pBestArray = position;
		}
		if(currFitness > nBestFitness) {
			nBestFitness = currFitness;
			nBestArray = position;
			updateNeighborhoodBest(currFitness, position);
		}
	}

	else if (testFunction.compare(ACKLEY_FUNCTION) == 0) {
		currFitness = evalAckley(position);
		if(currFitness > pBestFitness){
			pBestFitness = currFitness;
			pBestArray = position;
		}
		if(currFitness > nBestFitness) {
			nBestFitness = currFitness;
			nBestArray = position;
			updateNeighborhoodBest(currFitness, position);
		}
	}

	else if (testFunction.compare(RASTRIGIN_FUNCTION) == 0){
		currFitness = evalRastrigin(position);
		if(currFitness > pBestFitness) {
			pBestFitness = currFitness;
			pBestArray = position;
		}
		if(currFitness > nBestFitness) {
			nBestFitness = currFitness;
			nBestArray = position;
			updateNeighborhoodBest(currFitness, position);
		}
	}
	else {
		cerr << "Optimization Function does not exist" << endl;
	}
}

void Particle::updatePosition(){
	for(int i = 0; i < position.size(); i++) {
		position.at(i) = position.at(i) + velocity.at(i);
	}

}


// Possible here that we want a max velocity -- Discuss at next group meeting
void Particle::updateVelocity(){
	random_device seeder;
	mt19937 engine(seeder());
	uniform_real_distribution<double> randAcceleration(0.0, 3.0);
	for(int i = 0; i < position.size(); i++) {
		velocity.at(i) = velocity.at(i) + 
		(randAcceleration(engine) * (pBestArray.at(i) - position.at(i))) + 
		(randAcceleration(engine) * (nBestArray.at(i) - position.at(i)));
	}
}

void Particle::findNeighborhoodBest(){
	for(int i = 0; i < neighborsArray.size(); i++){
		if(neighborsArray[i]->pBestFitness < nBestFitness) {
			nBestArray = neighborsArray[i]->position;
			nBestFitness = neighborsArray[i]->pBestFitness;
		}
	}
}

void Particle::updateNeighborhoodBest(double bestFitness, vector<double> bestFitArray) { 

	for(int i = 0; i < neighborsArray.size(); i++ ) {
		neighborsArray.at(i)->nBestArray = bestFitArray;
		neighborsArray.at(i)->nBestFitness = bestFitness;
	}

}


/*
	Functions for swarm class
*/
void Swarm::initSwarm(int swarmSize, int numDimensions, 
			string neighborhoodTopology, string testFunction){
	this->swarmSize = swarmSize;
	gBestFitness = DBL_MAX;
	for(int i = 0; i < swarmSize; i++){
		shared_ptr<Particle> ptr(new Particle());
		ptr->initParticle(numDimensions, testFunction);
		swarm.push_back(ptr);
	}
}

void Swarm::findGlobalBest(){
	for (int i = 0; i < swarmSize; i++){
		if (swarm[i]->pBestFitness < gBestFitness){
			gBestArray = swarm[i]->pBestArray;
			gBestFitness = swarm[i]->pBestFitness;
		}
	}

}

//not sure we need this
void Swarm::globalTopology(){
	for(int i = 0; i < swarm(); i ++){
		for(int j = 0; j < swarm.size(); j++){
			swarm[i].neighbors.push_back(swarm[j]);

		}
		swarm[i].neighbors.push_back(swarm[i]);
	}
}


//This function bases neighborhoods on the 
void Swarm::ringTopology(){
	
	//takes care of all the elemetns between the first and last element
	for(int i = 0; i < swarmSize; i++){
		
		//takes care of first elements
		if (i == 0){
			swarm[i].neighbors.push_back(swarm[swarmSize-1]);
			swarm[i].neighbors.push_back(swarm[i]);
			swarm[i].neighbors.push_back(swarm[i+1])
		}

		//takes care of last element in the swarm
		else if (i == swarmSize - 1){
			swarm[i].neighbors.push_back(swarm[swarm.size()-2]);
			swarm[i].neighbors.push_back(swarm[swarm.size()-1]);
			swarm[i].neighbors.push_back(swarm[0]);
		}

		else {
			swarm[i].neighbors.push_back(swarm[i-1]);
			swarm[i].neighbors.push_back(swarm[i]);
			swarm[i].neighbors.push_back(swarm[i+1]);
		}
	}


}

void Swarm::vonNeumanTopology(){

}

void Swarm::randomTopology(){
	std::random_device seeder;
	std::mt19937 engine(seeder());
	std::uniform_int_distribution<int> gen(1, rankSum);


}



//returns distance between particles
double distance(Particle a, Particle b){

}


double evalAckley (vector<double> positions) {

    double firstSum = 0;
    double secondSum = 0;

    for(int i = 0; i < positions.size(); i++){
    	firstSum += positions[i];
    }
    firstSum = -0.2 * sqrt((1/positions.size()) * firstSum);

    for(int i = 0; i < positions.size(); i++){
    	secondSum += cos(2 * M_PI * positions[i]);
    }

    secondSum = exp(secondSum) + 20 + E;

    return -20 * firstSum - secondSum;
}  



 //evaluates rosenbrock for the specified number of dimensions
double evalRosenbrock (vector<double> position) {
	double sum = 0;

	for(int i = 0; i < position.size() -1; i++){
		sum+ = 100.0 * pow(position[i+1] - pow(position[i], 2), 2) + pow(position[i] -1, 2);
	}

	return sum;
}



 // returns the value of the Rastrigin Function at point (x, y)
 //   minimum is 0.0, which occurs at (0.0,...,0.0)
double evalRastrigin (vector<double> position) {

	double retVal = 0;

	for(int i = 0; i < position.size(); i++){
		retval += (pow(position[i], 2) - 10* cos(2* M_PI * position[i]) + 10);
	}
    return retVal;
}




void PSO(){

}



int main(int argc, char* argv[]){

	string neighborhoodTopology = string(argv[1]);
	int swarmSize = atoi(argv[2]);
	int numIterations = atoi(argv[3]);
	string testFunction = string(argv[4]);
	int numDimensions = atoi(argv[5]);

	Swarm swarm = new Swarm;
	swarm->initSwarm(neighborhoodTopology, swarmSize, testFunction, numDimensions);

	PSO(swarm)

	 	return 1;



}

