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

class Particle;

class Neighborhood {
	public:
		//update neighbors each time velocity update
		vector<Particle*> neighbors;
		double neighborhoodBest;
		void findNeighborhoodBest();
}


class Particle {
	public:
		vector<double> position;
		vector<double> pBest;
		vector<double> velocity;
		Neighborhood neighbors;

		double fitness;
		void calculateFitness();

		void updatePosition();
		void updateVelocity();

		void initParticle(int numDimensions, string testFunction);	
}


class Swarm {

	public:
		vector<*Particle> swarm;
		double swarmSize;


		vector<double> gBest;
		double globalFitness;
		void findGlobalBest();

		void globalTopology();
		void ringTopology();
		void vonNeumanTopology();
		void randomTopology();

		void initSwarm(int swarmSize, int numDimensions, 
			string neighborhoodTopology, string testFunction);
}


/*
	Functions for neighborhood class
*/
void Neighborhood::findNeighborhoodBest(){

}

/*
	Functions for particle class
*/
void Particle::initParticle(int numDimensions, string testFunction){

}

void Particle::calculateFitness(){

}

void Particle::updatePosition(){

}

void Particle::updateVelocity(){

}


/*
	Functions for swarm class
*/
void Swarm::initSwarm(int swarmSize, int numDimensions, 
			string neighborhoodTopology, string testFunction){
	this->swarmSize = swarmSize;
	this->globalFitness = DBL_MAX;
}

void Swarm::findGlobalBest(){
	for (int i = 0; i < swarmSize; i++){
		if (swarm[i]->fitness > globalFitness){
			gBest = swarm[i]->pBest;
			globalFitness = swarm[i]->fitness;
		}
	}

}

void Swarm::globalTopology(){

}

void Swarm::ringTopology(){

}

void Swarm::vonNeumanTopology(){

}

void Swarm::randomTopology(){


}

void Swarm::initSwarm(){

}



//returns distance between particles
double distance(particle a, particle b){

}


double evalAckley (vector<double> positions) {

    double firstSum = 0
    double secondSum = 0

    for(int i = 0; i < position.size(); i++){
    	firstSum += position[i];
    }
    firstSum = -0.2 * sqrt((1/positions.size()) * firstSum);

    for(int i = 0; i < position.size(); i++){
    	secondSum += cos(2 * M_PI * position[i])
    }

    secondSum = exp(secondSum) + 20 + E;

    return -20 * firstSum - secondSum;
}  



 //evaluates rosenbrock for the specified number of dimensions
public double evalRosenbrock (vector<double> position) {
	double sum = 0;

	for(int i = 0; i < position.size() -1; i++){
		sum+ = 100.0 * pow(position[i+1] - pow(position[i], 2), 2) + pow(position[i] -1, 2);
	}

	return sum;
}




  // returns the value of the Rastrigin Function at point (x, y)
  //   minimum is 0.0, which occurs at (0.0,...,0.0)
public double evalRastrigin (vector<double> position) {

	double retVal = 0;

	for(int i = 0; i < position.size(); i++){
		retval += (pow(position[i], 2) - 10* cos(2* M_PI * position[i]) + 10)
	}
    return retVal;
}




void PSO(Swarm swarm){

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

	 


}

	return 1;
}