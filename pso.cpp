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

const string GLOBAL_TOPOLOGY = "gl";
const string RING_TOPOLOGY = "ri";
const string VON_NEUMANN_TOPOLOGY = "vn";
const string RANDOM_TOPOLOGY = "ra";

const string ROSENBROCK_FUNCTION = "rok";
const string ACKLEY_FUNCTION = "ack";
const string RASTRIGIN_FUNCTION = "ras";


class Particle {
	public:
		vector<double> position;
		vector<double> pBest;
		vector<double> velocity;
		
		vector<Particle*> neighbors;
		double neighborhoodBest;

		double fitness;
		void calculateFitness();

		void updatePosition();
		void updateVelocity();
		void findNeighborhoodBest();

		void initParticle(int numDimensions, string testFunction);	
}


class Swarm {

	public:
		vector<shared_ptr<Particle> > swarm;
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
void particle::findNeighborhoodBest(){
	vector<double> minPos;
	double minFitness = 1000000;

	for(int i = 0; i < this->neighbors.size(); i++){
		if(this->neighbors[i].fitness < minFitness) {
			minPos = this->neighbors[i].position;
			minFitness = this->neighbors[i].fitness;
		}
	}
}

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
			this->position.push_back(genPosition(engine));
			this->velocity.push_back(genVelocity(engine));
		}
	}

	if(testFunction.compare(ACKLEY_FUNCTION)){
		uniform_real_distribution<double> genPosition(16.0, 32.0);
		uniform_real_distribution<double> genVelocity(-2.0, 4.0);
		for(int i = 0; i < numDimensions; i ++){
			this->position.push_back(genPosition(engine));
			this->velocity.push_back(genVelocity(engine));

		}
	}
	
	if(testFunction.compare(RASTRIGIN_FUNCTION)){
		uniform_real_distribution<double> genPosition(2.56, 5.12);
		uniform_real_distribution<double> genVelocity(-2.0, 4.0);
		for(int i = 0; i < numDimensions; i ++){
			this->position.push_back(genPosition(engine));
			this->velocity.push_back(genVelocity(engine));

		}
	}
	cerr << "Optimization Function does not exist" << endl;
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
	for(int i = 0; i < swarmSize; i++){
		shared_ptr<Particle> ptr(new Particle());
		ptr->initParticle(numDimensions, testFunction);
		swarm.push_back(ptr);
	}
}

void Swarm::findGlobalBest(){
	for (int i = 0; i < swarmSize; i++){
		if (swarm[i]->fitness < globalFitness){
			gBest = swarm[i]->pBest;
			globalFitness = swarm[i]->fitness;
		}
	}

}

void Swarm::globalTopology(vector<particle> swarm){
	for(int i = 0; i < swarm(); i ++){
		for(int j = 0; j < swarm.size(); j++){
			swarm[i].neighbors.push_back(swarm[j]);

		}
		swarm[i].neighbors.push_back(swarm[i]);
	}
}


//This function bases neighborhoods on the 
void Swarm::ringTopology(vector<particle> swarm){
	//takes care of first elements
	swarm[0].neighbors.push_back(swarm[swarm.size()]);
	swarm[0].neighbors.push_back(swarm[1]);
	swarm[0].neighbors.push_back(swarm[0])
	//takes care of all the elemetns between the first and last element
	for(int i = 1; i < swarm.size(); i ++){
		swarm[i].neighbors.push_back(swarm[i+1]);
		swarm[i].neighbors.push_back(swarm[i-1]);
		swarm[i].push_back(swarm[i])
	}

	//takes care of last element in the swarm
	swarm[swarm.size()].neighbors.push_back(swarm[swarm.size()]);
	swarm[swarm.size()].neighbors.push_back(swarm[swarm.size()-1]);
	swarm[swarm.size()].neighbors.push_back(swarm[0]);


}

void Swarm::vonNeumanTopology(){

}

void Swarm::randomTopology(){
	std::random_device seeder;
	std::mt19937 engine(seeder());
	std::uniform_int_distribution<int> gen(1, rankSum);


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

	 


}

	return 1;
}