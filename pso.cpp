// Authors Sean Cork, Kamaal Palmer, Luca Ostertag-Hill
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <stdint.h>
#include <random>
#include <string>
#include <vector>
#include <iostream>

# define M_PI    3.14159265358979323846



const double CONSTIRCTION_FACTOR = 0.7298;


class neighborhood{
	public:
		//update neighbors each time velocity update
		vector<*particle> neighbors;
		double neighborhoodBest;
		findNeighborhoodBest();
		particle Ring();
		particle vonNeuman();
}



class particle {
	public:
		vector<double> position;
		vector<double> pBest;
		vector<double> velocity;
		double fitness;
		void calculateFitness();
		void updateVelocity();
		void updatePosition();


}


class swarm {

	public:
		vector<neighborhood> swarm;
		double global bestl;

}

public initializeParticle(particle a){

}

//parses comand line and makes sure to enter the right things
void evaluate(){

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

    return -20.0 * exp(-0.2 * sqrt(firstSum/2.0)) - 
      exp(secondSum/2.0) + 20.0 + Math.E;
}  



 //evaluates rosenbrock for the specified number of dimenstions
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


void global(){

}

void ring(){

}

void vonNeuman(){

}

void random(){


}


void Pso(){

}



int main(int argc, char* argv[]){

	return 1;
}