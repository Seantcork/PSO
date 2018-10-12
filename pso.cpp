// Authors Sean Cork, Kamaal Palmer, Luca Ostertag-Hill
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <stdint.h>
#include <random>
#include <string>
#include <vector>
#include <iostream>

# define M_PI           3.14159265358979323846



const double CONSTIRCTION_FACTOR = 0.7298;


class neighborhood{
	public:
		//update neighbors each time velocity update
		vector<particle> neighbors;
		double neighborhoodBest;
		findNeighborhoodBest();
		particle Ring();
		particle vonNeuman();
}



class particle{
	public:
		double xPos;
		double yPos;
		double fitness;
		double velocity;
		//could be class of struct
		double personalBest;
		double calculateFitness(this->xPos, this->yPos);
		void updateVelocity();
		void updatePosition();


}


class swarm{

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


double evalAckley (double x, double y) {

    double firstSum = x*x + y*y;
    double secondSum = cos(2.0*M_PI*x) + cos(2.0*M_PI*y);

    return -20.0 * exp(-0.2 * sqrt(firstSum/2.0)) - 
      exp(secondSum/2.0) + 20.0 + Math.E;
}  

public double evalGriewank (vector<double> position) {


    double sumSquares = x*x + y*y;
    double productCos = cos(x/sqrt(1)) * cos(y/sqrt(2));

    return sumSquares/4000.0 - productCos + 1.0;

    double division = 1/4000.0

    double sum = 0;
    for(int i = 0; i < position.size(); i++){
    	
    }

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
public double evalRastrigin (double x, double y) {

	double retVal = 0;
    retVal += x*x - 10.0* cos(2.0*M_PI*x) + 10.0;
    retVal += y*y - 10.0* cos(2.0*M_PI*y) + 10.0;

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