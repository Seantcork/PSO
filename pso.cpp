// Authors Sean Cork, Kamaal Palmer, Luca Ostertag-Hill


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


double evalAckley (vector<double> dimensionVector) {

    double firstSum = x*x + y*y;
    double secondSum = Math.cos(2.0*Math.PI*x) + Math.cos(2.0*Math.PI*y);

    return -20.0 * Math.exp(-0.2 * Math.sqrt(firstSum/2.0)) - 
      Math.exp(secondSum/2.0) + 20.0 + Math.E;
}  

public double evalGriewank (vector<double> dimensionVector) {

    double sumSquares = x*x + y*y;
    double productCos = Math.cos(x/Math.sqrt(1)) * Math.cos(y/Math.sqrt(2));

    return sumSquares/4000.0 - productCos + 1.0;
 }  


public double evalRosenbrock (vector<double> dimensionVector) {

    return 100.0 * Math.pow(y - x*x, 2.0) + Math.pow(x-1.0, 2.0);
  }




  // returns the value of the Rastrigin Function at point (x, y)
  //   minimum is 0.0, which occurs at (0.0,...,0.0)
public double evalRastrigin (vector<double> dimensionVector) {

	double retVal = 0;
    retVal += x*x - 10.0*Math.cos(2.0*Math.PI*x) + 10.0;
    retVal += y*y - 10.0*Math.cos(2.0*Math.PI*y) + 10.0;

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