#include <math.h>
#include <fstream>
#include <cstdlib>
#include <stdint.h>
#include <random>
#include <string>
#include <vector>
#include <numeric>
#include <iostream>
#include <limits>
#include <set>

using namespace std;

const int ACK_VEC = 1;
const int ROK_VEC = 2;
const int RAS_VEC = 3;
const string GLOBAL_TOPOLOGY = "gl";
const string RING_TOPOLOGY = "ri";
const string VON_NEUMANN_TOPOLOGY = "vn";
const string RANDOM_TOPOLOGY = "ra";

const string ROSENBROCK_FUNCTION = "rok";
const string ACKLEY_FUNCTION = "ack";
const string RASTRIGIN_FUNCTION = "ras";


int main() {
	ifstream myfile;
	myfile.open("experimentResults.txt");
	string line;
	vector<vector<vector<vector<vector<float> > > > > resultsVec(3, vector<vector<vector<vector<float> > > >(4, vector<vector<vector<float> > >(3, vector<vector<float> >(11, vector<float>()))));
	int a, b, c, d, numLinesSince = 0;

	while(getline(myfile, line)) {
		if(numLinesSince == 0 && (line.find("Topology") != string::npos)) {

			if(line.substr(29,2) == "16") {
				c = 0;
			}
			else if(line.substr(29,2) == "30") {
				c = 1;
			}
			else if(line.substr(29,2) == "49") {
				c = 2;
			}
			if(line.substr(74,3) == ACKLEY_FUNCTION) {
				a = 0;
			}
			else if(line.substr(74,3) == ROSENBROCK_FUNCTION) {
				a = 1;
			}
			else if(line.substr(74,3) == RASTRIGIN_FUNCTION) {
				a = 2;
			}

			if(line.substr(15,2) == GLOBAL_TOPOLOGY) {
				b = 0;
			}
			else if(line.substr(15,2) == RING_TOPOLOGY) {
				b = 1;
			}
			else if(line.substr(15,2) == VON_NEUMANN_TOPOLOGY) {
				b = 2;
			}
			else if(line.substr(15,2) == RANDOM_TOPOLOGY) {
				b = 3;
			}
			numLinesSince++;
			continue;
		}

		if(line.substr(26, 1).compare("0") == 0){
			d = 0;
			resultsVec[a][b][c][d].push_back(stof(line.substr(31)));
		}
		if(line.substr(26,4).compare("1000") == 0) {
			d = 1;
			resultsVec[a][b][c][d].push_back(stof(line.substr(34)));
		}
		else if(line.substr(26,4).compare("2000") == 0) {
			d = 2;
			resultsVec[a][b][c][d].push_back(stof(line.substr(34)));
		}
		else if(line.substr(26,4).compare("3000") == 0) {
			d = 3;
			resultsVec[a][b][c][d].push_back(stof(line.substr(34)));
		}
		else if(line.substr(26,4).compare("4000") == 0) {
			d = 4;
			resultsVec[a][b][c][d].push_back(stof(line.substr(34)));
		}
		else if(line.substr(26,4).compare("5000") == 0) {
			d = 5;
			resultsVec[a][b][c][d].push_back(stof(line.substr(34)));
		}
		else if(line.substr(26,4).compare("6000") == 0) {
			d = 6;
			resultsVec[a][b][c][d].push_back(stof(line.substr(34)));
		}
		else if(line.substr(26,4).compare("7000") == 0) {
			d = 7;
			resultsVec[a][b][c][d].push_back(stof(line.substr(34)));
		}
		else if(line.substr(26,4).compare("8000") == 0) {
			d = 8;
			resultsVec[a][b][c][d].push_back(stof(line.substr(34)));
		}
		else if(line.substr(26,4).compare("9000") == 0) {
			d = 9;
			resultsVec[a][b][c][d].push_back(stof(line.substr(34)));
		}
		if(line.find("Best Overall Fitness") != string::npos ){
			resultsVec[a][b][c][10].push_back(stof(line.substr(28)));
			numLinesSince = 0;
		}
	}	// printing part

	for(int a = 0; a < 3; a++) {
		if(a == 0) {
			cout << "Ackley Function " << endl; 
		}
		else if (a == 1) {
			cout << "Rosenbrock Function"<< endl; 
		}
		else if (a == 2) {
			cout << "Rastrigin Function"<< endl; 
		}

		for(int b = 0; b < 4; b++ ) {
			if(b == 0 ) {
				cout <<" Global Topology" << endl; ;
			}
			else if (b == 1) {
				cout << " Ring Topology" << endl; ;
			}
			else if (b == 2) {
				cout << " Von Neumann Topology" << endl; ;
			}
			else if (b == 3 ) {
				cout << " Random Topology" << endl; ;
			}
			for(int c = 0; c < 3; c++ ) {
				if(c == 0) {
					cout << " Swarm Size: 16 " << endl; ;
				}
				else if(c == 1) {
					cout << " Swarm Size: 30 " << endl; ;
				}
				else if(c == 2) {
					cout << " Swarm Size: 49 " << endl; ;
				}
				for(int d = 0; d < 11; d++) {
					sort(resultsVec[a][b][c][d].begin(), resultsVec[a][b][c][d].end());
					if(d == 0){
						cout << "Iteration 0: ";
					}
					if(d == 1) {
						cout << "Iteration 1000: ";
					}
					else if(d == 2) {
						cout << "Iteration 2000: ";
					}
					else if(d == 3) {
						cout << "Iteration 3000: ";
					}
					else if(d == 4) {
						cout << "Iteration 4000: ";
					}
					else if(d == 5) {
						cout << "Iteration 5000: ";
					}
					else if(d == 6) {
						cout << "Iteration 6000: ";
					}
					else if(d == 7) {
						cout << "Iteration 7000: ";
					}
					else if(d == 8) {
						cout << "Iteration 8000: ";
					}
					else if(d == 9) {
						cout << "Iteration 9000: ";
					}
					else if(d == 10) {
						cout << "Iteration 10,000: ";
					}
					// cout << "resultsVec val: " << resultsVec[a][b][c][d][5] << endl;
					cout << "Median: " << resultsVec[a][b][c][d][9] << " Mean: " << (accumulate(resultsVec[a][b][c][d].begin(), resultsVec[a][b][c][d].end(), 0.0) / 20.0) << endl;
				}
			}

		}
	}


}
