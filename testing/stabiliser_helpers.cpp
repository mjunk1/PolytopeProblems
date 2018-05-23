#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <set>
#include <chrono>
#include <tclap/CmdLine.h>


#include "symplectic.h"
#include "utilities.h"
#include "stabiliser.h"


using namespace std;

// ----- Test of stabiliser generation using std::set in symplectic.h

int main(int argc, char** argv) {

unsigned n;
string infile;


cout << "Test product state projection" << endl;
cout << "-----------------------------" << endl;


vector<int> rho = {0,2};
vector<int> sigma = {-1};
vector<int> prod = pr_product_state(rho,sigma);

cout << "rho = ( ";
for(auto c : rho) {
	cout << c << ", ";
}
cout << " )" << endl;

cout << "sigma = ( ";
for(auto c : sigma) {
	cout << c << ", ";
}
cout << " )" << endl;


cout << "rho (x) sigma = ( ";
for(auto c : prod) {
	cout << c << ", ";
}
cout << " )" << endl;


}