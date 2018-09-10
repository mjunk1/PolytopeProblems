#include <iostream>
#include <fstream>
#include <streambuf>
#include <tclap/CmdLine.h>
#include <algorithm>
#include <chrono>

#include "symplectic.h"
#include "utilities.h"
#include "stabiliser.h"
#include "GLPKConvexSeparation.h"


using namespace std;



int main(int argc, char** argv) {

// ----------------------------------
// ----- parse comand-line parameters
// ----------------------------------

string infile;

try {

	TCLAP::CmdLine cmd("Program for generating the vertices of the projected stabiliser polytope", ' ', "0.1");

	// arguments
	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file name that contains the adjacency matrices for the non-isomorphic non-LC-equivalent graphs with n vertices in graph6 format", true, "in.dat", "string");
	cmd.add(input_arg);

	cmd.parse(argc, argv);

	infile = input_arg.getValue();

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}


// ----------------------------------
// ------ generation of states
// ----------------------------------


cout << "#-------------------------------------------------------" << endl;
cout << "# Select projected stabiliser states by purity criterion" << endl;
cout << "#-------------------------------------------------------" << endl;


// write program call to stdout
cout << endl;
cout << "# Program call" << endl;
cout << "#-------------" << endl;
cout << "   ";
for(unsigned i=0; i<argc; i++) {
	cout << argv[i] << " ";
}
cout << endl << endl;



// read states
vector<LabelledState> states;


cout << endl;
cout << "# Read states from file" << endl;
cout << "#----------------------" << endl;
cout << endl;

// test for already exisiting file
if(get_states(states, infile+".mat", infile+".lab") == 0) {
	cout << "   Read " << states.size() << " states." << endl;
}
else {
	cout << "   Failed." << endl;
}


// sort by decreasing purity
sort(states.begin(), states.end(), compare_by_purity2);

// compute mean purity
double mean_purity = 0;
unsigned nstates = states.size();
vector<double> purities (nstates, 0.);

for(unsigned i=0; i<nstates; i++) {
	purities.at(i) = purity(states.at(i));
	mean_purity += purities.at(i);
}
mean_purity /= nstates;

// delete all states with purity less than the mean
auto it = find_if(purities.begin(), purities.end(), [&](double val){ return val < mean_purity; });
unsigned d = distance(purities.begin(), it);
states.erase(states.begin()+d, states.end());

cout << "   Eliminated " << nstates - states.size() << " points using purity criterion." << endl << endl;

// ---- save states
write_states(states, infile+"_pure.mat", infile+"_pure.lab");


}