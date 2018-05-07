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

try {

	TCLAP::CmdLine cmd("Program for generating the vertices of the projected Clifford polytope", ' ', "0.1");

	// arguments
	TCLAP::ValueArg<unsigned> nqubits_arg ("q","qubits", "Number of qubits", true, 0, "Positive integer");
	cmd.add(nqubits_arg);

	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file name that contains the adjacency matrices for the non-isomorphic graphs with q vertices", true, "in.dat", "string");
	cmd.add(input_arg);

	cmd.parse(argc, argv);

	n = nqubits_arg.getValue();
	infile = input_arg.getValue();

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}


cout << "Test graph stuff" << endl;
cout << "-----------------" << endl;


vector<unsigned> graphs = read_graph_states(infile,n);

cout << "Read " << graphs.size() << " graphs !" << endl;


}