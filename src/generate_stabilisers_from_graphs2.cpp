#include <iostream>
#include <fstream>
#include <tclap/CmdLine.h>
#include <algorithm>
#include <chrono>

#include "symplectic.h"
#include "utilities.h"
#include "stabiliser.h"
#include "GLPKConvexSeparation.h"


using namespace std;


/* input to the program:

	* number of qubits

*/



int main(int argc, char** argv) {

// ----------------------------------
// ----- parse comand-line parameters
// ----------------------------------

unsigned n;
string outfile,infile;
bool verbose;

try {

	TCLAP::CmdLine cmd("Program for generating the vertices of the projected Clifford polytope", ' ', "0.1");

	// arguments
	TCLAP::ValueArg<unsigned> nqubits_arg ("q","qubits", "Number of qubits", true, 0, "Positive integer");
	cmd.add(nqubits_arg);

	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file name that contains the adjacency matrices for the non-isomorphic graphs with q vertices", true, "in.dat", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Output file name that will be used for writing the reduced constraint matrix", false, "", "string");
	cmd.add(output_arg);

	TCLAP::SwitchArg verb_arg ("v", "verbose", "Does what it promises.", false);
	cmd.add(verb_arg);


	cmd.parse(argc, argv);

	n = nqubits_arg.getValue();
	outfile = output_arg.getValue();
	infile = input_arg.getValue();
	verbose = verb_arg.getValue();

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}


// ----------------------------------
// ------ start generation
// ----------------------------------

cout << "------------------------------------------------------------" << endl;
cout << "Generate projected stabiliser states from graphs for n = " << n << endl;
cout << "------------------------------------------------------------" << endl;

// timing
auto t1 = chrono::high_resolution_clock::now();

GLPKConvexSeparation lp = generate_projected_stabiliser_states_from_graphs_conv(infile, n);
unsigned nvertices = lp.get_nvertices();

auto t2 = chrono::high_resolution_clock::now();	

chrono::duration<double, milli> fp_ms = t2 - t1;

cout << "Found " << nvertices << " images." << endl;
cout << "Generation took " <<  fp_ms.count() << " ms." << endl;

auto t3 = chrono::high_resolution_clock::now();	

lp.print_parameters();
lp.delete_redundant_points();
cout << "Deleted " << nvertices - lp.get_nvertices() << " points." << endl;
lp.print_parameters();

auto t4 = chrono::high_resolution_clock::now();	

fp_ms = t4 - t3;
chrono::duration<double, milli> fp_ms2 = t4 - t1;

cout << "Elimination took " <<  fp_ms.count() << " ms." << endl;
cout << "Total time " <<  fp_ms2.count() << " ms." << endl;

if(outfile != "") {

	// save states to file in sparse COO format
	lp.write_constraint_matrix(outfile+"_red.coo");

	// save states to file in dense format
	lp.write_dense_constraint_matrix(outfile+"_red.mat");

	// write the labels
	lp.write_labels(outfile+"_red_id.dat");

}

}