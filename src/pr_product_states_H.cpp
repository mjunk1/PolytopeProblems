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


/* Program to generate a sub-polytope of the projected stabiliser polytope which yields an approximate optimisation problem for RoM.

	The sub-polytope is spanned by product states made out of Bell-type factors and X eigenstates.

	Input to the program:
		(*) number of qubits 

*/



int main(int argc, char** argv) {

// ----------------------------------
// ----- parse comand-line parameters
// ----------------------------------

string outfile,infile;
unsigned n, max_factor;
bool elim,verbose,read_con=false;

try {

	TCLAP::CmdLine cmd("Program for generating the vertices of a sub-polytope of the projected stabiliser polytope spanned by projections of product states.", ' ', "0.1");

	// arguments
	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file name that contains the adjacency matrices for the non-isomorphic non-LC-equivalent graphs with n vertices in graph6 format", true, "in.dat", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Output file name that will be used for writing the reduced constraint matrix", true, "out", "string");
	cmd.add(output_arg);

	TCLAP::ValueArg<unsigned> n_arg ("q", "qubits", "Number of qubits", true, 1, "Unsigned integer");
	cmd.add(n_arg);

	TCLAP::ValueArg<unsigned> max_arg ("m", "max-factor", "Maximum number of qubit in a tensor factor", false, 2, "Unsigned integer");
	cmd.add(max_arg);

	TCLAP::SwitchArg elim_arg ("e", "eliminate", "Flag which enables elimination of redundant points from the set, i.e. those points which are convex combinations of the others.", false);
	cmd.add(elim_arg);

	TCLAP::SwitchArg verb_arg ("v", "verbose", "Does what it promises.", false);
	cmd.add(verb_arg);

	cmd.parse(argc, argv);

	outfile = output_arg.getValue();
	infile = input_arg.getValue();
	n = n_arg.getValue();
	max_factor = max_arg.getValue();
	elim = elim_arg.getValue();

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}


// ----------------------------------
// ------ generation of states
// ----------------------------------


// these will the states
vector<LabelledState> states;
vector<LabelledState> con_states;

cout << "----------------------------------------------------------------------" << endl;
cout << "Generate projected stabiliser sub-polytope for n = " << n << endl;
cout << "----------------------------------------------------------------------" << endl;


// write program call to stdout
cout << endl;
cout << "# Program call" << endl;
cout << "#-------------" << endl;
cout << "   ";
for(unsigned i=0; i<argc; i++) {
	cout << argv[i] << " ";
}
cout << endl << endl;


if(n <= max_factor) {
	// read connected vertices
	cout << endl;
	cout << "# Try to read connected vertices from file" << endl;
	cout << "#-----------------------------------------" << endl;
	cout << endl;

	// test for already exisiting file
	if(get_states(con_states, infile+to_string(n)+"_con.mat", infile+to_string(n)+"_con.lab") == 0) {
		read_con = true;

		cout << "   Read " << con_states.size() << " connected vertices." << endl;
	}
	else {
		cout << "   Failed." << endl;
	}

}


// timing
auto t1 = chrono::high_resolution_clock::now();

if( pr_product_states(n, max_factor, states, infile+"%u_con.mat", infile+"%u_con.lab") != 0 ) {
	cout << "Couldn't read vertices from lower-dimensional problems from files " << outfile+"%u_con.mat" << endl;
	cout << "Will now hold." << endl;
	return 1;
}

// timing
auto t2 = chrono::high_resolution_clock::now();
chrono::duration<double, milli> fp_ms = t2 - t1;

cout << "   Found " << states.size() << " images of product states." << endl;
cout << "   Generation took " <<  fp_ms.count() << " ms." << endl;


if(read_con) {
	// merge
	copy(con_states.begin(), con_states.end(), back_inserter(states));

	// free memory
	con_states.clear();
	con_states.resize(0);
	con_states.shrink_to_fit();
}


// ----------------------------------
// ----- final elimination round
// ----------------------------------


if(elim == true) {

	cout << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "Eliminate redundant points from the total set" << endl;
	cout << "-------------------------------------------------" << endl;
	cout << endl;

	sort(states.begin(), states.end());

	GLPKConvexSeparation lp (n);
	if(verbose == false) {
		lp.set_verbosity(1);
	}
	else {
		lp.set_verbosity(2);
	}

	lp.add_vertices(states);
	cout << "   Found " << lp.get_nvertices() << " extremal points." << endl;
	lp.print_parameters();

	lp.delete_redundant_points();
	cout << "   Found " << lp.get_nvertices() << " extremal points." << endl;
	lp.print_parameters();

 	states = lp.iget_labelled_vertices();

 	lp.write_constraint_matrix(outfile+to_string(n)+".coo");
}


// timing
auto t3 = chrono::high_resolution_clock::now();
chrono::duration<double, milli> fp_ms2 = t3 - t1;

cout << "Total runtime: " <<  fp_ms2.count() << " ms." << endl;


// ----------------------------------
// ----- write output
// ----------------------------------

write_states(states, outfile+to_string(n)+".mat", outfile+to_string(n)+".lab");

}