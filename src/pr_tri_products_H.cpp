#include <iostream>
#include <fstream>
#include <tclap/CmdLine.h>
#include <algorithm>
#include <chrono>

#include "symplectic.h"
#include "utilities.h"
#include "stabiliser.h"
#include "GLPKConvexSeparation.h"
#include "products.h"


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

string outfile;
unsigned n;
bool elim,verbose,mode;

try {

	TCLAP::CmdLine cmd("Program for generating the vertices of a sub-polytope of the projected stabiliser polytope spanned by projections of product states of Bell type.", ' ', "0.1");

	// arguments
	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Output file name that will be used for writing the reduced constraint matrix", true, "out", "string");
	cmd.add(output_arg);

	TCLAP::ValueArg<unsigned> n_arg ("q", "qubits", "Number of qubits", true, 1, "Unsigned integer");
	cmd.add(n_arg);

	TCLAP::SwitchArg elim_arg ("e", "eliminate", "Flag which enables elimination of redundant points from the set, i.e. those points which are convex combinations of the others.", false);
	cmd.add(elim_arg);

	TCLAP::SwitchArg verb_arg ("v", "verbose", "Does what it promises.", false);
	cmd.add(verb_arg);

	TCLAP::SwitchArg mode_arg ("m", "mode", "Change mode to generation of all possible Bell-type states", true);
	cmd.add(mode_arg);

	cmd.parse(argc, argv);

	outfile = output_arg.getValue();
	n = n_arg.getValue();
	elim = elim_arg.getValue();
	mode = mode_arg.getValue();

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}


// ----------------------------------
// ------ generation of states
// ----------------------------------


// these will hold connected and product states
vector<LabelledState> states;

cout << "----------------------------------------------------------------------" << endl;
cout << "Generate projected stabiliser sub-polytope for n = " << n << endl;
cout << "----------------------------------------------------------------------" << endl;

// timing
auto t1 = chrono::high_resolution_clock::now();

projected_tri_products(n, states);

// timing
auto t2 = chrono::high_resolution_clock::now();
chrono::duration<double, milli> fp_ms = t2 - t1;

cout << "   Found " << states.size() << " tri-partite product states." << endl;
cout << "   Generation took " <<  fp_ms.count() << " ms." << endl;



// ----------------------------------
// ----- final elimination round
// ----------------------------------

sort(states.begin(), states.end());

if(elim == true) {

	cout << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "Eliminate redundant points from the total set" << endl;
	cout << "-------------------------------------------------" << endl;
	cout << endl;

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
auto t5 = chrono::high_resolution_clock::now();
chrono::duration<double, milli> fp_ms3 = t5 - t1;

cout << "Total runtime: " <<  fp_ms3.count() << " ms." << endl;


// ----------------------------------
// ----- write output
// ----------------------------------

write_states(states, outfile+to_string(n)+".mat", outfile+to_string(n)+".lab");

}