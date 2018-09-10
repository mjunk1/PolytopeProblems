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
bool elim,verbose;

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

	cmd.parse(argc, argv);

	outfile = output_arg.getValue();
	n = n_arg.getValue();
	elim = elim_arg.getValue();

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}


// ----------------------------------
// ------ generation of states
// ----------------------------------


// these will hold connected and product states
vector<vector<int>> states;
vector<string> labels;

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


// timing
auto t1 = chrono::high_resolution_clock::now();


// get all partitions of n/2 of length 3: n/2 = k_1 + k_2 + k_3
// then use k_1 psi+, k_2 psi- and 2*k_3 (+1) +/- states

for(unsigned k1=0; k1<=n/2; k1++) {
	for(unsigned k2=0; k2<=n/2-k1; k2++) {
		unsigned k3 = n/2-k1-k2;

		vector<int> state;
		string label;

		for(unsigned i=0; i<k1; i++) {
			state = pr_product_state(pr_psiplus, state);
			label = pr_product_label(pr_psiplus_label, label);
		}
		for(unsigned i=0; i<k2; i++) {
			state = pr_product_state(pr_psiminus, state);
			label = pr_product_label(pr_psiminus_label, label);
		}

		vector<int> statem = state;
		string labelm = label;

		if(2*k3+n%2 != 0) {
			state = pr_product_state(pr_allplus(2*k3+n%2), state);
			label = pr_product_label(pr_allplus_label(2*k3+n%2), label);

			statem = pr_product_state(pr_allminus(2*k3+n%2), statem);
			labelm = pr_product_label(pr_allminus_label(2*k3+n%2), labelm);

			states.push_back(statem);
			labels.push_back(labelm);
		}

		// append
		states.push_back(state);
		labels.push_back(label);
		
	}
}


// timing
auto t2 = chrono::high_resolution_clock::now();
chrono::duration<double, milli> fp_ms = t2 - t1;

cout << "   Found " << states.size() << " Bell-type product states." << endl;
cout << "   Generation took " <<  fp_ms.count() << " ms." << endl;



// ----------------------------------
// ----- final elimination round
// ----------------------------------

GLPKConvexSeparation lp (states);
if(verbose == false) {
	lp.set_verbosity(1);
}
else {
	lp.set_verbosity(2);
}
lp.set_labels(labels);


if(elim == true) {

	cout << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "Eliminate redundant points from the total set" << endl;
	cout << "-------------------------------------------------" << endl;
	cout << endl;

	lp.print_parameters();
	int dels = lp.delete_redundant_points();
	cout << "   Deleted " << dels << " points." << endl;
	lp.print_parameters();

 	states = lp.iget_vertices();
	labels = lp.get_labels();
}

lp.write_constraint_matrix(outfile+to_string(n)+".coo");

// timing
auto t5 = chrono::high_resolution_clock::now();
chrono::duration<double, milli> fp_ms3 = t5 - t1;

cout << "Total runtime: " <<  fp_ms3.count() << " ms." << endl;


// ----------------------------------
// ----- write output
// ----------------------------------

write_states(states, outfile+to_string(n)+".mat");
write_labels(labels, outfile+to_string(n)+".lab");

}