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

string outfile,infile;
bool elim,verbose;

try {

	TCLAP::CmdLine cmd("Program for generating the vertices of the projected Clifford polytope", ' ', "0.1");

	// arguments
	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file name that contains the adjacency matrices for the non-isomorphic non-LC-equivalent graphs with q vertices in graph6 format", true, "in.dat", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Output file name that will be used for writing the reduced constraint matrix", true, "out", "string");
	cmd.add(output_arg);

	TCLAP::SwitchArg elim_arg ("e", "eliminate", "Flag which determines wether to eliminate redundant points from the set, i.e. those points which are convex combinations of the others.", false);
	cmd.add(elim_arg);

	TCLAP::SwitchArg verb_arg ("v", "verbose", "Does what it promises.", false);
	cmd.add(verb_arg);


	cmd.parse(argc, argv);

	outfile = output_arg.getValue();
	infile = input_arg.getValue();
	elim = elim_arg.getValue();
	verbose = verb_arg.getValue();

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}


// ----------------------------------
// ------ generation of states
// ----------------------------------

unsigned n = get_order_from_file(infile);

cout << "----------------------------------------------------------------------" << endl;
cout << "Generate projected stabiliser states from LC-orbits of graphs for n = " << n << endl;
cout << "----------------------------------------------------------------------" << endl;


// ---- generate orbits of connected graphs

cout << endl;
cout << "# Generate orbits of connected graphs" << endl;
cout << "# -----------------------------------" << endl;
cout << endl;

// timing
auto t1 = chrono::high_resolution_clock::now();

auto pr_con = pr_stabiliser_from_LCgraphs2(infile);

// timing
auto t2 = chrono::high_resolution_clock::now();
chrono::duration<double, milli> fp_ms = t2 - t1;

cout << "   Found " << pr_con.first.size() << " images of connected graphs." << endl;
cout << "   Generation took " <<  fp_ms.count() << " ms." << endl;


// ---- elimination of non-extremal connected states

if(elim == true) {

	cout << endl;
	cout << "# Eliminate redundant points from the connected set" << endl;
	cout << "# -------------------------------------------------" << endl;
	cout << endl;

	GLPKConvexSeparation lp (pr_con.first);
	if(verbose == false) {
		lp.set_verbosity(1);
	}
	lp.set_labels(pr_con.second);
	int ret_status;

	lp.print_parameters();
	unsigned nvertices = lp.get_nvertices();
	lp.delete_redundant_points();
	cout << "   Deleted " << nvertices - lp.get_nvertices() << " points." << endl;
	lp.print_parameters();

	pr_con = lp.iget_vertices();

}


// ---- generate product states of lower-dimensional extremal points, yielding disconnected states

if(n > 1) {

	cout << endl;
	cout << "# Generate product states of lower-dimensional extremal points" << endl;
	cout << "# ------------------------------------------------------------" << endl;
	cout << endl;

	// timing
	auto t3 = chrono::high_resolution_clock::now();

	auto pr_prod = pr_product_states(n, outfile+"%u.mat", outfile+"%u.lab");

	// timing
	auto t4 = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> fp_ms2 = t4 - t3;
	
	cout << "   Found " << pr_prod.first.size() << " images of product states." << endl;
	cout << "   Generation took " <<  fp_ms2.count() << " ms." << endl;

	// elimination of non-extremal disconnected states
	if(elim == true) {
		cout << endl;
		cout << "# Eliminate redundant points from the disconnected set" << endl;
		cout << "# ----------------------------------------------------" << endl;
		cout << endl;

		GLPKConvexSeparation lp (pr_prod.first);
		if(verbose == false) {
			lp.set_verbosity(1);
		}
		lp.set_labels(pr_prod.second);

		lp.print_parameters();
		unsigned nvertices = lp.get_nvertices();
		lp.delete_redundant_points();
		cout << "   Deleted " << nvertices - lp.get_nvertices() << " points." << endl;
		lp.print_parameters();

		pr_prod = lp.iget_vertices();
	}

	// merge disconnected and connected states
	copy(pr_prod.first.begin(), pr_prod.first.end(), back_inserter(pr_con.first));
	copy(pr_prod.second.begin(), pr_prod.second.end(), back_inserter(pr_con.second));

	// free memory
	pr_prod.first.clear();
	pr_prod.first.resize(0);
	pr_prod.first.shrink_to_fit();
}



// ----------------------------------
// ----- final elimination round
// ----------------------------------

if(elim == true && n > 1) {

	cout << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "Eliminate redundant points from the total set" << endl;
	cout << "-------------------------------------------------" << endl;
	cout << endl;

	GLPKConvexSeparation lp (pr_con.first);
	if(verbose == false) {
		lp.set_verbosity(1);
	}
	lp.set_labels(pr_con.second);

	lp.print_parameters();
	unsigned nvertices = lp.get_nvertices();
	lp.delete_redundant_points();
	cout << "   Deleted " << nvertices - lp.get_nvertices() << " points." << endl;
	lp.print_parameters();

	lp.write_constraint_matrix(outfile+to_string(n)+".coo");

	pr_con = lp.iget_vertices();
}

// timing
auto t5 = chrono::high_resolution_clock::now();
chrono::duration<double, milli> fp_ms3 = t5 - t1;

cout << "Total runtime: " <<  fp_ms3.count() << " ms." << endl;


// ----------------------------------
// ----- write output
// ----------------------------------

write_states(pr_con, outfile+to_string(n)+".mat", outfile+to_string(n)+".lab");

}