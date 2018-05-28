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
bool elim,verbose,con,dis,save_proj;

try {

	TCLAP::CmdLine cmd("Program for generating the vertices of the projected stabiliser polytope", ' ', "0.1");

	// arguments
	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file name that contains the adjacency matrices for the non-isomorphic non-LC-equivalent graphs with n vertices in graph6 format", true, "in.dat", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Output file name that will be used for writing the reduced constraint matrix", true, "out", "string");
	cmd.add(output_arg);

	TCLAP::SwitchArg elim_arg ("e", "eliminate", "Flag which disables eliminattion of redundant points from the set, i.e. those points which are convex combinations of the others.", true);
	cmd.add(elim_arg);

	TCLAP::SwitchArg verb_arg ("v", "verbose", "Does what it promises.", false);
	cmd.add(verb_arg);

	TCLAP::SwitchArg con_arg ("c", "connected", "Flag which disables the generation of connected states. If the flag is set, the program will check for the existence of a file containing these states and read them if is found.", true);
	cmd.add(con_arg);

	TCLAP::SwitchArg dis_arg ("d", "disconnected", "Flag which disables the generation of disconnected, i.e. product, states. If the flag is set, the program will check for the existence of a file containing these states and read them if is found.", true);
	cmd.add(dis_arg);

	TCLAP::SwitchArg proj_arg ("p", "save_projections", "Save all projections of stabiliser orbits to hard drive (only those from connected graphs).", false);
	cmd.add(proj_arg);


	cmd.parse(argc, argv);

	outfile = output_arg.getValue();
	infile = input_arg.getValue();
	elim = elim_arg.getValue();
	verbose = verb_arg.getValue();
	con = con_arg.getValue();
	dis = dis_arg.getValue();
	save_proj = proj_arg.getValue();


} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}


// ----------------------------------
// ------ generation of states
// ----------------------------------

unsigned n = get_order_from_file(infile);

// these will hold connected and product states
vector<LabelledState> con_states;
vector<LabelledState> prod_states;


cout << "----------------------------------------------------------------------" << endl;
cout << "Generate projected stabiliser states from LC-orbits of graphs for n = " << n << endl;
cout << "----------------------------------------------------------------------" << endl;

// timing
auto t1 = chrono::high_resolution_clock::now();

// ---- generate orbits of connected graphs
// ----------------------------------------

if(con == true) {
	cout << endl;
	cout << "# Generate orbits of connected graphs" << endl;
	cout << "# -----------------------------------" << endl;
	cout << endl;

	if(pr_stabiliser_from_LCgraphs(infile, con_states) != 0) {
		cout << "Couldn't read input file " << infile << endl;
		cout << "Will now hold." << endl;
		return 1;
	}

	// timing
	auto t2 = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> fp_ms = t2 - t1;

	cout << "   Found " << con_states.size() << " images of connected graphs." << endl;
	cout << "   Generation took " <<  fp_ms.count() << " ms." << endl;

	// ---- save orbits
	if(save_proj == true) {
		write_states(con_states, outfile+to_string(n)+"_con_orbits.mat", outfile+to_string(n)+"_con_orbits.lab");
	}

	// ---- elimination of non-extremal connected states

	if(elim == true) {

		cout << endl;
		cout << "# Eliminate redundant points from the connected set" << endl;
		cout << "# -------------------------------------------------" << endl;
		cout << endl;

		GLPKConvexSeparation lp (con_states);
		if(verbose == false) {
			lp.set_verbosity(1);
		}

		lp.print_parameters();
		unsigned nvertices = lp.get_nvertices();
		lp.delete_redundant_points();
		cout << "   Deleted " << nvertices - lp.get_nvertices() << " points." << endl;
		lp.print_parameters();

		con_states = lp.iget_labelled_vertices();

	}

	// write connected states
	write_states(con_states, outfile+to_string(n)+"_con.mat", outfile+to_string(n)+"_con.lab");
}
else { 
	cout << endl;
	cout << "# Try to read connected vertices from file" << endl;
	cout << "# -----------------------------------" << endl;
	cout << endl;

	// test for already exisiting file
	if(get_states(con_states, outfile+to_string(n)+"_con.mat", outfile+to_string(n)+"_con.lab") == 0) {

		cout << "   Read " << con_states.size() << " connected vertices." << endl;
	}
	else {
		cout << "   Failed." << endl;
	}


}


// ---- generate product states of lower-dimensional extremal points, yielding disconnected states
// -----------------------------------------------------------------------------------------------

if(n > 1 && dis == true) {

	cout << endl;
	cout << "# Generate product states of lower-dimensional extremal points" << endl;
	cout << "# ------------------------------------------------------------" << endl;
	cout << endl;

	// timing
	auto t3 = chrono::high_resolution_clock::now();

	if( pr_product_states(n, prod_states, outfile+"%u.mat", outfile+"%u.lab") != 0 ) {
		cout << "Couldn't read vertices from lower-dimensional problems from files " << outfile+"%u.mat" << endl;
		cout << "Will now hold." << endl;
		return 1;
	}

	// timing
	auto t4 = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> fp_ms2 = t4 - t3;
	
	cout << "   Found " << prod_states.size() << " images of product states." << endl;
	cout << "   Generation took " <<  fp_ms2.count() << " ms." << endl;

	// elimination of non-extremal disconnected states
	if(elim == true) {
		cout << endl;
		cout << "# Eliminate redundant points from the disconnected set" << endl;
		cout << "# ----------------------------------------------------" << endl;
		cout << endl;

		GLPKConvexSeparation lp (prod_states);
		if(verbose == false) {
			lp.set_verbosity(1);
		}

		lp.print_parameters();
		unsigned nvertices = lp.get_nvertices();
		lp.delete_redundant_points();
		cout << "   Deleted " << nvertices - lp.get_nvertices() << " points." << endl;
		lp.print_parameters();

		prod_states = lp.iget_labelled_vertices();
	}

	// merge disconnected and connected states
	copy(prod_states.begin(), prod_states.end(), back_inserter(con_states));

	// free memory
	prod_states.clear();
	prod_states.resize(0);
	prod_states.shrink_to_fit();



	// ----------------------------------
	// ----- final elimination round
	// ----------------------------------

	GLPKConvexSeparation lp (con_states);
	if(verbose == false) {
		lp.set_verbosity(1);
	}

	if(elim == true) {

		cout << endl;
		cout << "-------------------------------------------------" << endl;
		cout << "Eliminate redundant points from the total set" << endl;
		cout << "-------------------------------------------------" << endl;
		cout << endl;

		lp.print_parameters();
		unsigned nvertices = lp.get_nvertices();
		lp.delete_redundant_points();
		cout << "   Deleted " << nvertices - lp.get_nvertices() << " points." << endl;
		lp.print_parameters();

		con_states = lp.iget_labelled_vertices();
	}

	lp.write_constraint_matrix(outfile+to_string(n)+".coo");

	// timing
	auto t5 = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> fp_ms3 = t5 - t1;

	cout << "Total runtime: " <<  fp_ms3.count() << " ms." << endl;
}

if(dis == true) {

	// ----------------------------------
	// ----- write output
	// ----------------------------------

	write_states(con_states, outfile+to_string(n)+".mat", outfile+to_string(n)+".lab");

	// this is just for comparison with older data ...
	// sort(pr_con.first.begin(), pr_con.first.end());
	// write_states(pr_con, outfile+to_string(n)+"_sorted.mat");
}

}