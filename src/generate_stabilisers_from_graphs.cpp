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
	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file name that contains the adjacency matrices for the non-isomorphic graphs with q vertices", true, "in.dat", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Output file name that will be used for writing the reduced constraint matrix", false, "", "string");
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
// ------ start generation
// ----------------------------------

unsigned n = get_order_from_file(infile);

cout << "-------------------------------------------------------------" << endl;
cout << "Generate projected stabiliser states from graphs for n = " << n << endl;
cout << "-------------------------------------------------------------" << endl;

// timing
auto t1 = chrono::high_resolution_clock::now();

// vector<vector<double>> pr_states = pr_stabiliser_from_graphs(infile, n);
pair< vector<vector<int>>, vector<string> > pr_states = pr_stabiliser_from_graphs(infile);

auto t2 = chrono::high_resolution_clock::now();	

chrono::duration<double, milli> fp_ms = t2 - t1;


cout << "Found " << pr_states.first.size() << " images." << endl;
cout << "Generation took " <<  fp_ms.count() << " ms." << endl;

if(outfile != "") {

	// save states to file in sparse COO format
	fstream fout(outfile+".coo", ios::out);
	unsigned j = 1;

	if(fout.is_open()) {

		for(auto state : pr_states.first) {
			for(unsigned k=0; k<state.size(); k++) {
				if(state.at(k) != 0) {
					fout << j << " " << k+1 << " " << state.at(k) << endl;
				}
			}
			++j;
		}

	}
	fout.close();

	// save states to file in dense format
	fout.open(outfile+".mat", ios::out);

	if(fout.is_open()) {
		for(auto state : pr_states.first) {
			for(unsigned k=0; k<state.size(); k++) {
				fout  << state.at(k) << " ";
			}
			fout << endl;
		}
	}
	fout.close();

	// write the labels
	fout.open(outfile+"_id.dat", ios::out);

	if(fout.is_open()) {
		for(auto label : pr_states.second) {
			fout << label << endl;
		}
	}
	fout.close();

}

// ----------------------------------
// ----- start elimination
// ----------------------------------

if(elim == true && outfile != "") {
	// time it
	auto t3 = chrono::high_resolution_clock::now();	

	// free memory
	pr_states.first.clear();
	pr_states.first.resize(0);
	pr_states.first.shrink_to_fit();

	cout << endl;
	cout << "-------------------------------------------------" << endl;
	cout << "Eliminate redundant points from the generated set" << endl;
	cout << "-------------------------------------------------" << endl;
	cout << endl;

	GLPKConvexSeparation lp (outfile+".coo");
	if(verbose == false) {
		lp.set_verbosity(1);
	}
	lp.set_labels(pr_states.second);
	int ret_status;

	lp.print_parameters();
	unsigned nvertices = lp.get_nvertices();
	lp.delete_redundant_points();
	cout << "Deleted " << nvertices - lp.get_nvertices() << " points." << endl;
	lp.print_parameters();

	auto t4 = chrono::high_resolution_clock::now();	
	fp_ms = t4 - t3;
	chrono::duration<double, milli> fp_ms2 = t4 - t1;

	cout << "Elimination took " <<  fp_ms.count() << " ms." << endl;
	cout << "Total time " <<  fp_ms2.count() << " ms." << endl;


	lp.write_constraint_matrix(outfile+"_red.coo");
	lp.write_dense_constraint_matrix(outfile+"_red.mat");

	// write the labels
	fstream fout(outfile+"_red_id.dat", ios::out);

	if(fout.is_open()) {
		vector<string> labels = lp.get_labels();

		for(auto label : labels) {
			fout << label << endl;
		}

	}
	fout.close();
}

}