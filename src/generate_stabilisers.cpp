#include <iostream>
#include <fstream>
#include <tclap/CmdLine.h>
#include <algorithm>
#include <chrono>


#include "symplectic.h"
#include "utilities.h"
#include "stabiliser.h"


using namespace std;

/* input to the program:

	* number of qubits

*/



int main(int argc, char** argv) {

// ----------------------------------
// ----- parse comand-line parameters
// ----------------------------------

unsigned n;
string outfile;

try {

	TCLAP::CmdLine cmd("Program for generating the vertices of the projected Clifford polytope", ' ', "0.1");

	// arguments
	TCLAP::ValueArg<unsigned> nqubits_arg ("q","qubits", "Number of qubits", true, 0, "Positive integer");
	cmd.add(nqubits_arg);

	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Output file name that will be used for writing the reduced constraint matrix", true, "out.dat", "string");
	cmd.add(output_arg);

	cmd.parse(argc, argv);

	n = nqubits_arg.getValue();
	outfile = output_arg.getValue();

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}


// ----------------------------------
// ------ start computation
// ----------------------------------

cout << "Generate projected stabiliser states for n = " << n << endl;


// timing
auto t1 = chrono::high_resolution_clock::now();

// set<vector<double>> pr_states = generate_projected_stabiliser_states_set(n);
vector<vector<double>> pr_states = generate_projected_stabiliser_states(n);

auto t2 = chrono::high_resolution_clock::now();	

chrono::duration<double, milli> fp_ms = t2 - t1;


cout << "Found " << pr_states.size() << " images." << endl;
cout << "Generation took " <<  fp_ms.count() << " ms." << endl;

// save states to file
fstream fout(outfile, ios::out);
unsigned j = 1;

if(fout.is_open()) {

	for(auto state : pr_states) {
		for(unsigned k=0; k<n; k++) {
			if(state.at(k) != 0) {
				fout << j << " " << k+1 << " " << scientific << state.at(k) << endl;
			}
		}
		++j;
	}

}

fout.close();


}