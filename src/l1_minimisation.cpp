#include <iostream>
#include <fstream>
#include <tclap/CmdLine.h>

#include "GLPKL1Minimisation.h"
#include "stabiliser.h"

using namespace std;

/* Program for computing the minimum l1 norm ||x||_1 over decompositions of a point y = sum_i x_i v_i in vertices v_i.

	input to the program:

	* constraint matrix consisting of the coordinates of the polytope vertices, specified in some data file

*/



int main(int argc, char** argv) {

// ----------------------------------
// ----- parse comand-line parameters
// ----------------------------------

string infile,outfile,state;
bool Quiet;

try {

	TCLAP::CmdLine cmd("Program for solving the membership problem for the stabiliser polytope", ' ', "0.1");

	// arguments

	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file that contains the coordinates of the vertices as a sparse matrix with rows of the form i j a[i,j]", true, " ", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Name that will be used for output files", false, "solution.dat", "string");
	cmd.add(output_arg);

	TCLAP::SwitchArg Quiet_arg ("Q","Quiet","Suppress detailed output to files", cmd, false);

	vector<string> allowed_states = { "T", "H" };
	TCLAP::ValuesConstraint<string> state_con ( allowed_states );
	TCLAP::ValueArg<string> state_arg ("s", "state", "Which state to use: Either H or T", true, "H", &state_con);
	cmd.add(state_arg);

	cmd.parse(argc, argv);

	infile = input_arg.getValue();
	outfile = output_arg.getValue();
	Quiet = Quiet_arg.getValue();
	state = state_arg.getValue();

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}

// ----------------------------------
// ----- prepare computation
// ----------------------------------


cout << "#-------------------------------------------------" << endl;
cout << "# l1-minimisation over projected stabiliser states" << endl;
cout << "# for tensor powers of the " << state << " state" << endl;
cout << "#-------------------------------------------------" << endl;



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

unsigned nqubits = states.at(0).size();

// cout << nqubits << endl;
GLPKL1Minimisation lp (states);
lp.set_verbosity(3);
lp.print_parameters();

// ----------------------------------
// ------ start computation
// ----------------------------------


// solve & write solution
vector<double> y;
if(state == "H") {
	y = H_state_nb(nqubits);
}
else {
	y = T_state_nb(nqubits);
}

int ret_status = lp.check_point(y);
double ROM = lp.get_obj_value();

cout << "ROM(" << state << "^" << nqubits << ") = " << scientific << ROM << endl;

// write solution
lp.write_glpk_output(outfile);
lp.write_sol(outfile+".sol");

}