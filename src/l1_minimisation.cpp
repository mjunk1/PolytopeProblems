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

string cmatrix_file;
string output_prefix;
bool Quiet;
unsigned nqubits;
bool project;


int dimension; 

try {

	TCLAP::CmdLine cmd("Program for solving the membership problem for the stabiliser polytope", ' ', "0.1");

	// arguments

	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file that contains the coordinates of the vertices as a sparse matrix with rows of the form i j a[i,j]", true, " ", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "out_prefix", "Prefix that will be used for output files", false, "out_", "string");
	cmd.add(output_arg);


	TCLAP::SwitchArg Quiet_arg ("Q","Quiet","Suppress detailed output to files", cmd, false);

	TCLAP::ValueArg<bool> project_arg ("p","project", "Boolean variable which decides whether to project the states or not.", false, true, "Bool");
	cmd.add(project_arg);

	cmd.parse(argc, argv);

	cmatrix_file = input_arg.getValue();
	output_prefix = output_arg.getValue();
	Quiet = Quiet_arg.getValue();
	project = project_arg.getValue();

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}

// ----------------------------------
// ----- prepare computation
// ----------------------------------

GLPKL1Minimisation lp (cmatrix_file);
lp.set_verbosity(2);
lp.print_parameters();

// ----------------------------------
// ------ start computation
// ----------------------------------


// solve & write solution
vector<double> y;
if(project == true) {
	nqubits = lp.get_dimension();
	y = H_state_nb(nqubits);
}
else {
	nqubits = (unsigned)(log(lp.get_dimension())/log(4));
	y = H_state(nqubits);
}

int ret_status = lp.check_point(y);
double ROM = lp.get_obj_value();



cout << "ROM(H^" << nqubits << ")^(1/" << nqubits <<  ") = " << scientific << pow(ROM,1/(double)nqubits) << endl;


}