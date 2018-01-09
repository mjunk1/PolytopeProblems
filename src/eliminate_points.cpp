#include <iostream>
#include <fstream>
#include <tclap/CmdLine.h>

#include "GLPKConvexSeparation.h"

using namespace std;


int main(int argc, char** argv) {

// ----------------------------------
// ----- parse comand-line parameters
// ----------------------------------

string cmatrix_file;
string outfile;
string method;
bool Quiet;


try {

	TCLAP::CmdLine cmd("Program for eliminating redundant points from the convex hull of a set of points.", ' ', "0.1");

	// arguments

	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file that contains the coordinates of the vertices as a sparse matrix with rows of the form i j a[i,j]", true, " ", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Output file name that will be used for writing the reduced constraint matrix", false, "out.dat", "string");
	cmd.add(output_arg);


	TCLAP::SwitchArg Quiet_arg ("Q","Quiet","Suppress detailed output to files", cmd, false);

	vector<string> allowed2 = {"simplex", "interior_point"};
	TCLAP::ValuesConstraint<string> constraint2(allowed2);
	TCLAP::ValueArg<string> method_arg ("m","method", "LP solver to use", false, "simplex", &constraint2);
	cmd.add(method_arg);

	cmd.parse(argc, argv);

	cmatrix_file = input_arg.getValue();
	outfile = output_arg.getValue();
	method = method_arg.getValue();
	Quiet = Quiet_arg.getValue();



} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}

// ----------------------------------
// ----- elimination
// ----------------------------------

GLPKConvexSeparation lp (cmatrix_file);
lp.set_method(method);
lp.set_verbosity(1);
int ret_status;

lp.print_parameters();
unsigned nvertices = lp.get_nvertices();
lp.delete_redundant_points();
cout << "Deleted " << nvertices - lp.get_nvertices() << " points." << endl;
lp.print_parameters();

lp.write_constraint_matrix(outfile);



}