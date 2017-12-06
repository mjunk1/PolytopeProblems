#include <iostream>
#include <glpk.h>
#include <fstream>
#include <map>
#include <tclap/CmdLine.h>

#include "ConvexSeperation.h"
#include "noise.h"

using namespace std;

/* input to the program:

	* number of qubits
	* constraint matrix consisting of the coordinates of the polytope vertices, specified in some data file

*/




int main(int argc, char** argv) {

// ----------------------------------
// ----- parse comand-line parameters
// ----------------------------------

string cmatrix_file;
string outfile;
string method;
bool Quiet;

int dimension; // ambient dimension of the embedded polytope
int number_of_vertices;


try {

	TCLAP::CmdLine cmd("Program for solving the membership problem for the Clifford polytope", ' ', "0.1");

	// arguments

	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file that contains the coordinates of the vertices as a sparse matrix with rows of the form i j a[i,j]", true, " ", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Output file name that will be used for writing the reduced constraint matrix", false, "out.dat", "string");
	cmd.add(output_arg);

	TCLAP::ValueArg<int> vertices_arg ("v", "vertices", "Number of vertices of the polytope", false, 0, "Positive integer. By default, this will be deduced from the input file.");
	cmd.add(vertices_arg);

	TCLAP::ValueArg<int> dim_arg ("d", "dimension", "Ambient dimension of the polytope", true, 0, "Positive integer. By default, this will be set to 4^(nqbits+1)-1.");
	cmd.add(dim_arg);

	TCLAP::SwitchArg Quiet_arg ("Q","Quiet","Suppress detailed output to files", cmd, false);

	vector<string> allowed2 = {"simplex", "interior_point"};
	TCLAP::ValuesConstraint<string> constraint2(allowed2);
	TCLAP::ValueArg<string> method_arg ("m","method", "LP solver to use", false, "simplex", &constraint2);
	cmd.add(method_arg);

	cmd.parse(argc, argv);

	cmatrix_file = input_arg.getValue();
	outfile = output_arg.getValue();
	number_of_vertices = vertices_arg.getValue();
	dimension = dim_arg.getValue();
	method = method_arg.getValue();
	Quiet = Quiet_arg.getValue();


	if(number_of_vertices < 0) {
		cerr << "Error: Number of vertices is smaller than 0." << endl;
		return 1;
	} 


} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}

// ----------------------------------
// ----- prepare computation
// ----------------------------------

GLPKConvexSeperation lp (dimension, cmatrix_file);
lp.set_method(method);
int ret_status;

cout << "Initital number of vertices: " << lp.get_nvertices() << endl;

cout << "Number of non-redundant vertices: " << lp.delete_redundant_points() << endl;

lp.write_constraint_matrix(outfile);


// ----------------------------------
// ------ start computation
// ----------------------------------

// solve & write solution
// double pth = lp.check_family(y);
// cout << "Threshold value p_th = " << pth << endl;


}