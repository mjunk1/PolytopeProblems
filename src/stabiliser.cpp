#include <iostream>
#include <fstream>
#include <tclap/CmdLine.h>

#include "GLPKConvexSeparation.h"
#include "stabiliser.h"

using namespace std;

/* input to the program:

	* number of qubits
	* constraint matrix consisting of the coordinates of the polytope vertices, specified in some data file

*/

/* structure of the program:
	(1) read parameters, constraint matrix and set up solver
	(2) start iteration base on interval division method
		(i) generate point to check y = y(p)
		(ii) update p depending on the result
		(iii) restart at (i) until converged to the threshold value p = p_th 
	(3) output results
*/



int main(int argc, char** argv) {

// ----------------------------------
// ----- parse comand-line parameters
// ----------------------------------

string cmatrix_file;
string output_prefix;
double nlower, nupper, eps;
unsigned niter;
string method;
bool Quiet;

int dimension; // = pow(4, number_of_qubits + 1) - 1; // ambient dimension of the embedded polytope

try {

	TCLAP::CmdLine cmd("Program for solving the membership problem for the stabiliser polytope", ' ', "0.1");

	// arguments

	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file that contains the coordinates of the vertices as a sparse matrix with rows of the form i j a[i,j]", true, " ", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "out_prefix", "Prefix that will be used for output files", false, "out_", "string");
	cmd.add(output_arg);

	TCLAP::ValueArg<double> cnum1_arg ("n","noise", "Lower bound for the noise", false, 0.0, "Positive number between 0 and 1");
	cmd.add(cnum1_arg);

	TCLAP::ValueArg<double> cnum2_arg ("N","Noise", "Upper bound for the noise", false, 1.0, "Positive number between 0 and 1");
	cmd.add(cnum2_arg);

	TCLAP::ValueArg<double> eps_arg ("e","Error", "Error for the convergence of the noise threshold", false, 1e-6, "Positive number");
	cmd.add(eps_arg);

	TCLAP::ValueArg<unsigned> citer_arg ("i", "iterations", "Number of circuits the generate", false, 10, "Positive integer");
	cmd.add(citer_arg);

	TCLAP::SwitchArg Quiet_arg ("Q","Quiet","Suppress detailed output to files", cmd, false);

	vector<string> allowed2 = {"simplex", "interior_point"};
	TCLAP::ValuesConstraint<string> constraint2(allowed2);
	TCLAP::ValueArg<string> method_arg ("m","method", "LP solver to use", false, "simplex", &constraint2);
	cmd.add(method_arg);

	cmd.parse(argc, argv);

	cmatrix_file = input_arg.getValue();
	output_prefix = output_arg.getValue();
	nlower = cnum1_arg.getValue();
	nupper = cnum2_arg.getValue();
	niter = citer_arg.getValue();
	eps = eps_arg.getValue();
	method = method_arg.getValue();
	Quiet = Quiet_arg.getValue();


	if(nlower > nupper) {
		cerr << "Error: Lower bound for the noise strength is greater than the upper bound." << endl;
		return 1;
	}

	if(niter < 1) {
		cerr << "Error: Number of circuits to generate is less than 1." << endl;
		return 1;
	}

	if(eps <= 0) {
		cerr << "Error: Convergence error is less than or equal to 0." << endl;
		return 1;
	} 


} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}

// ----------------------------------
// ----- prepare computation
// ----------------------------------

GLPKConvexSeparation lp (cmatrix_file);
lp.set_method(method);
lp.set_verbosity(1);
lp.set_parameters(nlower,nupper,50,eps);
int ret_status;

// print parameters
// lp.print_parameters();


// set up point to check
NoisyHState y (lp.get_dimension());


// ----------------------------------
// ------ start computation
// ----------------------------------

// solve & write solution
double pth = lp.check_family(y);
cout << "Threshold value p_th = " << setprecision((int)(-log10(eps))) << pth << endl;


}