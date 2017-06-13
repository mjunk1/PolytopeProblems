#include <iostream>
#include <glpk.h>
#include <fstream>
#include <cmath>
#include <tclap/CmdLine.h>
using namespace std;

/* input to the program:

	* number of qubits
	* length of the circuit
	* density of CNOTs
	* constraint matrix consisting of the coordinates of the polytope vertices, specified in some data file

*/

/* structure of the program:
	(1) read parameters and constraint matrix etc.
	(2) create a random Clifford circuit								\
		(2a) propagate noise through the circuit						|	optional (later) or Mathematica?
		(2b) compute adjoint representation of final noisy T channel	/
	(3) for this channel, solve the membership problem as a function of the individual noise strength p
	(4) output data 
	(5) restart at (2)
*/

/* functions to write
	* channel representation as a function of final prob. distribution
	* random Clifford circuit generation
		- matrix-vector operations mod 2?
		- sparsity
*/

int get_number_of_lines(string filename) {
	int n;
	string line;
	fstream fin(filename, ios::in);
	if(fin.is_open()) 
		while(getline(fin, line))
			++n;
	fin.close();
	return n;
}

int get_dc_rep_1Q(double p, double *y) {
	y[1] = 1 - (4.*p)/3.;
	y[2] = 0;
	y[3] = 0; 
	y[4] = 0; 
	y[5] = (1 - p)/sqrt(2) - p/(3.*sqrt(2));
	y[6] = -((1 - p)/sqrt(2)) + p/(3.*sqrt(2));
	y[7] = 0;
	y[8] = (1 - p)/sqrt(2) - p/(3.*sqrt(2));
	y[9] = (1 - p)/sqrt(2) - p/(3.*sqrt(2));
}


int get_pauli_rep_1Q(double *p, double *y) {
	y[1] = p[1] - p[2] - p[3] + p[4];
	y[2] = y[3] = y[4] = 0;
	y[5] = p[1]/sqrt(2) + p[2]/sqrt(2) - p[3]/sqrt(2) - p[4]/sqrt(2);
	y[6] = -(p[1]/sqrt(2)) + p[2]/sqrt(2) - p[3]/sqrt(2) + p[4]/sqrt(2);
	y[7] = 0;
	y[8] = p[1]/sqrt(2) + p[2]/sqrt(2) - p[3]/sqrt(2) - p[4]/sqrt(2);
	y[9] = p[1]/sqrt(2) - p[2]/sqrt(2) + p[3]/sqrt(2) - p[4]/sqrt(2);
}

int convolve_mod2(double *a, double *b, double *c, int n) {
	for(int i=0; i<n; i++) {
		c[i] = 0;
		for(int j=0; j<n; j++) {
			c[i] += a[j] * b[i^j];
		} 
	}
	return 0;
}


int main(int argc, char** argv) {

// ----- parse comand-line parameters
int number_of_qubits;
string cmatrix_file;
string output_prefix;

try {

	TCLAP::CmdLine cmd("Program for solving the membership problem for the Clifford polytope", ' ', "0.1");

	// arguments
	vector<int> allowed = {1,2};
	TCLAP::ValuesConstraint<int> constraint(allowed);
	TCLAP::ValueArg<int> nqubits_arg ("q","qubits", "Number of qubits", true, 0, &constraint);
	cmd.add(nqubits_arg);

	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file that contains the coordinates of the vertices as a sparse matrix with rows of the form i j a[i,j]", true, " ", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "out_prefix", "Prefix that will be used for output files", false, "data/polytope_", "string");
	cmd.add(output_arg);

	cmd.parse(argc, argv);

	number_of_qubits = nqubits_arg.getValue();
	cmatrix_file = input_arg.getValue();
	output_prefix = output_arg.getValue();

} catch (TCLAP::ArgException &e)  // catch any exceptions
{ cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }

int dimension = pow(pow(2, 2*number_of_qubits)-1, 2);
int number_of_vertices = 24;

// int circuit_length;
// double cnot_density;



// ----- read sparse constraint matrix

// GLPK convention: data start at index 1
int number_of_entries = get_number_of_lines(cmatrix_file);
int *ia = new int[number_of_entries+1];
int *ja = new int[number_of_entries+1];
double *a = new double[number_of_entries+1];
int n = 1, i, j;
ia[0] = ja[0] = 0;
a[0] = 0;

// open file and read
fstream fin(cmatrix_file, ios::in);
if(fin.is_open()) 
	while(fin >> ia[n] >> ja[n] >> a[n])
		++n;

fin.close();

// ----- create GLPK problem 

glp_prob *lp = glp_create_prob();

// coordinates of the point to check
double *y = new double[dimension+2];
y[0] = 0.;
y[10] = -1.;

// set up problem parameters
glp_set_prob_name(lp, "Polytope membership");
glp_set_obj_dir(lp, GLP_MAX);

char s[5];

// setting up rows
glp_add_rows(lp, number_of_vertices+1);
for(i=1; i<=number_of_vertices; i++) {
	sprintf(s, "p%u", i);
	glp_set_row_name(lp, i, s);
	glp_set_row_bnds(lp, i, GLP_UP, 0.0, 0.0);
}
glp_set_row_name(lp, number_of_vertices+1, "q");
glp_set_row_bnds(lp, number_of_vertices+1, GLP_UP, 0.0, 1.0);

// setting up cols 
glp_add_cols(lp, dimension+1);
for(j=1; j<=dimension; j++) {
	sprintf(s, "x%u", j);
	glp_set_col_name(lp, j, s);
	glp_set_col_bnds(lp, j, GLP_FR, 0.0, 0.0);
}
sprintf(s, "x0");
glp_set_col_name(lp, dimension+1, s);
glp_set_col_bnds(lp, dimension+1, GLP_FR, 0.0, 0.0);

// setting up constraint matrix
glp_load_matrix(lp, number_of_entries, ia, ja, a);


// ----- start computation
int glp_ret;

// possible loops:
//  * random circuits
//  * noise strength 

// noise strength loop

int *ind = new int[dimension+2];
fin.open(output_prefix + "obj.dat", ios::out);
if(!fin.is_open()) {
	cerr << "Couldn't open output file " << output_prefix + "obj.dat" << " . Will now hold." << endl;
	return 1;
}

for(double p=0.; p<=1; p+=0.01) {

	// compute y
	get_dc_rep_1Q(p, y);

	// replacing the last row of the constraint matrix
	ind[0] = 0;
	for(i=1; i<=dimension+1; i++)
		ind[i] = i;
	glp_set_mat_row(lp, number_of_vertices+1, dimension+1, ind, y);

	// setting up objective function
	for(j=1; j<=dimension+1; j++) 
		glp_set_obj_coef(lp, j, y[j]);

	// solve
	glp_ret = glp_simplex(lp, NULL);

	// write output
	glp_print_sol(lp, (output_prefix + to_string(p) + ".dat").c_str());
	fin << p << " " << glp_get_obj_val(lp) << endl;
}

fin.close();

// free stuff
delete[] ia, ja, a;
delete[] y;
delete[] ind;
glp_delete_prob(lp);

}