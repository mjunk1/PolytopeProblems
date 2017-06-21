#include <iostream>
#include <glpk.h>
#include <fstream>
#include <random>
#include <map>
#include <tclap/CmdLine.h>
#include "utilities.h"

#ifndef POLYTOPEPROGRAM_H
#include "PolytopeProgram.h"
#endif

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
		(2b) compute adjoint representation of foutal noisy T channel	/
	(3) for this channel, solve the membership problem as a function of the individual noise strength p
	(4) output data 
	(5) restart at (2)
*/

int main(int argc, char** argv) {

// ----- parse comand-line parameters
unsigned number_of_qubits;
string cmatrix_file;
string output_prefix;
vector<string> circuit;
unsigned circuit_length;
double cnot_prob;
double nlower, nupper, ndelta;
unsigned niter;

try {

	TCLAP::CmdLine cmd("Program for solving the membership problem for the Clifford polytope", ' ', "0.1");

	// arguments
	vector<unsigned> allowed = {1,2};
	TCLAP::ValuesConstraint<unsigned> constraint(allowed);
	TCLAP::ValueArg<unsigned> nqubits_arg ("q","qubits", "Number of qubits", true, 0, &constraint);
	cmd.add(nqubits_arg);

	TCLAP::ValueArg<unsigned> clength_arg ("l","length", "Circuit length", false, 4, "positive number");
	cmd.add(clength_arg);

	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file that contains the coordinates of the vertices as a sparse matrix with rows of the form i j a[i,j]", true, " ", "string");
	cmd.add(input_arg);

	TCLAP::ValueArg<string> output_arg ("o", "out_prefix", "Prefix that will be used for output files", false, "out_", "string");
	cmd.add(output_arg);

	TCLAP::MultiArg<string> circuit_arg ("c", "circuit", "Specifies a specific order of gates that form a circuit", false, "string");
	cmd.add(circuit_arg);

	TCLAP::ValueArg<double> pc_arg ("p","prob", "CNOT probability in the circuit", false, 1./3., "positive number between 0 and 1");
	cmd.add(pc_arg);

	TCLAP::ValueArg<double> cnum1_arg ("n","noise", "Lower bound for the noise", false, 0.0, "Positive number between 0 and 1");
	cmd.add(cnum1_arg);

	TCLAP::ValueArg<double> cnum2_arg ("N","Noise", "Upper bound for the noise", false, 1.0, "Positive number between 0 and 1");
	cmd.add(cnum2_arg);

	TCLAP::ValueArg<double> cnumd_arg ("d","Delta", "Delta / stepsize for noise strength", false, 0.01, "Positive number between 0 and 1");
	cmd.add(cnumd_arg);

	TCLAP::ValueArg<unsigned> citer_arg ("i", "iterations", "Number of circuits the generate", false, 10, "Positive integer");
	cmd.add(citer_arg);

	cmd.parse(argc, argv);

	number_of_qubits = nqubits_arg.getValue();
	cmatrix_file = input_arg.getValue();
	output_prefix = output_arg.getValue();
	circuit_length = clength_arg.getValue();
	circuit = circuit_arg.getValue();
	cnot_prob = pc_arg.getValue();
	nlower = cnum1_arg.getValue();
	nupper = cnum2_arg.getValue();
	niter = citer_arg.getValue();
	ndelta = cnumd_arg.getValue();

	if(circuit.size() != 0) {
		circuit_length = circuit.size();
		niter = 1;
	}

	if(circuit_length < 1) {
		cerr << "Error: Circuit length is less than 1." << endl;
		return 1;
	}

	if(cnot_prob < 0 || cnot_prob > 1) {
		cerr << "Error: CNOT probability is less than 0 or greater than 1." << endl;
		return 1;
	}

	if(nlower > nupper) {
		cerr << "Error: Lower bound for the noise strength is greater than the upper bound." << endl;
		return 1;
	}

	if(niter < 1) {
		cerr << "Error: Number of circuits to generate is less than 1." << endl;
		return 1;
	}

	if(nupper - nlower < ndelta) {
		cerr << "Error: Stepsize for noise strength is greater than the specified interval." << endl;
		return 1;
	}

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}

int dimension = pow(pow(2, 2*number_of_qubits)-1, 2);
int number_of_vertices;
int len = pow(4,number_of_qubits);

if(number_of_qubits == 1)
	number_of_vertices = 24;
else
	number_of_vertices = 11520;



// ----- start computation

PolytopeProgram lp (number_of_vertices, dimension, cmatrix_file);
int ret_status;
string out;
fstream fout,fout2;

// probability vectors
double *dn = new double[len];
double *pdist0 = new double[len];
double *pdist1 = new double[len];

// coordinates of the noisy channel
double *y = new double[dimension];

// --- randomized circuit generation 

// symplectic representations of elementary channels
// unsigned H1[4] = { 0b1010, 0b0101, 0b0010, 0b0001 };
// unsigned H2[4] = { 0b1100, 0b0100, 0b0011, 0b0001 };
// unsigned S1[4] = { 0b0010, 0b0001, 0b1000, 0b0100 };
// unsigned S2[4] = { 0b0100, 0b1000, 0b0001, 0b0010 };
unsigned H1[4] = { 0b1100, 0b0100, 0b0010, 0b0001 };
unsigned H2[4] = { 0b1000, 0b0100, 0b0011, 0b0001 };
unsigned S1[4] = { 0b0100, 0b1000, 0b0010, 0b0001 };
unsigned S2[4] = { 0b1000, 0b0100, 0b0001, 0b0010 };
unsigned CNOT[4] = { 0b1010, 0b0110, 0b0010, 0b1101 };
unsigned ID[4] = { 0b1000, 0b0100, 0b0010, 0b0001 };

// making the gate access a bit easier
vector<unsigned*> gates = { H1, H2, S1, S2, CNOT, ID };
map<string, unsigned> gate_id = { {"H1", 0}, {"H2", 1}, {"S1", 2}, {"S2", 3}, {"CNOT", 4}, {"ID", 5} };
vector<unsigned> gate_order (circuit_length, 0);

// Prepare random circuit generation
// We use a Marsenne-Twister generator which is seeded by a random_device 
// cnot_prob is used to generate the distribution. It is assumed that H and S gates occur with the same probability
double hs_prob = (1.-cnot_prob)/4.;
random_device rd; 
mt19937 rd_gen(rd()); 
discrete_distribution<> rd_dist( { hs_prob, hs_prob, hs_prob, hs_prob, cnot_prob } ); 

// --- circuit loop

for(unsigned c=0; c < niter; c++) {

	if(circuit.size() != 0) {
		// this means an explicit circuit has been set
		for(unsigned i=0; i<circuit_length; i++) 
			gate_order.at(i) = gate_id[circuit.at(i)];

	} else {
		// generate a circuit of given length
		cout << "Generating circuit of length " << circuit_length << endl;
		cout << "---------------------------------------------------" << endl << endl;

		for(unsigned i=0; i<circuit_length; i++)			
			gate_order.at(i) = rd_dist(rd_gen);
	}


	// --- noise strength loop

	// the order of gates serves as identifier for the output file
	out = output_prefix;
	for(unsigned i : gate_order) {
		out += to_string(i);
		fout << i << " ";
	}
	fout << endl;

	cout << "Circuit file: " << out << endl << endl;

	// we will write the optimal value of the objective function for every value of the noise strength to this file
	fout.open(out + "_obj.dat", ios::out);
	if(!fout.is_open()) {
		cerr << "Couldn't open output file " << out + "_obj.dat" << " . Will now hold." << endl;
		return 1;
	}

	fout2.open(out + "_fdist.dat", ios::out);
	if(!fout2.is_open()) {
		cerr << "Couldn't open output file " << out + "_fdist.dat" << " . Will now hold." << endl;
		return 1;
	}

	for(double p=nlower; p<=nupper; p+=ndelta) {
		cout << "Starting optimization for noise strength p = " << p << endl;
		cout << "---------------------------------------------------" << endl << endl;

		// init pvec
		if(number_of_qubits == 1) {
			dn_pdist_1Q(p, dn);
			noisyT_1Q(dn, y);
		}
		else if(number_of_qubits == 2) {
			dn_pdist_2Q(p, dn);
			copy(dn, dn+len, pdist0);

			// write_pdist(out + "_" + to_string(p) + "_dist.dat", pdist0, number_of_qubits);

			// propagate through circuit
			for(unsigned i=0; i<circuit_length; i++) {
				symplectic_transform(pdist0, pdist1, gates[gate_order[i]], number_of_qubits);
				convolve_mod2(pdist1, dn, pdist0, 2*number_of_qubits);
			}

			fout2 << "Final distribution for p = " << p << endl;
			write_pdist(fout2, pdist0, number_of_qubits);

			// compute adjoint representation of the foutal channel
			noisyT_2Q(pdist0, y);
		}
		else {
			break;
		}

		
		// solve & write solution
		ret_status = lp.check_point(y);
		//lp.write_sol(output_prefix + to_string(p) + ".dat");
		fout << p << " " << lp.get_result() << " " << ret_status << " " << lp.get_status() << endl;

		cout << endl;
		cout << "---------------------------------------------------" << endl << endl;
	}

	fout.close();
	fout2.close();

}
// --- end of circuit loop


// ------ free stuff
delete[] y;
delete[] dn, pdist0, pdist1;

}