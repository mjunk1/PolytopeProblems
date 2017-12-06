#include <iostream>
#include <fstream>
#include <tclap/CmdLine.h>
#include <algorithm>

#include "symplectic.h"
#include "utilities.h"

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

	TCLAP::ValueArg<string> output_arg ("o", "outfile", "Output file name that will be used for writing the reduced constraint matrix", false, "out.dat", "string");
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

cout << "Generate the symplectic group for n = " << n << " ..." << endl;
cout << "Compute projections of Clifford group ... " << endl;

// vectors to store the symplectic matrix and Liouville representation
// we will store all Liouville at once to be able to delete duplicates
unsigned spn = sp_order(n);
unsigned paulin = pow(4,n);
vector<unsigned> S(2*n);
vector<vector<int>> L (paulin*spn, vector<int>(pow(4,n+1)-1) ); 

for(unsigned i=0; i<spn; i++) {
	// generate i-th symplectic matrix
	S = generate_symplectic_matrix(i,n);

	// compute the 1-qubit Pauli orbit and save the relevant entries of the Liouville matrix
	for(unsigned j=0; j<paulin; j++) {
		L.at(paulin*i+j) = projected_liouville_matrix(S,j);		
	}	
}

// sort and delete duplicates
sort(L.begin(),L.end());
auto last = unique(L.begin(),L.end());
int d = distance(last,L.end());
cout << "Deleted " << d << " duplicates. " << spn*paulin-d << " elements left." << endl;
L.erase(last, L.end());

// save L to file
fstream fout(outfile, ios::out);
unsigned j = 0;

for(auto l : L) {
	for(unsigned k=0; k<l.size(); k++) {
		if(l.at(k) != 0) {
			fout << j+1 << " " << k+1 << " " << l.at(k) << endl;
		}
	}
	++j;
}

fout.close();


}