#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <set>
#include <chrono>
#include <tclap/CmdLine.h>


#include "symplectic.h"
#include "utilities.h"
#include "stabiliser.h"


using namespace std;

// ----- Test of stabiliser generation using std::set in symplectic.h

int main(int argc, char** argv) {

string infile;

try {

	TCLAP::CmdLine cmd("Test blabla", ' ', "0.1");

	// arguments
	TCLAP::ValueArg<string> input_arg ("f", "file", "Input file name that contains the adjacency matrices for the non-isomorphic graphs with q vertices", true, "in.dat", "string");
	cmd.add(input_arg);

	cmd.parse(argc, argv);

	infile = input_arg.getValue();

} catch (TCLAP::ArgException &e) { 
	cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
}


cout << "Test graph stuff" << endl;
cout << "-----------------" << endl;


// open the file
fstream fin(infile, ios::in);

vector<string> graphs;
string row;

if(fin.is_open()) {
	while(fin >> row) {
		graphs.push_back(row);
	}

	cout << "Read " << graphs.size() << " graphs." << endl;

	// convert graph6 to adjacency matrices and back
	unsigned n = get_order(graphs.at(0));

	// print
	cout << "Check if conversion works ... " << endl;
	vector<binvec> adj_mat;
	for(unsigned k=0; k<graphs.size(); k++) {
		adj_mat = graph6_to_adj_mat(graphs.at(k));
		if(adj_mat_to_graph6(adj_mat) != graphs.at(k)) {
			cout << "Fail at graph #" << k << endl;
			cout << "Got " << adj_mat_to_graph6(adj_mat) << " instead of " << graphs.at(k) << endl;
		}

	}
	cout << "Done" << endl;

	// now check product labeling
	string label1 = graphs.at(0)+" ";
	for(unsigned i=0; i<n; i++)
		label1 += "0";
	label1 += " ";
	for(unsigned i=0; i<n; i++)
		label1 += "0";

	string label2 = graphs.at(1)+" ";
	for(unsigned i=0; i<n; i++)
		label2 += "0";
	label2 += " ";
	for(unsigned i=0; i<n; i++)
		label2 += "0";

	string label = pr_product_label(label1, label2);

	cout << "Product of " << endl;
	cout << label1 << endl;
	cout << "and" << endl;
	cout << label2 << endl;
	cout << "is" << endl;
	cout << label << endl;

} 
else {
	cout << "Couldn't open file " << infile << endl;
}
fin.close();


}