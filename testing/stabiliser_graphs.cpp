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
} 
else {
	cout << "Couldn't open file " << infile << endl;
}
fin.close();

cout << "Read " << graphs.size() << " graphs." << endl;

// convert graph6 to adjacency matrices and back
vector<binvec> adj_mat;
unsigned n = get_order(graphs.at(0));
for(auto g : graphs) {
	adj_mat.push_back( graph6_to_adj_mat(g) );
}

// print
cout << "Check if conversion works ... " << endl;
for(unsigned k=0; k<adj_mat.size(); k++) {
	if(adj_mat_to_graph6(adj_mat.at(k), n) != graphs.at(k)) {
		cout << "Fail at graph #" << k << endl;
	}
}
cout << "Done" << endl;



}