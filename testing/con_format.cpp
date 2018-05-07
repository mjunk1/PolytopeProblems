#include <iostream>
#include <fstream>
#include <vector>
#include <tclap/CmdLine.h>

#include "utilities.h"

using namespace std;

// ----- Test of stabiliser generation using std::set in symplectic.h

int main(int argc, char** argv) {

	string infile;
	string outfile;

	try {
		TCLAP::CmdLine cmd("Program for converting dense matrix into sparse matrix of COO format", ' ', "0.1");

		// arguments

		TCLAP::ValueArg<string> input_arg ("i", "in", "Input file that contains the rows of the dense matrix", true, " ", "string");
		cmd.add(input_arg);

		TCLAP::ValueArg<string> output_arg ("o", "out", "Output file where the program write the sparse matrix to. The rows will have the form i j a[i,j]", true, " ", "string");
		cmd.add(output_arg);

		cmd.parse(argc, argv);

		infile = input_arg.getValue();
		outfile = output_arg.getValue();

	} catch (TCLAP::ArgException &e) { 
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl; 
	}


	unsigned nrows = get_number_of_lines(infile);
	double val;
	vector<double> matrix;

	fstream fin (infile, ios::in);
	if(fin.is_open()) {
		while(fin >> val) {
			matrix.push_back(val);
		}
	}
	fin.close();
	
	unsigned ncols = matrix.size()/nrows;

	fstream fout (outfile, ios::out);

	if(fout.is_open()) {
		for(unsigned i=0; i<nrows; i++) {
			for(unsigned j=1; j<ncols; j++) {
				if(matrix.at(i*ncols+j) != 0) {
					fout << i+1 << " " << j << " " << matrix.at(i*ncols+j) << endl;;
				}
			}
		}
	}
	fout.close();
}