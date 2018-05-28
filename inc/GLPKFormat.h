#ifndef GLPKFORMAT_H
#define GLPKFORMAT_H
#include <fstream>
#include <iostream>
#include <vector>


#ifndef UTILITIES_H
#include "utilities.h"
#endif

using namespace std;

// ----- GLPK utilities

// data struct for GLPK
struct GLPKFormat {
	vector<int> rows;
	vector<int> cols;
	vector<double> values;
	unsigned non_zeros;
	unsigned nrows;
	unsigned ncols;
};

// reads sparse matrix in COO format and saves to struct
GLPKFormat to_GLPK_format(string cmatrix_file, bool transpose=false) {
	GLPKFormat ret;

	// number of non-zero values
	ret.non_zeros = get_number_of_lines(cmatrix_file);

	// GLPK convention: data start at index 1
	ret.rows = vector<int> (ret.non_zeros+1,0);
	ret.cols = vector<int> (ret.non_zeros+1,0);
	ret.values = vector<double> (ret.non_zeros+1,0.);

	if(transpose == true) {
		// open file and read
		fstream fin(cmatrix_file, ios::in);
		int n = 1;
		if(fin.is_open()) {
			while(fin >> ret.cols[n] >> ret.rows[n] >> ret.values[n]) {
				++n;
			}
		} else {
			cerr << "Error: Couldn't open file " << cmatrix_file << endl;
		}
		fin.close();

		// save number of rows and columns
		ret.nrows = *max_element(ret.rows.begin(),ret.rows.end());
		ret.ncols = *max_element(ret.cols.begin(),ret.cols.end());
	}
	else {
		// open file and read
		fstream fin(cmatrix_file, ios::in);
		int n = 1;
		if(fin.is_open()) {
			while(fin >> ret.rows[n] >> ret.cols[n] >> ret.values[n]) {
				++n;
			}
		} else {
			cerr << "Error: Couldn't open file " << cmatrix_file << endl;
		}
		fin.close();

		// save number of rows and columns
		ret.nrows = ret.rows.at(ret.non_zeros);
		ret.ncols = *max_element(ret.cols.begin(),ret.cols.end());
	}

	return ret;
}

// reads dense list of row vectors and saves to struct
GLPKFormat to_GLPK_format(vector<vector<int>> matrix) {
	GLPKFormat ret;

	// save number of rows and columns
	ret.nrows = matrix.size();
	ret.ncols = matrix.at(0).size();

	// number of non-zero values
	ret.non_zeros = ret.nrows * ret.ncols;

	// GLPK convention: data start at index 1
	ret.rows = vector<int> (ret.non_zeros+1,0);
	ret.cols = vector<int> (ret.non_zeros+1,0);
	ret.values = vector<double> (ret.non_zeros+1,0.);

	// convert data
	for(unsigned i=0; i<ret.nrows; i++) {
		for(unsigned j=0; j<ret.ncols; j++) {
			ret.rows.at( ret.ncols*i + j + 1 ) = i+1;
			ret.cols.at( ret.ncols*i + j + 1 ) = j+1;
			ret.values.at( ret.ncols*i + j + 1 ) = (double)(matrix.at(i).at(j));
		}
	}	

	return ret;
}

GLPKFormat to_GLPK_format(vector<LabelledPoint<int>> matrix) {
	GLPKFormat ret;

	// save number of rows and columns
	ret.nrows = matrix.size();
	ret.ncols = matrix.at(0).size();

	// number of non-zero values
	ret.non_zeros = ret.nrows * ret.ncols;

	// GLPK convention: data start at index 1
	ret.rows = vector<int> (ret.non_zeros+1,0);
	ret.cols = vector<int> (ret.non_zeros+1,0);
	ret.values = vector<double> (ret.non_zeros+1,0.);

	// convert data
	for(unsigned i=0; i<ret.nrows; i++) {
		for(unsigned j=0; j<ret.ncols; j++) {
			ret.rows.at( ret.ncols*i + j + 1 ) = i+1;
			ret.cols.at( ret.ncols*i + j + 1 ) = j+1;
			ret.values.at( ret.ncols*i + j + 1 ) = (double)(matrix.at(i).point.at(j));
		}
	}	

	return ret;
}

void print(GLPKFormat& data) {
	// GLPK = 1-based indexing
	unsigned k=1;
	for(unsigned i=1; i<=data.nrows; i++) {
		for(unsigned j=1; j<=data.ncols; j++) {
			if( (data.rows.at(k) == i) && (data.cols.at(k) == j) ) {
				cout << data.values.at(k) << " ";
				++k;
			}
			else{ 
				cout << 0 << " ";
			}
		}
		cout << endl;
	}
}

#endif