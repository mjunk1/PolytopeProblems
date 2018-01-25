#ifndef GLPKFORMAT_H
#define GLPKFORMAT_H
#include <fstream>
#include <iostream>
#include <vector>
#include <eigen3/Eigen/Sparse>


#ifndef UTILITIES_H
#include "utilities.h"
#endif

using namespace std;
using namespace Eigen;


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

// converts a Eigen::SparseMatrix to coordinate format
GLPKFormat to_GLPK_format(SparseMatrix<double,RowMajor>& other, bool transpose=false) {
	GLPKFormat ret;

	if(transpose==true) {
		ret.non_zeros = other.nonZeros();
		ret.ncols = other.cols();
		ret.nrows = other.rows();
	}
	else {
		ret.non_zeros = other.nonZeros();
		ret.nrows = other.rows();
		ret.ncols = other.cols();
	}

	// GLPK indexing starts at 1, we insert a 0 as 0th element.
	ret.rows = vector<int> (ret.non_zeros+1,0);
	ret.cols = vector<int> (ret.non_zeros+1,0);
	ret.values = vector<double> (ret.non_zeros+1,0.);

	if(transpose==true) {
		// copy data
		for(unsigned k=0; k<other.outerSize(); k++) {
			for (typename SparseMatrix<double,RowMajor>::InnerIterator it(other,k); it; ++it) {
				// GLPK indexing starts at 1, hence row and column indices have to be shifted
				ret.rows.at(k+1) = it.col()+1;
				ret.cols.at(k+1) = it.row()+1;
				ret.values.at(k+1) = it.value();
			}  				
		}	
	}
	else {
		// copy data
		for(unsigned k=0; k<other.outerSize(); k++) {
			for (typename SparseMatrix<double,RowMajor>::InnerIterator it(other,k); it; ++it) {
				// GLPK indexing starts at 1, hence row and column indices have to be shifted
				ret.rows.at(k+1) = it.row()+1;
				ret.cols.at(k+1) = it.col()+1;
				ret.values.at(k+1) = it.value();
			}  				
		}
	}

	return ret;
}

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

// // converts a Eigen::SparseMatrix to coordinate format
// GLPKFormat to_GLPK_format(SparseMatrix<double,RowMajor>& other) {
// 	GLPKFormat ret;

// 	ret.non_zeros = other.nonZeros();
// 	ret.nrows = other.rows();
// 	ret.ncols = other.cols();

// 	// GLPK indexing starts at 1, we insert a 0 as 0th element.
// 	ret.rows = vector<int> (ret.non_zeros+1,0);
// 	ret.cols = vector<int> (ret.non_zeros+1,0);
// 	ret.values = vector<double> (ret.non_zeros+1,0.);


// 	// copy data
// 	for(unsigned k=0; k<other.outerSize(); k++) {
// 		for (typename SparseMatrix<double,RowMajor>::InnerIterator it(other,k); it; ++it) {
// 			// GLPK indexing starts at 1, hence row and column indices have to be shifted
// 			ret.rows.at(k+1) = it.row()+1;
// 			ret.cols.at(k+1) = it.col()+1;
// 			ret.values.at(k+1) = it.value();
// 		}  				
// 	}

// 	return ret;
// }

// // reads sparse matrix in COO format and saves to struct
// GLPKFormat to_GLPK_format(string cmatrix_file, bool transpose=false) {
// 	GLPKFormat ret;

// 	// number of non-zero values
// 	ret.non_zeros = get_number_of_lines(cmatrix_file);

// 	// GLPK convention: data start at index 1
// 	ret.rows = vector<int> (ret.non_zeros+1,0);
// 	ret.cols = vector<int> (ret.non_zeros+1,0);
// 	ret.values = vector<double> (ret.non_zeros+1,0.);

	
// 	// open file and read
// 	fstream fin(cmatrix_file, ios::in);
// 	int n = 1;
// 	if(fin.is_open()) {
// 		while(fin >> ret.rows[n] >> ret.cols[n] >> ret.values[n]) {
// 			++n;
// 		}
// 	} else {
// 		cerr << "Error: Couldn't open file " << cmatrix_file << endl;
// 	}
// 	fin.close();

// 	// save number of rows and columns
// 	ret.nrows = ret.rows.at(ret.non_zeros);
// 	ret.ncols = *max_element(ret.cols.begin(),ret.cols.end());


// 	return ret;
// }

void print(GLPKFormat data) {
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