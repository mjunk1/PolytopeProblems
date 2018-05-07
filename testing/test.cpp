#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

// #include "symplectic.h"
#include "utilities.h"
#include "stabiliser.h"


using namespace std;

// ----- Test of functions 

int main(int argc, char** argv) {

	string file = "data/graphs/graph5.g6";
	fstream fin(file, ios::in);

	vector<string> graphs;
	string row;

	if(fin.is_open()) {
		while(fin >> row) {
			graphs.push_back(row);
		}
	} 
	else {
		cout << "Couldn't open file " << file << endl;
	}

	fin.close();

	// get number of qubits
	unsigned n = get_order( graphs.at(0) );
	binvec A;
	vector<binvec> B(n,0);
	unsigned N = n*(n-1)/2;
	unsigned i,j;
	unsigned counter = 1;

	for(auto g : graphs) {
		A = graph6_to_adj_mat(g);
		B.assign(n,0);

		cout << "Graph # " << counter << endl;
		// print A
		for(unsigned k=0; k<N; k++) {
			if( get_bit(A, k, N) == 1) {
				i = get_ut_row(k,n);
				j = get_ut_col(i,k,n);

				B[i] = set_bit(B[i], j, n);
				B[j] = set_bit(B[j], i, n);
			}
		}

		for(unsigned k=0; k<n; k++) {
			cout << write_bits(B[k],n) << endl;
		}

		cout << "Size of zero block: " << find_zero_block_size(A, n) << endl;

		cout << endl;

		counter++;
	}
}