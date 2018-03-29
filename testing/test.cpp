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

	unsigned n = 5;
	binvec A = graph6_to_adj_mat("DCo");

	// reserve memory
	vector<binvec> B(n,0);
	unsigned N = n*(n-1)/2;

	// fill B 
	unsigned i,j;

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
}