#include <iostream>
#include <vector>
#include <random>

#include "stabiliser.h"
#include "utilities.h"


using namespace std;

// ----- Test of functions in utilities.h

int main(int argc, char** argv) {

	unsigned N = 6;
	unsigned samples = 10;


	// ----- test of bit manipulation functions


	// cout << "=================================================" << endl;
	// cout << "===   Start testing of bin manipulation" << endl;
	// cout << "=================================================" << endl;
	
	// random_device rd; 
	// mt19937 gen(rd()); 
	// uniform_int_distribution<> dist (1, pow(2,N)-1);
	// unsigned x,y;
	
	// for(unsigned i=0; i<samples; i++) {
	// 	x = dist(gen);
	// 	y = invert_bits(x,N);

	// 	cout << write_bits(x,N) << " --> " << write_bits(y,N) << endl;	
	// }

	// cout << "------------------------" << endl;
	// for(unsigned i=0; i<samples; i++) {
	// 	x = dist(gen);
	// 	y = reset_bit(x,2,N);

	// 	cout << write_bits(x,N) << " --> " << write_bits(y,N) << endl;	
	// }


	// cout << "=================================================" << endl;
	// cout << "===   Start testing of index sets" << endl;
	// cout << "=================================================" << endl;

	// random_device rd; 
	// mt19937 gen(rd()); 
	// uniform_int_distribution<> dist (5, 15);

	// unsigned n = dist(gen);

	// cout << "Testing all index sets for n = " << n << " qubits for correctness ..." << endl;

	// for(unsigned k=0; k<=n; k++) {
	// 	for(unsigned kk=0; kk<=k; kk++) {
	// 		vector<unsigned> iset = index_set2(kk,k,n);
	// 		for(auto set : iset) {
	// 			vector<unsigned> weights = count_weights(set,n);
	// 			if(weights.at(1)+weights.at(3) != k || weights.at(3) != kk) {
	// 				cout << "Something went wrong with the following index set entry with (kk,k,n)=(" << kk << "," << k << "," << n << endl;
	// 				cout << write_bits(set,2*n) << endl;
	// 				break;
	// 			}
	// 		}
	// 	}
	// }

	cout << "=================================================" << endl;
	cout << "===   Start testing of partitions" << endl;
	cout << "=================================================" << endl;

	N = 5;
	cout << "All partitions for n = " << N << endl;

	vector<vector<unsigned>> parts = get_partitions_w_permutations(N,2);
	for(auto p : parts) {
		for(auto i : p) {
			cout << i << " ";
		}
		cout << endl;
	}
}