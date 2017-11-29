#include <iostream>
#include <vector>
#include <random>

#include "symplectic.h"
#include "utilities.h"


using namespace std;

// ----- Test of functions in symplectic.h

int main(int argc, char** argv) {

	unsigned n = 20;
	unsigned samples = 200;

	// ----- test of symplectic product

	// TBD


	// ----- test of transvections

	// list of transvection triples
	// TBD

	// random generation of pairs of points

	cout << "=================================================" << endl;
	cout << "===   Starting testing of transvections" << endl;
	cout << "=================================================" << endl;
	random_device rd; 
	mt19937 gen(rd()); 
	uniform_int_distribution<> dist (1, pow(4,n)-1);
	unsigned x,y,z;
	vector<unsigned> h;

	for(unsigned i=0; i<samples; i++) {
		x = dist(gen);
		y = dist(gen);
		h = find_transvection(x,y,n);

		z = transvection(h.at(1),x,n);
		z = transvection(h.at(0),z,n);


		cout << "  ** Transvection for random pair (x,y)=(" << write_bits(x,2*n) << "," << write_bits(y,2*n) << "):     h=(" << write_bits(h.at(0),2*n) << "," << write_bits(h.at(1),2*n) << ")  ...  ";

		if(z == y) {
			cout << "correct" << endl;
		}
		else {
			cout << "failed with status " << h.at(2) << endl;
			cerr << "Error: Transvection test not passed with random pair (x,y)=(" <<  write_bits(x,2*n) << "," << write_bits(y,2*n) << ")" << endl;
			cout << "Debug output: z=" << write_bits(h.at(3),2*n) << ", zz=" << write_bits(h.at(4),2*n) << endl;
			break;
		}
	}

}