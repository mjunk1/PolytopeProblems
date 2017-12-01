#include <iostream>
#include <vector>
#include <random>

// #include "symplectic.h"
#include "utilities.h"


using namespace std;

// ----- Test of functions in utilities.h

int main(int argc, char** argv) {

	unsigned N = 6;
	unsigned samples = 10;


	// ----- test of bit manipulation functions


	cout << "=================================================" << endl;
	cout << "===   Start testing of bin manipulation" << endl;
	cout << "=================================================" << endl;
	
	random_device rd; 
	mt19937 gen(rd()); 
	uniform_int_distribution<> dist (1, pow(2,N)-1);
	unsigned x,y;
	
	for(unsigned i=0; i<samples; i++) {
		x = dist(gen);
		y = invert_bits(x,N);

		cout << write_bits(x,N) << " --> " << write_bits(y,N) << endl;	
	}

	cout << "------------------------" << endl;
	for(unsigned i=0; i<samples; i++) {
		x = dist(gen);
		y = reset_bit(x,2,N);

		cout << write_bits(x,N) << " --> " << write_bits(y,N) << endl;	
	}
}