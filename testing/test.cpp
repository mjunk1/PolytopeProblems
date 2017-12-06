#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

// #include "symplectic.h"
#include "utilities.h"


using namespace std;

// ----- Test of functions in utilities.h

int main(int argc, char** argv) {

	unsigned n = 2;
	unsigned b;

	for(unsigned i=30; i<55; i++) {
		cout << write_bits(i,2*n+2) << endl;

		// get the last 2n bits of i by setting the 0th and 1st to 0
		b = reset_bit(i,0,2*n+2);
		b = reset_bit(b,1,2*n+2);

		cout << write_bits(get_bits(i,2,2*n+1,2*n+2),2*n) << endl;
	}
}