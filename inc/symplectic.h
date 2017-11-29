#ifndef SYMPLECTIC_H
#define SYMPLECTIC_H
#include <initializer_list>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>

#ifndef UTILITIES_H
#include "utilities.h"
#endif

using namespace std;



// ---- Phase space methods

unsigned eta(const unsigned a, const unsigned b, const unsigned n) {
	unsigned ret = 0;

	for(unsigned i=0; i<n; i++) {
		ret ^= ( get_bit(a,2U*i+1U,2U*n) * get_bit(b,2U*i,2U*n) );
		// ret ^= ( (a >> (2U*i+1U) ) & 1U ) & ( (b >> 2U*i ) & 1U );
	}

	return ret;
}

unsigned symplectic_form(const unsigned a, const unsigned b, const unsigned n) {
	return ( eta(a,b,n) ^ eta(b,a,n) );
}

// This is the convolution p3 = p1 * p2 in (Z_2)^N (addition mod 2)
int convolve_mod2(double *p1, double *p2, double *p3, const unsigned N) {
  unsigned len = pow(2U,N);
  unsigned i,j;

  for(i=0; i<len; i++) {
    p3[i] = 0;
    
    for(j=0; j<len; j++) {
      p3[i] += p1[j] * p2[i^j];
    } 
  }
  return 0;
}

// Computes the matrix vector product Ax over Z_2. It is assumed that a vector in Z_2^N is represented as a bitstring encoded as an unsigned integer. In this representation, every row of the matrix A corresponds to an unsigned integer, hence A is an unsigned array. The ouput is again an unsigned integer.
unsigned matrix_vector_prod_mod2(const unsigned *A, const unsigned x, const unsigned N) {
	unsigned b = 0;

	/* This performs the formula:
		b_i = sum_{j=1}^N A_ij x_j
		
		In the bitstring representation, this means that we have to multiply A[i] with x bitwise, afterwards adding up the resulting bits mod 2. The last step is just counting the number of ones and taking its modulus (here, called parity). Finally, we set the i-th bit of b to this result.
		Note that, in contrast to the usual convention, we count bits from the left since we see bitstrings as binary vectors. This means that the last bit is the least significant one.
	*/
	for (unsigned i=0; i<N; i++) 
		b ^= ((-parity( A[i]&x )) ^ b) & (1U << (N-i-1U)); 

	return b;	
}

void symplectic_transform(double *pin, double *pout, const unsigned *S, const unsigned n) {
	unsigned j;
	for(unsigned i=0; i<pow(4,n); i++) {
		j = matrix_vector_prod_mod2(S, i, 2U*n);
		pout[j] = pin[i];
	}
}

unsigned transvection(const unsigned h, const unsigned x, const unsigned n) {
	// returns the transvection of x by h, i.e.
	//    y = x + [x,h]h
	return ( x ^ (symplectic_form(x,h,n) * h) );
}

// ---- Symplectic group generation

unsigned phase_function(const unsigned a, const unsigned *S, const unsigned n) {
	// represents the phase function g of the trivial Clifford representative C with [C]=S in Sp(2n) = C(n)/P(n). It describes the action of C as follows:
	//		C W_a C^\dagger = (-1)^{g(a)} W_{Sa} 

	if(popcount(a) <= 1) {
		// hence, either a == 0 or a is a basis vector in phase space
		return 0;
	}

	// now we split a in two parts, a = aa + lsb, where lsb only contains the least significant bit of a and hence represents a basis vector. Then we use the summation rule for g on this decomposition

	// get least significant bit
	unsigned lsb = a & ~(a-1U);
	// reset it in a 
	unsigned aa = a & (a-1U);

	// write_bits(lsb,2*n);
	// write_bits(aa,2*n);

	// recursion. This will terminate if aa has only one set bit left
	return (( phase_function(aa,S,n) + eta(aa,lsb,n) + eta( matrix_vector_prod_mod2(S,aa,2*n), matrix_vector_prod_mod2(S,lsb,2*n), n ) ) % 2);
}

vector<unsigned> find_transvection(const unsigned x, const unsigned y, const unsigned n) {
	// Returns binary vectors h1, h2 (represented as a vector of 2 2n-bit unsigned integers)
	// such that 
	//    y = Z(h1)Z(h2)x
	// where Z(h) is the symplectic transvection given by
	//    Z(h)x = x + [x,h]h 		([,] being the symplectic product)
	// Note: for some pairs (x,y) only one transvection is needed and h1=0 is outputted
	//
	// Debug note: returns a third element in the vector which indicates at which part
	// of this function returned the vector

	if(x == 0 || y == 0) {
		// this doesn't work for zero vectors
		return vector<unsigned>(3,0U);
	}
	
	if(x == y) {
		// return 0
		return vector<unsigned>(3,0U);
	} 

	if(symplectic_form(x,y,n) == 1) {
		// return h = x + y (bitwise addition)
		return vector<unsigned>({0, x^y,1});
	}
	else {
		unsigned z = 0; 
		unsigned xx,yy;
		for(unsigned i=0; i<n; i++) {
			// extract bits 2i and 2i+1
			xx = get_bits(x,2U*i,2U*i+1U,2*n);
			yy = get_bits(y,2U*i,2U*i+1U,2*n);

			if( (xx != 0) && (yy != 0) ) {
				// found a pair, set z to bitwise addition of the pair
				z = xx ^ yy;
				if(z == 0) {	// i.e. xx == yy
					// set first bit of z
					z = set_bit(z,1,2); 

					if( xx != 3 ) {		// hence the bits are not equal 
						z = set_bit(z,0,2);
					}
				}
				// finish construction by shifting the bits in z to position 2i and 2i+1
				z = z << (2*n-2U*i-2U);
				return vector<unsigned>({x^z,y^z,2});
			}
		}

		// we didn't find a pair

		// look for a pair 00 in y
		for(unsigned i=0; i<n; i++) {
			// extract bits 2i and 2i+1
			xx = get_bits(x,2U*i,2U*i+1U,2*n);
			yy = get_bits(y,2U*i,2U*i+1U,2*n);

			if( (xx != 0) && (yy == 0) ) {
				// found it

				// xx is either 01, 10 or 11
				if(xx == 3) {
					// xx == 11
					z = set_bit(z,1,2);
				}
				else if(xx == 1) {			
					// set z to negation of that
					z = 2;
				}
				else {
					z = 1;
				}
				// shift bits in z again
				z = z << (2*n-2U*i-2U);
				break;
			}
		}

		// we need some temporary z
		unsigned zz = 0;

		// look for a pair 00 in x 
		for(unsigned i=0; i<n; i++) {
			// extract bits 2i and 2i+1
			xx = get_bits(x,2U*i,2U*i+1U,2*n);
			yy = get_bits(y,2U*i,2U*i+1U,2*n);

			if( (xx == 0) && (yy != 0) ) {
				// found it

				// yy is either 01, 10 or 11
				if(yy == 3) {
					// yy == 11
					zz = set_bit(zz,1,2);
				}
				else if(yy == 1) {			
					// set z to negation of that
					zz = 2;
				}
				else {
					zz = 1;
				}
				// shift bits in zz again
				zz = zz << (2*n-2U*i-2U);
				break;
			}
		}

		// finally we "add" z and zz (set bits shall not cancel, hence OR instead of XOR)
		z = z|zz;
		return vector<unsigned>({x^z,y^z,3});
	}

}

#endif