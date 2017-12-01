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

// ---------------------------------------------------
// ------ Header file for phase space functions ------
// ---------------------------------------------------
//
// Throughout this code, binary vectors in the n-qubit phase space Z_2^{2n} are represented as unsigned integers using their binary representation, e.g.
// 		(1,0,0,1) \in Z_2^4 ----> 1001 = 9
// We use the coordinate convention (z_1,x_1,z_2,x_2,...,z_n,x_n) in phase space such that Z_2^{2n} decomposes into a direct sum of 1-qubit phase spaces.
// Matrices are represented as a vector of unsigned integers, where every entry of this vector corresponds to a row of the matrix.
//
//
// Important note:
// ---------------
// Right now we use the datatype unsigned. Its size is 16 bit, meaning that in can at most represent a binary vector of length 2*8. Hence, this code is limited right now to at most 8 qubits! 
// In one of the next versions, we define a custom datatype which can be set as needed, e.g. to unsigned long
//


// ---- For working on finite fields

// actually returns the modulus instead of remainder
// e.g mod(-2,3) == 1 instead of (-2)%3 = -2
int mod(int x, int m) {
    return (x%m + m)%m;
}

// ---- Phase space methods

unsigned eta(const unsigned a, const unsigned b, const unsigned n) {
	// Computes the phase that is appearing in the composition law of Weyl operators:
	//		W(a)W(b) = i^{\eta(a,b)+\eta(b,a)} W(a+b)
	// \eta is actually a Z_4-valued function.
	unsigned ret = 0;

	for(unsigned i=0; i<n; i++) {
		ret += ( get_bit(a,2U*i+1U,2U*n) * get_bit(b,2U*i,2U*n) );
	}

	return ret;
}

unsigned eta(const vector<unsigned> a, const unsigned n) {
	// Compute the phase that is appearing in the product of many Weyl operators
	//		W(a_1)W(a_2)...W(a_m) = i^{-\eta(a_1,...,a_m)}W(a_1+...+a_m)
	// This phase is given by
	//		\eta(a_1,...,a_m) = \sum_{j=1}^{m-1} \eta(a_j,\sum_{k=j+1}^m a_k)
	//
	unsigned ret = 0;
	unsigned m = a.size();
	unsigned aa;
	for(unsigned j=0; j<m-1; j++) {
		aa = 0;
		// sum up vectors 
		for(unsigned k=j+1; k<m; k++) {
			aa ^= a.at(k);
		}
		// add eta function
		ret += eta(a.at(j),aa,n);
	}
	return ret;
}

int phi(const unsigned a, const unsigned b, const unsigned n) {
	// Computes the phase that is appearing in the composition law of Weyl operators:
	//		W(a)W(b) = i^{\phi(a,b)} W(a+b)
	// \phi is actually a Z_4-valued function.
	return ( 2*eta(a,b,n) + eta(a^b,a^b,n) - eta(a,a,n) - eta(b,b,n) );
}

int phi(const vector<unsigned> a, const unsigned n) {
	// Compute the phase that is appearing in the product of many Weyl operators
	//		W(a_1)W(a_2)...W(a_m) = i^{\phi(a_1,...,a_m)}W(a_1+...+a_m)
	// This phase is given by
	//		\phi(a_1,...,a_m) = \sum_{j=1}^{m-1}\sum_{k=j+1}^m \phi(a_j, a_k)
	//
	unsigned ret = 0;
	unsigned m = a.size();

	// sum up phases
	for(unsigned j=0; j<m-1; j++) {
		for(unsigned k=j+1; k<m; k++) {
			ret += phi(a.at(j),a.at(k),n);
		}
	}
	return ret;
}

unsigned symplectic_form(const unsigned a, const unsigned b, const unsigned n) {
	return ( eta(a,b,n) + eta(b,a,n) )%2;
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

unsigned matrix_vector_prod_mod2(const vector<unsigned> A, const unsigned x) {
	unsigned b = 0;

	/* This performs the formula:
		b_i = sum_{j=1}^N A_ij x_j
		
		In the bitstring representation, this means that we have to multiply A[i] with x bitwise, afterwards adding up the resulting bits mod 2. The last step is just counting the number of ones and taking its modulus (here, called parity). Finally, we set the i-th bit of b to this result.
		Note that, in contrast to the usual convention, we count bits from the left since we see bitstrings as binary vectors. This means that the last bit is the least significant one.
	*/
	unsigned N = A.size();

	for (unsigned i=0; i<N; i++) 
		b ^= ((-parity( A.at(i)&x )) ^ b) & (1U << (N-i-1U)); 

	return b;	
}

vector<unsigned> direct_sum(const vector<unsigned> A1, const vector<unsigned> A2) {
	// Returns the direct sum A1 + A2 of the square matrices A1 and A2

	unsigned n1 = A1.size();
	unsigned n2 = A2.size();

	// construct new vector
	vector<unsigned> R (n1+n2);

	for(unsigned j=0; j<n1; j++) {
		R.at(j) = A1.at(j) << n2; // fills in n2 zeros to the right
	}
	for(unsigned j=n1; j<n1+n2; j++) {
		R.at(j) = A2.at(j-n1);
	}

	return R;

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

unsigned phase_function(const unsigned a, const vector<unsigned> S) {
	// represents the phase function g of the trivial Clifford representative C with [C]=S in Sp(2n) = C(n)/P(n). It describes the action of C as follows:
	//		C W_a C^\dagger = (-1)^{g(a)} W_{Sa} 

	if(popcount(a) <= 1) {
		// hence, either a == 0 or a is a basis vector in phase space
		return 0;
	}

	unsigned n = S.size()/2;

	// now we split a in two parts, a = aa + lsb, where lsb only contains the least significant bit of a and hence represents a basis vector. Then we use the summation rule for g on this decomposition

	// get least significant bit
	unsigned lsb = a & ~(a-1U);
	// reset it in a 
	unsigned aa = a & (a-1U);

	// recursion. This will terminate if aa has only one set bit left
	return (( phase_function(aa,S) + eta(aa,lsb,n) + eta( matrix_vector_prod_mod2(S,aa), matrix_vector_prod_mod2(S,lsb), n ) ) % 2);
}

int phase_function2(const unsigned a, const vector<unsigned> S) {
	// represents the phase function g of the trivial Clifford representative C with [C]=S in Sp(2n) = C(n)/P(n). It describes the action of C as follows:
	//		C W_a C^\dagger = (-1)^{g(a)} W_{Sa} 

	if(popcount(a) <= 1) {
		// hence, either a == 0 or a is a basis vector in phase space
		return 0;
	}

	unsigned n = S.size()/2;

	// will be used for a basis expansion of a
	// we will actually only store set bits in a
	vector<unsigned> av;
	av.reserve(2*n);
	vector<unsigned> bv;
	bv.reserve(2*n);

	unsigned aa = a;
	while(aa != 0) {
		// push back least significant bit
		av.push_back(aa & ~(aa-1));
		bv.push_back(matrix_vector_prod_mod2(S,aa & ~(aa-1)));
		// reset it in aa
		aa = aa & (aa-1);
	}
	// compute phis
	return (phi(bv,n) - phi(av,n));

}

unsigned long sp_order(const unsigned n) {
	unsigned long r = 1;
	for(unsigned j=1; j<=n; j++) {
		r *= pow(4,j)-1;
	}
	r *= pow(2,n*n);
	return r;
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

vector<unsigned> generate_symplectic_matrix(const unsigned i, const unsigned n) {
	// Generate the i-th element in Sp(2n, Z_2) using transvections
	// Actually, returns the transposed element.

	unsigned nn = 2*n;
	unsigned s = (1<<nn)-1;

	// first and second basis vectors e1=(1,0,0,...,0), e2=(0,1,0,...,0)
	unsigned e1 = set_bit(0,0,nn);
	unsigned e2 = set_bit(0,1,nn);

	// compute phase space points
	unsigned f1 = invert_bits((i%s)+1,nn);
	unsigned b = invert_bits(i/s,nn-1);
	unsigned eprime = e1 ^ b;
	eprime = reset_bit(eprime,1,nn); // result is (1,0,b1,b2,b3,...,b2n-1)

	// compute and apply transvections
	vector<unsigned> T = find_transvection(e1,f1,n);
	
	unsigned h0 = transvection(T.at(1),eprime,n);
	h0 = transvection(T.at(0),h0,n);

	// check 0th bit of b
	if(get_bit(b,0,nn-1) == 1) {
		f1 = 0;
	}

	// unsigned f2 = transvection(T.at(1),e2,n);
	// f2 = transvection(T.at(0),f2,n);
	// f2 = transvection(Tprime.at(1),f2,n);
	// f2 = transvection(Tprime.at(0),f2,n);

	vector<unsigned> g;
	vector<unsigned> id = {2,1}; // 2x2 identity matrix

	if(n == 1) {
		g = id;
	}
	else {
		// recursive call
		g = direct_sum(id, generate_symplectic_matrix((i/s) >> (nn-1),n-1));
	}

	// finally, apply the transvections to each row in g
	for(unsigned j=0; j<nn; j++) {
		g[j] = transvection(T.at(1),g[j],n);
		g[j] = transvection(T.at(0),g[j],n);
		g[j] = transvection(h0,g[j],n);
		g[j] = transvection(f1,g[j],n);
	}

	return g;
}

vector<int> projected_liouville_matrix(const vector<unsigned> S, const unsigned a=0) {
	// Computes the Liouville matrix of the "trivial" Clifford representative corresponding to the symplectic matrix S, i.e. the one that acts Z/X Paulis as
	//		C:  W(e_j) ---> +W(Se_j)
	// Its action on an arbitrary Pauli operator W(b) is
 	//		C:  W(b) ---> (-1)^g(b) W(Sb)
 	// where g(b) is some phase function determined by S.
 	//
	// Optionally, if 0<=a<=3 is specified, we include a term (-1)^[c,b_1] in the phase, which corresponds to the action of the 1-qubit Pauli group P_1. b_1 are the first two entries in the binary vector b.

	unsigned n = S.size()/2;
	unsigned m = pow(4,n+1);
	unsigned b = 0;
	vector<int> L(m-1);

	for(unsigned i=1; i<m; i++) {
		// get the last 2n bits of i by setting the 0th and 1st to 0
		b = reset_bit(i,0,2*n+2);
		b = reset_bit(b,1,2*n+2);

		// check if the first two bits of Sb agree with the first two bits of i
		if(get_bits(i,0,1,2*n+2) == get_bits(matrix_vector_prod_mod2(S, b),0,1,2*n)) {
			// they do
			L.at(i-1) = pow(-1, phase_function(b,S) + symplectic_form(a,get_bits(b,0,1,2*n),1) );
		}
		else{
			L.at(i-1) = 0;
		}
	}
	return L;
}

vector<int> liouville_matrix(const vector<unsigned> S, const unsigned c=0) {
	// Computes the Liouville matrix of the "trivial" Clifford representative corresponding to the symplectic matrix S, i.e. the one that acts Z/X Paulis as
	//		C:  W(e_j) ---> +W(Se_j)
	// Its action on an arbitrary Pauli operator W(b) is
 	//		C:  W(b) ---> (-1)^g(b) W(Sb)
 	// where g(b) is some phase function determined by S.
 	//
	// Optionally, if 0<=c<2*n is specified, we include a term (-1)^[c,b] in the phase, which corresponds to the action of the n-qubit Pauli group P_n. 

	unsigned n = S.size()/2;
	unsigned m = 1 << 4*n;	// pow(16,n)
	unsigned a,b;
	vector<int> L(m);

	for(unsigned i=0; i<m; i++) {
		// get the first and last 2n bits of i
		a = get_bits(i,0,2*n-1,4*n);
		b = get_bits(i,2*n,4*n-1,4*n);

		// check if a = Sb
		if( a == matrix_vector_prod_mod2(S, b) ) {
			// they do
			L.at(i) = pow(-1, phase_function2(b,S) + symplectic_form(c,b,n) );
		}
		else{
			L.at(i) = 0;
		}
	}
	return L;
}

vector<int> liouville_matrix2(const vector<unsigned> S, const unsigned c=0) {
	// Computes the Liouville matrix of the "trivial" Clifford representative corresponding to the symplectic matrix S, i.e. the one that acts Z/X Paulis as
	//		C:  W(e_j) ---> +W(Se_j)
	// Its action on an arbitrary Pauli operator W(b) is
 	//		C:  W(b) ---> (-1)^g(b) W(Sb)
 	// where g(b) is some phase function determined by S.
 	//
	// Optionally, if 0<=c<2*n is specified, we include a term (-1)^[c,b] in the phase, which corresponds to the action of the n-qubit Pauli group P_n. 

	unsigned n = S.size()/2;
	unsigned m = 1 << 4*n;	// pow(16,n)
	unsigned a,b;
	vector<int> L(m);
	int ph;

	for(unsigned i=0; i<m; i++) {
		// get the first and last 2n bits of i
		a = get_bits(i,0,2*n-1,4*n);
		b = get_bits(i,2*n,4*n-1,4*n);

		// check if a = Sb
		if( a == matrix_vector_prod_mod2(S, b) ) {
			// they do
			L.at(i) = pow(-1, symplectic_form(c,b,n) );

			int ph = mod(phase_function2(b,S),4);
			assert(ph == 0 || ph == 2);
			if(ph == 2) {
				L.at(i) *= -1;
			}
			
			
		}
		else{
			L.at(i) = 0;
		}
	}
	return L;
}

#endif