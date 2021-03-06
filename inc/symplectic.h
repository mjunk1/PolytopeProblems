#ifndef SYMPLECTIC_H
#define SYMPLECTIC_H
#include <initializer_list>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <set>
#include <cstdint>

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


// ---- typedefs

// type that is used to represent binary vectors in phase space
// note that for n qubits, the binary vector corresponds to a sequence of 2n bits, so we use by default uint64_t which can represent 64 bits, so works definitely for up to n=32
typedef uint64_t binvec;
typedef vector<binvec> symplectic_matrix; 


// ----- bit manipulation for operating on binary vectors

unsigned popcount(unsigned i) {
	return __builtin_popcount(i);
}

unsigned parity(unsigned i) {
	return __builtin_parity(i);
}


// get j-th bit of b, counted from the end of the bitstring, i.e. from the left
inline binvec get_bit(binvec b, unsigned j, unsigned n) {
	return ((b >> ((n)-(j)-1U)) & 1U);
}

// get j-th bit of b, counted from the from the right
inline binvec get_bit2(binvec b, unsigned j) {
	return ((b >> j) & 1U);
}

inline binvec get_bits(binvec b, unsigned j, unsigned k, unsigned n) {
	return ((b >> (n-k-1U)) & ((1U << (k-j+1U))-1));
}


inline binvec set_bit(binvec b, unsigned j, unsigned n) {
	return (b | (1 << (n-j-1U)));
}

inline binvec reset_bit(binvec b, unsigned j, unsigned n) {
	return (b & (~(1 << (n-j-1U))));
}

string write_bits(binvec b, unsigned n) {
	string str;

	for(unsigned j=0; j<n; j++) {
		str += to_string(get_bit(b,j,n)); 
	}

	return str;
}

string write_trits(binvec t, unsigned n) {
	string str;

	vector<unsigned> trits(n,0);
	get_multi_index(n, 3, t, trits);

	for(unsigned j=0; j<n; j++) {
		str += to_string( trits.at(j) ); 
	}

	return str;
}

binvec invert_bits(const binvec b, const unsigned n) {
	// Returns a copy of b with inverted ordering of bits, e.g.
	// 		0010110 ---> 0110100
	// Note that this function "sees" only the first n bits, all the other bits are set to zero. I.e. calling the function for b = 0010110 with n = 6 will instead give
	// 		0010110 ---> 0011010

	unsigned r = 0;
	for(unsigned i=0; i<n; i++) {
		if(get_bit(b,i,n) == 1) {
			r |= 1<<i;
		}
	}
	return r;
}

// ---- For working on finite fields

// actually returns the modulus instead of remainder
// e.g mod(-2,3) == 1 instead of (-2)%3 = -2
int mod(int x, int m) {
    return (x%m + m)%m;
}

// ---- Phase space methods

unsigned eta(const binvec a, const binvec b, const unsigned n) {
	// Computes the phase that is appearing in definition of Weyl operators:
	//		W(a) = i^{-\eta(a,a)} Z(a_z)X(a_x)
	// \eta is actually a Z_4-valued function.
	unsigned ret = 0;

	for(unsigned i=0; i<n; i++) {
		ret += ( get_bit(a,2U*i+1U,2U*n) * get_bit(b,2U*i,2U*n) );
	}
	return ret;
}

int phi(const binvec a, const binvec b, const unsigned n) {
	// Computes the phase that is appearing in the composition law of Weyl operators:
	//		W(a)W(b) = i^{\phi(a,b)} W(a+b)
	// \phi is actually a Z_4-valued function.
	return ( 2*eta(a,b,n) + eta(a^b,a^b,n) - eta(a,a,n) - eta(b,b,n) );
}

int phi(const vector<binvec> &a, const unsigned n) {
	// Compute the phase that is appearing in the product of many Weyl operators
	//		W(a_1)W(a_2)...W(a_m) = i^{\phi(a_1,...,a_m)}W(a_1+...+a_m)
	// This phase is given by
	//		\phi(a_1,...,a_m) = \sum_{j=1}^{m-1} \phi(a_j, \sum_{k=j+1}^m a_k)
	//
	unsigned ret = 0;
	unsigned m = a.size();
	binvec aa;

	// sum up phases
	for(unsigned j=0; j<m-1; j++) {
		aa = 0;
		for(unsigned k=j+1; k<m; k++) {
			// sum up the vectors
			aa ^= a.at(k);
		}
		ret += phi(a.at(j),aa,n);
	}
	return ret;
}

unsigned symplectic_form(const binvec a, const binvec b, const unsigned n) {
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

binvec matrix_vector_prod_mod2(const binvec *A, const binvec x, const unsigned N) {
	binvec b = 0;

	/* This performs the formula:
		b_i = sum_{j=1}^N A_ij x_j
		
		In the bitstring representation, this means that we have to multiply A[i] with x bitwise, afterwards adding up the resulting bits mod 2. The last step is just counting the number of ones and taking its modulus (here, called parity). Finally, we set the i-th bit of b to this result.
		Note that, in contrast to the usual convention, we count bits from the left since we see bitstrings as binary vectors. This means that the last bit is the least significant one.
	*/
	for (unsigned i=0; i<N; i++) 
		b ^= ((-parity( A[i]&x )) ^ b) & (1U << (N-i-1U)); 

	return b;	
}

binvec matrix_vector_prod_mod2(const vector<binvec> &A, const binvec x) {
	binvec b = 0;

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

vector<binvec> direct_sum(const vector<binvec> &A1, const vector<binvec> &A2) {
	// Returns the direct sum A1 + A2 of the square matrices A1 and A2

	unsigned n1 = A1.size();
	unsigned n2 = A2.size();

	// construct new vector
	vector<binvec> R (n1+n2);

	for(unsigned j=0; j<n1; j++) {
		R.at(j) = A1.at(j) << n2; // fills in n2 zeros to the right
	}
	for(unsigned j=n1; j<n1+n2; j++) {
		R.at(j) = A2.at(j-n1);
	}

	return R;

}
vector<binvec> direct_sum(const vector<vector<binvec>> &Alist) {
	vector<binvec> temp = Alist.at(0);
	for(unsigned i=1; i<Alist.size(); i++) {
		temp = direct_sum(temp, Alist.at(i));
	}

	return temp;
}

vector<binvec> ut_to_dense_matrix(const binvec A, unsigned n) {
	vector<binvec> R (n,0);
	unsigned N = n*(n-1)/2;
	unsigned i,j;

	for(unsigned k=0; k<N; k++) {
		if(get_bit(A,k,N) == 1) {
			i = get_ut_row(k,n);
			j = get_ut_col(i,k,n);
			R.at(i) = set_bit(R.at(i), j, n);
			R.at(j) = set_bit(R.at(j), i, n);
		}
	}

	return R;
}

binvec dense_to_ut_matrix(vector<binvec> &A) {
	unsigned n = A.size();
	unsigned N = n*(n-1)/2;
	binvec R = 0;
	unsigned i,j;
	
	for(unsigned k=0; k<N; k++) {	
		i = get_ut_row(k,n);
		j = get_ut_col(i,k,n);
		if(get_bit(A.at(i),j,n) == 1) {
			R = set_bit(R, k, N);
		}
	}	

	return R;
}

// direct sum of two upper triangle matrices
binvec direct_sum_ut(const binvec A1, const unsigned n1, const binvec A2, const unsigned n2) {
	// Returns the direct sum A1 + A2 of the upper triangle matrices A1 and A2

	// construct dense matrices
	vector<binvec> a1 = ut_to_dense_matrix(A1, n1);
	vector<binvec> a2 = ut_to_dense_matrix(A2, n2);
	vector<binvec> sum = direct_sum(a1, a2);

	return dense_to_ut_matrix(sum);
}

void symplectic_transform(double *pin, double *pout, const binvec *S, const unsigned n) {
	unsigned j;
	for(unsigned i=0; i<pow(4,n); i++) {
		j = matrix_vector_prod_mod2(S, i, 2U*n);
		pout[j] = pin[i];
	}
}

binvec transvection(const binvec h, const binvec x, const binvec n) {
	// returns the transvection of x by h, i.e.
	//    y = x + [x,h]h
	return ( x ^ (symplectic_form(x,h,n) * h) );
}

vector<binvec> coord_matrix(const unsigned n) {
	vector<binvec> M (2*n,0);
	for(unsigned i=0; i<n; i++) {
		M[2*i] = set_bit(0, i, 2*n);
		M[2*i+1] = set_bit(0, n+i, 2*n);
	}
	return M;
}

binvec zx_to_product_coordinates(const binvec x, const unsigned n) {
	return matrix_vector_prod_mod2(coord_matrix(n), x);
}

// counts weights of the Pauli operator that corresponds to the phase space point a
// note that the order is 1,X,Z,Y
vector<unsigned> count_weights(const binvec a, const unsigned n) {
	binvec ai;
	vector<unsigned> weights(4,0);
	for(unsigned i=0; i<n; i++) {
		ai = get_bits(a, 2*i, 2*i+1, 2*n);
		weights.at(ai) += 1;
	}
	return weights;
}


// ---- Symplectic group generation

int phase_function(const binvec a, const binvec *S, const unsigned n) {
	// represents the phase function g of the trivial Clifford representative C with [C]=S in Sp(2n) = C(n)/P(n). It describes the action of C as follows:
	//		C W_a C^\dagger = (-1)^{g(a)} W_{Sa} 

	if(popcount(a) <= 1) {
		// hence, either a == 0 or a is a basis vector in phase space
		return 0;
	}

	// will be used for a basis expansion of a
	// we will actually only store set bits in a
	vector<binvec> av;
	av.reserve(2*n);
	vector<binvec> bv;
	bv.reserve(2*n);

	binvec aa = a;
	while(aa != 0) {
		// push back least significant bit
		av.push_back(aa & ~(aa-1));
		bv.push_back(matrix_vector_prod_mod2(S,aa & ~(aa-1),2*n));
		// reset it in aa
		aa = aa & (aa-1);
	}
	// compute phis
	return (phi(bv,n) - phi(av,n));

}

int phase_function(const binvec a, const vector<binvec> &S) {
	// represents the phase function g of the trivial Clifford representative C with [C]=S in Sp(2n) = C(n)/P(n). It describes the action of C as follows:
	//		C W_a C^\dagger = (-1)^{g(a)} W_{Sa} 

	if(popcount(a) <= 1) {
		// hence, either a == 0 or a is a basis vector in phase space
		return 0;
	}

	unsigned n = S.size()/2;

	// will be used for a basis expansion of a
	// we will actually only store set bits in a
	vector<binvec> av;
	av.reserve(2*n);
	vector<binvec> bv;
	bv.reserve(2*n);

	binvec aa = a;
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

vector<binvec> find_transvection(const binvec x, const binvec y, const unsigned n) {
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
		return vector<binvec>(3,0U);
	}
	
	if(x == y) {
		// return 0
		return vector<binvec>(3,0U);
	} 

	if(symplectic_form(x,y,n) == 1) {
		// return h = x + y (bitwise addition)
		return vector<binvec>({0, x^y,1});
	}
	else {
		binvec z = 0; 
		binvec xx,yy;
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
				return vector<binvec>({x^z,y^z,2});
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
		binvec zz = 0;

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
		return vector<binvec>({x^z,y^z,3});
	}

}

vector<binvec> generate_symplectic_matrix(const binvec i, const unsigned n) {
	// Generate the i-th element in Sp(2n, Z_2) using transvections
	// Actually, returns the transposed element.

	unsigned nn = 2*n;
	unsigned s = (1<<nn)-1;

	// first and second basis vectors e1=(1,0,0,...,0), e2=(0,1,0,...,0)
	binvec e1 = set_bit(0,0,nn);
	binvec e2 = set_bit(0,1,nn);

	// compute phase space points
	binvec f1 = invert_bits((i%s)+1,nn);
	binvec b = invert_bits(i/s,nn-1);
	binvec eprime = e1 ^ b;
	eprime = reset_bit(eprime,1,nn); // result is (1,0,b1,b2,b3,...,b2n-1)

	// compute and apply transvections
	vector<binvec> T = find_transvection(e1,f1,n);
	
	binvec h0 = transvection(T.at(1),eprime,n);
	h0 = transvection(T.at(0),h0,n);

	// check 0th bit of b
	if(get_bit(b,0,nn-1) == 1) {
		f1 = 0;
	}

	// unsigned f2 = transvection(T.at(1),e2,n);
	// f2 = transvection(T.at(0),f2,n);
	// f2 = transvection(Tprime.at(1),f2,n);
	// f2 = transvection(Tprime.at(0),f2,n);

	vector<binvec> g;
	vector<binvec> id = {2,1}; // 2x2 identity matrix

	if(n == 1) {
		g = id;
	}
	else {
		// recursive call
		g = direct_sum(id, generate_symplectic_matrix((i/s) >> (nn-1),n-1));
	}

	// finally, apply the transvections to each row in g
	for(binvec j=0; j<nn; j++) {
		g[j] = transvection(T.at(1),g[j],n);
		g[j] = transvection(T.at(0),g[j],n);
		g[j] = transvection(h0,g[j],n);
		g[j] = transvection(f1,g[j],n);
	}

	return g;
}

vector<int> projected_liouville_matrix(const vector<binvec> &S, const binvec a=0) {
	// Computes the Liouville matrix of the "canonical" Clifford representative corresponding to the symplectic matrix S, i.e. the one that acts on Z/X Paulis as
	//		C:  W(e_j) ---> +W(Se_j)
	// Its action on an arbitrary Pauli operator W(b) is
 	//		C:  W(b) ---> (-1)^g(b) W(Sb)
 	// where g(b) is some phase function determined by S.
 	//
 	// Note that this function returns the projection of the Liouville matrix
 	//
	// Optionally, if 0<=a<4^n is specified, we include a term (-1)^[a,b] in the phase, which corresponds to the action of the n-qubit Pauli group P_n. 

	unsigned n = S.size()/2;
	unsigned m = pow(4,n);
	unsigned i;
	binvec Sb;
	vector<int> L(4*m-1,0);
	int ph;

	for(unsigned b=1; b<m; b++) {
		// compute the action on b
		Sb = matrix_vector_prod_mod2(S, b);

		// check if the the last 2*n-2 bits of Sb agree with the last 2*n-2 bits of b
		if(get_bits(b,2,2*n-1,2*n) == get_bits(Sb,2,2*n-1,2*n)) {

			// they do, now get the first two bits in Sb and append them in front of b
			// this corresponds to the index in the projected Liouvillean
			i = b ^ ( get_bits(Sb,0,1,2*n) << 2*n );

			// this is the Pauli phase
			L.at(i-1) = pow(-1, symplectic_form(a,b,n) );

			// "canonical" phase
			int ph = mod(phase_function(b,S),4);
			assert(ph == 0 || ph == 2);
			if(ph == 2) {
				L.at(i-1) *= -1;
			}
		}
		// every other entry is zero 
	}
	return L;
}


vector<int> liouville_matrix(const vector<binvec> &S, const binvec c=0) {
	// Computes the Liouville matrix of the "trivial" Clifford representative corresponding to the symplectic matrix S, i.e. the one that acts Z/X Paulis as
	//		C:  W(e_j) ---> +W(Se_j)
	// Its action on an arbitrary Pauli operator W(b) is
 	//		C:  W(b) ---> (-1)^g(b) W(Sb)
 	// where g(b) is some phase function determined by S.
 	//
	// Optionally, if 0<=c<2*n is specified, we include a term (-1)^[c,b] in the phase, which corresponds to the action of the n-qubit Pauli group P_n. 

	unsigned n = S.size()/2;
	unsigned m = 1 << 4*n;	// pow(16,n)
	binvec a,b;
	vector<int> L(m,0);
	int ph;

	for(unsigned i=0; i<m; i++) {
		// get the first and last 2n bits of i
		a = get_bits(i,0,2*n-1,4*n);
		b = get_bits(i,2*n,4*n-1,4*n);

		// check if a = Sb
		if( a == matrix_vector_prod_mod2(S, b) ) {
			// they do
			L.at(i) = pow(-1, symplectic_form(c,b,n) );

			int ph = mod(phase_function(b,S),4);
			assert(ph == 0 || ph == 2);
			if(ph == 2) {
				L.at(i) *= -1;
			}			
		}
	}
	return L;
}


#endif