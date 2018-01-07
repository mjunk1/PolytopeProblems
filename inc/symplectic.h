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
	// Computes the phase that is appearing in definition of Weyl operators:
	//		W(a) = i^{-\eta(a,a)} Z(a_z)X(a_x)
	// \eta is actually a Z_4-valued function.
	unsigned ret = 0;

	for(unsigned i=0; i<n; i++) {
		ret += ( get_bit(a,2U*i+1U,2U*n) * get_bit(b,2U*i,2U*n) );
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
	//		\phi(a_1,...,a_m) = \sum_{j=1}^{m-1} \phi(a_j, \sum_{k=j+1}^m a_k)
	//
	unsigned ret = 0;
	unsigned m = a.size();
	unsigned aa;

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
vector<unsigned> direct_sum(const vector<vector<unsigned>> Alist) {
	vector<unsigned> temp = Alist.at(0);
	for(unsigned i=1; i<Alist.size(); i++) {
		temp = direct_sum(temp, Alist.at(i));
	}

	return temp;
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

vector<unsigned> coord_matrix(const unsigned n) {
	vector<unsigned> M (2*n,0);
	for(unsigned i=0; i<n; i++) {
		M[2*i] = set_bit(0, i, 2*n);
		M[2*i+1] = set_bit(0, n+i, 2*n);
	}
	return M;
}

unsigned zx_to_product_coordinates(const unsigned x, const unsigned n) {
	return matrix_vector_prod_mod2(coord_matrix(n), x);
}

// ---- Symplectic group generation

int phase_function(const unsigned a, const unsigned *S, const unsigned n) {
	// represents the phase function g of the trivial Clifford representative C with [C]=S in Sp(2n) = C(n)/P(n). It describes the action of C as follows:
	//		C W_a C^\dagger = (-1)^{g(a)} W_{Sa} 

	if(popcount(a) <= 1) {
		// hence, either a == 0 or a is a basis vector in phase space
		return 0;
	}

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
		bv.push_back(matrix_vector_prod_mod2(S,aa & ~(aa-1),2*n));
		// reset it in aa
		aa = aa & (aa-1);
	}
	// compute phis
	return (phi(bv,n) - phi(av,n));

}

int phase_function(const unsigned a, const vector<unsigned> S) {
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
	unsigned Sb;
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


// ---- Stabiliser state generation

// This generates the i-th symmetric nxn matrix and returns the representation of the corresponding graph state given in terms of a basis for the Lagrangian subspace
// The basis vectors are given by the column vectors of the matrix
// 			(  A  )
// 		B = ( --- )			(in (z_1,...,z_n,x_1,...,x_n) coordinates)
//			(  I  )
// where A is a nxn symmetric matrix and I is the nxn identity matrix. Note that since we are working in product coordinates (z_1,x_1,z_2,x_2,...,z_n,x_n), the function will return a permuted version of these vectors.
// The i-th vector from Z_2^{n(n+1)/2} is simply given by the binary representation of i
vector<unsigned> generate_gs_Lagrangian(const unsigned i, const unsigned n) {
	// reserve memory
	vector<unsigned> B(n,0);
	unsigned N = n*(n+1)/2;
	vector<unsigned> M = coord_matrix(n);

	for(unsigned j=0; j<n; j++) {
		// save the j-th basis vector
		for(unsigned k=0; k<n; k++) {
			if( get_bit(i, get_symmetric_index(k,j), N) == 1) {
				B[j] = set_bit(B[j], k, 2*n);
			}
		}

		// identity matrix
		B[j] = set_bit(B[j], n+j, 2*n);	

		// change to product coordinates
		B[j] = matrix_vector_prod_mod2(M, B[j]);	
	}

	return B;
}

bool in_Lagrangian(const unsigned x, const vector<unsigned> B, const unsigned n) {
	unsigned test = 0;
	for(unsigned v : B) {
		test += symplectic_form(x,v,n);
	} 
	return (test == 0);
}


// Generates the stabiliser state rho(B,s) that corresponds to the Lagrangian with basis B and a choice of signs (-1)^{s_i} where the bits s_i \in \Z_2 are specified in the length-n-bitstring s.
// The stabiliser state is given in the Pauli basis {W(a)|a\in\Z_2^{2n}} where the only non-vanishing components rho_a corresponds to points a in the Lagrangian subspace and can thus be written in the basis B = (b_1,...,b_n) as
//		a = \sum_{i=1}^n a_i b_i.
// These components are explicitly given as
//		rho_a = (-1)^{\sum_{i=1}^n a_i s_i} i^{\phi(a_1 b_1,...,a_n b_n)},
// where \phi is the phase function appearing in the composition law of Pauli operators (see above).
vector<double> generate_stabiliser_state(const vector<unsigned> B, const unsigned s) {
	unsigned n = B.size();
	vector<double> state (pow(4,n),0);
	vector<unsigned> avec (n);
	unsigned phase, idx;

	// loop over all points a in the Lagrangian
	for(unsigned a=0; a<pow(2,n); a++) {
		idx = 0;

		for(unsigned i=0; i<n; i++) {
			avec.at(i) = get_bit(a,i,n) * B.at(i);
			idx ^= avec.at(i); // this gives in the end the coordinates of a w.r.t. the canonical basis
		}

		// compute phi
		phase = mod(phi(avec,n),4);
		assert(phase==0 || phase==2);
		if(phase == 2) {
			phase = -1;
		}
		else {
			phase = 0;
		}

		phase += parity(a&s); // inner product of a and s

		state.at(idx) = pow(-1,phase);
	}

	return state;
}



vector<vector<double>> generate_stabiliser_states(const unsigned n) {
	// estimated size (overcounted!)
	unsigned N = pow(2,n*(n+1)/2);
	unsigned M = pow(6,n);
	unsigned S = pow(2,n);
	unsigned size = N*M*S;

	// needed later
	vector<vector<double>> states ( size, vector<double> (pow(4,n),0) );
	vector<unsigned> B;
	vector<unsigned> SB (n);

	// generate 1-qubit symplectic group
	vector<vector<unsigned>> symp_group (6);

	for(unsigned i=0; i<6; i++) {
		symp_group.at(i) = generate_symplectic_matrix(i,1);
	}

	// generate local n-qubit symplectic group
	vector<vector<unsigned>> loc_symp_group (M);
	vector<unsigned> indices (n);
	vector<vector<unsigned>> Slist (n);
	for(unsigned i=0; i<M; i++) {
		get_multi_index(n, 6, i, indices);
		for(unsigned j=0; j<n; j++) {
			Slist.at(j) = symp_group.at(indices.at(j));
		}
		loc_symp_group.at(i) = direct_sum(Slist);
	}

	// graph state loop
	for(unsigned i=0; i<N; i++) {
		B = generate_gs_Lagrangian(i,n);

		// local symplectic orbit
		for(unsigned j=0; j<M; j++) {
			// transform the basis vectors with the j-th matrix
			for(unsigned k=0; k<n; k++) {
				SB.at(k) = matrix_vector_prod_mod2(loc_symp_group.at(j), B.at(k));
			}

			// generate stabiliser states for all possible sign choices
			for(unsigned s=0; s<S; s++) {
				states.at(get_linear_index(3,{N,M,S},{i,j,s})) = generate_stabiliser_state(SB, s);
			}
		}
	}

	// eliminate duplicates
	sort(states.begin(), states.end());
	auto last = unique(states.begin(), states.end());
	int d = distance(last, states.end());
	states.erase(last,states.end());

	cout << "Found " << d << " duplicates. " << size-d << " elements left." << endl;

	return states;
}


#endif