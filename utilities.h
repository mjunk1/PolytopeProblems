#ifndef UTILITIES_H
#define UTILITIES_H
#include <initializer_list>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>


using namespace std;

/* Utilities header file 
	For usage with the noise propagation code. 

	Caution!
	--------
	Due to the use of __builtin_popcount(), this code can only be used with the GNU compiler suite!
	Furthermore, the code uses at least the C++14 standard.
*/




// ----- indexing routines

/* Computes the multi index for a given linear index assuming row-major order 
	The array dimensions are assumed to be N_1 x N_2 x ... x N_k and are specified by arr_dim = k and the array arr_ranges = {N_1, ..., N_k}. The result is stored in indices (should be allocated before the use of this function!)
*/
unsigned get_multi_index(const unsigned arr_dim, const unsigned* arr_ranges, const unsigned index, unsigned* indices) {

	unsigned mod;
	unsigned ind = index;
	if(arr_dim> 1) {
		for(unsigned i=arr_dim-1; i>0; i--) {
			mod = ind % arr_ranges[i];
			indices[i] = mod;
			ind = (ind - mod) / arr_ranges[i];
		}
	}
	indices[0] = ind;

	return 0;
}

unsigned get_multi_index(const std::vector<unsigned> &arr_ranges, const unsigned index, vector<unsigned> &indices) {
	indices.reserve(arr_ranges.size());
	return get_multi_index(arr_ranges.size(), arr_ranges.data(), index, indices.data());
}

// overload where all N_i are assumed to be the same
unsigned get_multi_index(const unsigned arr_dim, const unsigned arr_range, const unsigned index, unsigned* indices) {
	vector<unsigned> vec (arr_dim, arr_range);
	return get_multi_index(arr_dim, vec.data(), index, indices);
}

unsigned get_multi_index(const unsigned arr_dim, const unsigned arr_range, const unsigned index, vector<unsigned> &indices) {
	vector<unsigned> vec (arr_dim, arr_range);
	indices.reserve(arr_dim);
	return get_multi_index(arr_dim, vec.data(), index, indices.data());
}

// inverse stuff

unsigned get_linear_index(const unsigned arr_dim, const unsigned* arr_ranges, const unsigned* indices) {
	unsigned index = 0;
	unsigned jprod;
	for(unsigned i=0; i<arr_dim; i++) {//change order?
		jprod = 1;
		for(unsigned j=i+1; j<arr_dim; j++) {
			jprod *= arr_ranges[j];
		}
		jprod *= indices[i];
		index += jprod;
	}

	return index;
}

unsigned get_linear_index(const unsigned arr_dim, const unsigned* arr_ranges, const initializer_list<unsigned> indices){
	return get_linear_index(arr_dim, arr_ranges, const_cast<unsigned* >(indices.begin()));
}

// overloads where all N_i are assumed to be the same
unsigned get_linear_index(const unsigned arr_dim, const unsigned arr_range, const unsigned* indices) {
	vector<unsigned> vec (arr_dim, arr_range);
	return get_linear_index(arr_dim, vec.data(), indices);
}

unsigned get_linear_index(const unsigned arr_dim, const unsigned arr_range, const initializer_list<unsigned> indices) {
	vector<unsigned> vec (arr_dim, arr_range);
	return get_linear_index(arr_dim, vec.data(), const_cast<unsigned* >(indices.begin()));
}



// ----- computational routines

unsigned get_number_of_lines(string filename) {
	unsigned n=0;
	string line;
	fstream fin(filename, ios::in);
	if(fin.is_open())
		while(getline(fin, line))
			++n;
	fin.close();
	return n;
}

unsigned popcount(unsigned i) {
	return __builtin_popcount(i);
}

unsigned parity(unsigned i) {
	return __builtin_parity(i);
}



// ---- Phase space methods

// This is the convolution p3 = p1 * p2 in (Z_2)^n (addition mod 2)
int convolve_mod2(double *p1, double *p2, double *p3, const unsigned n) {
  unsigned len = pow(2U,n);
  unsigned i,j;

  for(i=0; i<len; i++) {
    p3[i] = 0;
    
    for(j=0; j<len; j++) {
      p3[i] += p1[j] * p2[i^j];
    } 
  }
  return 0;
}

// Computes the matrix vector product Ax over Z_2. It is assumed that a vector in Z_2^n is represented as a bitstring encoded as an unsigned integer. In this representation, every row of the matrix A corresponds to an unsigned integer, hence A is an unsigned array. The ouput is again an unsigned integer.
unsigned matrix_vector_prod_mod2(const unsigned *A, const unsigned x, const unsigned n) {
	unsigned b = 0;

	/* This performs the formula:
		b_i = sum_{j=1}^n A_ij x_j
		
		In the bitstring representation, this means that we have to multiply A[i] with x bitwise, afterwards adding up the resulting bits mod 2. The last step is just counting the number of ones and taking its modulus (here, called parity). Finally, we set the i-th bit of b to this result.
		Note that, in contrast to the usual convention, we count bits from the left since we see bitstrings as binary vectors. This means that the last bit is the least significant one.
	*/
	for (unsigned i=0; i<n; i++) 
		b ^= ((-parity( A[i]&x )) ^ b) & (1U << (n-i-1U)); 

	return b;	
}

int symplectic_transform(double *pin, double *pout, const unsigned *S, const unsigned nqubits) {
	unsigned j;
	for(unsigned i=0; i<pow(4,nqubits); i++) {
		j = matrix_vector_prod_mod2(S, i, 2U*nqubits);
		pout[j] = pin[i];
	}

	return 0;
}

int write_pdist(fstream &fout, double *p, unsigned nqubits) {
	
	unsigned len = pow(4,nqubits);
	unsigned bit;
	unsigned n = (unsigned)(2*nqubits);

	for(unsigned i=0; i<len; i++) {
		for(unsigned j=0; j<n; j++) {
			bit = (i >> (n-j-1U)) & 1U; // get the j-th bit of i
			fout << bit << " "; 
		}
		fout << p[i] << endl;		
	}

	return 0;
}


#endif