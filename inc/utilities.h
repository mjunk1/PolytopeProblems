#ifndef UTILITIES_H
#define UTILITIES_H
#include <initializer_list>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <eigen3/Eigen/Sparse>
#include <tuple>
#include <cblas.h>


using namespace std;
using namespace Eigen;

/* Utilities header file 
	For usage with the noise propagation code. 

	Caution!
	--------
	Due to the use of __builtin_popcount(), this code can only be used with the GNU compiler suite!
	Furthermore, the code uses at least the C++14 standard.
*/



// ----- helpers

// samples a random discrete probability distribution of length n
// corresponds to sampling a random point out of the (n-1)-simplex
vector<double> random_distribution(unsigned n) {
	random_device rd; 
	mt19937 gen(rd()); 
	uniform_real_distribution<> dist (0,1); // uniform distribution on [0,1]

	vector<double> p(n);
	double N = 0.;

	// corresponds to n+1 random exponentially distributed numbers
	for(unsigned i=0; i<n; i++) {
		p.at(i) = - log(dist(gen)); 
		N += p.at(i);
	}

	// normalise
	for(unsigned i=0; i<n; i++) {
		p.at(i) /= N;
	}

	return p;
}

// get all descending partitions of a positive number n, i.e. a list { (n), (n-1,1), ... } of a list of positive integers l_i such that 
// n = l_1 + ... + l_m (m=1,...,n), and
// n >= l_1 >= l_2 >= ... >= l_m >= 1
// The first element is n itself.
vector<vector<unsigned>> get_partitions(unsigned n) {
	if( n == 1 ) {
		// return the list { (1) }
		return vector<vector<unsigned>> (1, vector<unsigned>(1,1) );
	}

	vector<vector<unsigned>> ret = { vector<unsigned>(1,n) };
	vector<vector<unsigned>> ret2;
	
	for(unsigned k=1; k<n; k++) {
		ret2 = get_partitions(k);
		for(auto p : ret2) {
			if(p.at(0) <= n-k) {
				p.insert(p.begin(), n-k);
				ret.push_back(p);
			}
		}
	}

	return ret;
}


// linear algebra

// careful! This function applies the transposed matrix A^T to x ...
// in my case, this is convenient since ColMajor layout is faster
vector<double> my_matrix_vector_prod(const vector<double> &A, const vector<double> &x) {
	unsigned N = x.size(); // N = 4^n
	vector<double> y(N,0);

	// transform
	cblas_dgemv(CblasColMajor,CblasNoTrans,N,N,1.0,A.data(),N,x.data(),1,0.0,y.data(),1);

	return y;
}

unsigned factorial(const unsigned n) {
	unsigned f = 1;
	for(unsigned nn=n; nn>0; --nn) {
		f *= nn;
	}
	return f;
}

unsigned binomial_coeff(const unsigned n, const unsigned k) {
    vector<unsigned> C(k+1,0);
 
    C[0] = 1;
 
    for (unsigned i = 1; i <= n; i++) {
        // Compute next row of pascal triangle using
        // the previous row
        for (unsigned j = min(i, k); j > 0; j--)
            C[j] = C[j] + C[j-1];
    }
    return C[k];
}


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

unsigned get_linear_index(const unsigned arr_dim, const unsigned* arr_ranges, const initializer_list<unsigned> &indices){
	return get_linear_index(arr_dim, arr_ranges, const_cast<unsigned* >(indices.begin()));
}

unsigned get_linear_index(const unsigned arr_dim, const initializer_list<unsigned> &arr_ranges, const initializer_list<unsigned> &indices){
	return get_linear_index(arr_dim, const_cast<unsigned* >(arr_ranges.begin()), const_cast<unsigned* >(indices.begin()));
}

// overloads where all N_i are assumed to be the same
unsigned get_linear_index(const unsigned arr_dim, const unsigned arr_range, const unsigned* indices) {
	vector<unsigned> vec (arr_dim, arr_range);
	return get_linear_index(arr_dim, vec.data(), indices);
}

unsigned get_linear_index(const unsigned arr_dim, const unsigned arr_range, const initializer_list<unsigned> &indices) {
	vector<unsigned> vec (arr_dim, arr_range);
	return get_linear_index(arr_dim, vec.data(), const_cast<unsigned* >(indices.begin()));
}


// indexing for symmetric matrices, upper triangle storage
inline unsigned get_ut_index(const unsigned i, const unsigned j, const unsigned n) {
	return ( (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1 );
}

inline unsigned get_ut_row(const unsigned k, const unsigned n) {
	return n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5);
}

inline unsigned get_ut_col(const unsigned i, const unsigned k, const unsigned n) {
	return k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
}

// this is for column-wise storage
inline unsigned get_ut_index2(const unsigned i, const unsigned j, const unsigned n) {
	return ( i + j*(j-1)/2 );
}

unsigned get_ut_col2(const unsigned k, const unsigned n) {
	// return floor( 0.5 + 0.5*sqrt(8*k+1) );
	for(unsigned j=1; j<n; j++) {
		if(k >= j*(j-1)/2 && k < j*(j+1)/2) {
			return j;
		}
	}
}

inline unsigned get_ut_row2(const unsigned j, const unsigned k, const unsigned n) {
	return (k - j*(j-1)/2);
}

// conversion of upper-triangle storage to dense storage for symmetric matrices



// ----- bit manipulation

unsigned popcount(unsigned i) {
	return __builtin_popcount(i);
}

unsigned parity(unsigned i) {
	return __builtin_parity(i);
}


// get j-th bit of b, counted from the end of the bitstring, i.e. from the left
inline unsigned get_bit(unsigned b, unsigned j, unsigned n) {
	return ((b >> ((n)-(j)-1U)) & 1U);
}

// get j-th bit of b, counted from the from the right
inline unsigned get_bit2(unsigned b, unsigned j) {
	return ((b >> j) & 1U);
}

inline unsigned get_bits(unsigned b, unsigned j, unsigned k, unsigned n) {
	return ((b >> (n-k-1U)) & ((1U << (k-j+1U))-1));
}

inline unsigned set_bit(unsigned b, unsigned j, unsigned n) {
	return (b | (1 << (n-j-1U)));
}

inline unsigned reset_bit(unsigned b, unsigned j, unsigned n) {
	return (b & (~(1 << (n-j-1U))));
}

string write_bits(unsigned b, unsigned n) {
	string str;

	for(unsigned j=0; j<n; j++) {
		str += to_string(get_bit(b,j,n)); 
	}

	return str;
}

string write_trits(unsigned t, unsigned n) {
	string str;

	vector<unsigned> trits(n,0);
	get_multi_index(n, 3, t, trits);

	for(unsigned j=0; j<n; j++) {
		str += to_string( trits.at(j) ); 
	}

	return str;
}

unsigned invert_bits(const unsigned b, const unsigned n) {
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


// ----- input / output

unsigned get_number_of_lines(string filename) {
	unsigned n=0;
	string line;
	fstream fin(filename, ios::in);
	if(fin.is_open()) {
		while(getline(fin, line)) {
			++n;
		}
	}
	fin.close();
	return n;
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