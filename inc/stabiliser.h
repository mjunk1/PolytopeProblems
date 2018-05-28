#ifndef STABILISER_H
#define STABILISER_H
#include <initializer_list>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <set>
#include <fstream>
#include <utility>
#include <cstdio>
#include <sstream>

#ifndef SYMPLECTIC_H
#include "symplectic.h"
#endif

#ifndef GLPKCONVEXSEPARATION_H
#include "GLPKConvexSeparation.h"
#endif

// --- data type for points
typedef LabelledPoint<int> LabelledState;



// -----------------------------------
// ---- Helper functions
// -----------------------------------

// Liouville-Pauli representation of T gate
vector<double> Tmatrix = {1,0,0,0,0,1/sqrt(2),0,-1/sqrt(2),0,0,1,0,0,1/sqrt(2),0,1/sqrt(2)};


vector<double> rotation_matrix(const unsigned n) {
	// build transformation matrix A = (T^t)^{\otimes n}
	// actually the transposed ...
	unsigned N = pow(4,n);
	vector<double> A(N*N);
	vector<unsigned> ai(n);
	vector<unsigned> bi(n);

	for(unsigned a=0; a<N; a++) {
		// decompose indices	
		// get_multi_index(n, 4, a, ai);
		for(unsigned i=0; i<n; i++) {
			ai.at(i) = get_bits(a, 2*i, 2*i+1, 2*n);
		}

		for(unsigned b=0; b<N; b++) {
			// decompose indices	
			// get_multi_index(n, 4, b, bi);
			for(unsigned i=0; i<n; i++) {
				bi.at(i) = get_bits(b, 2*i, 2*i+1, 2*n);
			}

			A.at(a*N+b) = 1;
			for(unsigned i=0; i<n; i++) {
				A.at(a*N+b) *= Tmatrix.at(ai.at(i)*4+bi.at(i));
			}
		}
	}

	return A;
}

vector<unsigned> index_set(const unsigned k, const unsigned n) {
	// returns the index set corresponding to the subset of phase space points that have the form
	// 		a = (0,1,0,1,...,0,1,0,0,...,0,0)			((0,1) is appearing exactly k times)
	// and pairwise permutations thereof.
	assert(k<=n && n>0);

	if(k==0) {
		return vector<unsigned>({0});
	}

	if(n==k) {
		unsigned r = 0;
		for(unsigned i=0; i<n; i++) {
			// if r = 010101, then (r<<2)^1 = (01010100 )^1 = 01010101
			r = (r << 2)^1;
		}
		return vector<unsigned> ({r});
	}

	// first index (0,0)
	vector<unsigned> iset0 = index_set(k, n-1);

	// first index (0,1)
	vector<unsigned> iset1 = index_set(k-1, n-1);
	for(unsigned i=0; i<iset1.size(); i++) {
		// add 01 in front of every index
		iset1.at(i) = set_bit(iset1.at(i), 1, 2*n);
	}

	// merge the two sets
	iset0.reserve(iset0.size()+iset1.size());
	iset0.insert(iset0.end(), iset1.begin(), iset1.end());

	return iset0;
}

vector<unsigned> index_set_y(const unsigned k, const unsigned n, const bool use_zeros=true) {
	// same as index_set but with (1,1) instead of (0,1)
	// use_zeros is a switch of wether to fill the rest of the string with zeros or (0,1)'s
	assert(k<=n && n>0);

	if(k==0) {
		if(use_zeros == true) {
			return vector<unsigned>({0});
		}
		else {
			unsigned r = 0;
			for(unsigned i=0; i<n; i++) {
				// if r = 010101, then (r<<2)^1 = (01010100 )^1 = 01010101
				r = (r << 2)^1;
			}
			return vector<unsigned> ({r});
		}
	}

	if(n==k) {
		unsigned r = pow(2,2*n) - 1; // r = 1111...11
		return vector<unsigned> ({r});
	}

	// first index (0,0)
	vector<unsigned> iset0 = index_set_y(k, n-1, use_zeros);
	if(use_zeros == false) {
		for(unsigned i=0; i<iset0.size(); i++) {
			// add 01 in front of every index
			iset0.at(i) = set_bit(iset0.at(i), 1, 2*n);
		}
	}

	// first index (1,1)
	vector<unsigned> iset1 = index_set_y(k-1, n-1, use_zeros);
	for(unsigned i=0; i<iset1.size(); i++) {
		// add 11 in front of every index
		iset1.at(i) = set_bit(iset1.at(i), 0, 2*n);
		iset1.at(i) = set_bit(iset1.at(i), 1, 2*n);
	}

	// merge the two sets
	iset0.reserve(iset0.size()+iset1.size());
	iset0.insert(iset0.end(), iset1.begin(), iset1.end());

	return iset0;
}

vector<unsigned> index_set2(const unsigned kk, const unsigned k, const unsigned n) {
	// returns the index set corresponding to the subset of phase space points that have the form
	// 		a = (0,1,0,1,...,0,1,1,1,...,1,1,0,0,...,0,0)			((0,1) is appearing exactly k-kk times, (1,1) is appearing kk times)
	// and pairwise permutations thereof.
	assert(kk<=k && k<=n && n>0);

	if(kk==0) {
		// just the ordinary index set
		return index_set(k,n);
	}

	if(kk==k) {
		// now only (1,1) appearing, which is ordinary index set with replacing (0,1)->(1,1)
		return index_set_y(k,n);
	}

	if(n==k) {
		// no zeros at all
		return index_set_y(kk,n,false);
	}

	// first index (0,0)
	vector<unsigned> iset0 = index_set2(kk, k, n-1);

	// first index (0,1)
	vector<unsigned> iset1 = index_set2(kk, k-1, n-1);
	for(unsigned i=0; i<iset1.size(); i++) {
		// add 01 in front of every index
		iset1.at(i) = set_bit(iset1.at(i), 1, 2*n);
	}

	// first index (1,1)
	vector<unsigned> iset2 = index_set2(kk-1, k-1, n-1);
	for(unsigned i=0; i<iset2.size(); i++) {
		// add 11 in front of every index
		iset2.at(i) = set_bit(iset2.at(i), 0, 2*n);
		iset2.at(i) = set_bit(iset2.at(i), 1, 2*n);
	}

	// merge the three sets
	iset0.reserve(iset0.size()+iset1.size()+iset2.size());
	iset0.insert(iset0.end(), iset1.begin(), iset1.end());
	iset0.insert(iset0.end(), iset2.begin(), iset2.end());

	return iset0;
}

// -----------------------------------
// ---- Magic state representation
// -----------------------------------

// Representation of H state
vector<double> H = {1,1/sqrt(2),0,1/sqrt(2)};

vector<double> H_state(const unsigned n, const double p = 0) {
	vector<double> ret (pow(4,n), 0.);
	vector<unsigned> ind (n);

	ret.at(0) = 1;

	for(unsigned i=1; i<pow(4,n); i++) {
		get_multi_index(n, 4, i, ind);
		ret.at(i) = 1-p;
		for(unsigned j=0; j<n; j++) {
			ret.at(i) *= H.at(ind.at(j));
		}
	}

	return ret;
}

vector<double> H_state_rotated(const unsigned n) {
	vector<double> ret (pow(4,n), 0.);

	for(unsigned k=0; k<=n; k++) {
		vector<unsigned> aset = index_set(k,n);
		for(auto a : aset) {
			ret.at(a) = 1;
		}
	}

	return ret;
}

// H state in the (renormalised) number / Fock basis
vector<double> H_state_nb(const unsigned n, const double p = 0) {
	vector<double> ret (n);
	for(unsigned k=1; k<n+1; k++) {
		ret.at(k-1) = (1-p)*pow(sqrt(2),k)*binomial_coeff(n,k);
	}
	return ret;
}

class ProjectedNoisyHState: public PointGenerator {
protected:
	unsigned _n;

public:
	ProjectedNoisyHState(unsigned n) {
		_n = n;
	}

	vector<double> operator() (double p) {
		return H_state_nb(_n,p);
	}

};

class NoisyHState: public PointGenerator {
protected:
	unsigned _n;

public:
	NoisyHState(unsigned n) {
		_n = n;
	}

	vector<double> operator() (double p) {
		return H_state(_n,p);
	}

};


// -----------------------------------
// ---- Stabiliser generation
// -----------------------------------

// --- graph routines

// read graph6 format and returns the adjacency matrix in as a vector of bitstrings representing the rows
// IMPORTANT: this implementation is quick and dirty and hence only works for graphs with less than 63 vertices (fine for me)
vector<binvec> graph6_to_adj_mat(const string g6) {
	// The input string g6 is assumed to hold a graph6 representation of the graph in ASCII-encoded form, i.e. it's a sequence of ASCII characters

	// interprete g6 as integer
	unsigned n = (unsigned)g6.at(0)-63;
	unsigned l = g6.size();
	unsigned N = n*(n-1)/2;
	unsigned i,j;
	vector<binvec> ret (n,0);
	vector<binvec> x(l-1,0);

	if(N > 0) {
		for(i=1; i<l; i++) {
			x.at(i-1) = (binvec)(g6.at(i)-63);
		}

		// read out the bits and fill the adjacency matrix		
		for(unsigned k=0; k<N; k++) {
			if(get_bit(x.at(k/6), k%6, 6) == 1) {
				j = get_ut_col2(k,n);
				i = get_ut_row2(j,k,n);
				ret.at(i) = set_bit(ret.at(i), j, n);
				ret.at(j) = set_bit(ret.at(j), i, n);
			}
		}
	}

	return ret;
}

// inverse operation
string adj_mat_to_graph6(const vector<binvec> &A) {
	unsigned n = A.size();

	assert(n < 63);

	if(n == 0 || n == 1) {
		return string(1, 63+n);
	}

	// convert upper triangle of A to a column-wise ordered binary vector which is split in blocks of 6 bits each
	unsigned N = n*(n-1)/2;
	unsigned m = N/6;
	unsigned k;

	if(N%6 != 0) {
		++m;
	}
	vector<binvec> x(m,0);

	for(unsigned i=0; i<n; i++) {
		for(unsigned j=i+1; j<n; j++) {
			if(get_bit(A.at(i),j,n) == 1) {
				k = get_ut_index2(i,j,n);
				x.at(k/6) = set_bit(x.at(k/6), k%6, 6);
			}
		}
	}

	// now convert x to a string of ASCII characters
	// first character encodes n
	string ret (m+1, 63);	
	ret.at(0) = (char)(n+63);

	for(unsigned k=0; k<m; k++) {
		ret.at(k+1) = (char)( x.at(k) + 63 );
	}

	return ret;
}

unsigned get_order(const string g6) {
	return (unsigned)g6.at(0)-63;
}

unsigned get_order_from_file(const string file) {
	// open the file
	fstream fin(file, ios::in);
	string row;

	if(fin.is_open()) {
		fin >> row;
	} 
	else {
		cout << "Couldn't open file " << file << endl;
		return 1;
	}
	fin.close();

	return get_order(row);
}



// --- graph states

// represents generator matrix of a graph state
vector<binvec> graph_Lagrangian(const vector<binvec> &A, const unsigned n) {
	// This returns the representation of a graph state given in terms of a basis for the Lagrangian subspace. The graph of the graph state is specified by its adjacency matrix, given as a vector A of bitstring representing the row vectors in Z_2^n.
	// The basis vectors are given by the column vectors of the matrix
	// 			(  A  )
	// 		B = ( --- )			(in (z_1,...,z_n,x_1,...,x_n) coordinates)
	//			(  I  )
	// where A is the adjacency matrix, now interpreted as nxn symmetric matrix and I is the nxn identity matrix. Note that since we are working in product coordinates (z_1,x_1,z_2,x_2,...,z_n,x_n), the function will return a permuted version of these vectors.

	// reserve memory
	vector<binvec> B(n,0);
	unsigned N = n*(n-1)/2;
	vector<binvec> M = coord_matrix(n);

	// fill B 
	for(unsigned i=0; i<n; i++) {
		 // since A is symmetric, row vectors = column vectors, but we have to shift the bits by n because of the lower identity block
		B[i] = A[i] << n;
	}
	
	for(unsigned j=0; j<n; j++) {
		// identity matrix
		B[j] = set_bit(B[j], n+j, 2*n); 

		// change to product coordinates
		B[j] = matrix_vector_prod_mod2(M, B[j]);
	}

	return B;
}


// --- stabilisers 

bool in_Lagrangian(const binvec &x, const vector<binvec> &B, const unsigned n) {
	unsigned test = 0;
	for(binvec v : B) {
		test += symplectic_form(x,v,n);
	} 
	return (test == 0);
}

vector<double> stabiliser_state(const vector<binvec> &B, const unsigned s) {
	// Generates the stabiliser state rho(B,s) that corresponds to the Lagrangian with basis B and a choice of signs (-1)^{s_i} where the bits s_i \in \Z_2 are specified in the length-n-bitstring s.
	// The stabiliser state is given in the Pauli basis {W(a)|a\in\Z_2^{2n}} where the only non-vanishing components rho_a corresponds to points a in the Lagrangian subspace and can thus be written in the basis B = (b_1,...,b_n) as
	//		a = \sum_{i=1}^n a_i b_i.
	// These components are explicitly given as
	//		rho_a = (-1)^{\sum_{i=1}^n a_i s_i} i^{\phi(a_1 b_1,...,a_n b_n)},
	// where \phi is the phase function appearing in the composition law of Pauli operators (see above).

	unsigned n = B.size();
	vector<double> state (pow(4,n),0);
	vector<binvec> avec (n);
	unsigned phase, idx;

	// loop over all points a in the Lagrangian
	for(binvec a=0; a<pow(2,n); a++) {
		idx = 0;

		for(unsigned i=0; i<n; i++) {
			avec.at(i) = get_bit(a,i,n) * B.at(i);
			idx ^= avec.at(i); // this gives in the end the coordinates of a w.r.t. the canonical basis
		}

		// compute phi
		phase = mod(phi(avec,n),4);
		assert(phase==0 || phase==2);
		if(phase == 2) {
			phase = 1;
		}
		else {
			phase = 0;
		}

		phase += parity(a&s); // inner product of a and s

		state.at(idx) = pow(-1,phase);// / pow(2,n);
	}

	return state;
}


// ---- I/O of the generated states

int write_states(const vector<vector<int>> &states, string states_file) {
	// write the states
	fstream fout(states_file, ios::out);
	unsigned nstates = states.size();

	if(fout.is_open()) {
		for(unsigned i=0; i<nstates; i++) {
			for(auto c : states.at(i)) {
				fout << c << " ";
			}
			fout << endl;
		}	
	}
	else {
		cout << "Error in write_states : Couldn't open file " + states_file + " for writing" << endl;
		return 1;
	}

	return 0;
}

int write_labels(const vector<string> &labels, string labels_file) {
	fstream fout(labels_file, ios::out);
	if(fout.is_open()) {
		for(auto s : labels) {
			fout << s << endl;
		}	
	}
	else {
		cout << "Error in write_labels : Couldn't open file " + labels_file + " for writing" << endl;
		return 1;
	}

	return 0;
}

int get_states(vector<vector<int>> &states, string states_file) {
	unsigned nstates = get_number_of_lines(states_file);
	
	fstream fin(states_file, ios::in);
	string line;

	if(fin.is_open()) {
		states.resize(nstates);
		unsigned i = 0;
		while(getline(fin,line)) {
			istringstream sl (line);
			string buf;

			while(sl >> buf) {
				states.at(i).push_back(stoi(buf));
			}

			++i;
		}
	}
	else {
		cout << "Error in get_states : Couldn't open states file " + states_file + " for reading" << endl;
		return 1;
	}

	return 0;
}

int get_labels(vector<string> &labels, string labels_file) {
	unsigned nlabels = get_number_of_lines(labels_file);

	fstream fin(labels_file, ios::in);

	if(fin.is_open()) {
		labels.reserve(nlabels);
		string buf;
		while(getline(fin,buf)) {
			labels.push_back(buf);
		}
	}
	else {
		cout << "Error in get_labels : Couldn't open labels file " + labels_file + " for reading" << endl;
		return 1;
	}	

	return 0;
}

int write_states(const vector<LabelledState> &states, string states_file, string labels_file="") {
	// write the states
	fstream fout(states_file, ios::out);
	unsigned nstates = states.size();

	if(fout.is_open()) {
		for(unsigned i=0; i<nstates; i++) {
			for(auto c : states.at(i).point) {
				fout << c << " ";
			}
			fout << endl;
		}	
		fout.close();
	}
	else {
		cout << "Error in write_states : Couldn't open file " + states_file + " for writing" << endl;
		return 1;
	}

	if(labels_file != "") {
		fout.open(labels_file, ios::out);
		if(fout.is_open()) {
			for(auto s : states) {
				fout << s.label << endl;
			}	
		}
		else {
			cout << "Error in write_labels : Couldn't open file " + labels_file + " for writing" << endl;
			return 2;
		}
	}

	return 0;
}

int get_states(vector<LabelledState> &states, string states_file, string labels_file="") {
	unsigned nstates = get_number_of_lines(states_file);
	
	fstream fin(states_file, ios::in);
	string line;

	if(fin.is_open()) {
		states.resize(nstates);
		unsigned i = 0;
		while(getline(fin,line)) {
			istringstream sl (line);
			string buf;

			while(sl >> buf) {
				states.at(i).point.push_back(stoi(buf));
			}

			++i;
		}
		fin.close();
	}
	else {
		cout << "Error in get_states : Couldn't open states file " + states_file + " for reading" << endl;
		return 1;
	}

	if(labels_file != "") {

		fin.open(labels_file, ios::in);

		if(fin.is_open()) {
			string buf;
			unsigned i = 0;
			while(getline(fin,buf) && i < nstates) {
				states.at(i).label = buf;
				++i;
			}
		}
		else {
			cout << "Error in get_labels : Couldn't open labels file " + labels_file + " for reading" << endl;
			return 2;
		}	
	}

	return 0;
}

// ---- projections


vector<double> project_state(const vector<double> &state, const unsigned n) {
	// size
	vector<double> pr_state(n, 0.); // 0th component is ommited since it is forced to be 1 

	// project
	for(unsigned k=1; k<=n; k++) {
		// average over all components that correspond to k X or Y operators in total
		// here we loop over the number of Y's
		for(unsigned kk=0; kk<=k; kk++) {
			vector<unsigned> aset = index_set2(kk,k,n);

			for(auto a : aset) {
				pr_state.at(k-1) += state.at(a);
			}
		}
		// normalise
		// pr_state.at(k-1) /= pow(sqrt(2),k);
	}

	return pr_state;
}

vector<double> project_state(const vector<binvec> &B, const binvec s) {
	// --- variables
	unsigned n = B.size();
	vector<double> pr_state (n,0.); // 0th component is ommited since it is forced to be 1 
	vector<binvec> avec (n,0);
	vector<unsigned> weights (4,0);
	unsigned phase, k;
	binvec aa;

	// --- loop over all points a in the Lagrangian (except 0)
	for(binvec a=1; a<pow(2,n); a++) {
		aa = 0;

		// ----- determine the Fock component the Pauli component contributes to 

		// get coordinates w.r.t. the canonical basis
		for(unsigned i=0; i<n; i++) {
			avec.at(i) = get_bit(a,i,n) * B.at(i);
			aa ^= avec.at(i); 
		}

		// Compute the weights. Only points with zero Z weight have to be considered
		weights = count_weights(aa,n);
		if(weights.at(2) == 0) {
			// determine the Fock component to which the Pauli component conributes to
			k = weights.at(1) + weights.at(3);
			if(k > 0) {
				// compute Pauli component of the state corresponding to a, this is just given by a phase

				// compute the "Lagrangian" phase phi which is i^phi with phi=0,2,4,...
				// we will write it as i^phi = (-1)^(phi/2)
				phase = phi(avec,n); 
				phase /= 2; 

				// explicit stabiliser phase which is inner product of a and s 
				phase += parity(a&s); 

				// add it to the k-th component
				pr_state.at(k-1) += pow(-1,phase);
			}
		}

	}

	return pr_state;
}


struct proj_helper {
	binvec a;
	char fock_index;
	char pauli_component;

	proj_helper(binvec aa, char f, char p) : a(aa), fock_index(f), pauli_component(p) {
	}
};


// This tries to minimise the number of Hadmards
int pr_stabiliser_from_graphs(const string file, vector<vector<int>> &states, vector<string> &labels) {

	// open the file
	fstream fin(file, ios::in);

	vector<string> graphs;
	string row;

	if(fin.is_open()) {
		while(fin >> row) {
			graphs.push_back(row);
		}
	} 
	else {
		cout << "Couldn't open file " << file << endl;
		return 1;
	}
	fin.close();

	// get number of qubits
	unsigned n = get_order( graphs.at(0) );
	unsigned graph_len = graphs.size();

	// outer loop variables
	unsigned N = pow(2,n);
	states.reserve(pow(3,n)); // estimated size
	vector<binvec> B (n);
	vector<binvec> SB (n);

	// identifies of the projected states. labels the orbits
	labels.reserve(pow(3,n));

	// inner Lagrangian loop variables
	vector<int> pr_state (n,0.); // 0th component is ommited since it is forced to be 1 
	vector<proj_helper> Lagrangian_data;
	Lagrangian_data.reserve(N/2);
	vector<binvec> avec (n,0);
	vector<unsigned> weights (4,0);
	int phase;
	binvec aa;

	// representatives of the left cosets Sp(2,Z_2)/S (where S is the symplectic representation of the phase gate, S = (11,01))
	// note that it is sufficent to only consider Hadamard gates
	vector<vector<binvec>> cosets = { vector<binvec>({0b10,0b01}), vector<binvec>({0b01,0b10}) };

	vector<vector<binvec>> LC_cosets (N);
	vector<vector<binvec>> Slist (n);
	// Changed convention here: the bits of i encode wether a Hadamard is applied to the qubits; if the j-th bit from the right is 1, apply a Hadamard to the j-th qubit
	// Formerly: j-th bit from the left
	// This simplifies the reduced Hadamard loop, where Hadamards are only applied to the first r qubits
	for(binvec i=0; i<N; i++) {
		for(unsigned j=0; j<n; j++) {
			Slist.at(j) = cosets.at(get_bit2(i, j));
		}
		LC_cosets.at(i) = direct_sum(Slist);
	}
	// Note: This whole matrix thing can be replaced with bit manipulation

	// graph state loop
	int dist = 0;
	binvec b = 0;
	vector<binvec> A;

	cout << "  Compute projected stabiliser states for graph" << endl;

	// here, we take the Hadamard loop as outer loop to minimise the number of Hadamards.

	for(binvec h=0; h<N; h++) {

		// cout << "\r    #" << g << " / " << graph_len << flush;

		for(unsigned g=0; g<graph_len; g++) {
			
			A = graph6_to_adj_mat( graphs.at(g) );
			B = graph_Lagrangian(A,n);

			// we have to take the local Hadamard orbit of that graph

			for(unsigned k=0; k<n; k++) {
				SB.at(k) = matrix_vector_prod_mod2(LC_cosets.at(h), B.at(k));
			}		

			// loop over all points a in the Lagrangian SB (except 0) and fill the vector Lagrangian_data

			Lagrangian_data.clear();
			for(binvec a=1; a<N; a++) {
				aa = 0;

				// get coordinates w.r.t. the canonical basis
				for(unsigned i=0; i<n; i++) {
					avec.at(i) = get_bit(a,i,n) * SB.at(i);
					aa ^= avec.at(i); 
				}

				// Compute the weights. Only points with zero Z weight have to be considered
				weights = count_weights(aa,n);

				if(weights.at(2) == 0) {
					// compute Pauli component of the state corresponding to a, this is just given by a phase

					// compute the "Lagrangian" phase phi which is i^phi with phi=0,2,4,...
					// we will write it as i^phi = (-1)^(phi/2)
					phase = phi(avec,n); 
					phase /= 2; 

					// save for later
					Lagrangian_data.push_back(proj_helper(a, (char)(weights.at(1)+weights.at(3)), (char)pow(-1,phase)));				
				}
			}

			// --- now, loop over phases and use the Lagrangian data computed before
			for(binvec s=0; s<N; s++) {
				// clear pr_state
				pr_state.assign(n,0);

				// note that this loop will only involve non-trivial components
				for(binvec i=0; i<Lagrangian_data.size(); i++) {

					// explicit stabiliser phase which is inner product of a and s 
					phase = parity(Lagrangian_data.at(i).a & s); 

					// add it to the right Fock component
					pr_state.at(Lagrangian_data.at(i).fock_index - 1) += ( pow(-1,phase) * Lagrangian_data.at(i).pauli_component );		
				}

				// check if projected state pr_state already exists
				// to do that, we use binary search with std::lower_bound() which gives an iterator on the first element that is not less than pr_state
				auto it = lower_bound(states.begin(), states.end(), pr_state);
				if(it == states.end() || pr_state < *it) {
					// element was not redundant, so also add it to the list
					dist = distance(states.begin(),it);
					states.insert(it, pr_state);

					// add a label
					labels.insert(labels.begin()+dist, graphs.at(g)+" "+write_bits(invert_bits(h,n),n)+" 0 "+write_bits(s,n));
				}
				// note that procedures preserves the ordering in states ... 
			}
		}
	}

	cout << endl;

	return 0;
}

// ---- routines that rely on connectedness

// This assumes that the input are representatives of isomorphism and local complementation orbits of graphs. To get the full set of projected stabilisers, we have the following orbits of these representatives:
//	1. W.r.t. to the cosets of Sp(2,Z_2) / <S>  (where S is the symplectic representation of the phase gate, S = (11,01))
//  2. W.r.t. to the signs of the n generators
//
// The cosets are represented by {I,H,HS} and S can be effectively implemented by switching the representative diagonal entry of the adjacency matrix of the graph from 0 to 1.
int pr_stabiliser_from_LCgraphs(const string file, vector<vector<int>> &states, vector<string> &labels) {
	
	// open the file
	fstream fin(file, ios::in);

	vector<string> graphs;
	string row;

	if(fin.is_open()) {
		while(fin >> row) {
			graphs.push_back(row);
		}
	} 
	else {
		cout << "Couldn't open file " << file << endl;
		return 1;
	}

	fin.close();

	// get number of qubits
	unsigned n = get_order( graphs.at(0) );
	unsigned graph_len = graphs.size();

	// outer loop variables
	unsigned N = pow(2,n);
	unsigned M = pow(3,n);
	states.reserve(M); // estimated size
	vector<binvec> B (n);
	vector<binvec> SB (n);

	// identifies of the projected states. labels the orbits
	labels.reserve(M);

	// inner Lagrangian loop variables
	vector<int> pr_state (n,0); // 0th component is ommited since it is forced to be 1 
	vector<proj_helper> Lagrangian_data;
	Lagrangian_data.reserve(N/2);
	vector<binvec> avec (n,0);
	vector<unsigned> weights (4,0);
	int phase;
	binvec aa;

	// representatives of the left cosets Sp(2,Z_2)/S (where S is the symplectic representation of the phase gate, S = (11,01))
	vector<vector<binvec>> cosets = { vector<binvec>({0b10,0b01}), vector<binvec>({0b01,0b10}), vector<binvec>({0b01,0b11}) };

	vector<vector<binvec>> LC_cosets (M);
	vector<vector<binvec>> Slist (n);

	vector<unsigned> trits(n, 0);
	for(binvec i=0; i<M; i++) {
		get_multi_index(n, 3, i, trits);		
		for(unsigned j=0; j<n; j++) {
			Slist.at(j) = cosets.at(trits.at(j));
		}
		LC_cosets.at(i) = direct_sum(Slist);
	}

	// graph state loop
	int dist = 0;
	binvec b = 0;

	cout << "  Compute projected stabiliser states for graph" << endl;

	for(unsigned g=0; g<graph_len; g++) {

		cout << "\r    #" << g+1 << " / " << graph_len << flush;
		B = graph_Lagrangian( graph6_to_adj_mat( graphs.at(g) ),n);

		// we have to take the local coset orbit of that Lagrangian
		for(binvec h=0; h<M; h++) {

			for(unsigned k=0; k<n; k++) {
				SB.at(k) = matrix_vector_prod_mod2(LC_cosets.at(h), B.at(k));
			}

			// loop over all points a in the Lagrangian SB (except 0) and fill the vector Lagrangian_data

			Lagrangian_data.clear();
			for(binvec a=1; a<N; a++) {
				aa = 0;

				// get coordinates w.r.t. the canonical basis
				for(unsigned i=0; i<n; i++) {
					avec.at(i) = get_bit(a,i,n) * SB.at(i);
					aa ^= avec.at(i); 
				}

				// Compute the weights. Only points with zero Z weight have to be considered
				weights = count_weights(aa,n);

				if(weights.at(2) == 0) {
					// compute Pauli component of the state corresponding to a, this is just given by a phase

					// compute the "Lagrangian" phase phi which is i^phi with phi=0,2,4,...
					// we will write it as i^phi = (-1)^(phi/2)
					phase = phi(avec,n); 
					phase /= 2; 

					// save for later
					Lagrangian_data.push_back(proj_helper(a, (char)(weights.at(1)+weights.at(3)), (char)pow(-1,phase)));				
				}
			}

			// --- now, loop over phases and use the Lagrangian data computed before
			for(binvec s=0; s<N; s++) {
				// clear pr_state
				pr_state.assign(n,0);

				// note that this loop will only involve non-trivial components
				for(binvec i=0; i<Lagrangian_data.size(); i++) {

					// explicit stabiliser phase which is inner product of a and s 
					phase = parity(Lagrangian_data.at(i).a & s); 

					// add it to the right Fock component
					// cout << (unsigned)Lagrangian_data.at(i).fock_index << endl;
					pr_state.at(Lagrangian_data.at(i).fock_index - 1) += ( pow(-1,phase) * Lagrangian_data.at(i).pauli_component );		
				}

				// check if projected state pr_state already exists
				// to do that, we use binary search with std::lower_bound() which gives an iterator on the first element that is not less than pr_state
				auto it = lower_bound(states.begin(), states.end(), pr_state);
				if(it == states.end() || pr_state < *it) {
					// element was not redundant, so also add it to the list
					dist = distance(states.begin(),it);
					states.insert(it, pr_state);

					// add a label
					labels.insert(labels.begin()+dist, graphs.at(g)+" "+write_trits(h,n)+" "+write_bits(s,n));
				}
				// note that procedures preserves the ordering in states ... 
			}
		}
	}

	cout << endl;

	return 0;
}

int pr_stabiliser_from_LCgraphs(const string file, vector<LabelledState> &states) {
	
	// open the file
	fstream fin(file, ios::in);

	vector<string> graphs;
	string row;

	if(fin.is_open()) {
		while(fin >> row) {
			graphs.push_back(row);
		}
	} 
	else {
		cout << "Couldn't open file " << file << endl;
		return 1;
	}

	fin.close();

	// get number of qubits
	unsigned n = get_order( graphs.at(0) );
	unsigned graph_len = graphs.size();

	// outer loop variables
	unsigned N = pow(2,n);
	unsigned M = pow(3,n);
	states.reserve(M); // estimated size
	vector<binvec> B (n);
	vector<binvec> SB (n);

	// identifies of the projected states. labels the orbits
	// labels.reserve(M);

	// inner Lagrangian loop variables
	LabelledState pr_state (n); // 0th component is ommited since it is forced to be 1 
	vector<proj_helper> Lagrangian_data;
	Lagrangian_data.reserve(N/2);
	vector<binvec> avec (n,0);
	vector<unsigned> weights (4,0);
	int phase;
	binvec aa;

	// representatives of the left cosets Sp(2,Z_2)/S (where S is the symplectic representation of the phase gate, S = (11,01))
	vector<vector<binvec>> cosets = { vector<binvec>({0b10,0b01}), vector<binvec>({0b01,0b10}), vector<binvec>({0b01,0b11}) };

	vector<vector<binvec>> LC_cosets (M);
	vector<vector<binvec>> Slist (n);

	vector<unsigned> trits(n, 0);
	for(binvec i=0; i<M; i++) {
		get_multi_index(n, 3, i, trits);		
		for(unsigned j=0; j<n; j++) {
			Slist.at(j) = cosets.at(trits.at(j));
		}
		LC_cosets.at(i) = direct_sum(Slist);
	}

	// graph state loop
	cout << "  Compute projected stabiliser states for graph" << endl;

	for(unsigned g=0; g<graph_len; g++) {

		cout << "\r    #" << g+1 << " / " << graph_len << flush;
		B = graph_Lagrangian( graph6_to_adj_mat( graphs.at(g) ),n);

		// we have to take the local coset orbit of that Lagrangian
		for(binvec h=0; h<M; h++) {

			for(unsigned k=0; k<n; k++) {
				SB.at(k) = matrix_vector_prod_mod2(LC_cosets.at(h), B.at(k));
			}

			// loop over all points a in the Lagrangian SB (except 0) and fill the vector Lagrangian_data

			Lagrangian_data.clear();
			for(binvec a=1; a<N; a++) {
				aa = 0;

				// get coordinates w.r.t. the canonical basis
				for(unsigned i=0; i<n; i++) {
					avec.at(i) = get_bit(a,i,n) * SB.at(i);
					aa ^= avec.at(i); 
				}

				// Compute the weights. Only points with zero Z weight have to be considered
				weights = count_weights(aa,n);

				if(weights.at(2) == 0) {
					// compute Pauli component of the state corresponding to a, this is just given by a phase

					// compute the "Lagrangian" phase phi which is i^phi with phi=0,2,4,...
					// we will write it as i^phi = (-1)^(phi/2)
					phase = phi(avec,n); 
					phase /= 2; 

					// save for later
					Lagrangian_data.push_back(proj_helper(a, (char)(weights.at(1)+weights.at(3)), (char)pow(-1,phase)));				
				}
			}

			// --- now, loop over phases and use the Lagrangian data computed before
			for(binvec s=0; s<N; s++) {
				// clear pr_state
				pr_state.clear();

				// note that this loop will only involve non-trivial components
				for(binvec i=0; i<Lagrangian_data.size(); i++) {

					// explicit stabiliser phase which is inner product of a and s 
					phase = parity(Lagrangian_data.at(i).a & s); 

					// add it to the right Fock component
					// cout << (unsigned)Lagrangian_data.at(i).fock_index << endl;
					pr_state.point.at(Lagrangian_data.at(i).fock_index - 1) += ( pow(-1,phase) * Lagrangian_data.at(i).pauli_component );		
				}

				// add a label
				pr_state.label = graphs.at(g)+" "+write_trits(h,n)+" "+write_bits(s,n);

				// check if projected state pr_state already exists
				// to do that, we use binary search with std::lower_bound() which gives an iterator on the first element that is not less than pr_state
				// note that LabelledState has an overloaded < operator
				auto it = lower_bound(states.begin(), states.end(), pr_state);

				if(it == states.end() || pr_state < *it) {
					// element was not redundant, so add it to the list
					states.insert(it, pr_state);
				}
				// note that procedures preserves the ordering in states ... 
			}
		}
	}

	cout << endl;

	states.shrink_to_fit();

	return 0;
}


// ---- computes the projections of product states


// the coordinates of the projection of product state \rho \otimes \sigma can be expressed by their individual projections only. If A_k(.) denotes the k-th coordinate, i.e. the k-th signed XY weight enumerator, then
// A_k( \rho\otimes\sigma ) = \sum_{i=0}^k A_i(\rho) A_{k-i}(\sigma).
vector<int> pr_product_state(vector<int> rho, vector<int> sigma) {
	unsigned n = rho.size();
	unsigned m = sigma.size();

	vector<int> ret (n+m,0);
	vector<int> r (n+m+1,0);
	vector<int> s (n+m+1,0);

	// copy rho and sigma to r and s which are padded with zeros, to avoid headaches and bugs ;)
	copy(rho.begin(),rho.end(),r.begin()+1);
	copy(sigma.begin(),sigma.end(),s.begin()+1);
	r.at(0) = 1;
	s.at(0) = 1;

	for(unsigned k=1; k<=n+m; k++) {	
		for(unsigned i=0; i<=k; i++) {
			ret.at(k-1) += r.at(i)*s.at(k-i);
		}
	}

	return ret;
}

string pr_product_label(string label1, string label2) {
	// get the components of every label 
	vector<string> c1;
	vector<string> c2;

	istringstream ss (label1);
	string buf;

	while(ss >> buf) {
		c1.push_back(buf);
	}

	ss = istringstream(label2);
	while(ss >> buf) {
		c2.push_back(buf);
	}

	// get the graph6 representation of the product graph
	vector<binvec> A1 = graph6_to_adj_mat(c1.at(0));
	vector<binvec> A2 = graph6_to_adj_mat(c2.at(0));
	vector<binvec> A = direct_sum(A1, A2);

	string ret = adj_mat_to_graph6(A);

	// concatenate the operator and sign encoding
	ret += " "+c1.at(1)+c2.at(1);
	ret += " "+c1.at(2)+c2.at(2);

	return ret;
}

LabelledState pr_product_state(LabelledState rho, LabelledState sigma) {
	unsigned n = rho.size();
	unsigned m = sigma.size();

	LabelledState ret (n+m);
	vector<int> r (n+m+1,0);
	vector<int> s (n+m+1,0);

	// copy rho and sigma to r and s which are padded with zeros, to avoid headaches and bugs ;)
	copy(rho.point.begin(),rho.point.end(),r.begin()+1);
	copy(sigma.point.begin(),sigma.point.end(),s.begin()+1);
	r.at(0) = 1;
	s.at(0) = 1;

	for(unsigned k=1; k<=n+m; k++) {	
		for(unsigned i=0; i<=k; i++) {
			ret.point.at(k-1) += r.at(i)*s.at(k-i);
		}
	}

	ret.label = pr_product_label(rho.label, sigma.label);

	return ret;
}

// filename templates are assumed to contain placeholder like %d for system size
int pr_product_states(unsigned n, vector<vector<int>> &states, string states_file_tpl, vector<string> &labels, string labels_file_tpl="") {

	//  get partitions of n
	auto partitions = get_partitions(n);

	// clear containers
	states.clear();
	labels.clear();

	// read vertices 
	vector<vector<vector<int>>> vstates (n-1);
	vector<vector<string>> vlabels (n-1);
	vector<unsigned> nvertices (n-1);

	for(unsigned i=0; i<n-1; i++) {
		char *fstate = new char [2*states_file_tpl.size()];
		char *flabel = new char [2*labels_file_tpl.size()];

		snprintf(fstate, 2*states_file_tpl.size(), states_file_tpl.c_str(), i+1);

		if(get_states(vstates.at(i), string(fstate)) != 0) {
			return 1;
		}

		if(labels_file_tpl != "") {
			snprintf(flabel, 2*labels_file_tpl.size(), labels_file_tpl.c_str(), i+1);

			if(get_labels(vlabels.at(i), string(flabel)) != 0) {
				return 1;
			}

		}

		nvertices.at(i) = vstates.at(i).size();

		delete[] fstate;
		delete[] flabel;
	}

	for(auto p = partitions.begin()+1; p != partitions.end(); ++p) {
		
		unsigned L = p->size();
		unsigned N = 1;

		vector<unsigned> nv (L);

		// get total number of product states corresponding to that partition
		for(unsigned i=0; i<L; i++) {
			N *= nvertices.at( p->at(i)-1 );
			nv.at(i) = nvertices.at( p->at(i)-1 );
		}
		
		// build up all product states corresponding to tuples (i_1,...,i_L)
		vector<unsigned> ind (L);
		vector<int> state;
		string label;

		unsigned dist;

		for(unsigned i=0; i<N; i++) {
			get_multi_index(nv, i, ind);

			state = vstates.at( p->at(0)-1 ).at( ind.at(0) );

			for(unsigned k=1; k<L; k++) {
				state = pr_product_state( state, vstates.at( p->at(k)-1 ).at( ind.at(k) ) );
			}

			// check if state already exists in the list
			// to do that, we use binary search with std::lower_bound() which gives an iterator on the first element that is not less than state
			auto it = lower_bound(states.begin(), states.end(), state);
			if(it == states.end() || state < *it) {
				// element was not redundant, so also add it to the list
				dist = distance(states.begin(),it);
				states.insert(it, state);

				// add a label
				label = vlabels.at( p->at(0)-1 ).at( ind.at(0) );
				for(unsigned k=1; k<L; k++) {
					label = pr_product_label( label, vlabels.at( p->at(k)-1 ).at( ind.at(k) ) );
				}
				
				labels.insert(labels.begin()+dist, label);
			}
		}
	}

	return 0;
}

int pr_product_states(unsigned n, vector<LabelledState> &states, string states_file_tpl, string labels_file_tpl="") {

	//  get partitions of n
	auto partitions = get_partitions(n);

	// clear containers
	states.clear();

	// read vertices 
	vector<vector<LabelledState>> vertices (n-1);
	vector<unsigned> nvertices (n-1);

	for(unsigned i=0; i<n-1; i++) {
		char *fstate = new char [2*states_file_tpl.size()];
		char *flabel = new char [2*labels_file_tpl.size()];

		snprintf(fstate, 2*states_file_tpl.size(), states_file_tpl.c_str(), i+1);

		if(get_states(vertices.at(i), string(fstate)) != 0) {
			return 1;
		}

		if(labels_file_tpl != "") {
			snprintf(flabel, 2*labels_file_tpl.size(), labels_file_tpl.c_str(), i+1);
			snprintf(fstate, 2*states_file_tpl.size(), states_file_tpl.c_str(), i+1);

			if(get_states(vertices.at(i), string(fstate), string(flabel)) != 0) {
				return 1;
			}
		}
		else {
			snprintf(fstate, 2*states_file_tpl.size(), states_file_tpl.c_str(), i+1);

			if(get_states(vertices.at(i), string(fstate)) != 0) {
				return 1;
			}
		}

		nvertices.at(i) = vertices.at(i).size();

		delete[] fstate;
		delete[] flabel;
	}

	LabelledState pr_state;

	for(auto p = partitions.begin()+1; p != partitions.end(); ++p) {
		
		unsigned L = p->size();
		unsigned N = 1;

		vector<unsigned> nv (L);

		// get total number of product states corresponding to that partition
		for(unsigned i=0; i<L; i++) {
			N *= nvertices.at( p->at(i)-1 );
			nv.at(i) = nvertices.at( p->at(i)-1 );
		}
		
		// build up all product states corresponding to tuples (i_1,...,i_L)
		vector<unsigned> ind (L);

		for(unsigned i=0; i<N; i++) {
			get_multi_index(nv, i, ind);

			pr_state = vertices.at( p->at(0)-1 ).at( ind.at(0) );

			for(unsigned k=1; k<L; k++) {
				pr_state = pr_product_state( pr_state, vertices.at( p->at(k)-1 ).at( ind.at(k) ) );
			}

			// check if state already exists in the list
			// to do that, we use binary search with std::lower_bound() which gives an iterator on the first element that is not less than state
			auto it = lower_bound(states.begin(), states.end(), pr_state);
			if(it == states.end() || pr_state < *it) {
				// element was not redundant, so also add it to the list
				states.insert(it, pr_state);
			}
		}
	}

	return 0;
}


#endif