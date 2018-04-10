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

#ifndef SYMPLECTIC_H
#include "symplectic.h"
#endif

#ifndef GLPKCONVEXSEPARATION_H
#include "GLPKConvexSeparation.h"
#endif


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

// --- graph states

// read graph6 format and returns upper triangle of adjacency matrix in vectorised form
// IMPORTANT: this implementation is quick and dirty and hence only works for graphs with less than 63 vertices (fine for me)
binvec graph6_to_adj_mat(const string g6) {
	// The input string g6 is assumed to hold a graph6 representation of the graph in ASCII-encoded form, i.e. it's a sequence of ASCII characters

	// interprete g6 as integer
	binvec x = 0;
	unsigned n = (unsigned)g6.at(0)-63;
	unsigned l = g6.size();
	unsigned N = n*(n-1)/2;
	unsigned i,j;

	for(i=1; i<l; i++) {
		x ^= ((binvec)g6.at(i)-63) << (6*(l-i-1));
	}
	x = get_bits(x, 0, N-1, 6*(l-1));

	// the just created vector has column-wise ordering, convert it to row-wise ordering
	binvec ret = 0;
	for(unsigned k=0; k<N; k++) {
		if(get_bit(x,k,N) == 1) {
			j = get_ut_col2(k,n);
			i = get_ut_row2(j,k,n);
			ret = set_bit(ret, get_ut_index(i,j,n), N);
		}
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

// reads the geng output of adjacency matrices for all non-isomorphic graphs of n vertices
vector<binvec> read_graphs(const string file, const unsigned n) {
	// the data format is assumed to store the upper triangle of the adjacency matrices of the graphs as e.g.
	//    000
	//    00
	//    0
	// the matrices are proceded by some header line
	// 
	// The matrices are stored row-wise as a bitstring and returned as a vector<binvec_short>

	// open the file
	fstream fin(file, ios::in);

	vector<string> data (n-1);
	string row;
	unsigned counter=0;
	unsigned idx;

	vector<binvec> ret;

	if(fin.is_open()) {
		while(fin >> row) {
			idx = counter % n;
			if(idx != 0) {
				data.at(idx-1) = row;
			} 
			if(idx == n-1) {
				// end of block

				// concatenate all strings
				string tot;
				for(auto str : data) {
					tot += str;
				}
				// interprete as bitstring and convert to unsigned 
				ret.push_back( stoul(tot,0,2) );
			}
			++counter;
		}
	} 
	else {
		cout << "Couldn't open file " << file << endl;
	}

	fin.close();

	return ret;
}

vector<binvec> generate_gs_Lagrangian(const binvec &A, const unsigned n) {
	// This returns the representation of a graph state given in terms of a basis for the Lagrangian subspace. The graph of the graph state is specified by the upper triangle of its adjacency matrix, given as a binary vector A in Z_2^{n(n-1)/2}.
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
	unsigned i,j;

	for(unsigned k=0; k<N; k++) {
		if( get_bit(A, k, N) == 1) {
			i = get_ut_row(k,n);
			j = get_ut_col(i,k,n);

			B[j] = set_bit(B[j], i, 2*n);
			B[i] = set_bit(B[i], j, 2*n);
		}
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

vector<double> generate_stabiliser_state(const vector<binvec> &B, const unsigned s) {
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
	for(binvec_short a=0; a<pow(2,n); a++) {
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


// --- projection


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

vector<double> project_state(const vector<binvec> &B, const binvec_short s) {
	// --- variables
	unsigned n = B.size();
	vector<double> pr_state (n,0.); // 0th component is ommited since it is forced to be 1 
	vector<binvec> avec (n,0);
	vector<unsigned> weights (4,0);
	unsigned phase, k;
	binvec aa;

	// --- loop over all points a in the Lagrangian (except 0)
	for(binvec_short a=1; a<pow(2,n); a++) {
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
	binvec_short a;
	char fock_index;
	char pauli_component;

	proj_helper(binvec_short aa, char f, char p) : a(aa), fock_index(f), pauli_component(p) {
	}
};

pair< vector<vector<int>>, vector<string> > generate_projected_stabiliser_states_from_graphs(const string file) {

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
	}

	fin.close();

	// get number of qubits
	unsigned n = get_order( graphs.at(0) );
	unsigned graph_len = graphs.size();

	// outer loop variables
	unsigned N = pow(2,n);
	vector<vector<int>> states;
	states.reserve(pow(3,n)); // estimated size
	vector<binvec> B (n);
	vector<binvec> SB (n);

	// identifies of the projected states. labels the orbits
	vector<string> labels; 
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
	for(binvec i=0; i<N; i++) {
		for(unsigned j=0; j<n; j++) {
			Slist.at(j) = cosets.at(get_bit(i, j, n));
		}
		LC_cosets.at(i) = direct_sum(Slist);
	}

	// graph state loop
	int dist = 0;
	binvec b = 0;

	cout << "  Compute projected stabiliser states for graph" << endl;

	for(unsigned g=0; g<graph_len; g++) {
		cout << "\r    #" << g << " / " << graph_len << flush;
		B = generate_gs_Lagrangian( graph6_to_adj_mat( graphs.at(g) ),n);

		// we have to take the local Hadamard orbit of that graph, but only on nodes/qubits which form a maximal disconnected subgraph.
		for(binvec h=0; h<N; h++) {

			for(unsigned k=0; k<n; k++) {
				SB.at(k) = matrix_vector_prod_mod2(LC_cosets.at(h), B.at(k));
			}		

			// loop over all points a in the Lagrangian SB (except 0) and fill the vector Lagrangian_data

			Lagrangian_data.clear();
			for(binvec_short a=1; a<pow(2,n); a++) {
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
			for(binvec_short s=0; s<N; s++) {
				// clear pr_state
				pr_state.assign(n,0);

				// note that this loop will only involve non-trivial components
				for(binvec_short i=0; i<Lagrangian_data.size(); i++) {

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
					labels.insert(labels.begin()+dist, graphs.at(g)+" "+write_bits(h,n)+" 0 "+write_bits(s,n));
				}
				// note that procedures preserves the ordering in states ... 
			}
		}
	}

	cout << endl;

	return make_pair(states,labels);
}

pair< vector<vector<int>>, vector<string> > generate_projected_stabiliser_states_from_graphs_wloops(const string file) {

	
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
	}

	fin.close();

	// get number of qubits
	unsigned n = get_order( graphs.at(0) );
	unsigned graph_len = graphs.size();

	// outer loop variables
	unsigned N = pow(2,n);
	vector<vector<int>> states;
	states.reserve(pow(3,n)); // estimated size
	vector<binvec> B (n);
	vector<binvec> SB (n);

	// identifies of the projected states. labels the orbits
	vector<string> labels; 
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
	// note that HS = (01,11) is missing since it is equivalent to adding all possible loops to the graph + acting with H
	vector<vector<binvec>> cosets = { vector<binvec>({0b10,0b01}), vector<binvec>({0b01,0b10}) };

	vector<vector<binvec>> LC_cosets (N);
	vector<vector<binvec>> Slist (n);
	for(binvec i=0; i<N; i++) {
		for(unsigned j=0; j<n; j++) {
			Slist.at(j) = cosets.at(get_bit(i, j, n));
		}
		LC_cosets.at(i) = direct_sum(Slist);
	}

	// graph state loop
	int dist = 0;
	binvec b = 0;

	cout << "  Compute projected stabiliser states for graph" << endl;

	for(unsigned g=0; g<graph_len; g++) {
		cout << "\r    #" << g << " / " << graph_len << flush;
		B = generate_gs_Lagrangian( graph6_to_adj_mat( graphs.at(g) ),n);

		// we have to take the local Hadamard orbit of that Lagrangian + loops
		for(binvec h=0; h<N; h++) {
			// loops loop
			for(char l=0; l<=1; l++) {

				if(l==0) {
					for(unsigned k=0; k<n; k++) {
						SB.at(k) = matrix_vector_prod_mod2(LC_cosets.at(h), B.at(k));
					}
				}
				else {
					for(unsigned k=0; k<n; k++) {
						b = B.at(k);
						if(get_bit(h,k,n) == 1) {
							// change the k-th basis vector / stabiliser
							b = set_bit(b, 2*k, 2*n);
						}
						SB.at(k) = matrix_vector_prod_mod2(LC_cosets.at(h), b);
					}
				}

				// loop over all points a in the Lagrangian SB (except 0) and fill the vector Lagrangian_data

				Lagrangian_data.clear();
				for(binvec_short a=1; a<pow(2,n); a++) {
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
				for(binvec_short s=0; s<N; s++) {
					// clear pr_state
					pr_state.assign(n,0);

					// note that this loop will only involve non-trivial components
					for(binvec_short i=0; i<Lagrangian_data.size(); i++) {

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
						labels.insert(labels.begin()+dist, graphs.at(g)+" "+write_bits(h,n)+" "+to_string(l)+" "+write_bits(s,n));
					}
					// note that procedures preserves the ordering in states ... 
				}
			}
		}
	}

	cout << endl;

	return make_pair(states,labels);
}


vector<vector<int>> generate_projected_stabiliser_states_from_graph(const string file, const int gid, const unsigned n) {

	// outer loop variables
	unsigned N = pow(2,n);
	vector<vector<int>> states;
	states.reserve(pow(3,n)); // estimated size
	vector<binvec> B (n);
	vector<binvec> SB (n);
	int _gid;

	// inner Lagrangian loop variables
	vector<int> pr_state (n,0.); // 0th component is ommited since it is forced to be 1 
	vector<proj_helper> Lagrangian_data;
	Lagrangian_data.reserve(N/2);
	vector<binvec> avec (n,0);
	vector<unsigned> weights (4,0);
	unsigned phase;
	binvec aa;

	// generate group generated by local hadmards
	vector<vector<binvec>> H_group = { vector<binvec>({0b10,0b01}), vector<binvec>({0b01,0b10}) };

	vector<vector<binvec>> LH_group (N);
	vector<vector<binvec>> Slist (n);
	for(binvec i=0; i<N; i++) {
		for(unsigned j=0; j<n; j++) {
			Slist.at(j) = H_group.at(get_bit(i, j, n));
		}
		LH_group.at(i) = direct_sum(Slist);
	}

	// graph state loop
	vector<binvec> graphs = read_graphs(file,n);

	unsigned graph_len = graphs.size();
	cout << "  Compute projected stabiliser states for maximally connected graph" << endl;

	assert(abs(gid) < graph_len);
	if(gid < 0) {
		_gid = graph_len + gid;
	} else {
		_gid = gid;
	}
	auto th = graphs.at(_gid);

	B = generate_gs_Lagrangian(th,n);

	// we have to take the local Hadamard orbit of that Lagrangian
	for(auto S : LH_group) {
		for(unsigned k=0; k<n; k++) {
			SB.at(k) = matrix_vector_prod_mod2(S, B.at(k));
		}

		// --- loop over all points a in the Lagrangian SB (except 0) and fill the vector Lagrangian_data

		Lagrangian_data.clear();
		for(binvec_short a=1; a<pow(2,n); a++) {
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
		for(binvec_short s=0; s<N; s++) {
			// clear pr_state
			pr_state.assign(n,0);

			// note that this loop will only involve non-trivial components
			for(binvec_short i=0; i<Lagrangian_data.size(); i++) {

				// explicit stabiliser phase which is inner product of a and s 
				phase = parity(Lagrangian_data.at(i).a & s); 

				// add it to the right Fock component
				pr_state.at(Lagrangian_data.at(i).fock_index - 1) += pow(-1,phase) * Lagrangian_data.at(i).pauli_component;		
			}

			// check if projected state pr_state already exists
			// to do that, we use binary search with std::lower_bound() which gives an iterator on the first element that is not less than pr_state
			auto it = lower_bound(states.begin(), states.end(), pr_state);
			if(it == states.end() || pr_state < *it) {
				// element was not redundant, so also add it to the list
				states.insert(it, pr_state);
			}
			// note that procedures preserves the ordering in states ... 
		}
	}

	cout << endl;

	return states;
}

GLPKConvexSeparation generate_projected_stabiliser_states_from_graphs_conv(const string file, const unsigned n) {

	// outer loop variables
	unsigned N = pow(2,n);
	vector<vector<int>> states;
	states.reserve(pow(3,n)); // estimated size
	vector<binvec> B (n);
	vector<binvec> SB (n);

	// we keep a GLPKConvexSeparation object to check if elements that we want to add to the list are in the convex hull of the previous ones
	GLPKConvexSeparation lp (n);
	lp.set_verbosity(1);

	// inner Lagrangian loop variables
	vector<int> pr_state (n,0.); // 0th component is ommited since it is forced to be 1 
	vector<proj_helper> Lagrangian_data;
	Lagrangian_data.reserve(N/2);
	vector<binvec> avec (n,0);
	vector<unsigned> weights (4,0);
	unsigned phase;
	binvec aa;

	// generate group generated by local hadmards
	vector<vector<binvec>> H_group = { vector<binvec>({0b10,0b01}), vector<binvec>({0b01,0b10}) };

	vector<vector<binvec>> LH_group (N);
	vector<vector<binvec>> Slist (n);
	for(binvec i=0; i<N; i++) {
		for(unsigned j=0; j<n; j++) {
			Slist.at(j) = H_group.at(get_bit(i, j, n));
		}
		LH_group.at(i) = direct_sum(Slist);
	}

	// graph state loop
	vector<binvec> graphs = read_graphs(file,n);

	unsigned counter = 1;
	unsigned graph_len = graphs.size();
	cout << "  Compute projected stabiliser states for graph" << endl;

	for(unsigned g=0; g<graph_len; g++) {
		cout << "\r    #" << counter << " / " << graph_len << flush;
		B = generate_gs_Lagrangian(graphs.at(g),n);

		// we have to take the local Hadamard orbit of that Lagrangian
		for(binvec h=0; h<N; h++) {
			for(unsigned k=0; k<n; k++) {
				SB.at(k) = matrix_vector_prod_mod2(LH_group.at(h), B.at(k));
			}

			// --- loop over all points a in the Lagrangian SB (except 0) and fill the vector Lagrangian_data

			Lagrangian_data.clear();
			for(binvec_short a=1; a<pow(2,n); a++) {
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
			for(binvec_short s=0; s<N; s++) {
				// clear pr_state
				pr_state.assign(n,0);

				// note that this loop will only involve non-trivial components
				for(binvec_short i=0; i<Lagrangian_data.size(); i++) {

					// explicit stabiliser phase which is inner product of a and s 
					phase = parity(Lagrangian_data.at(i).a & s); 

					// add it to the right Fock component
					pr_state.at(Lagrangian_data.at(i).fock_index - 1) += pow(-1,phase) * Lagrangian_data.at(i).pauli_component;		
				}

				// check if projected state pr_state already exists
				// to do that, we use binary search with std::lower_bound() which gives an iterator on the first element that is not less than pr_state
				auto it = lower_bound(states.begin(), states.end(), pr_state);
				if(it == states.end() || pr_state < *it) {
					// the element is new, check now for convex dependence
					if(lp.add_vertex(pr_state, to_string(g)+" "+write_bits(h,n)+" "+write_bits(s,n)) == 0) {
						// element was not redundant, so also add it to the list
						states.insert(it, pr_state);
					}
				}
				// note that procedures preserves the ordering in states ... 
			}
		}
		++counter;
	}

	cout << endl;

	return lp;
}


#endif