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

#ifndef SYMPLECTIC_H
#include "symplectic.h"
#endif

#ifndef CONVEXSEPARATION_H
#include "ConvexSeparation.h"
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


// reads the geng output of non-isomorphic graphs 
// the data format is assumed to store the upper triangle of the adjacency matrices of the graphs as e.g.
//    000
//    00
//    0
// the matrices are proceded by some header line
// 
// The matrices are stored row-wise as a bitstring and returned as a vector<binvec_short>
vector<binvec_short> read_graph_states(const string file, const unsigned n) {
	// open the file
	fstream fin(file, ios::in);

	vector<string> data (n-1);
	string row;
	unsigned counter=0;
	unsigned idx;

	vector<binvec_short> ret;

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

vector<binvec> generate_gs_Lagrangian(const binvec_short th, const unsigned n) {
	// This generates the i-th symmetric nxn matrix and returns the representation of the corresponding graph state given in terms of a basis for the Lagrangian subspace
	// The basis vectors are given by the column vectors of the matrix
	// 			(  A  )
	// 		B = ( --- )			(in (z_1,...,z_n,x_1,...,x_n) coordinates)
	//			(  I  )
	// where A is a nxn symmetric matrix and I is the nxn identity matrix. Note that since we are working in product coordinates (z_1,x_1,z_2,x_2,...,z_n,x_n), the function will return a permuted version of these vectors.
	// The i-th vector from Z_2^{n(n+1)/2} is simply given by the binary representation of i

	// reserve memory
	vector<binvec> B(n,0);
	unsigned N = n*(n-1)/2;
	vector<binvec> M = coord_matrix(n);

	// fill B 
	unsigned i,j;

	for(unsigned k=0; k<N; k++) {
		if( get_bit(th, k, N) == 1) {
			i = get_symmetric_row(k,n);
			j = get_symmetric_col(i,k,n);

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

bool in_Lagrangian(const binvec x, const vector<binvec> &B, const unsigned n) {
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


vector<vector<double>> generate_projected_stabiliser_states_from_graphs(const string file, const unsigned n) {

	// needed later
	unsigned N = pow(2,n);
	vector<vector<double>> states;
	states.reserve(pow(3,n)); // estimated size
	vector<double> tmp_state (n,0);
	vector<binvec> B (n);
	vector<binvec> SB (n);


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
	vector<binvec_short> graphs = read_graph_states(file,n);

	unsigned counter = 1;
	unsigned graph_len = graphs.size();
	cout << "  Compute projected stabiliser states for graph" << endl;;
	for(auto th : graphs) {
		cout << "\r    #" << counter << " / " << graph_len << flush;
		B = generate_gs_Lagrangian(th,n);

		// we have to take the local Hadamard orbit of that Lagrangian
		for(auto S : LH_group) {
			for(unsigned k=0; k<n; k++) {
				SB.at(k) = matrix_vector_prod_mod2(S, B.at(k));
			}

			// generate stabiliser states for all possible sign choices and directly project them
			// note that the state is only added if it not already exists in the set states.
			for(binvec_short s=0; s<N; s++) {
				tmp_state = project_state(SB, s);

				// check if projected state tmp_state2 already exists
				// to do that, we use binary search with std::lower_bound() which gives an iterator on the first element that is not less than tmp_state2
				auto it = lower_bound(states.begin(), states.end(), tmp_state);
				if(it == states.end() || tmp_state < *it) {
					// the element is new
					states.insert(it, tmp_state);
				}
				// note that procedures preserves the ordering ... 
			}
		}
		++counter;
	}

	cout << endl;

	return states;
}

struct proj_helper {
	binvec_short a;
	char fock_index;
	char pauli_component;

	proj_helper(binvec_short aa, char f, char p) : a(aa), fock_index(f), pauli_component(p) {
	}
};

vector<vector<int>> generate_projected_stabiliser_states_from_graphs2(const string file, const unsigned n) {

	// outer loop variables
	unsigned N = pow(2,n);
	vector<vector<int>> states;
	states.reserve(pow(3,n)); // estimated size
	vector<binvec> B (n);
	vector<binvec> SB (n);

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
	vector<binvec_short> graphs = read_graph_states(file,n);

	unsigned counter = 1;
	unsigned graph_len = graphs.size();
	cout << "  Compute projected stabiliser states for graph" << endl;

	for(auto th : graphs) {
		cout << "\r    #" << counter << " / " << graph_len << flush;
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
					// the element is new
					states.insert(it, pr_state);
				}
				// note that procedures preserves the ordering ... 
			}
		}
		++counter;
	}

	cout << endl;

	return states;
}



#endif