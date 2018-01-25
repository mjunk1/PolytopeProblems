#ifndef STABILISER_H
#define STABILISER_H
#include <initializer_list>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <set>

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
	unsigned N = pow(4,n);
	vector<double> A(N*N);
	vector<unsigned> ai(n);
	vector<unsigned> bi(n);

	for(unsigned a=0; a<N; a++) {
		for(unsigned b=0; b<N; b++) {
			// decompose indices
			get_multi_index(n, 4, a, ai);
			get_multi_index(n, 4, b, bi);

			A.at(a*N+b) = 1;
			for(unsigned i=0; i<n; i++) {
				A.at(a*N+b) *= Tmatrix.at(bi.at(i)*4+ai.at(i));
			}
		}
	}

	return A;
}

vector<unsigned> index_set(unsigned k, unsigned n) {
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

vector<double> H_state_nb(const unsigned n, const double p = 0) {
	vector<double> ret (n);
	for(unsigned i=1; i<n+1; i++) {
		ret.at(i-1) = (1-p)*binomial_coeff(n,i);
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



// --- stabilisers 

vector<unsigned> generate_gs_Lagrangian(const unsigned i, const unsigned n) {
	// This generates the i-th symmetric nxn matrix and returns the representation of the corresponding graph state given in terms of a basis for the Lagrangian subspace
	// The basis vectors are given by the column vectors of the matrix
	// 			(  A  )
	// 		B = ( --- )			(in (z_1,...,z_n,x_1,...,x_n) coordinates)
	//			(  I  )
	// where A is a nxn symmetric matrix and I is the nxn identity matrix. Note that since we are working in product coordinates (z_1,x_1,z_2,x_2,...,z_n,x_n), the function will return a permuted version of these vectors.
	// The i-th vector from Z_2^{n(n+1)/2} is simply given by the binary representation of i

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

vector<double> generate_stabiliser_state(const vector<unsigned> B, const unsigned s) {
	// Generates the stabiliser state rho(B,s) that corresponds to the Lagrangian with basis B and a choice of signs (-1)^{s_i} where the bits s_i \in \Z_2 are specified in the length-n-bitstring s.
	// The stabiliser state is given in the Pauli basis {W(a)|a\in\Z_2^{2n}} where the only non-vanishing components rho_a corresponds to points a in the Lagrangian subspace and can thus be written in the basis B = (b_1,...,b_n) as
	//		a = \sum_{i=1}^n a_i b_i.
	// These components are explicitly given as
	//		rho_a = (-1)^{\sum_{i=1}^n a_i s_i} i^{\phi(a_1 b_1,...,a_n b_n)},
	// where \phi is the phase function appearing in the composition law of Pauli operators (see above).

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

		state.at(idx) = pow(-1,phase);// / pow(2,n);
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

	// cout << "Found " << d << " duplicates. " << size-d << " elements left." << endl;

	return states;
}



// -------- projection

vector<double> project_state(const vector<double> state, const unsigned n, vector<double> rot_matrix) {
	// size
	vector<double> pr_state(n, 0.); // 0th component is ommited since it is forced to be 1 
	vector<double> y = matrix_vector_prod(rot_matrix, state);

	// project
	for(unsigned k=1; k<n+1; k++) {
		vector<unsigned> aset = index_set(k,n);
		// transform the state
		for(auto a : aset) {
			pr_state.at(k-1) += y.at(a);
		}
		// pr_state.at(k) /= binomial_coeff(n,k);
	}

	return pr_state;
}

vector<double> project_state(const vector<double> state, const unsigned n) {
	// size
	vector<double> pr_state(n, 0.);

	// project
	for(unsigned k=1; k<n+1; k++) {
		vector<unsigned> aset = index_set(k,n);
		// transform the state
		for(auto a : aset) {
			pr_state.at(k-1) += state.at(a);
		}
		// pr_state.at(k) /= binomial_coeff(n,k);
	}

	return pr_state;
}

vector<double> project_state(vector<double> state, const unsigned n, bool rotate) {
	// size
	vector<double> pr_state(n, 0.);

	if(rotate == true) {
		// rotation matrix
		vector<double> A = rotation_matrix(n);

		pr_state = project_state(state, n, A);
	} 
	else {
		pr_state = project_state(state, n);
	}

	return pr_state;
}

vector<vector<double>> project_states(vector<vector<double>> states, const unsigned n) {
	// size
	unsigned nstates = states.size();
	vector<vector<double>> pr_states(nstates, vector<double>(n, 0.));

	// rotation matrix
	vector<double> A = rotation_matrix(n);

	for(unsigned i=0; i<nstates; i++) {
		// rotate
		vector<double> y = matrix_vector_prod(A, states.at(i));

		// project
		for(unsigned k=1; k<n+1; k++) {
			vector<unsigned> aset = index_set(k,n);
			// transform the i-th state
			for(auto a : aset) {
				pr_states.at(i).at(k-1) += y.at(a);
			}
			// pr_states.at(i).at(k) /= binomial_coeff(n,k);
		}
	}

	// eliminate duplicates
	sort(pr_states.begin(), pr_states.end());
	auto last = unique(pr_states.begin(), pr_states.end());
	pr_states.erase(last,pr_states.end());

	return pr_states;
}



// this function uses std::set which requires less memory for the states to store.
set<vector<double>> generate_projected_stabiliser_states_set(const unsigned n) {
	// estimated size (overcounted!)
	unsigned N = pow(2,n*(n+1)/2);
	unsigned M = pow(6,n);
	unsigned S = pow(2,n);
	unsigned size = N*M*S;

	// needed later
	set<vector<double>> states;
	vector<double> tmp_state (pow(4,n),0);
	vector<double> tmp_state2 (n,0);
	vector<unsigned> B;
	vector<unsigned> SB (n);
	vector<double> A = rotation_matrix(n);

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

			// generate stabiliser states for all possible sign choices and directly project them
			// note that the state is only added if it not already exists in the set states.
			for(unsigned s=0; s<S; s++) {
				tmp_state = generate_stabiliser_state(SB, s);
				tmp_state2 = project_state(tmp_state, n, A);
				states.insert(tmp_state2);
			}
		}
	}

	return states;
}

// this function uses a sorted std::vector which requires less memory for the states to store and should be faster than std::set
vector<vector<double>> generate_projected_stabiliser_states_vector(const unsigned n) {
	// estimated size (overcounted!)
	unsigned N = pow(2,n*(n+1)/2);
	unsigned M = pow(6,n);
	unsigned S = pow(2,n);
	unsigned size = N*M*S;

	// needed later
	vector<vector<double>> states;
	states.reserve(pow(3,n)); // estimated size
	vector<double> tmp_state (pow(4,n),0);
	vector<double> tmp_state2 (n,0);
	vector<unsigned> B;
	vector<unsigned> SB (n);
	vector<double> A = rotation_matrix(n);

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

			// generate stabiliser states for all possible sign choices and directly project them
			// note that the state is only added if it not already exists in the set states.
			for(unsigned s=0; s<S; s++) {
				tmp_state = generate_stabiliser_state(SB, s);
				tmp_state2 = project_state(tmp_state, n, A);

				// check if projected state tmp_state2 already exists
				// to do that, we use binary search with std::lower_bound() which gives an iterator on the first element that is not less than tmp_state2
				auto it = lower_bound(states.begin(), states.end(), tmp_state2);
				if(it == states.end() || tmp_state2 < *it) {
					// the element is new
					states.insert(it, tmp_state2);
				}
				// note that procedures preserves the ordering ... 
			}
		}
	}

	return states;
}

vector<vector<double>> generate_projected_stabiliser_states_vector_test(const unsigned n) {
	// estimated size (overcounted!)
	unsigned N = pow(2,n*(n+1)/2);
	unsigned M = pow(6,n);
	unsigned S = pow(2,n);
	unsigned size = N*M*S;

	// needed later
	vector<vector<double>> states;
	states.reserve(pow(3,n)); // estimated size
	vector<double> tmp_state (pow(4,n),0);
	vector<double> tmp_state2 (n,0);
	vector<unsigned> B;
	// vector<unsigned> SB (n);
	vector<double> A = rotation_matrix(n);

	// generate 1-qubit symplectic group
	// vector<vector<unsigned>> symp_group (6);

	// for(unsigned i=0; i<6; i++) {
	// 	symp_group.at(i) = generate_symplectic_matrix(i,1);
	// }

	// // generate local n-qubit symplectic group
	// vector<vector<unsigned>> loc_symp_group (M);
	// vector<unsigned> indices (n);
	// vector<vector<unsigned>> Slist (n);
	// for(unsigned i=0; i<M; i++) {
	// 	get_multi_index(n, 6, i, indices);
	// 	for(unsigned j=0; j<n; j++) {
	// 		Slist.at(j) = symp_group.at(indices.at(j));
	// 	}
	// 	loc_symp_group.at(i) = direct_sum(Slist);
	// }

	// graph state loop
	for(unsigned i=0; i<N; i++) {
		B = generate_gs_Lagrangian(i,n);

		// // local symplectic orbit
		// for(unsigned j=0; j<M; j++) {
		// 	// transform the basis vectors with the j-th matrix
		// 	for(unsigned k=0; k<n; k++) {
		// 		SB.at(k) = matrix_vector_prod_mod2(loc_symp_group.at(j), B.at(k));
		// 	}

			// generate stabiliser states for all possible sign choices and directly project them
			// note that the state is only added if it not already exists in the set states.
			for(unsigned s=0; s<S; s++) {
				tmp_state = generate_stabiliser_state(B, s);
				tmp_state2 = project_state(tmp_state, n, A);

				// check if projected state tmp_state2 already exists
				// to do that, we use binary search with std::lower_bound() which gives an iterator on the first element that is not less than tmp_state2
				auto it = lower_bound(states.begin(), states.end(), tmp_state2);
				if(it == states.end() || tmp_state2 < *it) {
					// the element is new
					states.insert(it, tmp_state2);
				}
				// note that procedures preserves the ordering ... 
			}
		// }
	}

	return states;
}


// This function tries to balance time and space complexity by using a std::vector that consumes about mem_size memory. Every couple of iterations, duplicates are deleted from that vector.
// mem_size should be given in MB
vector<vector<double>> generate_projected_stabiliser_states_vector2(const unsigned n, const unsigned mem_size) {
	// estimated size (overcounted!)
	unsigned N = pow(2,n*(n+1)/2);
	unsigned M = pow(6,n);
	unsigned S = pow(2,n);

	// estimate size of vector
	// mem_size*1024/8 = number of doubles that one can store in that memory
	unsigned size = mem_size*128 / n; 

	// needed later
	vector<vector<double>> states;
	states.reserve(size);
	vector<double> tmp_state (pow(4,n),0);
	vector<double> tmp_state2 (n,0);
	vector<unsigned> B;
	vector<unsigned> SB (n);
	vector<double> A = rotation_matrix(n);

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

			// generate stabiliser states for all possible sign choices and directly project them
			// note that the state is only added if it not already exists in the set states.
			for(unsigned s=0; s<S; s++) {
				tmp_state = generate_stabiliser_state(SB, s);
				tmp_state2 = project_state(tmp_state, n, A);

				// add state to the list
				states.push_back(tmp_state2);

				// if list is at maximum, delete duplicates
				if(states.size() == size) {
					// sort it and erase states
					sort(states.begin(), states.end());
					auto last = unique(states.begin(), states.end());
					states.erase(last,states.end());
				}
			}
		}
	}

	sort(states.begin(), states.end());
	auto last = unique(states.begin(), states.end());
	states.erase(last,states.end());

	states.shrink_to_fit();

	return states;
}

// this function is better suited since it needs exponentially less memory for the states to store.
vector<vector<double>> generate_projected_stabiliser_states_old(const unsigned n) {
	// estimated size (overcounted!)
	unsigned N = pow(2,n*(n+1)/2);
	unsigned M = pow(6,n);
	unsigned S = pow(2,n);
	unsigned size = N*M*S;

	// needed later
	vector<vector<double>> states ( size, vector<double> (n,0) );
	vector<double> tmp_state (pow(4,n),0);
	vector<unsigned> B;
	vector<unsigned> SB (n);
	vector<double> A = rotation_matrix(n);

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

			// generate stabiliser states for all possible sign choices and directly project them
			for(unsigned s=0; s<S; s++) {
				tmp_state = generate_stabiliser_state(SB, s);
				states.at(get_linear_index(3,{N,M,S},{i,j,s})) = project_state(tmp_state, n, A);
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


// this function uses a sorted std::vector which requires less memory for the states to store and should be faster than std::set
// however the advantage seems to diminuish for larger n
vector<vector<double>> generate_projected_stabiliser_states(const unsigned n) {
	return generate_projected_stabiliser_states_vector_test(n);
}

#endif