#ifndef PRODUCTS_H
#define PRODUCTS_H
#include <initializer_list>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <cstdio>

#ifndef SYMPLECTIC_H
#include "symplectic.h"
#endif


// ----- for building up product states for the H case
vector<int> pr_psiplus = {0, -2};
vector<int> pr_psiminus = {0, 2};
string pr_psiplus_label = "A_ 01 01";
string pr_psiminus_label = "A_ 01 11";

vector<LabelledState> tri_factors = { LabelledState({0, 3, -1}, "Bo 102 000 c"), LabelledState({0, 3, 1},  "Bo 102 100 c") };
vector<LabelledState> quad_factors = { LabelledState({0, 2, -4, 1}, "CR 0011 0000 c"), LabelledState({0, 2, 4, 1},  "CR 0011 0011 c"), LabelledState({0, 6, 0, 0},  "CF 0001 0001 c"), LabelledState({0, 6, 0, 2},  "CF 0001 0000 c") };


vector<int> pr_allplus (unsigned n) {
	vector<int> ret (n,0);
	for(unsigned k=1; k<=n; k++) {
		ret.at(k-1) = binomial_coeff(n,k);
	}
	return ret;
}

string pr_allplus_label (unsigned n) {
	if(n==0) {
		return "";
	}
	string ret = adj_mat_to_graph6(vector<binvec>(n,0));
	ret += " ";
	for(unsigned k=0; k<n; k++) {
		ret += "0";
	}
	ret += " ";
	for(unsigned k=0; k<n; k++) {
		ret += "0";
	}
	return ret;
}

vector<int> pr_allminus (unsigned n) {
	vector<int> ret (n,0);
	for(unsigned k=1; k<=n; k++) {
		ret.at(k-1) = binomial_coeff(n,k);
		if(k%2 != 0) {
			ret.at(k-1) *= -1;
		}
	}
	return ret;
}

string pr_allminus_label (unsigned n) {
	if(n==0) {
		return "";
	}
	string ret = adj_mat_to_graph6(vector<binvec>(n,0));
	ret += " ";
	for(unsigned k=0; k<n; k++) {
		ret += "0";
	}
	ret += " ";
	for(unsigned k=0; k<n; k++) {
		ret += "1";
	}
	return ret;
}

// products of n-qubit states containing up to k psi+ (and + states else)
// states and labels will be added to the vectors in the arguments
int pr_psiplus_products(unsigned n, unsigned k, vector<vector<int>> &states, vector<string> &labels) {
	// a n-qubit state can hold at most n/2 Bell pairs
	if(k > n/2) {
		cout << "Error in pr_psiplus_products(): number of requested Bell pairs k is greater than n/2. Set k = n/2." << endl;
		k = n/2;
	}

	// add the new states to the vectors
	states.push_back(pr_allplus(n));
	labels.push_back(pr_allplus_label(n));

	for(unsigned i=1; i<=k; i++) {
		// building up the state with i psi+ factors
		vector<int> st = pr_psiplus;
		string lab = pr_psiplus_label;

		for(unsigned j=2; j<=i; j++) {
			st = pr_product_state(st, pr_psiplus);
			lab = pr_product_label(lab, pr_psiplus_label);
		}

		// this is a 2i qubit state, so we add n-2i plus states
		if(n-2*i > 0) {
			st = pr_product_state(st, pr_allplus(n-2*i));
			lab = pr_product_label(lab, pr_allplus_label(n-2*i));
		}
		states.push_back(st);
		labels.push_back(lab);
	}

	return 0;
}

int pr_psiplus_products2(unsigned n, unsigned k, vector<vector<int>> &states, vector<string> &labels) {
	// a n-qubit state can hold at most n/2 Bell pairs
	if(k > n/2) {
		cout << "Error in pr_psiplus_products(): number of requested Bell pairs k is greater than n/2. Set k = n/2." << endl;
		k = n/2;
	}

	// add the new states to the vectors
	states.push_back(pr_allminus(n));
	labels.push_back(pr_allminus_label(n));

	for(unsigned i=1; i<=k; i++) {
		// building up the state with i psi+ factors
		vector<int> st = pr_psiplus;
		string lab = pr_psiplus_label;

		for(unsigned j=2; j<=i; j++) {
			st = pr_product_state(st, pr_psiplus);
			lab = pr_product_label(lab, pr_psiplus_label);
		}

		// this is a 2i qubit state, so we add n-2i minus states
		if(n-2*i > 0) {
			st = pr_product_state(st, pr_allminus(n-2*i));
			lab = pr_product_label(lab, pr_allminus_label(n-2*i));
		}
		states.push_back(st);
		labels.push_back(lab);
	}

	return 0;
}


// product state approximations
int projected_Bell_products(const unsigned n, vector<LabelledState> &states) {
	// get all partitions of n/2 of length 3: n/2 = k_1 + k_2 + k_3
	// then use k_1 psi+, k_2 psi- and 2*k_3 (+1) +/- states
	
	for(unsigned k1=0; k1<=n/2; k1++) {
		for(unsigned k2=0; k2<=n/2-k1; k2++) {
			unsigned k3 = n/2-k1-k2;

			LabelledState state;

			for(unsigned i=0; i<k1; i++) {
				state.object = pr_product_state(pr_psiplus, state.object);
				state.label = pr_product_label(pr_psiplus_label, state.label);
			}
			for(unsigned i=0; i<k2; i++) {
				state.object = pr_product_state(pr_psiminus, state.object);
				state.label = pr_product_label(pr_psiminus_label, state.label);
			}

			LabelledState statem;

			if(2*k3+n%2 != 0) {
				state.object = pr_product_state(pr_allplus(2*k3+n%2), state.object);
				state.label = pr_product_label(pr_allplus_label(2*k3+n%2), state.label);

				statem.object = pr_product_state(pr_allminus(2*k3+n%2), statem.object);
				statem.label = pr_product_label(pr_allminus_label(2*k3+n%2), statem.label);

				states.push_back(statem);
			}

			// append
			states.push_back(state);			
		}
	}
}

int projected_tri_products(const unsigned n, vector<LabelledState> &states) {
	// get all partitions of n/3 of length 3: n/2 = k_1 + k_2 + k_3
	// then use k_1 copies of the first, k_2 of the second and fill the rest with Bell products

	// would be better to loop over r in the beginning!
	
	for(unsigned k1=0; k1<=n/3; k1++) {
		for(unsigned k2=0; k2<=n/3-k1; k2++) {

			// remaining qubits
			unsigned r = n-3*k1-3*k2;

			// get all Bell products of r qubits
			vector<LabelledState> bell_products;
			projected_Bell_products(r, bell_products);

			// build 3-qubit product state
			LabelledState state;

			for(unsigned i=0; i<k1; i++) {
				state = pr_product_state(tri_factors.at(0), state);
			}
			for(unsigned i=0; i<k2; i++) {
				state = pr_product_state(tri_factors.at(1), state);
			}

			// now multiply and append
			for(unsigned i=0; i<bell_products.size(); i++) {
				states.push_back( pr_product_state(state, bell_products.at(i)) );
			}	
		}
	}
}



// ----- for building up product states for the T case
vector<int> pr_Bell1_T = {0, -3};
vector<int> pr_Bell2_T = {0, 3};
string pr_Bell1_label_T = "A_ 01 11 c";
string pr_Bell2_label_T = "A_ 00 00 c";




#endif