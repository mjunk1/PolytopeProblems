#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>

// #include "symplectic.h"
#include "utilities.h"
#include "GLPKConvexSeparation.h"
#include "noise.h"


using namespace std;
using namespace Eigen;

// ----- Test of ConvexSeparation class

int main(int argc, char** argv) {

	// parameters
	unsigned n = 5;
	unsigned samples = 100;
	vector<double> p (n);
	double eps = 1e-6;


	// -------------- Simplices test -----------------------

	cout << "-----------------------------------------------" << endl;
	cout << "Testing ConvexSeparation Class using n-simplex" << endl;
	cout << "-----------------------------------------------" << endl;


	// set up a vertex matrix for a n-simplex
	SparseMatrix<double,RowMajor> cons(n+1,n+1);
	vector<Triplet<double>> triplets;

	for(unsigned i=0; i<=n; i++) {
		triplets.push_back(Triplet<double>(i,i,1));
	}

	cons.setFromTriplets(triplets.begin(),triplets.end());

	// transform it to COO format

	auto glpk_cons = to_GLPK_format(cons);
	
	cout << "Constraint Matrix for a " << n << "-simplex:" << endl;
	print(glpk_cons);
	
	// set up a ConvexSeparation object
	GLPKConvexSeparation lp (cons);
	lp.set_verbosity(false);

	cout << endl;

	// --------- Checking random samples from n-simplex
	cout << "-----------------------------------------------" << endl;
	cout << "> Testing random samples from the " << n << "-simplex ... " << endl;
	cout << endl;

	for(unsigned k=0; k<samples; k++) {

		p = random_distribution(n+1);

		cout << "Check random probability distribution" << endl;
		cout << "   p = (";
		for(auto v : p) {
			cout << v << ",";
		}
		cout << ")\t\t\t";
		
		lp.check_point(p);

		if ( lp.get_obj_value() > 0 ) {
			cout << "Error" << endl;
			break;
		}
		else {
			cout << "ok" << endl;
		}
	}
	cout << endl;

	// --------- Checking altered random samples from n-simplex
	cout << "-----------------------------------------------" << endl;
	cout << "> Testing l1-close random samples ... " << endl;
	cout << endl;

	random_device rd; 
	mt19937 gen(rd()); 
	uniform_int_distribution<> dist (0,n-1); 
	unsigned i;

	for(unsigned k=0; k<samples; k++) {

		p = random_distribution(n+1);
		i = dist(gen);
		p.at(i) += eps;

		cout << "Check random point " << endl;
		cout << "   p = (";
		for(auto v : p) {
			cout << v << ",";
		}
		cout << ")\t\t\t";
		
		lp.check_point(p);

		if ( lp.get_obj_value() > 0 ) {
			cout << "ok" << endl;
		}
		else {
			cout << "Error" << endl;
			break;
		}
	}



	// -------------- Clifford test -----------------------

	cout << "-----------------------------------------------" << endl;
	cout << "Testing ConvexSeparation Class using Clifford polytope" << endl;
	cout << "-----------------------------------------------" << endl;

	// TBD
	// this should contain some testing scenarios using the 2 qubit Clifford polytope

	GLPKConvexSeparation lp2 ("testing/constraints/polytope_twirled_constraints2Q.dat");
	int ret_status;
	lp2.set_verbosity(false);

	lp2.print_parameters();

	// noise model to use
	DepolarisingNoise dp_noise(2);

	// configure circuit 
	CliffordCircuit2Q circuit_gen;
	circuit_gen.set_elemental_noise(&dp_noise);

	// representation of noisy channel
	// NoisyTChannel y(2, &circuit_gen);
	NoisyTChannel2Q y(&circuit_gen);

	// solve & write solution for different scenarios
	double pth;


	// D_p^1
	circuit_gen.set_circuit(vector<string>({"H1"}));
	pth = lp2.check_family(y);
	cout << "D_p^1: Threshold value p_th = " << pth << endl;

	// D_p^2
	circuit_gen.set_circuit(vector<string>({"H1","H1"}));
	pth = lp2.check_family(y);
	cout << "D_p^2: Threshold value p_th = " << pth << endl;

}