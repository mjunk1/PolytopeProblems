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

	cout << "----------------------------------" << endl;
	cout << "Testing ConvexSeparation methods" << endl;
	cout << "----------------------------------" << endl;


	// ------ set up a vertex matrix for a n-simplex
	SparseMatrix<double,RowMajor> cons(n+1,n+1);
	vector<Triplet<double>> triplets;

	for(unsigned i=0; i<=n; i++) {
		triplets.push_back(Triplet<double>(i,i,1));
	}
	cons.setFromTriplets(triplets.begin(),triplets.end());

	cout << "Constraint Matrix for a " << n << "-simplex:" << endl;
	cout << cons << endl;
	

	// ------ set up a ConvexSeparation object
	GLPKConvexSeparation lp;
	lp.set_verbosity(1);

	cout << "----------------------------------" << endl;
	cout << "Testing input methods ...." << endl;

	cout << "Reading input ... ";
	if(lp.read_vertex_matrix(cons) != 0) {
		cout << "failed" << endl;
	}
	else {
		cout << "ok" << endl;
	}

	cout << "----------------------------------" << endl;


	cout << "----------------------------------" << endl;
	cout << "Testing get methods ...." << endl;

	cout << "get_dimension() = " << lp.get_dimension() << " ... ";
	if(lp.get_dimension() == n+1) {
		cout << "ok" << endl;
	}
	else {
		cout << "failed" << endl;
	}

	cout << "get_nvertices() = " << lp.get_nvertices() << " ... ";
	if(lp.get_nvertices() == n+1) {
		cout << "ok" << endl;
	}
	else {
		cout << "failed" << endl;
	}

	cout << "get_nnz() = " << lp.get_nnz() << " ... ";
	if(lp.get_nnz() == n+1+n+2) {
		cout << "ok" << endl;
	}
	else {
		cout << "failed" << endl;
	}

	cout << "get_vertex() ... ";

	bool success = true;
	vector<double> v(n+1);
	for(unsigned i=0; i<n+1; i++) {
		v.at(i) = 1;
		if(lp.get_vertex(i) != v) {
			success = false;
		}
		v.at(i) = 0;
	}
	if(success == true) {
		cout << "ok" << endl;
	}
	else {
		cout << "failed" << endl;
	}


	cout << "----------------------------------" << endl;
	

	cout << "----------------------------------" << endl;
	cout << "Testing get methods ...." << endl;



	cout << "----------------------------------" << endl;

}