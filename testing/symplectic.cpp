#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include "symplectic.h"
#include "utilities.h"


using namespace std;

// ----- Test of functions in symplectic.h

int main(int argc, char** argv) {

	unsigned n = 2;
	unsigned m = n;
	unsigned samples = 200;
	random_device rd; 
	mt19937 gen(rd()); 
	uniform_int_distribution<> dist;
	unsigned x,y,z;

	// ----- test of matrix functions

	// cout << "=================================================" << endl;
	// cout << "===   Start testing of matrix functions" << endl;
	// cout << "=================================================" << endl;

	// random_device rd; 
	// mt19937 gen(rd()); 
	// dist = uniform_int_distribution<> (0, pow(4,n)-1);
	

	// vector<unsigned> A1 (2*n);
	// vector<unsigned> A2 (2*m);

	// for(unsigned j=0; j<2*n; j++) {
	// 	A1.at(j) = dist(gen);
	// }
	// for(unsigned j=0; j<2*m; j++) {
	// 	A2.at(j) = dist(gen);
	// }

	// cout << "Generated random matrices:" << endl;
	// for(unsigned j=0; j<2*n; j++) {
	// 	cout << "     " << write_bits(A1.at(j),2*n) << endl; 
	// }

	// for(unsigned j=0; j<2*m; j++) {
	// 	cout << "     ";
	// 	for(unsigned j=0; j<2*n; j++)
	// 		cout << " ";
	// 	cout << write_bits(A2.at(j),2*m) << endl; 
	// }
	// cout << "Direct sum is:" << endl;
	// vector<unsigned> S = direct_sum(A1,A2);

	// for(unsigned j=0; j<(2*n+2*m); j++) {
	// 	cout << "     " << write_bits(S.at(j),2*n+2*m) << endl; 
	// }


	// ----- test of transvections

	// random generation of pairs of points
	// compute transvections for them and test if the result is correct

	cout << "=================================================" << endl;
	cout << "===   Start testing of transvections" << endl;
	cout << "=================================================" << endl;

	n = 4;
	dist = uniform_int_distribution<> (1, pow(4,n)-1);
	vector<unsigned> h;

	for(unsigned i=0; i<samples; i++) {
		x = dist(gen);
		y = dist(gen);
		h = find_transvection(x,y,n);

		z = transvection(h.at(1),x,n);
		z = transvection(h.at(0),z,n);


		cout << "  ** Transvection for random pair (x,y)=(" << write_bits(x,2*n) << "," << write_bits(y,2*n) << "):     h=(" << write_bits(h.at(0),2*n) << "," << write_bits(h.at(1),2*n) << ")  ...  ";

		if(z == y) {
			cout << "correct" << endl;
		}
		else {
			cout << "failed with status " << h.at(2) << endl;
			cerr << "Error: Transvection test not passed with random pair (x,y)=(" <<  write_bits(x,2*n) << "," << write_bits(y,2*n) << ")" << endl;
			cout << "Debug output: z=" << write_bits(h.at(3),2*n) << ", zz=" << write_bits(h.at(4),2*n) << endl;
			break;
		}
	}


	// ----- test of symplectic matrix generation


	// cout << "=================================================" << endl;
	// cout << "===   Start testing of symplectic group generation" << endl;
	// cout << "=================================================" << endl;

	n = 3;
	samples = 20;
	dist = uniform_int_distribution<> (0, sp_order(n));
	uniform_int_distribution<> dist2 (0,pow(4,n));
	unsigned k;
	unsigned p1,p2;
	vector<unsigned> S;

	for(unsigned k=0; k<samples; k++) {	
		S = generate_symplectic_matrix(dist(gen),n);
		cout << "Generated random symplectic matrix S:" << endl;
		for(unsigned j=0; j<2*n; j++) {
			cout << "     " << write_bits(S.at(j),2*n) << endl; 
		}

		cout << "Test if it is really symplectic ... " << endl;

		for(unsigned i=0; i<10*samples; i++) {
			x = dist2(gen);
			y = dist2(gen);

			p1 = symplectic_form(x,y,n);
			p2 = symplectic_form(matrix_vector_prod_mod2(S,x),matrix_vector_prod_mod2(S,y),n);

			cout << "  ** Test for random pair (x,y)=(" << write_bits(x,2*n) << "," << write_bits(y,2*n) << "):     [x,y]=" << p1 << ",   [Sx,Sy]=" << p2 << "  ...  ";
			if(p1 == p2) {
				cout << "correct" << endl;
			}
			else {
				cout << "failed " << endl;
				cerr << "Error: Symplectic test not passed with random pair (x,y)=(" <<  write_bits(x,2*n) << "," << write_bits(y,2*n) << ")" << endl;
				break;
			}
		}
	}

	// cout << "Generate the whole group ..." << endl;
	// string filename = "sp_" + to_string(2*n) + ".dat";
	// fstream fout(filename, ios::out);

	// for(unsigned j=0; j<sp_order(n); j++) {
	// 	S = generate_symplectic_matrix(j,n);
	// 	for(unsigned k=0; k<2*n; k++) {
	// 		fout << S.at(k) << " ";
	// 	}
	// 	fout << endl;
	// }

	// fout.close();



	// ----- test of Liouville matrix


	cout << "=================================================" << endl;
	cout << "===   Start testing of Liouville representation" << endl;
	cout << "=================================================" << endl;

	// vector<unsigned> M = { 0b01, 0b10 }; // H
	vector<unsigned> M = { 0b11, 0b01 }; // S
	// vector<unsigned> M = { 0b0100, 0b1000, 0b0010, 0b0001 }; // H1
	// vector<unsigned> M = { 0b1000, 0b0100, 0b0001, 0b0010 }; // H2
	// vector<unsigned> M = { 0b1100, 0b0100, 0b0010, 0b0001 }; // S1
	// vector<unsigned> M = { 0b1000, 0b0100, 0b0011, 0b0001 }; // S2
	// vector<unsigned> M = { 0b1010, 0b0100, 0b0010, 0b0101 }; // CNOT

	n = M.size()/2;
	unsigned i = 0;
	vector<int> L = liouville_matrix(M);
	for(unsigned a=0; a<pow(4,n); a++) {
		// a-th row
		for(unsigned b=0; b<pow(4,n); b++) {
			i = a << 2*n;
			i ^= b;
			// cout << i << " " << write_bits(i,4*n) << " ";
			cout << L.at(i) << " ";
		}
		cout << endl;
	}

	// cout << phi(0b10,0b01,1) << " " << phi(matrix_vector_prod_mod2(M,0b10),matrix_vector_prod_mod2(M,0b01),1) << " " << phi(0b10,0b11,1) << endl;
 	// cout << phase_function2(0b11,M) << endl;
}