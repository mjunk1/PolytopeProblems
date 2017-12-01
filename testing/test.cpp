#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

// #include "symplectic.h"
#include "utilities.h"


using namespace std;

// ----- Test of functions in utilities.h

int main(int argc, char** argv) {

	unsigned N = 10;
	unsigned n = 4;
	vector<vector<int>> list (N,vector<int>(n)); // creates a vector of N vector<int> with length n

	// fill it with random numbers
	random_device rd; 
	mt19937 gen(rd()); 
	uniform_int_distribution<> dist (1, n);

	for(auto& v : list) {
		for(unsigned i=0; i<n; i++) {
			v.at(i) = dist(gen);
		}
	}

	// print it
	cout << "random list:" << endl;
	for(auto& v : list) {
		for(unsigned i=0; i<n; i++) {
			cout << v.at(i) << " ";
		}
		cout << endl;
	}

	// let's sort it 
	sort(list.begin(),list.end());

	// print it
	cout << "sorted list:" << endl;
	for(auto& v : list) {
		for(unsigned i=0; i<n; i++) {
			cout << v.at(i) << " ";
		}
		cout << endl;
	}
}