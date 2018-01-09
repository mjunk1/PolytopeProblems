#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <set>
#include <chrono>

#include "symplectic.h"
#include "utilities.h"


using namespace std;

// ----- Test of stabiliser generation using std::set in symplectic.h

int main(int argc, char** argv) {

	unsigned m;

	if(argc > 1) {
		m = atoi(argv[1]);
	}
	else {
		m = 1;
	}

	cout << m << endl;

	cout << "Test direct generation of projected states" << endl;
	cout << "------------------------------------------" << endl;

	cout << "Generate projected stabiliser states for n = " << m << endl;


	// timing
	auto t1 = chrono::high_resolution_clock::now();

	set<vector<double>> pr_states = generate_projected_stabiliser_states_set(m);

	auto t2 = chrono::high_resolution_clock::now();	

    chrono::duration<double, milli> fp_ms = t2 - t1;


    cout << "Found " << pr_states.size() << " images." << endl;
    cout << "Generation took " <<  fp_ms.count() << " ms." << endl;
	

	fstream fout ("pr_stab_states_"+to_string(m)+".dat", ios::out);

	if(fout.is_open()) {
		for(auto state : pr_states) {
			for(auto x : state) {
				fout << x << " ";
			}
			fout << endl;
		}
	}
	fout.close();
}