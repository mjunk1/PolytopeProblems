#ifndef UTILITIES_H
#define UTILITIES_H
#include <initializer_list>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

// ----- computational routines

int get_number_of_lines(string filename) {
	int n=0;
	string line;
	fstream fin(filename, ios::in);
	if(fin.is_open()) 
		while(getline(fin, line))
			++n;
	fin.close();
	return n;
}

// ----- indexing routines

/* Computes the multi index for a given linear index assuming row-major order 
	The array dimensions are assumed to be N_1 x N_2 x ... x N_k and are specified by arr_dim = k and the array arr_ranges = {N_1, ..., N_k}. The result is stored in indices (should be allocated before the use of this function!)
*/
int get_multi_index(const int arr_dim, const int* arr_ranges, const int index, int* indices) {

	int mod;
	int ind = index;
	if(arr_dim> 1) {
		for(int i=arr_dim-1; i>0; i--) {
			mod = ind % arr_ranges[i];
			indices[i] = mod;
			ind = (ind - mod) / arr_ranges[i];
		}
	}
	indices[0] = ind;

	return 0;
}

// overload where all N_i are assumed to be the same
int get_multi_index(const int arr_dim, const int arr_range, const int index, int* indices) {
	vector<int> vec (arr_dim, arr_range);
	return get_multi_index(arr_dim, vec.data(), index, indices);
}

int get_linear_index(const int arr_dim, const int* arr_ranges, const int* indices) {
	int index = 0;
	int jprod;
	for(int i=0; i<arr_dim; i++) {//change order?
		jprod = 1;
		for(int j=i+1; j<arr_dim; j++) {
			jprod *= arr_ranges[j];
		}
		jprod *= indices[i];
		index += jprod;
	}

	return index;
}

int get_linear_index(const int arr_dim, const int* arr_ranges, const initializer_list<int> indices){
	return get_linear_index(arr_dim, arr_ranges, const_cast<int* >(indices.begin()));
}

// overloads where all N_i are assumed to be the same
int get_linear_index(const int arr_dim, const int arr_range, const int* indices) {
	vector<int> vec (arr_dim, arr_range);
	return get_linear_index(arr_dim, vec.data(), indices);
}

int get_linear_index(const int arr_dim, const int arr_range, const initializer_list<int> indices) {
	vector<int> vec (arr_dim, arr_range);
	return get_linear_index(arr_dim, vec.data(), const_cast<int* >(indices.begin()));
}


#endif