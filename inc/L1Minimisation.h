#ifndef L1MINIMISATION_H
#define L1MINIMISATION_H
#include <string>
#include <vector>

using namespace std;


// ---------------------------------------
// ----- Abstract convex seperation class
// ---------------------------------------


class L1Minimisation {
protected:
	// size parameters
	// unsigned _nvertices = 0;
	unsigned _dim = 0;

public:

	// Constructors
	L1Minimisation() { }

	// Destructor
	virtual ~L1Minimisation() { }

	// abstract methods

	// reading & allocation method
	virtual int read_vertex_matrix(string cmatrix_file) = 0; 

	// single point checking method
	virtual int check_point(vector<double> &y) = 0;

	// write solution to outfile
	virtual void write_sol(string outfile) = 0;

	virtual double get_obj_value() = 0;

	virtual int get_status() = 0;
};


#endif