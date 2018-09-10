#ifndef CONVEXSEPERATION_H
#define CONVEXSEPERATION_H
#include <string>
#include <vector>

using namespace std;

#ifndef POINTGENERATOR_H
#include "PointGenerator.h"
#endif


// ---------------------------------------
// ----- Abstract convex seperation class
// ---------------------------------------


class ConvexSeparation {
protected:
	// size parameters
	// unsigned _nvertices = 0;
	unsigned _dim = 0;

public:

	// Constructors
	ConvexSeparation() { }

	// Destructor
	virtual ~ConvexSeparation() { }

	// abstract methods

	// reading & allocation method
	virtual int read_vertex_matrix(const string cmatrix_file) = 0; 

	// single point checking method
	virtual int check_point(const vector<double> &y) = 0;

	// given a 1-parameter family of points, y(p), compute the threshold value p = p_th for which the curve of points intersects the polytope and returns it
	virtual double check_family(const PointGenerator &y) = 0;

	// write solution to outfile
	virtual void write_sol(const string outfile) = 0;

	virtual double get_obj_value() = 0;

	virtual int get_status() = 0;
};


#endif