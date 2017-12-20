#ifndef CONVEXSEPERATION_H
#define CONVEXSEPERATION_H
#include <string>
#include <vector>

using namespace std;


// ---------------------------------
// ----- Abstract generator classes
// ---------------------------------

class PointGenerator {
public:
	virtual vector<double> operator() (double p) = 0;
};

class PSDistributionGenerator {
protected:
	unsigned _length;
public:
	virtual vector<double> operator() (double p) = 0;

	unsigned get_length() {
		return _length;
	}
};

class PSDistributionGeneratorReal : public PSDistributionGenerator {};
class PSDistributionGeneratorFourier : public PSDistributionGenerator {};

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
	virtual int read_vertex_matrix(string cmatrix_file) = 0; 

	// single point checking method
	virtual int check_point(vector<double> &y) = 0;

	// given a 1-parameter family of points, y(p), compute the threshold value p = p_th for which the curve of points intersects the polytope and returns it
	virtual double check_family(PointGenerator &y) = 0;

	// write solution to outfile
	virtual void write_sol(string outfile) = 0;

	virtual double get_obj_value() = 0;

	virtual int get_status() = 0;
};


#endif