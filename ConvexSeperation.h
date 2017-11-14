#ifndef CONVEXSEPERATION_H
#define CONVEXSEPERATION_H
#include <glpk.h>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>

#ifndef GENERATOR_H
#include "Generator.h"
#endif

using namespace std;


// ----- header file for solving convex seperation problems

class ConvexSeperation {
protected:
	// size parameters
	unsigned _nvertices = 0;
	unsigned _dim = 0;

public:

	// Constructors
	ConvexSeperation() { }

	// Destructor
	virtual ~ConvexSeperation() { }

	// abstract methods

	// reading & allocation method
	virtual int read_vertex_matrix(string cmatrix_file) = 0; 

	// single point checking method
	virtual int check_point(vector<double> &y) = 0;

	// given a 1-parameter family of points, y(p), compute the threshold value p = p_th for which the curve of points intersects the polytope and returns it
	virtual double check_family(PointGenerator &y) = 0;

	// write solution to outfile
	virtual void write_sol(string outfile) = 0;

	virtual double get_result() = 0;

	virtual int get_status() = 0;
};



class GLPKConvexSeperation: public ConvexSeperation {
protected:
	// size parameters
	unsigned _nentries = 0;

	// constraint matrix = coordinates of vertices
	vector<int> _ia;
	vector<int> _ja;
	vector<double> _a;	

	// coordinates of a point to check
	vector<double> _y;

	// GLPK variables
	glp_prob *_lp = nullptr;
	// glp_smcp _parm;
	int _glp_ret = 0;
	vector<int> _ind;
	string _method = "simplex";

	// interval division parameters
	double _lbnd = 0;
	double _ubnd = 1;
	unsigned _max_iter = 50;
	double _precision = 1e-6;

	// output
	bool _verbose = true;

public:

	// Constructors
	GLPKConvexSeperation() { }

	GLPKConvexSeperation(unsigned nvertices, unsigned dimension, string cmatrix_file) {

		_nvertices = nvertices;
		_dim = dimension;

		// preparing coordinates of the point to check
		_y.resize(_dim+2);
		_y.at(_dim+1) = -1.;

		// some helper stuff
		_ind.reserve(_dim+2);
		for(int i=0; i<=_dim+1; i++)
			_ind.push_back(i);

		// reading vertex coordinates
		if(read_vertex_matrix(cmatrix_file))
			exit (EXIT_FAILURE);
	}

	// Destructor
	~GLPKConvexSeperation() {
		glp_delete_prob(_lp);
	}

	// methods
	int read_vertex_matrix(string cmatrix_file) {

		// ----- read sparse constraint matrix

		// GLPK convention: data start at index 1
		_nentries = get_number_of_lines(cmatrix_file);
		_ia.resize(_nentries+1);
		_ja.resize(_nentries+1);
		_a.resize(_nentries+1);

		// open file and read
		fstream fin(cmatrix_file, ios::in);
		int n = 1;
		if(fin.is_open()) {
			while(fin >> _ia[n] >> _ja[n] >> _a[n])
				++n;
		} else {
			cerr << "Error: Couldn't open file " << cmatrix_file << endl;
			return 1;
		}

		fin.close();

		// ----- create GLPK problem 

		_lp = glp_create_prob();

		// set up some parameters


		// set up problem parameters
		glp_set_prob_name(_lp, "Polytope membership");
		glp_set_obj_dir(_lp, GLP_MAX);

		string s;

		// setting up rows
		glp_add_rows(_lp, _nvertices+1);
		for(int i=1; i<=_nvertices; i++) {
			s = "p" + to_string(i);
			glp_set_row_name(_lp, i, s.c_str());
			glp_set_row_bnds(_lp, i, GLP_UP, 0.0, 0.0);
		}
		glp_set_row_name(_lp, _nvertices+1, "q");
		glp_set_row_bnds(_lp, _nvertices+1, GLP_UP, 0.0, 1.0);

		// setting up cols 
		glp_add_cols(_lp, _dim+1);
		for(int j=1; j<=_dim; j++) {
			s = "x" + to_string(j);
			glp_set_col_name(_lp, j, s.c_str());
			glp_set_col_bnds(_lp, j, GLP_FR, 0.0, 0.0);
		}
		s = "x0";
		glp_set_col_name(_lp, _dim+1, s.c_str());
		glp_set_col_bnds(_lp, _dim+1, GLP_FR, 0.0, 0.0);

		// setting up constraint matrix
		glp_load_matrix(_lp, _nentries, _ia.data(), _ja.data(), _a.data());

		// setting last column
		vector<int> ind (_nvertices+2);
		vector<double> c (_nvertices+2);
		// ind.at(0) = 0;
		for(int i=1; i<=_nvertices+1; i++) {
			ind.at(i) = i;
			c.at(i) = -1.;
		}
		glp_set_mat_col(_lp, _dim+1, _nvertices+1, ind.data(), c.data());

		return 0;
	}

	void set_method(string s) {
		_method = s;
	}

	int check_point(vector<double> &y) {
		if(y.size() != _dim) {
			cout << "Error: check_point() received a vector of wrong size." << endl;
			return -1;
		}
		// copying y
		copy(y.begin(), y.end(), _y.begin()+1);

		// replacing the last row of the constraint matrix
		glp_set_mat_row(_lp, _nvertices+1, _dim+1, _ind.data(), _y.data());

		// setting up objective function
		for(int j=1; j<=_dim+1; j++) 
			glp_set_obj_coef(_lp, j, _y.at(j));


		// changes in the constraint matrix generally invalides the basis factorization
		//if(glp_bf_exists(_lp) == 0)
		//	glp_factorize(_lp);
		glp_warm_up(_lp);
		
		// solve
		if(_method == "simplex")
			_glp_ret = glp_simplex(_lp, NULL);
		else
			_glp_ret = glp_interior(_lp, NULL);

		return _glp_ret;
	}

	int set_parameters(double lbnd = 0, double ubnd = 1, unsigned max_iter = 50, double precision = 1e-6) {
		_lbnd = lbnd;
		_ubnd = ubnd;
		_max_iter = max_iter;
		_precision = precision;
	}

	double check_family(PointGenerator &y) {
		// initialize variables for interval division method
		double p1r = _ubnd;
		double p1l = _lbnd;
		double p0m = 0;
		double p1m = fabs(_ubnd + _lbnd)/2.;
		unsigned iter_counter = 0;
		vector<double> yy;

		// loop until desired accuracy is reached
		while( fabs(p1m-p0m) > _precision && iter_counter < _max_iter ) {
			
			if(_verbose == true)
				cout << "Solving convex seperation problem for y(p) with p = " << p1m << endl;

			// generate and check point
			yy = y(p1m);
			check_point(yy);

			// set new interval
			if(get_result() > 0) {
				// left of the optimal value
				p1l = p1m;
				p0m = p1m;
				p1m = fabs(p1r + p1l)/2.;

			} else {
				// right of the optimal value
				p1r = p1m;
				p0m = p1m;
				p1m = fabs(p1r + p1l)/2.;
			}

			++iter_counter;

			if(_verbose == true)
				cout << "-------------------------------------------------------------------" << endl;
		}

		if(iter_counter == _max_iter) {
			if(_verbose == true)
				cout << "Maximum number of iterations reached (max = " << _max_iter << ")" << endl;
		}
		else if(fabs(p1m-p0m) > _precision) {
			if(_verbose == true)
				cout << "Could not achieve required precision. Precision estimate = " << fabs(p1m-p0m) << endl;
		}
		else{
			if(_verbose == true)
				cout << "Converged to required precision." << endl;
		}

		return p1m;
	}

	void write_sol(string outfile) {
		glp_print_sol(_lp, outfile.c_str());
	}

	double get_result() {
		if(_method == "simplex")
			return glp_get_obj_val(_lp);
		else
			return glp_ipt_obj_val(_lp);
	}

	int get_status() {
		if(_method == "simplex")
			return glp_get_status(_lp);
		else
			return glp_ipt_status(_lp);
	}
};

// class SoPlexConvexSeperation : public ConvexSeperation {
	
// }

#endif