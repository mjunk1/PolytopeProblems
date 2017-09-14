#ifndef POLYTOPEPROGRAM_H
#define POLYTOPEPROGRAM_H
#include <glpk.h>
#include <fstream>
#include <string>
#include <cassert>

#ifndef UTILITIES_H
#include "utilities.h"
#endif

#define IPM

using namespace std;

class PolytopeProgram {
protected:
	// size parameters
	unsigned _nvertices = 0;
	unsigned _dim = 0;
	unsigned _nentries = 0;

	// constraint matrix = coordinates of vertices
	int *_ia = nullptr;
	int *_ja = nullptr;
	double *_a = nullptr;	

	// coordinates of a point to check
	double *_y = nullptr;

public:

	// Constructors
	PolytopeProgram() { }

	// Destructor
	virtual ~PolytopeProgram() {
		delete[] _ia, _ja, _a, _y;
	}

	// abstract methods

	// reading & allocation method
	virtual int read_vertex_matrix(string cmatrix_file) = 0; 

	// membership checking method
	virtual int check_point(double *y) = 0;

	// write solution to outfile
	virtual void write_sol(string outfile) = 0;

	virtual double get_result() = 0;

	virtual int get_status() = 0;
};



class GLPKPolytopeProgram : public PolytopeProgram {
protected:

	// GLPK variables
	glp_prob *_lp = nullptr;
	// glp_smcp _parm;
	int _glp_ret = 0;
	int *_ind;

	string _method = "simplex";

public:

	// Constructors
	GLPKPolytopeProgram() { }

	GLPKPolytopeProgram(unsigned nvertices, unsigned dimension, string cmatrix_file) {

		_nvertices = nvertices;
		_dim = dimension;

		// preparing coordinates of the point to check
		_y = new double[_dim+2];
		_y[0] = 0.;
		_y[_dim+1] = -1.;

		// some helper stuff
		_ind = new int[_dim+2];
		for(int i=0; i<=_dim+1; i++)
			_ind[i] = i;

		// reading vertex coordinates
		if(read_vertex_matrix(cmatrix_file))
			exit (EXIT_FAILURE);
	}

	// Destructor
	~GLPKPolytopeProgram() {
		delete[] _ind;
		glp_delete_prob(_lp);
	}

	// methods
	int read_vertex_matrix(string cmatrix_file) {

		// ----- read sparse constraint matrix

		// GLPK convention: data start at index 1
		_nentries = get_number_of_lines(cmatrix_file);
		_ia = new int[_nentries+1];
		_ja = new int[_nentries+1];
		_a = new double[_nentries+1];
		_ia[0] = _ja[0] = 0;
		_a[0] = 0;

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

		char s[50];

		// setting up rows
		glp_add_rows(_lp, _nvertices+1);
		for(int i=1; i<=_nvertices; i++) {
			sprintf(s, "p%u", i);
			glp_set_row_name(_lp, i, s);
			glp_set_row_bnds(_lp, i, GLP_UP, 0.0, 0.0);
		}
		glp_set_row_name(_lp, _nvertices+1, "q");
		glp_set_row_bnds(_lp, _nvertices+1, GLP_UP, 0.0, 1.0);

		// setting up cols 
		glp_add_cols(_lp, _dim+1);
		for(int j=1; j<=_dim; j++) {
			sprintf(s, "x%u", j);
			glp_set_col_name(_lp, j, s);
			glp_set_col_bnds(_lp, j, GLP_FR, 0.0, 0.0);
		}
		sprintf(s, "x0");
		glp_set_col_name(_lp, _dim+1, s);
		glp_set_col_bnds(_lp, _dim+1, GLP_FR, 0.0, 0.0);

		// setting up constraint matrix
		glp_load_matrix(_lp, _nentries, _ia, _ja, _a);

		// setting last column
		int *ind = new int[_nvertices+2];
		double *c = new double[_nvertices+2];
		ind[0] = 0;
		for(int i=1; i<=_nvertices+1; i++) {
			ind[i] = i;
			c[i] = -1.;
		}
		glp_set_mat_col(_lp, _dim+1, _nvertices+1, ind, c);

		delete[] ind;
		delete[] c;

		return 0;
	}

	void set_method(string s) {
		_method = s;
	}

	int check_point(double *y) {
		if(y == nullptr) {
			cout << "Error: check_point() received a null pointer." << endl;
			return -1;
		}
		// copying y
		copy(y,y+_dim+1,_y);

		// replacing the last row of the constraint matrix
		glp_set_mat_row(_lp, _nvertices+1, _dim+1, _ind, _y);

		// setting up objective function
		for(int j=1; j<=_dim+1; j++) 
			glp_set_obj_coef(_lp, j, _y[j]);


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

// class SoPlexPolytopeProgram : public PolytopeProgram {
	
// }

#endif