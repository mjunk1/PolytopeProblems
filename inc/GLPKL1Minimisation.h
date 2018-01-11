#ifndef GLPKL1MINIMISATION_H
#define GLPKL1MINIMISATION_H
#include <string>
#include <vector>

#ifndef GLPKFORMAT_H
#include "GLPKFormat.h"
#endif

#ifndef UTILITIES_H
#include "utilities.h"
#endif

#ifndef L1MINIMISATION_H
#include "L1Minimisation.h"
#endif

using namespace std;


// ---------------------------------------
// ----- GLPK specialisation
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


class GLPKL1Minimisation: public L1Minimisation {
protected:
	// size parameters
	// unsigned _nvertices; 
	unsigned _dim;
	
	// coordinates of a point to check
	vector<double> _y;

	// GLPK variables
	glp_prob *_lp = nullptr;
	glp_smcp _parm;
	int _glp_ret = 0;
	vector<int> _ind;
	string _method = "simplex";

	// output
	unsigned _verbose = 2;


	// this updates the GLPK problem using the GLPKFormat struct data
	void update_problem(GLPKFormat data) {
		// ----- create GLPK problem 

		unsigned nnz = data.non_zeros;
		unsigned nvertices = data.nrows;
		_dim = data.ncols;

		// create new problem
		if(_lp != nullptr) {
			glp_erase_prob(_lp);
		}
		else {
			_lp = glp_create_prob();
		}

		// set up problem parameters
		glp_set_prob_name(_lp, "Polytope membership");
		glp_set_obj_dir(_lp, GLP_MIN);

		string s;

		// setting up rows
		glp_add_rows(_lp, nvertices+1);
		for(int i=1; i<=nvertices; i++) {
			s = "p" + to_string(i);
			glp_set_row_name(_lp, i, s.c_str());
			glp_set_row_bnds(_lp, i, GLP_UP, 0.0, 0.0);
		}
		glp_set_row_name(_lp, nvertices+1, "q");
		glp_set_row_bnds(_lp, nvertices+1, GLP_UP, 0.0, 1.0);

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
		glp_load_matrix(_lp, nnz, data.rows.data(), data.cols.data(), data.values.data());

		// setting last column
		vector<int> ind (nvertices+2);
		vector<double> c (nvertices+2);
		// ind.at(0) = 0;
		for(int i=1; i<=nvertices+1; i++) {
			ind.at(i) = i;
			c.at(i) = -1.;
		}
		glp_set_mat_col(_lp, _dim+1, nvertices+1, ind.data(), c.data());



		// ----- additional preparation

		// preparing coordinates of the point to check
		_y.resize(_dim+2);
		_y.at(_dim+1) = -1.;

		// some helper stuff for GLPK
		_ind.reserve(_dim+2);
		for(int i=0; i<=_dim+1; i++)
			_ind.push_back(i);
	}

public:
	GLPKL1Minimisation() {
		// setting parameter struct to default values
		glp_init_smcp(&_parm);
	}

	GLPKL1Minimisation(string cmatrix_file) {
		// reading vertex coordinates and parameters
		if(read_vertex_matrix(cmatrix_file) != 0) {
			exit (EXIT_FAILURE);
		}

		// setting parameter struct to default values
		glp_init_smcp(&_parm);

	}

	GLPKL1Minimisation(SparseMatrix<double,RowMajor>& vertex_matrix) {
		// reading vertex coordinates and parameters
		if(read_vertex_matrix(vertex_matrix) != 0) {
			exit (EXIT_FAILURE);
		}

		// setting parameter struct to default values
		glp_init_smcp(&_parm);

	}

	// Destructor
	~GLPKL1Minimisation() {
		glp_delete_prob(_lp);
	}

	// input
	int read_vertex_matrix(string cmatrix_file) {
		GLPKFormat data = to_GLPK_format(cmatrix_file);
		update_problem(data);
		return 0;
	}

	int read_vertex_matrix(SparseMatrix<double,RowMajor>& vertex_matrix) {
		GLPKFormat data = to_GLPK_format(vertex_matrix);
		update_problem(data);
		return 0;
	}


	// operations
	int check_point(vector<double> &y) {
		assert(y.size() == _dim);
			
		// copying y
		copy(y.begin(), y.end(), _y.begin()+1);

		// make sure last element is -1
		_y.at(_dim+1) = -1.;

		// replacing the last row of the constraint matrix
		glp_set_mat_row(_lp, get_nvertices()+1, _dim+1, _ind.data(), _y.data());

		// setting up objective function
		for(int j=1; j<=_dim+1; j++) 
			glp_set_obj_coef(_lp, j, _y.at(j));


		// changes in the constraint matrix generally invalides the basis factorization
		//if(glp_bf_exists(_lp) == 0)
		//	glp_factorize(_lp);
		// glp_warm_up(_lp);
		glp_std_basis(_lp);
		
		// solve
		_glp_ret = glp_simplex(_lp, &_parm);

		return _glp_ret;
	}

	// get methods
	double get_obj_value() {
		return glp_get_obj_val(_lp);
		// if(_method == "simplex")
		// 	return glp_get_obj_val(_lp);
		// else
		// 	return glp_ipt_obj_val(_lp);
	}

	// vector<double> get_opt_solution() {
	// 	// TBD
	// }

	vector<double> get_vertex(unsigned i) {
		// returns a dense vector containing the coordinates of the i-th vertex
		assert(i < get_nvertices());

		vector<double> v (_dim, 0.);

		int *ind = new int[_dim+2];
		double *val = new double[_dim+2];
		int len = glp_get_mat_row(_lp, i+1, ind, val);

		for(unsigned j=1; j<=len; j++){
			if(ind[j] <= _dim)
				v.at(ind[j]-1) = val[j];
		}

		delete[] ind;
		delete[] val;

		return v;
	}

	int get_status() {
		if(_method == "simplex")
			return glp_get_status(_lp);
		else
			return glp_ipt_status(_lp);
	}

	unsigned get_nvertices() {
		return glp_get_num_rows(_lp)-1;
	}

	unsigned get_dimension() {
		return _dim;
	}

	unsigned get_nnz() {
		return glp_get_num_nz(_lp);
	}

	// set methods
	void set_verbosity(unsigned level) {
		// Sets verbosity level
		//   >= 2: All output
		//      1: Only errors
		//      0: No output 
		_verbose = level;
		if(_verbose == 0) {
			_parm.msg_lev = GLP_MSG_OFF;
		}
		else if(_verbose == 1) {
			_parm.msg_lev = GLP_MSG_ERR; // only error
		}
		else if(_verbose == 2) {
			_parm.msg_lev = GLP_MSG_ON;
		}
		else {
			_parm.msg_lev = GLP_MSG_ALL;
		}
	}

	void set_method(string s) {
		_method = s;
	}


	// output methods
	void write_sol(string outfile) {
		glp_print_sol(_lp, outfile.c_str());
	}

	void write_constraint_matrix(string outfile) {
		ofstream fout (outfile);

		int *ind = new int[_dim+2];
		double *val = new double[_dim+2];
		int len;

		if(fout.is_open()) {
			for(unsigned i = 1; i <= get_nvertices(); i++) {
				// get row
				len = glp_get_mat_row(_lp, i, ind, val);

				for(unsigned j=1; j<=len; j++) {
					if(ind[j] <= _dim) {
						fout << i << " " << ind[j] << " " << scientific << val[j] << endl;
					}
				}
			}
			fout.close();
		}
		else {
			cout << "Error in GLPKConvexSeparation::write_constraint_matrix : Couldn't open file " + outfile + " for writing" << endl;
		}

		delete[] ind;
		delete[] val;
	}

	void print_parameters() {
		cout << "---------------------------------------" << endl;
		cout << "GLPKL1Minimisation parameter output" << endl;
		cout << "---------------------------------------" << endl;

		cout << "Problem parameters:" << endl;
		cout << "  Dimension: " << get_dimension() << endl;
		cout << "  Vertices: " << get_nvertices() << endl;
		cout << "  Non-Zeros in the constraint matrix: " << get_nnz() << endl;
		cout << "---------------------------------------" << endl;
	}
};


#endif