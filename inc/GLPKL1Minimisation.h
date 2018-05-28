#ifndef GLPKL1MINIMISATION_H
#define GLPKL1MINIMISATION_H
#include <string>
#include <vector>
#include <glpk.h>
#include <cassert>

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



class GLPKL1Minimisation: public L1Minimisation {
protected:
	// size parameters
	// unsigned _nvertices; 
	unsigned _dim;
	

	// GLPK variables
	glp_prob *_lp = nullptr;
	glp_smcp _parm;
	int _glp_ret = 0;
	string _method = "simplex";

	// output
	unsigned _verbose = 2;


	// this updates the GLPK problem using the GLPKFormat struct data
	void update_problem(GLPKFormat& data) {
		// ----- create GLPK problem 

		unsigned nnz = data.non_zeros;
		unsigned nvertices = data.ncols;
		_dim = data.nrows;

		// create new problem
		if(_lp != nullptr) {
			glp_erase_prob(_lp);
		}
		else {
			_lp = glp_create_prob();
		}

		// set up problem parameters
		glp_set_prob_name(_lp, "L1 Minimisation");
		glp_set_obj_dir(_lp, GLP_MIN);

		string s;

		// setting up rows
		glp_add_rows(_lp, _dim + 2*nvertices + 1);
		for(int i=1; i<=_dim; i++) {
			s = "y" + to_string(i);
			glp_set_row_name(_lp, i, s.c_str());
			glp_set_row_bnds(_lp, i, GLP_FX, 0.0, 0.0); // fixed bound
		}
		for(int i=1; i<=nvertices; i++) {
			s = "r+" + to_string(i);
			glp_set_row_name(_lp, _dim+i, s.c_str());
			glp_set_row_bnds(_lp, _dim+i, GLP_UP, 0.0, 0.0); // upper bound 0
		}
		for(int i=1; i<=nvertices; i++) {
			s = "r-" + to_string(i);
			glp_set_row_name(_lp, _dim+nvertices+i, s.c_str());
			glp_set_row_bnds(_lp, _dim+nvertices+i, GLP_LO, 0.0, 0.0); // lower bound 0
		}

		// normalisation constraint
		s = "x norm.";
		glp_set_row_name(_lp, _dim+2*nvertices+1, s.c_str());
		glp_set_row_bnds(_lp, _dim+2*nvertices+1, GLP_FX, 1.0, 0.0); // fixed bound


		// setting up cols 
		glp_add_cols(_lp, 2*nvertices);
		for(int j=1; j<=nvertices; j++) {
			s = "x" + to_string(j);
			glp_set_col_name(_lp, j, s.c_str());
			glp_set_col_bnds(_lp, j, GLP_FR, 0.0, 0.0);
		}
		for(int j=1; j<=nvertices; j++) {
			s = "s" + to_string(j);
			glp_set_col_name(_lp, nvertices+j, s.c_str());
			glp_set_col_bnds(_lp, nvertices+j, GLP_FR, 0.0, 0.0);
		}

		// setting up objective function = 1^T.s
		vector<int> ind (nvertices+1, 1);
		for(int j=1; j<=nvertices; j++) 
			glp_set_obj_coef(_lp, nvertices+j, ind.at(j));


		// ----- setting up constraint matrix

		// load vertices
		glp_load_matrix(_lp, nnz, data.rows.data(), data.cols.data(), data.values.data());

		// setting up the rest
		ind = vector<int>(3,0);
		vector<double> val (3,0);
		val.at(1) = 1;

		for(int i=1; i<=nvertices; i++) {
			// cols to change
			ind.at(1) = i; 
			ind.at(2) = nvertices+i;

			// set row dim+i			
			val.at(2) = -1;
			glp_set_mat_row(_lp, _dim+i, 2, ind.data(), val.data());

			// set row dim+nvertices+i
			val.at(2) = 1;
			glp_set_mat_row(_lp, nvertices+_dim+i, 2, ind.data(), val.data());
		}

		// add another constraint such that the x's add up to 1
		ind = vector<int>(nvertices+1,0);
		val = vector<double>(nvertices+1,1);
		for(int i=1; i<=nvertices; i++) {
			ind.at(i) = i;
		}
		glp_set_mat_row(_lp, _dim+2*nvertices+1, nvertices, ind.data(), val.data());

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

	// GLPKL1Minimisation(SparseMatrix<double,RowMajor>& vertex_matrix) {
	// 	// reading vertex coordinates and parameters
	// 	if(read_vertex_matrix(vertex_matrix) != 0) {
	// 		exit (EXIT_FAILURE);
	// 	}

	// 	// setting parameter struct to default values
	// 	glp_init_smcp(&_parm);

	// }

	// Destructor
	~GLPKL1Minimisation() {
		glp_delete_prob(_lp);
	}

	// input
	int read_vertex_matrix(string cmatrix_file) {
		GLPKFormat data = to_GLPK_format(cmatrix_file, true);
		update_problem(data);
		return 0;
	}

	// int read_vertex_matrix(SparseMatrix<double,RowMajor>& vertex_matrix) {
	// 	GLPKFormat data = to_GLPK_format(vertex_matrix, true);
	// 	update_problem(data);
	// 	return 0;
	// }


	// operations
	int check_point(vector<double> &y) {
		assert(y.size() == _dim);

		// setting bounds
		for(int i=1; i<=_dim; i++) {
			glp_set_row_bnds(_lp, i, GLP_FX, y.at(i-1), 0.0); // fixed bound
		}

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

	int get_status() {
		if(_method == "simplex")
			return glp_get_status(_lp);
		else
			return glp_ipt_status(_lp);
	}

	unsigned get_nvertices() {
		return (glp_get_num_rows(_lp)-_dim)/2;
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
	void write_glpk_output(string outfile) {
		glp_print_sol(_lp, outfile.c_str());		
	}

	void write_sol(string outfile) {
		// write solution from GLPK object
		ofstream fout (outfile);
		unsigned M = get_nvertices();

		if(fout.is_open()) {
			for(unsigned i=1; i<=M; i++) {
				fout << glp_get_col_prim(_lp, i) << endl;
			}
			fout.close();
		}
		else {
			cout << "Error in GLPKL1Minimisation::write_sol : Couldn't open file " + outfile + " for writing" << endl;
		}

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

	void write_constraint_matrix(string outfile) {
		ofstream fout (outfile);

		int nrows = glp_get_num_rows(_lp);
		int nvertices = get_nvertices();

		int *ind = new int[nvertices+1];
		double *val = new double[nvertices+1];
		int len;

		if(fout.is_open()) {
			for(unsigned i = 1; i <= nrows; i++) {
				// get row
				len = glp_get_mat_row(_lp, i, ind, val);

				for(unsigned j=1; j<=len; j++) {
					fout << i << " " << ind[j] << " " << scientific << val[j] << endl;
				}
			}
			fout.close();
		}
		else {
			cout << "Error in GLPKL1Minimisation::write_constraint_matrix : Couldn't open file " + outfile + " for writing" << endl;
		}

		delete[] ind;
		delete[] val;
	}
};


#endif