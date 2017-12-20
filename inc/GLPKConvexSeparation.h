#ifndef GLPKCONVEXSEPERATION_H
#define GLPKCONVEXSEPERATION_H
#include <glpk.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <eigen3/Eigen/Sparse>

#ifndef UTILITIES_H
#include "utilities.h"
#endif

#ifndef CONVEXSEPERATION_H
#include "ConvexSeparation.h"
#endif

using namespace std;
using namespace Eigen;


// ----- GLPK utilities

// data struct for GLPK
struct GLPKFormat {
	vector<int> rows;
	vector<int> cols;
	vector<double> values;
	unsigned non_zeros;
	unsigned nrows;
	unsigned ncols;
};

// converts a Eigen::SparseMatrix to coordinate format
GLPKFormat to_GLPK_format(SparseMatrix<double,RowMajor>& other) {
	GLPKFormat ret;

	// number of non-zero values
	ret.non_zeros = other.nonZeros();
	ret.nrows = other.rows();
	ret.ncols = other.cols();

	// GLPK indexing starts at 1, we insert a 0 as 0th element.
	ret.rows = vector<int> (ret.non_zeros+1,0);
	ret.cols = vector<int> (ret.non_zeros+1,0);
	ret.values = vector<double> (ret.non_zeros+1,0.);

	// copy data
	for(unsigned k=0; k<other.outerSize(); k++) {
		for (typename SparseMatrix<double,RowMajor>::InnerIterator it(other,k); it; ++it) {
			// GLPK indexing starts at 1, hence row and column indices have to be shifted
			ret.rows.at(k+1) = it.row()+1;
			ret.cols.at(k+1) = it.col()+1;
			ret.values.at(k+1) = it.value();
		}  				
	}

	return ret;
}

// reads sparse matrix in COO format and saves to struct
GLPKFormat to_GLPK_format(string cmatrix_file) {
	GLPKFormat ret;

	// number of non-zero values
	ret.non_zeros = get_number_of_lines(cmatrix_file);

	// GLPK convention: data start at index 1
	ret.rows = vector<int> (ret.non_zeros+1,0);
	ret.cols = vector<int> (ret.non_zeros+1,0);
	ret.values = vector<double> (ret.non_zeros+1,0.);

	// open file and read
	fstream fin(cmatrix_file, ios::in);
	int n = 1;
	if(fin.is_open()) {
		while(fin >> ret.rows[n] >> ret.cols[n] >> ret.values[n]) {
			++n;
		}
	} else {
		cerr << "Error: Couldn't open file " << cmatrix_file << endl;
	}
	fin.close();

	// save number of rows and columns
	ret.nrows = ret.rows.at(ret.non_zeros);
	ret.ncols = *max_element(ret.cols.begin(),ret.cols.end());

	return ret;
}

void print(GLPKFormat data) {
	// GLPK = 1-based indexing
	unsigned k=1;
	for(unsigned i=1; i<=data.nrows; i++) {
		for(unsigned j=1; j<=data.ncols; j++) {
			if( (data.rows.at(k) == i) && (data.cols.at(k) == j) ) {
				cout << data.values.at(k) << " ";
				++k;
			}
			else{ 
				cout << 0 << " ";
			}
		}
		cout << endl;
	}
}


// ----- GLPK specialisation

class GLPKConvexSeparation: public ConvexSeparation {
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

	// interval division parameters
	double _lbnd = 0;
	double _ubnd = 1;
	unsigned _max_iter = 50;
	double _precision = 1e-6;

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
		glp_set_obj_dir(_lp, GLP_MAX);

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
	GLPKConvexSeparation() {
		// setting parameter struct to default values
		glp_init_smcp(&_parm);
	}

	GLPKConvexSeparation(string cmatrix_file) {
		// reading vertex coordinates and parameters
		if(read_vertex_matrix(cmatrix_file) != 0) {
			exit (EXIT_FAILURE);
		}

		// setting parameter struct to default values
		glp_init_smcp(&_parm);

	}

	GLPKConvexSeparation(SparseMatrix<double,RowMajor>& vertex_matrix) {
		// reading vertex coordinates and parameters
		if(read_vertex_matrix(vertex_matrix) != 0) {
			exit (EXIT_FAILURE);
		}

		// setting parameter struct to default values
		glp_init_smcp(&_parm);

	}

	// Destructor
	~GLPKConvexSeparation() {
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

	double check_family(PointGenerator &y) {
		// Note: This function should check for the correct "search direction" first. I.e. it has to identify which endpoint of the curve lies inside the polytope and which one lies outside. This is needed in order to correctly search of the intersection point parameter in the search interval [_lbnd, _ubnd].
		// Another note: There is certainly a smarter way to do this, too. This method can only handle a curve which starts outside the polytope, ends in the polytope and has only one intersection point with its boundary. That's enough for now, but for more general scenarios, this strategy has to be changed (e.g. starting from one endpoint and approaching the intersection point. As soon as it is passed, start interval division to find the precise location.)

		// initialize variables for interval division method
		double p1r = _ubnd;
		double p1l = _lbnd;
		double p0m = 0;
		double p1m = fabs(_ubnd + _lbnd)/2.;
		unsigned iter_counter = 0;
		vector<double> yy;

		// loop until desired accuracy is reached
		while( fabs(p1m-p0m) > _precision && iter_counter < _max_iter ) {
			
			if(_verbose >= 2)
				cout << "Solving convex seperation problem for y(p) with p = " << p1m << endl;

			// generate and check point
			yy = y(p1m);
			check_point(yy);

			// set new interval
			if(get_obj_value() > 0) {
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

			if(_verbose >= 2)
				cout << "-------------------------------------------------------------------" << endl;
		}

		if(iter_counter == _max_iter) {
			if(_verbose >= 1)
				cerr << "Maximum number of iterations reached (max = " << _max_iter << ")" << endl;
		}
		else if(fabs(p1m-p0m) > _precision) {
			if(_verbose >= 1)
				cerr << "Could not achieve required precision. Precision estimate = " << fabs(p1m-p0m) << endl;
		}
		else{
			if(_verbose >= 2)
				cout << "Converged to required precision." << endl;
		}

		return p1m;
	}

	void delete_point(unsigned number) {
		assert(number < get_nvertices());

		const int nums[] = {0, (int)(number+1)};
		glp_del_rows(_lp, 1, nums);
	}

	int delete_redundant_points(unsigned max_loops=1) {
		// Checks the generating set of the polytope for redundant points and deletes them. The remaining set consists of the vertices.
		//
		// If the parameter max_loops is given, the elimination procedure is repeated up to max_loops times until the number of points stabilises. This is just for testing purposes, since the algorithm should give the correct set of points after the first run.

		// indexing vars
		unsigned index;
		unsigned point;
		unsigned nvertices_tmp;
		unsigned nvertices_init = get_nvertices();
		unsigned cnt = 0;
		int ret;


		// used to temporarily save the coordinates of a point 
		vector<double> y (_dim, 0.);

		do {
			nvertices_tmp = get_nvertices();
			point = 0;

			while (point < get_nvertices()) {
				
				y = get_vertex(point);

				// temporarily set the i-th row to zero (remove the constraint corresponding to the vertex candidate)
				glp_set_mat_row(_lp, point+1, 0, NULL, NULL);


				// now check if our point is in the convex hull of all the other points
				ret = check_point(y);

				// note that the member _y now contains the row corresponding to y

				if(ret != 0) {
					if(_verbose >= 1) {
						cerr << "Error: LP solving process was not successful for point #" << point << endl;
					}

					// restore row
					glp_set_mat_row(_lp, point+1, _dim+1, _ind.data(), _y.data());

					return ret;
				}

				ret = get_status();
				if((ret != GLP_FEAS) && (ret != GLP_OPT)) {
					if(_verbose >= 1) {
						cerr << "Error: Couldn't find a feasible solution for point #" << point << endl;
					}

					// restore row
					glp_set_mat_row(_lp, point+1, _dim+1, _ind.data(), _y.data());

					return ret;
				}

				// check result and delete if necessary
				if(get_obj_value() > 0){
					// it is not redundant, restore row
					glp_set_mat_row(_lp, point+1, _dim+1, _ind.data(), _y.data());	

					// check the next point in the next loop
					++point;
				} 
				else {
					// it is redundant
					delete_point(point);

					// testing
					// cerr << "Deleted point #" << point << " with objective value " << scientific << get_obj_value() << endl;

					// do not increase point, since it refers now to the next index
					if(_verbose >= 2){
						cout << "Deleted point #" << point << endl;
					}
				}
			}

			++cnt;

			// cout << "Reduced number of vertices: " << get_nvertices() << endl;

		} while( (nvertices_tmp != get_nvertices()) && (cnt < max_loops) );

		if(cnt >= max_loops && max_loops > 1 && _verbose >= 1) {
			cerr << "Warning: Reached maximum number of iterations." << endl;
		}


		return 0;
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

	int set_parameters(double lbnd = 0, double ubnd = 1, unsigned max_iter = 50, double precision = 1e-6) {
		_lbnd = lbnd;
		_ubnd = ubnd;
		_max_iter = max_iter;
		_precision = precision;
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
						fout << i << " " << ind[j] << " " << val[j] << endl;
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
		cout << "GLPKConvexSeparation parameter output" << endl;
		cout << "---------------------------------------" << endl;

		cout << "Problem parameters:" << endl;
		cout << "  Dimension: " << get_dimension() << endl;
		cout << "  Vertices: " << get_nvertices() << endl;
		cout << "  Non-Zeros in the constraint matrix: " << get_nnz() << endl;

		cout << "Convex seperation parameters:" << endl;
		cout << "  Precision: " << _precision << endl;
		cout << "  Maximal iterations: " << _max_iter << endl;
		cout << "  Lower bound for interval division: " << _lbnd << endl;
		cout << "  Upper bound for interval division: " << _ubnd << endl;
		cout << "---------------------------------------" << endl;
	}
};






#endif