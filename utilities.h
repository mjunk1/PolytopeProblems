#ifndef UTILITIES_H
#define UTILITIES_H
#include <initializer_list>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>


using namespace std;

/* Utilities header file 
	For usage with the noise propagation code. 

	Caution!
	--------
	Due to the use of __builtin_popcount(), this code can only be used with the GNU compiler suite!
	Furthermore, the code uses at least the C++14 standard.
*/



// ----- indexing routines

/* Computes the multi index for a given linear index assuming row-major order 
	The array dimensions are assumed to be N_1 x N_2 x ... x N_k and are specified by arr_dim = k and the array arr_ranges = {N_1, ..., N_k}. The result is stored in indices (should be allocated before the use of this function!)
*/
unsigned get_multi_index(const unsigned arr_dim, const unsigned* arr_ranges, const unsigned index, unsigned* indices) {

	unsigned mod;
	unsigned ind = index;
	if(arr_dim> 1) {
		for(unsigned i=arr_dim-1; i>0; i--) {
			mod = ind % arr_ranges[i];
			indices[i] = mod;
			ind = (ind - mod) / arr_ranges[i];
		}
	}
	indices[0] = ind;

	return 0;
}

unsigned get_multi_index(const std::vector<unsigned> &arr_ranges, const unsigned index, vector<unsigned> &indices) {
	indices.reserve(arr_ranges.size());
	return get_multi_index(arr_ranges.size(), arr_ranges.data(), index, indices.data());
}

// overload where all N_i are assumed to be the same
unsigned get_multi_index(const unsigned arr_dim, const unsigned arr_range, const unsigned index, unsigned* indices) {
	vector<unsigned> vec (arr_dim, arr_range);
	return get_multi_index(arr_dim, vec.data(), index, indices);
}

unsigned get_multi_index(const unsigned arr_dim, const unsigned arr_range, const unsigned index, vector<unsigned> &indices) {
	vector<unsigned> vec (arr_dim, arr_range);
	indices.reserve(arr_dim);
	return get_multi_index(arr_dim, vec.data(), index, indices.data());
}

// inverse stuff

unsigned get_linear_index(const unsigned arr_dim, const unsigned* arr_ranges, const unsigned* indices) {
	unsigned index = 0;
	unsigned jprod;
	for(unsigned i=0; i<arr_dim; i++) {//change order?
		jprod = 1;
		for(unsigned j=i+1; j<arr_dim; j++) {
			jprod *= arr_ranges[j];
		}
		jprod *= indices[i];
		index += jprod;
	}

	return index;
}

unsigned get_linear_index(const unsigned arr_dim, const unsigned* arr_ranges, const initializer_list<unsigned> indices){
	return get_linear_index(arr_dim, arr_ranges, const_cast<unsigned* >(indices.begin()));
}

// overloads where all N_i are assumed to be the same
unsigned get_linear_index(const unsigned arr_dim, const unsigned arr_range, const unsigned* indices) {
	vector<unsigned> vec (arr_dim, arr_range);
	return get_linear_index(arr_dim, vec.data(), indices);
}

unsigned get_linear_index(const unsigned arr_dim, const unsigned arr_range, const initializer_list<unsigned> indices) {
	vector<unsigned> vec (arr_dim, arr_range);
	return get_linear_index(arr_dim, vec.data(), const_cast<unsigned* >(indices.begin()));
}



// ----- computational routines

unsigned get_number_of_lines(string filename) {
	unsigned n=0;
	string line;
	fstream fin(filename, ios::in);
	if(fin.is_open()) 
		while(getline(fin, line))
			++n;
	fin.close();
	return n;
}

unsigned popcount(unsigned i) {
	return __builtin_popcount(i);
}

unsigned parity(unsigned i) {
	return __builtin_parity(i);
}



// ----- Noise propagation

// probability distribution of depolarizing noise on 1 qubit
int dn_pdist_1Q(double p, double *pdist) {
  pdist[0] = (1.-p);
  pdist[1] = p/3.;
  pdist[2] = p/3.;
  pdist[3] = p/3.;
  return 0;
}

// probability distribution of depolarizing noise on 2 qubits
int dn_pdist_2Q(double p, double *pdist) {
  unsigned i,j,n;

  // i == 0, j == 0
  pdist[0] = (1.-p)*(1.-p);

  // i == 0, j > 0
  for(j=1; j<4; j++) {
    n = get_linear_index(2, 4, {0,j});
    pdist[n] = (1.-p)*p/3.;
  }   
  for(i=1; i<4; i++) {
    // i > 0, j == 0
    n = get_linear_index(2, 4, {i,0});
    pdist[n] = (1.-p)*p/3.;

    // i,j > 0
    for(j=1; j<4; j++) {
      n = get_linear_index(2, 4, {i,j});
      pdist[n] = p*p/9.;
    }     
  }
}

// This is the convolution in (Z_2)^n (addition mod 2)
int convolve_mod2(double *p1, double *p2, double *p3, const unsigned n) {
  unsigned len = pow(2U,n);
  unsigned i,j;

  for(i=0; i<len; i++) {
    p3[i] = 0;
    
    for(j=0; j<len; j++) {
      p3[i] += p1[j] * p2[i^j];
    } 
  }
  return 0;
}

// Computes the matrix vector product Ax over Z_2. It is assumed that a vector in Z_2^n is represented as a bitstring encoded as an unsigned integer. In this representation, every row of the matrix A corresponds to an unsigned integer, hence A is an unsigned array. The ouput is again an unsigned integer.
unsigned matrix_vector_prod_mod2(const unsigned *A, const unsigned x, const unsigned n) {
	unsigned b = 0;

	/* This performs the formula:
		b_i = sum_{j=1}^n A_ij x_j
		
		In the bitstring representation, this means that we have to multiply A[i] with x bitwise, afterwards adding up the resulting bits mod 2. The last step is just counting the number of ones and taking its modulus (here, called parity). Finally, we set the i-th bit of b to this result.
		Note that, in contrast to the usual convention, we count bits from the left since we see bitstrings as binary vectors. This means that the last bit is the least significant one.
	*/
	for (unsigned i=0; i<n; i++) 
		b ^= ((-parity( A[i]&x )) ^ b) & (1U << (n-i-1U)); 

	return b;	
}

int symplectic_transform(double *pin, double *pout, const unsigned *S, const unsigned nqubits) {
	unsigned j;
	for(unsigned i=0; i<pow(4,nqubits); i++) {
		j = matrix_vector_prod_mod2(S, i, 2U*nqubits);
		pout[j] = pin[i];
	}

	return 0;
}

int write_pdist(string outfile, double *p, unsigned nqubits) {
	
	unsigned len = pow(4,nqubits);
	unsigned bit;
	unsigned n = (unsigned)(2*nqubits);

	fstream fin (outfile, ios::out);

	if(fin.is_open()) {
		for(unsigned i=0; i<len; i++) {
			for(unsigned j=0; j<n; j++) {
				bit = (i >> (n-j-1U)) & 1U; // get the j-th bit of i
				fin << bit << " "; 
			}
			fin << p[i] << endl;		
		}
	} else {
		cerr << "Error: Couldn't open file " << outfile << endl;
		return 1;
	}

	fin.close();
	
	return 0;
}


// ----- Adjoint representation of noisy T channels


// T + random Pauli noise
int noisyT_1Q(double *p, double *y) {
	y[1] = p[0] - p[1] + p[2] - p[3];
	y[2] = y[3] = y[4] = 0;
	y[5] = p[0]/sqrt(2) + p[1]/sqrt(2) - p[2]/sqrt(2) - p[3]/sqrt(2);
	y[6] = -(p[0]/sqrt(2)) + p[1]/sqrt(2) + p[2]/sqrt(2) - p[3]/sqrt(2);
	y[7] = 0;
	y[8] = p[0]/sqrt(2) + p[1]/sqrt(2) - p[2]/sqrt(2) - p[3]/sqrt(2);
	y[9] = p[0]/sqrt(2) - p[1]/sqrt(2) - p[2]/sqrt(2) + p[3]/sqrt(2);
}

// T on first qubit + random Pauli noise on 2 qubits
int noisyT_2Q(double *p, double *y) {

	y[1] = p[0] - p[1] + p[12] - p[13] - p[14] + p[15] - p[2] + p[3]; 
	y[2] = 0;
	y[3] = 0;
	y[4] = 0;
	y[5] = 0;
	y[6] = p[10]/sqrt(3) - p[11]/sqrt(3) - p[4]/sqrt(3) + p[5]/sqrt(3) + p[6]/sqrt(3) - p[7]/sqrt(3) - p[8]/sqrt(3) + p[9]/sqrt(3);
	y[7] = 0;
	y[8] = 0;
	y[9] = 0;
	y[10] = 0;
	y[11] = -(sqrt(2./3.)*p[10]) + sqrt(2./3.)*p[11] + sqrt(2./3.)*p[4] - sqrt(2./3.)*p[5] - sqrt(2./3.)*p[6] + sqrt(2./3.)*p[7] + sqrt(2./3.)*p[8] - sqrt(2./3.)*p[9];
	y[12] = 0;
	y[13] = 0;
	y[14] = 0;
	y[15] = 0;
	y[16] = 0;
	y[17] = p[0] + p[1] + p[12] + p[13] - p[14] - p[15] - p[2] - p[3];
	y[18] = 0;
	y[19] = 0;
	y[20] = 0;
	y[21] = 0;
	y[22] = 0;
	y[23] = 0;
	y[24] = 0;
	y[25] = 0;
	y[26] = 0;
	y[27] = -p[10] - p[11] + p[4] + p[5] - p[6] - p[7] + p[8] + p[9];
	y[28] = 0;
	y[29] = 0;
	y[30] = 0;
	y[31] = 0;
	y[32] = 0;
	y[33] = p[0]/sqrt(2) - p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
	y[34] = 0;
	y[35] = 0;
	y[36] = 0;
	y[37] = 0;
	y[38] = p[1]/sqrt(2) - p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
	y[39] = -(p[0]/sqrt(2)) - p[11]/sqrt(2) + p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
	y[40] = 0;
	y[41] = 0;
	y[42] = 0;
	y[43] = 0;
	y[44] = -(p[1]/sqrt(2)) - p[10]/sqrt(2) + p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
	y[45] = 0;
	y[46] = 0;
	y[47] = 0;
	y[48] = 0;
	y[49] = p[0]/sqrt(2) + p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
	y[50] = 0;
	y[51] = 0;
	y[52] = p[1]/sqrt(2) + p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
	y[53] = 0;
	y[54] = 0;
	y[55] = -(p[1]/sqrt(2)) + p[11]/sqrt(2) + p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
	y[56] = 0;
	y[57] = 0;
	y[58] = -(p[0]/sqrt(2)) + p[10]/sqrt(2) + p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
	y[59] = 0;
	y[60] = 0;
	y[61] = 0;
	y[62] = 0;
	y[63] = 0;
	y[64] = 0;
	y[65] = p[0] - p[1] + p[12] - p[13] + p[14] - p[15] + p[2] - p[3];
	y[66] = 0;
	y[67] = 0;
	y[68] = 0;
	y[69] = 0;
	y[70] = 0;
	y[71] = 0;
	y[72] = 0;
	y[73] = 0;
	y[74] = 0;
	y[75] = p[10] - p[11] + p[4] - p[5] + p[6] - p[7] + p[8] - p[9];
	y[76] = p[10]/sqrt(3) - p[11]/sqrt(3) - p[4]/sqrt(3) + p[5]/sqrt(3) + p[6]/sqrt(3) - p[7]/sqrt(3) - p[8]/sqrt(3) + p[9]/sqrt(3);
	y[77] = 0;
	y[78] = 0;
	y[79] = 0;
	y[80] = 0;
	y[81] = p[0] + p[1]/3. - (2*p[10])/3. - (2*p[11])/3. + p[12] + p[13]/3. + p[14]/3. + p[15] + p[2]/3. + p[3] - (2*p[4])/3. - (2*p[5])/3. - (2*p[6])/3. - (2*p[7])/3. - (2*p[8])/3. - (2*p[9])/3.;
	y[82] = 0;
	y[83] = 0;
	y[84] = 0;
	y[85] = 0;
	y[86] = (2*sqrt(2)*p[1])/3. - (sqrt(2)*p[10])/3. - (sqrt(2)*p[11])/3. + (2*sqrt(2)*p[13])/3. + (2*sqrt(2)*p[14])/3. + (2*sqrt(2)*p[2])/3. - (sqrt(2)*p[4])/3. - (sqrt(2)*p[5])/3. - (sqrt(2)*p[6])/3. - (sqrt(2)*p[7])/3. - (sqrt(2)*p[8])/3. - (sqrt(2)*p[9])/3.;
	y[87] = 0;
	y[88] = 0;
	y[89] = 0;
	y[90] = 0;
	y[91] = 0;
	y[92] = 0;
	y[93] = 0;
	y[94] = p[1]/sqrt(2) + p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
	y[95] = 0;
	y[96] = 0;
	y[97] = p[0]/sqrt(2) + p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
	y[98] = 0;
	y[99] = 0;
	y[100] = -(p[0]/sqrt(2)) + p[10]/sqrt(2) + p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
	y[101] = 0;
	y[102] = 0;
	y[103] = -(p[1]/sqrt(2)) + p[11]/sqrt(2) + p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
	y[104] = 0;
	y[105] = 0;
	y[106] = 0;
	y[107] = 0;
	y[108] = p[1]/sqrt(2) - p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
	y[109] = 0;
	y[110] = 0;
	y[111] = 0;
	y[112] = 0;
	y[113] = p[0]/sqrt(2) - p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
	y[114] = -(p[1]/sqrt(2)) - p[10]/sqrt(2) + p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
	y[115] = 0;
	y[116] = 0;
	y[117] = 0;
	y[118] = 0;
	y[119] = -(p[0]/sqrt(2)) - p[11]/sqrt(2) + p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
	y[120] = 0;
	y[121] = 0;
	y[122] = 0;
	y[123] = p[0]/sqrt(2) - p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
	y[124] = 0;
	y[125] = 0;
	y[126] = 0;
	y[127] = 0;
	y[128] = p[1]/sqrt(2) - p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
	y[129] = p[0]/sqrt(2) + p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) - p[4]/sqrt(2) - p[7]/sqrt(2) + p[8]/sqrt(2);
	y[130] = 0;
	y[131] = 0;
	y[132] = 0;
	y[133] = 0;
	y[134] = p[1]/sqrt(2) + p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) - p[5]/sqrt(2) - p[6]/sqrt(2) + p[9]/sqrt(2);
	y[135] = 0;
	y[136] = 0;
	y[137] = 0;
	y[138] = 0;
	y[139] = p[1]/sqrt(2) + p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
	y[140] = 0;
	y[141] = 0;
	y[142] = p[0]/sqrt(2) + p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
	y[143] = 0;
	y[144] = 0;
	y[145] = p[0]/sqrt(2) - p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) - p[5]/sqrt(2) + p[6]/sqrt(2) + p[9]/sqrt(2);
	y[146] = 0;
	y[147] = 0;
	y[148] = p[1]/sqrt(2) - p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) - p[4]/sqrt(2) + p[7]/sqrt(2) + p[8]/sqrt(2);
	y[149] = 0;
	y[150] = 0;
	y[151] = -(sqrt(2./3.)*p[10]) + sqrt(2./3.)*p[11] + sqrt(2./3.)*p[4] - sqrt(2./3.)*p[5] - sqrt(2./3.)*p[6] + sqrt(2./3.)*p[7] + sqrt(2./3.)*p[8] - sqrt(2./3.)*p[9];
	y[152] = 0;
	y[153] = 0;
	y[154] = 0;
	y[155] = 0;
	y[156] = (2*sqrt(2)*p[1])/3. - (sqrt(2)*p[10])/3. - (sqrt(2)*p[11])/3. + (2*sqrt(2)*p[13])/3. + (2*sqrt(2)*p[14])/3. + (2*sqrt(2)*p[2])/3. - (sqrt(2)*p[4])/3. - (sqrt(2)*p[5])/3. - (sqrt(2)*p[6])/3. - (sqrt(2)*p[7])/3. - (sqrt(2)*p[8])/3. - (sqrt(2)*p[9])/3.;
	y[157] = 0;
	y[158] = 0;
	y[159] = 0;
	y[160] = 0;
	y[161] = p[0] - p[1]/3. - p[10]/3. - p[11]/3. + p[12] - p[13]/3. - p[14]/3. + p[15] - p[2]/3. + p[3] - p[4]/3. - p[5]/3. - p[6]/3. - p[7]/3. - p[8]/3. - p[9]/3.;
	y[162] = 0;
	y[163] = 0;
	y[164] = 0;
	y[165] = 0;
	y[166] = 0;
	y[167] = -p[10] - p[11] + p[4] + p[5] - p[6] - p[7] + p[8] + p[9];
	y[168] = 0;
	y[169] = 0;
	y[170] = 0;
	y[171] = 0;
	y[172] = 0;
	y[173] = 0;
	y[174] = 0;
	y[175] = 0;
	y[176] = 0;
	y[177] = p[0] + p[1] + p[12] + p[13] - p[14] - p[15] - p[2] - p[3];
	y[178] = 0;
	y[179] = 0;
	y[180] = 0;
	y[181] = 0;
	y[182] = 0;
	y[183] = 0;
	y[184] = p[0]/sqrt(2) + p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) + p[5]/sqrt(2) - p[6]/sqrt(2) - p[9]/sqrt(2);
	y[185] = 0;
	y[186] = 0;
	y[187] = p[1]/sqrt(2) + p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) + p[4]/sqrt(2) - p[7]/sqrt(2) - p[8]/sqrt(2);
	y[188] = 0;
	y[189] = 0;
	y[190] = p[1]/sqrt(2) - p[11]/sqrt(2) - p[13]/sqrt(2) + p[14]/sqrt(2) - p[2]/sqrt(2) - p[4]/sqrt(2) + p[7]/sqrt(2) + p[8]/sqrt(2);
	y[191] = 0;
	y[192] = 0;
	y[193] = p[0]/sqrt(2) - p[10]/sqrt(2) - p[12]/sqrt(2) + p[15]/sqrt(2) - p[3]/sqrt(2) - p[5]/sqrt(2) + p[6]/sqrt(2) + p[9]/sqrt(2);
	y[194] = 0;
	y[195] = 0;
	y[196] = 0;
	y[197] = 0;
	y[198] = p[1]/sqrt(2) - p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) + p[5]/sqrt(2) + p[6]/sqrt(2) - p[9]/sqrt(2);
	y[199] = 0;
	y[200] = 0;
	y[201] = 0;
	y[202] = 0;
	y[203] = p[0]/sqrt(2) - p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) + p[4]/sqrt(2) + p[7]/sqrt(2) - p[8]/sqrt(2);
	y[204] = p[1]/sqrt(2) + p[10]/sqrt(2) - p[13]/sqrt(2) - p[14]/sqrt(2) + p[2]/sqrt(2) - p[5]/sqrt(2) - p[6]/sqrt(2) + p[9]/sqrt(2);
	y[205] = 0;
	y[206] = 0;
	y[207] = 0;
	y[208] = 0;
	y[209] = p[0]/sqrt(2) + p[11]/sqrt(2) - p[12]/sqrt(2) - p[15]/sqrt(2) + p[3]/sqrt(2) - p[4]/sqrt(2) - p[7]/sqrt(2) + p[8]/sqrt(2);
	y[210] = 0;
	y[211] = 0;
	y[212] = 0;
	y[213] = 0;
	y[214] = 0;
	y[215] = p[10] - p[11] + p[4] - p[5] + p[6] - p[7] + p[8] - p[9];
	y[216] = 0;
	y[217] = 0;
	y[218] = 0;
	y[219] = 0;
	y[220] = 0;
	y[221] = 0;
	y[222] = 0;
	y[223] = 0;
	y[224] = 0;
	y[225] = p[0] - p[1] + p[12] - p[13] + p[14] - p[15] + p[2] - p[3];

	return 0;
}




#endif