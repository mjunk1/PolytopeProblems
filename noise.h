#ifndef NOISE_H
#define NOISE_H
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <random>
#include <cassert>
#include <stdexcept>

#ifndef UTILITIES_H
#include "utilities.h"
#endif

#ifndef CONVEXSEPERATION_H
#include "ConvexSeperation.h"
#endif
 


// ----- Noise models

// Depolarising noise phase space distribution
double dn_1Q(double p, unsigned a) {
	if(a == 0)
		return (1.-p);
	else
		return p/3.;
}

// Symplectic Fourier transform of depolarising noise distribution
double dn_hat_1Q(double p, unsigned a) {
	if(a == 0)
		return 1.;
	else
		return ( (1.-p) - p/3. );
}


// Liouville-Pauli representation of noisy T channels

vector<double> Tmatrix = {1,0,0,0,0,1/sqrt(2),0,-1/sqrt(2),0,0,1,0,0,1/sqrt(2),0,1/sqrt(2)};

// T + random Pauli noise
int noisyT_1Q(double *p, double *y) {
	y[0] = p[0] - p[1] + p[2] - p[3];
	y[1] = y[2] = y[3] = 0;
	y[4] = p[0]/sqrt(2) + p[1]/sqrt(2) - p[2]/sqrt(2) - p[3]/sqrt(2);
	y[5] = -(p[0]/sqrt(2)) + p[1]/sqrt(2) + p[2]/sqrt(2) - p[3]/sqrt(2);
	y[6] = 0;
	y[7] = p[0]/sqrt(2) + p[1]/sqrt(2) - p[2]/sqrt(2) - p[3]/sqrt(2);
	y[8] = p[0]/sqrt(2) - p[1]/sqrt(2) - p[2]/sqrt(2) + p[3]/sqrt(2);
}

// T on first qubit + random Pauli noise on 2 qubits
int noisyT_2Q(double *p, double *y) {
	y[0] = 0;
	y[1] = 0;
	y[2] = 0;
	y[3] = 0;
	y[4] = p[0] + p[1] - p[10] - p[11] + p[12] + p[13] - p[14] - p[15] - p[2] - p[3] + p[4] + p[5] - p[6] - p[7] + p[8] + p[9];
	y[5] = 0;
	y[6] = 0;
	y[7] = 0;
	y[8] = p[0] - p[1] + p[10] - p[11] + p[12] - p[13] + p[14] - p[15] + p[2] - p[3] + p[4] - p[5] + p[6] - p[7] + p[8] - p[9];
	y[9] = 0;
	y[10] = 0;
	y[11] = 0;
	y[12] = p[0] - p[1] - p[10] + p[11] + p[12] - p[13] - p[14] + p[15] - p[2] + p[3] + p[4] - p[5] - p[6] + p[7] + p[8] - p[9];
	y[13] = 0;
	y[14] = 0;
	y[15] = 0;
	y[16] = 0;
	y[17] = (p[0] + p[1] - p[10] - p[11] - p[12] - p[13] - p[14] - p[15] + p[2] + p[3] + p[4] + p[5] + p[6] + p[7] - p[8] - p[9])/sqrt(2);
	y[18] = (p[0] + p[1] - p[10] - p[11] - p[12] - p[13] - p[14] - p[15] + p[2] + p[3] + p[4] + p[5] + p[6] + p[7] - p[8] - p[9])/sqrt(2);
	y[19] = 0;
	y[20] = 0;
	y[21] = (p[0] + p[1] + p[10] + p[11] - p[12] - p[13] + p[14] + p[15] - p[2] - p[3] + p[4] + p[5] - p[6] - p[7] - p[8] - p[9])/sqrt(2);
	y[22] = (p[0] + p[1] + p[10] + p[11] - p[12] - p[13] + p[14] + p[15] - p[2] - p[3] + p[4] + p[5] - p[6] - p[7] - p[8] - p[9])/sqrt(2);
	y[23] = 0;
	y[24] = 0;
	y[25] = (p[0] - p[1] - p[10] + p[11] - p[12] + p[13] - p[14] + p[15] + p[2] - p[3] + p[4] - p[5] + p[6] - p[7] - p[8] + p[9])/sqrt(2);
	y[26] = (p[0] - p[1] - p[10] + p[11] - p[12] + p[13] - p[14] + p[15] + p[2] - p[3] + p[4] - p[5] + p[6] - p[7] - p[8] + p[9])/sqrt(2);
	y[27] = 0;
	y[28] = 0;
	y[29] = (p[0] - p[1] + p[10] - p[11] - p[12] + p[13] + p[14] - p[15] - p[2] + p[3] + p[4] - p[5] - p[6] + p[7] - p[8] + p[9])/sqrt(2);
	y[30] = (p[0] - p[1] + p[10] - p[11] - p[12] + p[13] + p[14] - p[15] - p[2] + p[3] + p[4] - p[5] - p[6] + p[7] - p[8] + p[9])/sqrt(2);
	y[31] = 0;
	y[32] = 0;
	y[33] = -((p[0] + p[1] + p[10] + p[11] - p[12] - p[13] - p[14] - p[15] + p[2] + p[3] - p[4] - p[5] - p[6] - p[7] + p[8] + p[9])/sqrt(2));
	y[34] = (p[0] + p[1] + p[10] + p[11] - p[12] - p[13] - p[14] - p[15] + p[2] + p[3] - p[4] - p[5] - p[6] - p[7] + p[8] + p[9])/sqrt(2);
	y[35] = 0;
	y[36] = 0;
	y[37] = -((p[0] + p[1] - p[10] - p[11] - p[12] - p[13] + p[14] + p[15] - p[2] - p[3] - p[4] - p[5] + p[6] + p[7] + p[8] + p[9])/sqrt(2));
	y[38] = (p[0] + p[1] - p[10] - p[11] - p[12] - p[13] + p[14] + p[15] - p[2] - p[3] - p[4] - p[5] + p[6] + p[7] + p[8] + p[9])/sqrt(2);
	y[39] = 0;
	y[40] = 0;
	y[41] = -((p[0] - p[1] + p[10] - p[11] - p[12] + p[13] - p[14] + p[15] + p[2] - p[3] - p[4] + p[5] - p[6] + p[7] + p[8] - p[9])/sqrt(2));
	y[42] = (p[0] - p[1] + p[10] - p[11] - p[12] + p[13] - p[14] + p[15] + p[2] - p[3] - p[4] + p[5] - p[6] + p[7] + p[8] - p[9])/sqrt(2);
	y[43] = 0;
	y[44] = 0;
	y[45] = -((p[0] - p[1] - p[10] + p[11] - p[12] + p[13] + p[14] - p[15] - p[2] + p[3] - p[4] + p[5] + p[6] - p[7] + p[8] - p[9])/sqrt(2));
	y[46] = (p[0] - p[1] - p[10] + p[11] - p[12] + p[13] + p[14] - p[15] - p[2] + p[3] - p[4] + p[5] + p[6] - p[7] + p[8] - p[9])/sqrt(2);
	y[47] = 0;
	y[48] = 0;
	y[49] = 0;
	y[50] = 0;
	y[51] = p[0] + p[1] - p[10] - p[11] + p[12] + p[13] + p[14] + p[15] + p[2] + p[3] - p[4] - p[5] - p[6] - p[7] - p[8] - p[9];
	y[52] = 0;
	y[53] = 0;
	y[54] = 0;
	y[55] = p[0] + p[1] + p[10] + p[11] + p[12] + p[13] - p[14] - p[15] - p[2] - p[3] - p[4] - p[5] + p[6] + p[7] - p[8] - p[9];
	y[56] = 0;
	y[57] = 0;
	y[58] = 0;
	y[59] = p[0] - p[1] - p[10] + p[11] + p[12] - p[13] + p[14] - p[15] + p[2] - p[3] - p[4] + p[5] - p[6] + p[7] - p[8] + p[9];
	y[60] = 0;
	y[61] = 0;
	y[62] = 0;
	y[63] = p[0] - p[1] + p[10] - p[11] + p[12] - p[13] - p[14] + p[15] - p[2] + p[3] - p[4] + p[5] + p[6] - p[7] - p[8] + p[9];

	return 0;
}

vector<double> noisyT_2Q(double *p) {
	vector<double> y(63);
	y[0] = 0;
	y[1] = 0;
	y[2] = 0;
	y[3] = p[0] + p[1] - p[10] - p[11] + p[12] + p[13] - p[14] - p[15] - p[2] - p[3] + p[4] + p[5] - p[6] - p[7] + p[8] + p[9];
	y[4] = 0;
	y[5] = 0;
	y[6] = 0;
	y[7] = p[0] - p[1] + p[10] - p[11] + p[12] - p[13] + p[14] - p[15] + p[2] - p[3] + p[4] - p[5] + p[6] - p[7] + p[8] - p[9];
	y[8] = 0;
	y[9] = 0;
	y[10] = 0;
	y[11] = p[0] - p[1] - p[10] + p[11] + p[12] - p[13] - p[14] + p[15] - p[2] + p[3] + p[4] - p[5] - p[6] + p[7] + p[8] - p[9];
	y[12] = 0;
	y[13] = 0;
	y[14] = 0;
	y[15] = 0;
	y[16] = (p[0] + p[1] - p[10] - p[11] - p[12] - p[13] - p[14] - p[15] + p[2] + p[3] + p[4] + p[5] + p[6] + p[7] - p[8] - p[9])/sqrt(2);
	y[17] = (p[0] + p[1] - p[10] - p[11] - p[12] - p[13] - p[14] - p[15] + p[2] + p[3] + p[4] + p[5] + p[6] + p[7] - p[8] - p[9])/sqrt(2);
	y[18] = 0;
	y[19] = 0;
	y[20] = (p[0] + p[1] + p[10] + p[11] - p[12] - p[13] + p[14] + p[15] - p[2] - p[3] + p[4] + p[5] - p[6] - p[7] - p[8] - p[9])/sqrt(2);
	y[21] = (p[0] + p[1] + p[10] + p[11] - p[12] - p[13] + p[14] + p[15] - p[2] - p[3] + p[4] + p[5] - p[6] - p[7] - p[8] - p[9])/sqrt(2);
	y[22] = 0;
	y[23] = 0;
	y[24] = (p[0] - p[1] - p[10] + p[11] - p[12] + p[13] - p[14] + p[15] + p[2] - p[3] + p[4] - p[5] + p[6] - p[7] - p[8] + p[9])/sqrt(2);
	y[25] = (p[0] - p[1] - p[10] + p[11] - p[12] + p[13] - p[14] + p[15] + p[2] - p[3] + p[4] - p[5] + p[6] - p[7] - p[8] + p[9])/sqrt(2);
	y[26] = 0;
	y[27] = 0;
	y[28] = (p[0] - p[1] + p[10] - p[11] - p[12] + p[13] + p[14] - p[15] - p[2] + p[3] + p[4] - p[5] - p[6] + p[7] - p[8] + p[9])/sqrt(2);
	y[29] = (p[0] - p[1] + p[10] - p[11] - p[12] + p[13] + p[14] - p[15] - p[2] + p[3] + p[4] - p[5] - p[6] + p[7] - p[8] + p[9])/sqrt(2);
	y[30] = 0;
	y[31] = 0;
	y[32] = -((p[0] + p[1] + p[10] + p[11] - p[12] - p[13] - p[14] - p[15] + p[2] + p[3] - p[4] - p[5] - p[6] - p[7] + p[8] + p[9])/sqrt(2));
	y[33] = (p[0] + p[1] + p[10] + p[11] - p[12] - p[13] - p[14] - p[15] + p[2] + p[3] - p[4] - p[5] - p[6] - p[7] + p[8] + p[9])/sqrt(2);
	y[34] = 0;
	y[35] = 0;
	y[36] = -((p[0] + p[1] - p[10] - p[11] - p[12] - p[13] + p[14] + p[15] - p[2] - p[3] - p[4] - p[5] + p[6] + p[7] + p[8] + p[9])/sqrt(2));
	y[37] = (p[0] + p[1] - p[10] - p[11] - p[12] - p[13] + p[14] + p[15] - p[2] - p[3] - p[4] - p[5] + p[6] + p[7] + p[8] + p[9])/sqrt(2);
	y[38] = 0;
	y[39] = 0;
	y[40] = -((p[0] - p[1] + p[10] - p[11] - p[12] + p[13] - p[14] + p[15] + p[2] - p[3] - p[4] + p[5] - p[6] + p[7] + p[8] - p[9])/sqrt(2));
	y[41] = (p[0] - p[1] + p[10] - p[11] - p[12] + p[13] - p[14] + p[15] + p[2] - p[3] - p[4] + p[5] - p[6] + p[7] + p[8] - p[9])/sqrt(2);
	y[42] = 0;
	y[43] = 0;
	y[44] = -((p[0] - p[1] - p[10] + p[11] - p[12] + p[13] + p[14] - p[15] - p[2] + p[3] - p[4] + p[5] + p[6] - p[7] + p[8] - p[9])/sqrt(2));
	y[45] = (p[0] - p[1] - p[10] + p[11] - p[12] + p[13] + p[14] - p[15] - p[2] + p[3] - p[4] + p[5] + p[6] - p[7] + p[8] - p[9])/sqrt(2);
	y[46] = 0;
	y[47] = 0;
	y[48] = 0;
	y[49] = 0;
	y[50] = p[0] + p[1] - p[10] - p[11] + p[12] + p[13] + p[14] + p[15] + p[2] + p[3] - p[4] - p[5] - p[6] - p[7] - p[8] - p[9];
	y[51] = 0;
	y[52] = 0;
	y[53] = 0;
	y[54] = p[0] + p[1] + p[10] + p[11] + p[12] + p[13] - p[14] - p[15] - p[2] - p[3] - p[4] - p[5] + p[6] + p[7] - p[8] - p[9];
	y[55] = 0;
	y[56] = 0;
	y[57] = 0;
	y[58] = p[0] - p[1] - p[10] + p[11] + p[12] - p[13] + p[14] - p[15] + p[2] - p[3] - p[4] + p[5] - p[6] + p[7] - p[8] + p[9];
	y[59] = 0;
	y[60] = 0;
	y[61] = 0;
	y[62] = p[0] - p[1] + p[10] - p[11] + p[12] - p[13] - p[14] + p[15] - p[2] + p[3] - p[4] + p[5] + p[6] - p[7] - p[8] + p[9];

	return y;
}

vector<double> noisyT(unsigned n, vector<double> &phat) {
	// note that the very first entry in the Liouville representation is trivial for unital channels
	unsigned len = pow(4,n+1);
	vector<double> ret( len - 1 );

	unsigned j,b;

	for(unsigned i=1; i < len; i++) {
		// extract first 4 digits of index
		j = i >> (2U*n-2U);

		// extract the last 2n digits

		// set bit 2n and 2n+1
		b = 0;
		b |= (1U << 2U*n); 
		b |= (1U << (2U*n+1U)); // b = ...001100..00
		// negate, b = ...110011..11
		b = ~b; 
		// bitwise and with i, which sets bit 2n and 2n+1 to 0
		b = i&b;

		// now set the array
		ret.at(i-1) = Tmatrix.at(j) * phat.at(b);
	}
	return ret;
}


// ----- PSDistributionGenerator specialisations

class DepolarisingNoise2Q: public PSDistributionGenerator {
public:
	DepolarisingNoise2Q() {
		this->_length = 16;
	}

	vector<double> operator() (double p) {
		vector<double> pdist(16);
		pdist[0] = pow(1 - p,2);
		pdist[1] = ((1 - p)*p)/3.;
		pdist[2] = ((1 - p)*p)/3.;
		pdist[3] = ((1 - p)*p)/3.;
		pdist[4] = ((1 - p)*p)/3.;
		pdist[5] = pow(p,2)/9.;
		pdist[6] = pow(p,2)/9.;
		pdist[7] = pow(p,2)/9.;
		pdist[8] = ((1 - p)*p)/3.;
		pdist[9] = pow(p,2)/9.;
		pdist[10] = pow(p,2)/9.;
		pdist[11] = pow(p,2)/9.;
		pdist[12] = ((1 - p)*p)/3.;
		pdist[13] = pow(p,2)/9.;
		pdist[14] = pow(p,2)/9.;
		pdist[15] = pow(p,2)/9.;
	
		return pdist;
	}
};

class DepolarisingNoise: public PSDistributionGenerator {
protected:
	unsigned _n; // number of qubits
public:
	DepolarisingNoise() {
		DepolarisingNoise(1);
	}

	DepolarisingNoise(unsigned n) {
		_n = n;
		this->_length = pow(4,_n);
	}

	vector<double> operator() (double p) {
		vector<double> pdist(this->_length, 1.);
		vector<unsigned> mindex(_n);

		for(unsigned i=0; i < this->_length; i++) {
			get_multi_index(_n, 4, i, mindex);
			
			for(unsigned m=0; m < _n; m ++)
				pdist.at(i) *= dn_1Q(p, mindex.at(m));
		}
	
		return pdist;
	}
};

class DepolarisingNoiseHat: public PSDistributionGenerator {
protected:
	unsigned _n; // number of qubits
public:
	DepolarisingNoiseHat() {
		DepolarisingNoiseHat(1);
	}

	DepolarisingNoiseHat(unsigned n) {
		_n = n;
		this->_length = pow(2,2*_n);
	}

	vector<double> operator() (double p) {
		vector<double> pdist(this->_length, 1.);
		vector<unsigned> mindex(_n);

		for(unsigned i=0; i < this->_length; i++) {
			get_multi_index(_n, 4, i, mindex);
			
			for(unsigned m=0; m < _n; m ++)
				pdist.at(i) *= dn_hat_1Q(p, mindex.at(m));
		}
	
		return pdist;
	}
};


class CliffordCircuit2Q:  public PSDistributionGenerator {
protected:
	// elemental noise at every position in the circuit
	PSDistributionGenerator *_elemental_noise = nullptr;

	// Parameters used for random circuit generation
	// We use a Marsenne-Twister generator which is seeded by a random_device 
	unsigned _circuit_length;
	double _cnot_prob;
	double _hs_prob;
	random_device _rd; 
	mt19937 _rd_gen; 
	discrete_distribution<> _rd_dist;

	// Clifford gates

	// symplectic representations of elementary channels
	unsigned _H1[4] = { 0b1100, 0b0100, 0b0010, 0b0001 };
	unsigned _H2[4] = { 0b1000, 0b0100, 0b0011, 0b0001 };
	unsigned _S1[4] = { 0b0100, 0b1000, 0b0010, 0b0001 };
	unsigned _S2[4] = { 0b1000, 0b0100, 0b0001, 0b0010 };
	unsigned _CNOT[4] = { 0b1010, 0b0110, 0b0010, 0b1101 };
	unsigned _ID[4] = { 0b1000, 0b0100, 0b0010, 0b0001 };

	// making the gate access a bit easier
	vector<unsigned*> _gates = { _H1, _H2, _S1, _S2, _CNOT, _ID };
	map<string, unsigned> _gate_id = { {"H1", 0}, {"H2", 1}, {"S1", 2}, {"S2", 3}, {"CNOT", 4}, {"ID", 5} };
	vector<unsigned> _gate_order;

public:
	// default CTOR
	CliffordCircuit2Q() {
		this->_length = 16;

		_circuit_length = 1;
		_cnot_prob = 0;
		_hs_prob = (1.-_cnot_prob)/4.;
		_rd_dist = discrete_distribution<>( { _hs_prob, _hs_prob, _hs_prob, _hs_prob, _cnot_prob } );
		_rd_gen = mt19937(_rd());

		_gate_order.resize(_circuit_length);
	}

	CliffordCircuit2Q(PSDistributionGenerator *gen): CliffordCircuit2Q() {
		set_elemental_noise(gen);
	}

	// CTOR which prepares for random circuit generation
	CliffordCircuit2Q(unsigned circuit_length, double cnot_prob) {
		this->_length = 16;

		// cnot_prob is used to generate the distribution. It is assumed that H and S gates occur with the same probability
		_circuit_length = circuit_length;
		_cnot_prob = cnot_prob;
		_hs_prob = (1.-_cnot_prob)/4.;
		_rd_dist = discrete_distribution<>( { _hs_prob, _hs_prob, _hs_prob, _hs_prob, _cnot_prob } );
		_rd_gen = mt19937(_rd());

		_gate_order.resize(_circuit_length);
	}

	CliffordCircuit2Q(unsigned circuit_length, double cnot_prob, PSDistributionGenerator *gen): CliffordCircuit2Q(circuit_length, cnot_prob) {
		set_elemental_noise(gen);
	}

	// set elemental noise
	void set_elemental_noise(PSDistributionGenerator *gen) {
		assert(gen->get_length() == 16);
		_elemental_noise = gen;
	}

	// explicitely set circuit
	void set_circuit(vector<string> circuit) {
		// check validity of circuit
		if(circuit.size() <= 0)
			throw invalid_argument("Argument is an empty vector.");

		for(string s : circuit) 
			if(_gate_id.find(s) == _gate_id.end())
				throw invalid_argument(s + " is not a valid gate.");

		// save
		_circuit_length = circuit.size();
		_gate_order.resize(_circuit_length);

		for(unsigned i = 0; i < _circuit_length; i++) 
			_gate_order.at(i) = _gate_id[circuit.at(i)];
	}

	// generate random circuit
	void randomise() {
		_gate_order.resize(_circuit_length);

		for(unsigned i = 0; i < _circuit_length; i++)			
			_gate_order.at(i) = _rd_dist(_rd_gen);
	}

	void randomise(unsigned circuit_length, double cnot_prob) {
		_circuit_length = circuit_length;
		_cnot_prob = cnot_prob;
		_hs_prob = (1.-_cnot_prob)/4.;
		_rd_dist = discrete_distribution<>( { _hs_prob, _hs_prob, _hs_prob, _hs_prob, _cnot_prob } );
		randomise();
	}


	// generate phase space distribution
	vector<double> operator() (double p) {
		assert(_gate_order.size() > 0);

		vector<double> el_noise = (*_elemental_noise)(p);
		vector<double> dist0 = el_noise;
		vector<double> dist1 (16);

		// propagate through circuit
		for(unsigned i = 1; i < _circuit_length; i++) {
			symplectic_transform(dist0.data(), dist1.data(), _gates[_gate_order[i]], 2);
			convolve_mod2(dist1.data(), el_noise.data(), dist0.data(), 4);
		}

		return dist0;
	}	
};

class CliffordCircuitHat2Q:  public PSDistributionGenerator {
protected:
	// elemental noise at every position in the circuit
	PSDistributionGenerator *_elemental_noise = nullptr;

	// Parameters used for random circuit generation
	// We use a Marsenne-Twister generator which is seeded by a random_device 
	unsigned _circuit_length;
	double _cnot_prob;
	double _hs_prob;
	random_device _rd; 
	mt19937 _rd_gen; 
	discrete_distribution<> _rd_dist;

	// Clifford gates

	// symplectic representations of elementary channels
	unsigned _H1[4] = { 0b1100, 0b0100, 0b0010, 0b0001 };
	unsigned _H2[4] = { 0b1000, 0b0100, 0b0011, 0b0001 };
	unsigned _S1[4] = { 0b0100, 0b1000, 0b0010, 0b0001 };
	unsigned _S2[4] = { 0b1000, 0b0100, 0b0001, 0b0010 };
	unsigned _CNOT[4] = { 0b1010, 0b0110, 0b0010, 0b1101 };
	unsigned _ID[4] = { 0b1000, 0b0100, 0b0010, 0b0001 };

	// making the gate access a bit easier
	vector<unsigned*> _gates = { _H1, _H2, _S1, _S2, _CNOT, _ID };
	map<string, unsigned> _gate_id = { {"H1", 0}, {"H2", 1}, {"S1", 2}, {"S2", 3}, {"CNOT", 4}, {"ID", 5} };
	vector<unsigned> _gate_order;

public:
	// default CTOR
	CliffordCircuitHat2Q() {
		this->_length = 16;

		_circuit_length = 1;
		_cnot_prob = 0;
		_hs_prob = (1.-_cnot_prob)/4.;
		_rd_dist = discrete_distribution<>( { _hs_prob, _hs_prob, _hs_prob, _hs_prob, _cnot_prob } );
		_rd_gen = mt19937(_rd());

		_gate_order.resize(_circuit_length);
	}

	CliffordCircuitHat2Q(PSDistributionGenerator *gen): CliffordCircuitHat2Q() {
		set_elemental_noise(gen);
	}

	// CTOR which prepares for random circuit generation
	CliffordCircuitHat2Q(unsigned circuit_length, double cnot_prob) {
		this->_length = 16;

		// cnot_prob is used to generate the distribution. It is assumed that H and S gates occur with the same probability
		_circuit_length = circuit_length;
		_cnot_prob = cnot_prob;
		_hs_prob = (1.-_cnot_prob)/4.;
		_rd_dist = discrete_distribution<>( { _hs_prob, _hs_prob, _hs_prob, _hs_prob, _cnot_prob } );
		_rd_gen = mt19937(_rd());

		_gate_order.resize(_circuit_length);
	}

	CliffordCircuitHat2Q(unsigned circuit_length, double cnot_prob, PSDistributionGenerator *gen): CliffordCircuitHat2Q(circuit_length, cnot_prob) {
		set_elemental_noise(gen);
	}

	// set elemental noise
	void set_elemental_noise(PSDistributionGenerator *gen) {
		assert(gen->get_length() == 16);
		_elemental_noise = gen;
	}

	// explicitely set circuit
	void set_circuit(vector<string> circuit) {
		// check validity of circuit
		if(circuit.size() <= 0)
			throw invalid_argument("Argument is an empty vector.");

		for(string s : circuit) 
			if(_gate_id.find(s) == _gate_id.end())
				throw invalid_argument(s + " is not a valid gate.");

		// save
		_circuit_length = circuit.size();
		_gate_order.resize(_circuit_length);

		for(unsigned i = 0; i < _circuit_length; i++) 
			_gate_order.at(i) = _gate_id[circuit.at(i)];
	}

	// generate random circuit
	void randomise() {
		_gate_order.resize(_circuit_length);

		for(unsigned i = 0; i < _circuit_length; i++)			
			_gate_order.at(i) = _rd_dist(_rd_gen);
	}

	void randomise(unsigned circuit_length, double cnot_prob) {
		_circuit_length = circuit_length;
		_cnot_prob = cnot_prob;
		_hs_prob = (1.-_cnot_prob)/4.;
		_rd_dist = discrete_distribution<>( { _hs_prob, _hs_prob, _hs_prob, _hs_prob, _cnot_prob } );
		randomise();
	}


	// generate phase space distribution
	vector<double> operator() (double p) {
		assert(_gate_order.size() > 0);

		vector<double> el_noise = (*_elemental_noise)(p);
		vector<double> dist0 = el_noise;
		vector<double> dist1 (16);

		// propagate through circuit
		for(unsigned i = 1; i < _circuit_length; i++) {
			symplectic_transform(dist0.data(), dist1.data(), _gates[_gate_order[i]], 2);
			//convolve_mod2(dist1.data(), el_noise.data(), dist0.data(), 4);
			for(unsigned j = 0; j < 16; j++)
				dist0.at(j) = dist1.at(j) * el_noise.at(j);
		}

		return dist0;
	}	
};

// ----- PointGenerator specialisations

class NoisyTChannel2Q: public PointGenerator {
protected:
	vector<double> _ps_dist;
	PSDistributionGenerator *_ps_gen;

public:
	NoisyTChannel2Q(PSDistributionGenerator *gen) {
		set_ps_generator(gen);
	}

	vector<double> operator() (double p) {
		_ps_dist = _ps_gen->operator()(p);
		return noisyT_2Q(_ps_dist.data());
	}

	void set_ps_generator(PSDistributionGenerator *gen) {
		assert(gen->get_length() == 16);
		_ps_gen = gen;
	}

};

class NoisyTChannel: public PointGenerator {
protected:
	vector<double> _ps_dist;
	PSDistributionGenerator *_ps_gen = nullptr;
	unsigned _n;

public:
	NoisyTChannel(unsigned n, PSDistributionGenerator *gen) {
		_n = n;
		set_ps_generator(gen);
	}

	vector<double> operator() (double p) {
		assert(_ps_gen != nullptr);
		_ps_dist = (*_ps_gen)(p);
		return noisyT(_n, _ps_dist);
	}

	void set_ps_generator(PSDistributionGenerator *gen) {
		assert(gen->get_length() == pow(4,_n));
		_ps_gen = gen;
	}

};



#endif