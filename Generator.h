#ifndef GENERATOR_H
#define GENERATOR_H
#include <glpk.h>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <cassert>
#include <stdexcept>

#ifndef UTILITIES_H
#include "utilities.h"
#endif

using namespace std;


// ----- Abstract Generator classes

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


// ----- PointGenerator specialisations

class noisyT_channel_2Q: public PointGenerator {
protected:
	vector<double> _ps_dist;
	PSDistributionGenerator *_ps_gen;

public:
	noisyT_channel_2Q(PSDistributionGenerator *gen) {
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

#endif