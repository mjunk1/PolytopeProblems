#ifndef POINTGENERATOR_H
#define POINTGENERATOR_H
#include <string>
#include <vector>

using namespace std;

// ---------------------------------
// ----- Abstract generator classes
// ---------------------------------

class PointGenerator {
public:
	virtual vector<double> operator() (double p) const = 0;
};

class PSDistributionGenerator {
protected:
	unsigned _length;
public:
	virtual vector<double> operator() (double p) const = 0;

	unsigned get_length() {
		return _length;
	}
};

class PSDistributionGeneratorReal : public PSDistributionGenerator {};
class PSDistributionGeneratorFourier : public PSDistributionGenerator {};

#endif