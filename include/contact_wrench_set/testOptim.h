#include <dwl/solver/OptimizationSolver.h>
#include <dwl/solver/IpoptNLP.h>
#include <contact_wrench_set/OptimizationProb.h>
#include <ros/package.h>

class testOptim {
public:
	void call_ipopt();
};