#include <dwl/solver/OptimizationSolver.h>
#include <dwl/solver/IpoptNLP.h>
#include <OptimizationProblem.h>

int main(int argc, char **argv)
{
	dwl::solver::OptimizationSolver* solver = new dwl::solver::IpoptNLP();
	OptimizationProblem sample_opt;
	solver->setOptimizationModel(&sample_opt);

	solver->setFromConfigFile("../config/ipopt_config.yaml");//should be in the file
	solver->init();
	solver->compute();
	return 0;
}
