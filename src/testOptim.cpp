#include <contact_wrench_set/testOptim.h>

void testOptim::call_ipopt()
{
	std::cout<<"hellooooo!"<<std::endl;
	dwl::solver::OptimizationSolver* solver = new dwl::solver::IpoptNLP();
	//std::shared_ptr< dwl::solver::OptimizationSolver> solver;
	OptimizationProb sample_opt;
	solver->setOptimizationModel(& sample_opt);
	//std::string pack_path = ros::package::getPath(std::string("contact_wrench_set"))+ "/config/ipopt_configLP.yaml";
    //solver->setFromConfigFile(pack_path);//should be in the folder
	solver->setFromConfigFile("/home/rorsolino/dls_ws/src/dls-distro/dls_utils/common_utilities/contact_wrench_set/config/ipopt_configLP.yaml");
	solver->init();
	solver->compute();
	Eigen::VectorXd solution = solver->getSolution();
	std::cout<<"The solution is: "<<solution.transpose()<<std::endl;
}
