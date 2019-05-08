 /*  Created on: October 12, 2017
 *      Author: Romeo Orsolino
 */
#include <contact_wrench_set/get_weights_lp.h>
#include <math_utils/utils.h>

using namespace Eigen;
using namespace iit;

OptWeights::OptWeights() {}


OptWeights::~OptWeights() {}


void OptWeights::init(bool soft_constraint)
{

	// Getting the dimension of the decision vector
	state_dimension_ = colonne+1; // the number of weights  plus the margin
	std::cout<<"state dimension: "<<state_dimension_<<std::endl;

	constraint_dimension_ = righe+1;
	std::cout<<"constraint dimension: "<<constraint_dimension_<<std::endl;


}


void OptWeights::getStartingPoint(double* decision, int decision_dim)
{
// Eigen interfacing to raw buffers
	Eigen::Map<Eigen::VectorXd> full_initial_point(decision, decision_dim);
	
	(full_initial_point.segment(0,righe)).setConstant(mean_coeff);
	full_initial_point(righe) = 0.0;
}

void OptWeights::evaluateBounds(double* decision_lbound, int decision_dim1,
										 double* decision_ubound, int decision_dim2,
										 double* constraint_lbound, int constraint_dim1,
										 double* constraint_ubound, int constraint_dim2)
{
	// Eigen interfacing to raw buffers
	Eigen::Map<Eigen::VectorXd> full_state_lower_bound(decision_lbound, decision_dim1);
	Eigen::Map<Eigen::VectorXd> full_state_upper_bound(decision_ubound, decision_dim2);
	Eigen::Map<Eigen::VectorXd> full_constraint_lower_bound(constraint_lbound, constraint_dim1);
	Eigen::Map<Eigen::VectorXd> full_constraint_upper_bound(constraint_ubound, constraint_dim2);

	double margin = 0.0;

	full_state_lower_bound.setConstant(0.0);
	full_state_upper_bound.setConstant(1.0);

	full_state_lower_bound(colonne) = -2e19;
	full_state_upper_bound(colonne) = 1.0;

	for(unsigned int j=0; j<righe; j++){
		full_constraint_upper_bound(j) = w_gi(j);
		full_constraint_lower_bound(j) = w_gi(j);
	}
	//constraint on sum of lamdas and s to 1
	full_constraint_upper_bound(righe) = 1.0;
	full_constraint_lower_bound(righe) = 1.0;

}

//void OptWeights::evaluateConstraints(Eigen::Ref<Eigen::VectorXd> full_constraint,
//                                              const Eigen::Ref<const Eigen::VectorXd>& decision_var)
	void OptWeights::evaluateConstraints(double* g, int m, const double* x, int n)
{
	Eigen::Map<Eigen::VectorXd> full_constraint(g, m); 
	const Eigen::Map<const Eigen::VectorXd> decision_var(x, n);
//	for(unsigned int j=0; j<righe; j++){
//		double tmp_constr = 0.0;
//
//		for(unsigned int i=0; i<colonne; i++){
//			 tmp_constr += A_v_description(j,i)*decision_var(i);
//		}
//
//		full_constraint(j) = tmp_constr;
//	}
    //more compact
    full_constraint.segment(0,righe) = A_v_description*decision_var.segment(0,colonne);
	full_constraint(righe) = decision_var.sum();
//	prt(full_constraint)
//	prt(decision_var)
//	std::cout<<"Decision variables in the constraint func"<<std::endl;
//	std::cout<<decision_var<<std::endl;
//	prt(w_gi)
//	prt(w_gi)
}


void OptWeights::evaluateCosts(double& cost,
									  const double* decision, int decision_dim)
{
	// Eigen interfacing to raw buffers
	const Eigen::Map<const Eigen::VectorXd> decision_var(decision, decision_dim);

	Eigen::VectorXd tmp_cost1;
	tmp_cost1.resize(colonne);
	tmp_cost1.setZero();

//	for (int i = 0; i <colonne; i++){
//		tmp_cost1(i) += (decision_var(i)- mean_coeff)*(decision_var(i)- mean_coeff);
//	}
//	std::cout<<"Decision variables"<<std::endl;
//	std::cout<<decision_var<<std::endl;
	cost = -decision_var(colonne);//tmp_cost1.sum();

}

dwl::WholeBodyTrajectory& OptWeights::evaluateSolution(const Eigen::Ref<const Eigen::VectorXd>& solution)
{
    std::cout << "=====================> the optimal margin is: " << solution(colonne) << std::endl;
//    std::cout << "the optimal weights are: " << solution.segment(0,colonne) << std::endl;

//    return motion_solution_;
}

void OptWeights::setInitialState(const Eigen::MatrixXd & A_v, const Eigen::VectorXd & wr)
{
A_v_description.resize(A_v.rows(),A_v.cols());
A_v_description = A_v;
w_gi.resize(A_v.rows());
w_gi.setZero();
w_gi = wr;
colonne = A_v_description.cols();
righe = A_v_description.rows();
mean_coeff = 1.0/((double)colonne);
std::cout<<"set initial state:"<<std::endl;
prt(w_gi.transpose())
prt(A_v_description.transpose())
}
