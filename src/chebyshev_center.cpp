 /*  Created on: October 12, 2017
 *      Author: Romeo Orsolino
 */
#include <contact_wrench_set/chebyshev_center.h>
#include <math_utils/utils.h>

using namespace Eigen;
using namespace iit;

ChebyshevCenter::ChebyshevCenter() {}

ChebyshevCenter::~ChebyshevCenter() {}

void ChebyshevCenter::init(bool soft_constraint)
{

	// Getting the dimension of the decision vector
	state_dimension_ = 1; // the number of weights  plus the margin
	std::cout<<"state dimension: "<<state_dimension_<<std::endl;

	constraint_dimension_ = rows;
	std::cout<<"constraint dimension: "<<constraint_dimension_<<std::endl;


}


void ChebyshevCenter::getStartingPoint(Eigen::Ref<Eigen::VectorXd> full_initial_point)
{

	full_initial_point.setZero();

}

void ChebyshevCenter::evaluateBounds(Eigen::Ref<Eigen::VectorXd> full_state_lower_bound,
		Eigen::Ref<Eigen::VectorXd> full_state_upper_bound,
		Eigen::Ref<Eigen::VectorXd> full_constraint_lower_bound,
		Eigen::Ref<Eigen::VectorXd> full_constraint_upper_bound)
{

	full_state_lower_bound.setConstant(-2e19);
	full_state_upper_bound.setConstant(2e19);

	full_constraint_upper_bound = b;
	full_constraint_lower_bound.setConstant(-2e19);

}

void ChebyshevCenter::evaluateConstraints(Eigen::Ref<Eigen::VectorXd> full_constraint,
                                              const Eigen::Ref<const Eigen::VectorXd>& decision_var)
{

    full_constraint = A*wrench_gi + decision_var(0)*A_norm;
//    prt(full_constraint)
//    prt(decision_var)

}


void ChebyshevCenter::evaluateCosts(double& cost,
                                        const Eigen::Ref<const Eigen::VectorXd>& decision_var)
{

	cost = - decision_var(0);

}

void ChebyshevCenter::evaluateSolution(const Eigen::Ref<const Eigen::VectorXd>& solution)
{
    std::cout << "the optimal margin is: " << solution(cols) << std::endl;
//    std::cout << "the optimal weights are: " << solution.segment(0,colonne) << std::endl;

}

void ChebyshevCenter::setInitialState(Eigen::MatrixXd A_hs_description, Eigen::VectorXd wrench)
{

	cols = A_hs_description.cols();
	rows = A_hs_description.rows();

	A_hs.resize(rows,cols);
	A_hs = A_hs_description;

	prt(A_hs)

	A.resize(rows, cols-1);
	A = A_hs.block(0,1,rows, cols-1);

	b.resize(rows);
	b = A_hs.block(0,0,rows,1);

	A_norm.resize(rows);
	for(int k=0; k<rows; k++){
		A_norm(k) = A.block(k,0,1,cols-1).norm();
	}

	wrench_gi.resize(cols-1);
	wrench_gi = wrench;

	prt(A)
	prt(b)
	prt(wrench_gi)

}



