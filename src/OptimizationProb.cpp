#include<contact_wrench_set/OptimizationProb.h>

OptimizationProb::OptimizationProb() {}


OptimizationProb::~OptimizationProb() {}


void OptimizationProb::init(bool soft_constraint)
{
        // Getting the dimension of the decision vector
//        state_dimension_ = 4;
//        constraint_dimension_ = 2;
        state_dimension_ = 2;
        constraint_dimension_ = 2;
}


void OptimizationProb::getStartingPoint(double* decision, int decision_dim)
{
  // Eigen interfacing to raw buffers
  Eigen::Map<Eigen::VectorXd> full_initial_point(decision, decision_dim);

//    full_initial_point(0) = 1.0;
//    full_initial_point(1) = 5.0;
//    full_initial_point(2) = 5.0;
//    full_initial_point(3) = 1.0;

    full_initial_point(0) = 2.0;
    full_initial_point(1) = 4.0;
}

void OptimizationProb::evaluateBounds(double* decision_lbound, int decision_dim1,
                     double* decision_ubound, int decision_dim2,
                     double* constraint_lbound, int constraint_dim1,
                     double* constraint_ubound, int constraint_dim2)
{
  // Eigen interfacing to raw buffers
  Eigen::Map<Eigen::VectorXd> full_state_lower_bound(decision_lbound, decision_dim1);
  Eigen::Map<Eigen::VectorXd> full_state_upper_bound(decision_ubound, decision_dim2);
  Eigen::Map<Eigen::VectorXd> full_constraint_lower_bound(constraint_lbound, constraint_dim1);
  Eigen::Map<Eigen::VectorXd> full_constraint_upper_bound(constraint_ubound, constraint_dim2);
    //bounds on state should be always defined (if they are not set they are set to 0,0)
//    full_state_lower_bound(0) = 1.0;
//    full_state_lower_bound(1) = 1.0;
//    full_state_lower_bound(2) = 1.0;
//    full_state_lower_bound(3) = 1.0;
//    full_state_upper_bound(0) = 5.0;
//    full_state_upper_bound(1) = 5.0;
//    full_state_upper_bound(2) = 5.0;
//    full_state_upper_bound(3) = 5.0;

    full_state_lower_bound(0) = -2e19;
    full_state_upper_bound(0) = 2e19;

    full_state_lower_bound(1) = -2e19;
    full_state_upper_bound(1) = 2e19;

    full_constraint_lower_bound(0) = -2e19;
    full_constraint_upper_bound(0) = 0.0;

    full_constraint_lower_bound(1) = 0.0;
    full_constraint_upper_bound(1) = 2e19;

}

void OptimizationProb::evaluateConstraints(double* constraint, int constraint_dim,
                        const double* decision, int decision_dim)
{
    Eigen::Map<Eigen::VectorXd> full_constraint(constraint, constraint_dim); 
    const Eigen::Map<const Eigen::VectorXd> decision_var(decision, decision_dim);
//    full_constraint = Eigen::VectorXd::Zero(constraint_dimension_);
//    full_constraint(0) = decision_var(0) * decision_var(1) *
//            decision_var(2) * decision_var(3);
//    full_constraint(1) = decision_var(0) * decision_var(0) +
//            decision_var(1) * decision_var(1) + decision_var(2) * decision_var(2) +
//            decision_var(3) * decision_var(3);
    full_constraint(0) = decision_var(0) -  decision_var(1);
    full_constraint(1) = decision_var(0) +  decision_var(1) - 2.0;


}


void OptimizationProb::evaluateCosts(double& cost,
                      const double* decision, int decision_dim)
{
  // Eigen interfacing to raw buffers
  const Eigen::Map<const Eigen::VectorXd> decision_var(decision, decision_dim);
//    cost = decision_var(0) * decision_var(3) * (decision_var(0) +
//            decision_var(1) + decision_var(2)) + decision_var(2);

//    cost = (decision_var - Eigen::Vector2d(1,2)).norm();
    cost = decision_var.norm();
}

dwl::WholeBodyTrajectory& OptimizationProb::evaluateSolution(const Eigen::Ref<const Eigen::VectorXd>& solution)
{
    std::cout << "the optimization solution is " << solution.transpose() << std::endl;

    //return motion_solution_;
}

