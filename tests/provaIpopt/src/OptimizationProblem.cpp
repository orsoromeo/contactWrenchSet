#include <OptimizationProblem.h>



OptimizationProblem::OptimizationProblem() {}


OptimizationProblem::~OptimizationProblem() {}


void OptimizationProblem::init(bool soft_constraint)
{
        // Getting the dimension of the decision vector
//        state_dimension_ = 4;
//        constraint_dimension_ = 2;
        state_dimension_ = 2;
        constraint_dimension_ = 1;
}


void OptimizationProblem::getStartingPoint(Eigen::Ref<Eigen::VectorXd> full_initial_point)
{
//    full_initial_point(0) = 1.0;
//    full_initial_point(1) = 5.0;
//    full_initial_point(2) = 5.0;
//    full_initial_point(3) = 1.0;

    full_initial_point(0) = 4.0;
    full_initial_point(1) = 4.0;
}

void OptimizationProblem::evaluateBounds(Eigen::Ref<Eigen::VectorXd> full_state_lower_bound,
                                         Eigen::Ref<Eigen::VectorXd> full_state_upper_bound,
                                         Eigen::Ref<Eigen::VectorXd> full_constraint_lower_bound,
                                         Eigen::Ref<Eigen::VectorXd> full_constraint_upper_bound)
{

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
    full_state_lower_bound(1) = -2e19;
    full_state_upper_bound(0) = 2e19;
    full_state_upper_bound(1) = 2e19;

    full_constraint_lower_bound(0) = 0.0;
//    full_constraint_lower_bound(1) = 40.0;
   full_constraint_upper_bound(0) = 0.0;
//    full_constraint_upper_bound(1) = 40.0;


}

void OptimizationProblem::evaluateConstraints(Eigen::Ref<Eigen::VectorXd> full_constraint,
                                              const Eigen::Ref<const Eigen::VectorXd>& decision_var)
{
//    full_constraint = Eigen::VectorXd::Zero(constraint_dimension_);
//    full_constraint(0) = decision_var(0) * decision_var(1) *
//            decision_var(2) * decision_var(3);
//    full_constraint(1) = decision_var(0) * decision_var(0) +
//            decision_var(1) * decision_var(1) + decision_var(2) * decision_var(2) +
//            decision_var(3) * decision_var(3);
    full_constraint(0) = decision_var(0) +  decision_var(1) ;


}


void OptimizationProblem::evaluateCosts(double& cost,
                                        const Eigen::Ref<const Eigen::VectorXd>& decision_var)
{
//    cost = decision_var(0) * decision_var(3) * (decision_var(0) +
//            decision_var(1) + decision_var(2)) + decision_var(2);

    cost = (decision_var - Eigen::Vector2d(1,2)).norm();
}

dwl::WholeBodyTrajectory& OptimizationProblem::evaluateSolution(const Eigen::Ref<const Eigen::VectorXd>& solution)
{
    std::cout << "the optimization solution is " << solution.transpose() << std::endl;

    //return motion_solution_;
}

