#include <dwl/model/OptimizationModel.h>


class OptimizationProblem : public dwl::model::OptimizationModel
{

public:
    /** @brief Constructor function */
    OptimizationProblem();

    /** @brief Destructor function */
    ~OptimizationProblem();

    /** @brief Initializes the optimization model, i.e. the dimensions
     * of the optimization vectors */
    void init(bool soft_constraint);
    /**
     * @brief Gets the starting point of the problem
     * @param Eigen::Ref<Eigen::VectorXd> Full initial point
     */
    void getStartingPoint(Eigen::Ref<Eigen::VectorXd> full_initial_point);

    /**
     * @brief Evaluates the bounds of the problem
     * @param Eigen::Ref<Eigen::VectorXd> Full state lower bound
     * @param Eigen::Ref<Eigen::VectorXd> Full state upper bound
     * @param Eigen::Ref<Eigen::VectorXd> Full constraint lower bound
     * @param Eigen::Ref<Eigen::VectorXd> Full constraint upper bound
     */
    void evaluateBounds(Eigen::Ref<Eigen::VectorXd> full_state_lower_bound,
                        Eigen::Ref<Eigen::VectorXd> full_state_upper_bound,
                        Eigen::Ref<Eigen::VectorXd> full_constraint_lower_bound,
                        Eigen::Ref<Eigen::VectorXd> full_constraint_upper_bound);

    /**
     * @brief Evaluates the constraint function given a current decision state
     * @param Eigen::Ref<Eigen::VectorXd> Constraint vector
     * @param const Eigen::Ref<const Eigen:VectorXd>& Decision vector
     */
    void evaluateConstraints(Eigen::Ref<Eigen::VectorXd> full_constraint,
                             const Eigen::Ref<const Eigen::VectorXd>& decision_var);

    /**
     * @brief Evaluates the cost function given a current decision state
     * @param double& Cost value
     * @param const Eigen::Ref<const Eigen:VectorXd>& Decision vector
     */
    void evaluateCosts(double& cost,
                       const Eigen::Ref<const Eigen::VectorXd>& decision_var);

    /**
     * @brief Evaluates the solution from an optimizer
     * @param const Eigen::Ref<const Eigen::VectorXd>& Solution vector
     * @return WholeBodyTrajectory& Returns the whole-body trajectory solution
     */
    dwl::WholeBodyTrajectory& evaluateSolution(const Eigen::Ref<const Eigen::VectorXd>& solution);

};

