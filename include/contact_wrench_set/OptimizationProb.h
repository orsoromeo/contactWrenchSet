#include <dwl/model/OptimizationModel.h>


class OptimizationProb : public dwl::model::OptimizationModel
{

public:
    /** @brief Constructor function */
    OptimizationProb();

    /** @brief Destructor function */
    ~OptimizationProb();

    /** @brief Initializes the optimization model, i.e. the dimensions
     * of the optimization vectors */
    void init(bool soft_constraint);
    /**
     * @brief Gets the starting point of the problem
     * @param double* Initial values for the decision variables, $x$
     * @param int Number of the decision variables
     */
    void getStartingPoint(double* decision, int decision_dim);
    /**
     * @brief Evaluates the bounds of the problem
     * @param double* Lower bounds $x^L$ for $x$
     * @param double* Upper bounds $x^L$ for $x$
     * @param int Number of decision variables (dimension of $x$)
     * @param double* Lower bounds of the constraints $g^L$ for $x$
     * @param double* Upper bounds of the constraints $g^L$ for $x$
     * @param int Number of constraints (dimension of $g(x)$)
     */
    void evaluateBounds(double* decision_lbound, int decision_dim1,
                        double* decision_ubound, int decision_dim2,
                        double* constraint_lbound, int constraint_dim1,
                        double* constraint_ubound, int constraint_dim2);
    /**
     * @brief Evaluates the constraint function given a current decision state
     * @param double* Array of constraint function values, $g(x)$
     * @param int Number of constraint variables (dimension of $g(x)$)
     * @param const double* Array of the decision variables, $x$, at which the constraint functions,
     * $g(x)$, are evaluated
     * @param int Number of decision variables (dimension of $x$)
     */
    void evaluateConstraints(double* constraint, int constraint_dim,
                             const double* decision, int decision_dim);

    /**
     * @brief Evaluates the cost function given a current decision state
     * @param double& Value of the objective function ($f(x)$).
     * @param const double* Array of the decision variables, $x$, at which the cost functions,
     * $f(x)$, is evaluated
     * @param int Number of decision variables (dimension of $x$)
     */
    void evaluateCosts(double& cost,
                           const double* decision, int decision_dim);
     
    dwl::WholeBodyTrajectory& evaluateSolution(const Eigen::Ref<const Eigen::VectorXd>& solution);

};

