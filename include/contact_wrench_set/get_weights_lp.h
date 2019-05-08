 /*  Created on: October 12, 2017
 *      Author: Romeo Orsolino
 */
#include <dwl/model/OptimizationModel.h>
#include <iit/rbd/rbd.h>
#include <iit/rbd/utils.h>
#include <iit/commons/dog/leg_data_map.h>
#include <dwl/WholeBodyState.h>
#include <dwl/model/WholeBodyKinematics.h>
#include <dwl/model/WholeBodyDynamics.h>
#include <ros/ros.h>
#include <string>

class OptWeights : public dwl::model::OptimizationModel
{

public:
    /** @brief Constructor function */
    OptWeights();

    /** @brief Destructor function */
    ~OptWeights();

    /** @brief Initializes the optimization model, i.e. the dimensions
     * of the optimization vectors */
    void init(bool soft_constraint);
    /**
     * @brief Gets the starting point of the problem
     * @param Eigen::Ref<Eigen::VectorXd> Full initial point
     */
    //    void getStartingPoint(Eigen::Ref<Eigen::VectorXd> full_initial_point);

    void getStartingPoint(double* decision, int decision_dim);

    /**
     * @brief Evaluates the bounds of the problem
     * @param Eigen::Ref<Eigen::VectorXd> Full state lower bound
     * @param Eigen::Ref<Eigen::VectorXd> Full state upper bound
     * @param Eigen::Ref<Eigen::VectorXd> Full constraint lower bound
     * @param Eigen::Ref<Eigen::VectorXd> Full constraint upper bound
     */

    //void evaluateBounds(Eigen::Ref<Eigen::VectorXd> full_state_lower_bound,
    //                    Eigen::Ref<Eigen::VectorXd> full_state_upper_bound,
    //                    Eigen::Ref<Eigen::VectorXd> full_constraint_lower_bound,
    //                    Eigen::Ref<Eigen::VectorXd> full_constraint_upper_bound);

    void evaluateBounds(double* decision_lbound, int decision_dim1,
                            double* decision_ubound, int decision_dim2,
                            double* constraint_lbound, int constraint_dim1,
                            double* constraint_ubound, int constraint_dim2);

    /**
     * @brief Evaluates the constraint function given a current decision state
     * @param Eigen::Ref<Eigen::VectorXd> Constraint vector
     * @param const Eigen::Ref<const Eigen:VectorXd>& Decision vector
     */
    //void evaluateConstraints(Eigen::Ref<Eigen::VectorXd> full_constraint,
    //                         const Eigen::Ref<const Eigen::VectorXd>& decision_var);
    void evaluateConstraints(double* g, int m, const double* x, int n);

    /**
     * @brief Evaluates the jacobian of the constraint function given a current decision
     * state
     * @param Eigen::MatrixXd& Jacobian of the constraint function
     * @param const Eigen::VectorXd& Decision vector
     */
//    void evaluateConstraintJacobian(Eigen::MatrixXd& jacobian,
//                                                   const Eigen::VectorXd& decision_var);

    /**
     * @brief Evaluates the cost function given a current decision state
     * @param double& Cost value
     * @param const Eigen::Ref<const Eigen:VectorXd>& Decision vector
     */
    //void evaluateCosts(double& cost,
    //                   const Eigen::Ref<const Eigen::VectorXd>& decision_var);

    void evaluateCosts(double& cost, const double* decision, int decision_dim);


//    void evaluateCostGradient(Eigen::MatrixXd& gradient,
//                                              const Eigen::VectorXd& decision_var);
    /**
     * @brief Evaluates the solution from an optimizer
     * @param const Eigen::Ref<const Eigen::VectorXd>& Solution vector
     * @return WholeBodyTrajectory& Returns the whole-body trajectory solution
     */
    dwl::WholeBodyTrajectory& evaluateSolution(const Eigen::Ref<const Eigen::VectorXd>& solution);

    void setInitialState(const Eigen::MatrixXd & A_v, const Eigen::VectorXd & wr);


private:

    iit::dog::LegDataMap<Eigen::Vector3d> footPosW;
    iit::dog::LegDataMap<bool> stanceCWS;
    iit::rbd::Matrix66d Ic;
    iit::rbd::Vector6D basePoseW0;

    bool weights_as_decision_var = true;
    bool print_all = false;
    Eigen::MatrixXd A_hs, A_v_description, reduced_A_v_descr, A_hs_slack, pseudoInverse;
	Eigen::VectorXd weights;
	iit::rbd::Matrix66d base_inertia_mat;
	Eigen::VectorXd w_gi; // Tau_x, Tau_y, Tau_z, Fx, Fy, Fz
	Eigen::Vector3d com_position; //com_x, com_z
	Eigen::VectorXd biasTerm;
	double mean_coeff, des_height;
    int righe = 0;
    int colonne = 0;
    iit::rbd::ForceVector gravitoInertialWrench;
    Eigen::Vector3d comPosWF, comAccWF, old_des_base_xd, initial_com_pos, final_heuristics_com_pos;
	unsigned int nodes_num = 3;
	unsigned int states_per_knot = 1;
	unsigned int constraints_per_knot = 1;
	double horizon_time = 1.0;
	double sampling_time = horizon_time/(double)nodes_num;

    /** @brief Whole-body current state */
    dwl::WholeBodyState current_wbs;
    /** @brief Kinematic model */
    dwl::model::WholeBodyKinematics wb_kin;
    /** @brief Dynamic model */
    dwl::model::WholeBodyDynamics wb_dyn;
    /** Floating-base system */
    dwl::model::FloatingBaseSystem fbs;

	/** @brief Joint and foot id mapping between urdf to robcogen */
	std::vector<unsigned int> joint_id_map_;
	std::vector<unsigned int> foot_id_map_;


};

