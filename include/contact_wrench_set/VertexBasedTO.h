#include <dwl/model/OptimizationModel.h>
#include <iit/commons/dog/leg_data_map.h>
#include <contact_wrench_set/FeasibleWrenchPolytope_API.h>
#include <dwl/WholeBodyState.h>
#include <dwl/model/WholeBodyKinematics.h>
#include <dwl/model/WholeBodyDynamics.h>
#include <ros/ros.h>
#include <string>

class VertexBasedTO : public dwl::model::OptimizationModel
{

public:
    /** @brief Constructor function */
    VertexBasedTO();

    /** @brief Destructor function */
    ~VertexBasedTO();

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
    void evaluateCosts(double& cost,
                       const Eigen::Ref<const Eigen::VectorXd>& decision_var);

    void checkPoint(Eigen::Vector3d & inputPoint);

    bool optimize_next_com_pos(const Eigen::Vector3d final_com_pos,
    							const Eigen::Vector3d actual_com_pos,
    							const Eigen::Vector3d actual_com_orient,
								const iit::dog::LegDataMap<double> mu,
								const iit::dog::LegDataMap<Eigen::Vector3d> vec_incl,
								const iit::dog::LegDataMap<Eigen::Vector3d> footPosBF,
								const iit::dog::LegID swing_leg_index,
								const iit::rbd::Matrix66d inertia_mat,
								const double heuristic_move_base,
								const bool use_weights_as_dv,
								const unsigned int nodes_n,
								const FeasibleWrenchPolytope_API::CWSData cwc,
								const FeasibleWrenchPolytope_API::TLSData tls,
								const FeasibleWrenchPolytope_API::FWSData fwp_option,
								Eigen::VectorXd & opt_com_traj);

//    void evaluateCostGradient(Eigen::MatrixXd& gradient,
//                                              const Eigen::VectorXd& decision_var);
    /**
     * @brief Evaluates the solution from an optimizer
     * @param const Eigen::Ref<const Eigen::VectorXd>& Solution vector
     * @return WholeBodyTrajectory& Returns the whole-body trajectory solution
     */
    void evaluateSolution(const Eigen::Ref<const Eigen::VectorXd>& solution);

    bool setInitialState(const Eigen::Vector3d base_posW_inital_pos,
			const Eigen::Vector3d base_posW_heuristic_guess,
            const Eigen::Vector3d base_orient,
            const iit::dog::LegDataMap<double> mu,
            const iit::dog::LegDataMap<Eigen::Vector3d> normals,
			const iit::dog::LegDataMap<double> max_normal_grf,
            const iit::dog::LegDataMap<Eigen::Vector3d> footPos,
            const iit::dog::LegID swing_index,
			const iit::rbd::Matrix66d inertia_mat,
			const bool use_weights_as_dv,
			const unsigned int nodes_n,
			const  FeasibleWrenchPolytope_API::CWSData cwc,
			const FeasibleWrenchPolytope_API::TLSData tls,
			const FeasibleWrenchPolytope_API::FWSData fwp_option);

    void setPrintAll(bool print_all_flag)
    {
        print_all = print_all_flag;
    }

    bool getPrintAllState(){
        return print_all;;
    }

private:

    iit::dog::LegDataMap<Eigen::Vector3d> footPosW;
    iit::dog::LegDataMap<bool> stanceCWS;
    iit::rbd::Matrix66d Ic;
    iit::rbd::Vector6D basePoseW0;
    FeasibleWrenchPolytope_API::FWSData fws_data;
    FeasibleWrenchPolytope_API::CWSData cws_data;
    FeasibleWrenchPolytope_API::TLSData tls_data;
//    FeasibleWrenchPolytope_API::FWSData fwp_opt;
    FeasibleWrenchPolytope_API cws_margin;

    bool weights_as_decision_var;
    bool print_all = false;
    bool use_pseudo_inverse = true;
    Eigen::MatrixXd A_hs, A_v_description, A_v_description_avg, reduced_A_v_descr, A_hs_slack, pseudoInverse, pseudoInverse_avg;
	Eigen::VectorXd weights, vertex_avg;
	iit::rbd::Matrix66d base_inertia_mat;
	Eigen::VectorXd wrench; // Tau_x, Tau_y, Tau_z, Fx, Fy, Fz
	Eigen::Vector3d com_position; //com_x, com_z
	Eigen::VectorXd biasTerm, biasTerm_avg;
	double mean_coeff, des_height, margin;
    int nRows = 0;
    int nCols = 0;
    iit::rbd::ForceVector gravitoInertialWrench;
    Eigen::Vector3d comPosWF, comAccWF, old_des_base_xd, initial_com_pos, final_heuristics_com_pos;
    double heuristic_move_base_duration = 0.5;
    double min_move_base_duration = 0.0;
    double max_move_base_duration = 0.0;
    unsigned int nodes_num;
    bool optimize_time_duration = false;
    double horizon_time = 0.5;
    unsigned int states_per_knot = 1;
    unsigned int constraints_per_knot = 1;
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
	std::vector<unsigned int> foot_id_map_;


};

