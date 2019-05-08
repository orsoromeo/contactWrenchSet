/*
 * CCWS_margin.h
 *
 *  Created on: April 12, 2017
 *      Author: Romeo Orsolino
 */

#ifndef FeasibleWrenchPolytope_API_H_
#define FeasibleWrenchPolytope_API_H_

#include <Eigen/Dense>
#include <math.h>
#include <ros/package.h>
#include <iit/rbd/rbd.h>
#include <iit/rbd/utils.h>
#include <iit/commons/dog/leg_data_map.h>
#include <iit/commons/dog/joint_id_tricks.h>
#include <rbdl/rbdl.h>
#include <stdio.h>
#include <string.h>

//DWL
#include <dwl/WholeBodyState.h>
#include <dwl/model/WholeBodyKinematics.h>
#include <dwl/model/WholeBodyDynamics.h>
#include <dwl/utils/RigidBodyDynamics.h>
#include <dwl/utils/Orientation.h>
#include <dwl/solver/OptimizationSolver.h>
#include <dwl/solver/IpoptNLP.h>
#include <dwl_rviz_plugin/DisplayInterface.h>

//CDD
//#include "setoper.h"
//#include "cdd_f.h"

//FWP class
#include <contact_wrench_set/mink_glpk.h>
#include <contact_wrench_set/get_weights_lp.h>
#include <contact_wrench_set/chebyshev_center.h>
#include <contact_wrench_set/testOptim.h>

// Politopix
#include "../thirdparty/politopix/trunk/politopixAPI.h"
#include "IpIpoptApplication.hpp"


#include <cassert>
#include <iostream>
#include <cstdlib>
#include <random>
#include <thread>
#include <mutex>              // std::mutex, std::unique_lock
#include <condition_variable> // std::condition_variable
#include <realtime_tools/realtime_buffer.h>

//for config files // to read in config files - include before SL includes to avoid name collisions
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using namespace iit;

class FeasibleWrenchPolytope_API {
public:
	typedef Eigen::Matrix<double, 4,1> Vector4D;
public:
	FeasibleWrenchPolytope_API();
	virtual ~FeasibleWrenchPolytope_API();

	void init();

	struct CWSData {
		iit::dog::LegDataMap<rbd::Vector3d> stance_feet_WF;
		iit::dog::LegDataMap<bool> stance_legs_flag;
		iit::dog::LegDataMap<rbd::Vector3d> normal;
		iit::dog::LegDataMap<double> friction;
		iit::dog::LegDataMap<double> max_normal_force;
		double robot_mass = 92.0;
	};

	struct TLSData {
		iit::dog::LegDataMap<rbd::Vector3d> stance_feet_WF;
		iit::dog::LegDataMap<rbd::Vector3d> max_torque;
		iit::dog::LegDataMap<rbd::Vector3d> legs_grav_torque;
		Eigen::Vector3d comAWPApproxWF;
		Eigen::MatrixXd fixed_base_jac;
	};

	enum constraint{no_constraint = 0, only_friction = 1, only_actuation = 2, friction_and_actuation = 3};

	enum wrenchType{full_6D = 0 , angular_3D, linear_3D};

	struct FWSData	{
		iit::dog::LegDataMap<rbd::Vector3d> stance_feet_WF;
		iit::dog::LegDataMap<bool> stance_legs_flag;
		iit::dog::LegDataMap<rbd::Vector3d> normal;
		iit::dog::LegDataMap<double> friction;
		iit::dog::LegDataMap<double> max_normal_force;
		iit::dog::LegDataMap<rbd::Vector3d> max_torque;
		iit::dog::LegDataMap<rbd::Vector3d> legs_grav_torque;
		// FWP optimization, default = 0 (no optimization), CWC constraints only = 1,  AWP constraints only = 2, full FWP constraints = 3
//		int constraint_type = 3;
		constraint constraint_type;
		wrenchType wrench_type;
		bool optimize_time;
	};

	typedef Eigen::MatrixXd hs_description;
	typedef Eigen::MatrixXd v_description;
	typedef Eigen::MatrixXi topology_map;
	typedef Eigen::VectorXi topology_vec;

	struct polytope {
		v_description v_rep;
		hs_description hs_rep;
		topology_map top_map;
		topology_vec top_vec;
	};

	void ZMP(Eigen::Vector3d r,
	                        Eigen::Vector3d r_dd,
						   Eigen::Vector3d & ZMP);

	void ZMPstability(iit::dog::LegDataMap<rbd::Vector3d> stance_feet_WF,
							    iit::dog::LegDataMap<bool> stance_legs_flag,
								Eigen::Vector3d r,
	                            Eigen::Vector3d r_dd,
								Eigen::Vector3d & ZMP,
								double & ZMPdistance);

	void distPointToLine(Eigen::Vector3d P1,
			Eigen::Vector3d P2,
			Eigen::Vector3d X,
			double & dist);

	void IntertialWrench_BF(Eigen::Vector3d r_d,
		                        Eigen::Vector3d r_dd,
		                        Eigen::Vector3d orient_d,
		                        Eigen::Vector3d orient_dd,
							   iit::rbd::Matrix66d InertiaMatrix,
		                        iit::rbd::ForceVector & hdot_BF);

	void IntertialWrench_WF(Eigen::Vector3d r,
		                        Eigen::Vector3d r_d,
		                        Eigen::Vector3d r_dd,
		                        Eigen::Vector3d orient,
		                        Eigen::Vector3d orient_d,
		                        Eigen::Vector3d orient_dd,
							   iit::rbd::Matrix66d InertiaMatrix,
							   iit::rbd::ForceVector & hdot_WF);

	void GravitationalWrench_WF(Eigen::Vector3d r, CWSData cws_struct);

	void GravitationalWrench_WF(Eigen::Vector3d r, CWSData cws_struct, Eigen::VectorXd & grav_wrench_W);

//	void compute_cddlib(iit::dog::LegDataMap<rbd::Vector3d> stance_feet_WF,
//        			iit::dog::LegDataMap<bool> stance_legs_flag,
//					iit::dog::LegDataMap<double> friction_coeffs,
//					iit::dog::LegDataMap<rbd::Vector3d> terrain_normals,
//					Eigen::MatrixXd & A_hs_description);
//
//	void compute_cddlib(CWSData cws_struct, Eigen::MatrixXd & A_hs_description);

	void eigen_hs_description(const boost::shared_ptr<Polytope_Rn> primal_polytope, Eigen::MatrixXd & A_hs_description);

	void eigen_v_description(const boost::shared_ptr<Polytope_Rn> primal_polytope, Eigen::MatrixXd & A_v_description);

	void contact_wrench_set(const CWSData cws_struct, Eigen::MatrixXd & A_v_description);

	void get_hs_description_politopix(Eigen::MatrixXd & A_v_description,
										boost::shared_ptr<Polytope_Rn> & primal_polytope);

	void get_hs_description_politopix(Eigen::MatrixXd & A_v_description,
										Eigen::MatrixXd & A_hs_description);

	void force_polygon_analytic_hs_rep(const v_description v_rep, hs_description & hs_rep);

	void normalize_half_spaces(Eigen::MatrixXd & A_hs);

	bool compute_intersection(FeasibleWrenchPolytope_API::v_description friction_cone,
											FeasibleWrenchPolytope_API::v_description force_polytope,
											const int leg_id,
											FeasibleWrenchPolytope_API::v_description & intersection);

	void CreateEdges(const iit::dog::LegDataMap<rbd::Vector3d> stance_feet_WF,
									const iit::dog::LegDataMap<double> friction_coeff,
									const iit::dog::LegDataMap<rbd::Vector3d> terrain_normal,
									const iit::dog::LegDataMap<double> max_normal_grf,
									iit::dog::LegDataMap< Eigen::MatrixXd > & legs_cone);

	void FeasibilityMargin(const FeasibleWrenchPolytope_API::wrenchType wrench_type,
							const Eigen::MatrixXd M,
							const Eigen::VectorXd wrench_gi,
							double & cws_feasibility);

	void hs_based_margin(const FeasibleWrenchPolytope_API::wrenchType Wrench_type,
										const Eigen::MatrixXd M1,
										const Eigen::VectorXd wrench_gi,
										double & fwp_hs_margin);

	void NormalToEdgesRotation(const double half_cone_angle,
			const unsigned int edges_per_cone,
			iit::dog::LegDataMap<Eigen::Matrix3d> & R);

	double GetVectorLenght(Eigen::Vector3d vec);

	void remove_repeated_vertices(boost::shared_ptr<Polytope_Rn> primal_polytope);

	void remove_hs_with_null_norm(const Eigen::MatrixXd A_hs, Eigen::MatrixXd & A_hs_reduced);

	void fill_angular_component(const Eigen::Vector3d stance_feet_pos,
			const Eigen::MatrixXd polytope_3d,
			Eigen::MatrixXd & polytope_6d);

	void legs_torque_limits_3d(	const FeasibleWrenchPolytope_API::TLSData tls,
								iit::dog::LegDataMap< polytope > & max_lin_grf);

	void legs_friction_cone_3d(const iit::dog::LegDataMap<double> friction_coeff,
			const iit::dog::LegDataMap<rbd::Vector3d> terrain_normal,
			const iit::dog::LegDataMap<double> max_normal_grf,
			iit::dog::LegDataMap< Eigen::MatrixXd > & legs_set);

	bool feasibility_wrench_set(const TLSData tlws_struct,
								const CWSData cws_struct,
								const FWSData fwp_options,
								Eigen::MatrixXd & FWS_v_description);

	bool bounded_friction_polytopes(const TLSData tlws_struct,
			const CWSData cws_struct,
			const FWSData fwp_options,
			iit::dog::LegDataMap<Eigen::MatrixXd > & bounded_fp);

	void vertex_based_margin_LP(const wrenchType & wrench_type,
														const v_description & A_v, 
														const Eigen::VectorXd & wrench_gi, 
														double & fwp_margin);

	void residual_radius_LP(const hs_description A_hs, const Eigen::VectorXd wrench_gi, double & residual_radius);

	void get_jacobian(const int constraint_type,
						const dog::LegDataMap<Eigen::Vector3d> footPos_BF,
						Eigen::MatrixXd & fixed_base_jacobian);

//	void eigen2politopix(const iit::dog::LegDataMap<rbd::ForceVector> edges,
//			                 boost::shared_ptr<Polytope_Rn> & polytope);

//	ddf_MatrixPtr eigen2cdd(iit::dog::LegDataMap<rbd::ForceVector> edge);
//
//	void cdd2eigen(ddf_MatrixPtr A, Eigen::MatrixXd & A_eig);

//	void cdd2politopix(ddf_MatrixPtr A, boost::shared_ptr<Polytope_Rn> & poly);

//	void debugPolytope(boost::shared_ptr<Polytope_Rn> polytope);

	iit::rbd::ForceVector getGravitoInertialWrench();

	void torque_limits_wrench_set(const TLSData tlws_struct,
			const Eigen::MatrixXd jacobian,
			const Eigen::Vector3d com_orientation,
			Eigen::MatrixXd & TLWS_v_description);

	void joint_space_zonotope(const Eigen::Vector3d joints_lim, Eigen::Matrix<double, 3, 8> & joint_space_zono);

	void get_bounded_friction_polytopes(iit::dog::LegDataMap< polytope > & bfp);

	void get_linear_friction_cones(iit::dog::LegDataMap<polytope > & lfc);

	void get_force_polytopes(iit::dog::LegDataMap<polytope> & fp);

	void get_fwp(polytope & fw_poly);

	void clock_wise_sort(const FeasibleWrenchPolytope_API::v_description v_rep,
						const FeasibleWrenchPolytope_API::topology_vec top_vec,
						FeasibleWrenchPolytope_API::topology_map & top_map);

	void build_topology(const boost::shared_ptr<Polytope_Rn> primal_polytope,
						topology_map & topology,
						topology_vec & top_vec);

	void draw_polygon_vertices(iit::dog::LegDataMap<polytope> poly,
						iit::dog::LegDataMap<Eigen::Vector3d > foot_pos,
						dwl::Color color,
						double scaling_factor,
						std::shared_ptr<dwl_rviz_plugin::DisplayInterface> display_);

	void draw_polygon_vertices(polytope poly,
						Eigen::Vector3d application_point,
						dwl::Color color,
						double scaling_factor,
						std::shared_ptr<dwl_rviz_plugin::DisplayInterface> display_);

	void draw_polygon_edges(iit::dog::LegDataMap< polytope > poly,
			iit::dog::LegDataMap<Eigen::Vector3d > foot_pos,
			dwl::Color color,
			double scaling_factor,
			std::shared_ptr<dwl_rviz_plugin::DisplayInterface> display_);

	void draw_friction_cone_edges(iit::dog::LegDataMap< polytope > poly,
			iit::dog::LegDataMap<Eigen::Vector3d > foot_pos,
			dwl::Color color,
			double scaling_factor,
			std::shared_ptr<dwl_rviz_plugin::DisplayInterface> display_);

	void check_friction_cone_topology(boost::shared_ptr<Polytope_Rn> primal_polytope,
			topology_map & topology,
			topology_vec & top_vec);

	void reorder_v_description(const boost::shared_ptr<Polytope_Rn> & primal_polytope,
								Eigen::MatrixXd & A_v_description);

	void set_stance_change(const iit::dog::LegDataMap<bool> previous_stance_legs,
										const iit::dog::LegDataMap<bool> new_stance_legs,
										bool & stance_flag);

	void reset(CWSData & cwc,
			TLSData & awp,
			FWSData & fwp_opt,
			polytope & fwp);

	void remove_average_point(const FWSData fwp_options,
			const v_description & polytope,
			const Eigen::VectorXd & wrench_gi,
			v_description & avg_polytope,
			Eigen::VectorXd & wrench_gi_avg);

	void kill();

	void setPrintAll(bool print_all_flag){
	    this->print_all = print_all_flag;
	}

	Eigen::MatrixXd A_matrix;

	typedef realtime_tools::RealtimeBuffer<Eigen::VectorXd> cws_buffer;
	typedef realtime_tools::RealtimeBuffer<FeasibleWrenchPolytope_API::CWSData> cws_buffer_read;
	typedef std::mutex cws_mtx;
	typedef std::condition_variable cws_condv;

//	void getPointsNum(unsigned int & pointsNumber);

private:
	double g0;
	double mass;

        iit::dog::JointState q;
	double TOL = 10^-3;

//	std::string param_file = ros::package::getPath(std::string("contact_wrench_set")) + "/config/cws_options.ini";
	//config files
//	boost::property_tree::ptree config_;
//	bool print_friction_edges = config_.get<double>("Options.print_friction_cone_edges");
	bool print_friction_edges = false;
	// in case that only a 3D dynamics is used (therefore for use_6d_dynamics = false) you can
	// chose between using the only the linear wrench or only the angular wrench
	bool use_linear_wrench = true;
	bool print_all = false;
	bool use_torque_set;
	bool use_simplified_torque_limits_set = false;

//	Generator_Rn vx = Generator_Rn(cardinality);
	iit::dog::LegDataMap<Eigen::Vector3d> footPos_WF;
	RigidBodyDynamics::Math::SpatialVector vel = RigidBodyDynamics::Math::SpatialVector::Zero();
	RigidBodyDynamics::Math::SpatialVector acc = RigidBodyDynamics::Math::SpatialVector::Zero();
	iit::dog::LegDataMap<double> max_normal_force;
	iit::dog::LegDataMap<double> normal_force_projection;
	iit::dog::LegDataMap<double> half_cone_angle;
	iit::dog::LegDataMap<Eigen::Matrix3d> R;
	iit::dog::LegDataMap< polytope > bounded_friction_poly,force_polytopes;
	iit::dog::LegDataMap<Eigen::MatrixXd > linear_friction_cones;
	iit::dog::LegDataMap<iit::dog::LegDataMap<rbd::Vector3d> > friction_cone_verteces;
	iit::dog::LegDataMap<iit::dog::LegDataMap<rbd::Vector3d> > friction_cone_edges;
	iit::dog::LegDataMap<iit::dog::LegDataMap<rbd::ForceVector> > cws_pyramid;

	Eigen::MatrixXd A_hs;
	topology_map lfc_top_map, bfc_top_map, fp_top_map, fwp_top_map;
	topology_vec lfc_top_vec, bfc_top_vec, fp_top_vec, fwp_top_vec;
	//std::shared_ptr< dwl::solver::OptimizationSolver> solver2;
	dwl::solver::OptimizationSolver* margin_solver = new dwl::solver::IpoptNLP();
    OptWeights opt_weights;
	ChebyshevCenter cheb_center;
	std::shared_ptr< dwl::solver::OptimizationSolver> solver;
	polytope cwc, awp, fwp;

	int stance_num;
	unsigned int cardinality = 6;
	unsigned int pointsNum = 5;
	double dist, num, denum;
	double margin, old_cws_feasibility;
	rbd::ForceVector inertial_wrench_WF, grav_wrench_WF;
//	ddf_rowindex newpos;
//	ddf_rowset impl_lin,redset;

	/** @brief Kinematical model */
	dwl::model::WholeBodyKinematics wb_kin;
	/** @brief Kinematical model */
	dwl::model::WholeBodyDynamics wb_dyn;
	/** @brief Actual whole-body state information */
	dwl::WholeBodyState current_wbs;
	/** Floating-base system */
	dwl::model::FloatingBaseSystem fbs;

	// Resetting the system from the hyq urdf file
	std::string model_;
	string robot_;
//	dwl::rbd::BodyVector3d ik_pos;
//	dwl::rbd::BodySelector feet_names;
//	std::vector<unsigned int> foot_id_map_;
};

#endif /* FeasibleWrenchPolytope_API_H_ */





