#include <contact_wrench_set/VertexBasedTO.h>
//#include <math_utils/utils.h>


using namespace Eigen;
using namespace iit;

VertexBasedTO::VertexBasedTO() {}


VertexBasedTO::~VertexBasedTO() {}


void VertexBasedTO::init(bool soft_constraint)
{
    std::cout<<"started init function"<<std::endl;
	if (weights_as_decision_var){
        // Getting the dimension of the decision vector
		states_per_knot = com_position.size() + nCols;
        state_dimension_ = nodes_num * states_per_knot; // number of dimensions of the COM position + the number of weights
        std::cout<<"state dimension: "<<state_dimension_<<std::endl;
        constraints_per_knot = 1 + 3;
        constraint_dimension_ = nodes_num * constraints_per_knot;
        std::cout<<"constraint dimension: "<<constraint_dimension_<<std::endl;
	}else{
        // Getting the dimension of the decision vector
		if(optimize_time_duration){
	        state_dimension_ = nodes_num * 4 + 1; // com_x, com_x_dot, com_y, com_y_dot for each time know plus the total time horizon
		}else{
	        state_dimension_ = nodes_num * 4; // com_x, com_x_dot, com_y, com_y_dot for each time know plus the total time horizon
		}
        constraint_dimension_ = nodes_num * 2 + 4;
        if (print_all){
            std::cout<<"state dimension: "<<state_dimension_<<std::endl;
        	std::cout<<"constraint dimension: "<<constraint_dimension_<<std::endl;
        }
	}
	cws_margin.setPrintAll(print_all);
}


void VertexBasedTO::getStartingPoint(Eigen::Ref<Eigen::VectorXd> full_initial_point)
{

	if (weights_as_decision_var){
		for(unsigned int j=0; j<nodes_num; j++){
			(full_initial_point.segment(j*states_per_knot,com_position.size()))= final_heuristics_com_pos;
//			full_initial_point(2) = des_height;
			(full_initial_point.segment(j*states_per_knot + com_position.size(),nCols)).setConstant(mean_coeff);
		}
	}else{
		Eigen::Vector3d mean_coeff_vec;
		mean_coeff_vec.resize(nCols-1);
		mean_coeff_vec.setConstant(mean_coeff);
		Eigen::VectorXd h_dot_des, static_com_wrench;
		cws_margin.GravitationalWrench_WF(com_position, cws_data, static_com_wrench);

		for(unsigned int j=0; j<nodes_num; j++){
			full_initial_point(4*j) = final_heuristics_com_pos(0);
			full_initial_point(4*j+1) = 0.0;
			full_initial_point(4*j+2) = final_heuristics_com_pos(1);
			full_initial_point(4*j+3) = 0.0;
		}
		if(optimize_time_duration) full_initial_point(4*nodes_num) = heuristic_move_base_duration;



	}

}

void VertexBasedTO::evaluateBounds(Eigen::Ref<Eigen::VectorXd> full_state_lower_bound,
                                         Eigen::Ref<Eigen::VectorXd> full_state_upper_bound,
                                         Eigen::Ref<Eigen::VectorXd> full_constraint_lower_bound,
                                         Eigen::Ref<Eigen::VectorXd> full_constraint_upper_bound)
{

	if (weights_as_decision_var){
		double margin = 0.01;
		for(unsigned int j=0; j<nodes_num; j++){
			(full_state_lower_bound.segment(j*states_per_knot,com_position.size())).setConstant(-2e19);
			(full_state_upper_bound.segment(j*states_per_knot,com_position.size())).setConstant(2e19);
			(full_state_lower_bound.segment(j*states_per_knot+com_position.size(),nCols)).setConstant(0.0);
			(full_state_upper_bound.segment(j*states_per_knot+com_position.size(),nCols)).setConstant(1.0);

			(full_constraint_upper_bound.segment(j*constraints_per_knot, 1)).setConstant(1.0-margin);
			(full_constraint_lower_bound.segment(j*constraints_per_knot, 1)).setConstant(1.0-margin);

			(full_constraint_upper_bound.segment(j*constraints_per_knot+1, 3)).setZero();
			(full_constraint_lower_bound.segment(j*constraints_per_knot+1, 3)).setZero();
		}
	}else{

		full_state_lower_bound.setConstant(-2e19);
		full_state_upper_bound.setConstant(2e19);
		if(optimize_time_duration){
			full_state_lower_bound(nodes_num*4) = min_move_base_duration;
			full_state_upper_bound(nodes_num*4) = max_move_base_duration;
		}

		full_constraint_upper_bound.setZero();
		full_constraint_lower_bound.setZero();
//
//		full_constraint_upper_bound(1) = 1.0;
//		full_constraint_lower_bound(1) = 0.0;
//
//		full_constraint_upper_bound(2) = 1.0;
//		full_constraint_lower_bound(2) = 0.0;

	}

}

void VertexBasedTO::evaluateConstraints(Eigen::Ref<Eigen::VectorXd> full_constraint,
                                              const Eigen::Ref<const Eigen::VectorXd>& decision_var)
{

	if (weights_as_decision_var){
		// TODO: 1) add speed as decision variable, 2) solve the problem in the trajectory, 3) set the constraint w = V * lambda rather than tau = c x f
		Eigen::VectorXd sum_weights;
		sum_weights.resize(nodes_num);

		for(unsigned int j=0; j<nodes_num; j++){
			sum_weights(j) = (decision_var.segment(j*states_per_knot + com_position.size(),nCols)).sum();
			full_constraint(j*constraints_per_knot) = sum_weights(j);
			Eigen::VectorXd lambdas(nCols);

			if(nRows == 6){
				// wrench is in this order: wrench = (Tau x, Tau y, Tau z, Fx, Fy, Fz)^T
				wrench = A_v_description.block(0,0,nRows,nCols-1)*decision_var.segment(com_position.size(),nCols-1);
				Eigen::Vector3d lin_force = iit::rbd::linearPart(wrench);
				Eigen::Vector3d torque = iit::rbd::angularPart(wrench);
				com_position = decision_var.segment(0,3);
				full_constraint.segment(j*constraints_per_knot+1,3) = com_position.cross(lin_force) - torque;

//				Eigen::Vector3d com_pos, com_vel, com_acc, com_orient, com_orient_d, com_orient_dd;
//				com_orient.setZero();
//				com_orient_d.setZero();
//				com_orient_dd.setZero();
//
//				com_pos(0) = decision_var(constraints_per_knot*j);
//				com_pos(1) = decision_var(constraints_per_knot*j+1);
//				com_pos(2) = 1.0;
//
//				com_vel(0) = decision_var(4*j+1);
//				com_vel(1) = decision_var(4*j+3);
//				com_vel(2) = 0.0;
//
//				if(j == 0){
//					com_acc(0) = (decision_var(1))/sampling_time;
//					com_acc(1) = (decision_var(3))/sampling_time;
//					com_acc(2) = 0.0;
//				}else{
//					com_acc(0) = (decision_var(4*j+1) - decision_var(4*(j-1)+1))/sampling_time;
//					com_acc(1) = (decision_var(4*j+3) - decision_var(4*(j-1)+3))/sampling_time;
//					com_acc(2) = 0.0;
//				}
//
//				cws_margin.IntertialWrench_WF(com_pos,
//												com_vel,
//												com_acc,
//												com_orient,
//												com_orient_d,
//												com_orient_dd,
//												base_inertia_mat,
//												inertial_wrench);

			}else{
				lambdas.setZero();
				lambdas = decision_var.segment(j*states_per_knot + com_position.size(),nCols);
				Eigen::VectorXd static_com_wrench;
				Eigen::Vector3d torque, static_com_torque;
				Eigen::Vector3d com_pos;
				com_pos(0) = decision_var(states_per_knot*j);
				com_pos(1) = decision_var(states_per_knot*j+1);
				com_pos(2) = decision_var(states_per_knot*j+2);

				torque = reduced_A_v_descr*lambdas.segment(1,nCols-1);
				cws_margin.GravitationalWrench_WF(com_pos, cws_data, static_com_wrench);
//				prt(torque.transpose())
//				prt(static_com_wrench.transpose())
				full_constraint.segment(j*constraints_per_knot+1,3) = torque - static_com_wrench.segment(0,3);

			}
		}
		//    prt(full_constraint);

	}else{
		Eigen::VectorXd h_dot_des, static_com_wrench, intertial_wrench;
		com_position(0) = decision_var(0);
		com_position(1) = decision_var(1);
		com_position(2) = 2.0;
		cws_margin.GravitationalWrench_WF(com_position, cws_data, static_com_wrench);

//		// TODO calcolare omega e omega_d (che è diverso orient_d e orient_dd)
//		cws_margin.IntertialWrench_WF(com_pos,
//				com_vel,
//				com_acc,
//				com_orient,
//				com_orient_d,
//				com_orient_dd,
//				base_inertia_mat,
//				inertial_wrench);
//
//		h_dot_des = intertial_wrench - static_com_wrench;
//		wrench = h_dot_des;

//		weights = pseudoInverse * wrench + biasTerm;
		//	prt(static_hdot)
		//	prt(weights)
		//	Eigen::VectorXd check_wrench = reduced_A_v_descr*weights;
		//	prt(check_wrench)
		//    full_constraint.segment(0,nCols) = weights;
		double time_horizon;
		if(optimize_time_duration){
			time_horizon = decision_var(4*nodes_num);
		}else{
			time_horizon = heuristic_move_base_duration;
		}
		sampling_time = time_horizon/(double)nodes_num;
		full_constraint(0) = (decision_var(0)-initial_com_pos(0))/sampling_time - decision_var(1);
		full_constraint(1) = (decision_var(2)-initial_com_pos(1))/sampling_time - decision_var(3);
		for(unsigned int j=1; j<nodes_num; j++){
			full_constraint(2*j) = (decision_var(4*j)-decision_var(4*(j-1)))/sampling_time - decision_var(4*j+1);
			full_constraint(2*j+1) = (decision_var(4*j+2)-decision_var(4*(j-1)+2))/sampling_time - decision_var(4*j+3);
		}
		// Initial and final velocity constraits
		full_constraint(nodes_num*2) = decision_var(1);
		full_constraint(nodes_num*2+1) = decision_var(3);
		full_constraint(nodes_num*2+2) = decision_var((nodes_num-1)*4 + 1);
		full_constraint(nodes_num*2+3) = decision_var((nodes_num-1)*4 + 3);

//		full_constraint(0) = weights.sum();
//		full_constraint(1) = weights.minCoeff();
//		full_constraint(2) = weights.maxCoeff();
		if (print_all) prt(full_constraint);

	}


}


void VertexBasedTO::evaluateCosts(double& cost,
                                        const Eigen::Ref<const Eigen::VectorXd>& decision_var)
{

	if (weights_as_decision_var){
			Eigen::VectorXd lambdas(nCols);
			Eigen::VectorXd tmp_cost1, tmp_cost2;
			tmp_cost1.resize(nodes_num);
			tmp_cost1.setZero();
			tmp_cost2.resize(nodes_num);
			tmp_cost2.setZero();

			for(unsigned int j=0; j<nodes_num; j++){
				lambdas = decision_var.segment(j*states_per_knot + com_position.size(),nCols);
				//		std::cout<<"weights: "<<lambdas<<std::endl;
				for (int i = 0; i <nCols; i++){
					tmp_cost1(j) += (lambdas(i)- mean_coeff)*(lambdas(i)- mean_coeff);
				}

	//			if(nRows==3){
	//				Eigen::VectorXd static_com_wrench;
	//				Eigen::Vector3d torque, static_com_torque;
	//				Eigen::Vector3d com_pos;
	//				com_pos(0) = decision_var(states_per_knot*j);
	//				com_pos(1) = decision_var(states_per_knot*j+1);
	//				com_pos(2) = decision_var(states_per_knot*j+2);
	//
	//				torque = reduced_A_v_descr*lambdas.segment(1,nCols-1);
	//				cws_margin.GravitationalWrench_WF(com_pos, static_com_wrench);
	//				static_com_torque = static_com_wrench.segment(0,3);
	////				std::cout<<"static torque: "<<static_com_torque.transpose()<<std::endl;
	////				std::cout<<"torque: "<<torque.transpose()<<std::endl;
	//				tmp_cost2(j) = (torque-static_com_torque).transpose()*(torque-static_com_torque);
	//			}else{
	//				tmp_cost2(j) = 0.0;
	//			}
			}
	//		cost = -lambdas.minCoeff();// + 0.001*(wrench(5)-9.81)*(wrench(5)-9.81);
	//		cost = tmp_cost + (decision_var(1)-des_height)*(decision_var(1)-des_height);
			cost = tmp_cost1.sum();// + (com_position(2)-des_height)*(com_position(2)-des_height);
	}else{

		Eigen::VectorXd h_dot_des, static_com_wrench;
		iit::rbd::ForceVector inertial_wrench;
		Eigen::VectorXd tmp_cost1, tmp_cost2;
		tmp_cost1.resize(nodes_num);
		tmp_cost1.setZero();
		tmp_cost2.resize(nodes_num);
		tmp_cost2.setZero();
		Eigen::Vector3d com_pos, com_vel, com_acc, com_orient, com_orient_d, com_orient_dd;
		com_orient.setZero();
		com_orient_d.setZero();
		com_orient_dd.setZero();

		double time_horizon;
		if(optimize_time_duration){
			time_horizon = decision_var(4*nodes_num);
		}else{
			time_horizon = heuristic_move_base_duration;
		}

		sampling_time = time_horizon/(double)nodes_num;

		for(unsigned int j=0; j<nodes_num; j++){
			com_pos(0) = decision_var(4*j);
			com_pos(1) = decision_var(4*j+2);
			com_pos(2) = 1.0;

			com_vel(0) = decision_var(4*j+1);
			com_vel(1) = decision_var(4*j+3);
			com_vel(2) = 0.0;

			if(j == 0){
				com_acc(0) = (decision_var(1))/sampling_time;
				com_acc(1) = (decision_var(3))/sampling_time;
				com_acc(2) = 0.0;
			}else{
				com_acc(0) = (decision_var(4*j+1) - decision_var(4*(j-1)+1))/sampling_time;
				com_acc(1) = (decision_var(4*j+3) - decision_var(4*(j-1)+3))/sampling_time;
				com_acc(2) = 0.0;
			}

			cws_margin.GravitationalWrench_WF(com_pos, cws_data, static_com_wrench);

// TODO calcolare omega e omega_d (che è diverso orient_d e orient_dd)
			cws_margin.IntertialWrench_WF(com_pos,
											com_vel,
											com_acc,
											com_orient,
											com_orient_d,
											com_orient_dd,
											base_inertia_mat,
											inertial_wrench);

			h_dot_des = inertial_wrench - static_com_wrench;

			wrench = h_dot_des;
            //shift everything in the origin
            Eigen::VectorXd wrench_gi_avg;
            wrench_gi_avg = wrench - vertex_avg;


			if(use_pseudo_inverse){
				weights = pseudoInverse * wrench + biasTerm;
//				prt(weights.transpose())
				//alternative way removing the avg
				//VectorXd weights_avg;weights_avg.resize(nCols);
				//weights= pseudoInverse_avg * wrench_gi_avg + biasTerm_avg;
//				prt(vertex_avg.transpose())
//				prt(wrench_gi_avg.transpose())
//				prt(biasTerm_avg.transpose())
//				prt(pseudoInverse_avg.transpose())
//				prt(weights_avg.transpose())

			}else{
				std::cout<<"-------------- Starting LP program -----------"<<std::endl;
//				std::cout<<"rows:"<<A_v_description.rows()<<" cols: "<<A_v_description.cols()<<std::endl;
				Eigen::MatrixXd V_test(A_v_description.rows(),A_v_description.cols());
				V_test = A_v_description;
//				V_test << 1.0, 2.0, 2.0, 1.0,
//							1.0, 1.0, 2.0, 2.0;
			    dwl::solver::OptimizationSolver* solver2 = new dwl::solver::IpoptNLP();
			    OptWeights opt_weights;
			    solver2->setOptimizationModel(& opt_weights);
			    Eigen::VectorXd avg, tmp_avg;
			    Eigen::MatrixXd A_avg;
				Eigen::VectorXd wrench_gi;
				wrench_gi.resize(6);
				wrench_gi = wrench;
			    avg.resize(V_test.rows());
			    tmp_avg.resize(V_test.rows());
			    A_avg.resize(V_test.rows(),V_test.cols());
				for(unsigned int j=0; j<V_test.rows(); j++){
					tmp_avg(j) = (V_test.block(j,0,1,V_test.cols())).sum();
					avg(j) = tmp_avg(j)/(double)V_test.cols();
					for(unsigned int i=0; i<V_test.cols(); i++){
						A_avg(j,i) = V_test(j,i)-avg(j);
					}
					wrench_gi(j) = wrench_gi(j)- avg(j);
				}
			    opt_weights.setInitialState(A_avg, wrench_gi);
			    std::cout<<" i am here"<<std::endl;
//			    solver2->setFromConfigFile("config/ipopt_config.yaml");//should be in the folder
                solver2->init();
			    solver2->compute(10.0);
			    Eigen::VectorXd solution = solver2->getSolution();
			    margin = solution(V_test.cols());
			    weights = solution.segment(0,V_test.cols());
//			    weights = pseudoInverse * wrench + biasTerm;
			}
			//		prt(weights)

			for (int i = 0; i <nCols; i++){
				tmp_cost1(j) += (weights(i)- mean_coeff)*(weights(i)- mean_coeff);
			}

			tmp_cost2(j) = margin;
		}
//		wrench = reduced_A_v_descr*weights;
		//if (print_all) prt(decision_var)

//		prt(h_dot)
//		prt(wrench)
//		cost = (decision_var(2)-des_height)*(decision_var(2)-des_height);
//		cost = (wrench(5) - static_com_wrench(5))*(wrench(5) - static_com_wrench(5));
//		Eigen::Vector3d torque_diff = wrench.segment(0,3) - h_dot.segment(0,3);
//		Eigen::Vector3d force_diff = wrench.segment(3,3) - h_dot.segment(3,3);
//		cost = 0.5*(weights.minCoeff() - mean_coeff)*(weights.minCoeff() - mean_coeff) + decision_var(3)*decision_var(3) + decision_var(4)*decision_var(4) + (decision_var(5) - h_dot_des(5))*(decision_var(5) - h_dot_des(5));
//		cost = (weights.maxCoeff() - mean_coeff)*(weights.maxCoeff() - mean_coeff);
		//		cost = tmp_cost;
		if(optimize_time_duration){
			cost = tmp_cost1.sum();// + 0.1*decision_var(4*nodes_num)*decision_var(4*nodes_num);
		}else{
			cost = tmp_cost1.sum();
			//		cost = (weights.minCoeff() - mean_coeff)*(weights.minCoeff() - mean_coeff);
			//		cost = (weights.minCoeff() - mean_coeff)*(weights.minCoeff() - mean_coeff); // works good! (It's convex)
		}


		if(print_all){
			prt(mean_coeff)
				prt(wrench)
				prt(weights.minCoeff())
				prt(weights.maxCoeff())
				prt(weights.sum())
				prt(cost)
		}

	}

}

//void VertexBasedTO::evaluateCostGradient(Eigen::MatrixXd& gradient,
//                                          const Eigen::VectorXd& decision_var){
//    gradient(0,0) = 0.0;
//    gradient(0,1) = 0.0;
//    gradient(0,2) = -1.0;
//
//}

/**
 * @brief Evaluates the jacobian of the constraint function given a current decision
 * state
 * @param Eigen::MatrixXd& Jacobian of the constraint function
 * @param const Eigen::VectorXd& Decision vector
 */

//void VertexBasedTO::evaluateConstraintJacobian(Eigen::MatrixXd& jacobian,
//                                               const Eigen::VectorXd& decision_var){
////    for(int j = 0;j<nRows;j++){
////        jacobian(j,0) = -A_hs(j,1)*gravitoInertialWrench(rbd::LZ)+A_hs(j,2)*gravitoInertialWrench(rbd::LY);
////        jacobian(j,1) =  A_hs(j,0)*gravitoInertialWrench(rbd::LZ)-A_hs(j,2)*gravitoInertialWrench(rbd::LX);
////        jacobian(j,2) = -A_hs(j,0)*gravitoInertialWrench(rbd::LY)+A_hs(j,1)*gravitoInertialWrench(rbd::LX);
////        jacobian(j,3) = 1.0;
////    }
//}

void VertexBasedTO::evaluateSolution(const Eigen::Ref<const Eigen::VectorXd>& solution)
{
    if(print_all) std::cout << "the optimization solution is " << solution.transpose() << std::endl;

}

bool VertexBasedTO::setInitialState(const Eigen::Vector3d base_posW_inital_pos,
		const Eigen::Vector3d base_posW_heuristic_guess,
		const Eigen::Vector3d base_orient,
		const dog::LegDataMap<double> mu,
		const dog::LegDataMap<Vector3d> normals,
		const dog::LegDataMap<double> max_normal_grf,
		const dog::LegDataMap<Vector3d> footPosW,
		const dog::LegID swing_index,
		const iit::rbd::Matrix66d inertia_mat,
		const bool use_weights_as_dv,
		const unsigned int nodes_n,
		const FeasibleWrenchPolytope_API::CWSData cwc,
		const FeasibleWrenchPolytope_API::TLSData tls,
		const FeasibleWrenchPolytope_API::FWSData fwp_option)
{
	std::cout<<"set initial states"<<std::endl;
	optimize_time_duration = fwp_option.optimize_time;
	nodes_num = nodes_n;
	weights_as_decision_var = use_weights_as_dv;
	bool flag = false;
	//this could come from heuristics
	rbd::linearPart(this->basePoseW0) = base_posW_inital_pos;
	rbd::angularPart(this->basePoseW0) = base_orient;

	for (int leg = iit::dog::LF; leg <= iit::dog::RH; leg++){
		stanceCWS[leg] = true;
	}
	stanceCWS[swing_index] = false;

	Matrix3d R = commons::rpyToRot(base_orient);
	dog::LegDataMap<Vector3d> footPosB;
	for (int leg = iit::dog::LF; leg <= iit::dog::RH; leg++){
		footPosB[leg] = R*(footPosW[leg] - tls.comAWPApproxWF);
	}
	if (print_all){
		std::cout<<"feet positions in the WF:"<<std::endl;
		std::cout<<footPosW<<std::endl;
	}

	Eigen::MatrixXd fixed_base_jacobian;
	fixed_base_jacobian.resize(12,12);
	fixed_base_jacobian.setZero();

	cws_margin.get_jacobian(fwp_option.constraint_type, footPosB, fixed_base_jacobian);

//	prt(fixed_base_jacobian)
	//prt(stanceCWS)
	//prt(normals)
	//prt(mu)
	//prt(R)
	//prt(base_posW_inital_pos)

	Ic.setIdentity();

	cws_data.friction = mu;
	cws_data.stance_legs_flag = stanceCWS;
	cws_data.normal = normals;
	cws_data.stance_feet_WF = footPosW;
	cws_data.max_normal_force = max_normal_grf;
	cws_data.robot_mass = cwc.robot_mass;

	tls_data.stance_feet_WF = footPosW;
	tls_data.max_torque = tls.max_torque;
	tls_data.legs_grav_torque = tls.legs_grav_torque;
	tls_data.fixed_base_jac = fixed_base_jacobian;

	//compute the wrench set (either fwp or awp)
	Eigen::MatrixXd feasibility_lims_v_descr;

	flag = cws_margin.feasibility_wrench_set(tls_data, cws_data, fwp_option, feasibility_lims_v_descr);

	A_v_description.resize(feasibility_lims_v_descr.rows(), feasibility_lims_v_descr.cols());
	Eigen::VectorXd tmp_avg;
	tmp_avg.resize(6);
	cws_margin.remove_average_point(fwp_option, feasibility_lims_v_descr, tmp_avg, A_v_description, tmp_avg);
	A_v_description *= 0.9;
    if(print_all) prt(A_v_description); // each column is a vertex of the polytope
	Eigen::MatrixXd tmp, I, I_null;
	nRows = A_v_description.rows();
	nCols = A_v_description.cols();
	if(nCols!=0) mean_coeff = 1.0/(double)nCols; else mean_coeff = 0.0;
	std::cout<<"cols: "<<nCols<<std::endl;
	std::cout<<"rows: "<<nRows<<std::endl;
//	std::cout<<"mean coeff: "<<mean_coeff<<std::endl;
	wrench.resize(nRows);
	weights.resize(nCols);
	weights.setConstant(mean_coeff);
	I.resize(nRows, nRows);
	I.setIdentity();
	I *= 0.001;
	reduced_A_v_descr = A_v_description.block(0,1,nRows,nCols-1);
	if(print_all) prt(reduced_A_v_descr);
	tmp = reduced_A_v_descr*reduced_A_v_descr.transpose() + I;
	tmp = tmp.inverse();
	//prt(tmp)
	I_null.resize(nCols-1,nCols-1);
	I_null.setIdentity();
	Eigen::MatrixXd reducedPseudoInverse = reduced_A_v_descr.transpose()*tmp;
	//pseudoInverse = fullPseudoInverse.block(1,0,nCols-1,nRows);
	Eigen::MatrixXd nullSpaceProjector = (I_null - reducedPseudoInverse*reduced_A_v_descr);
	Eigen::VectorXd mean_coeff_vec;
	mean_coeff_vec.resize(nCols-1);
	mean_coeff_vec.setConstant(mean_coeff);
	//std::cout<<mean_coeff_vec.transpose()<<std::endl;
	biasTerm = nullSpaceProjector * mean_coeff_vec;
	des_height = base_posW_inital_pos(2);
	//std::cout<<pseudoInverse<<std::endl;
	//cws_margin.getHSdescription(A_v_description, A_hs);
	//prt(reducedPseudoInverse)
	pseudoInverse = reducedPseudoInverse;
	initial_com_pos = base_posW_inital_pos;
	final_heuristics_com_pos = base_posW_heuristic_guess;
	if(print_all) prt(pseudoInverse);

	//compute the centroid of the polytope
	vertex_avg.setZero();vertex_avg.resize(nRows);
	for (int i=0; i<nCols;i++)
	{
		vertex_avg+=A_v_description.col(i);
	}
	vertex_avg/=(double)nCols;
	//shift the all vertices in the matrix
	A_v_description_avg = A_v_description.colwise() - vertex_avg;
	pseudoInverse_avg.resize(nCols, nRows);
	//compute the pseudoinverse
	pseudoInverse_avg = commons::psdInv(A_v_description_avg, 1E-06);//A_v_description_avg.transpose()*(A_v_description_avg*A_v_description_avg.transpose() + I).inverse();
	//compute the nullspace projector
	I_null.resize(nCols,nCols);
	I_null.setIdentity();
	Eigen::MatrixXd nullSpaceProjector_avg = (I_null - pseudoInverse_avg*A_v_description_avg);
	//compute the bias term
	Eigen::VectorXd mean_coeff_vec_avg;
	mean_coeff_vec_avg.resize(nCols);
	mean_coeff_vec_avg.setConstant(mean_coeff);
	biasTerm_avg = nullSpaceProjector_avg * mean_coeff_vec_avg;

	return flag;
}

void VertexBasedTO::checkPoint(Eigen::Vector3d & inputPoint)
{
    cws_margin.GravitationalWrench_WF(inputPoint, cws_data);
    //prt(A_hs_slack.block(0,0,nRows,6))
//    prt(cws_margin.getGravitoInertialWrench().transpose())
    VectorXd output(nRows);
    output = A_hs_slack.block(0,0,nRows,6)*cws_margin.getGravitoInertialWrench();
    prt(output.transpose())
}

bool VertexBasedTO::optimize_next_com_pos(const Eigen::Vector3d final_pos_CoM_W0,
										const Eigen::Vector3d actual_pos_CoM_W0,
										const Eigen::Vector3d actual_orient,
										const iit::dog::LegDataMap<double> mu,
										const iit::dog::LegDataMap<Eigen::Vector3d> terrain_normal,
										const iit::dog::LegDataMap<Eigen::Vector3d> footPosWF,
										const iit::dog::LegID swing_leg_index,
										const iit::rbd::Matrix66d inertia_mat,
										const double heuristic_move_base,
										const bool use_weights_as_dv,
										const unsigned int nodes_n,
										const FeasibleWrenchPolytope_API::CWSData cwc,
										const FeasibleWrenchPolytope_API::TLSData tls,
										const FeasibleWrenchPolytope_API::FWSData fwp_option,
										Eigen::VectorXd & sol){

    if(print_all){
    	prt(actual_pos_CoM_W0)
		prt(actual_orient)
		prt(mu)
		prt(terrain_normal)
		prt(footPosWF)
		prt(swing_leg_index)
    }

	iit::dog::LegDataMap<double> max_normal_grf;
    max_normal_grf = cwc.max_normal_force;
    double duty_fac = 0.5;
    double min_swing_time = 0.2;
    heuristic_move_base_duration = heuristic_move_base;
    min_move_base_duration = 0.7;
    max_move_base_duration = heuristic_move_base + 2.0;

	bool flag = false;
//  Vertices based optimization problem
    dwl::solver::OptimizationSolver* solver = new dwl::solver::IpoptNLP();
    solver->setOptimizationModel(this);
    flag = setInitialState(actual_pos_CoM_W0,
    							final_pos_CoM_W0,
    							actual_orient,
								mu,
								terrain_normal,
								max_normal_grf,
								footPosWF,
								swing_leg_index,
								inertia_mat,
								use_weights_as_dv,
								nodes_n,
								cwc,
								tls,
								fwp_option);

    std::string pack_path = ros::package::getPath(std::string("contact_wrench_set")) +
						  "/config/ipopt_configNL.yaml";

	std::cout<<"I am here!"<<std::endl;
    solver->setFromConfigFile(pack_path);//should be in the folder
	std::cout<<"I am here 2!"<<std::endl;
    clock_t start = clock();
    std::cout<<"I am here 2.5!"<<std::endl;
    solver->init();
   	std::cout<<"I am here 3!"<<std::endl;
    solver->compute();
	std::cout<<"I am here 4!"<<std::endl;
	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf ("seconds needed for the trajectory optimization: %f  \n", seconds);

    //sol = solver->getSolution();


    return flag;

}
