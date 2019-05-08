#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <Eigen/Dense>
#include <glpk.h>

class GLPK_wrapper{

public:

GLPK_wrapper(){};

~GLPK_wrapper(){};

void call_glpk();

void solve_lp(const Eigen::VectorXd & obj_fun,
				const Eigen::VectorXd & b, 
				const Eigen::MatrixXd & A, 
				Eigen::VectorXd & sol);

private:

void new_lp();

void delete_lp();

void set_constraints(const Eigen::MatrixXd & A);

void set_coefficients(const Eigen::VectorXd & b);

void set_cost_function(const Eigen::VectorXd & obj_fun);

void solve();

void get_solution(Eigen::VectorXd & sol);

glp_prob * lin_prog;

unsigned int states_number = 0, constraints_number = 0;

};
