// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: MyNLP.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "MyNLP.hpp"

#include <cassert>

using namespace Ipopt;

/* Constructor. */
MyNLP::MyNLP()
{}

MyNLP::~MyNLP()
{}

bool MyNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // The problem described in MyNLP.hpp has 2 variables, x1, & x2,
  n = 5;

  // one equality constraint,
  m = 1;

  // 2 nonzeros in the jacobian (one for x1, and one for x2),
  nnz_jac_g = 0;

  // and 2 nonzeros in the hessian of the lagrangian
  // (one in the hessian of the objective for x2,
  //  and one in the hessian of the constraints for x1)
  nnz_h_lag = 0;

  // We use the standard fortran index style for row/col entries
  index_style = FORTRAN_STYLE;

  return true;
}

bool MyNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == 5);
  assert(m == 1);

  // x2 has no upper or lower bound, so we set them to
  // a large negative and a large positive number.
  // The value that is interpretted as -/+infinity can be
  // set in the options, but it defaults to -/+1e19
  x_l[0] = -10000.0;
  x_l[1] = -10000.0;
  x_l[2] = -10000.0;
  x_l[3] = -10000.0;
  x_l[4] = -10000.0;

  x_u[0] = +10000.0;
  x_u[1] = +10000.0;
  x_u[2] = +10000.0;
  x_u[3] = +10000.0;
  x_u[4] = +10000.0;

  // we have one equality constraint, so we set the bounds on this constraint
  // to be equal (and zero).
  g_l[0] = g_u[0] = 1.0;
//  g_l[1] = g_u[1] = 0.0;
//  g_l[2] = g_u[2] = 0.0;

  return true;
}

//bool MyNLP::get_scaling_parameters(Number& obj_scaling,
//									bool& use_x_scaling, Index n,
//									Number* x_scaling,
//									bool& use_g_scaling, Index m,
//									Number* g_scaling)
//{
//	obj_scaling = 1.0;
//	}

bool MyNLP::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // we initialize x in bounds, in the upper right quadrant
  x[0] = 0.2;
  x[1] = 0.2;
  x[2] = 0.2;
  x[3] = 0.2;
  x[4] = 0.2;

  return true;
}

bool MyNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  // return the value of the objective function

  obj_value = (-1.0 + x[4]);
  return true;
}

bool MyNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // return the gradient of the objective function grad_{x} f(x)
	assert(n == 5);
  // grad_{x1} f(x): x1 is not in the objective
  grad_f[0] = 0.0;
  grad_f[1] = 0.0;
  grad_f[2] = 0.0;
  grad_f[3] = 0.0;
  grad_f[4] = 1.0;

  // grad_{x2} f(x):
//  Number x2 = x[1];
//  grad_f[1] = -2.0*(x2 - 2.0);

  return true;
}

bool MyNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  // return the value of the constraints: g(x)
  double P0[2] = {1.0,1.0};
  double P1[2] = {1.0,-1.0};
  double P2[2] = {-1.0,1.0};
  double P3[2] = {-1.0,-1.0};
  double P4[2] = {0.0,0.0};
  g[0] = 1.0*x[0] + 1.0*x[1] + 1.0*x[2] + 1.0*x[3] + 1.0*x[4];
//  g[1] = P0[1]*x[0] + P1[1]*x[1] + P2[1]*x[2] + P3[1]*x[3] + P4[1]*x[4];
//  g[2] = x[0] + x[1] + x[2] + x[3] + x[4] - 1.0;

  return true;
}

bool MyNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{

	  double P0[2] = {1.0,1.0};
	  double P1[2] = {1.0,-1.0};
	  double P2[2] = {-1.0,1.0};
	  double P3[2] = {-1.0,-1.0};
	  double P4[2] = {0.0,0.0};


	  if (values == NULL) {
	    // return the structure of the Jacobian

	    // this particular Jacobian is dense
//	    iRow[0] = 0; jCol[0] = 0;
//	    iRow[1] = 0; jCol[1] = 1;
//	    iRow[2] = 0; jCol[2] = 2;
//	    iRow[3] = 0; jCol[3] = 3;
//	    iRow[4] = 0; jCol[4] = 4;

//	    iRow[5] = 1; jCol[5] = 0;
//	    iRow[6] = 1; jCol[6] = 1;
//	    iRow[7] = 1; jCol[7] = 2;
//	    iRow[8] = 1; jCol[8] = 3;
//	    iRow[9] = 1; jCol[9] = 4;

//	    iRow[10] = 2; jCol[10] = 0;
//	    iRow[11] = 2; jCol[11] = 1;
//	    iRow[12] = 2; jCol[12] = 2;
//	    iRow[13] = 2; jCol[13] = 3;
//	    iRow[14] = 2; jCol[13] = 4;
	  }
	  else {
	    // return the values of the Jacobian of the constraints

//	    values[0] = 1.0; // 0,0
//	    values[1] = 1.0; // 0,1
//	    values[2] = 1.0; // 0,2
//	    values[3] = 1.0; // 0,3
//	    values[4] = 1.0; // 0,4

//	    values[5] = P0[1]; // 1,0
//	    values[6] = P1[1]; // 1,1
//	    values[7] = P2[1]; // 1,2
//	    values[8] = P3[1]; // 1,3
//	    values[9] = P4[1]; // 1,4

//	    values[10] = 1.0; // 2,0
//	    values[11] = 1.0; // 2,1
//	    values[12] = 1.0; // 2,2
//	    values[13] = 1.0; // 2,3
//	    values[14] = 1.0; // 2,4
	  }

	  return true;
}

bool MyNLP::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
//  if (values == NULL) {
//    // return the structure. This is a symmetric matrix, fill the lower left
//    // triangle only.
//
//    // element at 1,1: grad^2_{x1,x1} L(x,lambda)
//    iRow[0] = 1;
//    jCol[0] = 1;
//
//    // element at 2,2: grad^2_{x2,x2} L(x,lambda)
//    iRow[1] = 2;
//    jCol[1] = 2;
//
//    // Note: off-diagonal elements are zero for this problem
//  }
//  else {
//    // return the values
//
//    // element at 1,1: grad^2_{x1,x1} L(x,lambda)
//    values[0] = -2.0 * lambda[0];
//
//    // element at 2,2: grad^2_{x2,x2} L(x,lambda)
//    values[1] = -2.0 * obj_factor;
//
//    // Note: off-diagonal elements are zero for this problem
//  }

  return true;
}

void MyNLP::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
			      const IpoptData* ip_data,
			      IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // For this example, we write the solution to the console
  printf("\n\nSolution of the primal variables, x\n");
  for (Index i=0; i<n; i++) {
    printf("x[%d] = %e\n", i, x[i]);
  }

  printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
  for (Index i=0; i<n; i++) {
    printf("z_L[%d] = %e\n", i, z_L[i]);
  }
  for (Index i=0; i<n; i++) {
    printf("z_U[%d] = %e\n", i, z_U[i]);
  }

  printf("\n\nObjective value\n");
  printf("f(x*) = %e\n", obj_value);
}
