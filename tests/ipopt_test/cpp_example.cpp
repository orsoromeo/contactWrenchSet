// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: cpp_example.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

//#include "IpIpoptApplication.hpp"
//#include "IpSolveStatistics.hpp"
//#include "MyNLP.hpp"
#include <mink_glpk.h>
#include <iostream>

//using namespace Ipopt;


int main() {

  int call_ipopt();
//  void distancePointHyperPlane(const Eigen::VectorXd);

  clock_t start = clock();
  int i;
  for (i=0; i<1;i++){
	  call_ipopt();
  }
  clock_t end = clock();
  float seconds = (float)(end - start) / CLOCKS_PER_SEC;
  printf ("seconds: %f  \n", seconds);
  return 0;
}

int call_ipopt()
{
////////////////////////////////////
//Computing the vertices of the CWS
////////////////////////////////////
//	const matrix< double > m1(3,3);
//    const matrix< double > m2(3,3);
//    unsigned int u = 0;
//    unsigned int v = 1;
//    double tol=0.000001;

	MS_LP ms_lp0;
	ms_lp0.compute_politopix();

///////////////////////////////////////////
//	solving the linear program with ipopt
///////////////////////////////////////////

//  // Create an instance of your nlp...
//  SmartPtr<TNLP> mynlp = new MyNLP();
//
//  // Create an instance of the IpoptApplication
//  //
//  // We are using the factory, since this allows us to compile this
//  // example with an Ipopt Windows DLL
//  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
//
//  // Initialize the IpoptApplication and process the options
//  ApplicationReturnStatus status;
//  status = app->Initialize();
//  if (status != Solve_Succeeded) {
//    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
//    return (int) status;
//  }
//
//  status = app->OptimizeTNLP(mynlp);
//
//  if (status == Solve_Succeeded) {
//    // Retrieve some statistics about the solve
//    Index iter_count = app->Statistics()->IterationCount();
//    std::cout << std::endl << std::endl << "*** The problem solved in " << iter_count << " iterations!" << std::endl;
//
//    Number final_obj = app->Statistics()->FinalObjective();
//    std::cout << std::endl << std::endl << "*** The final value of the objective function is " << final_obj << '.' << std::endl;
//
//  }

  return (int) 0;
}
