// mink_glpk_algo2 creates files to be processed by glpk
//     Copyright (C) 2014 : Delos Vincent
//
//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//       \file mink_glpk_algo2.cpp
//       \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)

// g++ mink_glpk_algo2.cpp -L /usr/local/lib -l glpk -I /cygdrive/C/Users/vdelos/Documents/CPP/glpk-4.54/src -o mink_glpk_algo2
// mink_glpk_algo2 -f1 v_desc1.ptop -f2 v_desc2.ptop -d 6 -o v_desc3.ptop
// ./mink_glpk_algo2.exe -f1 o2.ptop -f2 o1.ptop -d 6 -o res2.txt > /dev/null 2>&1
//
// File format for v_desc1.ptop
//L0: Comment
//L1: NumberOfVertices
//L2: V00 V01 V02 ...
//L3: V10 V11 V12 ...
//...


//  g++ mink_glpk_algo2.cpp -L /usr/local/lib -l glpk -o mink_glpk_algo2
//  ./mink_glpk_algo2 -p1 cube1.ptop -p2 cube2.ptop -d 6 -o cube3.ptop

#include <glpk.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <boost/timer.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <Eigen/Dense>


//using namespace boost::numeric::ublas;

/// \class Dedicated to finding Minkowski vertices with linear programming techniques for V-polytopes.
class MS_LP {

public:
  MS_LP(){}

  ~MS_LP(){}

  void setup_lp(const boost::numeric::ublas::matrix< double >& m1,
      const boost::numeric::ublas::matrix< double >& m2,
      unsigned int u,
      unsigned int v,
      double tol=0.000001);

  void kill_lp();

  int mink_sum_with_vertices(const Eigen::MatrixXd & first_summand,
			const Eigen::MatrixXd & second_summand,
			Eigen::MatrixXd & A_v_description);

  double getObjectiveValue() const {return _objVal;}

  double getSolutionVectorValue(unsigned int i) const {return glp_get_col_prim(_linearProblem, i);}

  /// Write problem data in CPLEX LP format 
  void dumpCPLEXFile(const std::string& fileOut) const {glp_write_lp(_linearProblem, NULL, fileOut.c_str());}

  /// Write problem data in MPS format
  void dumpMPSFile(const std::string& fileOut) const {
	  //glp_write_mps(_linearProblem, GLP_MPS_FILE, NULL, fileOut.c_str());}
	  glp_write_mps(_linearProblem, GLP_MPS_DECK, NULL, fileOut.c_str());
  }
  void setPrintAll(bool print_all_flag){
      this->print_all = print_all_flag;
  }

  int terminal_output_flag = 0;

protected:
  unsigned int  _numberOfGenerators1;
  unsigned int  _numberOfGenerators2;
  unsigned int  _dimension;
  double        _objVal;
  double        _tolerance;
  glp_prob*     _linearProblem;
  const glp_smcp *smcp;
  bool print_all = false;

};

