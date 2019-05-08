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
/// \file mink_glpk_algo2.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)

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


//  g++ glpk_simple_test.cpp -L /usr/local/lib -l glpk -o glpk_simple_test
//  ./glpk_simple_test


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


using namespace boost::numeric::ublas;

/// \class Dedicated to finding Minkowski vertices with linear programming techniques for V-polytopes.
class MS_LP {

public:
  MS_LP(
      const matrix< double >& m1,
      const matrix< double >& m2,
      unsigned int u,
      unsigned int v,
      double tol=0.000001)
  {
    _tolerance = tol;
    _dimension = m1.size1()-1;
    _numberOfGenerators1 = m1.size2()-1;
    _numberOfGenerators2 = m2.size2()-1;

//	std::cout<<"i am here"<<std::endl;
//	std::cout<<"u "<<u<<std::endl;
//	std::cout<<"v "<<v<<std::endl;
//	std::cout<<"n gen 1 "<<_numberOfGenerators1<<std::endl;
//	std::cout<<"n gen 2 "<<_numberOfGenerators1<<std::endl;
//	std::cout<<"tolerance "<<_tolerance<<std::endl;

    // Declare variables
    //int ia[1+dimension+2], ja[1+numberOfGenerators1+numberOfGenerators2];
    int* ia = new int[(_dimension+3)*(1+_numberOfGenerators1+_numberOfGenerators2)];
    int* ja = new int[(_dimension+3)*(1+_numberOfGenerators1+_numberOfGenerators2)];
    double* ar = new double[(_dimension+3)*(1+_numberOfGenerators1+_numberOfGenerators2)];

    // Create problem
    _linearProblemg = glp_create_prob();
    glp_set_prob_name(_linearProblem, "Minkowski vertices");
    glp_set_obj_dir(_linearProblem, GLP_MAX);
    glp_add_rows(_linearProblem, _dimension+2);
    glp_add_cols(_linearProblem, _numberOfGenerators1+_numberOfGenerators2);
    // Deal with the objective function f* = max(2 - x_u - x_v) and variables x1, ... , xn
    {for (unsigned int col_number=1; col_number<=_numberOfGenerators1+_numberOfGenerators2; ++col_number) {
      if (col_number != u && col_number != v+_numberOfGenerators1)
        glp_set_obj_coef(_linearProblem, col_number,  0.);
      else
        glp_set_obj_coef(_linearProblem, col_number, -1.);
      // Variables are positive
      glp_set_col_bnds(_linearProblem, col_number, GLP_LO, 0., 0.);
    }}
   // The main loop.
    unsigned int counter=0;
    //{for (unsigned int row_number=1; row_number<=_dimension; ++row_number) {
    //  // Build the row name.
    //  std::ostringstream stream_;
    //  stream_ << "r";
    //  stream_ << row_number;
    //  std::string row_number_string = stream_.str();
    //  glp_set_row_name(_linearProblem, row_number, row_number_string.c_str());
    //  double sec_mem =  m1(row_number,u)+m2(row_number,v);
    //  // Insert the second member as an equality constraint. It gives under global notations :
    //  // r_i : a_{i,1} x_1 + ... + a_{i,k} x_k + a_{i,k+1} x_{k+1} + ... + a_{i,k+l} x_{k+l} = a_{i,u}+a_{i,k+v}
    //  glp_set_row_bnds(_linearProblem, row_number, GLP_FX, sec_mem, sec_mem);
//
    //  {for (unsigned int col_number=1; col_number<=_numberOfGenerators1; ++col_number) {
    //    counter++;
    //    ia[counter] = row_number;
    //    ja[counter] = col_number;
    //    ar[counter] = m1(row_number,col_number);
    //  }}
    //  {for (unsigned int col_number=_numberOfGenerators1+1; col_number<=_numberOfGenerators1+_numberOfGenerators2; ++col_number) {
    //    counter++;
    //    ia[counter] = row_number;
    //    ja[counter] = col_number;
    //    ar[counter] = m2(row_number,col_number-_numberOfGenerators1);
    //  }}
    //}}
    //{for (unsigned int row_number=_dimension+1; row_number<=_dimension+2; ++row_number) {
    //  // Build the row name.
    //  std::ostringstream stream_;
    //  stream_ << "r";
    //  stream_ << row_number;
    //  std::string row_number_string = stream_.str();
    //  glp_set_row_name(_linearProblem, row_number, row_number_string.c_str());
    //  // Insert the second member as a constraint equal to 1.
    //  glp_set_row_bnds(_linearProblem, row_number, GLP_FX, 1., 1.);
//
    //  {for (unsigned int col_number=1; col_number<=_numberOfGenerators1; ++col_number) {
    //    counter++;
    //    ia[counter] = row_number;
    //    ja[counter] = col_number;
    //    if (row_number == _dimension+1)
    //      ar[counter] = 1.;
    //    else
    //      ar[counter] = 0.;
    //  }}
    //  {for (unsigned int col_number=_numberOfGenerators1+1; col_number<=_numberOfGenerators1+_numberOfGenerators2; ++col_number) {
    //    counter++;
    //    ia[counter] = row_number;
    //    ja[counter] = col_number;
    //    if (row_number == _dimension+1)
    //      ar[counter] = 0.;
    //    else
    //      ar[counter] = 1.;
    //  }}
    //}}
    // Load the problem.
    glp_load_matrix(_linearProblem, ((_dimension+2)*(_numberOfGenerators1+_numberOfGenerators2)), ia, ja, ar);

    // Solve
    int retval = glp_simplex(_linearProblem, NULL);

    delete[] ia;
    delete[] ja;
    delete[] ar;
    // Print solutions
    //glp_print_sol(P, valString.c_str());

    _objVal = 2 + glp_get_obj_val(_linearProblem);
  }

  ~MS_LP() {
    // Housekeeping
    glp_delete_prob(_linearProblem);
    glp_free_env();
  }

  double getObjectiveValue() const {return _objVal;}

  double getSolutionVectorValue(unsigned int i) const {return glp_get_col_prim(_linearProblem, i);}

  /// Write problem data in CPLEX LP format 
  void dumpCPLEXFile(const std::string& fileOut) const {glp_write_lp(_linearProblem, NULL, fileOut.c_str());}

  /// Write problem data in MPS format
  void dumpMPSFile(const std::string& fileOut) const {
    //glp_write_mps(_linearProblem, GLP_MPS_FILE, NULL, fileOut.c_str());}
    glp_write_mps(_linearProblem, GLP_MPS_DECK, NULL, fileOut.c_str());}

protected:
  unsigned int  _numberOfGenerators1;
  unsigned int  _numberOfGenerators2;
  unsigned int  _dimension;
  double        _objVal;
  double        _tolerance;
  glp_prob*     _linearProblem;
};


int main() {

  // The first argument is the ptop file, the second argument the qhull one.
  std::string qhullFile1;
  std::string qhullFile2;
  std::string fOut;
  std::string version("Version 3.3.0");
  int dimension=6, numberOfGenerators1, numberOfGenerators2, numberOfHalfSpaces1, numberOfHalfSpaces2;
  bool output=false;
  double TOL=0.000001;


  std::string file3 = "my_file.ptop";
  fOut = file3;
  output = true;

  std::string line;

  // To store the results.
  std::vector< unsigned int > allMinkowskiVertices_i;
  std::vector< unsigned int > allMinkowskiVertices_j;
  std::vector< std::vector< double > > allMinkowskiVertices;


  ////////////////////////////////////// My test  /////////////
  dimension=2;
  numberOfGenerators1 = 4;
  numberOfGenerators2 = 4;
  // By convention we do not want not to use the 0 column and line as glpk does not make use of them.
  matrix< double > m1(dimension+1, numberOfGenerators1+1);
  matrix< double > m2(dimension+1, numberOfGenerators2+1);

  m1(1,1) = 0.0;	  m1(2,1) = 0.0;
  m1(1,2) = 1.0;	  m1(2,2) = 0.0;
  m1(1,3) = 0.0;	  m1(2,3) = 1.0;
  m1(1,4) = 1.0;	  m1(2,4) = 1.0;

  m2(1,1) = 0.0;	  m2(2,1) = 0.0;
  m2(1,2) = 1.0;	  m2(2,2) = 0.0;
  m2(1,3) = 0.0;	  m2(2,3) = 1.0;
  m2(1,4) = 1.0;	  m2(2,4) = 1.0;

  boost::timer thisTimer;
  std::cerr.precision(15);
  std::cout.precision(15);
  //std::cerr << "(" << u << "," << v << ") ";
  {for (unsigned int u=1; u<=numberOfGenerators1; ++u) {
    {for (unsigned int v=1; v<=numberOfGenerators2; ++v) {
      // Create a linear problem
      MS_LP ms_lp(m1, m2, u, v);

      //std::cout << "(u,v) = (" << u << "," << v << ") ";
      double z = ms_lp.getObjectiveValue();
      bool minkVertex=false;
      boost::numeric::ublas::vector<double> lp_sol1(dimension+1);
      lp_sol1.clear();
      {for (unsigned int i=1; i<=numberOfGenerators1; ++i) {
         if (ms_lp.getSolutionVectorValue(i) != 0) {
	   boost::numeric::ublas::matrix_column< boost::numeric::ublas::matrix<double> > matCol1(m1, i);
	   lp_sol1 += ms_lp.getSolutionVectorValue(i)*matCol1;
           //std::cout << "x" << i << "=" << ms_lp.getSolutionVectorValue(u) << " ";
         }
      }}
      boost::numeric::ublas::matrix_column< boost::numeric::ublas::matrix<double> > ref_point1(m1, u);
      boost::numeric::ublas::vector<double> diff_point = lp_sol1-ref_point1;
      double norm = norm_2(diff_point);
      //std::cout << "(u,v)=(" << u << "," << v << ")" << std::endl;
      //std::cout << "norm=" << norm << std::endl;
      if (norm < TOL) {
	// We can go on.
	boost::numeric::ublas::vector<double> lp_sol2(dimension+1);
	lp_sol2.clear();
	{for (unsigned int i=numberOfGenerators1+1; i<=numberOfGenerators1+numberOfGenerators2; ++i) {
	  if (ms_lp.getSolutionVectorValue(i) != 0) {
	    boost::numeric::ublas::matrix_column< boost::numeric::ublas::matrix<double> > matCol2(m2, i-numberOfGenerators1);
	    lp_sol2 += ms_lp.getSolutionVectorValue(i)*matCol2;
	    //std::cout << "x" << i << "=" << ms_lp.getSolutionVectorValue(i) << " ";
	  }
	}}
	//std::cout << "OK=" << std::endl;
	boost::numeric::ublas::matrix_column< boost::numeric::ublas::matrix<double> > ref_point2(m2, v);
	boost::numeric::ublas::vector<double> diff_point2 = lp_sol2-ref_point2;
	double norm_b = norm_2(diff_point2);
	//std::cout << "norm_b=" << norm_b << std::endl;
	if (norm_b < TOL)
	  minkVertex = true;
      }
      if (minkVertex == true) {
        std::vector< double > minkVert;
        // We have a Minkowski vertex.
        for (unsigned int coord_count=1; coord_count<=dimension; ++coord_count) {
          minkVert.push_back( m1(coord_count,u) + m2(coord_count,v) );
        }
        allMinkowskiVertices.push_back( minkVert );
	allMinkowskiVertices_i.push_back(u);
	allMinkowskiVertices_j.push_back(v);
 	///
	//std::cout << "Minkowski vertex: ";
	//std::copy(minkVert.begin(), minkVert.end(), std::ostream_iterator<double>(std::cout, " ") );
 	//std::cout << std::endl << "***";
	//std::cout << std::endl << std::endl;
     }
    }} // {for (unsigned int v=1; v<=numberOfGenerators2; ++v) {
  }} // {for (unsigned int u=1; u<=numberOfGenerators1; ++u) {
  std::cout << "TIME=" << thisTimer.elapsed() << std::endl;

  std::ofstream fileOUT(fOut.c_str(), std::ifstream::out);
  if (!fileOUT) {
    std::string s("Unable to open ");
    s += fOut;
    s += "\n";
    throw std::ios_base::failure(s);
  }
  fileOUT.precision(15);

  // Write the results into the file.
  //fileOUT << "# GENERATORS : V = (v1, ..., vn)" << std::endl;
  //fileOUT << allMinkowskiVertices.size() << std::endl;
  std::vector< unsigned int >::const_iterator iteMV_i = allMinkowskiVertices_i.begin();
  std::vector< unsigned int >::const_iterator iteMV_j = allMinkowskiVertices_j.begin();
  std::vector< std::vector< double > >::const_iterator iteMV = allMinkowskiVertices.begin();
  {for ( ; iteMV != allMinkowskiVertices.end(); ++iteMV,++iteMV_i,++iteMV_j ) {
      //fileOUT << *iteMV_i << "+" << *iteMV_j << " ";
    std::vector< double >::const_iterator iteCoord = iteMV->begin();
    {for ( ; iteCoord != iteMV->end(); ++iteCoord ) {
      //fileOUT << *iteCoord << " ";
      std::cout << *iteCoord << " ";
    }}
    //fileOUT << "(";
    //{for (unsigned int coord_count=1; coord_count<=dimension; ++coord_count) {
	//fileOUT << m1(coord_count,*iteMV_i) << " ";
    //}}
    //fileOUT << ")+(";
    //{for (unsigned int coord_count=1; coord_count<=dimension; ++coord_count) {
	//fileOUT << m2(coord_count,*iteMV_j) << " ";
    //}}
    //fileOUT << ")";
    fileOUT << std::endl;
    std::cout << std::endl;
  }}

  fileOUT.close();

  return 0;
}
