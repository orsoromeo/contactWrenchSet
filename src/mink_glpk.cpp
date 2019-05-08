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


//  g++ mink_glpk_algo2.cpp -L /usr/local/lib -l glpk -o mink_glpk_algo2
//  ./mink_glpk_algo2 -p1 cube1.ptop -p2 cube2.ptop -d 6 -o cube3.ptop


#include "contact_wrench_set/mink_glpk.h"

void glp_init_smcp(glp_smcp *parm)
{     parm->msg_lev = GLP_MSG_ALL;
      parm->meth = GLP_PRIMAL;
      parm->pricing = GLP_PT_PSE;
      parm->r_test = GLP_RT_HAR;
      parm->tol_bnd = 1e-7;
      parm->tol_dj = 1e-7;
      parm->tol_piv = 1e-10;
//      parm->obj_ll = -DBL_MAX;
//      parm->obj_ul = +DBL_MAX;
//      parm->it_lim = INT_MAX;
//      parm->tm_lim = INT_MAX;
      parm->it_lim = 1000;
      parm->tm_lim = 2000;
//      std::cout<<"max iter num GLPK: "<<parm->it_lim<<std::endl;
//      std::cout<<"max time limit GLPK: "<<parm->tm_lim<<" [ms]"<<std::endl;
      parm->out_frq = 200;
      parm->out_dly = 0;
      parm->presolve = GLP_OFF;
      return;
}

void MS_LP::setup_lp(const boost::numeric::ublas::matrix< double >& m1,
                     const boost::numeric::ublas::matrix< double >& m2,
                     unsigned int u,
                     unsigned int v,
                     double tol)
{
	_tolerance = tol;
	_dimension = m1.size1()-1;
	_numberOfGenerators1 = m1.size2()-1;
	_numberOfGenerators2 = m2.size2()-1;

	glp_term_out(terminal_output_flag);
	// Declare variables
	int* ia = new int[(_dimension+3)*(1+_numberOfGenerators1+_numberOfGenerators2)];
	int* ja = new int[(_dimension+3)*(1+_numberOfGenerators1+_numberOfGenerators2)];
	double* ar = new double[(_dimension+3)*(1+_numberOfGenerators1+_numberOfGenerators2)];

	// Create problem
	_linearProblem = glp_create_prob();
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
	{for (unsigned int row_number=1; row_number<=_dimension; ++row_number) {
		// Build the row name.
		std::ostringstream stream_;
		stream_ << "r";
		stream_ << row_number;
		std::string row_number_string = stream_.str();
		glp_set_row_name(_linearProblem, row_number, row_number_string.c_str());
		double sec_mem =  m1(row_number,u)+m2(row_number,v);
		// Insert the second member as an equality constraint. It gives under global notations :
		// r_i : a_{i,1} x_1 + ... + a_{i,k} x_k + a_{i,k+1} x_{k+1} + ... + a_{i,k+l} x_{k+l} = a_{i,u}+a_{i,k+v}
		glp_set_row_bnds(_linearProblem, row_number, GLP_FX, sec_mem, sec_mem);

		{for (unsigned int col_number=1; col_number<=_numberOfGenerators1; ++col_number) {
			counter++;
			ia[counter] = row_number;
			ja[counter] = col_number;
			ar[counter] = m1(row_number,col_number);
		}}
		{for (unsigned int col_number=_numberOfGenerators1+1; col_number<=_numberOfGenerators1+_numberOfGenerators2; ++col_number) {
			counter++;
			ia[counter] = row_number;
			ja[counter] = col_number;
			ar[counter] = m2(row_number,col_number-_numberOfGenerators1);
		}}
	}}
	{for (unsigned int row_number=_dimension+1; row_number<=_dimension+2; ++row_number) {
		// Build the row name.
		std::ostringstream stream_;
		stream_ << "r";
		stream_ << row_number;
		std::string row_number_string = stream_.str();
		glp_set_row_name(_linearProblem, row_number, row_number_string.c_str());
		// Insert the second member as a constraint equal to 1.
		glp_set_row_bnds(_linearProblem, row_number, GLP_FX, 1., 1.);

		{for (unsigned int col_number=1; col_number<=_numberOfGenerators1; ++col_number) {
			counter++;
			ia[counter] = row_number;
			ja[counter] = col_number;
			if (row_number == _dimension+1)
				ar[counter] = 1.;
			else
				ar[counter] = 0.;
		}}
		{for (unsigned int col_number=_numberOfGenerators1+1; col_number<=_numberOfGenerators1+_numberOfGenerators2; ++col_number) {
			counter++;
			ia[counter] = row_number;
			ja[counter] = col_number;
			if (row_number == _dimension+1)
				ar[counter] = 0.;
			else
				ar[counter] = 1.;
		}}
	}}
	// Load the problem.
	glp_load_matrix(_linearProblem, ((_dimension+2)*(_numberOfGenerators1+_numberOfGenerators2)), ia, ja, ar);

	// Set control parameters of the simplex method
//	glp_init_smcp(_smcp);
	glp_smcp _smcp;
	smcp = & _smcp;
	glp_init_smcp((glp_smcp *)smcp);
	// Solve
    glp_simplex(_linearProblem, smcp);

	delete[] ia;
	delete[] ja;
	delete[] ar;
	// Print solutions
	//glp_print_sol(P, valString.c_str());

	_objVal = 2 + glp_get_obj_val(_linearProblem);
}

void MS_LP::kill_lp(){
    // Housekeeping
    glp_delete_prob(_linearProblem);
    glp_free_env();
  }

int MS_LP::mink_sum_with_vertices(const Eigen::MatrixXd & first_summand,
		const Eigen::MatrixXd & second_summand,
		Eigen::MatrixXd & A_v_description) {

      if (print_all){
		std::cout<<"first: "<<std::endl;
		std::cout<<first_summand.transpose()<<std::endl;
		std::cout<<"second: "<<std::endl;
		std::cout<<second_summand.transpose()<<std::endl;
      }
	  // The first argument is the ptop file, the second argument the qhull one.
	  std::string qhullFile1;
	  std::string qhullFile2;
//	  std::string fOut;
//	  std::string version("Version 3.3.0");
      unsigned int dimension=6, numberOfGenerators1, numberOfGenerators2;
	  bool output=false;
      //bool ptopformat1=true,ptopformat2=true;
	  double TOL=0.000001;


	  std::string file3 = "my_file.ptop";
//	  fOut = file3;
	  output = false;
	  std::string t1 = "cube4.ptop";
	  qhullFile1 = t1;
      //ptopformat1 = false;

	  std::string t2 = "cube5.ptop";
	  qhullFile2 = t2;
      //ptopformat2 = false;

	  std::ifstream fileIN1(qhullFile1.c_str(), std::ifstream::in);
	  if (!fileIN1) {
		//std::cout<<"Unable to open the first file"<<std::endl;
//	    s += qhullFile1;
//	    s += "\n";
//	    throw std::ios_base::failure(s);
	  }
	  std::ifstream fileIN2(qhullFile2.c_str(), std::ifstream::in);
	  if (!fileIN2) {
		//	std::cout<<"Unable to open the second file"<<std::endl;
//	    s += qhullFile2;
//	    s += "\n";
//	    throw std::ios_base::failure(s);
	  }
	  fileIN1.precision(15);
	  fileIN2.precision(15);


	  // To store the results.
	  std::vector< unsigned int > allMinkowskiVertices_i;
	  std::vector< unsigned int > allMinkowskiVertices_j;
	  std::vector< std::vector< double > > allMinkowskiVertices;
	//
	//  ///////////////////////
	//  // Generators block //
	//  /////////////////////
	//  double val;
	//  if (numberOfGenerators1 != 0) {
	//    for (unsigned int vtx_count=1; vtx_count<=numberOfGenerators1; vtx_count++) {
	//      for (unsigned int coord_count=1; coord_count<=dimension; coord_count++) {
	//        fileIN1 >> val;
	//        m1(coord_count, vtx_count) = val;
	//      }
	//    }
	//  }
	//  if (numberOfGenerators2 != 0) {
	//    for (unsigned int vtx_count=1; vtx_count<=numberOfGenerators2; vtx_count++) {
	//      for (unsigned int coord_count=1; coord_count<=dimension; coord_count++) {
	//        fileIN2 >> val;
	//        m2(coord_count, vtx_count) = val;
	//      }
	//    }
	//  }
	//  fileIN1.close();
	//  fileIN2.close();


	  ////////////////////////////////////// My test  /////////////
	  dimension= first_summand.rows();
	  numberOfGenerators1 = first_summand.cols();
	  numberOfGenerators2 = second_summand.cols();

	  if (print_all)  std::cout<<"dim: "<<dimension<<" numberOfGenerators1: "<<numberOfGenerators1<<" numberOfGenerators2: "<<numberOfGenerators2<<std::endl;
	  // By convention we do not want not to use the 0 column and line as glpk does not make use of them.
      boost::numeric::ublas::matrix< double > m1(dimension+1, numberOfGenerators1+1);
      boost::numeric::ublas::matrix< double > m2(dimension+1, numberOfGenerators2+1);

      for (unsigned int iter1=1; iter1<=numberOfGenerators1; iter1++){
          for(unsigned int j = 1; j<= dimension; j++){
			  m1(j,iter1) = first_summand(j-1,iter1-1);
		  }
	  }
//	  m1(1,1) = 0.0;	  m1(2,1) = 0.0;      m1(3,1) = 0.0;	  m1(4,1) = 0.0;	  m1(5,1) = 0.0;      m1(6,1) = 0.0;
//	  m1(1,2) = 2.0;	  m1(2,2) = 0.0;      m1(3,2) = 0.0;	  m1(4,2) = 0.0;	  m1(5,2) = 0.0;      m1(6,2) = 0.0;
//	  m1(1,3) = 0.0;	  m1(2,3) = 2.0;      m1(3,3) = 0.0;	  m1(4,3) = 0.0;	  m1(5,3) = 0.0;      m1(6,3) = 0.0;
//	  m1(1,4) = 2.0;	  m1(2,4) = 2.0;      m1(3,4) = 0.0;	  m1(4,4) = 0.0;	  m1(5,4) = 0.0;      m1(6,4) = 0.0;
//	  m1(1,5) = 0.0;	  m1(2,5) = 0.0;      m1(3,5) = 1.0;	  m1(4,5) = 0.0;	  m1(5,5) = 0.0;      m1(6,5) = 0.0;
//	  m1(1,6) = 2.0;	  m1(2,6) = 0.0;      m1(3,6) = 1.0;	  m1(4,6) = 0.0;	  m1(5,6) = 0.0;      m1(6,6) = 0.0;
//	  m1(1,7) = 0.0;	  m1(2,7) = 2.0;      m1(3,7) = 1.0;	  m1(4,7) = 0.0;	  m1(5,7) = 0.0;      m1(6,7) = 0.0;
//	  m1(1,8) = 2.0;	  m1(2,8) = 2.0;      m1(3,8) = 1.0;	  m1(4,8) = 0.0;	  m1(5,8) = 0.0;      m1(6,8) = 0.0;

      for (unsigned int iter2=1; iter2<=numberOfGenerators2; iter2++){
          for(unsigned int j = 1; j<= dimension; j++){
			  m2(j,iter2) = second_summand(j-1,iter2-1);
		  }
	  }

//	  m2(1,1) = 0.0;	  m2(2,1) = 0.0;      m2(3,1) = 0.0;	  m2(4,1) = 0.0;	  m2(5,1) = 0.0;      m2(6,1) = 0.0;
//	  m2(1,2) = 1.0;	  m2(2,2) = 2.0;      m2(3,2) = 0.0;	  m2(4,2) = 0.0;	  m2(5,2) = 0.0;      m2(6,2) = 0.0;
//	  m2(1,3) = 2.0;	  m2(2,3) = 0.0;      m2(3,3) = 0.0;	  m2(4,3) = 0.0;	  m2(5,3) = 0.0;      m2(6,3) = 0.0;
//	  m2(1,4) = 0.0;	  m2(2,4) = 0.0;      m2(3,4) = 1.0;	  m2(4,4) = 0.0;	  m2(5,4) = 0.0;      m2(6,4) = 0.0;
//	  m2(1,5) = 1.0;	  m2(2,5) = 2.0;      m2(3,5) = 1.0;	  m2(4,5) = 0.0;	  m2(5,5) = 0.0;      m2(6,5) = 0.0;
//	  m2(1,6) = 2.0;	  m2(2,6) = 0.0;      m2(3,6) = 1.0;	  m2(4,6) = 0.0;	  m2(5,6) = 0.0;      m2(6,6) = 0.0;


	  boost::timer thisTimer;
	  std::cerr.precision(15);
	  std::cout.precision(15);
	  //std::cerr << "(" << u << "," << v << ") ";
	  {for (unsigned int u=1; u<=numberOfGenerators1; ++u) {
	    {for (unsigned int v=1; v<=numberOfGenerators2; ++v) {
	      // Create a linear problem
	      setup_lp(m1, m2, u, v);
	      //std::cout << "(u,v) = (" << u << "," << v << ") ";
          //double z = getObjectiveValue();
	      bool minkVertex=false;
	      boost::numeric::ublas::vector<double> lp_sol1(dimension+1);
	      lp_sol1.clear();
	      {for (unsigned int i=1; i<=numberOfGenerators1; ++i) {
	         if (getSolutionVectorValue(i) != 0) {
		   boost::numeric::ublas::matrix_column< boost::numeric::ublas::matrix<double> > matCol1(m1, i);
		   lp_sol1 += getSolutionVectorValue(i)*matCol1;
	           //std::cout << "x" << i << "=" << getSolutionVectorValue(u) << " ";
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
		  if (getSolutionVectorValue(i) != 0) {
		    boost::numeric::ublas::matrix_column< boost::numeric::ublas::matrix<double> > matCol2(m2, i-numberOfGenerators1);
		    lp_sol2 += getSolutionVectorValue(i)*matCol2;
		    //std::cout << "x" << i << "=" << getSolutionVectorValue(i) << " ";
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

		//std::cout << "Minkowski vertex: ";
		//std::copy(minkVert.begin(), minkVert.end(), std::ostream_iterator<double>(std::cout, " ") );
	 	//std::cout << std::endl << "***";
		//std::cout << std::endl << std::endl;
	     }
	      //else if (z > TOL) {
	        // Interior point or just frontier point not being a vertex
	      //}
	      //else {
	        //std::cerr << "Inconsistent solution for u=" << u << " and v=" << v;
	        //std::cerr << " z=" << z << " ; (x_u, x_v) = (" << getSolutionVectorValue(u);
	        //std::cerr << ", " << getSolutionVectorValue(numberOfGenerators1+v) << ")" << std::endl;
	        // Inconsistent value so print solutions
	        //std::ostringstream stream_;
	        //stream_ << "_lp_";
	        //stream_ << u;
	        //stream_ << "_";
	        //stream_ << v;
	        //stream_ << ".txt";
	        //std::string valString = stream_.str();
	        //dumpFile(valString);
	      //}
	    }} // {for (unsigned int v=1; v<=numberOfGenerators2; ++v) {
	  }} // {for (unsigned int u=1; u<=numberOfGenerators1; ++u) {

	  //std::cout << "TIME=" << thisTimer.elapsed() << std::endl;


	//   The result file containing the Minkowski vertices.
	  if (output == false) {
	    // First remove the .ptop extension file.
	    qhullFile1.resize(qhullFile1.size()-5);

	  }
	  else {
	    // Nothing to do as the option -o provided the file name.
	  }
//	  std::ofstream fileOUT(fOut.c_str(), std::ifstream::out);
//	  if (!fileOUT) {
//	    std::string s("Unable to open ");
//	    s += fOut;
//	    s += "\n";
//	    throw std::ios_base::failure(s);
//	  }
//	  fileOUT.precision(15);

	  // Write the results into the file.
	  A_v_description.resize(dimension,allMinkowskiVertices.size());
	  if (print_all)  std::cout<<"size of the result: "<<dimension<<" x "<<allMinkowskiVertices.size()<<std::endl;
//	  fileOUT << "# GENERATORS : V = (v1, ..., vn)" << std::endl;
//	  fileOUT << allMinkowskiVertices.size() << std::endl;
	  std::vector< unsigned int >::const_iterator iteMV_i = allMinkowskiVertices_i.begin();
	  std::vector< unsigned int >::const_iterator iteMV_j = allMinkowskiVertices_j.begin();
	  std::vector< std::vector< double > >::const_iterator iteMV = allMinkowskiVertices.begin();
	  unsigned int row_iter = 0;
	  {for ( ; iteMV != allMinkowskiVertices.end(); ++iteMV,++iteMV_i,++iteMV_j ) {
	      //fileOUT << *iteMV_i << "+" << *iteMV_j << " ";
	    std::vector< double >::const_iterator iteCoord = iteMV->begin();
	    unsigned int col_iter = 0;
	    {for ( ; iteCoord != iteMV->end(); ++iteCoord ) {
//	      fileOUT << *iteCoord << " ";
	      if (print_all) std::cout << *iteCoord << " ";
	      A_v_description(col_iter, row_iter) = *iteCoord;
	      A_v_description(col_iter, row_iter) = round(A_v_description(col_iter, row_iter));
	      col_iter++;
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
//	    fileOUT << std::endl;
	    if (print_all) std::cout << std::endl;
	    row_iter++;
	  }}

//	  fileOUT.close();
	  return 1;
}
