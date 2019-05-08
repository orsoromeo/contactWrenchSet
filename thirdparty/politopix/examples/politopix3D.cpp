// politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2017 : Delos Vincent
//
//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU Lesser General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU Lesser General Public License for more details.
//
//     You should have received a copy of the GNU Lesser General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <boost/shared_ptr.hpp>
#include <boost/timer.hpp>
#include <iostream>
#include <string.h>
#include <vector>
#include <cmath>
#include <time.h>
#include "/home/user/Downloads/politopix/trunk/politopixAPI.h"
#include "/home/user/Downloads/politopix/trunk/Generator_Rn.h"
#include "/home/user/Downloads/politopix/trunk/Rn.h"
// export PATH=$PATH:$pol/bin
// linux:  g++ -W  -Wall  -Wextra  -Wswitch  -Wformat  -Wchar-subscripts  -Wparentheses  -Wmultichar  -Wtrigraphs  -Wpointer-arith  -Wcast-align  -Wreturn-type  -Wshadow  -Wundef  -Woverloaded-virtual  -Wno-unused-function -O3 -fno-strict-aliasing test_politopix.cpp  -I /usr/include/boost/ -I $pol/trunk $pol/bin/libpolito.a -o test_politopix
// cygwin: g++ -W  -Wall  -Wextra  -Wswitch  -Wformat  -Wchar-subscripts  -Wparentheses  -Wmultichar  -Wtrigraphs  -Wpointer-arith  -Wcast-align  -Wreturn-type  -Wshadow  -Wundef  -Woverloaded-virtual  -Wno-unused-function -O3 -fno-strict-aliasing test_politopix.cpp  -I /usr/include/boost/ -I $pol/trunk $pol/bin/libpolito.dll.a -o test_politopix
// time ./test_politopix.exe 

int main() {


  void callPolitopix();

  clock_t start = clock();
  for (int i=0; i<100;i++){
	  callPolitopix();
  }
  clock_t end = clock();
  float seconds = (float)(end - start) / CLOCKS_PER_SEC;
  std::cout<< "seconds: "<<seconds <<std::endl;
  return 0;
}

void callPolitopix(){
	  unsigned int dimension=3;
	  Rn::setDimension(dimension);
	  Rn::setTolerance(0.0001);
	  boost::shared_ptr<Polytope_Rn> LF_cone(new Polytope_Rn());
	  boost::shared_ptr<Polytope_Rn> LH_cone(new Polytope_Rn());
	  boost::shared_ptr<Polytope_Rn> RF_cone(new Polytope_Rn());
	  boost::shared_ptr<Polytope_Rn> RH_cone(new Polytope_Rn());
	  boost::shared_ptr<Generator_Rn> VX0(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX1(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX2(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX3(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX4(new Generator_Rn(dimension));

	  //////////////////////////////////////////////////
	  //// Define the friction cones of each leg
	  //////////////////////////////////////////////////

	  //	  First friction cone, Left-Front leg (LF)
	  VX0->setCoordinate(0, 1.0); VX0->setCoordinate(1,1.0);  VX0->setCoordinate(2, 0.6);
	  VX1->setCoordinate(0, 1.0); VX1->setCoordinate(1,-1.0); VX1->setCoordinate(2, 0.6);
	  VX2->setCoordinate(0, -1.0); VX2->setCoordinate(1,1.0); VX2->setCoordinate(2, 1.4);
	  VX3->setCoordinate(0, -1.0); VX3->setCoordinate(1,-1.0);VX3->setCoordinate(2, 1.4);
	  VX4->setCoordinate(0, 0.0);  VX4->setCoordinate(1,0.0);  VX4->setCoordinate(2, 0.0);

	  LF_cone->addGenerator(VX0); LF_cone->addGenerator(VX1); LF_cone->addGenerator(VX2);
	  LF_cone->addGenerator(VX3); LF_cone->addGenerator(VX4);

	  // Compute the double description removing the interior point.
	  politopixAPI::computeDoubleDescription(LF_cone,1000.);

	  //	  Second friction cone, Left-hind leg (LH)
	  VX0->setCoordinate(0, 1.0); VX0->setCoordinate(1,1.0);  VX0->setCoordinate(2, 1.0);
	  VX1->setCoordinate(0, 1.0); VX1->setCoordinate(1,-1.0); VX1->setCoordinate(2, 1.0);
	  VX2->setCoordinate(0, -1.0); VX2->setCoordinate(1,1.0); VX2->setCoordinate(2, 1.0);
	  VX3->setCoordinate(0, -1.0); VX3->setCoordinate(1,-1.0);VX3->setCoordinate(2, 1.0);
	  VX4->setCoordinate(0, 0.0);  VX4->setCoordinate(1,0.0);  VX4->setCoordinate(2, 0.0);

	  LH_cone->addGenerator(VX0); LH_cone->addGenerator(VX1); LH_cone->addGenerator(VX2);
	  LH_cone->addGenerator(VX3); LH_cone->addGenerator(VX4);

	  // Compute the double description removing the interior point.
	  politopixAPI::computeDoubleDescription(LH_cone,1000.);

	  //	  Third friction cone, Right-Front leg (RF)
	  VX0->setCoordinate(0, 1.0); VX0->setCoordinate(1,1.0);  VX0->setCoordinate(2, 1.0);
	  VX1->setCoordinate(0, 1.0); VX1->setCoordinate(1,-1.0); VX1->setCoordinate(2, 1.0);
	  VX2->setCoordinate(0, -1.0); VX2->setCoordinate(1,1.0); VX2->setCoordinate(2, 1.0);
	  VX3->setCoordinate(0, -1.0); VX3->setCoordinate(1,-1.0);VX3->setCoordinate(2, 1.0);
	  VX4->setCoordinate(0, 0.0);  VX4->setCoordinate(1,0.0);  VX4->setCoordinate(2, 0.0);

	  RF_cone->addGenerator(VX0); RF_cone->addGenerator(VX1); RF_cone->addGenerator(VX2);
	  RF_cone->addGenerator(VX3); RF_cone->addGenerator(VX4);

	  // Compute the double description removing the interior point.
	  politopixAPI::computeDoubleDescription(RF_cone,1000.);

	  //	  Second friction cone, Right-hind leg (RH)
	  VX0->setCoordinate(0, 1.0); VX0->setCoordinate(1,1.0);  VX0->setCoordinate(2, 1.0);
	  VX1->setCoordinate(0, 1.0); VX1->setCoordinate(1,-1.0); VX1->setCoordinate(2, 1.0);
	  VX2->setCoordinate(0, -1.0); VX2->setCoordinate(1,1.0); VX2->setCoordinate(2, 1.0);
	  VX3->setCoordinate(0, -1.0); VX3->setCoordinate(1,-1.0);VX3->setCoordinate(2, 1.0);
	  VX4->setCoordinate(0, 0.0);  VX4->setCoordinate(1,0.0);  VX4->setCoordinate(2, 0.0);

	  RH_cone->addGenerator(VX0); RH_cone->addGenerator(VX1); RH_cone->addGenerator(VX2);
	  RH_cone->addGenerator(VX3); RH_cone->addGenerator(VX4);

	  // Compute the double description removing the interior point.
	  politopixAPI::computeDoubleDescription(RH_cone,1000.);

	  //////////////////////////////////////////////////////
	  //	  Compute the Contact Wrench Set (CWS)   ///////
	  //////////////////////////////////////////////////////
	  boost::shared_ptr<Polytope_Rn> CWS(new Polytope_Rn());
	  boost::shared_ptr<Polytope_Rn> CWS_tmp1(new Polytope_Rn());
	  boost::shared_ptr<Polytope_Rn> CWS_tmp2(new Polytope_Rn());
	  politopixAPI::computeMinkowskiSumOfPolytopes(LF_cone, LH_cone, CWS_tmp1);
	  politopixAPI::computeMinkowskiSumOfPolytopes(RF_cone, RH_cone, CWS_tmp2);
	  politopixAPI::computeMinkowskiSumOfPolytopes(CWS_tmp1, CWS_tmp2, CWS);
	  // Dump the result on the standard output.
	  //	  LF_cone->dump(std::cout);

	  unsigned int hs_size = CWS->numberOfHalfSpaces();
	  std::cout<< "num of half spaces: "<<hs_size <<std::endl;
 	  unsigned int gen_size = CWS->numberOfGenerators();
  	  std::cout<< "num of gens: "<<gen_size <<std::endl;
  	  double minimal_distance = 10000.0;
  	  double dist;

  	  //////////////////////////////////////////////////////
  	  //////// set the state of the robot
  	  ////////////////////////////////////////////////////
	  boost::shared_ptr<Generator_Rn> check_point(new Generator_Rn(dimension));
      check_point->setCoordinate(0, 0.0);
      check_point->setCoordinate(1, 0.0);
      check_point->setCoordinate(2, 4.0);

      /////////////////////////////////////////////////////////////////////////////////////////
      ///// check the margin between the state of the robot and the contact wrench set (CWS)
      ///////////////////////////////////////////////////////////////////////////////////////////

  	  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > hs_iterator(CWS->getListOfHalfSpaces());
	  for (hs_iterator.begin(); hs_iterator.end()!=true; hs_iterator.next()) {
		  const boost::shared_ptr<HalfSpace_Rn>& currentHalfSpace = hs_iterator.current();
		  hs_size = currentHalfSpace->dimension();
		  boost::numeric::ublas::vector<double> HScoord = currentHalfSpace->vect();
		  std::cout<<"coord:"<<HScoord <<std::endl;

		  dist = currentHalfSpace->computeDistancePointHyperplane(check_point->vect());
		  std::cout<<"dist: "<<dist <<std::endl;
		  if(dist < minimal_distance){
			  minimal_distance = dist;
		  }
	  }
	  std::cout<<"margin: "<< minimal_distance <<std::endl;
	  if (minimal_distance<0.0){
		  std::cout<<" the robots state is out of the contact wrench set!!!" <<std::endl;
	  }
}
