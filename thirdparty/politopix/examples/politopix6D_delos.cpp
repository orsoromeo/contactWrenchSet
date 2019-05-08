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
#include "politopixAPI.h"
#include "Generator_Rn.h"
#include "DoubleDescription_Rn.h"
#include "Rn.h"
// export PATH=$PATH:$pol/bin
//I run this code with:
//g++ -W  -Wall  -Wextra  -Wswitch  -Wformat  -Wchar-subscripts  -Wparentheses  -Wmultichar  -Wtrigraphs -Wpointer-arith  -Wcast-align  -Wreturn-type  -Wshadow  -Wundef  -Woverloaded-virtual -Wno-unused-function -O3 -fno-strict-aliasing politopix6D.cpp -I /usr/include/boost/ -I $pol/trunk $pol/bin/libpolito.dll.a -o politopix6D
//
// time ./test_politopix6D

int main() {


  void callPolitopix();

  clock_t start = clock();
  for (int i=0; i<1;i++){
    callPolitopix();
  }
  clock_t end = clock();
  float seconds = (float)(end - start) / CLOCKS_PER_SEC;

  return 0;
}

void callPolitopix(){
  unsigned int dimension=6;
  double margin_tol = -10e-4;
  Rn::setDimension(dimension);
  Rn::setTolerance(0.000001);
  boost::shared_ptr<Polytope_Rn> primal_cone(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> primal_polytope(new Polytope_Rn());
  boost::shared_ptr<Generator_Rn> VX00(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX10(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX20(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX30(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX40(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX01(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX11(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX21(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX31(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX41(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX02(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX12(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX22(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX32(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX42(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX03(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX13(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX23(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX33(new Generator_Rn(dimension));
  boost::shared_ptr<Generator_Rn> VX43(new Generator_Rn(dimension));

  //////////////////////////////////////////////////
  //// Define the friction cones of each leg
  //////////////////////////////////////////////////
  VX00->setCoordinate(0, 1.0); VX00->setCoordinate(1,1.0);  VX00->setCoordinate(2, 1.5);  VX00->setCoordinate(3, 1.0);  VX00->setCoordinate(4,1.0);  VX00->setCoordinate(5, 0.5);
  VX10->setCoordinate(0, 1.0); VX10->setCoordinate(1,-1.0); VX10->setCoordinate(2, 2.0);  VX10->setCoordinate(3, 1.0);  VX10->setCoordinate(4,-1.0); VX10->setCoordinate(5, 0.5);
  VX20->setCoordinate(0, -1.0); VX20->setCoordinate(1,1.0); VX20->setCoordinate(2, 2.5);  VX20->setCoordinate(3, -1.0); VX20->setCoordinate(4,1.0);  VX20->setCoordinate(5, 0.5);
  VX30->setCoordinate(0, -1.0); VX30->setCoordinate(1,-1.0);VX30->setCoordinate(2, 1.5);  VX30->setCoordinate(3, -1.0); VX30->setCoordinate(4,-1.0); VX30->setCoordinate(5, 0.5);
  VX01->setCoordinate(0, 1.0); VX01->setCoordinate(1,1.0);  VX01->setCoordinate(2, 3.0);  VX01->setCoordinate(3, 0.0);  VX01->setCoordinate(4,0.0);  VX01->setCoordinate(5, 0.0);
  VX11->setCoordinate(0, 1.0); VX11->setCoordinate(1,-1.0); VX11->setCoordinate(2, 3.0);  VX11->setCoordinate(3, 0.0);  VX11->setCoordinate(4,0.0);  VX11->setCoordinate(5, 0.0);
  VX21->setCoordinate(0, -1.0); VX21->setCoordinate(1,1.0); VX21->setCoordinate(2, 3.7);  VX21->setCoordinate(3, 0.0);  VX21->setCoordinate(4,0.0);  VX21->setCoordinate(5, 0.0);
  VX31->setCoordinate(0, -1.0); VX31->setCoordinate(1,-1.0);VX31->setCoordinate(2, 3.0);  VX31->setCoordinate(3, 0.0);  VX31->setCoordinate(4,0.0);  VX31->setCoordinate(5, 0.0);
  VX02->setCoordinate(0, 1.0); VX02->setCoordinate(1,1.0);  VX02->setCoordinate(2, 1.0);  VX02->setCoordinate(3, 1.0);  VX02->setCoordinate(4,1.0);  VX02->setCoordinate(5, 1.5);
  VX12->setCoordinate(0, 1.0); VX12->setCoordinate(1,-1.0); VX12->setCoordinate(2, 1.5);  VX12->setCoordinate(3, 1.0);  VX12->setCoordinate(4,-1.0); VX12->setCoordinate(5, 1.5);
  VX22->setCoordinate(0, -1.0); VX22->setCoordinate(1,1.0); VX22->setCoordinate(2, 1.0);  VX22->setCoordinate(3, -1.0); VX22->setCoordinate(4,1.0);  VX22->setCoordinate(5, 1.5);
  VX32->setCoordinate(0, -1.0); VX32->setCoordinate(1,-1.0);VX32->setCoordinate(2, 1.2);  VX32->setCoordinate(3, -1.0); VX32->setCoordinate(4,-1.0); VX32->setCoordinate(5, 1.5);
  VX03->setCoordinate(0, 1.0); VX03->setCoordinate(1,1.5);  VX03->setCoordinate(2, 0.5);  VX03->setCoordinate(3, 1.0);  VX03->setCoordinate(4,1.0);  VX03->setCoordinate(5, 1.0);
  VX13->setCoordinate(0, 1.0); VX13->setCoordinate(1,-1.5); VX13->setCoordinate(2, 0.1);  VX13->setCoordinate(3, 1.0);  VX13->setCoordinate(4,-1.0); VX13->setCoordinate(5, 1.0);
  VX23->setCoordinate(0, -1.0); VX23->setCoordinate(1,1.0); VX23->setCoordinate(2, 0.5);  VX23->setCoordinate(3, -1.0); VX23->setCoordinate(4,1.0);  VX23->setCoordinate(5, 1.0);
  VX33->setCoordinate(0, -1.0); VX33->setCoordinate(1,-1.0);VX33->setCoordinate(2, 0.5);  VX33->setCoordinate(3, -1.0); VX33->setCoordinate(4,-1.0); VX33->setCoordinate(5, 1.0);

  boost::shared_ptr<Generator_Rn> origin(new Generator_Rn(dimension));
  origin->setCoordinate(0,0.0);
  origin->setCoordinate(1,0.0);
  origin->setCoordinate(2,0.0);
  origin->setCoordinate(3,0.0);
  origin->setCoordinate(4,0.0);
  origin->setCoordinate(5,0.0);

  primal_polytope->addGenerator(origin);
  primal_polytope->addGenerator(VX00); primal_polytope->addGenerator(VX10); primal_polytope->addGenerator(VX20); primal_polytope->addGenerator(VX30);
  primal_polytope->addGenerator(VX01); primal_polytope->addGenerator(VX11); primal_polytope->addGenerator(VX21); primal_polytope->addGenerator(VX31);
  primal_polytope->addGenerator(VX02); primal_polytope->addGenerator(VX12); primal_polytope->addGenerator(VX22); primal_polytope->addGenerator(VX32);
  primal_polytope->addGenerator(VX03); primal_polytope->addGenerator(VX13); primal_polytope->addGenerator(VX23); primal_polytope->addGenerator(VX33);
	  
  // Compute the double description removing the interior point.
  politopixAPI::computeDoubleDescription(primal_polytope,1000.);

//  primal_polytope->dump(std::cout);

  string path_test1("./");
  boost::shared_ptr<Polytope_Rn> _polytopeTest(new Polytope_Rn());

//  politopixAPI::savePolytope(path_test1 + string("primal_polytope.ptop"), primal_polytope);
//  politopixAPI::loadPolytope( std::string("./primal_polytope.ptop"), primal_polytope);

  //////// set the state of the robot
  ////////////////////////////////////////////////////
//  boost::shared_ptr<Generator_Rn> check_point(new Generator_Rn(dimension));
//  check_point->setCoordinate(0,0.0);
//  check_point->setCoordinate(1,0.0);
//  check_point->setCoordinate(2,9.81);
//  check_point->setCoordinate(3,0.0);
//  check_point->setCoordinate(4,0.0);
//  check_point->setCoordinate(5,9.81);

  /////////////////////////////////////////////////////////////////////////////////////////
  ///// check the margin between the state of the robot and the contact wrench set (CWS)
  ///////////////////////////////////////////////////////////////////////////////////////////
  unsigned int hs_size = primal_polytope->numberOfHalfSpaces();
  std::cout<< "num of half spaces: "<<hs_size <<std::endl;
	  unsigned int gen_size = primal_polytope->numberOfGenerators();
	  std::cout<< "num of gens: "<<gen_size <<std::endl;
	  double minimal_distance = 10e10;
	  double dist;
	  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > hs_iterator(primal_polytope->getListOfHalfSpaces());
  for (hs_iterator.begin(); hs_iterator.end()!=true; hs_iterator.next()) {
	  const boost::shared_ptr<HalfSpace_Rn>& currentHalfSpace = hs_iterator.current();
	  hs_size = currentHalfSpace->dimension();
	  boost::numeric::ublas::vector<double> HScoord = currentHalfSpace->vect();
	  std::cout<<"coord:"<<HScoord <<std::endl;

//	  dist = currentHalfSpace->computeDistancePointHyperplane(check_point->vect());
//	  std::cout<<"dist: "<<dist <<std::endl;
//	  if(dist < minimal_distance){
//		  minimal_distance = dist;
//	  }
  }
  std::cout<<"margin: "<< minimal_distance <<std::endl;
  if (minimal_distance<margin_tol){
	  std::cout<<" the robots state is out of the contact wrench set!!!" <<std::endl;
  }

}
