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

  return 0;
}

void callPolitopix(){
	  unsigned int dimension=4;
	  Rn::setDimension(dimension);
	  Rn::setTolerance(0.00001);
	  boost::shared_ptr<Polytope_Rn> Cube(new Polytope_Rn());
	  boost::shared_ptr<Generator_Rn> VX0(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX1(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX2(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX3(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX4(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX5(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX6(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX7(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX8(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX9(new Generator_Rn(dimension));
	  boost::shared_ptr<Generator_Rn> VX10(new Generator_Rn(dimension));
	  VX0->setCoordinate(0, 1.0); VX0->setCoordinate(1,1.0); VX0->setCoordinate(2, 1.0);
	  VX1->setCoordinate(0, 0.0); VX1->setCoordinate(1,1.0); VX1->setCoordinate(2, 1.0);
	  VX2->setCoordinate(0, 0.0); VX2->setCoordinate(1,0.0); VX2->setCoordinate(2, 1.0);
	  VX3->setCoordinate(0, 0.0); VX3->setCoordinate(1,0.0); VX3->setCoordinate(2, 1.0);
	  VX4->setCoordinate(0, 0.0); VX4->setCoordinate(1,0.0); VX4->setCoordinate(2, 0.0);
	  VX5->setCoordinate(0, 0.0); VX5->setCoordinate(1,0.0); VX5->setCoordinate(2, 0.0);
	  VX6->setCoordinate(0, 0.0); VX6->setCoordinate(1,0.0); VX6->setCoordinate(2, 0.0);
	  VX7->setCoordinate(0, 0.0); VX7->setCoordinate(1,0.0); VX7->setCoordinate(2, 0.0);
	  VX8->setCoordinate(0, 0.0); VX8->setCoordinate(1,0.0); VX8->setCoordinate(2, 0.0);
	  VX9->setCoordinate(0, 0.5); VX9->setCoordinate(1,0.5); VX9->setCoordinate(2, 0.5);
	  VX10->setCoordinate(0, 0.5); VX10->setCoordinate(1,0.5); VX10->setCoordinate(2, 0.5);

	  VX0->setCoordinate(3, 1.0);
	  VX1->setCoordinate(3, 1.0);
	  VX2->setCoordinate(3, 1.0);
	  VX3->setCoordinate(3, 1.0);
	  VX4->setCoordinate(3, 1.0);
	  VX5->setCoordinate(3, 0.0);
	  VX6->setCoordinate(3, 0.0);
	  VX7->setCoordinate(3, 0.0);
	  VX8->setCoordinate(3, 0.0);
	  VX9->setCoordinate(3, 0.5);
	  VX10->setCoordinate(3, 0.5);

	  Cube->addGenerator(VX0); Cube->addGenerator(VX1); Cube->addGenerator(VX2);
	  Cube->addGenerator(VX3); Cube->addGenerator(VX4); Cube->addGenerator(VX5);
	  Cube->addGenerator(VX6); Cube->addGenerator(VX7); Cube->addGenerator(VX8);
	  Cube->addGenerator(VX9); Cube->addGenerator(VX10);

	  // Compute the double description removing the interior point.
	  politopixAPI::computeDoubleDescription(Cube,1000.);
	  // Dump the result on the standard output.
//	  Cube->dump(std::cout);

	  unsigned int gen_size = Cube->numberOfGenerators();
	  std::cout<< "num of gens: "<<gen_size <<std::endl;
	  listOfGeometricObjects< boost::shared_ptr<Generator_Rn> > listOfEdges = Cube->getListOfGenerators();
	  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > gen_iterator(listOfEdges);
	  for (gen_iterator.begin(); gen_iterator.end()!=true; gen_iterator.next()) {
		  const boost::shared_ptr<Generator_Rn>& currentGenerator = gen_iterator.current();
		  gen_size = currentGenerator->dimension();
		  boost::numeric::ublas::vector<double> GenCoord = currentGenerator->vect();
		  std::cout<< "gen size: "<<gen_size<<" coord:"<<GenCoord <<std::endl;
	  }

	  unsigned int hs_size = Cube->numberOfHalfSpaces();
	  std::cout<< "num of half spaces: "<<hs_size <<std::endl;
	  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > hs_iterator(Cube->getListOfHalfSpaces());
	  for (hs_iterator.begin(); hs_iterator.end()!=true; hs_iterator.next()) {
		  const boost::shared_ptr<HalfSpace_Rn>& currentHalfSpace = hs_iterator.current();
		  hs_size = currentHalfSpace->dimension();
		  boost::numeric::ublas::vector<double> HScoord = currentHalfSpace->vect();
		  std::cout<<"coord:"<<HScoord <<std::endl;
	  }
//
//	    
//	  // Save on file.
//	  politopixAPI::savePolytope(std::string("thisCube.ptop"), Cube);
}

  // {
  //   boost::timer this_timer100;
  //   std::vector< boost::shared_ptr<Polytope_Rn> > allVoronoiCells_Angers100;
  //   politopixAPI::computeVoronoiDiagram(inputSpaceAngers, listOfPoints3DAngers, allVoronoiCells_Angers100, Voronoi_Rn::CellByCellDistNgb);
  //   cout << "TIME100=" << this_timer100.elapsed() << std::endl;
  //   for (unsigned int i=0; i<allVoronoiCells_Angers100.size(); ++i) {
  //     std::ostringstream stream_;
  //     stream_ << "test/CellByCellDistNgb_cellVor";
  //     stream_ << i;
  //     stream_ << ".ptop";
  //     std::string valString = stream_.str();
  //     //std::cout << valString << std::endl;
  //     politopixAPI::savePolytope(vorPathLascaux + valString, allVoronoiCells_Angers100[i]);
  //   }
  // }

  // {
  //   boost::timer this_timer40;
  //   std::vector< boost::shared_ptr<Polytope_Rn> > allVoronoiCells_Angers40;
  //   politopixAPI::computeVoronoiDiagram(inputSpaceAngers, listOfPoints3DAngers, allVoronoiCells_Angers40, Voronoi_Rn::CellByCellDistNgb_40pc);
  //   cout << "TIME40  =" << this_timer40.elapsed() << std::endl;
  //   for (unsigned int i=0; i<allVoronoiCells_Angers40.size(); ++i) {
  //     std::ostringstream stream_;
  //     stream_ << "test/CellByCellDistNgb_40pc_cellVor";
  //     stream_ << i;
  //     stream_ << ".ptop";
  //     std::string valString = stream_.str();
  //     //std::cout << valString << std::endl;
  //     politopixAPI::savePolytope(vorPathLascaux + valString, allVoronoiCells_Angers40[i]);
  //   }
  // }

  // {
  //   boost::timer this_timer30;
  //   std::vector< boost::shared_ptr<Polytope_Rn> > allVoronoiCells_Angers30;
  //   politopixAPI::computeVoronoiDiagram(inputSpaceAngers, listOfPoints3DAngers, allVoronoiCells_Angers30, Voronoi_Rn::CellByCellDistNgb_30pc);
  //   cout << "TIME30 =" << this_timer30.elapsed() << std::endl;
  //   for (unsigned int i=0; i<allVoronoiCells_Angers30.size(); ++i) {
  //     std::ostringstream stream_;
  //     stream_ << "test/CellByCellDistNgb_30pc_cellVor";
  //     stream_ << i;
  //     stream_ << ".ptop";
  //     std::string valString = stream_.str();
  //     //std::cout << valString << std::endl;
  //     politopixAPI::savePolytope(vorPathLascaux + valString, allVoronoiCells_Angers30[i]);
  //   }
  // }

  // {
  //   boost::timer this_timer10;
  //   std::vector< boost::shared_ptr<Polytope_Rn> > allVoronoiCells_Angers10;
  //   politopixAPI::computeVoronoiDiagram(inputSpaceAngers, listOfPoints3DAngers, allVoronoiCells_Angers10, Voronoi_Rn::CellByCellDistNgb_10pc);
  //   cout << "TIME10 =" << this_timer10.elapsed() << std::endl;
  //   for (unsigned int i=0; i<allVoronoiCells_Angers10.size(); ++i) {
  //     std::ostringstream stream_;
  //     stream_ << "test/CellByCellDistNgb_10pc_cellVor";
  //     stream_ << i;
  //     stream_ << ".ptop";
  //     std::string valString = stream_.str();
  //     //std::cout << valString << std::endl;
  //     politopixAPI::savePolytope(vorPathLascaux + valString, allVoronoiCells_Angers10[i]);
  //   }
  // }

