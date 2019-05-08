// test_politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2014-2016 : Delos Vincent
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
//
/// \file test_politopix.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)

// g++ -W  -Wall  -Wextra  -Wswitch  -Wformat  -Wchar-subscripts  -Wparentheses  -Wmultichar  -Wtrigraphs  -Wpointer-arith  -Wcast-align  -Wreturn-type  -Wshadow  -Wundef  -Woverloaded-virtual  -Wno-unused-function -O3 -fno-strict-aliasing test_politopix.cpp  -I /usr/include/boost/ -I $pol/trunk $pol/build/libpolito.a -o test_politopix
// To run the tests (the last one takes far longer):
// ./test_politopix.exe
// ./test_politopix.exe piston

#include <boost/test/minimal.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/timer.hpp>
#include <iostream>
#include <string.h>
#include <vector>
#include <cmath>
#include "politopixAPI.h"
#include "NormalFan_Rn.h"

using namespace std;


int test_main(int argc, char* argv[]) {

  Rn::setDimension(6);
  Rn::setTolerance(0.000001);

  cout << endl;
  cout << "////////////" << endl;
  cout << "// DATA1 //" << endl;
  cout << "//////////" << endl;
  string path_test1("test/DATA1/");
  boost::shared_ptr<Polytope_Rn> _polytopeTest(new Polytope_Rn());
  cout << "######################" << endl;
  cout << "# DATA1: truncations #" << endl;
  cout << "######################" << endl;
  Rn::setDimension(3);
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("lmp1.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 6 );
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("lmp2.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 6 );
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("lmp3.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 12 );
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("SharirCube.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 8 );
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("CUT3-3-4.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 4 );
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("CR0-3-6.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 6 );
  Rn::setDimension(4);
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("neighborly_4_8.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 8 );
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("cyclic.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 10000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 8 );
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("BIR3-4-6.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 6 );
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("p4_7_simplicial_1.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 7 );
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("p4_7_nonsimplicial_01.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 7 );
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("120-cell.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 600 );
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("not-cubical4.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 32 );
  Rn::setDimension(6);
  _polytopeTest.reset(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> _polytope2(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("8D6.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("8D6_v.ptop"), _polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(_polytopeTest, _polytope2) == TEST_OK );
  _polytope2.reset(new Polytope_Rn());
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("10D6.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("10D6_v.ptop"), _polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(_polytopeTest, _polytope2) == TEST_OK );
  _polytope2.reset(new Polytope_Rn());
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("20D6.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("20D6_v.ptop"), _polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(_polytopeTest, _polytope2) == TEST_OK );
  Rn::setDimension(10);
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("CF-10-11.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  cout << "########################" << endl;
  cout << "# DATA1: intersections #" << endl;
  cout << "########################" << endl;
  Rn::setDimension(3);
  _polytope2.reset(new Polytope_Rn());
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("lmp1_v.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("SharirCube_v.ptop"), _polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeIntersection(_polytopeTest, _polytope2) == TEST_OK );
  BOOST_REQUIRE( _polytopeTest->numberOfGenerators() == 8 );
  cout << "###############" << endl;
  cout << "# DATA1: sums #" << endl;
  cout << "###############" << endl;
  Rn::setDimension(3);
  Rn::setTolerance(0.0001);
  boost::shared_ptr<Polytope_Rn> _sum;
  _sum.reset(new Polytope_Rn());
  _polytope2.reset(new Polytope_Rn());
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("lmp1_v.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("cube3D.ptop"), _polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytopeTest, _polytope2, _sum) == TEST_OK );
  BOOST_REQUIRE( _sum->numberOfGenerators() == 24 );
  BOOST_REQUIRE( _sum->numberOfHalfSpaces() == 26 );
  _sum.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytope2, _polytopeTest, _sum) == TEST_OK );
  BOOST_REQUIRE( _sum->numberOfGenerators() == 24 );
  BOOST_REQUIRE( _sum->numberOfHalfSpaces() == 26 );
  _sum.reset(new Polytope_Rn());
  _polytope2.reset(new Polytope_Rn());
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("cyclic_v.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("120-cell_v.ptop"), _polytope2) == TEST_OK );
  //BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytopeTest, _polytope2, _sum) == TEST_OK );
  _polytope2.reset(new Polytope_Rn());
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("cube_1.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("cube_2.ptop"), _polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytope2, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(_polytope2, _polytopeTest, true) == TEST_OK );
  //BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytope2, _polytopeTest, _sum) == TEST_OK );
  Rn::setTolerance(0.000001);
  _sum.reset(new Polytope_Rn());
  _polytope2.reset(new Polytope_Rn());
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("tetra3D_v_1.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("tetra3D_v_2.ptop"), _polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytopeTest, _polytope2, _sum) == TEST_OK );
  BOOST_REQUIRE( _sum->numberOfGenerators() == 9 );
  BOOST_REQUIRE( _sum->numberOfHalfSpaces() == 7 );
  _sum.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytope2, _polytopeTest, _sum) == TEST_OK );
  BOOST_REQUIRE( _sum->numberOfGenerators() == 9 );
  BOOST_REQUIRE( _sum->numberOfHalfSpaces() == 7 );
  _sum.reset(new Polytope_Rn());
  _polytope2.reset(new Polytope_Rn());
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("tetra3D_v_3.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("tetra3D_v_4.ptop"), _polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytopeTest, _polytope2, _sum) == TEST_OK );
  BOOST_REQUIRE( _sum->numberOfGenerators() == 17 );
  BOOST_REQUIRE( _sum->numberOfHalfSpaces() == 17 );
  _sum.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytope2, _polytopeTest, _sum) == TEST_OK );
  BOOST_REQUIRE( _sum->numberOfGenerators() == 17 );
  BOOST_REQUIRE( _sum->numberOfHalfSpaces() == 17 );
  Rn::setDimension(4);
  Rn::setTolerance(0.0001);
  _sum.reset(new Polytope_Rn());
  _polytope2.reset(new Polytope_Rn());
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("_p4s.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("_p4ns.ptop"), _polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytopeTest, _polytope2, _sum) == TEST_OK );
  BOOST_REQUIRE( _sum->numberOfGenerators() == 25 );
  BOOST_REQUIRE( _sum->numberOfHalfSpaces() == 33 );
  _sum.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytope2, _polytopeTest, _sum) == TEST_OK );
  BOOST_REQUIRE( _sum->numberOfGenerators() == 25 );
  BOOST_REQUIRE( _sum->numberOfHalfSpaces() == 33 );
  Rn::setTolerance(0.000001);
  _sum.reset(new Polytope_Rn());
  _polytope2.reset(new Polytope_Rn());
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("neighborly_4_8_v.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("BIR3-4-6_v.ptop"), _polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytopeTest, _polytope2, _sum) == TEST_OK );
  BOOST_REQUIRE( _sum->numberOfGenerators() == 34 );
  BOOST_REQUIRE( _sum->numberOfHalfSpaces() == 40 );
  _sum.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytope2, _polytopeTest, _sum) == TEST_OK );
  BOOST_REQUIRE( _sum->numberOfGenerators() == 34 );
  BOOST_REQUIRE( _sum->numberOfHalfSpaces() == 40 );
  Rn::setDimension(6);
  Rn::setTolerance(0.000001);
  _sum.reset(new Polytope_Rn());
  _polytope2.reset(new Polytope_Rn());
  _polytopeTest.reset(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> _pol2compare;
  _pol2compare.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("cube15pts.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("sph10D6.ptop"), _polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(_polytopeTest, _polytope2, _sum) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("cube15-sph10_D6.ptop"), _pol2compare) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(_sum, _pol2compare) == TEST_OK );

  cout << endl;
  cout << "/////////////" << endl;
  cout << "// DATA 2 //" << endl;
  cout << "///////////" << endl;
  Rn::setTolerance(0.000001);
  Rn::setDimension(3);
  path_test1 = string("test/DATA2/");
  boost::shared_ptr<Polytope_Rn> v2h(new Polytope_Rn()), comp_v2h(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("v_lmp1.ptop"), v2h) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(v2h, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("comp_v_lmp1.ptop"), comp_v2h) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(v2h, comp_v2h) == TEST_OK );
  v2h.reset(new Polytope_Rn()); comp_v2h.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("v_lmp2.ptop"), v2h) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(v2h, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("comp_v_lmp2.ptop"), comp_v2h) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(v2h, comp_v2h) == TEST_OK );
  v2h.reset(new Polytope_Rn()); comp_v2h.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("v_lmp3.ptop"), v2h) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(v2h, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("comp_v_lmp3.ptop"), comp_v2h) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(v2h, comp_v2h) == TEST_OK );

  Rn::setDimension(6);
  Rn::setTolerance(0.000001);
  std::vector< boost::shared_ptr<Polytope_Rn> > arrayOfPolytopes;
  arrayOfPolytopes.resize(7);
  path_test1 = string("test/DATA2/T1/");
  string path_test2("test/DATA2/T2/");
  string path_test3("test/DATA2/T3/");
  string path_test4("test/DATA2/T4/");
  string path_test5("test/DATA2/T5/");
  std::vector< string > Ppath1;
  std::vector< string > Ppath2;
  std::vector< string > Ppath3;
  std::vector< string > Ppath4;
  std::vector< string > Ppath5;
  //cout << "ok" << endl;
  Ppath1.push_back( string("outputcthp_0_1-ctl_848_0.ptop") );
  Ppath1.push_back( string("outputcthp_10_0-ctl_524_1.ptop") );
  Ppath1.push_back( string("outputcthp_102_1-ctl_221_0.ptop") );
  Ppath1.push_back( string("outputcthp_103_0-ctl_122_2.ptop") );
  Ppath1.push_back( string("outputcthp_107_3-ctl_145_1.ptop") );
  Ppath1.push_back( string("outputcthp_109_2-ctl_452_2.ptop") );
  Ppath1.push_back( string("outputcthp_113_0-ctl_989_2.ptop") );
  Ppath2.push_back( string("outputcthp_0_2-ctl_82_2.ptop") );
  Ppath2.push_back( string("outputcthp_10_3-ctl_99_3.ptop") );
  Ppath2.push_back( string("outputcthp_100_0-ctl_914_1.ptop") );
  Ppath2.push_back( string("outputcthp_104_2-ctl_857_3.ptop") );
  Ppath2.push_back( string("outputcthp_108_3-ctl_196_1.ptop") );
  Ppath2.push_back( string("outputcthp_109_2-ctl_452_2.ptop") );
  Ppath2.push_back( string("outputcthp_110_1-ctl_697_0.ptop") );
  Ppath3.push_back( string("outputcthp_111_1-ctl_111_1.ptop") );
  Ppath3.push_back( string("outputcthp_162_0-ctl_42_1.ptop") );
  Ppath3.push_back( string("outputcthp_366_1-ctl_201_0.ptop") );
  Ppath3.push_back( string("outputcthp_57_3-ctl_330_2.ptop") );
  Ppath3.push_back( string("outputcthp_586_0-ctl_932_2.ptop") );
  Ppath3.push_back( string("outputcthp_741_3-ctl_568_0.ptop") );
  Ppath3.push_back( string("outputcthp_915_3-ctl_803_1.ptop") );
  Ppath4.push_back( string("outputcthp_177_1-ctl_83_2.ptop") );
  Ppath4.push_back( string("outputcthp_231_1-ctl_830_2.ptop") );
  Ppath4.push_back( string("outputcthp_361_1-ctl_578_2.ptop") );
  Ppath4.push_back( string("outputcthp_44_1-ctl_51_2.ptop") );
  Ppath4.push_back( string("outputcthp_482_2-ctl_880_0.ptop") );
  Ppath4.push_back( string("outputcthp_712_1-ctl_693_1.ptop") );
  Ppath4.push_back( string("outputcthp_899_3-ctl_233_3.ptop") );
  Ppath5.push_back( string("outputcthp_0_2-ctl_82_2.ptop") );
  Ppath5.push_back( string("outputcthp_10_0-ctl_524_1.ptop") );
  Ppath5.push_back( string("outputcthp_101_0-ctl_74_0.ptop") );
  Ppath5.push_back( string("outputcthp_104_1-ctl_357_2.ptop") );
  Ppath5.push_back( string("outputcthp_107_2-ctl_186_0.ptop") );
  Ppath5.push_back( string("outputcthp_112_2-ctl_527_3.ptop") );
  Ppath5.push_back( string("outputcthp_113_2-ctl_424_1.ptop") );
  //cout << "ok" << endl;
  {for (unsigned int i=0; i<7; ++i) {
    //cout << "ok " << path_test1+Ppath1[i] << endl;
    boost::shared_ptr<Polytope_Rn> T1(new Polytope_Rn());
    boost::shared_ptr<Polytope_Rn> T2(new Polytope_Rn());
    boost::shared_ptr<Polytope_Rn> T3(new Polytope_Rn());
    boost::shared_ptr<Polytope_Rn> T4(new Polytope_Rn());
    boost::shared_ptr<Polytope_Rn> T5(new Polytope_Rn());
    BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + Ppath1[i], T1) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::loadPolytope(path_test2 + Ppath2[i], T2) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + Ppath3[i], T3) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::loadPolytope(path_test4 + Ppath4[i], T4) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::loadPolytope(path_test5 + Ppath5[i], T5) == TEST_OK );
    cout << "##################" << endl;
    cout << "# DATA 2: sums ";
    cout << i+1 << " #" << endl;
    cout << "##################" << endl;
    boost::shared_ptr<Polytope_Rn> T12(new Polytope_Rn());
    boost::shared_ptr<Polytope_Rn> T123(new Polytope_Rn());
    boost::shared_ptr<Polytope_Rn> T1234(new Polytope_Rn());
    boost::shared_ptr<Polytope_Rn> T12345(new Polytope_Rn());
    BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(T1, T2, T12) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(T12, T3, T123) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(T123, T4, T1234) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(T1234, T5, T12345) == TEST_OK );
    boost::shared_ptr<Polytope_Rn> T54(new Polytope_Rn());
    boost::shared_ptr<Polytope_Rn> T543(new Polytope_Rn());
    boost::shared_ptr<Polytope_Rn> T5432(new Polytope_Rn());
    boost::shared_ptr<Polytope_Rn> T54321(new Polytope_Rn());
    BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(T5, T4, T54) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(T54, T3, T543) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(T543, T2, T5432) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(T5432, T1, T54321) == TEST_OK );
    //BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(T12345, T54321) == TEST_OK );
    arrayOfPolytopes[i] = T12345;
  }}


  Rn::setDimension(6);
  cout << endl;
  cout << "////////////////" << endl;
  cout << "// SSS DISC1 //" << endl;
  cout << "//////////////" << endl;
  path_test1 = string("test/TM/Sss1_disc6/");
  arrayOfPolytopes.clear();
  arrayOfPolytopes.resize(12);
  cout << "##########################" << endl;
  cout << "# SSS DISC1: truncations #" << endl;
  cout << "##########################" << endl;
  {for (unsigned int i=0; i<12; ++i) {
    ostringstream ss;
    ss << (i+1);
    string p_name = string("Polytope_") + ss.str() + string(".ptop");
    arrayOfPolytopes[i].reset(new Polytope_Rn());
    BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + p_name, arrayOfPolytopes[i]) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::computeDoubleDescription(arrayOfPolytopes[i], 1000.) == TEST_OK );
  }}
  cout << "############################" << endl;
  cout << "# SSS DISC1: intersections #" << endl;
  cout << "############################" << endl;
  {for (unsigned int i=1; i<5 ; ++i) {
    BOOST_REQUIRE( politopixAPI::computeIntersection(arrayOfPolytopes[0], arrayOfPolytopes[i]) == TEST_OK );
  }}
  {for (unsigned int i=6; i<10; ++i) {
    BOOST_REQUIRE( politopixAPI::computeIntersection(arrayOfPolytopes[5], arrayOfPolytopes[i]) == TEST_OK );
  }}
  cout << "# intersections 11 & 12 #" << endl;
  BOOST_REQUIRE( politopixAPI::computeIntersection(arrayOfPolytopes[0], arrayOfPolytopes[10]) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeIntersection(arrayOfPolytopes[5], arrayOfPolytopes[11]) == TEST_OK );
  cout << "###################" << endl;
  cout << "# SSS DISC1: sums #" << endl;
  cout << "###################" << endl;
  boost::shared_ptr<Polytope_Rn> C05(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> C50(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(arrayOfPolytopes[0], arrayOfPolytopes[5], C05) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(arrayOfPolytopes[5], arrayOfPolytopes[0], C50) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(C05, C50) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::savePolytope(path_test1 + string("test_resultat12.ptop"), C05) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::savePolytope(path_test1 + string("test_resultat21.ptop"), C50) == TEST_OK );


  cout << endl;
  cout << "//////////////////" << endl;
  cout << "// BBB_D_5.966 //" << endl;
  cout << "////////////////" << endl;
  path_test2 = string("test/TM/Interference_Bbb/BBB_D_5.966/");
  arrayOfPolytopes.clear();
  arrayOfPolytopes.resize(13);
  {for (unsigned int i=1; i<=13; ++i) {
    cout << "############################" << endl;
    cout << "# BBB_D_5.966: truncations #" << endl;
    cout << "############################" << endl;
    ostringstream ss1;
    ss1 << i;
    std::vector< boost::shared_ptr<Polytope_Rn> > tmpArrayOfPtop(5);
    {for (unsigned int j=1; j<=5; ++j) {
      ostringstream ss2;
      ss2 << j;
      string p_name = string("ctl_") + ss1.str();
      p_name += string("_") + ss2.str() + string(".ptop");
      tmpArrayOfPtop[j-1].reset(new Polytope_Rn());
      BOOST_REQUIRE( politopixAPI::loadPolytope(path_test2 + p_name,  tmpArrayOfPtop[j-1]) == TEST_OK );
      BOOST_REQUIRE( politopixAPI::computeDoubleDescription(tmpArrayOfPtop[j-1], 1000.) == TEST_OK );
    }}
    cout << "##############################" << endl;
    cout << "# BBB_D_5.966: intersections #" << endl;
    cout << "##############################" << endl;
    {for (unsigned int j=1; j<=4; ++j) {
      BOOST_REQUIRE( politopixAPI::computeIntersection(tmpArrayOfPtop[0], tmpArrayOfPtop[j]) == TEST_OK );
    }}
    arrayOfPolytopes.push_back( tmpArrayOfPtop[0] );
    string checkPtop;
    checkPtop += path_test2 + string("qhull_ctl_");
    checkPtop += ss1.str() + string("_x.ptop");
    boost::shared_ptr<Polytope_Rn> CP(new Polytope_Rn());
    BOOST_REQUIRE( politopixAPI::loadPolytope(checkPtop, CP) == TEST_OK );
    BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(tmpArrayOfPtop[0], CP) == TEST_OK );
  }}


  cout << endl;
  cout << "//////////////////" << endl;
  cout << "// BBB_D_5.976 //" << endl;
  cout << "////////////////" << endl;
  path_test2 = string("test/TM/Interference_Bbb/BBB_D_5.976/");
  arrayOfPolytopes.clear();
  arrayOfPolytopes.resize(13);
  {for (unsigned int i=1; i<=13; ++i) {
    cout << "############################" << endl;
    cout << "# BBB_D_5.976: truncations #" << endl;
    cout << "############################" << endl;
    ostringstream ss1;
    ss1 << i;
    std::vector< boost::shared_ptr<Polytope_Rn> > tmpArrayOfPtop(5);
    {for (unsigned int j=1; j<=5; ++j) {
      ostringstream ss2;
      ss2 << j;
      string p_name = string("ctl_") + ss1.str();
      p_name += string("_") + ss2.str() + string(".ptop");
      tmpArrayOfPtop[j-1].reset(new Polytope_Rn());
      BOOST_REQUIRE( politopixAPI::loadPolytope(path_test2 + p_name, tmpArrayOfPtop[j-1]) == TEST_OK );
      BOOST_REQUIRE( politopixAPI::computeDoubleDescription(tmpArrayOfPtop[j-1], 1000.) == TEST_OK );
    }}
    cout << "##############################" << endl;
    cout << "# BBB_D_5.976: intersections #" << endl;
    cout << "##############################" << endl;
    {for (unsigned int j=1; j<=4; ++j) {
      BOOST_REQUIRE( politopixAPI::computeIntersection(tmpArrayOfPtop[0], tmpArrayOfPtop[j]) == TEST_OK );
    }}
    arrayOfPolytopes.push_back( tmpArrayOfPtop[0] );
    if (i != 10) {
      string checkPtop;
      checkPtop += path_test2 + string("qhull_ctl_");
      checkPtop += ss1.str() + string("_x.ptop");
      boost::shared_ptr<Polytope_Rn> CP(new Polytope_Rn());
      BOOST_REQUIRE( politopixAPI::loadPolytope(checkPtop, CP) == TEST_OK );
      BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(tmpArrayOfPtop[0], CP) == TEST_OK );
    }
  }}


  cout << endl;
  cout << "/////////////////////" << endl;
  cout << "// DTL_BPTL/D3-R6 //" << endl;
  cout << "///////////////////" << endl;
  path_test3 = string("test/TM/Interference_Bbb/DTL_BPTL/D3-R6/");
  arrayOfPolytopes.resize(4);
  arrayOfPolytopes[0].reset(new Polytope_Rn()); arrayOfPolytopes[1].reset(new Polytope_Rn());
  arrayOfPolytopes[2].reset(new Polytope_Rn()); arrayOfPolytopes[3].reset(new Polytope_Rn());
  cout << "###############################" << endl;
  cout << "# DTL_BPTL/D3-R6: truncations #" << endl;
  cout << "###############################" << endl;
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("BPHP_BPTL.ptop"), arrayOfPolytopes[0]) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(arrayOfPolytopes[0], 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("SPAR_BPHP.ptop"), arrayOfPolytopes[1]) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(arrayOfPolytopes[1], 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("CCCTL_CCCTHP_SPAR.ptop"), arrayOfPolytopes[2]) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(arrayOfPolytopes[2], 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("DTL2_CCCTL.ptop"), arrayOfPolytopes[3]) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(arrayOfPolytopes[3], 1000.) == TEST_OK );
  boost::shared_ptr<Polytope_Rn> BPHP_BPTL_v(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> SPAR_BPHP_v(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> CCCTL_CCCTHP_SPAR_v(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> DTL2_CCCTL_v(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("BPHP_BPTL_v.ptop"), BPHP_BPTL_v) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("SPAR_BPHP_v.ptop"), SPAR_BPHP_v) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("CCCTL_CCCTHP_SPAR_v.ptop"), CCCTL_CCCTHP_SPAR_v) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("DTL2_CCCTL_v.ptop"), DTL2_CCCTL_v) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(arrayOfPolytopes[0], BPHP_BPTL_v) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(arrayOfPolytopes[1], SPAR_BPHP_v) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(arrayOfPolytopes[2], CCCTL_CCCTHP_SPAR_v) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(arrayOfPolytopes[3], DTL2_CCCTL_v) == TEST_OK );
  cout << "#####################################" << endl;
  cout << "# DTL_BPTL/D3-R6: polar truncations #" << endl;
  cout << "#####################################" << endl;
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("BPHP_BPTL.ptop"), arrayOfPolytopes[0]) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("SPAR_BPHP.ptop"), arrayOfPolytopes[1]) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("CCCTL_CCCTHP_SPAR.ptop"), arrayOfPolytopes[2]) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("DTL2_CCCTL.ptop"), arrayOfPolytopes[3]) == TEST_OK );
  {for (unsigned int i=0; i<4; ++i) {
    boost::shared_ptr<Polytope_Rn> thisVPol(new Polytope_Rn());
    // Get the V-description.
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(arrayOfPolytopes[i]->getListOfGenerators());
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
      thisVPol->addGenerator(iteGN.current());
    }}
    BOOST_REQUIRE( politopixAPI::computeDoubleDescription(thisVPol, 1000.) == TEST_OK );
    arrayOfPolytopes[i] = thisVPol;
  }}
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(arrayOfPolytopes[0], BPHP_BPTL_v) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(arrayOfPolytopes[1], SPAR_BPHP_v) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(arrayOfPolytopes[2], CCCTL_CCCTHP_SPAR_v) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(arrayOfPolytopes[3], DTL2_CCCTL_v) == TEST_OK );
  cout << "########################" << endl;
  cout << "# DTL_BPTL/D3-R6: sums #" << endl;
  cout << "########################" << endl;
  boost::shared_ptr<Polytope_Rn> BPTL_SPAR12(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> BPTL_SPAR21(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> BPTL_CCCTL12(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> BPTL_CCCTL21(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> BPTL_DTL212(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> BPTL_DTL221(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(arrayOfPolytopes[0], arrayOfPolytopes[1], BPTL_SPAR12) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(arrayOfPolytopes[1], arrayOfPolytopes[0], BPTL_SPAR21) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(arrayOfPolytopes[2], BPTL_SPAR12, BPTL_CCCTL12) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(BPTL_SPAR21, arrayOfPolytopes[2], BPTL_CCCTL21) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(arrayOfPolytopes[3], BPTL_CCCTL12, BPTL_DTL212) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(BPTL_CCCTL21, arrayOfPolytopes[3], BPTL_DTL221) == TEST_OK );
  boost::shared_ptr<Polytope_Rn> qhl_BPTL_SPAR12(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> qhl_BPTL_CCCTL21(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> qhl_BPTL_DTL221(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("qhl_BPTL_SPAR12.ptop"), qhl_BPTL_SPAR12) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("qhl_BPTL_CCCTL21.ptop"), qhl_BPTL_CCCTL21) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test3 + string("qhl_BPTL_DTL221.ptop"), qhl_BPTL_DTL221) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(BPTL_SPAR12, BPTL_SPAR21) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(BPTL_CCCTL12, BPTL_CCCTL21) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(BPTL_DTL212, BPTL_DTL221) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(BPTL_SPAR12, qhl_BPTL_SPAR12) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(BPTL_CCCTL12, qhl_BPTL_CCCTL21) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(BPTL_DTL221, qhl_BPTL_DTL221) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::savePolytope(path_test3 + string("test_BPTL_DTL221.ptop"), BPTL_DTL221) == TEST_OK );
  Rn::setDimension(3);
  boost::shared_ptr<Polytope_Rn> pol2polarize2(new Polytope_Rn()), polarizedPol2(new Polytope_Rn()), polarizedPolComp2(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope("test/DATA1/outputpolar.ptop", pol2polarize2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope("test/DATA1/newpol.ptop", polarizedPolComp2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::PolarPolytope(pol2polarize2, polarizedPol2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkTopologyAndGeometry(polarizedPol2) == TEST_OK );
  //BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(polarizedPol2, polarizedPolComp2) == TEST_OK );

  cout << endl;
  cout << "//////////////////" << endl;
  cout << "// SPECTOMETER //" << endl;
  cout << "////////////////" << endl;
  cout << "#########" << endl;
  cout << "# SUM 1 #" << endl;
  cout << "#########" << endl;
  string SPEC_path_test = string("test/SPECTRO/");
  Rn::setDimension(6);
  Rn::setTolerance(0.000001);
  boost::shared_ptr<Polytope_Rn> SPEC_sum12(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> SPEC_sum21(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> SPEC_polytope1(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> SPEC_polytope2(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> SPEC_polytope2test(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(SPEC_path_test + string("P_10_20_FC2.ptop"), SPEC_polytope1) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(SPEC_path_test + string("P_20_30_FC2.ptop"), SPEC_polytope2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(SPEC_path_test + string("comp_10_20_30.ptop"), SPEC_polytope2test) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(SPEC_polytope1, SPEC_polytope2, SPEC_sum12) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(SPEC_polytope2, SPEC_polytope1, SPEC_sum21) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(SPEC_sum12, SPEC_polytope2test) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(SPEC_sum21, SPEC_polytope2test) == TEST_OK );

  cout << endl;
  cout << "//////////////////////" << endl;
  cout << "// COMPUTE VOLUMES //" << endl;
  cout << "////////////////////" << endl;
  string Pvol_path_test = string("test/VOL/");
  Rn::setDimension(3);
  Rn::setTolerance(0.000001);
  boost::shared_ptr<Polytope_Rn> Pvol1(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol2(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol3(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol4(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol5(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol6(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol7(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol8(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol9(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol10(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol11(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol_cube3d(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol_cube4d(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol_cube5d(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol_cube6d(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("teis1_3d.ptop"), Pvol1) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("teis2_3d.ptop"), Pvol2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("teis3_3d.ptop"), Pvol3) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("teis4_3d.ptop"), Pvol4) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("teis5_3d.ptop"), Pvol5) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("teis6_3d.ptop"), Pvol6) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("teis7_3d.ptop"), Pvol7) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("teis8_3d.ptop"), Pvol8) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("teis9_3d.ptop"), Pvol9) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("teis10_3d.ptop"), Pvol10) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("teis11_3d.ptop"), Pvol11) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("cube3D.ptop"), Pvol_cube3d) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("cube4D.ptop"), Pvol_cube4d) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("cube5D.ptop"), Pvol_cube5d) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("cube6D.ptop"), Pvol_cube6d) == TEST_OK );
  double volume1 = politopixAPI::computeVolume(Pvol1);
  double volume2 = politopixAPI::computeVolume(Pvol2);
  double volume3 = politopixAPI::computeVolume(Pvol3);
  double volume4 = politopixAPI::computeVolume(Pvol4);
  double volume5 = politopixAPI::computeVolume(Pvol5);
  double volume6 = politopixAPI::computeVolume(Pvol6);
  double volume7 = politopixAPI::computeVolume(Pvol7);
  double volume8 = politopixAPI::computeVolume(Pvol8);
  double volume9 = politopixAPI::computeVolume(Pvol9);
  double volume10 = politopixAPI::computeVolume(Pvol10);
  double volume11 = politopixAPI::computeVolume(Pvol11);
  BOOST_REQUIRE(6.65110e-07  < volume1 && volume1 < 6.65112e-07);
  BOOST_REQUIRE(6.07855e-07  < volume2 && volume2 < 6.07857e-07);
  BOOST_REQUIRE(1.81130e-06  < volume3 && volume3 < 1.81132e-06);
  BOOST_REQUIRE(7.37121e-07  < volume4 && volume4 < 7.37123e-07);
  BOOST_REQUIRE(5.73780e-07  < volume5 && volume5 < 5.73782e-07);
  BOOST_REQUIRE(3.08791e-06  < volume6 && volume6 < 3.08793e-06);
  BOOST_REQUIRE(1.11362e-06  < volume7 && volume7 < 1.11364e-06);
  BOOST_REQUIRE(1.54424e-06  < volume8 && volume8 < 1.54426e-06);
  BOOST_REQUIRE(1.24113e-05  < volume9 && volume9 < 1.24116e-05);
  BOOST_REQUIRE(3.21623e-06  < volume10 && volume10 < 3.21625e-06);
  BOOST_REQUIRE(2.55674e-05  < volume11 && volume11 < 2.55676e-05);
  double volumeC1 = politopixAPI::computeVolume(Pvol_cube3d);
  Rn::setDimension(4);
  double volumeC2 = politopixAPI::computeVolume(Pvol_cube4d);
  Rn::setDimension(5);
  double volumeC3 = politopixAPI::computeVolume(Pvol_cube5d);
  Rn::setDimension(6);
  double volumeC4 = politopixAPI::computeVolume(Pvol_cube6d);
  BOOST_REQUIRE(volumeC1 == 1000);
  BOOST_REQUIRE(volumeC2 == 10000);
  BOOST_REQUIRE(volumeC3 == 100000);
  BOOST_REQUIRE(volumeC4 == 1000000);

  boost::shared_ptr<Polytope_Rn> Pvol_A1(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol_A2(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol_A3(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pvol_A4(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("BIR3-4-6_v.ptop"), Pvol_A1) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("20D6_v.ptop"), Pvol_A2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("DG1.ptop"), Pvol_A3) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(Pvol_path_test + string("DG2.ptop"), Pvol_A4) == TEST_OK );
  Rn::setDimension(4);
  double volumeA1 = politopixAPI::computeVolume(Pvol_A1);
  Rn::setDimension(6);
  double volumeA2 = politopixAPI::computeVolume(Pvol_A2);
  double volumeA3 = politopixAPI::computeVolume(Pvol_A3);
  double volumeA4 = politopixAPI::computeVolume(Pvol_A4);
  BOOST_REQUIRE(0.124  < volumeA1 && volumeA1 < 0.126);
  BOOST_REQUIRE(0.0114730  < volumeA2 && volumeA2 < 0.0114733);
  BOOST_REQUIRE(0.015947725  < volumeA3 && volumeA3 < 0.015947727);
  BOOST_REQUIRE(38.8410  < volumeA4 && volumeA4 < 38.8412);

  cout << endl;
  cout << "/////////////////" << endl;
  cout << "// TEST CUBES //" << endl;
  cout << "///////////////" << endl;
  Rn::setDimension(6);
  Rn::setTolerance(0.000001);
  boost::shared_ptr<Polytope_Rn> cube_0, cube_1, cube_2;
  boost::shared_ptr<Polytope_Rn> res_cube_a(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> res_cube_b(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::makeCube(cube_0, 1.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::makeCube(cube_1, 10.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::makeCube(cube_2, 100.) == TEST_OK );
  cube_0->checkTopologyAndGeometry();
  cube_1->checkTopologyAndGeometry();
  BOOST_REQUIRE( politopixAPI::computeIntersection(cube_0, cube_1, res_cube_a) == TEST_OK );
  res_cube_a->checkTopologyAndGeometry();
  BOOST_REQUIRE( politopixAPI::computeIntersection(res_cube_a, cube_2, res_cube_b) == TEST_OK );

  cout << endl;
  cout << "/////////////////////////////" << endl;
  cout << "// REMOVE CAP HALF-SPACES //" << endl;
  cout << "///////////////////////////" << endl;
  Rn::setDimension(5);
  Rn::setTolerance(0.000001);
  string ppiston = string("test/CAPS/PISTON/5D/");
  /// PSEUDO SUM METHOD ///
  boost::shared_ptr<Polytope_Rn> P1_5D(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> P2_5D(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> P3_5D(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> P4_5D(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(ppiston + string("P1.ptop"), P1_5D) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(P1_5D, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(ppiston + string("P2.ptop"), P2_5D) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(P2_5D, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(ppiston + string("P3.ptop"), P3_5D) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(P3_5D, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(ppiston + string("P4.ptop"), P4_5D) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(P4_5D, 1000.) == TEST_OK );
  boost::shared_ptr<Polytope_Rn> ps_P12_5D_S1(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> ps_P12_5D_S2(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> ps_P123_5D_S1(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> ps_P123_5D_S2(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> ps_P1234_5D_S1(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> ps_P1234_5D_S2(new Polytope_Rn());
  std::set< unsigned int > P1Caps, P2Caps, P3Caps, P4Caps;
  std::set< unsigned int > Sum12Caps, Sum123Caps, Sum1234Caps;
  {for (unsigned int i=0; i<2; ++i) {
    P2Caps.insert(i);
    P3Caps.insert(i);
    P4Caps.insert(i);
  }}
  {for (unsigned int i=0; i<4; ++i) {
    P1Caps.insert(i);
  }}
  cout << "############################" << endl;
  cout << "# pseudosum: P1 + P2 = P12 #" << endl;
  cout << "############################" << endl;
  BOOST_REQUIRE( politopixAPI::pseudoSum(P1_5D, P2_5D, ps_P12_5D_S1, P1Caps, P2Caps, Sum12Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoSum(P2_5D, P1_5D, ps_P12_5D_S2, P2Caps, P1Caps, Sum12Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(ps_P12_5D_S1, ps_P12_5D_S2) == TEST_OK );
  cout << "##############################" << endl;
  cout << "# pseudosum: P12 + P3 = P123 #" << endl;
  cout << "##############################" << endl;
  BOOST_REQUIRE( politopixAPI::pseudoSum(ps_P12_5D_S1, P3_5D, ps_P123_5D_S1, Sum12Caps, P3Caps, Sum123Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoSum(P3_5D, ps_P12_5D_S1, ps_P123_5D_S2, P3Caps, Sum12Caps, Sum123Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(ps_P123_5D_S1, ps_P123_5D_S2) == TEST_OK );
  cout << "################################" << endl;
  cout << "# pseudosum: P123 + P4 = P1234 #" << endl;
  cout << "################################" << endl;
  BOOST_REQUIRE( politopixAPI::pseudoSum(ps_P123_5D_S1, P4_5D, ps_P1234_5D_S1, Sum123Caps, P4Caps, Sum1234Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoSum(P4_5D, ps_P123_5D_S1, ps_P1234_5D_S2, P4Caps, Sum123Caps, Sum1234Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(ps_P1234_5D_S1, ps_P1234_5D_S2) == TEST_OK );
  /// FULL METHOD ///
  boost::shared_ptr<Polytope_Rn> polP12_5D(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> shgP12_5D(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> polP123_5D(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> shgP123_5D(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> polP1234_5D(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> shgP1234_5D(new Polytope_Rn());
  cout << "############################" << endl;
  cout << "# with caps: P1 + P2 = P12 #" << endl;
  cout << "############################" << endl;
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(P1_5D, P2_5D, polP12_5D) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(ppiston + string("shgP12.ptop"), shgP12_5D) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfVertices(polP12_5D, shgP12_5D) == TEST_OK );
  cout << "##############################" << endl;
  cout << "# with caps: P12 + P3 = P123 #" << endl;
  cout << "##############################" << endl;
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(polP12_5D, P3_5D, polP123_5D) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(ppiston + string("shgP123.ptop"), shgP123_5D) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfVertices(polP123_5D, shgP123_5D) == TEST_OK );
  cout << "################################" << endl;
  cout << "# with caps: P123 + P4 = P1234 #" << endl;
  cout << "################################" << endl;
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(polP123_5D, P4_5D, polP1234_5D) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(ppiston + string("shgP1234.ptop"), shgP1234_5D) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfVertices(polP1234_5D, shgP1234_5D) == TEST_OK );
  ///
  std::set< unsigned int > hyperplanes2project2D;
  hyperplanes2project2D.insert(3);
  hyperplanes2project2D.insert(4);
  boost::shared_ptr<Polytope_Rn> proj_polP1234_5D(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> proj_ps_P1234_5D_S1(new Polytope_Rn());
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project2D, ps_P1234_5D_S1, proj_ps_P1234_5D_S1);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project2D, polP1234_5D, proj_polP1234_5D);
  Rn::setDimension(2);
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(proj_ps_P1234_5D_S1, proj_polP1234_5D) == TEST_OK );
  boost::shared_ptr<Polytope_Rn> finalP2D(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(ppiston + string("P1234_2D.ptop"), finalP2D) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(proj_ps_P1234_5D_S1, finalP2D) == TEST_OK );
  /////////////////////////////////////////////
  string path_test_caps = string("test/CAPS/");
  {for (int i=1; i<argc; ++i) {
    if (strcmp(argv[i], "piston") == 0) {
      Rn::setDimension(6);
      Rn::setTolerance(0.000001);
      cout << "##########" << endl;
      cout << "# Piston #" << endl;
      cout << "##########" << endl;
      string path_test_caps_piston = path_test_caps+string("PISTON/");
      boost::shared_ptr<Polytope_Rn> pist1_TestCaps(new Polytope_Rn());
      boost::shared_ptr<Polytope_Rn> pist2_TestCaps(new Polytope_Rn());
      boost::shared_ptr<Polytope_Rn> pist_Sum_TestCaps(new Polytope_Rn());
      boost::shared_ptr<Polytope_Rn> pist_Sum_TestCaps_check(new Polytope_Rn());
      BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps_piston + string("P1_6D_HV.ptop"), pist1_TestCaps) == TEST_OK );
      BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps_piston + string("P2_6D_HV.ptop"), pist2_TestCaps) == TEST_OK );
      std::set< unsigned int > PO1Caps, PO2Caps, PSum_1Caps_;
      PO1Caps.insert(0); PO1Caps.insert(1); PO1Caps.insert(2); PO1Caps.insert(3);
      PO2Caps.insert(0); PO2Caps.insert(1); PO2Caps.insert(2); PO2Caps.insert(3);
      BOOST_REQUIRE( politopixAPI::pseudoSum(pist1_TestCaps, pist2_TestCaps, pist_Sum_TestCaps, PO1Caps, PO2Caps, PSum_1Caps_) == TEST_OK );
      BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps_piston + string("Sum_6D_proj_HV.ptop"), pist_Sum_TestCaps_check) == TEST_OK );
      BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(pist_Sum_TestCaps, pist_Sum_TestCaps_check) == TEST_OK );
    }
  }}
  Rn::setDimension(3);
  Rn::setTolerance(0.0001);
  cout << "#########" << endl;
  cout << "# R3 Ex #" << endl;
  cout << "#########" << endl;
  boost::shared_ptr<Polytope_Rn> p1_TestCaps(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> p2_TestCaps(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Sum_TestCaps(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Sum_TestCaps_check(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps + string("cube3D.ptop"), p1_TestCaps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps + string("lmp1_v.ptop"), p2_TestCaps) == TEST_OK );
  std::set< unsigned int > O1Caps, O2Caps, Sum_1Caps_;
  O1Caps.insert(0); O1Caps.insert(1); O1Caps.insert(3); O1Caps.insert(4);
  BOOST_REQUIRE( politopixAPI::pseudoSum(p1_TestCaps, p2_TestCaps, Sum_TestCaps, O1Caps, O2Caps, Sum_1Caps_) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps + string("lmp1_v_cube3D_testCaps.ptop"), Sum_TestCaps_check) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Sum_TestCaps, Sum_TestCaps_check) == TEST_OK );

  cout << "###########" << endl;
  cout << "# CIRP R6 #" << endl;
  cout << "###########" << endl;
  Rn::setDimension(6);
  Rn::setTolerance(0.000001);
  boost::shared_ptr<Polytope_Rn> p1_TestCaps_b(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> p2_TestCaps_b(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Sum_TestCaps_b(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Sum_TestCaps_check_b(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps + string("P1_6D_HV.ptop"), p1_TestCaps_b) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps + string("P2_6D_HV.ptop"), p2_TestCaps_b) == TEST_OK );
  std::set< unsigned int > O1Caps_b, O2Caps_b, Sum_1Caps_b;
  O1Caps_b.insert(0); O1Caps_b.insert(1); O1Caps_b.insert(2); O1Caps_b.insert(3); O1Caps_b.insert(4); O1Caps_b.insert(5);
  O2Caps_b.insert(0); O2Caps_b.insert(1); O2Caps_b.insert(2); O2Caps_b.insert(3); O2Caps_b.insert(4); O2Caps_b.insert(5);
  BOOST_REQUIRE( politopixAPI::pseudoSum(p1_TestCaps_b, p2_TestCaps_b, Sum_TestCaps_b, O1Caps_b, O2Caps_b, Sum_1Caps_b) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps + string("P1_P2_testCaps.ptop"), Sum_TestCaps_check_b) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Sum_TestCaps_b, Sum_TestCaps_check_b) == TEST_OK );

  cout << "############" << endl;
  cout << "# SPEC CAP #" << endl;
  cout << "############" << endl;
  Rn::setDimension(6);
  Rn::setTolerance(0.000001);
  string path_test_caps2 = string("test/CAPS/spec_small_disc/");
  boost::shared_ptr<Polytope_Rn> Pc_21(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pc_22(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pc_23(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pc_24(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pc_25(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pc_26(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pg_37(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pg_38(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Pg_14(new Polytope_Rn());
  cout << "#########################" << endl;
  cout << "# SPEC CAP: truncations #" << endl;
  cout << "#########################" << endl;
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps2 + string("Polytope_10.ptop"), Pc_21) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(Pc_21, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps2 + string("Polytope_11.ptop"), Pc_22) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(Pc_22, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps2 + string("Polytope_3.ptop"),  Pc_23) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(Pc_23, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps2 + string("Polytope_4.ptop"),  Pc_24) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(Pc_24, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps2 + string("Polytope_12.ptop"), Pc_25) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(Pc_25, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps2 + string("Polytope_13.ptop"), Pc_26) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(Pc_26, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps2 + string("Polytope_7.ptop"),  Pg_37) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(Pg_37, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps2 + string("Polytope_8.ptop"),  Pg_38) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(Pg_38, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps2 + string("Polytope_9.ptop"),  Pg_14) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(Pg_14, 1000.) == TEST_OK );
  std::set< unsigned int > Int_1Caps, Int_2Caps, Sum_1Caps, Sum_1_s2Caps, P_11_20_FC1Caps, P_20_30_FC1Caps;
  std::set< unsigned int > Pc_21Caps, Pc_22Caps, Pc_23Caps, Pc_24Caps, Pc_25Caps, Pc_26Caps, Pg_37Caps;
  {for (unsigned int i=0; i<8; ++i) {
    Pc_21Caps.insert(i);Pc_22Caps.insert(i);
  }}
  {for (unsigned int i=0; i<6; ++i) {
    Pc_23Caps.insert(i);Pc_24Caps.insert(i);
  }}
  {for (unsigned int i=0; i<4; ++i) {
    Pg_37Caps.insert(i);
  }}
  Pc_25Caps.insert(0);Pc_26Caps.insert(0);Pc_25Caps.insert(1);Pc_26Caps.insert(1);Pc_25Caps.insert(6);Pc_26Caps.insert(6);Pc_25Caps.insert(7);Pc_26Caps.insert(7);
  cout << endl;
  cout << "######################" << endl;
  cout << "# Operations for FC1 #" << endl;
  cout << "######################" << endl;
  boost::shared_ptr<Polytope_Rn> Int_1, Int_2, P_11_20_FC1, P_20_30_FC1;
  boost::shared_ptr<Polytope_Rn> Sum_1(new Polytope_Rn()); boost::shared_ptr<Polytope_Rn> Sum_1_s2(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::pseudoIntersection(Pc_21, Pc_22, Int_1, Pc_21Caps, Pc_22Caps, Int_1Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoSum(Pc_23, Int_1, Sum_1, Pc_23Caps, Int_1Caps, Sum_1Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoSum(Int_1, Pc_23, Sum_1_s2, Int_1Caps, Pc_23Caps, Sum_1_s2Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Sum_1, Sum_1_s2) == TEST_OK );

  Int_2.reset(new Polytope_Rn()); P_11_20_FC1.reset(new Polytope_Rn()); P_20_30_FC1.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::pseudoIntersection(Sum_1, Int_1, P_11_20_FC1, Sum_1Caps, Int_1Caps, P_11_20_FC1Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoIntersection(Pc_24, Pc_25, Int_2, Pc_24Caps, Pc_25Caps, Int_2Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoIntersection(Int_2, Pc_26, P_20_30_FC1, Int_2Caps, Pc_26Caps, P_20_30_FC1Caps) == TEST_OK );

  cout << "#########################################" << endl;
  cout << "# Polytope from node 1.0 to 3.0 for FC1 #" << endl;
  cout << "#########################################" << endl;
  boost::shared_ptr<Polytope_Rn> P_11_30_FC1_s1(new Polytope_Rn()); boost::shared_ptr<Polytope_Rn> P_11_30_FC1_s2(new Polytope_Rn());
  std::set< unsigned int > P_11_30_FC1_s2Caps, P_11_30_FC1_s1Caps;
  BOOST_REQUIRE( politopixAPI::pseudoSum(P_20_30_FC1, P_11_20_FC1, P_11_30_FC1_s2, P_20_30_FC1Caps, P_11_20_FC1Caps, P_11_30_FC1_s2Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoSum(P_11_20_FC1, P_20_30_FC1, P_11_30_FC1_s1, P_11_20_FC1Caps, P_20_30_FC1Caps, P_11_30_FC1_s1Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(P_11_30_FC1_s1, P_11_30_FC1_s2) == TEST_OK );

  cout << "#####################" << endl;
  cout << "# Final sum for FC1 #" << endl;
  cout << "#####################" << endl;
  boost::shared_ptr<Polytope_Rn> Pfin_FC1_s1(new Polytope_Rn()); boost::shared_ptr<Polytope_Rn> Pfin_FC1_s2(new Polytope_Rn());
  std::set< unsigned int > Pfin_FC1_s1Caps, Pfin_FC1_s2Caps;
  BOOST_REQUIRE( politopixAPI::pseudoSum(P_11_30_FC1_s1, Pg_37, Pfin_FC1_s1, P_11_30_FC1_s1Caps, Pg_37Caps, Pfin_FC1_s1Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoSum(Pg_37, P_11_30_FC1_s1, Pfin_FC1_s2, Pg_37Caps, P_11_30_FC1_s1Caps, Pfin_FC1_s2Caps) == TEST_OK );
  //std::copy(Pfin_FC1_s1Caps.begin(), Pfin_FC1_s1Caps.end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pfin_FC1_s1, Pfin_FC1_s2) == TEST_OK );
  boost::shared_ptr<Polytope_Rn> final_sant(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_caps2 + string("final_sant.ptop"), final_sant) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pfin_FC1_s1, final_sant) == TEST_OK );
  cout << endl;



  cout << "#############################" << endl;
  cout << "# Simulation cutting system #" << endl;
  cout << "#############################" << endl;
  Rn::setDimension(6);
  Rn::setTolerance(0.000001);
  string path_test_cutting = string("test/AFILADO/");
  boost::shared_ptr<Polytope_Rn> P11(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> P33_43(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> P34_44(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> P41(new Polytope_Rn());
  cout << "########################" << endl;
  cout << "# CUTTING: truncations #" << endl;
  cout << "########################" << endl;
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_cutting + string("Afilado_Polytope_1.ptop"),  P11) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_cutting + string("Afilado_Polytope_6.ptop"),  P33_43) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_cutting + string("Afilado_Polytope_10.ptop"), P34_44) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_cutting + string("Afilado_Polytope_11.ptop"), P41) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(P11, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(P33_43, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(P34_44, 1000.) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(P41, 1000.) == TEST_OK );
  boost::shared_ptr<Polytope_Rn> P30_40(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> P11_41S2(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> P11_41S1(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> PfinS1(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> PfinS2(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> PfinS3(new Polytope_Rn());
  std::set< unsigned int > P11Caps, P33_43Caps, P34_44Caps, P41Caps;
  std::set< unsigned int > P30_40Caps, P11_41Caps, PfinCaps;
  {for (unsigned int i=0; i<10; ++i) {
    P34_44Caps.insert(i);
  }}
  {for (unsigned int i=0; i<4; ++i) {
    P11Caps.insert(i);
    P41Caps.insert(i);
    P33_43Caps.insert(i);
  }}
	
  cout << "##############" << endl;
  cout << "# Operations #" << endl;
  cout << "##############" << endl;
  BOOST_REQUIRE( politopixAPI::pseudoIntersection(P33_43, P34_44, P30_40, P33_43Caps, P34_44Caps, P30_40Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoSum(P11, P41, P11_41S1, P11Caps, P41Caps, P11_41Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoSum(P41, P11, P11_41S2, P41Caps, P11Caps, P11_41Caps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(P11_41S1, P11_41S2) == TEST_OK );
  cout << "#############" << endl;
  cout << "# Final sum #" << endl;
  cout << "#############" << endl;
  BOOST_REQUIRE( politopixAPI::pseudoSum(P11_41S1, P30_40, PfinS1, P11_41Caps, P30_40Caps, PfinCaps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::pseudoSum(P30_40, P11_41S1, PfinS2, P30_40Caps, P11_41Caps, PfinCaps) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(PfinS1, PfinS2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test_cutting + string("Pfin.ptop"), PfinS3) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(PfinS1, PfinS3) == TEST_OK );
  cout << endl;

  cout << "###################" << endl;
  cout << "# DATA1: Lattices #" << endl;
  cout << "###################" << endl;
  Rn::setDimension(3);
  path_test1 = string("./test/DATA1/");
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("lmp1.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  FaceEnumeration FaceEnumLmp1(_polytopeTest);
  FaceEnumeration::Compute(_polytopeTest, FaceEnumLmp1);
  //FaceEnumLmp1.printFacesWithFacets(std::cout);
  //FaceEnumLmp1.printFacesWithVertices(std::cout);
  std::vector< std::vector< ListOfFaces > > Lattice2Compare;
  FaceEnumeration::load(path_test1 + string("lmp1.latt"), Lattice2Compare);
  BOOST_REQUIRE(FaceEnumLmp1.getFacesWithVertices() == Lattice2Compare);
 
  Rn::setDimension(3);
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("cube3D.ptop"), _polytopeTest) == TEST_OK );
  FaceEnumeration FaceEnumCube(_polytopeTest);
  FaceEnumeration::Compute(_polytopeTest, FaceEnumCube);
  //FaceEnumCube.printFacesWithFacets(std::cout);
  //FaceEnumCube.printFacesWithVertices(std::cout);
  Lattice2Compare.clear();
  FaceEnumeration::load(path_test1 + string("cube3d.latt"), Lattice2Compare);
  BOOST_REQUIRE(FaceEnumCube.getFacesWithVertices() == Lattice2Compare);

  Rn::setDimension(6);
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("8D6.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  FaceEnumeration FaceEnum(_polytopeTest);
  FaceEnumeration::Compute(_polytopeTest, FaceEnum);
  //FaceEnum.printFacesWithFacets(std::cout);
  //FaceEnum.printFacesWithVertices(std::cout);
  Lattice2Compare.clear();
  FaceEnumeration::load(path_test1 + string("8D6.latt"), Lattice2Compare);
  BOOST_REQUIRE(FaceEnum.getFacesWithVertices() == Lattice2Compare);

  Rn::setDimension(6);
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("10D6.ptop"), _polytopeTest) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::computeDoubleDescription(_polytopeTest, 1000.) == TEST_OK );
  FaceEnumeration FaceEnum2(_polytopeTest);
  FaceEnumeration::Compute(_polytopeTest, FaceEnum2);
  //FaceEnum2.printFacesWithFacets(std::cout);
  //FaceEnum2.printFacesWithVertices(std::cout);
  Lattice2Compare.clear();
  FaceEnumeration::load(path_test1 + string("10D6.latt"), Lattice2Compare);
  BOOST_REQUIRE(FaceEnum2.getFacesWithVertices() == Lattice2Compare);

  Rn::setDimension(6);
  _polytopeTest.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(path_test1 + string("P9.ptop"), _polytopeTest) == TEST_OK );
  FaceEnumeration FaceEnum_9(_polytopeTest);
  FaceEnumeration::Compute(_polytopeTest, FaceEnum_9);
  //FaceEnum_9.printFacesWithFacets(std::cout);
  //FaceEnum_9.printFacesWithVerticesToSage(std::cout);
  //FaceEnum_9.save(std::cout);
  Lattice2Compare.clear();
  FaceEnumeration::load(path_test1 + string("P9.latt"), Lattice2Compare);
  BOOST_REQUIRE(FaceEnum_9.getFacesWithVertices() == Lattice2Compare);
  cout << endl;

  cout << "###############" << endl;
  cout << "# PROJECTIONS #" << endl;
  cout << "###############" << endl;
  Rn::setDimension(6);
  Rn::setTolerance(0.00001);
  boost::shared_ptr<HalfSpace_Rn> H1(new HalfSpace_Rn(6)), H2(new HalfSpace_Rn(6)), H3(new HalfSpace_Rn(6));
  {for (unsigned int i=0; i<6; ++i) {
    H1->setCoefficient(i,0);
    H2->setCoefficient(i,0);
    H3->setCoefficient(i,0);
  }}
  H1->setCoefficient(0,1);
  H2->setCoefficient(1,1);
  H3->setCoefficient(2,1);
  H1->setConstant(0);
  H2->setConstant(0);
  H3->setConstant(0);
  std::vector< boost::shared_ptr<HalfSpace_Rn> > arrayOfHS;
  arrayOfHS.push_back(H1);
  arrayOfHS.push_back(H2);
  arrayOfHS.push_back(H3);
  //{for (unsigned int i=0; i<arrayOfHS.size(); ++i) {
    //arrayOfHS[i]->dump(cout);
    //cout << endl;
  //}}
  string pcaps = string("test/PROJ/");
  boost::shared_ptr<Polytope_Rn> P2Test(new Polytope_Rn()), Pproj(new Polytope_Rn()), P2Compare(new Polytope_Rn()), transPol(new Polytope_Rn());

  P2Test.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("P11_41.ptop"), P2Test) == TEST_OK );
  std::set< unsigned int > hyperplanes2project;
  hyperplanes2project.insert(1);hyperplanes2project.insert(2);hyperplanes2project.insert(3);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, P2Test, Pproj);
  P2Compare.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("_P11_41_3D_H.ptop"), P2Compare) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );
  P2Compare.reset(new Polytope_Rn());
  transPol.reset(new Polytope_Rn());
  NormalFan_Rn NF(P2Test);
  NF.computeHyperplanesSeparationForProjection(arrayOfHS, transPol);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, transPol, P2Compare);
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );

  cout << endl;
  P2Test.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("P30_40.ptop"), P2Test) == TEST_OK );
  hyperplanes2project.clear();
  hyperplanes2project.insert(1);hyperplanes2project.insert(2);hyperplanes2project.insert(3);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, P2Test, Pproj);
  P2Compare.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("_P30_40_3D_H.ptop"), P2Compare) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );
  P2Compare.reset(new Polytope_Rn());
  transPol.reset(new Polytope_Rn());
  NormalFan_Rn NF_P30_40(P2Test);
  NF_P30_40.computeHyperplanesSeparationForProjection(arrayOfHS, transPol);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, transPol, P2Compare);
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );

  cout << endl;
  P2Test.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("Polytope_1.ptop"), P2Test) == TEST_OK );
  hyperplanes2project.clear();
  hyperplanes2project.insert(1);hyperplanes2project.insert(2);hyperplanes2project.insert(3);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, P2Test, Pproj);
  P2Compare.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("_Polytope_1_3D_H.ptop"), P2Compare) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );
  P2Compare.reset(new Polytope_Rn());
  transPol.reset(new Polytope_Rn());
  NormalFan_Rn NF_Polytope_1(P2Test);
  NF_Polytope_1.computeHyperplanesSeparationForProjection(arrayOfHS, transPol);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, transPol, P2Compare);
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );

  cout << endl;
  P2Test.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("Polytope_6.ptop"), P2Test) == TEST_OK );
  hyperplanes2project.clear();
  hyperplanes2project.insert(1);hyperplanes2project.insert(2);hyperplanes2project.insert(3);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, P2Test, Pproj);
  P2Compare.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("_Polytope_6_3D_H.ptop"), P2Compare) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );
  P2Compare.reset(new Polytope_Rn());
  transPol.reset(new Polytope_Rn());
  NormalFan_Rn NF_Polytope_6(P2Test);
  NF_Polytope_6.computeHyperplanesSeparationForProjection(arrayOfHS, transPol);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, transPol, P2Compare);
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );

  cout << endl;
  P2Test.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("Polytope_10.ptop"), P2Test) == TEST_OK );
  hyperplanes2project.clear();
  hyperplanes2project.insert(1);hyperplanes2project.insert(2);hyperplanes2project.insert(3);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, P2Test, Pproj);
  P2Compare.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("_Polytope_10_3D_H.ptop"), P2Compare) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );
  P2Compare.reset(new Polytope_Rn());
  transPol.reset(new Polytope_Rn());
  NormalFan_Rn NF_Polytope_10(P2Test);
  //NF_Polytope_10.dump(cout);
  NF_Polytope_10.computeHyperplanesSeparationForProjection(arrayOfHS, transPol);
  BOOST_REQUIRE( politopixAPI::savePolytope(string("step1.ptop"), transPol) == TEST_OK );
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, transPol, P2Compare);
  BOOST_REQUIRE( politopixAPI::savePolytope(string("step2.ptop"), P2Compare) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );

  cout << endl;
  P2Test.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("Polytope_11.ptop"), P2Test) == TEST_OK );
  hyperplanes2project.clear();
  hyperplanes2project.insert(1);hyperplanes2project.insert(2);hyperplanes2project.insert(3);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, P2Test, Pproj);
  P2Compare.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(pcaps + string("_Polytope_11_3D_H.ptop"), P2Compare) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );
  P2Compare.reset(new Polytope_Rn());
  transPol.reset(new Polytope_Rn());
  NormalFan_Rn NF_Polytope_11(P2Test);
  NF_Polytope_11.computeHyperplanesSeparationForProjection(arrayOfHS, transPol);
  TopGeomTools::projectPolytopeOnCanonicalHyperplanes(hyperplanes2project, transPol, P2Compare);
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Pproj, P2Compare) == TEST_OK );

  cout << "#####################" << endl;
  cout << "# TRACING GENERATORS #" << endl;
  cout << "#####################" << endl;
  Rn::setDimension(6);
  Rn::setTolerance(0.000001);
  string tal_path("test/TALADRO/");
  boost::shared_ptr<Polytope_Rn> Tal_CP1(new Polytope_Rn()),Tal_CP2(new Polytope_Rn()),Tal_SUM12(new Polytope_Rn());
  boost::shared_ptr<Polytope_Rn> Tal_checkSUM12(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(tal_path + string("CP1.ptop"), Tal_CP1) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(tal_path + string("CP2.ptop"), Tal_CP2) == TEST_OK );
  std::vector< std::vector<int> > genitorsOfTal_CP1(Tal_CP1->numberOfGenerators());
  std::vector< std::vector<int> > genitorsOfTal_CP2(Tal_CP2->numberOfGenerators());
  std::vector< std::vector<int> > genitorsOfTal_CP_1_2;
  {
    for (unsigned int i=0; i<genitorsOfTal_CP1.size(); ++i)
      genitorsOfTal_CP1[i].push_back(i);
    for (unsigned int j=0; j<genitorsOfTal_CP2.size(); ++j)
      genitorsOfTal_CP2[j].push_back(j);
  }
  BOOST_REQUIRE( politopixAPI::computeMinkowskiSumOfPolytopes(Tal_CP1, Tal_CP2, Tal_SUM12, genitorsOfTal_CP1, genitorsOfTal_CP2, genitorsOfTal_CP_1_2) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::loadPolytope(tal_path + string("CP_1_2.ptop"), Tal_checkSUM12) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(Tal_SUM12, Tal_checkSUM12) == TEST_OK );
  {
    unsigned int thisDim = Rn::getDimension();
    double tol_2 = Rn::getTolerance()*Rn::getTolerance();
    for (unsigned int i=0; i<genitorsOfTal_CP_1_2.size(); ++i) {
      boost::shared_ptr<Generator_Rn> VX(new Generator_Rn(thisDim));
      VX->makeSum(Tal_CP1->getGenerator(genitorsOfTal_CP_1_2[i][0]), Tal_CP2->getGenerator(genitorsOfTal_CP_1_2[i][1]));
      BOOST_REQUIRE( VX->isEqual2(Tal_SUM12->getGenerator(i), thisDim, tol_2) == true );
    }
  }

  cout << "####################" << endl;
  cout << "# VORONOI DIAGRAMS #" << endl;
  cout << "####################" << endl;
  Rn::setDimension(2);
  Rn::setTolerance(0.000001);
  boost::shared_ptr<Polytope_Rn> inputSpace(new Polytope_Rn());
  inputSpace->createBoundingBox(10);
  Point_Rn P1(2), P2(2), P3(2), P4(2), P5(2);
  P1.setCoordinate(0,1);
  P1.setCoordinate(1,1);
  P2.setCoordinate(0,8);
  P2.setCoordinate(1,4);
  P3.setCoordinate(0,-7);
  P3.setCoordinate(1,7);
  P4.setCoordinate(0,-4);
  P4.setCoordinate(1,-2);
  P5.setCoordinate(0,3);
  P5.setCoordinate(1,-3);
  std::vector<Point_Rn> listOfPoints;
  listOfPoints.push_back(P1);listOfPoints.push_back(P2);listOfPoints.push_back(P3);listOfPoints.push_back(P4);listOfPoints.push_back(P5);
  std::vector< boost::shared_ptr<Polytope_Rn> > allVoronoiCells;
  BOOST_REQUIRE( politopixAPI::computeVoronoiDiagram(inputSpace, listOfPoints, allVoronoiCells) == TEST_OK );
  string vorPath = string("./test/VORONOI/2D/");
  boost::shared_ptr<Polytope_Rn> Vor_Test(new Polytope_Rn());
  Vor_Test.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(vorPath + string("cell0.ptop"), Vor_Test) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(allVoronoiCells[0], Vor_Test) == TEST_OK );
  Vor_Test.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(vorPath + string("cell1.ptop"), Vor_Test) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(allVoronoiCells[1], Vor_Test) == TEST_OK );
  Vor_Test.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(vorPath + string("cell2.ptop"), Vor_Test) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(allVoronoiCells[2], Vor_Test) == TEST_OK );
  Vor_Test.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(vorPath + string("cell3.ptop"), Vor_Test) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(allVoronoiCells[3], Vor_Test) == TEST_OK );
  Vor_Test.reset(new Polytope_Rn());
  BOOST_REQUIRE( politopixAPI::loadPolytope(vorPath + string("cell4.ptop"), Vor_Test) == TEST_OK );
  BOOST_REQUIRE( politopixAPI::checkEqualityOfPolytopes(allVoronoiCells[4], Vor_Test) == TEST_OK );


  return TEST_OK;

}
