// politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2015-2016 : Delos Vincent
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
/// \file politopixAPI.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include <stdio.h>
#include <boost/timer.hpp>
#include "politopixAPI.h"


int politopixAPI::savePolytope(const string& pathA, boost::shared_ptr<Polytope_Rn>& A)  throw (ios_base::failure)
{
  try {
    IO_Polytope::save(pathA, A);
  }
  catch(ios_base::failure& except) {
    cerr << "In/out exception in politopixAPI::save_polytope(const string&, boost::shared_ptr<PolyhedralCone_Rn>&): " << except.what() << endl;
    throw except;
  }
  return TEST_OK;
}

int politopixAPI::loadPolytope(const string& pathA, boost::shared_ptr<Polytope_Rn>& A)  throw (ios_base::failure)
{
  try {
    IO_Polytope::load(pathA, A);
  }
  catch(ios_base::failure& except) {
    cerr << "In/out exception in politopixAPI::load_polytope(const string&, boost::shared_ptr<PolyhedralCone_Rn>&): " << except.what() << endl;
    throw except;
  }
  return TEST_OK;
}

int politopixAPI::addHalfspace(boost::shared_ptr<Polytope_Rn>& A, const boost::shared_ptr<HalfSpace_Rn>& HS) {
  A->addHalfSpace(HS);
  return TEST_OK;
}

int politopixAPI::addGenerator(boost::shared_ptr<Polytope_Rn>& A, const boost::shared_ptr<Generator_Rn>& GN) {
  A->addGenerator(GN);
  return TEST_OK;
}

int politopixAPI::computeDoubleDescriptionWithoutCheck(boost::shared_ptr<Polytope_Rn>& A, double bb_size) throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  boost::timer this_timer;
  try {
    // Do we compute the generators or the half-spaces?
    if (A->numberOfGenerators() == 0) {
      // Create the bounding box centered on the origin including the H-polytope A.
      A->createBoundingBox(bb_size);
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > lexmin_ite(A->getListOfHalfSpaces());
      StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
      // The number of facets for a cube.
      unsigned int truncationStep = 2*A->dimension();
      DoubleDescription<
        boost::shared_ptr<PolyhedralCone_Rn>,
        constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
        StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
          DD(A, lexmin_ite, NRP, truncationStep);
    }
    else if (A->numberOfHalfSpaces() == 0) {
      DoubleDescriptionFromGenerators::Compute(A, bb_size);
    }
  }
  catch(invalid_argument& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Invalid argument exception in politopixAPI::computeDoubleDescriptionWithoutCheck() " << except.what() << endl;
    return TEST_KO;
  }
  catch(out_of_range& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Out of range exception in politopixAPI::computeDoubleDescriptionWithoutCheck() " << except.what() << endl;
    return TEST_KO;
  }
  catch(ios_base::failure& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "In/out exception in politopixAPI::computeDoubleDescriptionWithoutCheck() " << except.what() << endl;
    return TEST_KO;
  }
  catch(logic_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Logic error exception in politopixAPI::computeDoubleDescriptionWithoutCheck() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::computeDoubleDescriptionWithoutCheck() !" << endl;
    return TEST_KO;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;

  return TEST_OK;
}

int politopixAPI::computeDoubleDescription(boost::shared_ptr<Polytope_Rn>& A, double bb_size) throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  int res = computeDoubleDescriptionWithoutCheck(A, bb_size);
  if (res == TEST_KO)
    return TEST_KO;
  return A->checkTopologyAndGeometry() == true ? TEST_OK : TEST_KO;
}

int politopixAPI::computeIntersection(boost::shared_ptr<Polytope_Rn>& A,
				  const boost::shared_ptr<Polytope_Rn>& B) throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHSB(B->getListOfHalfSpaces());
  unsigned int truncationStep = A->numberOfHalfSpaces();
  for (iteHSB.begin(); iteHSB.end()!=true; iteHSB.next()) {
    A->addHalfSpace(iteHSB.current());
  }
  boost::timer this_timer;
  try {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > lexmin_ite(A->getListOfHalfSpaces());
    StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
    DoubleDescription<
      boost::shared_ptr<PolyhedralCone_Rn>,
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
      StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
        DD(A, lexmin_ite, NRP, truncationStep);
  }
  catch(invalid_argument& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Invalid argument exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(out_of_range& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Out of range exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(ios_base::failure& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "In/out exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(logic_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Logic error exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::computeIntersection() !" << endl;
    return TEST_KO;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;
    
  if (A->numberOfGenerators()==0 && A->numberOfHalfSpaces()==0) {
    cout << "The result is empty." << std::endl;
    return TEST_OK;
  }

  return A->checkTopologyAndGeometry() == true ? TEST_OK : TEST_KO;
}

int politopixAPI::computeIntersection(const boost::shared_ptr<Polytope_Rn>& A, const boost::shared_ptr<Polytope_Rn>& B,
          boost::shared_ptr<Polytope_Rn>& C) throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  // Copy the generators of A.
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteHS0(A->getListOfGenerators());
  for (iteHS0.begin(); iteHS0.end()!=true; iteHS0.next()) {
    C->addGenerator(iteHS0.current());
  }
  // Copy the half-spaces of A.
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS1(A->getListOfHalfSpaces());
  for (iteHS1.begin(); iteHS1.end()!=true; iteHS1.next()) {
    C->addHalfSpace(iteHS1.current());
  }
  // Copy the half-spaces of B.
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS2(B->getListOfHalfSpaces());
  for (iteHS2.begin(); iteHS2.end()!=true; iteHS2.next()) {
    C->addHalfSpace(iteHS2.current());
  }

  boost::timer this_timer;
  try {
    // Now the truncation.
    unsigned int truncationStep = A->numberOfHalfSpaces();
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > lexmin_ite(C->getListOfHalfSpaces());
    StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
    // The number of facets for a cube.
    DoubleDescription<
      boost::shared_ptr<PolyhedralCone_Rn>,
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
      StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
      DD(C, lexmin_ite, NRP, truncationStep);
  }
  catch(invalid_argument& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Invalid argument exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(out_of_range& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Out of range exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(ios_base::failure& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "In/out exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(logic_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Logic error exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::computeIntersection() !" << endl;
    return TEST_KO;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;

  if (C->numberOfGenerators()==0 && C->numberOfHalfSpaces()==0) {
    cout << "The result is empty." << std::endl;
    return TEST_OK;
  }

  return C->checkTopologyAndGeometry() == true ? TEST_OK : TEST_KO;
}

bool politopixAPI::isIncluded(const boost::shared_ptr<Polytope_Rn>& A,
          const boost::shared_ptr<Polytope_Rn>& B) throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  boost::timer this_timer;
  bool isInside = false;
  try {
    isInside = A->isIncluded(B);
  }
  catch(invalid_argument& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Invalid argument exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(out_of_range& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Out of range exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(ios_base::failure& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "In/out exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(logic_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Logic error exception in politopixAPI::computeIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::computeIntersection() !" << endl;
    return TEST_KO;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;

  return isInside;
}

int politopixAPI::computeIntersectionWithoutCheck(boost::shared_ptr<Polytope_Rn>& A,
					      const boost::shared_ptr<Polytope_Rn>& B) throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHSB(B->getListOfHalfSpaces());
  unsigned int truncationStep = A->numberOfHalfSpaces();
  for (iteHSB.begin(); iteHSB.end()!=true; iteHSB.next()) {
    A->addHalfSpace(iteHSB.current());
  }
  boost::timer this_timer;
  try {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > lexmin_ite(A->getListOfHalfSpaces());
    StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
    DoubleDescription<
      boost::shared_ptr<PolyhedralCone_Rn>,
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
      StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
      DD(A, lexmin_ite, NRP, truncationStep);
  }
  catch(invalid_argument& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Invalid argument exception in politopixAPI::computeIntersectionWithoutCheck() " << except.what() << endl;
    return TEST_KO;
  }
  catch(out_of_range& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Out of range exception in politopixAPI::computeIntersectionWithoutCheck() " << except.what() << endl;
    return TEST_KO;
  }
  catch(ios_base::failure& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "In/out exception in politopixAPI::computeIntersectionWithoutCheck() " << except.what() << endl;
    return TEST_KO;
  }
  catch(logic_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Logic error exception in politopixAPI::computeIntersectionWithoutCheck() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::computeIntersectionWithoutCheck() !" << endl;
    return TEST_KO;
  }
  //cout << "TIME=" << this_timer.elapsed() << endl;
    
  if (A->numberOfGenerators()==0 && A->numberOfHalfSpaces()==0) {
    //cout << "The result is empty." << std::endl;
    return TEST_OK;
  }

  return TEST_OK;
}
  
int politopixAPI::computeMinkowskiSumOfPolytopes(const boost::shared_ptr<Polytope_Rn>& A,
    const boost::shared_ptr<Polytope_Rn>& B, boost::shared_ptr<Polytope_Rn>& C)
      throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  boost::timer this_timer;
  try {
    MinkowskiSum Ope(A,B,C);
  }
  catch(invalid_argument& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Invalid argument exception in politopixAPI::computeMinkowskiSumOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(out_of_range& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Out of range exception in politopixAPI::computeMinkowskiSumOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(ios_base::failure& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "In/out exception in politopixAPI::computeMinkowskiSumOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(logic_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Logic error exception in politopixAPI::computeMinkowskiSumOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::computeMinkowskiSumOfPolytopes() !" << endl;
    return TEST_KO;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;

  return C->checkTopologyAndGeometry() == true ? TEST_OK : TEST_KO;
}

int politopixAPI::computeMinkowskiSumOfPolytopes(const boost::shared_ptr<Polytope_Rn>& A,
  const boost::shared_ptr<Polytope_Rn>& B,
  boost::shared_ptr<Polytope_Rn>& C,
  const std::vector< std::vector<int> >& genitorsOfGeneratorsA,
  const std::vector< std::vector<int> >& genitorsOfGeneratorsB,
  std::vector< std::vector<int> >& traceGenerators)
      throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  boost::timer this_timer;
  try {
    MinkowskiSum Ope(A,B,C,genitorsOfGeneratorsA,genitorsOfGeneratorsB,traceGenerators);
  }
  catch(invalid_argument& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Invalid argument exception in politopixAPI::computeMinkowskiSumOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(out_of_range& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Out of range exception in politopixAPI::computeMinkowskiSumOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(ios_base::failure& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "In/out exception in politopixAPI::computeMinkowskiSumOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(logic_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Logic error exception in politopixAPI::computeMinkowskiSumOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::computeMinkowskiSumOfPolytopes() !" << endl;
    return TEST_KO;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;

  return C->checkTopologyAndGeometry() == true ? TEST_OK : TEST_KO;
}

int politopixAPI::checkEqualityOfPolytopes(const boost::shared_ptr<Polytope_Rn>& A,
    const boost::shared_ptr<Polytope_Rn>& B, bool getFaceMapping) throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  bool ret_val;
  boost::timer this_timer;
  try {
    ret_val = A->checkEquality(B, getFaceMapping);
  }
  catch(invalid_argument& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Invalid argument exception in politopixAPI::checkEqualityOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(out_of_range& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Out of range exception in politopixAPI::checkEqualityOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(ios_base::failure& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "In/out exception in politopixAPI::checkEqualityOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(logic_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Logic error exception in politopixAPI::checkEqualityOfPolytopes() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::checkEqualityOfPolytopes() " << endl;
    return TEST_KO;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;

  return ret_val == true ? TEST_OK : TEST_KO;
}

bool politopixAPI::checkEqualityOfVertices(const boost::shared_ptr<Polytope_Rn>& A,
    const boost::shared_ptr<Polytope_Rn>& B) throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  bool ret_val;
  boost::timer this_timer;
  try {
    ret_val = A->checkEqualityOfVertices(B);
  }
  catch(invalid_argument& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Invalid argument exception in politopixAPI::checkEqualityOfVertices() " << except.what() << endl;
    return TEST_KO;
  }
  catch(out_of_range& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Out of range exception in politopixAPI::checkEqualityOfVertices() " << except.what() << endl;
    return TEST_KO;
  }
  catch(ios_base::failure& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "In/out exception in politopixAPI::checkEqualityOfVertices() " << except.what() << endl;
    return TEST_KO;
  }
  catch(logic_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Logic error exception in politopixAPI::checkEqualityOfVertices() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::checkEqualityOfVertices() " << endl;
    return TEST_KO;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;

  return ret_val == true ? TEST_OK : TEST_KO;
}

double politopixAPI::computeVolume(const boost::shared_ptr<Polytope_Rn> Pol)  throw (invalid_argument, out_of_range, ios_base::failure)
{
  double vol = -1.;
  boost::timer this_timer;
  try {
    vol = VolumeOfPolytopes_Rn::compute(Pol);
  }
  catch(std::invalid_argument& e) {
    std::cout << "politopixAPI::computeVolume() : invalid argument exception " << e.what() << std::endl;
    return -1;
  }
  catch(std::out_of_range& e) {
    std::cout << "politopixAPI::computeVolume() : out of range exception " << e.what() << std::endl;
    return -1;
  }
  catch(std::ios_base::failure& e) {
    std::cout << "politopixAPI::computeVolume() : in/out exception " << e.what() << std::endl;
    return -1;
  }
  catch(...) {
    std::cout << "politopixAPI::computeVolume() : unexpected exception caught !" << std::endl;
    return -1;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;

  return vol;
}

int politopixAPI::checkTopologyAndGeometry(const boost::shared_ptr<PolyhedralCone_Rn>& A) {
  return A->checkTopologyAndGeometry() == true ? TEST_OK : TEST_KO;
}

int politopixAPI::makeCube(boost::shared_ptr<Polytope_Rn>& A, double M) {
  A.reset(new Polytope_Rn());
  A->createBoundingBox(M);
  return TEST_OK;
}

int politopixAPI::PolarPolytope(const boost::shared_ptr<Polytope_Rn>& original_pol, boost::shared_ptr<Polytope_Rn>& polar_pol) {
  return TopGeomTools::PolarPolytope(original_pol, polar_pol);
}

int politopixAPI::Translate(boost::shared_ptr<Polytope_Rn>& pol, const boost::numeric::ublas::vector<double>& v2t) {
  return TopGeomTools::Translate(pol, v2t);
}

int politopixAPI::pseudoIntersection(
    const boost::shared_ptr<Polytope_Rn>& A,
    const boost::shared_ptr<Polytope_Rn>& B,
    boost::shared_ptr<Polytope_Rn>& C,
    const std::set< unsigned int >& firstOperandCaps,
    const std::set< unsigned int >& secondOperandCaps,
    std::set< unsigned int >& newCaps,
    double bb_size) throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  boost::timer this_timer;
  try {
    PseudoIntersectionWithoutCaps PIWC(A, B, C, firstOperandCaps, secondOperandCaps, newCaps, bb_size);
  }
  catch(invalid_argument& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Invalid argument exception in politopixAPI::PseudoIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(out_of_range& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Out of range exception in politopixAPI::PseudoIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(ios_base::failure& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "In/out exception in politopixAPI::PseudoIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(logic_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Logic error exception in politopixAPI::PseudoIntersection() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::PseudoIntersection() !" << endl;
    return TEST_KO;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;

  if (C->numberOfGenerators()==0 && C->numberOfHalfSpaces()==0) {
    cout << "The result is empty." << std::endl;
    return TEST_OK;
  }

  return C->checkTopologyAndGeometry() == true ? TEST_OK : TEST_KO;
}

int politopixAPI::pseudoSum(
    const boost::shared_ptr<Polytope_Rn>& A,
    const boost::shared_ptr<Polytope_Rn>& B,
    boost::shared_ptr<Polytope_Rn>& C,
    const std::set< unsigned int >& firstOperandCaps,
    const std::set< unsigned int >& secondOperandCaps,
    std::set< unsigned int >& newCaps,
    double bb_size) throw (invalid_argument,out_of_range,ios_base::failure,logic_error)
{
  boost::timer this_timer;
  try {
    PseudoSumWithoutCaps PSWC(A, B, C, firstOperandCaps, secondOperandCaps, newCaps, bb_size);
  }
  catch(invalid_argument& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Invalid argument exception in politopixAPI::pseudoSum() " << except.what() << endl;
    return TEST_KO;
  }
  catch(out_of_range& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Out of range exception in politopixAPI::pseudoSum() " << except.what() << endl;
    return TEST_KO;
  }
  catch(ios_base::failure& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "In/out exception in politopixAPI::pseudoSum() " << except.what() << endl;
    return TEST_KO;
  }
  catch(logic_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Logic error exception in politopixAPI::pseudoSum() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::pseudoSum() !" << endl;
    return TEST_KO;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;

  if (C->numberOfGenerators()==0 && C->numberOfHalfSpaces()==0) {
    cout << "The result is empty." << std::endl;
    return TEST_OK;
  }

  return C->checkTopologyAndGeometry() == true ? TEST_OK : TEST_KO;
}

int politopixAPI::computeVoronoiDiagram(
   const boost::shared_ptr<Polytope_Rn>& inputSpace,
   const std::vector<Point_Rn>& listOfSeeds,
   std::vector< boost::shared_ptr<Polytope_Rn> >& VoronoiCells) throw (std::length_error)
{
  boost::timer this_timer;
  Voronoi_Rn VD(inputSpace, listOfSeeds);
  try {
    VD.compute();
    VoronoiCells = VD.getVoronoiCells();
  }
  catch(length_error& except) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Length error exception in politopixAPI::computeVoronoiDiagram() " << except.what() << endl;
    return TEST_KO;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << endl;
    cerr << "Unexpected exception caught in politopixAPI::computeVoronoiDiagram() !" << endl;
    return TEST_KO;
  }
  cout << "TIME=" << this_timer.elapsed() << endl;

  return VD.checkTopologyAndGeometry() == true ? TEST_OK : TEST_KO;
}

