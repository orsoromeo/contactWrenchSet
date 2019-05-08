// politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2011-2016 : Delos Vincent
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
/// \file PolyhedralCone_Rn.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include <math.h>
#include <sstream>
#include <iostream>
#include "Rn.h"
#include "PolyhedralCone_Rn.h"


// HALF-SPACES
HalfSpace_Rn::State PolyhedralCone_Rn::checkPoint(const Point_Rn& thisPoint) const {
  HalfSpace_Rn::State thisState = HalfSpace_Rn::hs_IN;
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHSB(_listOfHalfSpaces);
  for (iteHSB.begin(); iteHSB.end()!=true; iteHSB.next()) {
    double halfSpaceNorm =
      std::inner_product(iteHSB.current()->begin(), iteHSB.current()->end(), iteHSB.current()->begin(), 0.);
    halfSpaceNorm = sqrt(halfSpaceNorm);
    double scalarProduct =
      std::inner_product(thisPoint.begin(), thisPoint.end(), iteHSB.current()->begin(), 0.);
    double distanceToHyperplane = (scalarProduct+iteHSB.current()->getConstant()) / halfSpaceNorm;
    if (distanceToHyperplane < -Rn::getTolerance())
      return HalfSpace_Rn::hs_OUT;
    else if (distanceToHyperplane < Rn::getTolerance())
      thisState = HalfSpace_Rn::hs_ON;
  }
  return thisState;
}

HalfSpace_Rn::State PolyhedralCone_Rn::checkPoint(
    const boost::shared_ptr<Generator_Rn>& point,
    const boost::shared_ptr<HalfSpace_Rn>& halfSpace,
    double halfSpaceNorm) const {
  if (!point)
    throw std::invalid_argument(std::string("Invalid point in PolyhedralCone_Rn::checkPoint()"));
  if (!halfSpace)
    throw std::invalid_argument(std::string("Invalid half space in PolyhedralCone_Rn::checkPoint()"));

  double scalarProduct =
    std::inner_product(point->begin(), point->end(), halfSpace->begin(), 0.);
  double distanceToHyperplane = (scalarProduct+halfSpace->getConstant()) / halfSpaceNorm;
  if (distanceToHyperplane > Rn::getTolerance()) {
    return HalfSpace_Rn::hs_IN;
  }
  else if (distanceToHyperplane < -Rn::getTolerance()) {
    return HalfSpace_Rn::hs_OUT;
  }
  else {
    return HalfSpace_Rn::hs_ON;
  }
}

bool PolyhedralCone_Rn::checkDuplicateGenerators(unsigned int& a, unsigned int& b) {
  bool dup=false;
  double TOL=Rn::getTolerance();
  unsigned int RnDIM=dimension();
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> >
    iteGN1(_listOfGenerators);
  {for (iteGN1.begin(); iteGN1.end()!=true; iteGN1.next()) {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> >
      iteGN2(_listOfGenerators);
    {for (iteGN2.begin(); iteGN2.end()!=true; iteGN2.next()) {
      if (iteGN1.current() != iteGN2.current()) {
        double dist=0.;
        bool goNext=false;
        {for (unsigned int j=0; j<RnDIM; j++) {
       	  double distj = fabs(iteGN1.current()->getCoordinate(j)-iteGN2.current()->getCoordinate(j));
          if (distj > TOL) {
            goNext = true;
            break;
          }
          dist = dist + distj*distj;
        }}
        if (dist<TOL && goNext==false) {
        	dup = true;
        	a = iteGN1.currentIteratorNumber();
        	b = iteGN2.currentIteratorNumber();
          std::cout << "@ same vertices (a=" << a << ", b=" << b <<")" << std::endl;
          std::cout << "### V" << a << " = (";
          {for (unsigned int ii=0; ii<RnDIM; ii++) {
            std::cout << iteGN1.current()->getCoordinate(ii);
            if (ii != RnDIM-1)
              std::cout << ", ";
          }}
          std::cout << ")" << std::endl << "{ ";
          {for (unsigned int ii=0; ii<iteGN1.current()->numberOfFacets(); ii++) {
            std::cout << "(";
            {for (unsigned int jj=0; jj<RnDIM; jj++) {
              std::cout << iteGN1.current()->getFacet(ii)->getCoefficient(jj) << " ";
              if (jj == RnDIM-1)
               std::cout << iteGN1.current()->getFacet(ii)->getSideAsText() << " " << -iteGN1.current()->getFacet(ii)->getConstant();
            }}
            std::cout << ") ";
          }}
          std::cout << "}" << std::endl;
          std::cout << "### V" << b << " = (";
          {for (unsigned int ii=0; ii<RnDIM; ii++) {
            std::cout << iteGN2.current()->getCoordinate(ii);
            if (ii != RnDIM-1)
              std::cout << ", ";
          }}
          std::cout << ")" << std::endl << "{ ";
          {for (unsigned int ii=0; ii<iteGN2.current()->numberOfFacets(); ii++) {
            std::cout << "(";
            {for (unsigned int jj=0; jj<RnDIM; jj++) {
              std::cout << iteGN2.current()->getFacet(ii)->getCoefficient(jj) << " ";
              if (jj == RnDIM-1)
                std::cout << iteGN2.current()->getFacet(ii)->getSideAsText() << " " << -iteGN2.current()->getFacet(ii)->getConstant();
            }}
            std::cout << ") ";
          }}
          std::cout << "}" << std::endl;
        }
      }
    }}
  }}
  return dup;
}

boost::shared_ptr<HalfSpace_Rn> PolyhedralCone_Rn::addHalfSpace(boost::shared_ptr<HalfSpace_Rn> hs, bool check) {
  if (check == false) {
    _listOfHalfSpaces.push_back(hs);
    return hs;
  }
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS(_listOfHalfSpaces);
  for (iteHS.begin(); iteHS.end()!=true; iteHS.next()) {
    unsigned int i=0;
    bool equal=true;
    while (i<dimension() && equal==true) {
      if (fabs(iteHS.current()->getCoefficient(i) - hs->getCoefficient(i)) > Rn::getTolerance())
        equal = false;
      i++;
    }
    if (equal == true)
      return iteHS.current();
  }
  // The half-space has not been found in list so insert it.
  _listOfHalfSpaces.push_back(hs);
  return hs;
}

const boost::shared_ptr<HalfSpace_Rn>& PolyhedralCone_Rn::getHalfSpace(unsigned int i) const throw (std::out_of_range) {
  if (i < numberOfHalfSpaces())
    return _listOfHalfSpaces[i];
  else {
    std::string errorMessage = Point_Rn::concatStrings(i, "PolyhedralCone_Rn::getHalfSpace_Rn");
    throw std::out_of_range(errorMessage);
  }
}


// CONVEX HULL

const boost::shared_ptr<Generator_Rn>& PolyhedralCone_Rn::getGenerator(unsigned int i) const throw (std::out_of_range) {
  if (i < numberOfGenerators())
    return _listOfGenerators[i];
  else {
    std::string errorMessage = Point_Rn::concatStrings(i, "PolyhedralCone_Rn::getGenerator");
    throw std::out_of_range(errorMessage);
  }
}

unsigned int PolyhedralCone_Rn::getGeneratorNumber(boost::shared_ptr<Generator_Rn> G) const throw (std::out_of_range,std::invalid_argument) {
  if (!G)
    throw std::invalid_argument(std::string("Invalid generator G in PolyhedralCone_Rn::getGeneratorNumber(G)"));
  if (_listOfGenerators.size() == 0)
    throw std::out_of_range(std::string("Non allocated array of generators in PolyhedralCone_Rn::getGeneratorNumber(G)"));

  unsigned int i=0;
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(_listOfGenerators);
  for (iteGN.begin(); iteGN.end()!=true && iteGN.current()!=G; iteGN.next()) {
    i++;
  }
  if (iteGN.end()==true) {
    std::ostringstream stream_;
    stream_ << "Generator G ";
    G->dump(stream_);
    stream_ << " is not stored in array in PolyhedralCone_Rn::getGeneratorNumber(G)";
    std::string valString = stream_.str();
    throw std::out_of_range(valString);
  }
  return i;
}

// CHECK POLYHEDRON

void PolyhedralCone_Rn::dump(std::ostream &this_ostream) const {

  if (_listOfHalfSpaces.size() == 0 && _listOfGenerators.size() == 0)
    return;
  this_ostream << std::endl << "#POLYTOPE" << std::endl;
  this_ostream << dimension() << " ";
  this_ostream << _listOfHalfSpaces.size() << " ";
  this_ostream << _listOfGenerators.size() << std::endl;
  this_ostream << std::endl;

  if (_listOfHalfSpaces.size() != 0) {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS(_listOfHalfSpaces);
    {for (iteHS.begin(); iteHS.end()!=true; iteHS.next()) {
      this_ostream << iteHS.current()->getConstant() << "\t";
      {for (unsigned int j=0; j<dimension(); j++) {
        this_ostream << iteHS.current()->getCoefficient(j) << "\t";
      }}
      this_ostream << iteHS.current()->getSideAsText() << "\t0." << std::endl;
    }}
    this_ostream << std::endl;
  }

  if (_listOfGenerators.size() != 0) {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(_listOfGenerators);
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
      {for (unsigned int j=0; j<dimension(); j++) {
        this_ostream << iteGN.current()->getCoordinate(j) << "\t";
      }}
      this_ostream << std::endl;
    }}
    this_ostream << std::endl;
  }

  if (_listOfHalfSpaces.size() != 0 && _listOfGenerators.size() != 0) {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(_listOfGenerators);
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
      {for (unsigned int j=0; j<_listOfGenerators[iteGN.currentIteratorNumber()]->numberOfFacets(); j++) {
        //this_ostream << getGeneratorNumber(_listOfGenerators[i]->getNeighbourGenerator(j)) << "\t";
        boost::shared_ptr<HalfSpace_Rn> F = iteGN.current()->getFacet(j);
        this_ostream << F << " " << std::endl;
        //this_ostream << getFacetNumber(F) << "\t";
        //this_ostream << "(";
        //{for (unsigned int jj=0; jj<Rn::getDimension(); jj++) {
        //this_ostream << F->getCoefficient(jj) << " ";
        //if (jj == Rn::getDimension()-1)
        //this_ostream << F->getSideAsText() << " " << -F->getConstant();
        //}}
        //this_ostream << ") " << std::endl;
      }}
      this_ostream << std::endl;
    }}
    this_ostream << std::endl;
  }
}

// ALGORITHM

void PolyhedralCone_Rn::computeMinkowskiSum(
      const boost::shared_ptr<PolyhedralCone_Rn>&,
      const boost::shared_ptr<PolyhedralCone_Rn>&) {
}

boost::shared_ptr<PolyhedralCone_Rn> PolyhedralCone_Rn::computeDualPolyhedralCone() const {
  boost::shared_ptr<PolyhedralCone_Rn> dualCone;
  dualCone.reset(new PolyhedralCone_Rn());

  // Temp array to store constraints.
  std::vector< boost::shared_ptr<HalfSpace_Rn> > HSvect;
  ///////////////////
  // Facets block //
  /////////////////
  for (unsigned int i=0; i<numberOfGenerators(); i++) {
    boost::shared_ptr<HalfSpace_Rn> HS;
    HS.reset(new HalfSpace_Rn(dimension()));
    HS->setConstant(0.);
    // The normal vector
    for (unsigned int coord_count=0; coord_count<dimension(); coord_count++) {
      HS->setCoefficient(coord_count, getGenerator(i)->getCoordinate(coord_count));
    }
    dualCone->addHalfSpace(HS);
    HSvect.push_back(HS);
  }

  ///////////////////////
  // Generators block //
  /////////////////////
  boost::shared_ptr<Generator_Rn> VX;
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHSA(_listOfHalfSpaces);
  for (iteHSA.begin(); iteHSA.end()!=true; iteHSA.next()) {
    VX.reset(new Generator_Rn(dimension()));
    for (unsigned int coord_count=0; coord_count<dimension(); coord_count++) {
      VX->setCoordinate(coord_count, iteHSA.current()->getCoefficient(coord_count));
    }
    dualCone->addGenerator(VX);
  }

  //////////////////////////////
  // Facets per vertex block //
  ////////////////////////////
  // {Fi1, Fi2, ... }
  //for (unsigned int i=0; i<dualCone->numberOfGenerators(); i++) {
  for (iteHSA.begin(); iteHSA.end()!=true; iteHSA.next()) {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> >
      iteGN(_listOfGenerators);
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
    //for (unsigned int vtx_count=0; vtx_count<numberOfGenerators(); vtx_count++) {
      // The first generator in the dual corresponds to the first facet in the primal.
      for (unsigned int j=0; j<iteGN.current()->numberOfFacets(); j++) {
        boost::shared_ptr<HalfSpace_Rn> Fj = iteGN.current()->getFacet(j);
        if (iteHSA.current() == Fj)
          dualCone->getGenerator(iteHSA.currentIteratorNumber())->setFacet(HSvect[iteGN.currentIteratorNumber()]);
      }
    }}
  }
  //}

  return dualCone;
}

void PolyhedralCone_Rn::createTruncatedGenerator(
  const boost::shared_ptr<Generator_Rn_SD>& currentGeneratorOut,
  const boost::shared_ptr<Generator_Rn_SD>& currentGeneratorIn,
  boost::shared_ptr<Generator_Rn_SD> newV, double ay, double az, double) const {
  // For polyhedral cones.
  //for (unsigned int k=0; k<Rn::getDimension(); k++) {
    //double yi = currentGeneratorOut->getCoordinate(k);
    //double zi =  currentGeneratorIn->getCoordinate(k);
    //newV->setCoordinate(k, az*yi/(az-ay) - ay*zi/(az-ay));
  //}
  // Final operation to compare connection facets in Minkowski sums.
  //newV->normalize();
  double coef1 = az/(az-ay);
  double coef2 =-ay/(az-ay);
  newV->makeCoefSum(currentGeneratorOut, currentGeneratorIn, coef1, coef2);
}

bool PolyhedralCone_Rn::isIncluded(const boost::shared_ptr<PolyhedralCone_Rn>& B) const {
  //std::cout << Rn::getDimension() << " " << numberOfHalfSpaces() << " " << numberOfGenerators() << std::endl;
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS_B(B->_listOfHalfSpaces);
  {for (iteHS_B.begin(); iteHS_B.end()!=true; iteHS_B.next()) {
    double halfSpaceNorm = std::inner_product(iteHS_B.current()->begin(), iteHS_B.current()->end(), iteHS_B.current()->begin(), 0.);
    halfSpaceNorm = sqrt(halfSpaceNorm);
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_A(_listOfGenerators);
    {for (iteGN_A.begin(); iteGN_A.end()!=true; iteGN_A.next()) {
      HalfSpace_Rn::State currentState = checkPoint(iteGN_A.current(), iteHS_B.current(), halfSpaceNorm);
      if (currentState == HalfSpace_Rn::hs_OUT)
        return false;
    }}
  }}
  return true;
}

void PolyhedralCone_Rn::relocateGenerators() {
  double TOL = Rn::getTolerance();
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(_listOfGenerators);
  {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
    boost::numeric::ublas::vector<double> averagePoint(dimension());
    for (unsigned i=0; i<averagePoint.size(); ++i)
      averagePoint(i) = 0.;
    bool isVeryClose = true;
    {for (unsigned int j=0; j<iteGN.current()->numberOfFacets(); j++) {
      boost::shared_ptr<HalfSpace_Rn> HS = iteGN.current()->getFacet(j);
      boost::numeric::ublas::vector<double> projectedPoint;
      double halfSpaceNorm = norm_2(HS->vect());
      //halfSpaceNorm = sqrt(halfSpaceNorm);
      double disPoint2Hyp = HS->computeDistancePointHyperplane(iteGN.current()->vect(), projectedPoint, halfSpaceNorm);
      //std::cout << "d" << iteGN.currentIteratorNumber() << " = " << disPoint2Hyp << std::endl;
      if (disPoint2Hyp > 0.25*TOL || disPoint2Hyp < -0.25*TOL)
        isVeryClose = false;
      averagePoint += projectedPoint;
    }}
    if (isVeryClose == false) {
      averagePoint /= iteGN.current()->numberOfFacets();
      iteGN.current()->setCoordinates(averagePoint);
    }
  }}
}

bool PolyhedralCone_Rn::checkEquality(const boost::shared_ptr<PolyhedralCone_Rn>& B, bool getFaceMapping) const {
   bool res1, res2;
   std::cout << "Check generators of A inside the half-spaces of B ..... ";
   res1 = checkGenerators(   _listOfGenerators, B->_listOfHalfSpaces, true);
   std::cout << "Check generators of B inside the half-spaces of A ..... ";
   res2 = checkGenerators(B->_listOfGenerators,    _listOfHalfSpaces, true);
   if (getFaceMapping && res1 && res2) {
     //# Dimension NumberOfHalfspaces NumberOfGenerators       # Dimension NumberOfHalfspaces NumberOfGenerators
     //3 6 8                                                   3 6 8
     //# HALFSPACES : a0 + a1.x1 + ... + an.xn >= 0.           # HALFSPACES : a0 + a1.x1 + ... + an.xn >= 0.
     //1 1 0 0                                                 1 -1 0 0
     //1 0 1 0                                                 1 0 0 1
     //1 0 0 1                                                 1 0 -1 0
     //1 -1 0 0                                                1 0 1 0
     //1 0 -1 0                                                1 0 0 -1
     //1 0 0 -1                                                1 1 0 0
     //
     //Mapping Va <-> Vb
     //0 <-> 5
     //1 <-> 1
     //2 <-> 4
     //3 <-> 0
     //4 <-> 7
     //5 <-> 3
     //6 <-> 6
     //7 <-> 2
     //Mapping Fa <=> Fb
     //0 <=> 5
     //1 <=> 3
     //2 <=> 1
     //3 <=> 0
     //4 <=> 2
     //5 <=> 4
     double TOL = Rn::getTolerance();
     std::vector<int> VaVb(   _listOfGenerators.size());
     std::vector<int> VbVa(B->_listOfGenerators.size());
     std::cout << "Mapping Va <-> Vb" << std::endl;
     {
       constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_A_(_listOfGenerators);
       {for (iteGN_A_.begin(); iteGN_A_.end()!=true; iteGN_A_.next()) {
         bool eq = false;
         constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_B_(B->_listOfGenerators);
         {for (iteGN_B_.begin(); iteGN_B_.end()!=true && eq==false; iteGN_B_.next()) {
           double dist = iteGN_B_.current()->distanceFrom(*iteGN_A_.current());
           if (dist < TOL) {
             eq = true;
             VaVb[iteGN_A_.currentIteratorNumber()] = iteGN_B_.currentIteratorNumber();
             VbVa[iteGN_B_.currentIteratorNumber()] = iteGN_A_.currentIteratorNumber();
             std::cout << iteGN_A_.currentIteratorNumber() << " <-> " << iteGN_B_.currentIteratorNumber() << std::endl;
           }
         }}
       }}
     }
     std::cout << "Mapping Fa <=> Fb" << std::endl;
     std::vector< std::vector<int> > FacetsOfB_WithGeneratorOfB(B->_listOfHalfSpaces.size());
     constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS_B(B->_listOfHalfSpaces);
     {for (iteHS_B.begin(); iteHS_B.end()!=true; iteHS_B.next()) {
       double halfSpaceNorm =
         std::inner_product(iteHS_B.current()->begin(), iteHS_B.current()->end(), iteHS_B.current()->begin(), 0.);
       halfSpaceNorm = sqrt(halfSpaceNorm);
       constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_B_(B->_listOfGenerators);
       {for (iteGN_B_.begin(); iteGN_B_.end()!=true; iteGN_B_.next()) {
         HalfSpace_Rn::State currentState = checkPoint(iteGN_B_.current(), iteHS_B.current(), halfSpaceNorm);
         if (currentState == HalfSpace_Rn::hs_ON)
           FacetsOfB_WithGeneratorOfB[ iteHS_B.currentIteratorNumber() ].push_back( iteGN_B_.currentIteratorNumber() );
       }}
     }}
     std::vector< std::vector<int> > FacetsOfA_WithGeneratorOfB(_listOfHalfSpaces.size());
     constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS_A(_listOfHalfSpaces);
     {for (iteHS_A.begin(); iteHS_A.end()!=true; iteHS_A.next()) {
       double halfSpaceNorm =
         std::inner_product(iteHS_A.current()->begin(), iteHS_A.current()->end(), iteHS_A.current()->begin(), 0.);
       halfSpaceNorm = sqrt(halfSpaceNorm);
       constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_B_(B->_listOfGenerators);
       {for (iteGN_B_.begin(); iteGN_B_.end()!=true; iteGN_B_.next()) {
         HalfSpace_Rn::State currentState = checkPoint(iteGN_B_.current(), iteHS_A.current(), halfSpaceNorm);
         if (currentState == HalfSpace_Rn::hs_ON)
           FacetsOfA_WithGeneratorOfB[ iteHS_A.currentIteratorNumber() ].push_back( iteGN_B_.currentIteratorNumber() );
       }}
     }}
     //Sort all arrays.
     std::vector< std::vector<int> >::iterator itV = FacetsOfB_WithGeneratorOfB.begin();
     {for (; itV!=FacetsOfB_WithGeneratorOfB.end(); ++itV) {
       std::sort(itV->begin(), itV->end());
     }}
     itV = FacetsOfA_WithGeneratorOfB.begin();
     {for (; itV!=FacetsOfA_WithGeneratorOfB.end(); ++itV) {
       std::sort(itV->begin(), itV->end());
     }}
     //std::cout << "Size of FacetsOfA_WithGeneratorOfB = " << FacetsOfA_WithGeneratorOfB.size() << std::endl;//@
     //{for (std::vector< std::vector<int> >::const_iterator iteFacetA = FacetsOfA_WithGeneratorOfB.begin();
       //    iteFacetA != FacetsOfA_WithGeneratorOfB.end(); ++iteFacetA) {
       //std::copy(iteFacetA->begin(), iteFacetA->end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
       //std::cout << std::endl;
     //}}
     //std::cout << "Size of FacetsOfB_WithGeneratorOfB = " << FacetsOfB_WithGeneratorOfB.size() << std::endl;//@
     //{for (std::vector< std::vector<int> >::const_iterator iteFacetB = FacetsOfB_WithGeneratorOfB.begin();
       //  iteFacetB != FacetsOfB_WithGeneratorOfB.end(); ++iteFacetB) {
       //std::copy(iteFacetB->begin(), iteFacetB->end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
       //std::cout << std::endl;
     //}}

     {for (std::vector< std::vector<int> >::const_iterator iteFacetA = FacetsOfA_WithGeneratorOfB.begin();
         iteFacetA != FacetsOfA_WithGeneratorOfB.end(); ++iteFacetA) {
       const std::vector<int>& currentGeneratorsOfB = *iteFacetA;
       std::vector<int> equivalentGeneratorsOfA(currentGeneratorsOfB.size());
       //std::vector<int>::const_iterator iteGN;
       {for (unsigned int i=0; i<currentGeneratorsOfB.size(); ++i) {
         equivalentGeneratorsOfA[i] = VbVa[currentGeneratorsOfB[i]];
       }}
       std::sort(equivalentGeneratorsOfA.begin(), equivalentGeneratorsOfA.end());
       //std::copy(equivalentGeneratorsOfA.begin(), equivalentGeneratorsOfA.end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
       std::vector< std::vector<int> >::const_iterator iteFacetB =
           std::find(FacetsOfB_WithGeneratorOfB.begin(), FacetsOfB_WithGeneratorOfB.end(), currentGeneratorsOfB);
       int d1 = iteFacetA - FacetsOfA_WithGeneratorOfB.begin();
       int d2 = iteFacetB - FacetsOfB_WithGeneratorOfB.begin();
       std::cout << d1 << " <=> " << d2 << std::endl;
     }}
   }
   return (res1 || res2);
 }

