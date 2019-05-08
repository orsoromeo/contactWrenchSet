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
/// \file Polytope_Rn.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "Rn.h"
#include "Polytope_Rn.h"
#include "IO_Polytope.h"
#include "DoubleDescription_Rn.h"


/// \brief Class dedicated to degeneration processing when looking for neighbors.
class RealNeighbours {

public:
  RealNeighbours():_iterator(0) {}

  /// \brief Tell whether a pseudo neighbor is a genuine one comparing set of half-spaces.
  /// \param commonFacets the set of common half-spaces pointers between <i>this</i> and <i>gen</i>
  /// \param gen the other generator candidate to be a genuine neighbor of <i>this</i>
  /// \param number the other generator number, if provided
  /// \return true if SO FAR none match has been found.
  bool addNeighbour(
    const std::vector< HalfSpace_Rn* >& commonFacets,
    const boost::shared_ptr<Generator_Rn>& gen,
    unsigned int number=0) {
    {for (unsigned int i=0; i<_HSPerPseudoNeighbours.size(); ++i) {
      unsigned int nbCommonFacet=0;
      {for (unsigned int j=0; j<_HSPerPseudoNeighbours[i].size(); ++j) {
        {for (unsigned int k=0; k<commonFacets.size(); ++k) {
          if (_HSPerPseudoNeighbours[i][j] == commonFacets[k]) {
            ++nbCommonFacet;
            if (nbCommonFacet == commonFacets.size()) {
              // No need to do anything as the current generator set of
              // half-spaces is included in the set number j.
              return false;
            }
            else if (nbCommonFacet == _HSPerPseudoNeighbours[i].size()) {
              // Substitute the old generator to be removed with the current one.
              _HSPerPseudoNeighbours[i] = commonFacets;
              _pseudoNeighboursNumber[i] = number;
              _pseudoNeighbours[i] = gen;
              return false;
            }
          }
        }}
      }}
    }}
    _pseudoNeighbours.push_back(gen);
    _pseudoNeighboursNumber.push_back(number);
    _HSPerPseudoNeighbours.push_back(commonFacets);
    return true;
  }

  /// \brief Iterator function.
  void begin() {_iterator=0;}

  /// \brief Iterator function.
  void next() {++_iterator;}

  /// \brief Iterator function.
  bool end() {return (_iterator==_pseudoNeighbours.size());}

  /// \brief Iterator function.
  unsigned int currentGeneratorNumber() {return _pseudoNeighboursNumber[_iterator];}

  /// \brief Iterator function.
  boost::shared_ptr<Generator_Rn> currentNeighbour() {return _pseudoNeighbours[_iterator];}

  /// \brief Display the content on the stream passed as an argument.
  void dump(std::ostream &ofs) {
    ofs << "Real ngb:" << std::endl;
    std::vector< std::vector< HalfSpace_Rn* > >::const_iterator ite;
    for (ite=_HSPerPseudoNeighbours.begin(); ite!=_HSPerPseudoNeighbours.end(); ++ite) {
      std::copy((*ite).begin(), (*ite).end(), std::ostream_iterator< HalfSpace_Rn* >(ofs, " ") );
      ofs << std::endl;
    }
  }

protected:
  /// A runner to iterate through the list of genuine neighbors.
  unsigned int _iterator;
  /// The generator numbers in a global list.
  std::vector< unsigned int > _pseudoNeighboursNumber;
  /// The generator smart pointer.
  std::vector< boost::shared_ptr<Generator_Rn> > _pseudoNeighbours;
  /// For each generator, store all raw pointers on their corresponding half-spaces.
  std::vector< std::vector< HalfSpace_Rn* > > _HSPerPseudoNeighbours;
};


bool Polytope_Rn::checkEdges() const {
  // Check all edges, built from neighbor generators and topological data, are on frontier and not inside the polytope.
  std::cout << "Check edges.......... ";
  unsigned int RnDIM= this->dimension();
  bool checkEdgesOK = true;
  {for (unsigned int i=0; i<numberOfGenerators(); i++) {
    boost::shared_ptr<Generator_Rn> V1 = _listOfGenerators[i];
    for (unsigned int j=0; j<numberOfGenerators(); j++) {
      boost::shared_ptr<Generator_Rn> V2 = _listOfGenerators[j];
      std::vector< boost::shared_ptr<HalfSpace_Rn> > commonFacets;
      if (V1 != V2 && checkNeighbours(V1, V2, commonFacets) == true) {
        // Build the edge with V1 & V2 and check it is inside all facets.
        std::vector<double> edge(RnDIM);
        for (unsigned int k=0; k<RnDIM; k++)
          edge[k] = V1->getCoordinate(k) - V2->getCoordinate(k);
        unsigned int countFacets = 0;
        double scalarProduct = 0.;
        for (unsigned int l=0; l<commonFacets.size(); l++) {
          for (unsigned int n=0; n<RnDIM; n++)
            scalarProduct = scalarProduct + edge[n]*commonFacets[l]->getCoefficient(n);
          if (scalarProduct < -Rn::getTolerance() || scalarProduct > Rn::getTolerance()) {
            std::cout << "\t### Edge for generators " << i << " and " << j << " not on facet (";
            {for (unsigned int jj=0; jj<RnDIM; jj++) {
               std::cout << commonFacets[l]->getCoefficient(jj) << " ";
               if (jj == RnDIM-1)
                 std::cout << commonFacets[l]->getSideAsText() << " " << commonFacets[l]->getConstant();
            }}
            std::cout << ") " << std::endl;
            checkEdgesOK = false;
          }
          else
            countFacets++;
        }
        // For polytopes only.
        if (countFacets < RnDIM-1) {
          std::cout << "\t### Edge for generators " << i << " and " << j << " has not enough facets" << std::endl;
          checkEdgesOK = false;
        }
      }
    }
  }}
  if (checkEdgesOK == true)
    std::cout << "OK" << std::endl;
  return checkEdgesOK;
}

/// When only the polytope facets are given, include it in a huge simplex that we are going to truncate.
void Polytope_Rn::createBoundingSimplex(double M) {
  double MIN =-M;
  double MAX = (2*dimension()-1)*M;
  unsigned int dim = dimension();
  unsigned int currentNumberOfFacets = (unsigned int)dimension()+1;
  unsigned int currentNumberOfGenerators = (unsigned int)dimension()+1;
  //std::cout << "Number of generators = " << currentNumberOfGenerators << std::endl;

  // Fill the data structure with blank Generators, 4 in dim 2, 8 in dim 3, ...
  {for (unsigned int vtx_count=0; vtx_count<currentNumberOfGenerators; vtx_count++) {
    boost::shared_ptr<Generator_Rn> VX;
    VX.reset(new Generator_Rn(dim));
    addGenerator(VX);
  }}
  // Create the simplex generators (MIN, MIN, ..., MAX, ..., MIN) where MAX
  // moves from 0 to (n-1) and the last vertex is (MIN, ..., MIN).
  {for (unsigned int coord_count=0; coord_count<dim; coord_count++) {
    {for (unsigned int vtx_count=0; vtx_count<currentNumberOfGenerators; vtx_count++) {
      //vtx_count==currentNumberOfGenerators-1 => this vertex is the corner i.e. (MIN, ..., MIN)
      if (vtx_count==coord_count && vtx_count<currentNumberOfGenerators-1)
      	getGenerator(vtx_count)->setCoordinate(coord_count, MAX);
      else
      	getGenerator(vtx_count)->setCoordinate(coord_count, MIN);
    }}
  }}
  // Create the simplex array of facets excepted the oblic one, for the moment.
  std::vector< boost::shared_ptr<HalfSpace_Rn> > supportFacets;
  // M + xi >= 0.
  {for (unsigned int facet_count=0; facet_count<currentNumberOfFacets-1; facet_count++) {
    boost::shared_ptr<HalfSpace_Rn> HalfSp(new HalfSpace_Rn(dim));
    {for (unsigned int coord_count=0; coord_count<dim; coord_count++) {
      if (facet_count == coord_count)
      	HalfSp->setCoefficient(coord_count, 1.);
      else
      	HalfSp->setCoefficient(coord_count, 0.);
    }}
    HalfSp->setConstant(M);
    supportFacets.push_back(HalfSp);
    // The first generator is the corner so it has
    // to contain all the perpendicular half-spaces.
    getGenerator(currentNumberOfGenerators-1)->setFacet(HalfSp);
  }}
  // Oblic half-space (Rn::getDimension()*M, -1., ..., -1.)
  boost::shared_ptr<HalfSpace_Rn> oblicHalfSp(new HalfSpace_Rn(dim));
  {for (unsigned int coord_count=0; coord_count<dim; coord_count++) {
    oblicHalfSp->setCoefficient(coord_count, -1.);
  }}
  oblicHalfSp->setConstant(dimension()*M);
  supportFacets.push_back(oblicHalfSp);
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteH(getListOfHalfSpaces());
  {for (iteH.begin(); iteH.end()!=true; iteH.next()) {
    supportFacets.push_back(iteH.current());
  }}
  _listOfHalfSpaces.clear();
  std::vector< boost::shared_ptr<HalfSpace_Rn> >::const_iterator itF;
  {for (itF=supportFacets.begin(); itF!=supportFacets.end();++itF) {
      _listOfHalfSpaces.push_back(*itF);
  }}
  // Make connections between the generators excepted the corner (already treated).
  {for (unsigned int vtx_count1=0; vtx_count1<currentNumberOfGenerators-1; vtx_count1++) {
    boost::shared_ptr<Generator_Rn> VX1 = getGenerator(vtx_count1);
    VX1->setFacet(oblicHalfSp);
    {for (unsigned int coord_count=0; coord_count<dim; coord_count++) {
      if (VX1->getCoordinate(coord_count) == MIN) {
        VX1->setFacet(_listOfHalfSpaces[coord_count]);
      }
    }}
  }}
  //std::string FileOut("simplex.ptop");
  //boost::shared_ptr<Polytope_Rn> A(this);
  //IO_Polytope::save(FileOut, A);
  //exit(0);
}

/// When only the polytope facets are given, include it in a huge cube that we are going to truncate.
void Polytope_Rn::createBoundingBox(double M) {
  double MAX = M;
  unsigned int dim = dimension();
  unsigned int currentNumberOfGenerators = (unsigned int)pow(static_cast<double>(2), static_cast<int>(dim));
  //std::cout << "Number of generators = " << currentNumberOfGenerators << std::endl;

  // Fill the data structure with blank Generators, 4 in dim 2, 8 in dim 3, ...
  for (unsigned int vtx_count=0; vtx_count<currentNumberOfGenerators; vtx_count++) {
    boost::shared_ptr<Generator_Rn> VX;
    VX.reset(new Generator_Rn(dim));
    addGenerator(VX);
  }
  // Create the cube generators.
  for (unsigned int coord_count=0; coord_count<dim; coord_count++) {
    int counter = 0;
    for (unsigned int vtx_count=0; vtx_count<currentNumberOfGenerators; vtx_count++) {
      if (counter >= pow(static_cast<double>(2), static_cast<int>(coord_count))) {
      	MAX = -1. * MAX;
      	counter = 0;
      }
      getGenerator(vtx_count)->setCoordinate(coord_count, MAX);
      counter++;
    }
  }
  // Create the cube array of facets.
  std::vector< boost::shared_ptr<HalfSpace_Rn> > supportFacets;
  // xi < M
  for (unsigned int facet_count=0; facet_count<dim; facet_count++) {
    boost::shared_ptr<HalfSpace_Rn> HalfSp(new HalfSpace_Rn(dim));
    for (unsigned int coord_count=0; coord_count<dim; coord_count++) {
      if (facet_count == coord_count)
      	HalfSp->setCoefficient(coord_count, -1.);
      else
      	HalfSp->setCoefficient(coord_count, 0.);
    }
    HalfSp->setConstant(M);
    supportFacets.push_back(HalfSp);
  }
  // xi > -M
  for (unsigned int facet_count=0; facet_count<dim; facet_count++) {
    boost::shared_ptr<HalfSpace_Rn> HalfSp(new HalfSpace_Rn(dim));
    for (unsigned int coord_count=0; coord_count<dim; coord_count++) {
      if (facet_count == coord_count)
      	HalfSp->setCoefficient(coord_count, 1.);
      else
      	HalfSp->setCoefficient(coord_count, 0.);
    }
    HalfSp->setConstant(M);
    supportFacets.push_back(HalfSp);
  }
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteH(getListOfHalfSpaces());
  {for (iteH.begin(); iteH.end()!=true; iteH.next()) {
    supportFacets.push_back(iteH.current());
  }}
  _listOfHalfSpaces.clear();
  std::vector< boost::shared_ptr<HalfSpace_Rn> >::const_iterator itF;
  {for (itF=supportFacets.begin(); itF!=supportFacets.end();++itF) {
      _listOfHalfSpaces.push_back(*itF);
  }}
  // Make connections between the generators.
  for (unsigned int vtx_count1=0; vtx_count1<currentNumberOfGenerators; vtx_count1++) {
    boost::shared_ptr<Generator_Rn> VX1 = getGenerator(vtx_count1);
    for (unsigned int coord_count=0; coord_count<dim; coord_count++) {
      if (VX1->getCoordinate(coord_count) == M)
      	VX1->setFacet(_listOfHalfSpaces[coord_count]);
      else if (VX1->getCoordinate(coord_count) == -M)
      	VX1->setFacet(_listOfHalfSpaces[coord_count+dim]);
    }
  }
}

/// This is the intersection vertex in the truncating algorithm, defined by the intersection between an edge and an hyperplane.
void Polytope_Rn::createTruncatedGenerator(
  const boost::shared_ptr<Generator_Rn_SD>& out,
  const boost::shared_ptr<Generator_Rn_SD>& in,
  boost::shared_ptr<Generator_Rn_SD> newV, double ay, double az, double b) const {
  // Formula for polytopes.
  double t = (-b-az) / (ay-az);
  //for (unsigned int k=0; k<Rn::getDimension(); k++) {
    //double yi = out->getCoordinate(k);
    //double zi =  in->getCoordinate(k);
    //newV->setCoordinate(k, zi+(yi-zi)*t);
  //}
  newV->makeDiff(out, in);
  newV->makeCoefSum(in, newV, 1., t);
}

bool Polytope_Rn::checkEqualityOfVertices(const boost::shared_ptr<Polytope_Rn>& B, bool printOnScreen) const {
  double RnDIM = Rn::getDimension();
  double TOL2 = Rn::getTolerance()*Rn::getTolerance();
  std::vector<bool> vtxOfA(numberOfGenerators());
  std::vector<bool> vtxOfB(B->numberOfGenerators());
  for (unsigned int i=0; i<vtxOfA.size(); ++i)
    vtxOfA[i] = false;
  for (unsigned int i=0; i<vtxOfB.size(); ++i)
    vtxOfB[i] = false;
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_A(getListOfGenerators());
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_B(B->getListOfGenerators());
  {for (iteGN_A.begin(); iteGN_A.end()!=true; iteGN_A.next()) {
    {for (iteGN_B.begin(); iteGN_B.end()!=true; iteGN_B.next()) {
      if (iteGN_A.current()->isEqual2(iteGN_B.current(), RnDIM, TOL2)) {
        vtxOfA[ iteGN_A.currentIteratorNumber() ] = true;
        vtxOfB[ iteGN_B.currentIteratorNumber() ] = true;
      }
    }}
  }}
  bool AinsideB = true, BinsideA = true;
  {for (unsigned int i=0; i<vtxOfA.size(); ++i) {
    if (vtxOfA[i] == false) {
      AinsideB = false;
      break;
    }
  }}
  {for (unsigned int i=0; i<vtxOfB.size(); ++i) {
    if (vtxOfB[i] == false) {
      BinsideA = false;
      break;
    }
  }}
  if (printOnScreen == true) {
    if (AinsideB == true)
      std::cout << "The set of vertices of A in included in the set of vertices of B." << std::endl;
    if (BinsideA == true)
      std::cout << "The set of vertices of B in included in the set of vertices of A." << std::endl;
  }

  return (AinsideB && BinsideA);
}


boost::shared_ptr<PolyhedralCone_Rn> Polytope_Rn::getPrimalCone(unsigned int u) const throw (std::out_of_range) {
  if (u >= numberOfGenerators())
    throw std::out_of_range("Polytope_Rn::getPrimalCone() index out of range");

  boost::shared_ptr<PolyhedralCone_Rn> primalCone;
  boost::shared_ptr<Generator_Rn> vx1 = getGenerator(u);
  //////////////////////////////////////////////////////////
  // Split the polytope into its primal polyhedral cones //
  ////////////////////////////////////////////////////////
  primalCone.reset(new PolyhedralCone_Rn());
  // Insert all half-spaces connected to the current vertex into the primal cone.
  for (unsigned int i=0; i<vx1->numberOfFacets(); i++) {
    primalCone->addHalfSpace(vx1->getFacet(i));
  }
  for (unsigned int v=0; v<numberOfGenerators(); v++) {
    if (u != v) {
      boost::shared_ptr<Generator_Rn> vx2 = getGenerator(v);
      std::vector< boost::shared_ptr<HalfSpace_Rn> > commonFacets;
      if (checkNeighbours(vx1, vx2, commonFacets) == true) {
        // Build the edge.
        boost::shared_ptr<Generator_Rn> edge12(new Generator_Rn(dimension()));
        for (unsigned int k=0; k<dimension(); k++) {
          edge12->setCoordinate(k, vx2->getCoordinate(k)-vx1->getCoordinate(k));
        }
        // Build the facet list.
        for (unsigned int j=0; j<commonFacets.size(); j++) {
          boost::shared_ptr<HalfSpace_Rn> Fj = commonFacets[j];
          edge12->setFacet(Fj);
        }
        // Insert the edge into the primal cone.
        primalCone->addGenerator(edge12);
      }
    }
  }
  return primalCone;
}

boost::shared_ptr<PolyhedralCone_Rn> Polytope_Rn::getPrimalCone(const boost::shared_ptr<Generator_Rn>& vx1) const {
  unsigned int RnDIM=dimension();
  boost::shared_ptr<PolyhedralCone_Rn> primalCone;
  //////////////////////////////////////////////////////////
  // Split the polytope into its primal polyhedral cones //
  ////////////////////////////////////////////////////////
  primalCone.reset(new PolyhedralCone_Rn());
  // Insert all half-spaces connected to the current vertex into the primal cone.
  for (unsigned int i=0; i<vx1->numberOfFacets(); i++) {
    primalCone->addHalfSpace(vx1->getFacet(i));
  }
  RealNeighbours rn;
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_A(getListOfGenerators());
  {for (iteGN_A.begin(); iteGN_A.end()!=true; iteGN_A.next()) {
    boost::shared_ptr<Generator_Rn> vx2 = iteGN_A.current();
    if (vx1 != vx2) {
      std::vector< HalfSpace_Rn* > commonPFacets;
      // Create an empty set for the signature.
      std::set< boost::shared_ptr<HalfSpace_Rn> > listOfRedundantHS;
      if (checkNeighboursWithHSnumbers(vx1, vx2, commonPFacets, listOfRedundantHS) == true) {
        rn.addNeighbour(commonPFacets, vx2);
      }
    }
  }}
  //rn.dump(std::cout);
  for (rn.begin(); rn.end()!=true; rn.next()) {
    boost::shared_ptr<Generator_Rn> vx2 = rn.currentNeighbour();
    // Build the edge.
    boost::shared_ptr<Generator_Rn> edge12(new Generator_Rn(RnDIM));
    edge12->makeDiff(vx2, vx1);
    // Build the facet list.
    std::vector< boost::shared_ptr<HalfSpace_Rn> > commonFacets;
    // Get facets.
    checkNeighbours(vx1, vx2, commonFacets);
    for (unsigned int j=0; j<commonFacets.size(); ++j) {
      boost::shared_ptr<HalfSpace_Rn> Fj = commonFacets[j];
      edge12->setFacet(Fj);
    }
    // Insert the edge into the primal cone.
    primalCone->addGenerator(edge12);
  }
  return primalCone;
}
