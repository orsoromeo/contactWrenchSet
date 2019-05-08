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
/// \file PolyhedralCone_Rn.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef POLYHEDRALCONE_Rn
#define POLYHEDRALCONE_Rn

#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <exception>
#include <iterator>
#include <numeric>
#include <vector>
#include <set>
#include "polito_Export.h"
#include "GeometricObjectIterator_Rn.h"
#include "Neighbours_Rn.h"
#include "Generator_Rn.h"
#include "HalfSpace_Rn.h"
#include "Rn.h"

using namespace boost::numeric::ublas;


/// \class PolyhedralCone_Rn
/// \brief Model a polyhedral cone using its two equivalent definitions : the convex hull and the half-space intersection.
/// We store its edges in _listOfHS and the positive combination of these vectors generates the polyhedral cone.
class polito_EXPORT PolyhedralCone_Rn {
  friend class lexIteratorOfListOfHalfSpaces;
  friend class constIteratorOfListOfHalfSpaces;

 public:
  /// Constructor
  PolyhedralCone_Rn() {}

  /// Constructor
  PolyhedralCone_Rn(const PolyhedralCone_Rn& A) {
    _listOfHalfSpaces = A._listOfHalfSpaces;
    _listOfGenerators = A._listOfGenerators;
  }

  /// Destructor
  virtual ~PolyhedralCone_Rn() {}

  /// Return the space dimension
  virtual unsigned int dimension() const {
    if (numberOfHalfSpaces()==0 && numberOfGenerators()==0)
      return Rn::getDimension();
    else if (numberOfGenerators()==0) {
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS(getListOfHalfSpaces());
      const boost::shared_ptr<HalfSpace_Rn>& currentHalfSpace = iteHS.current();
      return currentHalfSpace->dimension();
    }
    else {
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(getListOfGenerators());
      const boost::shared_ptr<Generator_Rn>& currentGenerator = iteGN.current();
      return currentGenerator->dimension();
    }
  }

  /// Tell whether this polyhedron is bounded or not, polyhedral cones are not.
  virtual bool isBounded() const {return false;}

  /// Two edges are neighbours in a polyhedral cone <=> they share at least (n-2) facets.
  virtual unsigned int neigbourhoodCondition() const {return (dimension()-2);}

  /// Each facet in a polyhedral cone has got (n-1) edges.
  virtual unsigned int numberOfGeneratorsPerFacet() const {return (dimension()-1);}

  /// Remove all half-spaces and generators.
  void reset() {_listOfHalfSpaces.clear(); _listOfGenerators.clear();}

  // HALF-SPACES

  /// Get the total number of half-spaces.
  unsigned int numberOfHalfSpaces() const {return _listOfHalfSpaces.size();}

  /// Add the current half-space in its list.
  boost::shared_ptr<HalfSpace_Rn> addHalfSpace(boost::shared_ptr<HalfSpace_Rn> hs, bool check=false);

  /// Remove the half-space number j from its list.
  void removeHalfSpace(unsigned int j) {_listOfHalfSpaces.removeGeometricObject(j);}

  /// Remove the half-space number j from its list.
  void removeHalfSpaces(const std::set< boost::shared_ptr<HalfSpace_Rn> >& setOfRedundantHS) {
    _listOfHalfSpaces.removeGeometricObjects( setOfRedundantHS );
  }

  /// For a given half-space, return its list index.
  unsigned int getHalfSpaceNumber(const boost::shared_ptr<HalfSpace_Rn>& F) const throw (std::out_of_range) {
    unsigned int index;
    try {
      index = _listOfHalfSpaces.find(F);
    }
    catch(std::out_of_range& excep) {
      std::string ex(excep.what());
      ex += std::string(" : Catch in PolyhedralCone_Rn::getHalfSpaceNumber()");
      throw std::out_of_range(ex);
    }
    return index;
  }

  /// Return the i-th generator.
  const boost::shared_ptr<HalfSpace_Rn>& getHalfSpace(unsigned int i) const throw (std::out_of_range);

  /// Return the list of half-spaces.
  listOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >& getListOfHalfSpaces() {return _listOfHalfSpaces;}

  /// Return the list of half-spaces.
  const listOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >& getListOfHalfSpaces() const {return _listOfHalfSpaces;}


  // CONVEX HULL

  /// Give the total number of generators.
  unsigned int numberOfGenerators() const {return _listOfGenerators.size();}

  /// Add the given generator.
  void addGenerator(boost::shared_ptr<Generator_Rn> vx) {_listOfGenerators.push_back(vx);}

  /// Return the i-th generator.
  const boost::shared_ptr<Generator_Rn>& getGenerator(unsigned int i) const throw (std::out_of_range);

  /// For a given generator, return its list index.
  unsigned int getGeneratorNumber(boost::shared_ptr<Generator_Rn> G) const throw (std::out_of_range,std::invalid_argument);

  /// Return the list of generators.
  const listOfGeometricObjects< boost::shared_ptr<Generator_Rn> >& getListOfGenerators() const {return _listOfGenerators;}

  unsigned int getListOfGeneratorsSD(std::vector< boost::shared_ptr<Generator_Rn_SD> >& currentListOfGeneratorsSD) {
    unsigned int generatorNumber=0;
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(getListOfGenerators());
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
      boost::shared_ptr<Generator_Rn_SD>
        currentGeneratorSD(new Generator_Rn_SD(*iteGN.current(), generatorNumber, Generator_Rn_SD::UNCHANGED));
      {for (unsigned int ii=0; ii<iteGN.current()->numberOfFacets(); ii++) {
        unsigned int fctNumber = getHalfSpaceNumber( iteGN.current()->getFacet(ii) );
        currentGeneratorSD->setFacet(fctNumber);
      }}
      currentGeneratorSD->orderFacets();
      currentListOfGeneratorsSD.push_back(currentGeneratorSD);
      ++generatorNumber;
    }}
    return generatorNumber;
  }

  /// Set a new list of generators. The list of half-spaces should have been previously set.
  void setListOfGeneratorsSD(const std::vector< boost::shared_ptr<Generator_Rn_SD> >& gnList) {
    listOfGeometricObjects< boost::shared_ptr<Generator_Rn> > listOfGenerators;
    std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN;
    {for (iteGN=gnList.begin(); iteGN!=gnList.end(); ++iteGN) {
      if ((*iteGN)->getStatus() ==  Generator_Rn_SD::CREATED ||
          (*iteGN)->getStatus() ==  Generator_Rn_SD::CREATED_AND_MODIFIED) {
        boost::shared_ptr<Generator_Rn> currentGenerator = (*iteGN)->makeGenerator_Rn();
        {for (unsigned int ii=0; ii<(*iteGN)->numberOfFacets(); ii++) {
          currentGenerator->setFacet( _listOfHalfSpaces[(*iteGN)->getFacet(ii)] );
        }}
        listOfGenerators.push_back( currentGenerator );
      }
      else if ((*iteGN)->getStatus() ==  Generator_Rn_SD::MODIFIED) {
        // Substitute all old half-spaces with the new ones.
        unsigned int nbGN = (*iteGN)->getGeneratorNumber();
        boost::shared_ptr<Generator_Rn> currentGenerator = _listOfGenerators[nbGN];
        currentGenerator->clearFacets();
        {for (unsigned int ii=0; ii<(*iteGN)->numberOfFacets(); ii++) {
          currentGenerator->setFacet( _listOfHalfSpaces[(*iteGN)->getFacet(ii)] );
        }}
        listOfGenerators.push_back(_listOfGenerators[nbGN]);
      }
      else if ((*iteGN)->getStatus() ==  Generator_Rn_SD::UNCHANGED) {
        // Substitute all old half-spaces with the new ones.
        unsigned int nbGN = (*iteGN)->getGeneratorNumber();
        listOfGenerators.push_back(_listOfGenerators[nbGN]);
      }
      // If it was deleted we obviously don't need it.
    }}
    _listOfGenerators.clear();
    _listOfGenerators = listOfGenerators;
  }

  void relocateGenerators();

  /// Set a new list of generators.
  void setListOfGenerators(const listOfGeometricObjects< boost::shared_ptr<Generator_Rn> >& gnList) {_listOfGenerators = gnList;}


  // CHECK POLYHEDRON

  /// Check for polytopes that vertices share <i>(n-1)</i> facets. For polyhedral cones, it must check that vectors share <i>(n-2)</i> facets.
  /// \param V1 the first  vertex to check
  /// \param V2 the second vertex to check
  /// \param commonFacets the set of common facets between V1 and V2
  /// \param listOfRedundantHS the set of half-spaces that will not be taken into account to compute commonFacets
  /// \return false if they are not neighbours, true otherwise with the list of common facets.
  bool checkNeighbours(const boost::shared_ptr<Generator_Rn>& V1,
                       const boost::shared_ptr<Generator_Rn>& V2,
                       std::vector< boost::shared_ptr<HalfSpace_Rn> >& commonFacets,
                       const std::set< boost::shared_ptr<HalfSpace_Rn> >& listOfRedundantHS) const throw (std::invalid_argument) {
    unsigned int sameHyperplane=0;
    for (unsigned int i=0; i<V1->numberOfFacets(); i++) {
      for (unsigned int j=0; j<V2->numberOfFacets(); j++) {
        if (V1->getRawFacet(i)==V2->getRawFacet(j) && listOfRedundantHS.find(V1->getFacet(i))==listOfRedundantHS.end()) {
          // Check the facet is not declared redundant !
          sameHyperplane++;
          commonFacets.push_back(V1->getFacet(i));
        }
      }
    }
    // For polyhedral cones it is Rn::getDimension()-1, for polytopes it is Rn::getDimension()-2.
    if (sameHyperplane >= neigbourhoodCondition())
      return true;
    return false;
  }

  /// \brief Test whether the current polytope V-description is inside the polytope B H-description.
  /// \param B The H-polytope
  /// \return true if \f[ this \subset B \f], false otherwise.
  bool isIncluded(const boost::shared_ptr<PolyhedralCone_Rn>& B) const;

  /// Check for polytopes that vertices share <i>(n-1)</i> facets. For polyhedral cones, it must check that vectors share <i>(n-2)</i> facets.
  /// \param V1 the first  vertex to check
  /// \param V2 the second vertex to check
  /// \param commonFacets the set of common facets between V1 and V2
  /// \param topologicalCode model the level of neighborhood: 1 for an edge, ..., (n-1) for a facet in a n-dimensional space
  /// \return false if they are not neighbours, true otherwise with the list of common facets.
  bool checkNeighbours(const boost::shared_ptr<Generator_Rn>& V1,
                       const boost::shared_ptr<Generator_Rn>& V2,
                       std::vector< boost::shared_ptr<HalfSpace_Rn> >& commonFacets,
                       unsigned int topologicalCode=1) const throw (std::invalid_argument) {
    unsigned int sameHyperplane=0;
    for (unsigned int i=0; i<V1->numberOfFacets(); i++) {
      for (unsigned int j=0; j<V2->numberOfFacets(); j++) {
        if (V1->getRawFacet(i) == V2->getRawFacet(j)) {
          // Check the facet is not declared redundant !
          sameHyperplane++;
          commonFacets.push_back(V1->getFacet(i));
        }
      }
    }
    // For polyhedral cones neigbourhoodCondition() is Rn::getDimension()-1, for polytopes it is Rn::getDimension()-2.
    if (sameHyperplane >= neigbourhoodCondition()-topologicalCode+1) {
      return true;
    }
    return false;
  }

  /// \brief Check for polytopes that vertices share <i>(n-1)</i> facets. For polyhedral cones, it must check that vectors share <i>(n-2)</i> facets.
  /// \param V1 the first  vertex to check
  /// \param V2 the second vertex to check
  /// \param commonFacets the list of common facets between V1 and V2
  /// \param topologicalCode model the level of neighborhood: 1 for an edge, ..., (n-1) for a facet in a n-dimensional space
  /// \return false if they are not neighbours, true otherwise with the list of common facets.
  bool checkNeighbours(
      const boost::shared_ptr<Generator_Rn>& V1,
      const boost::shared_ptr<Generator_Rn>& V2,
      std::vector< HalfSpace_Rn* >& commonFacets,
      unsigned int topologicalCode=1) const throw (std::invalid_argument) {
    unsigned int sameHyperplane=0;
    for (unsigned int i=0; i<V1->numberOfFacets(); ++i) {
      for (unsigned int j=0; j<V2->numberOfFacets(); ++j) {
        if (V1->getRawFacet(i)==V2->getRawFacet(j)) {
          // Check the facet is not declared redundant !
          sameHyperplane++;
          commonFacets.push_back(V1->getRawFacet(i));
       }
      }
    }
    // For polyhedral cones neigbourhoodCondition() is Rn::getDimension()-1, for polytopes it is Rn::getDimension()-2.
    if (sameHyperplane >= neigbourhoodCondition()-topologicalCode+1) {
      return true;
    }
    return false;
  }

  /// \brief Check for polytopes that vertices share <i>(n-1)</i> facets. For polyhedral cones, it must check that vectors share <i>(n-2)</i> facets.
  /// \param genIn the first  vertex to check
  /// \param genOut the second vertex to check
  /// \param commonFacets the list of common facets indices between genIn and genOut
  bool checkNeighbours(
    const boost::shared_ptr<Generator_Rn_SD>& genIn,
    const boost::shared_ptr<Generator_Rn_SD>& genOut,
    std::vector< unsigned int >& commonFacets) {
    // Compute the intersection on sorted lists of facets.
    std::set_intersection(genIn->facetsBegin(), genIn->facetsEnd(), genOut->facetsBegin(), genOut->facetsEnd(),
      std::inserter(commonFacets, commonFacets.end()));
    // For polyhedral cones neigbourhoodCondition() is Rn::getDimension()-1, for polytopes it is Rn::getDimension()-2.
    if (commonFacets.size() >= neigbourhoodCondition())
      return true;

    return false;
  }


  /// \brief Check for polytopes that vertices share <i>(n-1)</i> facets. For polyhedral cones, it must check that vectors share <i>(n-2)</i> facets.
  /// \param V1 the first  vertex to check
  /// \param V2 the second vertex to check
  /// \param commonFacets the list of common facets between V1 and V2
  /// \param listOfRedundantHS the list of facets to ignore when we process V1 and V2
  /// \param topologicalCode model the level of neighborhood: 1 for an edge, ..., (n-1) for a facet in a n-dimensional space
  /// \return false if they are not neighbours, true otherwise with the list of common facets.
  bool checkNeighboursWithHSnumbers(
      const boost::shared_ptr<Generator_Rn>& V1,
      const boost::shared_ptr<Generator_Rn>& V2,
      std::vector< HalfSpace_Rn* >& commonFacets,
      const std::set< boost::shared_ptr<HalfSpace_Rn> >& listOfRedundantHS,
      unsigned int topologicalCode=1) const throw (std::invalid_argument) {
    unsigned int sameHyperplane=0;
    for (unsigned int i=0; i<V1->numberOfFacets(); ++i) {
      for (unsigned int j=0; j<V2->numberOfFacets(); ++j) {
        if (V1->getRawFacet(i)==V2->getRawFacet(j) &&
            listOfRedundantHS.find(V1->getFacet(i))==listOfRedundantHS.end()) {
          // Check the facet is not declared redundant !
          sameHyperplane++;
          commonFacets.push_back(V1->getRawFacet(i));
       }
      }
    }
    // For polyhedral cones neigbourhoodCondition() is Rn::getDimension()-1, for polytopes it is Rn::getDimension()-2.
    if (sameHyperplane >= neigbourhoodCondition()-topologicalCode+1) {
      return true;
    }
    return false;
  }

  /// Compute and store all neighborhood relations between generators.
  /// \param neighboursA a sparse matrix where each line models a generator
  /// \param topologicalCode model the level of neighborhood: 1 for an edge, ..., (n-1) for a facet in a n-dimensional space
  void fillNeighbourMatrix(std::vector< std::vector<unsigned int> >& neighboursA, unsigned int topologicalCode=1) const throw (std::out_of_range) {
    if (neighboursA.size() != numberOfGenerators()) {
      std::string errorMessage("PolyhedralCone_Rn::fillNeighbourMatrix wrong matrix size");
      throw std::out_of_range(errorMessage);
    }
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_A1(getListOfGenerators());
    //Neighbours_Rn newGeneratorsTool;
    {for (iteGN_A1.begin(); iteGN_A1.end()!=true; iteGN_A1.next()) {
      const boost::shared_ptr<Generator_Rn>& v1_A = iteGN_A1.current();
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_A2(getListOfGenerators());
      {for (iteGN_A2.begin(); iteGN_A2.end()!=true; iteGN_A2.next()) {
        const boost::shared_ptr<Generator_Rn>& v2_A = iteGN_A2.current();
        std::vector< HalfSpace_Rn* > commonFacets;
        if (v1_A != v2_A) {
          //neighboursA(iteGN_A1.currentIteratorNumber(), iteGN_A2.currentIteratorNumber()) =
          if (checkNeighbours(v1_A, v2_A, commonFacets, topologicalCode) == true) {
            neighboursA[iteGN_A1.currentIteratorNumber()].push_back(iteGN_A2.currentIteratorNumber());
            //newGeneratorsTool.addGenerator(
            //    commonFacets,
            //    v1_A,
            //    v2_A,
            //    iteGN_A1.currentIteratorNumber(),
            //    iteGN_A2.currentIteratorNumber());
          }
        }
        else {
          //std::cout << "v1=" << iteGN_A1.currentIteratorNumber() << ", v2=" << iteGN_A2.currentIteratorNumber() << std::endl;
          neighboursA[iteGN_A1.currentIteratorNumber()].push_back(iteGN_A2.currentIteratorNumber());
        }
      }}
    }}
    //@@
    //newGeneratorsTool.dump(std::cout);
    //for (newGeneratorsTool.begin(); newGeneratorsTool.end()!=true; newGeneratorsTool.next()) {
    //  neighboursA(newGeneratorsTool.currentGenInNumber(),  newGeneratorsTool.currentGenOutNumber()) = 1;
    //  neighboursA(newGeneratorsTool.currentGenOutNumber(), newGeneratorsTool.currentGenInNumber()) = 1;
    //}
  }

  /// \brief As stated by <b>Komei Fukuda</b> <i> "the complexity of the polyhedral verification problem is unknown.
  /// Is it in P or in coNP-complete?" </i>
  /// So only 3 verifications are made in R<sup>n</sup> :
  /// <ul> <li> check all the generators are inside the H-representation </li>
  /// <li> check that all generators have at least <i>(n-1)</i> facets in the case of a polyhedral cone (<i>n</i> for a polytope) </li>
  /// <li> check that all facets have at least <i>(n-1)</i> generators in the case of a polyhedral cone (<i>n</i> for a polytope). </li> </ul>
  virtual bool checkTopologyAndGeometry() const {
    std::cout << "# Dimension NumberOfHalfspaces NumberOfGenerators" << std::endl;
    std::cout << dimension() << " " << numberOfHalfSpaces() << " " << numberOfGenerators() << std::endl;

    // Check that all generators have at least (n-1) facets in the case of a polyhedral cone and are inside the H-representation.
    bool checkGeneratorsOK = checkGenerators(_listOfGenerators, _listOfHalfSpaces);

    // Check that all facets have at least (n-1) generators in the case of a polyhedral cone.
    bool checkFacetsOK = checkFacets();

    return (checkGeneratorsOK==true && checkFacetsOK==true);
  }

  /// \brief Check a point state against the whole polyhedron.
  /// \return HalfSpace_Rn::OUT, HalfSpace_Rn::IN or HalfSpace_Rn::ON.
  HalfSpace_Rn::State checkPoint(const Point_Rn& thisPoint) const;

  /// \brief Check a point state against a half-space.
  /// \return HalfSpace_Rn::OUT, HalfSpace_Rn::IN or HalfSpace_Rn::ON.
  HalfSpace_Rn::State checkPoint(
      const boost::shared_ptr<Generator_Rn>& point,
      const boost::shared_ptr<HalfSpace_Rn>& halfSpace,
      double halfSpaceNorm) const;

  /// \brief Make sure no duplicate generators are stored.
  /// \param a in case of equality store in this parameter the index of the  first equal generator
  /// \param b in case of equality store in this parameter the index of the second equal generator
  /// \return true if there are equal generators, false otherwise.
  bool checkDuplicateGenerators(unsigned int& a, unsigned int& b);

  /// \brief Compute  all distance from the current point to all half-spaces frontiers.
  void checkGenerator(unsigned int vtxNumber, std::ostream &this_ostream) const throw (std::out_of_range) {
    if (vtxNumber >= _listOfGenerators.size())
      throw std::out_of_range("PolyhedralCone_Rn::checkGenerator(unsigned int fctNumber, std::ostream &this_ostream) inadequate generator number!");
    this_ostream.precision(15);
    double TOL=Rn::getTolerance();
    this_ostream << "Check generator " << vtxNumber << std::endl;
    _listOfGenerators[vtxNumber]->dump(std::cout);
    this_ostream << std::endl;
    {for (unsigned int i=0; i< _listOfGenerators[vtxNumber]->numberOfFacets(); ++i) {
      this_ostream << getHalfSpaceNumber( _listOfGenerators[vtxNumber]->getFacet(i) ) << " ";
    }}
    this_ostream << std::endl;
    this_ostream << std::endl;
    {for (unsigned int i=0; i< _listOfGenerators[vtxNumber]->numberOfFacets(); ++i) {
      this_ostream << "H" << getHalfSpaceNumber( _listOfGenerators[vtxNumber]->getFacet(i) ) << std::endl;
      _listOfGenerators[vtxNumber]->getFacet(i)->dump(this_ostream);
      this_ostream << std::endl;
    }}
    this_ostream << std::endl;
    this_ostream << "Check generator " << vtxNumber << std::endl;
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS_B(_listOfHalfSpaces);
    {for (iteHS_B.begin(); iteHS_B.end()!=true; iteHS_B.next()) {
      double halfSpaceNorm = std::inner_product(iteHS_B.current()->begin(), iteHS_B.current()->end(), iteHS_B.current()->begin(), 0.);
      halfSpaceNorm = sqrt(halfSpaceNorm);
      double scalarProduct = std::inner_product(_listOfGenerators[vtxNumber]->begin(), _listOfGenerators[vtxNumber]->end(), iteHS_B.current()->begin(), 0.);
      double distanceToHyperplane = (scalarProduct+iteHS_B.current()->getConstant()) / halfSpaceNorm;
      this_ostream << "H" << iteHS_B.currentIteratorNumber() << "\t";
      this_ostream << ": " << distanceToHyperplane;
      if (distanceToHyperplane<TOL && distanceToHyperplane>-TOL)
        this_ostream << " (*)";
      this_ostream << std::endl;
    }}
  }

  /// \brief Check the number of facets per generator and make sure it is compliant with its current constraints.
  /// It must verify the following property :
  /// \f[ \forall G=(g_1,...,g_n) \in \mathbb{R}^n, \exists \, (n-1) \, \bar{H}_i = \left\{ x : \sum_{j=1}^n a_{ij}x_j \geq 0 \right\}  such \, as \, G \in \bar{H}_i, i \in \{1,...,n-1\} \f]
  /// \return true if everything's OK, false otherwise.
  bool checkGenerators(
      const listOfGeometricObjects< boost::shared_ptr<Generator_Rn> >& listGenA,
      const listOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >& listHSB, bool check_all=false) const {
    std::cout << "Check generators..... ";
    std::vector< std::vector<int> > listOfFacetsPerGenerator(listGenA.size());
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS_B(listHSB);
    {for (iteHS_B.begin(); iteHS_B.end()!=true; iteHS_B.next()) {
      double halfSpaceNorm =
        std::inner_product(iteHS_B.current()->begin(), iteHS_B.current()->end(), iteHS_B.current()->begin(), 0.);
      halfSpaceNorm = sqrt(halfSpaceNorm);
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_A_(listGenA);
      {for (iteGN_A_.begin(); iteGN_A_.end()!=true; iteGN_A_.next()) {
        HalfSpace_Rn::State currentState = checkPoint(iteGN_A_.current(), iteHS_B.current(), halfSpaceNorm);
        if (currentState == HalfSpace_Rn::hs_ON)
          listOfFacetsPerGenerator[ iteGN_A_.currentIteratorNumber() ].push_back( iteHS_B.currentIteratorNumber() );
        else if (currentState == HalfSpace_Rn::hs_OUT)
          std::cout << "\t### G" << iteGN_A_.currentIteratorNumber() << " is out for half-space "
            << iteHS_B.currentIteratorNumber() << std::endl;
      }}
    }}
    bool checkGeneratorsOK=true;
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_A(listGenA);
    {for (iteGN_A.begin(); iteGN_A.end()!=true; iteGN_A.next()) {
      if (listOfFacetsPerGenerator[ iteGN_A.currentIteratorNumber() ].empty()==true ||
          listOfFacetsPerGenerator[ iteGN_A.currentIteratorNumber() ].size() < numberOfGeneratorsPerFacet())
        checkGeneratorsOK = false;
      if (check_all &&
          listOfFacetsPerGenerator[ iteGN_A.currentIteratorNumber() ].size() != iteGN_A.current()->numberOfFacets())
        checkGeneratorsOK = false;
    }}
    if (checkGeneratorsOK == true)
      std::cout << "OK" << std::endl;
    else {
      std::cout << "KO" << std::endl;
      {for (iteGN_A.begin(); iteGN_A.end()!=true; iteGN_A.next()) {
        if (listOfFacetsPerGenerator[ iteGN_A.currentIteratorNumber() ].empty()==true ||
            listOfFacetsPerGenerator[ iteGN_A.currentIteratorNumber() ].size() < numberOfGeneratorsPerFacet())
          std::cout << "\t### G"  <<  iteGN_A.currentIteratorNumber() << " has only "
            << listOfFacetsPerGenerator[ iteGN_A.currentIteratorNumber() ].size() << " facets." << std::endl;
        if (listOfFacetsPerGenerator[ iteGN_A.currentIteratorNumber() ].size() != iteGN_A.current()->numberOfFacets()) {
          std::cout << "\t### G" << iteGN_A.currentIteratorNumber() << " has different facet size "
            << listOfFacetsPerGenerator[ iteGN_A.currentIteratorNumber() ].size() << " and "
            << iteGN_A.current()->numberOfFacets() << std::endl;
          //std::copy(
            //listOfFacetsPerGenerator[ iteGN_A.currentIteratorNumber() ].begin(),
            //listOfFacetsPerGenerator[ iteGN_A.currentIteratorNumber() ].end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
        }
      }}
    }
    return checkGeneratorsOK;
  }

  /// Remove the generator number j from its list.
  void removeGenerator(unsigned int j) {_listOfGenerators.removeGeometricObject(j);}

  /// Always true in the polyhedral cone case.
  virtual bool checkEdges() const {return true;}

  void checkFacet(unsigned int fctNumber, std::ostream &this_ostream) const throw (std::out_of_range) {
    if (fctNumber >= _listOfHalfSpaces.size())
      throw std::out_of_range("PolyhedralCone_Rn::checkFacet(unsigned int fctNumber, std::ostream &this_ostream) inadequate facet number!");
    this_ostream.precision(15);
    this_ostream << "Check facet " << fctNumber << std::endl;
    _listOfHalfSpaces[fctNumber]->dump(this_ostream);
    this_ostream << std::endl;
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(_listOfGenerators);
    double halfSpaceNorm = std::inner_product(
        _listOfHalfSpaces[fctNumber]->begin(), _listOfHalfSpaces[fctNumber]->end(), _listOfHalfSpaces[fctNumber]->begin(), 0.);
    halfSpaceNorm = sqrt(halfSpaceNorm);
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
      HalfSpace_Rn::State currentState = checkPoint(iteGN.current(), _listOfHalfSpaces[fctNumber], halfSpaceNorm);
      if (currentState == HalfSpace_Rn::hs_ON)
        this_ostream << iteGN.currentIteratorNumber() << " ";
    }}
    ////////////////////
    //std::vector< boost::shared_ptr<Generator_Rn> > arrayOfGen;
    //this_ostream << "Maximize" << std::endl << "obj:";
    //this_ostream << "x_1" << std::endl << std::endl << "Subject To" << std::endl;
    //{for (unsigned int j=0; j<Rn::getDimension(); ++j) {
      //this_ostream << "r" << j << ":" << std::endl;
      //{for (unsigned int i=0; i<arrayOfGen.size(); ++i) {
        //if (arrayOfGen[i]->getCoordinate(j) > 0.)
          //this_ostream << " + ";
        //else
          //this_ostream << " - ";
        //this_ostream << fabs(arrayOfGen[i]->getCoordinate(j)) << "x_" << i << " ";
      //}}
      //this_ostream << " = " << std::endl;
    //}}
    //this_ostream << "r" << Rn::getDimension() << ":" << std::endl;
    //{for (unsigned int i=0; i<arrayOfGen.size(); ++i) {
      //this_ostream << " + x_" << i;
    //}}
    //this_ostream << " = 1." << std::endl << std::endl;
    //this_ostream << "Bounds" << std::endl;
    //{for (unsigned int i=0; i<arrayOfGen.size(); ++i) {
      //this_ostream << 0. << " <= x_" << i << std::endl;
    //}}
    //this_ostream << std::endl;
    //this_ostream << std::endl << "End" << std::endl;
    this_ostream << std::endl;
    double TOL=Rn::getTolerance();
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
      double scalarProduct = std::inner_product(iteGN.current()->begin(),iteGN.current()->end(), _listOfHalfSpaces[fctNumber]->begin(), 0.);
      double distanceToHyperplane = (scalarProduct+ _listOfHalfSpaces[fctNumber]->getConstant()) / halfSpaceNorm;
      this_ostream << "V" << iteGN.currentIteratorNumber() << "\t";
      this_ostream << ": " << distanceToHyperplane;
      if (distanceToHyperplane<TOL && distanceToHyperplane>-TOL)
        this_ostream << " (*)";
      this_ostream << std::endl;
    }}
  }

  /// \brief Detect redundant half spaces and make sure each facet has the right number of generators.
  /// \f[ A = Conv(G_1, ..., G_k) = \displaystyle{ \bigcap_{i=1}^m \bar{H}_i^+ } \f]
  /// Check that the polytope does not have generators violating a constraint defined by a half-space.<br>
  /// \f[ \forall G=(g_1,...,g_n) \in \mathbb{R}^n, \nexists \, \bar{H}_i^+ = \left\{ x : \sum_{j=1}^n a_{ij}x_j \geq 0 \right\} such \, as \, G \not\in \bar{H}_i^+ \f]
  /// Check that the polytope does have at least (n-1) generators per half-space frontier.<br>
  /// \f[ \forall \bar{H}_i, \exists \, (n-1) \, G=(g_1,...,g_n) \in \mathbb{R}^n, such \, as \, G \in \bar{H}_i \f]
  /// \return true if everything's OK, false otherwise.
  bool checkFacets() const {
    std::cout << "Check facets......... ";
    bool checkFacetsOK=true;
    unsigned int j=0;
    std::vector<int> listOfGeneratorsPerFacetON(numberOfGenerators());
    std::vector<int> listOfGeneratorsPerFacetOUT(numberOfGenerators());
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS(_listOfHalfSpaces);
    {for (iteHS.begin(); iteHS.end()!=true; iteHS.next()) {
      listOfGeneratorsPerFacetON.clear();
      listOfGeneratorsPerFacetOUT.clear();
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(_listOfGenerators);
      double halfSpaceNorm = std::inner_product(iteHS.current()->begin(), iteHS.current()->end(), iteHS.current()->begin(), 0.);
      halfSpaceNorm = sqrt(halfSpaceNorm);
      {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
        HalfSpace_Rn::State currentState = checkPoint(iteGN.current(), iteHS.current(), halfSpaceNorm);
        unsigned int i=iteGN.currentIteratorNumber();
        if (currentState == HalfSpace_Rn::hs_ON)
          listOfGeneratorsPerFacetON.push_back(i);
        else if (currentState == HalfSpace_Rn::hs_OUT)
          listOfGeneratorsPerFacetOUT.push_back(i);
      }}
      if (listOfGeneratorsPerFacetON.empty() == true || listOfGeneratorsPerFacetON.size() < numberOfGeneratorsPerFacet())
        std::cout<< std::endl << "\t=== Facet " << j << " has only " << listOfGeneratorsPerFacetON.size() << " generators."<< std::endl;
      if (listOfGeneratorsPerFacetOUT.empty() != true) {
        {for (unsigned int i=0; i<listOfGeneratorsPerFacetOUT.size(); i++) {
          checkFacetsOK = false;
          std::cout<< std::endl;
          std::cout << "\t### Facet " << j << " excludes generator " << listOfGeneratorsPerFacetOUT[i] << std::endl;
        }}
      }
      j++;
    }}
    if (checkFacetsOK == true)
      std::cout << "OK" << std::endl;
    return checkFacetsOK;
  }

  /// Check whether the current polyhedral cones and B are equal or not.
  /// If the last variable is true, print the mapping between the generators and faces of both polyhedral cones.
  bool checkEquality(const boost::shared_ptr<PolyhedralCone_Rn>& B, bool getFaceMapping=false) const;


  // ALGORITHMS

  /// Compute the symmetrical polyhedral cone.
  void negate() {_listOfHalfSpaces.negate(); _listOfGenerators.negate();}

  /// \brief Create the intersection edge in the truncating algorithm.
  /// It is defined by the intersection between a 2-face and a hyperplane, i.e. a (n-1)-face.
  /// The new egde is given by this formula where H is the current half space :
  /// \f[ newG = \displaystyle{\frac{\langle a, z \rangle}{\langle a, z-y \rangle} y - \frac{\langle a, y \rangle}{\langle a, z-y \rangle} z, H = \left\{ x : \langle a, x \rangle = \sum_{j=1}^n a_{j}x_j \geq 0 \right\}} \f]
  /// \param y The generator outside the current half-space \f$ \langle a, y \rangle < 0\f$
  /// \param z The generator  inside the current half-space \f$ \langle a, z \rangle > 0\f$
  /// \param newG The new generator is the intersection between the 2-face created by y and z, and the current halfspace hyperplane
  /// \param ay Scalar product \f$ \langle a, y \rangle \f$
  /// \param az Scalar product \f$ \langle a, z \rangle \f$
  /// \param b Half-space constant, null in the polyhedral cone case.
  virtual void createTruncatedGenerator(
      const boost::shared_ptr<Generator_Rn_SD>& y,
      const boost::shared_ptr<Generator_Rn_SD>& z,
      boost::shared_ptr<Generator_Rn_SD> newG, double ay, double az, double b=0.) const;

  /// \brief Build the dual polyhedral cone of the current one whose edge are orthogonal to the primal and vice-versa.
  boost::shared_ptr<PolyhedralCone_Rn> computeDualPolyhedralCone() const;

  /// Compute the Minkowski sum of two polyhedral cones.
  virtual void computeMinkowskiSum(
      const boost::shared_ptr<PolyhedralCone_Rn>& A, 
      const boost::shared_ptr<PolyhedralCone_Rn>& B);

  /// \brief Dump the polyhedral structure on std::cout.
  void dump(std::ostream &this_ostream) const;

  /// At the moment this function is useless in the case of polyhedral cones.
  virtual void createBoundingBox(double) {}

  /// At the moment this function is useless in the case of polyhedral cones.
  virtual void createBoundingSimplex(double) {}

 protected:
  /// The list of half-spaces defining the polytope.
  listOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > _listOfHalfSpaces;
  /// The convex hull of connected points.
  listOfGeometricObjects< boost::shared_ptr<Generator_Rn> > _listOfGenerators;
};


#endif // POLYHEDRALCONE_Rn
