// politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2014-2015 : Delos Vincent
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
/// \file VolumeOfPolytopes_Rn.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef VOLUME_POLYTOPES_Rn
#define VOLUME_POLYTOPES_Rn

#include <stdexcept>
#include <exception>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <set>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "polito_Export.h"
#include "Polytope_Rn.h"
#include "Rn.h"

using namespace boost::numeric::ublas;

class PolytopeToSimplexes {

public:
  PolytopeToSimplexes(const std::vector< unsigned int >& splitPol, unsigned int dim) {
    _dimension = dim;
    //_listOfProcessedVertices.push_back( splitPol[0] );
    std::vector< unsigned int >::const_iterator itEnd = splitPol.end();
    std::vector< unsigned int >::const_iterator itBeg = splitPol.begin();
    //++itBeg;
    // Insert all elements except the first one.
    _listOfVerticesToSplit.insert(_listOfVerticesToSplit.begin(), itBeg, itEnd);
  }

  PolytopeToSimplexes(const std::vector< unsigned int >& vtx2plit, unsigned int apexNb, unsigned int dim) {
    _dimension = dim;
    _listOfProcessedVertices.push_back( apexNb );
    _listOfVerticesToSplit = vtx2plit;
  }

  PolytopeToSimplexes(const std::vector< unsigned int >& vtx2plit, const std::vector< unsigned int >& prcVtx, unsigned int apexNb, unsigned int dim) {
    _dimension = dim;
    _listOfProcessedVertices = prcVtx;
    _listOfProcessedVertices.push_back( apexNb );
    _listOfVerticesToSplit = vtx2plit;
  }

  const std::vector< unsigned int >& getListOfProcessedVertices() const {
    return _listOfProcessedVertices;
  }

  void setListOfProcessedVertices(const std::vector< unsigned int >& LPV) {
    _listOfProcessedVertices = LPV;
  }

  const std::vector< unsigned int >& getListOfVerticesToSplit() const {
    return _listOfVerticesToSplit;
  }

  void setListOfVerticesToSplit(const std::vector< unsigned int >& VTS) {
    _listOfVerticesToSplit = VTS;
  }

  void addProcessedVertex(unsigned int ap) {
    _listOfProcessedVertices.push_back(ap);
  }

  bool checkSimplex() const {
    if (_dimension+1 == _listOfVerticesToSplit.size())
      return true;
    return false;
  }

  std::vector< unsigned int > buildPolytope() const {
    std::vector< unsigned int > pol;
    std::vector< unsigned int >::const_iterator it;
    for (it=_listOfProcessedVertices.begin(); it!=_listOfProcessedVertices.end(); ++it)
      pol.push_back(*it);
    for (it=_listOfVerticesToSplit.begin(); it!=_listOfVerticesToSplit.end(); ++it)
      pol.push_back(*it);
    return pol;
  }

  void switchToFullDimension() {
    _dimension += _listOfProcessedVertices.size();
    _listOfProcessedVertices.insert(_listOfProcessedVertices.end(), _listOfVerticesToSplit.begin(), _listOfVerticesToSplit.end());
    _listOfVerticesToSplit.clear();
  }

  unsigned int dimension() const {
    return _dimension;
  }

  void dump(std::ostream &this_ostream, unsigned int shift=0) const {
    for (unsigned int i=0; i<shift; ++i) this_ostream << "\t";
    this_ostream << "{" << _dimension << ":";
    std::copy(getListOfProcessedVertices().begin(), getListOfProcessedVertices().end(), std::ostream_iterator<unsigned int>(this_ostream, " "));
    this_ostream << std::endl;
    for (unsigned int i=0; i<shift; ++i) std::cout << "\t";
    std::copy(getListOfVerticesToSplit().begin(), getListOfVerticesToSplit().end(), std::ostream_iterator<unsigned int>(this_ostream, " "));
    this_ostream << "}" << std::endl;
  }


protected:
  /// The ordered list of vertices as we go down in smaller dimensions spaces.
  std::vector< unsigned int > _listOfProcessedVertices;
  /// The ordered list of vertices to be split in lower dimensions.
  std::vector< unsigned int > _listOfVerticesToSplit;
  /// The current dimension space we work in.
  unsigned int _dimension;
};

/// \brief Split a polytope into simplices to compute its volume. <br>
/// <i> Two Algorithms for Determining Volumes of Convex Polyhedra </i> (1979) by <b>Jacques Cohen</b> and <b>Timothy Hickey</b> <br>
/// Journal of the ACM (JACM) JACM Homepage archive <br>
/// Volume 26 Issue 3, july 1979 <br>
/// Pages 401-414  <br>
class polito_EXPORT VolumeOfPolytopes_Rn {

 public:

  /// Constructor.
  VolumeOfPolytopes_Rn(const boost::shared_ptr<Polytope_Rn> P);

  ~VolumeOfPolytopes_Rn() {}

  /// Build all simplices to partition the polytope.
  /// \param DIM the dimension of the space we work in
  void splitCloudOfVertices(unsigned int DIM);

  /// Return the volume of the given polytope P.
  static double compute(const boost::shared_ptr<Polytope_Rn> P) {
    VolumeOfPolytopes_Rn volComp(P);
    unsigned int currentDIM = Rn::getDimension();
    // Put all other generators in a set (in a logical way as they are numbers)
    std::vector< unsigned int > attempt2Simplex;
    std::vector< unsigned int > listOfOthers;
    for (unsigned int i=1; i<P->numberOfGenerators(); ++i)
      listOfOthers.push_back(i);
    try {
      // Do the hard job.
      volComp.splitCloudOfVertices(currentDIM);
#ifdef DEBUG
      volComp.dump(std::cout);
#endif
      return volComp.volume();
    }
    catch(std::invalid_argument& e) {
      std::cout << "VolumeOfPolytopes_Rn::compute() : invalid argument exception " << e.what() << std::endl;
      return -1;
    }
    catch(std::out_of_range& e) {
      std::cout << "VolumeOfPolytopes_Rn::compute() : out of range exception " << e.what() << std::endl;
      return -1;
    }
    catch(std::ios_base::failure& e) {
      std::cout << "VolumeOfPolytopes_Rn::compute() : in/out exception " << e.what() << std::endl;
      return -1;
    }
    catch(...) {
      std::cout << "VolumeOfPolytopes_Rn::compute() : unexpected exception caught !" << std::endl;
      return -1;
    }
  }

  void check() const throw (std::domain_error);

  /// Sum the volumes of all simplices partitionning the polytope.
  double volume();

  /// Compute the volume of a simplex making use of the following formula :
  /// \f$ \displaystyle Vol \big( conv(v_0, ..., v_k) \big) = \frac{ \vert det( v_1-v_0, ..., v_k-v_0 ) \vert }{n!} \f$
  double computeSimplexVolume(const std::set< boost::shared_ptr<Generator_Rn> >& listOfSimplexVertices) const;

  /// Called by computeSimplexVolume() to compute a square matrix determinant.
  /// As we run it on small matrices we just use the minors method.
  double determinant(boost::numeric::ublas::matrix<double> a) const;

  void dump(std::ostream &this_ostream) {
    std::cout << "List of simplices:" << std::endl;
    unsigned int counter=0;
    std::vector< PolytopeToSimplexes >::iterator iteSet;
    {for (iteSet=_allSimplices.begin(); iteSet!=_allSimplices.end(); ++iteSet) {
      this_ostream << "Simplex number " << counter << std::endl;
      std::vector< unsigned int >::const_iterator iteVtx;
      {for (iteVtx=iteSet->getListOfVerticesToSplit().begin(); iteVtx!=iteSet->getListOfVerticesToSplit().end(); ++iteVtx) {
        boost::shared_ptr<Generator_Rn> gn = _polytope->getGenerator(*iteVtx);
        gn->dump( this_ostream );
      }}
      this_ostream << std::endl;
      counter++;
    }}
  }

  void dumpDS(std::ostream &this_ostream) const {
    std::cout << "Vertices By Facets:" << std::endl;
    std::vector< std::vector< unsigned int > >::const_iterator iteSet;
    unsigned int counter=0;
    {for (iteSet=_verticesByFacets.begin(); iteSet!=_verticesByFacets.end(); ++iteSet) {
      this_ostream << counter << ": ";
      std::copy(iteSet->begin(), iteSet->end(), std::ostream_iterator<unsigned int>(std::cout, " "));
      this_ostream << std::endl;
      counter++;
    }}
    std::cout << "Facets By Vertices:" << std::endl;
    counter=0;
    {for (iteSet=_facetsByVertices.begin(); iteSet!=_facetsByVertices.end(); ++iteSet) {
      this_ostream << counter << ": ";
      std::copy(iteSet->begin(), iteSet->end(), std::ostream_iterator<unsigned int>(std::cout, " "));
      this_ostream << std::endl;
      counter++;
    }}
  }

  void dumpAllSimplices(std::ostream &this_ostream) const {
    std::vector< PolytopeToSimplexes >::const_iterator itS;
    {for (itS=_allSimplices.begin(); itS!=_allSimplices.end(); ++itS) {
      itS->dump(this_ostream);
      this_ostream << std::endl;
    }}
  }


 protected:
  /// The ordered list of all vertices stored by facets
  std::vector< std::vector< unsigned int > > _verticesByFacets;
  /// The ordered list of all facets stored by vertices
  std::vector< std::vector< unsigned int > > _facetsByVertices;
  /// List to store all the simplices partitioning the polytope
  std::vector< PolytopeToSimplexes > _allSimplices;
  /// The current polytope we are working on.
  boost::shared_ptr<Polytope_Rn> _polytope;
  /// The volume of the polytope
  double _volume;
  /// As the algorithm goes down in lower dimensions, we want to store the starting space dimension.
  unsigned int _dimension;
  /// The number of facets of the current polytope
  unsigned int _numberOfFacets;
  /// The number of vertices of the current polytope
  unsigned int _numberOfVertices;

};

#endif // VOLUME_POLYTOPES_Rn
