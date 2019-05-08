// politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2011-2015 : Delos Vincent
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
/// \file NormalFan_Rn.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef NORMALFAN_Rn
#define NORMALFAN_Rn

#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <vector>
#include "Rn.h"
#include "Generator_Rn.h"
#include "HalfSpace_Rn.h"
#include "Polytope_Rn.h"


/// \class NormalFan_Rn
/// \brief Model a normal fan.
class polito_EXPORT NormalFan_Rn {
  friend class constIteratorOfListOfPolyhedralCones;

 public:
  /// Constructor
  NormalFan_Rn() {}

  /// Constructor
  NormalFan_Rn(const boost::shared_ptr<Polytope_Rn>& A);

  /// Destructor
  ~NormalFan_Rn() {}

  /// Get the total number of polyhedral cones.
  unsigned int numberOfPolyhedralCones() const {return _listOfPolyhedralCones.size();}

  /// Add the current half-space in its list.
  void addPolyhedralCone(boost::shared_ptr<PolyhedralCone_Rn> hs) {_listOfPolyhedralCones.push_back(hs);}

  /// Add the current vertex in its list.
  void addVertex(boost::shared_ptr<Generator_Rn> vx) {_listOfVertices.push_back(vx);}

  const std::vector< boost::shared_ptr<Generator_Rn> >& getListOfGenerators() const {return _listOfVertices;}

  const std::vector< boost::shared_ptr<PolyhedralCone_Rn> >& getListOfPolyhedralCones() const {return _listOfPolyhedralCones;}


  // CHECK POLYHEDRON

  bool checkTopologyAndGeometry() const throw (std::domain_error);


  // ALGORITHMS

  /// \brief Compute \f$ N(A+B) = N(A) \wedge N(B) \f$
  ///
  /// Compute the intersection of all polyhedral cones from the first normal fan N(A) with all polyhedral cones from the second normal fan N(B).<br>
  /// \f[ N(A+B) = N(A) \wedge N(B) = \left\{ C_{a_i} \bigcap C_{b_j}, \forall \, C_{a_i} \in N(A), \forall \, C_{b_j} \in N(B) \right\} \f]
  /// \param NA The first  normal fan computed from polytope A
  /// \param NB The second normal fan computed from polytope B
  void computeCommonRefinement(const NormalFan_Rn& NA, const NormalFan_Rn& NB);

  void computeHyperplanesSeparationForProjection(const std::vector< boost::shared_ptr<HalfSpace_Rn> >&, boost::shared_ptr<Polytope_Rn>& );

  /// Dump the polyhedral structure on std::cout.
  void dump(std::ostream& out) const;

 protected:
  /// The list of polyhedral cones partitioning the whole space.
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> > _listOfPolyhedralCones;
  /// The list of vertices attached to their respective dual polyhedral cones.
  std::vector< boost::shared_ptr<Generator_Rn> > _listOfVertices;
};

#endif // NORMALFAN_Rn
