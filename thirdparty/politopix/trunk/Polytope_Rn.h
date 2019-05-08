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
/// \file Polytope_Rn.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef POLYTOPE_Rn
#define POLYTOPE_Rn

#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <exception>
#include <vector>
#include "polito_Export.h"
#include "PolyhedralCone_Rn.h"
#include "Generator_Rn.h"
#include "HalfSpace_Rn.h"

/// Model a polytope using its two equivalent definitions : the convex hull and the half-space intersection.
class polito_EXPORT Polytope_Rn : public PolyhedralCone_Rn {

 public:
  /// Constructor for polytopes i.e. bounded convex polyhedra.
  Polytope_Rn():PolyhedralCone_Rn() {}

  /// Destructor
  virtual ~Polytope_Rn() {}

  /// Tell whether this polyhedron is bounded or not, polytopes are bounded.
  virtual bool isBounded() const {return true;}

  /// Two vertices are neighbours in a polytope <=> they share at least (n-1) facets.
  virtual unsigned int neigbourhoodCondition() const {return (dimension()-1);}

  /// Each facet in a polytope has got n vertices.
  virtual unsigned int numberOfGeneratorsPerFacet() const {return (dimension());}

  /// Initialize the truncating algorithm building a M-sized simplex around the polytope.
  virtual void createBoundingSimplex(double M);

  /// Initialize the truncating algorithm building a M-sized bounding box around the polytope.
  virtual void createBoundingBox(double M);

  virtual bool checkEdges() const;

  /// \brief Check whether two V-polytopes are identical
  /// Check whether the sets of vertices of A and B are equal.
  /// \param A The V-polytope
  /// \return true if \f[ \mathcal{V}_A = \mathcal{V}_B \f], false otherwise.
  bool checkEqualityOfVertices(const boost::shared_ptr<Polytope_Rn>& B, bool printOnScreen = false) const;

  /// This is the intersection vertex in the truncating algorithm.
  /// It is defined by the intersection between an edge, i.e. a 1-face, and an hyperplane, i.e. a (n-1)-face.
  virtual void createTruncatedGenerator(
      const boost::shared_ptr<Generator_Rn_SD>& in,
      const boost::shared_ptr<Generator_Rn_SD>& out,
      boost::shared_ptr<Generator_Rn_SD> newV, double ay, double az, double b=0.) const;

  /// Return the i-th primal cone C(i) of the polytope A.
  /// If A has k vertices then \f$ A = \displaystyle{ \bigcap_{i=1}^k C(i) }\f$ where \f$ C(i) = \displaystyle{ \bigcap_{i=1}^l \bar{H}_i^+ } \f$
  boost::shared_ptr<PolyhedralCone_Rn> getPrimalCone(unsigned int i) const throw (std::out_of_range);

  /// Return the primal cone C(vx) of the polytope A associated to the vertex vx.
  boost::shared_ptr<PolyhedralCone_Rn> getPrimalCone(const boost::shared_ptr<Generator_Rn>& vx) const;

};

#endif // POLYTOPE_Rn
