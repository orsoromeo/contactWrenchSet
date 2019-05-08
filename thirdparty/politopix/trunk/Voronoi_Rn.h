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
/// \file Voronoi_Rn.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef VORONOI_Rn
#define VORONOI_Rn

#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <vector>
#include "Rn.h"
#include "Generator_Rn.h"
#include "HalfSpace_Rn.h"
#include "Polytope_Rn.h"


/// \class Voronoi_Rn
/// \brief Compute a n-dimensional Voronoi diagram.
/// It is a partitioning of a space into regions based on distance to points. Both the space and the list of points are provided as input.
class Voronoi_Rn {

 public:
  /// Constructor
  //Voronoi_Rn() { _inputSpace=0; }

  /// Constructor
  Voronoi_Rn(const boost::shared_ptr<Polytope_Rn>& inputSpace, const std::vector<Point_Rn>& listOfPoints);

  /// Destructor
  ~Voronoi_Rn() {}

  /// Run the whole algorithm
  bool compute() throw (std::length_error);

  /// Compute the half-space containing seed1, in between seed1 and seed2, according to the growing seed property.
  boost::shared_ptr<HalfSpace_Rn> computeMidPlane(std::vector<Point_Rn>::const_iterator seed1, std::vector<Point_Rn>::const_iterator seed2);

  const std::vector< boost::shared_ptr<Polytope_Rn> >& getVoronoiCells() const {return _listOfVoronoiCells;}
  std::vector< boost::shared_ptr<Polytope_Rn> > getVoronoiCells() {return _listOfVoronoiCells;}

  // CHECK POLYHEDRON
  bool checkTopologyAndGeometry() const throw (std::domain_error);

  /// Dump the cell structure on the given output.
  void dump(std::ostream& out) const;

  void gnuplot(std::ostream& out) const throw (std::domain_error);

 protected:
  /// The original space to be divided.
  const boost::shared_ptr<Polytope_Rn>& _inputSpace;
  /// The list of input points.
  const std::vector< Point_Rn >& _listOfSeeds;
  /// The list of polytopes partitioning the whole space.
  std::vector< boost::shared_ptr<Polytope_Rn> > _listOfVoronoiCells;
};

#endif // VORONOI_Rn
