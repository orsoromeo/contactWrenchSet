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
/// \file Point_Rn.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef POINT_Rn
#define POINT_Rn

#include <algorithm>
#include <stdexcept>
#include <exception>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "polito_Export.h"

using namespace boost::numeric::ublas;

/// Creation of a n-coordinate geometric point designed to be shared by its neighbour faces.
class polito_EXPORT Point_Rn {

 public:
  // Think about in the future to an empty constructor
  // reading the dimension stored in a static variable.
  /// Create a n-coordinates point.
  Point_Rn(unsigned int n);

  /// Fill the n-dimensional point with the constant u.
  Point_Rn(unsigned int n, double u);

  /// Create a 3-dimensional point.
  Point_Rn(double u1, double u2, double u3);

  virtual ~Point_Rn();

  double normalize();

  double distanceFrom(const Point_Rn&);

  void setCoordinate(unsigned int i, double val) throw (std::out_of_range);

  double getCoordinate(unsigned int i) const throw (std::out_of_range);

  int dimension() const {return _coordinates.size();}

  void load(std::istream &this_istream);

  void save(std::ostream &this_ostream) const;

  vector<double>::const_iterator begin() const {return _coordinates.begin();}

  vector<double>::const_iterator  end()  const {return _coordinates.end();}

  const vector<double>& vect() const {return _coordinates;}

  void negate() { _coordinates = -1.*_coordinates; }

  /// Useful function to provide error message to the exception mechanism.
  static std::string concatStrings(int i, const std::string& functionName);

  /// Useful function to provide error message to the exception mechanism.
  static std::string concatStrings(int i, double val, const std::string& functionName);

 protected:
  // Create an empty set of coordinates.
  vector<double> _coordinates;

};


#endif // POINT_Rn
