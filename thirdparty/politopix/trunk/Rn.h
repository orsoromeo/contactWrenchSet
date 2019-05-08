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
/// \file Rn.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef INC_Rn
#define INC_Rn

#include "polito_Export.h"

/// \brief This class stores static function that dispatch the main geometric values we use.
class Rn {

 public:
  Rn() {}

  /// Set the dimension for the cartesian space we work in.
  static polito_EXPORT void setDimension(unsigned int dim);

  /// Return the dimension of the cartesian space we work in.
  static polito_EXPORT unsigned int getDimension();

  /// Give the minimum distance between two points.
  static polito_EXPORT double getTolerance();

  /// Give the minimum distance between two points.
  static polito_EXPORT void setTolerance(double t);

 protected:
   /// Rn dimension
   static unsigned int _dimension;
   /// Rn dimension
   static double _tolerance;

};

#endif // INC_Rn
