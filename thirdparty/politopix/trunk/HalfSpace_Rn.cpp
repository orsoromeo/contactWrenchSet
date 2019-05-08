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
/// \file HalfSpace_Rn.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include "Point_Rn.h"
#include "HalfSpace_Rn.h"

HalfSpace_Rn::HalfSpace_Rn(unsigned int n):_coefficients(n),_constant(0.) {
}

HalfSpace_Rn::~HalfSpace_Rn() {
}

double HalfSpace_Rn::getCoefficient(unsigned int i) const throw (std::out_of_range) {
  if (i < _coefficients.size()) {
    return _coefficients[i];
  }
  else {
    std::string newErrorMessage = Point_Rn::concatStrings(i, "HalfSpace_Rn::getCoefficient");
    newErrorMessage += "\n";
    throw std::out_of_range(newErrorMessage);
  }
}

void HalfSpace_Rn::setCoefficient(unsigned int i, double c) throw (std::out_of_range) {
  if (i < _coefficients.size()) {
    _coefficients[i] = c;
  }
  else {
    std::string newErrorMessage = Point_Rn::concatStrings(i, "HalfSpace_Rn::setCoefficient");
    newErrorMessage += "\n";
    throw std::out_of_range(newErrorMessage);
  }
}

std::string HalfSpace_Rn::getStateAsText(const HalfSpace_Rn::State& state) {
  if (state == HalfSpace_Rn::hs_ON)
      return std::string("ON");
  if (state == HalfSpace_Rn::hs_OUT)
      return std::string("OUT");
  if (state == HalfSpace_Rn::hs_IN)
      return std::string("IN");
  if (state == HalfSpace_Rn::hs_IN_OR_OUT)
      return std::string("IN_OR_OUT");
  return std::string("UNKNOWN");
}
