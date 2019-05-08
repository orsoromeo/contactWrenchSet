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
/// \file Generator_Rn.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include "Generator_Rn.h"

Generator_Rn::Generator_Rn(unsigned int n) {
  _coordinates.resize(n);
}

Generator_Rn::~Generator_Rn() {
}

void Generator_Rn::removeFacet(unsigned int j) throw (std::out_of_range,std::domain_error) {
  if (numberOfFacets() == 0) {
    std::string errorMessage("List of neighbour vertices not allocated in Generator_Rn::removeNeighbourVertex");
    throw std::domain_error(errorMessage);
  }
  else if (j >= numberOfFacets()) {
    std::string errorMessage = Point_Rn::concatStrings(j, "Generator_Rn::removeNeighbourVertex");
    throw std::out_of_range(errorMessage);
  }
  std::vector< boost::shared_ptr<HalfSpace_Rn> >::iterator itRemove = _supportFacets.begin() + j;
  // Now, actually remove this element from the array
  _supportFacets.erase(itRemove);

}

bool Generator_Rn::isFacetInside(boost::shared_ptr<HalfSpace_Rn> F) const {
  std::vector< boost::shared_ptr<HalfSpace_Rn> >::const_iterator it = _supportFacets.begin();
  {for (; it!=_supportFacets.end(); ++it) {
    if (*it == F)
      return true;
  }}
  return false;
}
