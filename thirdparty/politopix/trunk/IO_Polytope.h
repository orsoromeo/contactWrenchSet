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
/// \file IO_Polytope.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef IO_POLYTOPES
#define IO_POLYTOPES

#include <stdexcept>
#include <exception>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "polito_Export.h"
#include "Polytope_Rn.h"
#include "PolyhedralCone_Rn.h"

/// \brief Read/write polytopes. <br>
/// The way we store polytopes : <br>
/// 1st  line : comments = "# Dimension NumberOfHalfspaces NumberOfGenerators" <br>
/// 2nd  line : cartesian_space_dimension number_of_facets number_of_generators <br>
/// 3rd  line : comments = "# HALFSPACES : a0 + a1.x1 + ... + an.xn >= 0." <br>
/// 4th  line : a00 a10 ... an0 <br>
/// k-th line : a0k a1k ... ank <br>
/// (k+1)-th line : comments = "# GENERATORS : V = (v1, ..., vn)" <br>
/// (k+2)-th line : v11 ... v1n <br>
/// l-th line : vl1 ... vln <br>
/// (l+1)-th line : comments = "# FACETS PER GENERATOR : {Fi1, Fi2, ...}" <br>
/// (l+2)-th line the neigh : Fr ... Fs <br>
/// m-th line : Fs ... Ft <br>
/// If (number_of_vertices == 0) then compute the vertices from the <br>
/// facets including the polytope into a huge cube containing it. <br>
/// In this case the blocks "GENERATORS" and "FACETS PER GENERATOR" are ignored.
class IO_Polytope {

 public:
  IO_Polytope() {}

  ~IO_Polytope() {}

  /// Load the main data format to store polytopes.
  static polito_EXPORT void load(const std::string& filename, boost::shared_ptr<PolyhedralCone_Rn> POLY)
    throw (std::ios_base::failure,std::out_of_range);

  /// Save the polytope to the main data format.
  static polito_EXPORT void save(const std::string& filename, boost::shared_ptr<PolyhedralCone_Rn> POLY)
    throw (std::ios_base::failure);

};

#endif // IO_POLYTOPES
