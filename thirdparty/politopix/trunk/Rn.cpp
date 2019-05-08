// politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2015 : Delos Vincent
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
/// \file Rn.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include "Rn.h"


unsigned int Rn::_dimension = 6;
double Rn::_tolerance = 1.e-06;

void Rn::setDimension(unsigned int dim) {_dimension = dim;}

unsigned int Rn::getDimension() {return _dimension;}

double Rn::getTolerance() {return _tolerance;}

void Rn::setTolerance(double t) {_tolerance = t;}

