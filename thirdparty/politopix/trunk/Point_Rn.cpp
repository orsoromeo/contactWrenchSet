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
/// \file Point_Rn.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cmath>
#include "Point_Rn.h"
#include "Rn.h"

Point_Rn::Point_Rn(unsigned int n) {
  _coordinates.resize(n);
}

Point_Rn::Point_Rn(unsigned int n, double u) {
  _coordinates.resize(n, u);
}

Point_Rn::Point_Rn(double u1, double u2, double u3) {
  _coordinates.resize(3);
  _coordinates[0] = u1;
  _coordinates[1] = u2;
  _coordinates[2] = u3;
}

Point_Rn::~Point_Rn() {
}

double Point_Rn::normalize() {
  double norm = norm_2(_coordinates);
  _coordinates = _coordinates / norm;
  return norm;
}

double Point_Rn::distanceFrom(const Point_Rn& P) {
  vector<double> distV = _coordinates;
  distV -= P._coordinates;
  double dist = norm_2(distV);
  return dist;
}

void Point_Rn::setCoordinate(unsigned int i, double val) throw (std::out_of_range) {
  if (i < _coordinates.size())
    _coordinates[i] = val;
  else {
    std::string errorMessage = Point_Rn::concatStrings(i, val, "Point_Rn::setCoordinate");
    throw std::out_of_range(errorMessage);
  }
}

double Point_Rn::getCoordinate(unsigned int i) const throw (std::out_of_range) {
  if (i < _coordinates.size())
    return _coordinates[i];
  else {
    std::string errorMessage = Point_Rn::concatStrings(i, "Point_Rn::getCoordinate");
    throw std::out_of_range(errorMessage);
  }
}

std::string Point_Rn::concatStrings(int i, const std::string& functionName) {
  std::ostringstream stream_;
  stream_ << "index=";
  stream_ << i;
  std::string valString = stream_.str();

  std::string errorMessage(functionName);
  errorMessage += " ";
  errorMessage += valString;
  errorMessage += " ";

  return errorMessage;
}

std::string Point_Rn::concatStrings(int i, double val, const std::string& functionName) {
  std::ostringstream stream_;
  stream_ << "index=";
  stream_ << i;
  stream_ << ", value=";
  stream_ << val;
  std::string valString = stream_.str();

  std::string errorMessage(functionName);
  errorMessage += " ";
  errorMessage += valString;
  errorMessage += " ";

  return errorMessage;
}

void Point_Rn::load(std::istream &this_istream) {
  double val;
  for (unsigned int coord_count=0; coord_count<_coordinates.size(); coord_count++) {
    this_istream >> val;
    setCoordinate(coord_count, val);
  }
}

void Point_Rn::save(std::ostream &this_ostream) const {
  for (unsigned int coord_count=0; coord_count<_coordinates.size(); coord_count++) {
    this_ostream << getCoordinate(coord_count) << " ";
  }
  //this_ostream << std::endl;
}

