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
/// \file HalfSpace_Rn.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef HALFSPACE_Rn
#define HALFSPACE_Rn

#include <numeric>
#include <stdexcept>
#include <exception>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "Rn.h"

using namespace boost::numeric::ublas;


/// \brief A half-space whose frontier is a linear (n-1) dimension space. <br>
/// _constant + _coefficients[0].x1 + ... + _coefficients[n-1].xn >= 0.
class polito_EXPORT HalfSpace_Rn {

 public:

  enum State {
    hs_ON = 0,
    hs_IN = 1,
    hs_OUT = 2,
    hs_UNKNOWN = 3,
    hs_IN_OR_OUT = 4};

  /// Constructor.
  HalfSpace_Rn(unsigned int n);

  ~HalfSpace_Rn();

  double getCoefficient(unsigned int i) const throw (std::out_of_range);

   void  setCoefficient(unsigned int i, double c) throw (std::out_of_range);

   void  setConstant(double c) {_constant = c;}

  double getConstant() const {return _constant;}

  void negate() {_coefficients *= -1.;}

  int dimension() const {return _coefficients.size();}

  std::string getSideAsText() const {return std::string(">=");}

  boost::numeric::ublas::vector<double>::const_iterator begin() const {return _coefficients.begin();}

  boost::numeric::ublas::vector<double>::const_iterator  end()  const {return _coefficients.end();}

  const boost::numeric::ublas::vector<double>& vect() const {return _coefficients;}

  static std::string getStateAsText(const HalfSpace_Rn::State&);

  double computeDistancePointHyperplane(const boost::numeric::ublas::vector<double>& thisPoint) const {
    double halfSpaceNorm = std::inner_product(_coefficients.begin(), _coefficients.end(), _coefficients.begin(), 0.);
    halfSpaceNorm = sqrt(halfSpaceNorm);
    double scalarProduct = std::inner_product(thisPoint.begin(), thisPoint.end(), _coefficients.begin(), 0.);
    double distanceToHyperplane = (scalarProduct+_constant) / halfSpaceNorm;
    return distanceToHyperplane;
  }

  double computeDistancePointHyperplane(const boost::numeric::ublas::vector<double>& thisPoint, double halfSpaceNorm) const {
    double scalarProduct = std::inner_product(thisPoint.begin(), thisPoint.end(), _coefficients.begin(), 0.);
    double distanceToHyperplane = (scalarProduct+_constant) / halfSpaceNorm;
    return distanceToHyperplane;
  }

  double computeDistancePointHyperplane(
      const boost::numeric::ublas::vector<double>& thisPoint, 
      boost::numeric::ublas::vector<double>& projectedPoint, 
      double halfSpaceNorm) const {
    double scalarProduct = std::inner_product(thisPoint.begin(), thisPoint.end(), _coefficients.begin(), 0.);
    double distanceToHyperplane = (scalarProduct+_constant) / halfSpaceNorm;
    projectedPoint = thisPoint - distanceToHyperplane*_coefficients;
    return distanceToHyperplane;
  }

  void dump(std::ostream &this_ostream) const {
    this_ostream << "(" << getConstant() << ", ";
    unsigned int RnDIM = _coefficients.size();
    {for (unsigned int ii=0; ii<RnDIM; ++ii) {
      this_ostream << getCoefficient(ii);
      if (ii != RnDIM-1)
        this_ostream << ", ";
    }}
    this_ostream << ")";
  }

 protected:
  /// The normal vector
  boost::numeric::ublas::vector<double> _coefficients;
  /// The second member constant
  double _constant;

};

#endif // HALFSPACE_Rn
