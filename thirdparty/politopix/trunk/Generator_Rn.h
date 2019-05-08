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
/// \file Generator_Rn.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef GENERATOR_Rn
#define GENERATOR_Rn

#include <stdexcept>
#include <exception>
#include <vector>
#include <cmath>
#include <set>
#include "polito_Export.h"
#include "Rn.h"
#include "Tracking.h"
#include "Point_Rn.h"
#include "HalfSpace_Rn.h"


/// \brief A n-coordinates generator, which can be a vertex or an edge whether
/// it is contained by a polytope or a polyhedral cone. It contains all of its support facets.
class polito_EXPORT Generator_Rn {
  friend class Generator_Rn_SD;

 public:
  /// Creates a n-coordinates generator.
  Generator_Rn(unsigned int n);

  /// Copy constructor
  Generator_Rn(const Generator_Rn& gn) {
    _coordinates = gn._coordinates;
    _supportFacets = gn._supportFacets;
  }

  /// Destructor.
  ~Generator_Rn();

  int dimension() const {return _coordinates.size();}

  void setCoordinate(unsigned int i, double val) {_coordinates[i] = val;}

  double getCoordinate(unsigned int i) const {return _coordinates[i];}

  vector<double>::const_iterator begin() const {return _coordinates.begin();}

  vector<double>::const_iterator  end()  const {return _coordinates.end();}

  const vector<double>& vect() const {return _coordinates;}

  void setCoordinates(const vector<double>& vec) {_coordinates = vec;}

  void negate() { _coordinates = -1.*_coordinates; }

  bool isEqual1(const boost::shared_ptr<Generator_Rn>& gn, unsigned int RnDIM, double TOL2) {
    unsigned int i=0;
    double sum=0.;
    while (i < RnDIM) {
      sum = sum + fabs(getCoordinate(i) - gn->getCoordinate(i));
      if (sum > TOL2)
        return false;
      i++;
    }
    return true;
  }

  bool isEqual2(const boost::shared_ptr<Generator_Rn>& gn, unsigned int RnDIM, double TOL2) {
    unsigned int i=0;
    double sum=0.;
    while (i < RnDIM) {
      sum = sum + fabs(getCoordinate(i) - gn->getCoordinate(i)) * fabs(getCoordinate(i) - gn->getCoordinate(i));
      if (sum > TOL2)
        return false;
      i++;
    }
    return true;
  }

  /// Clear the list of facets.
  void clearFacets() {_supportFacets.clear();}

  /// Clear the list of facets.
  void switchFacets(const std::vector< boost::shared_ptr<HalfSpace_Rn> >& tab) {
    _supportFacets.clear();
    _supportFacets = tab;
  }

  /// Insert a new support facet for the current generator.
  void setFacet(boost::shared_ptr<HalfSpace_Rn> F) {_supportFacets.push_back(F);}

  /// Insert all facets stored in the argument.
  void importFacets(const std::set< boost::shared_ptr<HalfSpace_Rn> >& setOfFacets) {
    _supportFacets.clear();
    std::set< boost::shared_ptr<HalfSpace_Rn> >::const_iterator it;
    {for (it=setOfFacets.begin(); it!=setOfFacets.end(); ++it) {
      _supportFacets.push_back(*it);
    }}
  }

  /// Store all facets in a set.
  void exportFacets(std::set< boost::shared_ptr<HalfSpace_Rn> >& setOfFacets) const {
    std::vector< boost::shared_ptr<HalfSpace_Rn> >::const_iterator it;
    {for (it=_supportFacets.begin(); it!=_supportFacets.end(); ++it) {
      setOfFacets.insert(*it);
    }}
  }

  /// Remove the i-th facet in list.
  void removeFacet(unsigned int i) throw (std::out_of_range,std::domain_error);

  /// Return the i-th facet.
  /// The user has to check validity of the returned smart pointer.
  boost::shared_ptr<HalfSpace_Rn> getFacet(unsigned int i) const throw (std::out_of_range) {
    if (i < _supportFacets.size()) {
      return _supportFacets[i];
    }
    else {
      std::string errorMessage = Point_Rn::concatStrings(i, "Generator_Rn::getFacet");
      throw std::out_of_range(errorMessage);
    }
  }

  /// Return the i-th facet as a pointer for very fast comparisons. No check is performed!
  HalfSpace_Rn* getRawFacet(unsigned int i) {return _supportFacets[i].get();}

  /// Check whether the given half-space is inside the generator's list.
  bool isFacetInside(boost::shared_ptr<HalfSpace_Rn> F) const;

  /// Return the total number of support faces.
  unsigned int numberOfFacets() const {return _supportFacets.size();}

  void makeDiff(const boost::shared_ptr<Generator_Rn>& gn1, const boost::shared_ptr<Generator_Rn>& gn2) {
    _coordinates = gn1->_coordinates - gn2->_coordinates;
  }

  void makeSum( const boost::shared_ptr<Generator_Rn>& gn1, const boost::shared_ptr<Generator_Rn>& gn2) {
    _coordinates = gn1->_coordinates + gn2->_coordinates;
  }

  void makeCoefSum(
      const boost::shared_ptr<Generator_Rn>& gn1,
      const boost::shared_ptr<Generator_Rn>& gn2,
      double coef1,
      double coef2) {
    _coordinates = coef1*gn1->_coordinates + coef2*gn2->_coordinates;
  }

  /// Return the square distance of the generator gn1 to the straight line
  /// defined by _coordinates and passing through the origin.
  double getNormalDistance(const boost::shared_ptr<Generator_Rn>& gn1, double coef, unsigned int RnDIM) {
    double sum=0.;
    {for (unsigned int ii=0; ii<RnDIM; ++ii) {
      sum += (gn1->_coordinates[ii] - _coordinates[ii]*coef)*(gn1->_coordinates[ii] - _coordinates[ii]*coef);
    }}
    //std::cout << " (sum=" << sum << ") ";
    return sum;
  }

  double normalize() {
    double norm = norm_2(_coordinates);
    _coordinates = _coordinates / norm;
    return norm;
  }

  double distanceFrom(const Generator_Rn& P) {
    vector<double> distV = _coordinates;
    distV -= P._coordinates;
    double dist = norm_2(distV);
    return dist;
  }

  void dump(std::ostream &this_ostream) const {
    this_ostream << "[";
    //unsigned int RnDIM=Rn::getDimension();
    {for (unsigned int ii=0; ii<_coordinates.size(); ++ii) {
      this_ostream << getCoordinate(ii);
      if (ii != _coordinates.size()-1)
        this_ostream << ", ";
    }}
    this_ostream << "]";
  }

  void load(std::istream &this_istream) {
    double val;
    for (unsigned int coord_count=0; coord_count<_coordinates.size(); coord_count++) {
      this_istream >> val;
      setCoordinate(coord_count, val);
    }
  }

  void save(std::ostream &this_ostream) const {
    for (unsigned int coord_count=0; coord_count<_coordinates.size(); coord_count++) {
      this_ostream << getCoordinate(coord_count) << " ";
    }
    //this_ostream << std::endl;
  }


 protected:
  /// The set of coordinates.
  vector<double> _coordinates;
  /// Contain the list of all support facets.
  std::vector< boost::shared_ptr<HalfSpace_Rn> > _supportFacets;

};

/// \brief A n-coordinates generator for internal data structure.
/// It can be a vertex or an edge whether it is embedded in a polytope or
/// a polyhedral cone. It contains all of its support facets.
class Generator_Rn_SD {
  friend class Generator_Rn;

 public:

  enum Status {
    // When just made as a copy of a Generator_Rn
    UNCHANGED,
    // When some half-spaces have been modified
    MODIFIED,
    // When it was created as the result of an operation
    CREATED,
    // When it was created as the result of an operation and then some half-spaces have been modified
    CREATED_AND_MODIFIED,
    // When it was deleted by an operation
    DELETED,
    UNKNOWN
  };

  /// Creates a n-coordinates generator.
  Generator_Rn_SD(unsigned int n, unsigned int nb, Status st):_generatorNumber(nb),_status(st) {_coordinates.resize(n);}

  /// Copy constructor
  Generator_Rn_SD(const Generator_Rn_SD& gn):_generatorNumber(gn._generatorNumber),_status(gn._status) {
    _coordinates = gn._coordinates;
    _supportIntFacets = gn._supportIntFacets;
  }

  /// Constructor with a Generator_Rn
  Generator_Rn_SD(const Generator_Rn& gn, unsigned int nb, Status st):_generatorNumber(nb),_status(st) {
    _coordinates = gn._coordinates;
  }

  /// Destructor
  ~Generator_Rn_SD() {}

  int dimension() const {return _coordinates.size();}

  /// To make a Generator_Rn out of a Generator_Rn_SD
  boost::shared_ptr<Generator_Rn> makeGenerator_Rn() const {
    boost::shared_ptr<Generator_Rn> gn(new Generator_Rn( _coordinates.size() ));
    gn->_coordinates = _coordinates;
    return gn;
  }

  void setCoordinate(unsigned int i, double val) {_coordinates[i] = val;}

  double getCoordinate(unsigned int i) const {return _coordinates[i];}

  void setGeneratorNumber(unsigned int gn) {_generatorNumber = gn;}

  unsigned int getGeneratorNumber() const {return _generatorNumber;}

  void setStatus(Status st) {_status = st;}

  Status getStatus() const {return _status;}

  vector<double>::const_iterator begin() const {return _coordinates.begin();}

  vector<double>::const_iterator  end()  const {return _coordinates.end();}

  const vector<double>& vect() const {return _coordinates;}

  void negate() {_coordinates = -1.*_coordinates;}

  bool isEqual1(const boost::shared_ptr<Generator_Rn_SD>& gn, unsigned int RnDIM, double TOL2) {
    unsigned int i=0;
    double sum=0.;
    while (i < RnDIM) {
      sum = sum + fabs(getCoordinate(i) - gn->getCoordinate(i));
      if (sum > TOL2)
        return false;
      i++;
    }
    return true;
  }

  bool isEqual2(const boost::shared_ptr<Generator_Rn_SD>& gn, unsigned int RnDIM, double TOL2) {
    unsigned int i=0;
    double sum=0.;
    while (i < RnDIM) {
      sum = sum + fabs(getCoordinate(i) - gn->getCoordinate(i)) * fabs(getCoordinate(i) - gn->getCoordinate(i));
      if (sum > TOL2)
        return false;
      i++;
    }
    return true;
  }

  /// Insert a new support facet for the current generator.
  void setFacet(unsigned int F) {_supportIntFacets.push_back(F);}

  /// Insert a new support facet for the current generator.
  void setAllFacets(const std::vector<unsigned int>& AF) {_supportIntFacets = AF;}

  /// Insert all facets stored in the argument.
  void importFacets(const std::set< unsigned int >& setOfFacets) {
    _supportIntFacets.clear();
    std::set< unsigned int >::const_iterator it;
    {for (it=setOfFacets.begin(); it!=setOfFacets.end(); ++it) {
      _supportIntFacets.push_back(*it);
    }}
  }

  /// Store all facets in a set.
  void exportFacets(std::set< unsigned int >& setOfFacets) const {
    std::vector< unsigned int >::const_iterator it;
    {for (it=_supportIntFacets.begin(); it!=_supportIntFacets.end(); ++it) {
      setOfFacets.insert(*it);
    }}
  }

  /// Remove the i-th facet in list.
  void removeFacet(unsigned int i) throw (std::out_of_range,std::domain_error) {
    std::vector< unsigned int >::iterator itRemove = _supportIntFacets.begin() + i;
    // Now, actually remove this element from the array
    _supportIntFacets.erase(itRemove);
  }

  /// Return the i-th facet number.
  unsigned int getFacet(unsigned int i) const throw (std::out_of_range) {
    if (i < _supportIntFacets.size()) {
      return _supportIntFacets[i];
    }
    else {
      std::string errorMessage = Point_Rn::concatStrings(i, "Generator_Rn::getFacet");
      throw std::out_of_range(errorMessage);
    }
  }

  /// Return the i-th facet. No check is performed!
  unsigned int getRawFacet(unsigned int i) {return _supportIntFacets[i];}

  /// Check whether the given half-space is inside the generator's list.
  bool isFacetInside(unsigned int F) const {
    std::vector< unsigned int >::const_iterator it = _supportIntFacets.begin();
    {for (; it!=_supportIntFacets.end(); ++it) {
      if (*it == F)
        return true;
    }}
    return false;
  }

  void orderFacets() { std::sort(_supportIntFacets.begin(), _supportIntFacets.end());}

  std::vector<unsigned int>::const_iterator facetsBegin() const {return _supportIntFacets.begin();}

  std::vector<unsigned int>::const_iterator  facetsEnd()  const {return _supportIntFacets.end();}

  /// Return the total number of support faces.
  unsigned int numberOfFacets() const {return _supportIntFacets.size();}

  void makeDiff(const boost::shared_ptr<Generator_Rn_SD>& gn1, const boost::shared_ptr<Generator_Rn_SD>& gn2) {
    _coordinates = gn1->_coordinates - gn2->_coordinates;
  }

  void makeSum( const boost::shared_ptr<Generator_Rn_SD>& gn1, const boost::shared_ptr<Generator_Rn_SD>& gn2) {
    _coordinates = gn1->_coordinates + gn2->_coordinates;
  }

  void makeCoefSum(
      const boost::shared_ptr<Generator_Rn_SD>& gn1,
      const boost::shared_ptr<Generator_Rn_SD>& gn2,
      double coef1,
      double coef2) {
    _coordinates = coef1*gn1->_coordinates + coef2*gn2->_coordinates;
  }

  /// Return the square distance of the generator gn1 to the straight line
  /// defined by _coordinates and passing through the origin.
  double getNormalDistance(const boost::shared_ptr<Generator_Rn_SD>& gn1, double coef, unsigned int RnDIM) {
    double sum=0.;
    {for (unsigned int ii=0; ii<RnDIM; ++ii) {
      sum += (gn1->_coordinates[ii] - _coordinates[ii]*coef)*(gn1->_coordinates[ii] - _coordinates[ii]*coef);
    }}
    //std::cout << " (sum=" << sum << ") ";
    return sum;
  }

  double normalize() {
    double norm = norm_2(_coordinates);
    _coordinates = _coordinates / norm;
    return norm;
  }

  double distanceFrom(const Generator_Rn_SD& P) {
    vector<double> distV = _coordinates;
    distV -= P._coordinates;
    double dist = norm_2(distV);
    return dist;
  }

  void dump(std::ostream &this_ostream) const {
    this_ostream << "[";
    //unsigned int RnDIM=Rn::getDimension();
    {for (unsigned int ii=0; ii<_coordinates.size(); ++ii) {
      this_ostream << getCoordinate(ii);
      if (ii != _coordinates.size()-1)
        this_ostream << ", ";
    }}
    this_ostream << "]";
  }

  void load(std::istream &this_istream) {
    double val;
    for (unsigned int coord_count=0; coord_count<_coordinates.size(); coord_count++) {
      this_istream >> val;
      setCoordinate(coord_count, val);
    }
  }

  void save(std::ostream &this_ostream) const {
    for (unsigned int coord_count=0; coord_count<_coordinates.size(); coord_count++) {
      this_ostream << getCoordinate(coord_count) << " ";
    }
    //this_ostream << std::endl;
  }


 protected:
  /// The set of coordinates.
  vector<double> _coordinates;
  /// Contain the list of all support facets.
  std::vector< unsigned int > _supportIntFacets;
  /// The SD generator embeds its own number
  unsigned int _generatorNumber;
  /// The SD generator embeds its status to trace the operations
  Status _status;

};


#endif // GENERATOR_Rn
