// politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2012-2015 : Delos Vincent
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
/// \file GeometricObjectIterator_Rn.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef GEOMETRIC_OBJECT_ITERATOR_Rn
#define GEOMETRIC_OBJECT_ITERATOR_Rn

#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <exception>
#include <vector>
#include <stack>
#include <set>
#include "Generator_Rn.h"
#include "HalfSpace_Rn.h"
#include "Rn.h"

using namespace boost::numeric::ublas;


/// \brief This class is designed to contain the list of all generators
/// or half-spaces representing a polytope or a polyhedral cone.
template< class GEOMETRIC_OBJECT > class listOfGeometricObjects {

 public:
  /// Constructor.
  listOfGeometricObjects() {}

  /// Include a new half space in the list.
  void push_back(const GEOMETRIC_OBJECT& gn) {_GO.push_back(gn);}

  /// Return the i-th generator.
  const GEOMETRIC_OBJECT& operator [](unsigned int i) const {return _GO[i];}

  /// Copies all elements from listOfGN to _GN
  void operator=(const listOfGeometricObjects< GEOMETRIC_OBJECT >& listOfGN)
  {_GO.clear(); _GO = listOfGN._GO;}

  /// Copies all elements from listOfGN to _GN
  void assign(const listOfGeometricObjects< GEOMETRIC_OBJECT >& listOfGN)
  {_GO.clear(); _GO = listOfGN._GO;}

  /// Check whether the set is empty or not.
  bool empty() const {return _GO.empty();}

  /// Get the total number of genuine facets.
  unsigned int size() const {return _GO.size();}

  /// Clear the whole list.
  void clear() {_GO.clear();}

  /// Multiply all generators or half-spaces by -1.
  void negate() {
    for (typename std::vector< GEOMETRIC_OBJECT >::iterator it=_GO.begin(); it!=_GO.end(); ++it)
      (*it)->negate();
  }

  /// Find a given object in list..
  unsigned int find(const GEOMETRIC_OBJECT& GO) const throw (std::out_of_range) {
    unsigned int index=0;
    typename std::vector< GEOMETRIC_OBJECT >::const_iterator it;
    for (it=_GO.begin(); it!=_GO.end() && *it!=GO; ++it)
      index++;
    if (it == _GO.end()) {
      std::ostringstream stream_;
      stream_ << "listOfGeometricObjects::find() The current object is not stored in array";
      std::string valString = stream_.str();
      throw std::out_of_range(valString);
    }
    return index;
  }

  /// Get rid of all the objects stored in the set.
  void removeGeometricObjects(const std::set< GEOMETRIC_OBJECT >& setToRemove) {
    for (typename std::vector< GEOMETRIC_OBJECT >::iterator it=_GO.begin(); it!=_GO.end(); ) {
      if (setToRemove.find(*it) != setToRemove.end())
        it = _GO.erase(it);
      else
        ++it;
    }
  }

  /// Remove the geometric object number <i> j </i> from the list.
  void removeGeometricObject(unsigned int j) {_GO.erase(_GO.begin()+j);}

  /// Tell whether a given object is declared inferior to another one.
  static bool inferior(const GEOMETRIC_OBJECT& HS1, const GEOMETRIC_OBJECT& HS2) {
    unsigned int i=0;
    unsigned int n=HS1->dimension();
    do {
      if (HS1->getCoefficient(i) < HS2->getCoefficient(i))
        return true;
      else if (HS1->getCoefficient(i) > HS2->getCoefficient(i))
        return false;
      else
        i++;
    } while (i < n);
    return true;
  }

  /// The opposite of the function inferior(HS1, HS2)
  static bool superior(const GEOMETRIC_OBJECT& HS1, const GEOMETRIC_OBJECT& HS2) {
    return !inferior(HS1, HS2);
  }

  void lexminSort(unsigned int step) {
    typename std::vector< GEOMETRIC_OBJECT >::iterator it = _GO.begin();
    std::advance(it, step);
    std::sort( it, _GO.end(), &inferior );
  }

  void lexmaxSort(unsigned int step) {
    typename std::vector< GEOMETRIC_OBJECT >::iterator it = _GO.begin();
    std::advance(it, step);
    std::sort( it, _GO.end(), &superior );
  }

 protected:
  /// The full list of half spaces or generators for example.
  std::vector< GEOMETRIC_OBJECT > _GO;
};


/// \brief This class is designed to run the list of all geometric objects representing a polytope.
template< class GEOMETRIC_OBJECT > class constIteratorOfListOfGeometricObjects {

 public:
  /// Constructor
  constIteratorOfListOfGeometricObjects(const listOfGeometricObjects<GEOMETRIC_OBJECT>& l):_list(l)
  {_iterator=0;_step=0;}

  /// Move the iterator at the beginning of the list.
  void begin() {_iterator=_step;}

  /// Move the iterator one step forward.
  void next() {_iterator++;}

  /// Step forward in the list geometric elements.
  void setStep(unsigned int n) {_step=n; advance(n);}

  /// Step forward in the list geometric elements.
  void advance(unsigned int n) {
    if (n > _list.size()) return;
    _iterator = _iterator+n;
    //while (_iterator < _list._fullSize && counter < n) {
      //_iterator++;
      //if (_list._activeFacets[_iterator] == true) {
        //counter++;
      //}
    //}
  }

  /// Tell whether we have reached the end of the list.
  bool end() const {return (_iterator==_list.size());}

  /// Return the current geometric element.
  const GEOMETRIC_OBJECT current() {
    /// Return the current geometric element.
    if (_iterator < _list.size())
      return _list[_iterator];
    else
      return GEOMETRIC_OBJECT();
  }

  /// Return the current position in the list.
  int currentIteratorNumber() const {return _iterator;}

 protected:
  /// The current position in the list.
  unsigned int _iterator;
  /// The actual list of geometric elements.
  const listOfGeometricObjects<GEOMETRIC_OBJECT>& _list;
  /// To perform a step.
  unsigned int _step;
};

/// Insert the half-spaces in the list in a lexicographically order, whether min or max.
template< class GEOMETRIC_OBJECT > class lexIteratorOfListOfGeometricObjects {

 public:
  /// Constructor
  lexIteratorOfListOfGeometricObjects(listOfGeometricObjects<GEOMETRIC_OBJECT>& l):
    _iterator(0),_list(l),_alreadySorted(false),_step(0)
    {}

  /// Move the iterator at the beginning of the list.
  void begin() {_iterator=_step;}

  /// Move the iterator one step forward.
  void next() {_iterator++;}

  /// Step forward in the list geometric elements.
  void setStep(unsigned int n) {_step=n;}

  /// Tell whether we have reached the end of the list.
  bool end() const {return (_iterator==_list.size());}

  /// Return the current geometric element.
  const GEOMETRIC_OBJECT current() {
    if (_iterator < _list.size())
      return _list[_iterator];
    else
      return GEOMETRIC_OBJECT();
  }

  /// Return the current position in the list.
  int currentIteratorNumber() const {return _iterator;}

 protected:
  /// The current position in the list.
  unsigned int _iterator;
  /// The actual list of geometric elements.
  listOfGeometricObjects<GEOMETRIC_OBJECT>& _list;
  /// Not to do the same job twice.
  bool _alreadySorted;
  /// Sort after a step.
  unsigned int _step;

};

/// Insert the half-spaces in the list in lexicographically increasing order.
template< class GEOMETRIC_OBJECT > class lexminIteratorOfListOfGeometricObjects :
  public lexIteratorOfListOfGeometricObjects< GEOMETRIC_OBJECT > {

 public:
  /// Constructor
  lexminIteratorOfListOfGeometricObjects(listOfGeometricObjects<GEOMETRIC_OBJECT>& l):
    lexIteratorOfListOfGeometricObjects< GEOMETRIC_OBJECT >(l) {}

  /// Step forward in the list geometric elements.
  void setStep(unsigned int n) {
    this->_step = n;
    this->_iterator = 0;
    if (this->_alreadySorted == false) {
      this->_list.lexminSort(this->_step);
      this->_alreadySorted = true;
    }
    this->_iterator = this->_step;
  }

  /// Move the iterator at the beginning of the list.
  void begin() {this->_iterator = this->_step;}

};

/// Insert the half-spaces in the list in lexicographically decreasing order.
template< class GEOMETRIC_OBJECT > class lexmaxIteratorOfListOfGeometricObjects :
    public lexIteratorOfListOfGeometricObjects< GEOMETRIC_OBJECT > {

 public:
  /// Constructor
  lexmaxIteratorOfListOfGeometricObjects(listOfGeometricObjects<GEOMETRIC_OBJECT>& l):
    lexIteratorOfListOfGeometricObjects< GEOMETRIC_OBJECT >(l) {}

  /// Move the iterator at the beginning of the list.
  void setStep(unsigned int n) {
    this->_step = n;
    this->_iterator = 0;
    if (this->_alreadySorted == false) {
      this->_list.lexmaxSort(this->_step);
      this->_alreadySorted = true;
    }
    this->_iterator = this->_step;
  }

  /// Move the iterator at the beginning of the list.
  void begin() {this->_iterator = this->_step;}

};


#endif // GEOMETRIC_OBJECT_ITERATOR_Rn
