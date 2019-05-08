// politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2014-2015 : Delos Vincent
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
/// \file Tracking.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef INC_Tracking
#define INC_Tracking


enum StatusBefore {
  Operator_UNCHANGED = 0,
  Operator_MODIFIED = 1,
  Operator_DELETED = 2,
  Operator_UNKNOWN = 4
};

enum StatusAfter {
  Result_UNCHANGED = 0,
  Result_MODIFIED = 1,
  Result_CREATED = 2,
  Result_UNKNOWN = 4
};


/// \brief This class stores static function that dispatch the main geometric values we use.
class TrackingOperatorToResult {

 public:

  TrackingOperatorToResult() {}

  void setNumbersOfEntities(unsigned int nbEntBefore, unsigned int nbEntAfter) {
    _operator.resize(nbEntBefore);
    _result.resize(nbEntAfter);
    // Set DELETED by default as we will change the status further on.
    for (unsigned int i=0;i<_operator.size();++i)
      _operator[i].first = Operator_DELETED;
    for (unsigned int i=0;i<_result.size();++i)
      _result[i].first = Result_UNCHANGED;
  }

  /// Mark as Operator_UNCHANGED the nb-th entity before the operation
  void setOperatorEntityAsUnchanged(unsigned int nb) {_operator[nb].first = Operator_UNCHANGED;}

  /// Mark as Operator_MODIFIED the nb-th entity before the operation
  void  setOperatorEntityAsModified(unsigned int nb) {_operator[nb].first = Operator_MODIFIED;}

  /// Mark as Operator_DELETED the nb-th entity before the operation
  void   setOperatorEntityAsDeleted(unsigned int nb) {_operator[nb].first = Operator_DELETED;}

  /// Mark as Result_UNCHANGED the nb-th entity before the operation
  void setResultEntityAsUnchanged(unsigned int nb) {_result[nb].first = Result_UNCHANGED;}

  /// Mark as Result_MODIFIED the nb-th entity after the operation
  void  setResultEntityAsModified(unsigned int nb) {_result[nb].first = Result_MODIFIED;}

  /// Mark as Result_UNKNOWN the nb-th entity after the operation
  void   setResultEntityAsUnknown(unsigned int nb) {_result[nb].first = Result_UNKNOWN;}

  /// Mark as Result_CREATED the nb-th entity after the operation
  void   setResultEntityAsCreated(unsigned int nb) {_result[nb].first = Result_CREATED;}

  /// Get the nb-th entity status after the operation
  StatusAfter    getResultEntityStatus(unsigned int nb) {return   _result[nb].first;}

  /// Get the nb-th entity status after the operation
  StatusBefore getOperatorEntityStatus(unsigned int nb) {return _operator[nb].first;}


  /// Set the link between the nb-th entity of operator1 and the nbRes-th entity of the result.
  void setOperator_Result(unsigned int nb, int nbRes) {_operator[nb].second = nbRes;}

  /// Set the link between the nb-th entity of operator1 and the nbRes-th entity of the result.
  void setResult_Operator(unsigned int nbRes, int nb) {_result[nbRes].second = nb;}

  ///
  const std::vector< std::pair< StatusBefore, int > >& getOperatorToResult() const {return _operator;}

  ///
  const std::vector< std::pair< StatusAfter,  int > >& getResultToOperator() const {return _result;}

 protected:

  /// List of the entities before the given operation with their status and the entity number it is connected to after the operation (if need be).
  std::vector< std::pair< StatusBefore, int > > _operator;
  /// List of the entities after the given operation with their status and the entity number it is connected to before the operation (if need be).
  std::vector< std::pair< StatusAfter,  int > > _result;

};


typedef std::pair< int, int > operator1_operator2;

/// \brief This class stores static function that dispatch the main geometric values we use.
class TrackingBinaryOperation {

 public:

  TrackingBinaryOperation() {}

  void setNumbersOfEntities(unsigned int nbEntBefore1, unsigned int nbEntBefore2, unsigned int nbEntAfter) {
    _operator1.resize(nbEntBefore1);
    _operator2.resize(nbEntBefore2);
    _result.resize(nbEntAfter);
    {for (unsigned int i=0;i<_operator1.size();++i) {
      _operator1[i].first = Operator_DELETED;
      // Put -1 at the moment as we don't know the corresponding element.
      _operator1[i].second = -1;
    }}
    {for (unsigned int i=0;i<_operator2.size();++i) {
      _operator2[i].first = Operator_DELETED;
      _operator2[i].second = -1;
    }}
    {for (unsigned int i=0;i<_result.size();++i) {
      _result[i].first = Result_UNCHANGED;
      _result[i].second.first  = -1;
      _result[i].second.second = -1;
    }}
  }

  /// Set the link between the nbOp1-th entity of operator1 and the nbRes-th entity of the result.
  void setOperator1_Result(unsigned int nbOp1, int nbRes) {_operator1[nbOp1].second = nbRes;}

  /// Set the link between the nbOp2-th entity of operator2 and the nbRes-th entity of the result.
  void setOperator2_Result(unsigned int nbOp2, int nbRes) {_operator2[nbOp2].second = nbRes;}

  /// Set the link between the nbRes-th entity of the result and the nbOp1-th entity of operator1 and the nbOp2-th entity of operator2.
  void setResult_Operator1Operator2(unsigned int nbRes, int nbOp1, int nbOp2) {
    _result[nbRes].second.first  = nbOp1;
    _result[nbRes].second.second = nbOp2;
  }

  ///
  const std::vector< std::pair< StatusAfter, operator1_operator2 > >& getResultStatus_Operator1Operator2() const { return _result;}

  /// Mark Operator_UNCHANGED the nbOp1-th entity of operator1
  void setOperator1EntityAsUnchanged(unsigned int nbOp1) {_operator1[nbOp1].first = Operator_UNCHANGED;}

  /// Mark Operator_MODIFIED the nbOp1-th entity of operator1
  void setOperator1EntityAsModified( unsigned int nbOp1) {_operator1[nbOp1].first = Operator_MODIFIED;}

  /// Mark Operator_DELETED the nbOp1-th entity of operator1
  void setOperator1EntityAsDeleted(  unsigned int nbOp1) {_operator1[nbOp1].first = Operator_DELETED;}

  /// Mark Operator_UNCHANGED the nbOp2-th entity of operator2
  void setOperator2EntityAsUnchanged(unsigned int nbOp2) {_operator2[nbOp2].first = Operator_UNCHANGED;}

  /// Mark Operator_MODIFIED the nbOp2-th entity of operator2
  void setOperator2EntityAsModified( unsigned int nbOp2) {_operator2[nbOp2].first = Operator_MODIFIED;}

  /// Mark Operator_DELETED the nbOp2-th entity of operator2
  void  setOperator2EntityAsDeleted( unsigned int nbOp2) {_operator2[nbOp2].first = Operator_DELETED;}

  /// Mark Result_UNCHANGED the nbRes-th entity of the result
  void setResultEntityAsUnchanged(unsigned int nbRes) { _result[nbRes].first = Result_UNCHANGED;}

  /// Mark Result_MODIFIED the nbRes-th entity of the result
  void setResultEntityAsModified( unsigned int nbRes) { _result[nbRes].first = Result_MODIFIED;}

  /// Mark Result_UNKNOWN the nbRes-th entity of the result
  void  setResultEntityAsUnknown( unsigned int nbRes) { _result[nbRes].first = Result_UNKNOWN;}

  /// Mark Result_CREATED the nbRes-th entity of the result
  void  setResultEntityAsCreated( unsigned int nbRes) { _result[nbRes].first = Result_CREATED;}


 protected:
  /// List of the entities before the given operation with their status and the entity number it is connected to after the operation (if need be).
  std::vector< std::pair< StatusBefore, int > > _operator1;
  /// List of the entities before the given operation with their status and the entity number it is connected to after the operation (if need be).
  std::vector< std::pair< StatusBefore, int > > _operator2;
  /// List of the entities after the given operation with their status and the entity numbers it is connected to before the operation (if need be).
  std::vector< std::pair< StatusAfter, operator1_operator2 > > _result;

};

#endif // INC_Tracking
