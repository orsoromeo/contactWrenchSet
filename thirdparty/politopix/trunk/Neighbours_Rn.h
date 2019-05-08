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
/// \file Neighbours_Rn.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef NEIGHBOURS_Rn
#define NEIGHBOURS_Rn
#include "HalfSpace_Rn.h"


/// \brief Class dedicated to degeneration processing when looking for neighbours.<br>
/// Let <i>A</i> be a polytope of \f$ \mathbb{R}^n, A = H_1^+ \cap H_2^+ \cap ... H_{n-1}^+ \cap ... H_k^+ \f$
/// where <i>k>n</i>.<br>
/// Let <i> [v<sub>i</sub>, v<sub>j</sub>] </i> be a segment of two vertices of <i>A</i> such as:
/// <ul>
/// <li> they share (n-1) hyperplanes in common: \f$ H_{ij} = \{ H_1, H_2, ... H_{n-1} \} \f$ </li>
/// <li> there exists an hyperplane \f$ \mathcal{H} \f$ separating <i> [v<sub>j</sub>, v<sub>k</sub>] </i> </li>
/// </ul>
/// We call <i> [v<sub>i</sub>, v<sub>j</sub>] </i> a pseudo-edge if it respects the first assumption.
/// The question is: to which condition <i> [v<sub>i</sub>, v<sub>j</sub>] </i> is a genuine edge?
/// Can we answer the question only by processing the pseudo-edges separated by the hyperplane \f$ \mathcal{H} \f$?
/// The straight line \f$ (v_i, v_j) \subset H_1 \cap H_2 \cap ... H_{n-1} \text{~with~} \mathcal{H} \neq H_u \f$,
/// \f$ u \in \{ 1, ..., n-1 \} \f$. So \f$ H_1 \cap H_2 \cap ... H_{n-1} \cap H_n^+ \cap ... H_k^+ = F_A \f$
/// is a face of <i>A</i> of dimension at least 1.
/// Let's assume <i> [v<sub>i</sub>, v<sub>j</sub>] </i> is not an edge of <i>A</i>,
/// then it is not an edge of <i>F<sub>A</sub></i> and we cannot have <i>F<sub>A</sub></i> included neither in
/// \f$ \mathcal{H}^+ \f$ nor in \f$ \mathcal{H}^- \f$ as \f$ \mathcal{H} \f$
/// separates <i> [v<sub>i</sub>, v<sub>j</sub>] </i>.
/// So \f$ \mathcal{H} \f$ separates <i>F<sub>A</sub></i> and we can
/// find an edge <i> [v<sub>a</sub>, v<sub>b</sub>] </i> of
/// <i>F<sub>A</sub></i> such as:
/// <ul>
/// <li> \f$ \mathcal{H} \f$ separates <i> [v<sub>a</sub>, v<sub>b</sub>] </i> or </li>
/// <li> \f$ \mathcal{H} \f$ passes through <i> v<sub>a</sub> </i> or </li>
/// <li> \f$ \mathcal{H} \f$ passes through <i> v<sub>b</sub> </i> </li>
/// </ul>
/// As <i> [v<sub>a</sub>, v<sub>b</sub>] </i> is an edge of <i>F<sub>A</sub></i>,
/// <i> v<sub>a</sub> </i> and <i> v<sub>b</sub> </i> share in common the list of hyperplanes <i> H<sub>ab</sub> </i>.
/// <i> H<sub>ab</sub> </i> contains the (n-1) half-spaces
/// \f$ H_1, H_2, ... H_{n-1} \f$ and others because the intersection \f$ H_1 \cap H_2 \cap ... H_{n-1} \f$
/// does not define a straight line. So if <i> [v<sub>i</sub>, v<sub>j</sub>] </i> is not an edge of <i>A</i>, we can
/// find a genuine edge <i> [v<sub>a</sub>, v<sub>b</sub>] </i>, intersecting with \f$ \mathcal{H} \f$, such as
/// \f$ H_{ij} \subset H_{ab} \f$
class Neighbours_Rn {

public:
  Neighbours_Rn():_iterator(0) {}

  /// \brief Tell whether a pseudo neighbor is a genuine one comparing set of half-spaces.
  /// \param commonFacets the set of common half-spaces pointers between <i>this</i> and <i>gen</i>
  /// \param numbergenIN  the generator number candidate to be a genuine end of edge
  /// \param numbergenOUT the generator number candidate to be the other genuine end of edge
  /// \param state equal to HalfSpace_Rn::hs_IN_OR_OUT or HalfSpace_Rn::hs_ON according to the edge property
  void addGenerator(
    const std::vector< unsigned int >& commonFacets,
    unsigned int numbergenIN,
    unsigned int numbergenOUT,
    HalfSpace_Rn::State state) {

    //std::cout << std::endl << "DUMP CANDIDATE:" << std::endl;
    //std::cout << "Gin candidate: " << numbergenIN << " ";
    //std::cout << ", Gout candidate: " << numbergenOUT << std::endl;
    //dump(std::cout);
    bool alreadyInsertedPair=false;
    {for (unsigned int i=0; i<_HSPerNewGenerators.size(); ++i) {
      if (_HSPerNewGenerators[i].size() > commonFacets.size()) {
        // Larger set first.
        if (std::includes(_HSPerNewGenerators[i].begin(), _HSPerNewGenerators[i].end(), commonFacets.begin(), commonFacets.end())) {
          // No need to do anything as the current generator set of
          // half-spaces is included in the set number j.
          //std::cout << "false1 (size=" << i << ")" << std::endl;
          return;
        }
      }
      else {
        if (std::includes(commonFacets.begin(), commonFacets.end(), _HSPerNewGenerators[i].begin(), _HSPerNewGenerators[i].end())) {
          // Substitute the old generator to be removed with the current one.
          if (alreadyInsertedPair == false) {
            _HSPerNewGenerators[i]  = commonFacets;
            _GeneratorsInNumber[i]  = numbergenIN;
            _GeneratorsOutNumber[i] = numbergenOUT;
            _GeneratorsState[i] = state;
            alreadyInsertedPair = true;
          }
          else {
            // This case is rare but can occur, when we need to erase two generators already
            // in the list. Take care of not inserting twice the same so shift the arrays.
            {for (unsigned int ii=i; ii<_HSPerNewGenerators.size()-1; ++ii) {
              _GeneratorsOutNumber[ii]= _GeneratorsOutNumber[ii+1];
              _HSPerNewGenerators[ii] = _HSPerNewGenerators[ii+1];
              _GeneratorsInNumber[ii] = _GeneratorsInNumber[ii+1];
              _GeneratorsState[ii] = _GeneratorsState[ii+1];
            }}
            _GeneratorsOutNumber.pop_back();
            _GeneratorsInNumber.pop_back();
            _HSPerNewGenerators.pop_back();
            _GeneratorsState.pop_back();
          }
        }
      }
    }}
    if (alreadyInsertedPair == false) {
      _GeneratorsInNumber.push_back(numbergenIN);
      _GeneratorsOutNumber.push_back(numbergenOUT);
      _HSPerNewGenerators.push_back(commonFacets);
      _GeneratorsState.push_back(state);
      //std::cout << "true (total=" << _HSPerNewGenerators.size() << std::endl;
    }
    //dump(std::cout);
  }

  /// \brief Iterator function.
  void begin() {_iterator=0; checkIterator();}

  /// \brief Iterator function.
  void next() {++_iterator; checkIterator();}

  /// \brief  Make sure we don't point on a generator with state ON.
  void checkIterator() {
    while (_iterator!=_HSPerNewGenerators.size() && _GeneratorsState[ _iterator]==HalfSpace_Rn::hs_ON)
      ++_iterator;
  }

  /// \brief Iterator function.
  bool end() {return (_iterator==_HSPerNewGenerators.size());}

  /// \brief Iterator function.
  unsigned int currentGenInNumber()  {return _GeneratorsInNumber[ _iterator];}

  /// \brief Iterator function.
  unsigned int currentGenOutNumber() {return _GeneratorsOutNumber[_iterator];}

  /// \brief Display the content on the stream passed as an argument.
  void dump(std::ostream &ofs) {
    ofs << "New gen:" << std::endl;
    unsigned int i=0;
    std::vector< std::vector< unsigned int > >::const_iterator ite;
    for (ite=_HSPerNewGenerators.begin(); ite!=_HSPerNewGenerators.end(); ++ite) {
      ofs << "(Gin=  " << _GeneratorsInNumber[i]  << ",";
      ofs << " Gout= " << _GeneratorsOutNumber[i] << ") ";
      ofs << "(State=" << HalfSpace_Rn::getStateAsText(_GeneratorsState[i]) << ")" << std::endl;
      i++;
      std::copy((*ite).begin(), (*ite).end(), std::ostream_iterator< unsigned int >(ofs, " ") );
      ofs << std::endl;
    }
  }

protected:
  /// A runner to iterate through the list of genuine neighbors.
  unsigned int _iterator;
  /// The pair of generators state
  std::vector< HalfSpace_Rn::State > _GeneratorsState;
  /// The generator numbers  IN in a global list.
  std::vector< unsigned int > _GeneratorsInNumber;
  /// The generator numbers OUT in a global list.
  std::vector< unsigned int > _GeneratorsOutNumber;
  /// For each generator, store all raw pointers on their corresponding half-spaces.
  std::vector< std::vector< unsigned int > > _HSPerNewGenerators;
};

#endif
