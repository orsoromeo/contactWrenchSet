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
/// \file VolumeOfPolytopes_Rn.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include <string>
#include <iostream>
#include <boost/math/special_functions/factorials.hpp>
#include "VolumeOfPolytopes_Rn.h"

VolumeOfPolytopes_Rn::VolumeOfPolytopes_Rn(const boost::shared_ptr<Polytope_Rn> P) {
  _volume = -1.;
  _polytope = P;
  _dimension = Rn::getDimension();
  _numberOfFacets = P->numberOfHalfSpaces();
  _numberOfVertices = P->numberOfGenerators();

  // Fill the ordered list of vertices by facet, by example :
  // Facet0 : vertex_i, vertex_j, vertex_k, ...
  // Facet1 : vertex_l, vertex_m, ...
  // ...
  _verticesByFacets = std::vector< std::vector< unsigned int > >(_numberOfFacets);
  {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS(P->getListOfHalfSpaces());
    {for (iteHS.begin(); iteHS.end()!=true; iteHS.next()) {
      unsigned int facetNumber = iteHS.currentIteratorNumber();
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(P->getListOfGenerators());
      for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
        unsigned int vertexNumber = iteGN.currentIteratorNumber();
        {for (unsigned int j=0; j<iteGN.current()->numberOfFacets(); j++) {
          if (iteGN.current()->getFacet(j) == iteHS.current())
            _verticesByFacets[facetNumber].push_back(vertexNumber);
        }}
      }
    }}
  }

  // Fill the ordered list of facets by vertices
  // Vertex0 : facet_u, facet_v, facet_w, ...
  // Vertex1 : facet_x, facet_y, ...
  // ...
  _facetsByVertices = std::vector< std::vector< unsigned int > >(_numberOfVertices);
  {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(P->getListOfGenerators());
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
      unsigned int vertexNumber = iteGN.currentIteratorNumber();
      {for (unsigned int j=0; j<iteGN.current()->numberOfFacets(); j++) {
        unsigned int facetNumber = P->getHalfSpaceNumber( iteGN.current()->getFacet(j) );
        _facetsByVertices[vertexNumber].push_back(facetNumber);
      }}
    }}
  }

#ifdef DEBUG
  dumpDS(std::cout);
#endif
}

void VolumeOfPolytopes_Rn::check() const throw (std::domain_error) {
  //if (_listOfGN.size() == 0)
    //throw std::domain_error("VolumeOfPolytopes_Rn::runComputations() empty list of generators  !");
  //if (_listOfHS.size() == 0)
    //throw std::domain_error("VolumeOfPolytopes_Rn::runComputations() empty list of half-spaces !");
}

void VolumeOfPolytopes_Rn::splitCloudOfVertices(unsigned int DIM) {
  // The algorithm main data structure to store the computations
  std::stack< PolytopeToSimplexes > stackOfPolToSimp;
  std::vector< unsigned int > listOfVtx;
  // Build the first element to initialize the algorithm
  for (unsigned int i=0; i<_polytope->numberOfGenerators(); ++i)
    listOfVtx.push_back(i);
  PolytopeToSimplexes P2S(listOfVtx, DIM);
  stackOfPolToSimp.push(P2S);

  while (stackOfPolToSimp.size() != 0) {
#ifdef DEBUG
    unsigned int shift = _dimension-stackOfPolToSimp.top().dimension();
    std::cout << "The top of the stack is:" << std::endl;
    stackOfPolToSimp.top().dump(std::cout, shift);
#endif
    const PolytopeToSimplexes& topOfStack = stackOfPolToSimp.top();
    if (topOfStack.checkSimplex() == true) {
#ifdef DEBUG
      std::cout << "[*** The current polytope is a simplex ***]" << std::endl;
#endif
      // Here we have two possibilities: full dimension or low dimension.
      if (topOfStack.dimension() < _dimension) {
        stackOfPolToSimp.top().switchToFullDimension();
      }
      // Now it's a fully dimensional simplex store it on the side and remove it from the stack.
      _allSimplices.push_back(stackOfPolToSimp.top());
#ifdef DEBUG
      std::cout << "After switchToFullDimension():" << std::endl;
      stackOfPolToSimp.top().dump(std::cout, shift);
#endif
      stackOfPolToSimp.pop();
    }
    else {
      // 1. Split the current polytope into sub-polytopes
      // 2. Insert the sub-polytopes into the list
      // 3. Remove the father polytope
#ifdef DEBUG
      std::cout << "[*** The current polytope is to be split ***]" << std::endl;
#endif
      // Work on a copy
      PolytopeToSimplexes topOfStack_cp = stackOfPolToSimp.top();
      const std::vector< unsigned int >& vtx2plit = topOfStack_cp.getListOfVerticesToSplit();
      unsigned int apexNumber = vtx2plit[0];
      std::vector< unsigned int >::const_iterator itEnd = vtx2plit.end();
      std::vector< unsigned int >::const_iterator itBeg = vtx2plit.begin();
      ++itBeg;
      // Insert all elements except the first one, this list should be ordered.
      std::vector< unsigned int > newVtx;
      newVtx.insert(newVtx.begin(), itBeg, itEnd);

      // Now remove the old polytope that has been split.
      stackOfPolToSimp.pop();

      {for (unsigned int facetNumber=0; facetNumber<_numberOfFacets; ++facetNumber) {
        // Exclude all facets belonging to the apex.
#ifdef DEBUG
        std::cout << "Search " << facetNumber << " in : ";
        std::copy(_facetsByVertices[apexNumber].begin(), _facetsByVertices[apexNumber].end(), std::ostream_iterator<unsigned int>(std::cout, " "));
        std::cout << std::endl;
#endif
        if (std::binary_search(_facetsByVertices[apexNumber].begin(), _facetsByVertices[apexNumber].end(), facetNumber) == false) {
          // This is not a facet of the apex so split the base now.
#ifdef DEBUG
          for (unsigned int i=0; i<_dimension-topOfStack_cp.dimension(); ++i) std::cout << "\t";
          std::cout << "Not found facet number : " << facetNumber << std::endl;
#endif
          // Split the cloud of points in baseToSimplices according to their facets.
          // verticesOnGivenFacet = vertices of baseToSimplices which are on the facetNumber-th half-space.
          std::vector< unsigned int > verticesOnGivenFacet;
          std::set_intersection(
              newVtx.begin(),
              newVtx.end(),
              _verticesByFacets[facetNumber].begin(),
              _verticesByFacets[facetNumber].end(),
              std::inserter( verticesOnGivenFacet, verticesOnGivenFacet.begin() ) );
          if (verticesOnGivenFacet.size() >= topOfStack_cp.dimension()) {
            PolytopeToSimplexes newP2S(verticesOnGivenFacet, topOfStack_cp.getListOfProcessedVertices(), apexNumber, topOfStack_cp.dimension()-1);
            stackOfPolToSimp.push(newP2S);
#ifdef DEBUG
            std::cout << "{ ";
            std::copy(newVtx.begin(), newVtx.end(), std::ostream_iterator<unsigned int>(std::cout, " "));
            std::cout << "} INTER { ";
            std::copy(_verticesByFacets[facetNumber].begin(), _verticesByFacets[facetNumber].end(), std::ostream_iterator<unsigned int>(std::cout, " "));
            std::cout << "}" << std::endl;
            for (unsigned int i=0; i<_dimension-newP2S.dimension(); ++i) std::cout << "\t";
            std::cout << "{" << newP2S.dimension() << ":";
            std::copy(newP2S.getListOfProcessedVertices().begin(), newP2S.getListOfProcessedVertices().end(), std::ostream_iterator<unsigned int>(std::cout, " "));
            std::cout << std::endl;
            for (unsigned int i=0; i<_dimension-newP2S.dimension(); ++i) std::cout << "\t";
            std::copy(newP2S.getListOfVerticesToSplit().begin(), newP2S.getListOfVerticesToSplit().end(), std::ostream_iterator<unsigned int>(std::cout, " "));
            std::cout << "}" << std::endl;
#endif
          }
        }
      }}
    }
  }


#ifdef DEBUG
  std::cout << "_allSimplices : " << std::endl;
  for (std::vector< PolytopeToSimplexes >::iterator it=_allSimplices.begin(); it!=_allSimplices.end(); ++it) {
    std::cout << "{" << it->dimension() << ":";
    std::copy(it->getListOfProcessedVertices().begin(), it->getListOfProcessedVertices().end(), std::ostream_iterator<unsigned int>(std::cout, " "));
    std::cout << std::endl;
    std::copy(it->getListOfVerticesToSplit().begin(), it->getListOfVerticesToSplit().end(), std::ostream_iterator<unsigned int>(std::cout, " "));
    std::cout << "}" << std::endl;
    std::cout << std::endl;
  }
#endif
}

double VolumeOfPolytopes_Rn::volume() {
  _volume = 0.;
#ifdef DEBUG
  std::cout << "VolumeOfPolytopes_Rn::volume() : " << _allSimplices.size() << std::endl;
#endif
  std::set< std::set< boost::shared_ptr<Generator_Rn> > > listOfCoordSimplices;
  std::vector< PolytopeToSimplexes >::iterator it1;
  {for (it1=_allSimplices.begin(); it1!=_allSimplices.end(); ++it1) {
    if (it1->dimension() == _dimension) {
      std::set< boost::shared_ptr<Generator_Rn> > setOfGN;
      std::vector< unsigned int >::const_iterator it2;
#ifdef DEBUG
      std::cout << "Current simplex : ";
      std::copy(it1->getListOfProcessedVertices().begin(), it1->getListOfProcessedVertices().end(), std::ostream_iterator<unsigned int>(std::cout, " "));
      std::cout << std::endl;
#endif
      {for (it2=it1->getListOfProcessedVertices().begin(); it2!=it1->getListOfProcessedVertices().end(); ++it2) {
        setOfGN.insert(_polytope->getGenerator( *it2 ));
      }}
      listOfCoordSimplices. insert(setOfGN);
    }
  }}

  std::set< std::set< boost::shared_ptr<Generator_Rn> > >::const_iterator constITER_setGN;
  {for (constITER_setGN=listOfCoordSimplices.begin(); constITER_setGN!=listOfCoordSimplices.end(); ++constITER_setGN) {
    double current_volume = fabs( computeSimplexVolume(*constITER_setGN) );
#ifdef DEBUG
    std::cout << "volumeOfSimplex=" << current_volume << std::endl;
#endif
    _volume += current_volume;
  }}
  return _volume / (double) boost::math::factorial<double>(_dimension);
}

double VolumeOfPolytopes_Rn::computeSimplexVolume(const std::set< boost::shared_ptr<Generator_Rn> >& listOfSimplexVertices) const
{
  unsigned int RnDIM=Rn::getDimension();
  boost::shared_ptr<Generator_Rn> vx0 = *(listOfSimplexVertices.begin());
  boost::numeric::ublas::matrix<double> mat2computeDET(RnDIM, RnDIM);
  std::set< boost::shared_ptr<Generator_Rn> >::const_iterator constITER_GN = listOfSimplexVertices.begin();
  constITER_GN++;
  unsigned int i=0;
  {for ( ; constITER_GN!=listOfSimplexVertices.end(); ++constITER_GN) {
    for (unsigned int j=0; j<mat2computeDET.size2(); ++j)
      mat2computeDET(i,j) = (*constITER_GN)->getCoordinate(j) - vx0->getCoordinate(j);
    i++;
  }}
#ifdef DEBUG
  std::cout << "mat2computeDET=" << mat2computeDET << std::endl;
#endif
  return determinant(mat2computeDET);
}

double VolumeOfPolytopes_Rn::determinant(boost::numeric::ublas::matrix<double> mat) const
{
  double det = 0;
  unsigned int i,j,j1,j2;
  unsigned int n = mat.size1();

  if (n < 1) { // Error
  }
  else if (n == 1) { // Shouldn't get used
    det = mat(0,0);
  }
  else if (n == 2) {
    det = mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1);
  }
  else {
    det = 0;
    for (j1=0; j1<n; j1++) {
      boost::numeric::ublas::matrix<double> smaller_mat(n-1, n-1);
      for (i=1; i<n; i++) {
        j2 = 0;
        for (j=0; j<n; j++) {
          if (j == j1)
            continue;
          smaller_mat(i-1,j2) = mat(i,j);
          j2++;
        }
      }
      det += pow(-1.0, 1.0+j1+1.0) * mat(0,j1) * VolumeOfPolytopes_Rn::determinant(smaller_mat);
    }
  }
  return(det);
}
