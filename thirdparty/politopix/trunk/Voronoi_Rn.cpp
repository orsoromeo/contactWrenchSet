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
/// \file Voronoi_Rn.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include "GeometricObjectIterator_Rn.h"
#include "PolyhedralAlgorithms_Rn.h"
#include "DoubleDescription_Rn.h"
#include "Voronoi_Rn.h"
#include "IO_Polytope.h"
#include <math.h>
#include <numeric>
#include <sstream>


// Constructor
Voronoi_Rn::Voronoi_Rn(const boost::shared_ptr<Polytope_Rn>& inputSpace, 
		       const std::vector<Point_Rn>& listOfSeeds):_inputSpace(inputSpace),_listOfSeeds(listOfSeeds),_listOfVoronoiCells(0)
{}

bool Voronoi_Rn::compute() throw (std::length_error) {

  if (_listOfSeeds.empty() == true)
    throw std::length_error("Voronoi_Rn::compute() needs seeds as input!");

  if (_listOfSeeds.size() == 1) {
    // Make a deep copy of the current input space
    boost::shared_ptr<Polytope_Rn> is(new Polytope_Rn( *_inputSpace ));
    _listOfVoronoiCells.push_back(is);
    return true;
  }

  // The first cell is the whole input space.
  boost::shared_ptr<Polytope_Rn> is(new Polytope_Rn( *_inputSpace ));
  _listOfVoronoiCells.push_back(is);
  //unsigned int alreadyProcessedSeeds = 1;
  {for (std::vector<Point_Rn>::const_iterator newVoronoiSeed = (_listOfSeeds.begin()+1); newVoronoiSeed!=_listOfSeeds.end(); ++newVoronoiSeed) {

#ifdef DEBUG
    std::cout << std::endl;
    std::cout << "Seed " << newVoronoiSeed-_listOfSeeds.begin() << ": ";
    (newVoronoiSeed)->save(std::cout);
    std::cout << std::endl;
#endif
    ///////////////
    // NEW CELL //
    /////////////
    // Construction of the current cell as a first copy of the input space, then to be chopped.
    boost::shared_ptr<Polytope_Rn> newVoronoiCell(new Polytope_Rn(*_inputSpace));
    unsigned int newTruncationStep = newVoronoiCell->numberOfHalfSpaces();
    {for (std::vector<Point_Rn>::const_iterator processedSeed = _listOfSeeds.begin(); processedSeed!=newVoronoiSeed; ++processedSeed) {
      // Compute all the half-spaces of the current cell associated to newVoronoiSeed.
      boost::shared_ptr<HalfSpace_Rn> HS = computeMidPlane(newVoronoiSeed, processedSeed);
      newVoronoiCell->addHalfSpace(HS);
    }}
    // Now the truncation of the new cell
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > lexmin_ite(newVoronoiCell->getListOfHalfSpaces());
    StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
    DoubleDescription<
      boost::shared_ptr<PolyhedralCone_Rn>,
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
      StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
        DD(newVoronoiCell, lexmin_ite, NRP, newTruncationStep);

    ////////////////
    // OLD CELLS //
    //////////////
    // Chopping of the already created cells.
    std::vector<Point_Rn>::const_iterator processedSeed;
    std::vector< boost::shared_ptr<Polytope_Rn> >::iterator oldVoronoiCell;
    {for (oldVoronoiCell = _listOfVoronoiCells.begin(), processedSeed = _listOfSeeds.begin();
        oldVoronoiCell!=_listOfVoronoiCells.end(); ++oldVoronoiCell, ++processedSeed) {
#ifdef DEBUG
    //gnuplot(std::cout);
    //(newVoronoiSeed)->save(std::cout);
    //std::cout << "### Voronoi cell ###";
    //(*oldVoronoiCell)->dump(std::cout);
    //std::cout << std::endl;
#endif
      // Build the half-spaces the other way
      boost::shared_ptr<HalfSpace_Rn> HS = computeMidPlane(processedSeed, newVoronoiSeed);
      unsigned int oldTruncationStep = (*oldVoronoiCell)->numberOfHalfSpaces();
      (*oldVoronoiCell)->addHalfSpace(HS);
      // Now the truncation of the old cells
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > lexmin_ite2((*oldVoronoiCell)->getListOfHalfSpaces());
      StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP2;
      DoubleDescription<
        boost::shared_ptr<PolyhedralCone_Rn>,
        constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
        StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
          DD2(*oldVoronoiCell, lexmin_ite2, NRP2, oldTruncationStep);
    }}
    _listOfVoronoiCells.push_back(newVoronoiCell);
  }}
#ifdef DEBUG
    //gnuplot(std::cout);
    //(newVoronoiSeed)->save(std::cout);
    //std::cout << "### Voronoi cell ###";
    //(*oldVoronoiCell)->dump(std::cout);
    //std::cout << std::endl;
    //dump(std::cout);
#endif

  return true;
}

boost::shared_ptr<HalfSpace_Rn> Voronoi_Rn::computeMidPlane(std::vector<Point_Rn>::const_iterator iteSeed1, std::vector<Point_Rn>::const_iterator iteSeed2) {
  unsigned int dimension = iteSeed1->dimension();
  boost::shared_ptr<HalfSpace_Rn> HS;
  HS.reset(new HalfSpace_Rn(dimension));
  double sum_square = 0.;
  for (unsigned int coord_count=0; coord_count<dimension; coord_count++) {
    HS->setCoefficient(coord_count, iteSeed1->getCoordinate(coord_count)-iteSeed2->getCoordinate(coord_count));
    double sq = iteSeed1->getCoordinate(coord_count)*iteSeed1->getCoordinate(coord_count);
    sq -= iteSeed2->getCoordinate(coord_count)*iteSeed2->getCoordinate(coord_count);
    sq /= 2.;
    sum_square += sq;
  }
  HS->setConstant(-sum_square);
  return HS;
}

// CHECK POLYHEDRON
bool Voronoi_Rn::checkTopologyAndGeometry() const throw (std::domain_error) {
  return true;
}

/// Dump the cell structure on the given output.
void Voronoi_Rn::dump(std::ostream& out) const {
  out << std::endl << "#VORONOI DIAGRAM" << std::endl;
  unsigned int nb=0;
  std::vector< boost::shared_ptr<Polytope_Rn> >::const_iterator itePC=_listOfVoronoiCells.begin();
  {for (; itePC!=_listOfVoronoiCells.end(); ++itePC) {
    out << "# Cell "<< nb << std::endl;
    (*itePC)->dump(out);
    out << std::endl;
#ifdef DEBUG
    //std::ostringstream stream_C;
    //stream_C << "cell";
    //stream_C << nb;
    //stream_C << ".ptop";
    //std::string cellfile = stream_C.str();
    //IO_Polytope::save(cellfile, *itePC);
#endif
   ++nb;
  }}
}

void Voronoi_Rn::gnuplot(std::ostream& out)  const throw (std::domain_error) {
  if (Rn::getDimension() != 2)
    throw std::domain_error("Voronoi_Rn::gnuplot(std::ostream& out) dimension is not 2 ");
  unsigned int nb = 1;
  // To code an unique color for each polygon.
  std::vector< boost::shared_ptr<Polytope_Rn> >::const_iterator iteVC=_listOfVoronoiCells.begin();
  {for (; iteVC!=_listOfVoronoiCells.end(); ++iteVC) {
    double frac = (double)nb / (double)(*iteVC)->numberOfGenerators();
    std::ostringstream stream_;
    stream_ << nb;
    std::string valString = stream_.str();
    Visualization::gnuplot2D(*iteVC, valString, frac, out);
    ++nb;
  }}
}

