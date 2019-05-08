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
/// \file IO_Polytope.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include <cmath>
#include "Rn.h"
#include "IO_Polytope.h"


void IO_Polytope::load(const std::string& filename, boost::shared_ptr<PolyhedralCone_Rn> POLY)
  throw (std::ios_base::failure,std::out_of_range) {

  int dimension, numberOfGenerators, numberOfHalfSpaces;
  std::ifstream file(filename.c_str(), std::ifstream::in);

  if (!file) {
    std::string s("Unable to open ");
    s += filename;
    s += "\n";
    throw std::ios_base::failure(s);
  }

  // Temp array to store constraints.
  std::vector< boost::shared_ptr<HalfSpace_Rn> > HSvect;
  std::string line;
  // Read the comment line.
  std::getline(file,line);
  // Get the 3 integers : space dimension, number of generators and facets.
  std::getline(file,line);
  std::istringstream iline(line);
  iline >> dimension;
  iline >> numberOfHalfSpaces;
  iline >> numberOfGenerators;
  double val;
  ///////////////////
  // Facets block //
  /////////////////
  if (numberOfHalfSpaces != 0) {
    // Read the comment line.
    std::getline(file,line);
    // H : a0 + a1.x1 + ... + an.xn >= 0.
    boost::shared_ptr<HalfSpace_Rn> HS;
    for (int hyp_count=0; hyp_count<numberOfHalfSpaces; hyp_count++) {
      HS.reset(new HalfSpace_Rn(dimension));
      std::getline(file,line);
      std::istringstream iline3(line);
      // The second member b
      iline3 >> val;
      HS->setConstant(val);
      // The normal vector
      for (int coord_count=0; coord_count<dimension; coord_count++) {
        iline3 >> val;
        HS->setCoefficient(coord_count, val);
      }
      POLY->addHalfSpace(HS);
      HSvect.push_back(HS);
    }
  }
  // When we only have the half-spaces, use the truncating algorithm to compute the generators.
  if (numberOfGenerators == 0) {
    file.close();
    return;
  }
  else {
    ///////////////////////
    // Generators block //
    /////////////////////
    // Read the comment line.
    std::getline(file,line);
    // V = (v1, ..., vn)
    boost::shared_ptr<Generator_Rn> VX;
    for (int vtx_count=0; vtx_count<numberOfGenerators; vtx_count++) {
      VX.reset(new Generator_Rn(dimension));
      std::getline(file,line);
      std::istringstream iline2(line);
      VX->load(iline2);
      //for (int coord_count=0; coord_count<dimension; coord_count++) {
      //iline2 >> val;
      //std::cout << "val=" << val << " ";
      //VX->setCoordinate(coord_count, val);
      //}
      POLY->addGenerator(VX);
    }
    //////////////////////////////
    // Facets per vertex block //
    ////////////////////////////
    if (numberOfHalfSpaces != 0 && numberOfGenerators != 0) {
      // Read the comment line.
      std::getline(file,line);
      // {Fi1, Fi2, ... }
      for (int vtx_count=0; vtx_count<numberOfGenerators; vtx_count++) {
        VX = POLY->getGenerator(vtx_count);
        std::string lineOfFacets;
        std::getline(file, lineOfFacets);
        std::istringstream stream(lineOfFacets);
        int facetForGenerator;
        while (stream >> facetForGenerator) {
          if (facetForGenerator < 0 || facetForGenerator >= numberOfHalfSpaces) {
            std::ostringstream stream_;
            stream_ << "Facet number ";
            stream_ << facetForGenerator;
            stream_ << " is outside the range [0, ";
            stream_ << numberOfHalfSpaces;
            stream_ << "] in function IO_Polytope::load() ";
            std::string valString = stream_.str();
            throw std::out_of_range(valString);
          }
          boost::shared_ptr<HalfSpace_Rn> Fi = HSvect[facetForGenerator];
          VX->setFacet(Fi);
        }
      }
    }
  }
  file.close();
}

void IO_Polytope::save(const std::string& filename, boost::shared_ptr<PolyhedralCone_Rn> POLY) throw (std::ios_base::failure) {
  std::ofstream ofs(filename.c_str());

  if (!ofs) {
    std::string s("Unable to open ");
    s += filename;
    s += "\n";
    throw std::ios_base::failure(s);
  }

  // Make sure we print on the stream values consistent with the tolerance in use.
  ofs.precision(15);

  ofs << "# Dimension NumberOfHalfspaces NumberOfGenerators" << std::endl;
  ofs << POLY->dimension() << " ";
  ofs << POLY->numberOfHalfSpaces() << " ";
  ofs << POLY->numberOfGenerators() << std::endl;
  //////////////////
  // Facet block //
  ////////////////
  if (POLY->numberOfHalfSpaces() != 0) {
    ofs << "# HALFSPACES : a0 + a1.x1 + ... + an.xn >= 0." << std::endl;
    // H : a0 + a1.x1 + ... + an.xn >= 0.
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >
      iteHS(POLY->getListOfHalfSpaces());
    for (iteHS.begin(); iteHS.end()!=true; iteHS.next()) {
      boost::shared_ptr<HalfSpace_Rn> HS = iteHS.current();
      // The second member b
      ofs << HS->getConstant() << " ";
      for (unsigned int coord_count=0; coord_count<POLY->dimension(); coord_count++) {
        ofs << HS->getCoefficient(coord_count) << " ";
      }
      ofs << std::endl;
    }
  }
  ofs << "# GENERATORS : V = (v1, ..., vn)" << std::endl;
  //////////////////////
  // Generator block //
  ////////////////////
  // V = (v1, ..., vn)
  boost::shared_ptr<Generator_Rn> VX;
  for (unsigned int vtx_count=0; vtx_count<POLY->numberOfGenerators(); vtx_count++) {
    VX = POLY->getGenerator(vtx_count);
    VX->save(ofs);
    //for (unsigned int coord_count=0; coord_count<Rn::getDimension(); coord_count++) {
    //ofs << VX->getCoordinate(coord_count) << " ";
    //}
    ofs << std::endl;
  }
  //////////////////////////////
  // Facets per vertex block //
  ////////////////////////////
  if (POLY->numberOfHalfSpaces() != 0) {
    ofs << "# GENERATOR CONNECTION : Ha, Hb, ..." << std::endl;
    // {Fi1, Fi2, ... }
    for (unsigned int vtx_count=0; vtx_count<POLY->numberOfGenerators(); vtx_count++) {
      VX = POLY->getGenerator(vtx_count);
      for (unsigned int ngb_count=0; ngb_count<VX->numberOfFacets(); ngb_count++) {
        boost::shared_ptr<HalfSpace_Rn> Fi = VX->getFacet(ngb_count);
        ofs << POLY->getHalfSpaceNumber(Fi) << " ";
      }
      ofs << std::endl;
    }
  }
  ofs.close();
}
