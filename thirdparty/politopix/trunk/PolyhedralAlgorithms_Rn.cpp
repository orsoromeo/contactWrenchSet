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
/// \file PolyhedralAlgorithms.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <numeric>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>
#include "Rn.h"
#include "Polytope_Rn.h"
#include "IO_Polytope.h"
#include "NormalFan_Rn.h"
#include "PolyhedralAlgorithms_Rn.h"

//typedef std::vector< unsigned int > ListOfFaces;


void FaceEnumeration::Compute(const boost::shared_ptr<Polytope_Rn>& thisPolytope) {
  FaceEnumeration FE(thisPolytope);
  FaceEnumeration::ComputeWithFacets(thisPolytope, FE);
  FaceEnumeration::ComputeWithVertices(thisPolytope, FE);
}

void FaceEnumeration::Compute(const boost::shared_ptr<Polytope_Rn>& thisPolytope, FaceEnumeration& FaceEnum) {
  FaceEnumeration::ComputeWithFacets(thisPolytope, FaceEnum);
  FaceEnumeration::ComputeWithVertices(thisPolytope, FaceEnum);
}

void FaceEnumeration::ComputeWithFacets(const boost::shared_ptr<Polytope_Rn>& thisPolytope, FaceEnumeration& FaceEnum) {
  unsigned int RnDIM = Rn::getDimension();
  FaceEnum._allFacesWithFacets.clear();
  std::set< ListOfFaces > allPotentialFaces;
  std::vector< ListOfFaces > thisListOfFacetsPerVertices;
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > constITER_GN1(thisPolytope->getListOfGenerators());
  {for (constITER_GN1.begin(); constITER_GN1.end()!=true; constITER_GN1.next()) {
    unsigned int NbFacets = constITER_GN1.current()->numberOfFacets();
    ListOfFaces thisListOfFaces;
    {for (unsigned int j=0; j<NbFacets; j++) {
      unsigned int NbFj = thisPolytope->getHalfSpaceNumber( constITER_GN1.current()->getFacet(j) );
      // GenPerFacet: number(HalfSpace_Rn) => {Generator_Rn_0, Generator_Rn_1, Generator_Rn_2, ...}
      thisListOfFaces.push_back( NbFj );
    }}
    thisListOfFacetsPerVertices.push_back( thisListOfFaces );
    allPotentialFaces.insert(thisListOfFaces);
  }}
  FaceEnum._allFacesWithFacets.push_back( thisListOfFacetsPerVertices );
  //std::cout << "RnDIM = " << RnDIM << std::endl;
  {for (unsigned int dimensionOfFace=1; dimensionOfFace<=RnDIM; ++dimensionOfFace) {
    //std::cout << std::endl << "D=" << dimensionOfFace << std::endl;
    std::vector< ListOfFaces > kp1_Faces;
    std::set< ListOfFaces > kp1_FacesSet;
    std::vector< ListOfFaces >& k_Faces = FaceEnum._allFacesWithFacets[dimensionOfFace-1];
    std::vector< ListOfFaces >::iterator it1 = k_Faces.begin(), it2 = k_Faces.begin();
    for (it1 = k_Faces.begin(); it1 != k_Faces.end(); ) {
      // To know whether we had to move forward.
      bool it1StepForward = false;
      for (it2 = it1+1; it2 != k_Faces.end(); ) {
        std::vector< unsigned int > interFace;
        std::set_intersection(it1->begin(), it1->end(), it2->begin(), it2->end(), std::inserter(interFace, interFace.end()));
        // Check whether this face as already been processed.
        //std::cout << std::endl << "| it1 = { ";
        //std::copy(it1->begin(), it1->end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
        //std::cout << "} INTER it2 = { ";
        //std::copy(it2->begin(), it2->end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
        //std::cout << "} = { ";
        //std::copy(interFace.begin(), interFace.end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
        //std::cout << "}  ";
        //std::cout << "inserted ";
        if (interFace.empty() == false) {
          if (interFace == *it2) {
            it2 = k_Faces.erase(it2);
            //std::cout << "(it2 erased) ";
          }
          else
            ++it2;
          if (interFace == *it1) {
            it1 = k_Faces.erase(it1);
            //std::cout << "(it1 erased) ";
            it1StepForward = true;
            if (it1 != k_Faces.end())
              it2 = it1+1;
          }
          //std::cout << "in kp1 ";
          // Check it was not already in.
          if (kp1_FacesSet.insert(interFace).second == true)
            kp1_Faces.push_back(interFace);
        }
        else {
          //std::cout << "empty intersection ";
          ++it2;
        }
      } // for (it2 = k_Faces.begin(); it2 != k_Faces.end(); ) {
      // Only move forward iif we didn't do it before.
      if (it1StepForward == false)
        ++it1;
    }
    FaceEnum._allFacesWithFacets.push_back(kp1_Faces);
  }}
  //FaceEnum.printFacesWithFacets(std::cout);
}

void FaceEnumeration::ComputeWithVertices(const boost::shared_ptr<Polytope_Rn>& thisPolytope, FaceEnumeration& FaceEnum) {
  unsigned int nbHS = thisPolytope->numberOfHalfSpaces();
  FaceEnum._allFacesWithVertices.clear();
  FaceEnum._allFacesWithVertices.resize(Rn::getDimension()+1);
  // number(HalfSpace_Rn0) => {Generator_Rn0*, Generator_Rn1*, Generator_Rn2*, ...}
  std::vector< std::vector< unsigned int > > allGenPerHS;
  {for (unsigned int j=0; j<nbHS; ++j) {
    std::vector< unsigned int > genSet;
    allGenPerHS.push_back(genSet);
  }}

  unsigned int generatorNumber=0;
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(thisPolytope->getListOfGenerators());
  {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
    {for (unsigned int ii=0; ii<iteGN.current()->numberOfFacets(); ii++) {
      unsigned int fctNumber = thisPolytope->getHalfSpaceNumber( iteGN.current()->getFacet(ii) );
      allGenPerHS[fctNumber].push_back( iteGN.currentIteratorNumber() );
    }}
    ++generatorNumber;
  }}

  unsigned int dim = 0;
  std::vector< std::vector< ListOfFaces > >::const_iterator iteF;
  for (iteF=FaceEnum._allFacesWithFacets.begin(); iteF!=FaceEnum._allFacesWithFacets.end(); ++iteF) {
    //this_ostream << "F" << dim << ": ";
    std::vector< ListOfFaces >::const_iterator iteFF;
    for (iteFF=iteF->begin(); iteFF!=iteF->end(); ++iteFF) {
      ListOfFaces::const_iterator iteFFF;
      //this_ostream << " { ";
      iteFFF=iteFF->begin();
      std::vector< unsigned int > INTER = allGenPerHS[*iteFFF];
      ++iteFFF;
      for (; iteFFF!=iteFF->end(); ++iteFFF) {
        //this_ostream << *iteFFF << " ";
        std::vector< unsigned int > partial_INTER;
        std::set_intersection(INTER.begin(), INTER.end(),
			      allGenPerHS[*iteFFF].begin(), allGenPerHS[*iteFFF].end(),
			      std::inserter(partial_INTER, partial_INTER.end()));
        INTER = partial_INTER;
      }
      FaceEnum._allFacesWithVertices[dim].push_back(INTER);
      //std::copy(INTER.begin(), INTER.end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
      //this_ostream << "} ";
    }
    //this_ostream << std::endl;
    ++dim;
  }
}

void FaceEnumeration::printFacesWithFacets(std::ostream &this_ostream) const {
  unsigned int dim = 0, allSizes = 0;
  std::vector< std::vector< ListOfFaces > >::const_iterator iteF;
  for (iteF=_allFacesWithFacets.begin(); iteF!=_allFacesWithFacets.end(); ++iteF) {
    allSizes = allSizes + iteF->size();
  }
  // Add the empty polytope and the total polytope.
  allSizes = allSizes + 2;
  this_ostream << "FaceEnumeration::printFacesWithFacets, total number of elements = " << allSizes << std::endl;
  for (iteF=_allFacesWithFacets.begin(); iteF!=_allFacesWithFacets.end(); ++iteF) {
    this_ostream << "F" << dim << ": ";
    std::vector< ListOfFaces >::const_iterator iteFF;
    for (iteFF=iteF->begin(); iteFF!=iteF->end(); ++iteFF) {
      ListOfFaces::const_iterator iteFFF;
      this_ostream << " { ";
      for (iteFFF=iteFF->begin(); iteFFF!=iteFF->end(); ++iteFFF) {
        this_ostream << *iteFFF << " ";
      }
      this_ostream << "} ";
    }
    this_ostream << std::endl;
    ++dim;
  }
}

void FaceEnumeration::printFacesWithVertices(std::ostream &this_ostream) const {
  unsigned int dim = 0, allSizes = 0;
  std::vector< std::vector< ListOfFaces > >::const_iterator iteF;
  for (iteF=_allFacesWithVertices.begin(); iteF!=_allFacesWithVertices.end(); ++iteF) {
    allSizes = allSizes + iteF->size();
  }
  // Add the empty polytope and the total polytope.
  allSizes = allSizes + 2;
  this_ostream << "FaceEnumeration::printFacesWithVertices, total number of elements = " << allSizes << std::endl;
  for (iteF=_allFacesWithVertices.begin(); iteF!=_allFacesWithVertices.end(); ++iteF) {
    this_ostream << "F" << dim << ": ";
    std::vector< ListOfFaces >::const_iterator iteFF;
    for (iteFF=iteF->begin(); iteFF!=iteF->end(); ++iteFF) {
      ListOfFaces::const_iterator iteFFF;
      this_ostream << " { ";
      for (iteFFF=iteFF->begin(); iteFFF!=iteFF->end(); ++iteFFF) {
        this_ostream << *iteFFF << " ";
      }
      this_ostream << "} ";
    }
    this_ostream << std::endl;
    ++dim;
  }
}

void FaceEnumeration::printFacesWithVerticesToSage(std::ostream &this_ostream) const {
  unsigned int dim = 0, allSizes = 0;
  std::vector< std::vector< ListOfFaces > >::const_iterator iteF;
  for (iteF=_allFacesWithVertices.begin(); iteF!=_allFacesWithVertices.end(); ++iteF) {
    allSizes = allSizes + iteF->size();
  }
  // Add the empty polytope and the total polytope.
  allSizes = allSizes + 2;
  this_ostream << "FaceEnumeration::printFacesWithVerticesToSage, total number of elements = " << allSizes << std::endl;
  for (iteF=_allFacesWithVertices.begin(); iteF!=_allFacesWithVertices.end(); ++iteF) {
    this_ostream << "F" << dim << ": " << std::endl;
    std::vector< ListOfFaces >::const_iterator iteFF;
    //std::cout << "ok "; //@
    std::vector< ListOfFaces >::const_iterator stopFF=iteF->end()-1;
    //std::cout << ", stopFF=" << std::endl; //@
    for (iteFF=iteF->begin(); iteFF!=stopFF && iteFF!=iteF->end(); ++iteFF) {
      ListOfFaces::const_iterator iteFFF, stopFFF=iteFF->end()-1;
      this_ostream << "(";
      for (iteFFF=iteFF->begin(); iteFFF!=stopFFF; ++iteFFF) {
        this_ostream << *iteFFF << ", ";
      }
      this_ostream << *iteFFF << "), ";
    }
    //std::cout << "ok2 " << std::endl; //@
    // Print the last one.
    if (iteFF!=iteF->end()) {
      ListOfFaces::const_iterator iteFFF, stopFFF=iteFF->end()-1;
      this_ostream << "(";
      for (iteFFF=iteFF->begin(); iteFFF!=stopFFF; ++iteFFF) {
        this_ostream << *iteFFF << ", ";
      }
      this_ostream << *iteFFF << ") ";
    }
    //std::cout << "ok3 " << std::endl; //@
    this_ostream << std::endl;
    ++dim;
  }
}

void FaceEnumeration::load(const std::string& filename, std::vector< std::vector< ListOfFaces > >& latt) 
  throw (std::ios_base::failure) {
  std::ifstream file(filename.c_str(), std::ifstream::in);
  if (!file) {
    std::string s("Unable to open ");
    s += filename;
    s += "\n";
    throw std::ios_base::failure(s);
  }
  FaceEnumeration::load(file, latt);
  file.close();
}

void FaceEnumeration::save(const std::string& filename, const std::vector< std::vector< ListOfFaces > >& latt) {
  std::ofstream file(filename.c_str());
  if (!file) {
    std::string s("Unable to open ");
    s += filename;
    s += "\n";
    throw std::ios_base::failure(s);
  }
  FaceEnumeration::save(file, latt);
  file.close();
}

void FaceEnumeration::load(std::istream& this_stream, std::vector< std::vector< ListOfFaces > >& latt) 
  throw (std::out_of_range){
  std::string line;
  // Read the comment line.
  std::getline(this_stream, line);
  // Get the 2 integers : space dimension, number of faces.
  std::getline(this_stream, line);
  std::istringstream iline(line);
  unsigned int spaceDimension, totalNumberOfFaces;
  iline >> spaceDimension;
  iline >> totalNumberOfFaces;
  latt.resize(spaceDimension+1);
  unsigned nbReadFaces = 0;
  while (nbReadFaces != totalNumberOfFaces) {
    ListOfFaces thisOne;
    std::getline(this_stream, line);
    std::istringstream iline3(line);
    unsigned int dimensionOfFace, numberOfVtx, val;
    iline3 >> dimensionOfFace;
    iline3 >> numberOfVtx;
    if (dimensionOfFace > spaceDimension) {
      std::string errorMsg("FaceEnumeration::load() wrong face dimension");
      throw std::out_of_range(errorMsg);
    }
    // dimensionOfFace numberOfVtx a0 a1 ... an
    for (unsigned int count=0; count<numberOfVtx; ++count) {
      iline3 >> val;
      thisOne.push_back(val);
    }
    latt[dimensionOfFace].push_back(thisOne);
    ++nbReadFaces;
  }
}

void FaceEnumeration::save(std::ostream& this_stream, const std::vector< std::vector< ListOfFaces > >& latt) {
  unsigned int allSizes = 0;
  std::vector< std::vector< ListOfFaces > >::const_iterator iteF;
  for (iteF=latt.begin(); iteF!=latt.end(); ++iteF) {
    allSizes = allSizes + iteF->size();
  }
  // Don't add the empty polytope and the total polytope.
  //allSizes = allSizes + 2;
  //this_ostream << "FaceEnumeration::printFacesWithFacets, total number of elements = " << allSizes << std::endl;
  unsigned int spaceDimension = Rn::getDimension();
  unsigned int totalNumberOfFaces = allSizes;
  this_stream << "# Dimension LatticeSize" << std::endl;
  this_stream << spaceDimension << " ";
  this_stream << totalNumberOfFaces << std::endl;
  unsigned dimReadFaces = 0;
  for (iteF=latt.begin(); iteF!=latt.end(); ++iteF) {
    std::vector< ListOfFaces >::const_iterator iteFF;
    for (iteFF=iteF->begin(); iteFF!=iteF->end(); ++iteFF) {
      this_stream << dimReadFaces << " " << iteFF->size() << " ";
      ListOfFaces::const_iterator iteFFF;
      for (iteFFF=iteFF->begin(); iteFFF!=iteFF->end(); ++iteFFF)
        this_stream << *iteFFF << " ";
      this_stream << std::endl;
    }
    ++dimReadFaces;
  }
}

boost::shared_ptr<Polytope_Rn> MinkowskiSum::compute() throw (std::domain_error) {
  unsigned int RnDIM=Rn::getDimension();

  if (_firstOperand->numberOfGenerators()==0 ||  _firstOperand->numberOfHalfSpaces()==0 ||
     _secondOperand->numberOfGenerators()==0 || _secondOperand->numberOfHalfSpaces()==0)
    throw std::domain_error("MinkowskiSum::compute() needs two double-description polytopes i.e. with both vertices and half-spaces.");

  // minkowskiVertices(i,j)==1 <=> a_i + b_j is a Minkowski vertex.
  //  MinkowskiDecomposition[k] = (a_i,b_j) <=> c_k = a_i + b_j is a Minkowski vertex.
  //   _neighboursA(a_i,a_j)==1 <=> a_i R a_j
  //   _neighboursB(b_i,b_j)==1 <=> b_u R b_v
  _neighboursA.resize(_firstOperand->numberOfGenerators());
  _neighboursB.resize(_secondOperand->numberOfGenerators());
  _firstOperand->fillNeighbourMatrix(_neighboursA, RnDIM-1);
  _secondOperand->fillNeighbourMatrix(_neighboursB, RnDIM-1);
  //std::cout << "_neighboursA: " << std::endl;
  //{for (unsigned int i=0; i<_neighboursA.size(); ++i) {
    //std::cout << i << ": ";
    //std::copy(_neighboursA[i].begin(), _neighboursA[i].end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
    //std::cout << std::endl;
  //}}
  //std::cout << "_neighboursB: " << std::endl;
  //{for (unsigned int i=0; i<_neighboursB.size(); ++i) {
    //std::cout << i << ": ";
    //std::copy(_neighboursB[i].begin(), _neighboursB[i].end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
    //std::cout << std::endl;
  //}}
  // Give the relation between a vertex of A and all of its associated Minkowski vertices i.e. its polyhedrical cap
  // _A2C[i] = [ c_u, c_v, ... ]  }
  //                             } => c_u = a_i + b_j
  // _B2C[j] = [ c_u, c_w, ... ]  }
  _A2C.resize(_firstOperand->numberOfGenerators());
  _B2C.resize(_secondOperand->numberOfGenerators());

  listOfGeometricObjects< boost::shared_ptr<PolyhedralCone_Rn> > listOfDualCones_A;
  listOfGeometricObjects< boost::shared_ptr<PolyhedralCone_Rn> > listOfDualCones_B;
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_A(_firstOperand->getListOfGenerators());
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN_B(_secondOperand->getListOfGenerators());
  //std::cout << "_neighboursA_:" << std::endl;
  {for (iteGN_A.begin(); iteGN_A.end()!=true; iteGN_A.next()) {
    const boost::shared_ptr<PolyhedralCone_Rn>& cone_a = _firstOperand->getPrimalCone( iteGN_A.current() );
    listOfDualCones_A.push_back( cone_a->computeDualPolyhedralCone() );
  }}
  //std::cout << "_neighboursB_:" << std::endl;
  {for (iteGN_B.begin(); iteGN_B.end()!=true; iteGN_B.next()) {
    const boost::shared_ptr<PolyhedralCone_Rn>& cone_b = _secondOperand->getPrimalCone( iteGN_B.current() );
    listOfDualCones_B.push_back( cone_b->computeDualPolyhedralCone() );
  }}
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > LVX_A(_firstOperand->getListOfGenerators());
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > LVX_B(_secondOperand->getListOfGenerators());
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<PolyhedralCone_Rn> > LPC_A(listOfDualCones_A);
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<PolyhedralCone_Rn> > LPC_B(listOfDualCones_B);
  unsigned int minkowskiVertexNumber=0;
  {for (LPC_B.begin(), LVX_B.begin(); LPC_B.end()!=true; LPC_B.next(), LVX_B.next()) {
    boost::shared_ptr<PolyhedralCone_Rn> Cb = LPC_B.current();
    {for (LPC_A.begin(), LVX_A.begin(); LPC_A.end()!=true; LPC_A.next(), LVX_A.next()) {
      boost::shared_ptr<PolyhedralCone_Rn> Ca = LPC_A.current();
#ifdef DEBUG
      //std::cout << "=== Cone from A number " << LPC_A.currentIteratorNumber() << std::endl;
      //std::cout << "=== Cone from B number " << LPC_B.currentIteratorNumber() << std::endl;
      //std::ostringstream stream_A_;
      //stream_A_ << "EX/ConeA";
      //stream_A_ << LPC_A.currentIteratorNumber();
      //stream_A_ << ".pcon";
      //std::string coneAfile_ = stream_A_.str();
      // Primal cones.
      //IO_Polytope::save(coneAfile_, LPC_A.current());
      //LPC_A.current()->checkTopologyAndGeometry();
      //std::ostringstream stream_B_;
      //stream_B_ << "EX/ConeB";
      //stream_B_ << LPC_B.currentIteratorNumber();
      //stream_B_ << ".pcon";
      //std::string coneBfile_ = stream_B_.str();
      //IO_Polytope::save(coneBfile_, LPC_B.current());
      //LPC_B.current()->checkTopologyAndGeometry();
#endif
      boost::shared_ptr<PolyhedralCone_Rn> intersectionCone;
      intersectionCone.reset(new PolyhedralCone_Rn());
#ifdef DEBUG
      //unsigned int i=LPC_A.currentIteratorNumber();
      //unsigned int j=LPC_B.currentIteratorNumber();
      //std::cout << "=== Cone from A number " << i << std::endl;
      //std::cout << "=== Cone from B number " << j << std::endl;
      //std::ostringstream stream_A;
      //stream_A << "EX/ConeA";
      //stream_A << i;
      //stream_A << ".pcon";
      //std::string coneAfile = stream_A.str();
      //IO_Polytope::save(coneAfile, Ca);
      //Ca->checkTopologyAndGeometry();
      //std::ostringstream stream_B;
      //stream_B << "EX/ConeB";
      //stream_B << j;
      //stream_B << ".pcon";
      //std::string coneBfile = stream_B.str();
      //IO_Polytope::save(coneBfile, Cb);
      //Cb->checkTopologyAndGeometry();
#endif
      // Fill the data structures.
      int truncationStep = Ca->numberOfHalfSpaces();
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHSA(Ca->getListOfHalfSpaces());
      for (iteHSA.begin(); iteHSA.end()!=true; iteHSA.next())
        intersectionCone->addHalfSpace(iteHSA.current());
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGNA(Ca->getListOfGenerators());
      for (iteGNA.begin(); iteGNA.end()!=true; iteGNA.next()) {
        // Make a deep copy of each generator
        boost::shared_ptr<Generator_Rn> gn(new Generator_Rn( *(iteGNA.current()) ));
        intersectionCone->addGenerator(gn);
      }
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHSB(Cb->getListOfHalfSpaces());
      for (iteHSB.begin(); iteHSB.end()!=true; iteHSB.next())
        intersectionCone->addHalfSpace(iteHSB.current());
      // Compute the intersection between the two polyhedral cones.
      //bool notEmpty = intersectionCone->truncate(truncationStep);
      //std::cout << iteHS.current()->getConstant() << "\t";
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >
        lexmin_ite(intersectionCone->getListOfHalfSpaces());
      //NoRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
      StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
      DoubleDescription<
        boost::shared_ptr<PolyhedralCone_Rn>,
        constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
        //NoRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
        StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
          DD(intersectionCone, lexmin_ite, NRP, truncationStep);
      bool notEmpty = !DD.getIsEmpty();
#ifdef DEBUG
      std::cout << "@@@ nbEDG=" << intersectionCone->numberOfGenerators() << "( ";
#endif
      // The second test is necessary when dealing with parallelism.
      if (notEmpty == true && intersectionCone->numberOfGenerators() >= RnDIM) {
        _NF_Cones.push_back(intersectionCone);
        boost::shared_ptr<Generator_Rn> VX(new Generator_Rn(RnDIM));
        VX->makeSum(LVX_A.current(), LVX_B.current());
#ifdef DEBUG
        for (unsigned int k=0; k<RnDIM; k++) {
          std::cout << LVX_A.current()->getCoordinate(k) << "+" << LVX_B.current()->getCoordinate(k)
              << "=" << VX->getCoordinate(k) << "\t";
        }
        //unsigned int u=LPC_A.currentIteratorNumber();
        //unsigned int v=LPC_B.currentIteratorNumber();
        //std::ostringstream stream_C;
        //stream_C << "EX_2/ConeC";
        //stream_C << minkowskiVertexNumber;
        //stream_C << "--";
        //stream_C << u;
        //stream_C << "_";
        //stream_C << v;
        //stream_C << ".pcon";
        //std::string coneCfile = stream_C.str();
        //intersectionCone->checkTopologyAndGeometry();
        //IO_Polytope::save(coneCfile, intersectionCone);
#endif
        _NF_Vertices.push_back(VX);
        // Now store the bijection (a_i, b_j) <-> c_k
        _MinkowskiDecomposition.push_back( std::make_pair(LPC_A.currentIteratorNumber(), LPC_B.currentIteratorNumber()) );
        _MinkowskiDecompositionOK.push_back(true);
        // Now write the connection between a_i and c_k and also between b_j and c_k.
        _A2C[LPC_A.currentIteratorNumber()].push_back(minkowskiVertexNumber);
        _B2C[LPC_B.currentIteratorNumber()].push_back(minkowskiVertexNumber);
        ++minkowskiVertexNumber;
      }
#ifdef DEBUG
      std::cout << ")" << std::endl;
#endif
    }}
  }}

  processNormalFan2();

  //_sum->checkTopologyAndGeometry();
  return _sum;
}

void MinkowskiSum::processNormalFan2()
{
  unsigned int RnDIM=Rn::getDimension();
  double TOL=Rn::getTolerance();

  boost::numeric::ublas::matrix<double> matOfVtx(_NF_Vertices.size(), RnDIM);
  {
    unsigned int vtxNumber=0;
    std::vector< boost::shared_ptr<Generator_Rn> >::iterator iteVX2;
    {for (iteVX2 = _NF_Vertices.begin(); iteVX2!=_NF_Vertices.end(); ++iteVX2) {
      boost::numeric::ublas::matrix_row< boost::numeric::ublas::matrix<double> > matRow(matOfVtx, vtxNumber); 
      std::copy( (*iteVX2)->begin(), (*iteVX2)->end(), matRow.begin());
      ++vtxNumber;
    }}
  }
  std::vector< std::vector< unsigned int > > VerticesPerFacet;
  // To count the total number of edges
  std::vector< boost::shared_ptr<HalfSpace_Rn> > setOfAllHalfSpaces;
  std::vector< boost::shared_ptr<Generator_Rn> >::iterator iteVX = _NF_Vertices.begin();
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> >::iterator itePC = _NF_Cones.begin();
  {for (; itePC!=_NF_Cones.end(); ++itePC, ++iteVX) {
    // The big job begins
    {for (unsigned int u=0; u<(*itePC)->numberOfGenerators(); ++u) {
      // Build the corresponding half-space.
      boost::shared_ptr<HalfSpace_Rn> HS;
      HS.reset(new HalfSpace_Rn(RnDIM));
      double sum=0.;
      // The normal vector
      for (unsigned int coord_count=0; coord_count<RnDIM; coord_count++) {
        double this_coord = (*itePC)->getGenerator(u)->getCoordinate(coord_count);
        //HYP->setCoefficient(coord_count, iteEDG.current()->getCoordinate(coord_count));
        HS->setCoefficient(coord_count, this_coord);
        sum = sum + this_coord*(*iteVX)->getCoordinate(coord_count);
      }
      HS->setConstant(-sum);
      // Prepare the array of vertices IN/ON.
      std::vector< unsigned int > vtxStates;
      double halfSpaceNorm = std::inner_product(HS->begin(), HS->end(), HS->begin(), 0.);
      halfSpaceNorm = sqrt(halfSpaceNorm);
      // Now iterate on the list of vertices to classify them.
      boost::numeric::ublas::vector<double> dist2Hyp = prod(matOfVtx, HS->vect());
      // Possible future optimization in high dimensions
      //boost::numeric::ublas::vector<double> dist2Hyp(NF_Vertices.size());
      //axpy_prod(matOfVtx, HS->vect(), dist2Hyp, true);
      boost::numeric::ublas::scalar_vector<double> scalsum(_NF_Vertices.size(),sum);
      dist2Hyp = dist2Hyp - scalsum;
      dist2Hyp = dist2Hyp / halfSpaceNorm;;
      {for (unsigned int vtxNumber2=0; vtxNumber2<dist2Hyp.size(); ++vtxNumber2) {
        if (dist2Hyp(vtxNumber2) > TOL) {
          // Nothing to do !
        }
        //else if (distanceToHyperplane < -TOL) {
          //std::cerr.precision(15);
          //std::cerr << "distanceToHyperplane=" << distanceToHyperplane;
          //std::cerr << std::endl;
          //HS->dump(std::cerr);
          //std::cerr << std::endl;
          //(*iteVX2)->dump(std::cerr);
          //std::cerr << std::endl;
          //throw std::domain_error("Polytope_Rn::processNormalFan2() current vertex is out !");
        //}
        else {
          vtxStates.push_back(vtxNumber2);
        }
      }}
      // Now all the states are set for this current half-space.

      bool inserted=false;
      {for (unsigned int i=0; i<VerticesPerFacet.size() && inserted==false; ++i) {
        if (VerticesPerFacet[i].size() == vtxStates.size()) {
          if (VerticesPerFacet[i] == vtxStates)
            inserted = true;
        }
        else if (VerticesPerFacet[i].size() > vtxStates.size()) {
          if (std::includes(VerticesPerFacet[i].begin(), VerticesPerFacet[i].end(), vtxStates.begin(), vtxStates.end()) == true) {
            // No need to do anything as the current set of
            // generators is included in the set number j.
            //std::cout << "false1 (size=" << i << ")" << std::endl;
            inserted = true;
          }
        }
        else {
          if (std::includes(vtxStates.begin(), vtxStates.end(), VerticesPerFacet[i].begin(), VerticesPerFacet[i].end()) == true) {
            VerticesPerFacet[i] = vtxStates;
            setOfAllHalfSpaces[i] = HS;
            inserted = true;
            // Do not return yet as we could need to erase another cell in the array.
          }
        }
      }}
      if (inserted == false) {
        VerticesPerFacet.push_back(vtxStates);
        setOfAllHalfSpaces.push_back(HS);
        //std::cout << "inserted at the end (" << VerticesPerFacet.size() << ")" << std::endl;
      }

    }} // {for (unsigned int u=0; u<(*itePC)->numberOfGenerators(); ++u) {
  }} // {for (; itePC!=NF_Cones.end(); ++itePC, ++iteVX) {

  // Last job: reconstruct the polytope
  unsigned int vtxNb=0;
  iteVX=_NF_Vertices.begin();
  while (iteVX != _NF_Vertices.end()) {
    // We insert a vertex only if it has the correct number of half-spaces.
    // This is due to the fact that somes dual cones have the right number of edges
    // but in fact are "flat" numerically speaking so in the reduction process they
    // "lose" edges and furthermore do not have at least RnDIM edges.
    bool foundEnoughHS=false;
    unsigned int nbHS=0;
    {for (unsigned int i=0; i<VerticesPerFacet.size() && foundEnoughHS==false; ++i) {
      {for (unsigned int j=0; j<VerticesPerFacet[i].size(); ++j) {
        if (VerticesPerFacet[i][j] == vtxNb)
          ++nbHS;
        if (nbHS >= RnDIM)
          foundEnoughHS = true;
      }}
    }}
    if (foundEnoughHS == true) {
      _sum->addGenerator( (*iteVX) );
    }
    else {
      // A very important step if we want to use the array _MinkowskiDecomposition: some dual cones are flat
      // so they lose edges during the process and have foundEnoughHS to false. So we need to mark them in
      // the corresponding array.
      _MinkowskiDecompositionOK[vtxNb] = false;
    }
    ++vtxNb;
    ++iteVX;
  }
  {for (unsigned int i=0; i<VerticesPerFacet.size(); ++i) {
    {for (unsigned int j=0; j<VerticesPerFacet[i].size(); ++j) {
      // Give the vertex its corresponding facets.
      _NF_Vertices[ VerticesPerFacet[i][j] ]->setFacet( setOfAllHalfSpaces[i] );
    }}
    _sum->addHalfSpace( setOfAllHalfSpaces[i], false);
  }}
}

void MinkowskiSum::processNormalFan1()
{
  unsigned int RnDIM=Rn::getDimension();
  double TOL=Rn::getTolerance();
  double TOL2=TOL*TOL;
  // Give the relation between a vertex of A and all of its associated Minkowski vertices i.e. its polyhedrical cap
  // _A2C[i] = [ c_u, c_v, ... ]  }
  //                             } => c_u = a_i + b_j
  // _B2C[j] = [ c_u, c_w, ... ]  }
  // minkowskiVertices(i,j)==1 <=> a_i + b_j is a Minkowski vertex.
  //  MinkowskiDecomposition[k] = (a_i,b_j) <=> c_k = a_i + b_j is a Minkowski vertex.
  //   _neighboursA(a_i,a_j)==1 <=> a_i R a_j
  //   _neighboursB(b_i,b_j)==1 <=> b_u R b_v
  unsigned int currentNumber=0;
  std::vector< boost::shared_ptr<HalfSpace_Rn> > setOfAllHalfSpaces;
  std::vector< boost::shared_ptr<Generator_Rn> >::iterator iteVX = _NF_Vertices.begin();
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> >::iterator itePC = _NF_Cones.begin();
  // Map between an edge and its norm
  std::map< boost::shared_ptr<Generator_Rn>, double > allEdgesNorms;
  {for (; itePC!=_NF_Cones.end(); ++itePC, ++iteVX) {
    // Build a half-space from the dual polyhedral cone generator.
    // So iterate on the list of its edges.

    // First of all get the ancestors of the current Minkowski vertex
    unsigned int a_i = _MinkowskiDecomposition[ currentNumber ].first;
    unsigned int b_j = _MinkowskiDecomposition[ currentNumber ].second;
    // Now get all of their facet neighbours
    // Now inspect all the possible combinations to know whether they provide a Minkowski vertex.
    std::vector< unsigned int > c_kMinkowskiVertices;
    std::vector< unsigned int >::const_iterator iteNgb;
    // Now get all a_i neighbours and find their contributions in terms of Minkowski vertices
    {for (unsigned int j=0; j<_neighboursA[a_i].size(); ++j) {
      unsigned int a_i_ngb = _neighboursA[a_i][j];
      {for (unsigned int k=0; k<_A2C[a_i_ngb].size(); ++k) {
        // Make sure we will not process the current normal fan with itself.
        if (_A2C[a_i_ngb][k] != currentNumber)
          c_kMinkowskiVertices.push_back( _A2C[a_i_ngb][k] );
      }}
    }}
    {for (unsigned int i=0; i<_neighboursB[b_j].size(); ++i) {
      unsigned int b_j_ngb=_neighboursB[b_j][i];
      {for (unsigned int k=0; k<_B2C[b_j_ngb].size(); ++k) {
        if (_B2C[b_j_ngb][k] != currentNumber)
          c_kMinkowskiVertices.push_back( _B2C[b_j_ngb][k] );
      }}
    }}
    //std::copy(c_kMinkowskiVertices.begin(), c_kMinkowskiVertices.end(), std::ostream_iterator<unsigned int>(std::cout, " ") );

    // The big job begins
    {for (unsigned int u=0; u<(*itePC)->numberOfGenerators(); ++u) {
      // Build the corresponding half-space.
      boost::shared_ptr<HalfSpace_Rn> HS;
      HS.reset(new HalfSpace_Rn(RnDIM));
      double sum=0.;
      // The normal vector
      for (unsigned int coord_count=0; coord_count<RnDIM; coord_count++) {
        double this_coord = (*itePC)->getGenerator(u)->getCoordinate(coord_count);
        //HYP->setCoefficient(coord_count, iteEDG.current()->getCoordinate(coord_count));
        HS->setCoefficient(coord_count, this_coord);
        sum = sum + this_coord*(*iteVX)->getCoordinate(coord_count);
      }
      HS->setConstant(-sum);
      const boost::shared_ptr<Generator_Rn>& gn1=(*itePC)->getGenerator(u);
      double gn1_sum2 = std::inner_product(gn1->begin(), gn1->end(), gn1->begin(), 0.);
      {for (iteNgb=c_kMinkowskiVertices.begin(); iteNgb!=c_kMinkowskiVertices.end(); ++iteNgb) {
        // Check all generators of the neighbour cone against the current one.
        bool check=true;
        {for (unsigned int v=0; v<_NF_Cones[*iteNgb]->numberOfGenerators() && check==true; ++v) {
          const boost::shared_ptr<Generator_Rn>& gn2=_NF_Cones[*iteNgb]->getGenerator(v);
          double scal_prod = std::inner_product(gn1->begin(), gn1->end(), gn2->begin(), 0.);
          if (scal_prod > 0.) {
            double coef1 = scal_prod / gn1_sum2;
            // Check whether the current edge is in the "tube" of the second one, if not switch them.
            if (gn1->getNormalDistance(gn2, coef1, RnDIM) < TOL2) {
              // We found a copy of the current edge so remove it from the explored polyhedrical cone.
              _NF_Cones[*iteNgb]->removeGenerator(v);
              _NF_Vertices[*iteNgb]->setFacet(HS);
              check = false;
            }
            else if (scal_prod > 0.) {
              double gn2_sum2 = 0.;
              // First try to get the norm if it was previously computed
              std::map< boost::shared_ptr<Generator_Rn>, double >::iterator iteNorm = allEdgesNorms.find( gn2 );
              if (iteNorm == allEdgesNorms.end())
                gn2_sum2 = std::inner_product(gn2->begin(), gn2->end(), gn2->begin(), 0.);
              else
                gn2_sum2 = iteNorm->second;
              double coef2 = scal_prod / gn2_sum2;
              if (gn2->getNormalDistance(gn1, coef2, RnDIM) < TOL2) {
                _NF_Cones[*iteNgb]->removeGenerator(v);
                _NF_Vertices[*iteNgb]->setFacet(HS);
                check = false;
                if (iteNorm != allEdgesNorms.end())
                  allEdgesNorms.erase(iteNorm);
              }
            }
          }
        }}
      }}
      // Insert the half-space in the polytope list, if not already in,
      // and insert it into the vertex facet list.
      (*iteVX)->setFacet(HS);
      setOfAllHalfSpaces.push_back( HS );
      // Now remove the current generator as it has become useless.
      (*itePC)->removeGenerator(u);
    }} // {for (unsigned int u=0; u<(*itePC)->numberOfGenerators() ; ++u) {

    // Insert the vertex into the polytope vertex list.
    _sum->addGenerator( (*iteVX) );
    ++currentNumber;
  }} // {for (; itePC!=NF_Cones.end(); ++itePC, ++iteVX) {

  // Check everything's fine
  {for (itePC=_NF_Cones.begin(); itePC!=_NF_Cones.end(); ++itePC) {
    if ((*itePC)->numberOfGenerators() != 0)
      throw std::domain_error("MinkowskiSum::processNormalFan1() dual cones reduction not operational !");
  }}

  // Transfer all half-spaces into the main list.
  std::vector< boost::shared_ptr<HalfSpace_Rn> >::const_iterator iteHS=setOfAllHalfSpaces.begin();
  {for (iteHS=setOfAllHalfSpaces.begin(); iteHS!=setOfAllHalfSpaces.end(); ++iteHS) {
    // Thanks to the previous work we do not need to check.
    _sum->addHalfSpace( *iteHS, false);
  }}

}

void MinkowskiSum::processNormalFan0()
{
  unsigned int RnDIM=Rn::getDimension();

  std::vector< boost::shared_ptr<Generator_Rn> >::iterator iteVX=_NF_Vertices.begin();
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> >::iterator itePC=_NF_Cones.begin();
  for (; itePC!=_NF_Cones.end(); ++itePC, ++iteVX) {
    // Build a half-space from the dual polyhedral cone generator.
    // So iterate on the list of its edges.

    //constIteratorOfListOfGenerators iteEDG(itePC.current()->getListOfGenerators());
    //for (iteEDG.begin(); iteEDG.end()!=true; iteEDG.next()) {
    for (unsigned int u=0; u<(*itePC)->numberOfGenerators(); u++) {
      boost::shared_ptr<HalfSpace_Rn> HS;
      HS.reset(new HalfSpace_Rn(RnDIM));
      double sum=0.;
      // The normal vector
      for (unsigned int coord_count=0; coord_count<RnDIM; coord_count++) {
        double this_coord = (*itePC)->getGenerator(u)->getCoordinate(coord_count);
        //HYP->setCoefficient(coord_count, iteEDG.current()->getCoordinate(coord_count));
        HS->setCoefficient(coord_count, this_coord);
        sum = sum + this_coord*(*iteVX)->getCoordinate(coord_count);
      }
      HS->setConstant(-sum);
      // Insert the half-space in the polytope list, if not already in,
      // and insert it into the vertex facet list.
      (*iteVX)->setFacet(_sum->addHalfSpace(HS, true));
    }
    // Insert the vertex into the polytope vertex list.
    _sum->addGenerator( (*iteVX) );
  }

}

//========================================================================================

/// For each facet of the sum F_C, build F_A and F_B such as F_C = F_A + F_B
void PseudoSumWithoutCaps::computeCapHalfSpaces(
    const std::set< unsigned int >& firstOperandCaps,
    const std::set< unsigned int >& secondOperandCaps,
    std::set< unsigned int >& sumCaps) throw (std::domain_error)
{
#ifdef DEBUG
  std::cout << "PseudoSumWithoutCaps::computeCapHalfSpaces()" << endl;
#endif
  if (_sum->numberOfGenerators()==0)
    throw std::domain_error("PseudoSumWithoutCaps::computeCapHalfSpaces() We need to compute the Minkowski sum first.");

  // For each facet of the sum, store its vertices.
  unsigned int facetNumber=0;
  std::vector< std::vector<unsigned int> > listOfVerticesPerFacet(_sum->numberOfHalfSpaces());
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS(_sum->getListOfHalfSpaces());
  {for (iteHS.begin(); iteHS.end()!=true; iteHS.next()) {
#ifdef DEBUG
    std::cout << "facetNumber=" << facetNumber << ":";
#endif
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(_sum->getListOfGenerators());
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
      for (unsigned int i=0; i<iteGN.current()->numberOfFacets(); ++i) {
        if (iteGN.current()->getFacet(i) == iteHS.current()) {
          listOfVerticesPerFacet[facetNumber].push_back(iteGN.currentIteratorNumber());
#ifdef DEBUG
          std::cout << " " << iteGN.currentIteratorNumber();
#endif
        }
      }
    }}
    ++facetNumber;
#ifdef DEBUG
    std::cout << endl;
#endif
  }}

  // Reorder the data structures in the father class.
  std::vector< std::pair< unsigned int, unsigned int > > thisMinkDecomposition;
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> > theseNF_Cones;
  std::vector< boost::shared_ptr<Generator_Rn> > theseNF_Vertices;
  {for (unsigned int i=0; i<_MinkowskiDecompositionOK.size(); ++i) {
    if (_MinkowskiDecompositionOK[i] == true) {
      thisMinkDecomposition.push_back( _MinkowskiDecomposition[i] );
      theseNF_Vertices.push_back( _NF_Vertices[i] );
      theseNF_Cones.push_back( _NF_Cones[i] );
    }
  }}
  _MinkowskiDecomposition = thisMinkDecomposition;
  _NF_Vertices = theseNF_Vertices;
  _NF_Cones = theseNF_Cones;

  // Now facet by facet, decompose the Minkowski vertices into vertices of A and B.
  unsigned int sumFacetNb = 0;
  std::vector< std::vector<unsigned int> >::const_iterator itVtx = listOfVerticesPerFacet.begin();
  {for (; itVtx!=listOfVerticesPerFacet.end(); ++itVtx) {
    std::set< unsigned int > listOfVtxA, listOfVtxB;
    std::vector<unsigned int>::const_iterator MinkVtx = itVtx->begin();
#ifdef DEBUG
    std::cout << endl << "Compute F_A and F_B for facet " << sumFacetNb << endl;
#endif
    {for (; MinkVtx!=itVtx->end(); ++MinkVtx) {
      unsigned int MinkNumber = *MinkVtx;
      listOfVtxA.insert( _MinkowskiDecomposition[MinkNumber].first );
      listOfVtxB.insert( _MinkowskiDecomposition[MinkNumber].second);
    }}

    std::set<unsigned int> listOfFacetsOfVtxA;
    std::set<unsigned int> tmpInterResForA, tmp2InterResForA;
    std::set<unsigned int>::const_iterator itVtxNbOperand;
    itVtxNbOperand = listOfVtxA.begin();
    for (unsigned int ngb_count=0; ngb_count<_firstOperand->getGenerator(*itVtxNbOperand)->numberOfFacets(); ngb_count++) {
      boost::shared_ptr<HalfSpace_Rn> Fi = _firstOperand->getGenerator(*itVtxNbOperand)->getFacet(ngb_count);
      listOfFacetsOfVtxA.insert(_firstOperand->getHalfSpaceNumber(Fi));
    }
    //_firstOperand->getGenerator(*itVtxNbOperand)->exportFacets(listOfFacetsOfVtxA);
#ifdef DEBUG
    std::cout << "F_A:";
    std::cout << " V = {";
    std::cout << " " << *itVtxNbOperand;
#endif
    tmpInterResForA = listOfFacetsOfVtxA;
    itVtxNbOperand++;
    {for (; itVtxNbOperand!=listOfVtxA.end(); ++itVtxNbOperand) {
#ifdef DEBUG
      std::cout << " " << *itVtxNbOperand;
#endif
      std::set<unsigned int> curListOfFacetsOfVtxA;
      for (unsigned int ngb_count=0; ngb_count<_firstOperand->getGenerator(*itVtxNbOperand)->numberOfFacets(); ngb_count++) {
        boost::shared_ptr<HalfSpace_Rn> Fi = _firstOperand->getGenerator(*itVtxNbOperand)->getFacet(ngb_count);
        curListOfFacetsOfVtxA.insert(_firstOperand->getHalfSpaceNumber(Fi));
      }
      tmp2InterResForA = tmpInterResForA;
      tmpInterResForA.clear();
      std::set_intersection(
          tmp2InterResForA.begin(),
          tmp2InterResForA.end(),
          curListOfFacetsOfVtxA.begin(),
          curListOfFacetsOfVtxA.end(),
          std::inserter(tmpInterResForA, tmpInterResForA.end()));
      if (listOfFacetsOfVtxA.empty() == true)
        throw std::domain_error("PseudoSumWithoutCaps::computeCapHalfSpaces() the i-Face of polytope A has 0 half-space.");
    }}
#ifdef DEBUG
    std::cout << " }";
#endif
    // At this step the half-spaces support of F_A are computed.
    std::set<unsigned int> InterResForA;
    std::set_intersection(
        tmpInterResForA.begin(),
        tmpInterResForA.end(),
        firstOperandCaps.begin(),
        firstOperandCaps.end(),
        std::inserter(InterResForA, InterResForA.end()));
#ifdef DEBUG
    std::cout << " ; H = { ";
    std::copy(InterResForA.begin(), InterResForA.end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
    //std::cout << listOfFacetsOfVtxA.size();
    std::cout << "}";
    std::cout << endl;
#endif

    if (InterResForA.empty() != true)
      sumCaps.insert( sumFacetNb );
    else {
      ///////////////
      // 2nd part //
      /////////////
      // Get the facets of all the vertices of the j-face in B and intersect them to extract the support.
      std::set<unsigned int> listOfFacetsOfVtxB;
      std::set<unsigned int> tmpInterResForB, tmp2InterResForB;
      itVtxNbOperand = listOfVtxB.begin();
      for (unsigned int ngb_count=0; ngb_count<_secondOperand->getGenerator(*itVtxNbOperand)->numberOfFacets(); ngb_count++) {
        boost::shared_ptr<HalfSpace_Rn> Fi = _secondOperand->getGenerator(*itVtxNbOperand)->getFacet(ngb_count);
        listOfFacetsOfVtxB.insert(_secondOperand->getHalfSpaceNumber(Fi));
      }
#ifdef DEBUG
      std::cout << "F_B:";
      std::cout << " V = {";
      std::cout << " " << *itVtxNbOperand;
#endif
      tmpInterResForB = listOfFacetsOfVtxB;
      itVtxNbOperand++;
      {for (; itVtxNbOperand!=listOfVtxB.end(); ++itVtxNbOperand) {
#ifdef DEBUG
        std::cout << " " << *itVtxNbOperand;
#endif
        std::set<unsigned int> curListOfFacetsOfVtxB;
        for (unsigned int ngb_count=0; ngb_count<_secondOperand->getGenerator(*itVtxNbOperand)->numberOfFacets(); ngb_count++) {
          boost::shared_ptr<HalfSpace_Rn> Fi = _secondOperand->getGenerator(*itVtxNbOperand)->getFacet(ngb_count);
          curListOfFacetsOfVtxB.insert(_secondOperand->getHalfSpaceNumber(Fi));
        }
        tmp2InterResForB = tmpInterResForB;
        tmpInterResForB.clear();
        std::set_intersection(
            tmp2InterResForB.begin(),
            tmp2InterResForB.end(),
            curListOfFacetsOfVtxB.begin(),
            curListOfFacetsOfVtxB.end(),
            std::inserter(tmpInterResForB, tmpInterResForB.end()));
        if (listOfFacetsOfVtxB.empty() == true)
          throw std::domain_error("PseudoSumWithoutCaps::computeCapHalfSpaces() the j-Face of polytope B has 0 half-space.");
      }}
#ifdef DEBUG
      std::cout << " }";
#endif
      // Now we have all the geometric supports of F_B, see whether there is a cap half-space inside.
      std::set<unsigned int> InterResForB;
      std::set_intersection(
          tmpInterResForB.begin(),
          tmpInterResForB.end(),
          secondOperandCaps.begin(),
          secondOperandCaps.end(),
          std::inserter(InterResForB, InterResForB.end()));
#ifdef DEBUG
      std::cout << " ; H = { ";
      std::copy(InterResForB.begin(), InterResForB.end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
      //std::cout << listOfFacetsOfVtxB.size();
      std::cout << "}";
      std::cout << endl;
#endif

      // If at least one cap half-space is in listOfFacetsOfVtxA or listOfFacetsOfVtxB then the current half-space of the sum is capped too.
      if (InterResForB.empty() != true)
        sumCaps.insert( sumFacetNb );
    }

    sumFacetNb++;
  }}

}

// Remove the cap half-spaces stored in sets and then truncate again.
boost::shared_ptr<Polytope_Rn> PseudoSumWithoutCaps::rebuildSum(
    const std::set<unsigned int>& firstOperandCaps,
    const std::set<unsigned int>& secondOperandCaps,
    std::set<unsigned int>& newCaps,
    double bb_size)
{
  if (_sum->numberOfHalfSpaces() == 0 || _sum->numberOfGeneratorsPerFacet() == 0)
    throw std::domain_error("PseudoSumWithoutCaps::rebuildSum() _sum is not computed.");

  std::set< unsigned int > sumCaps;
  computeCapHalfSpaces(firstOperandCaps, secondOperandCaps, sumCaps);
#ifdef DEBUG
  std::cout << "firstOperandCaps=" << firstOperandCaps.size() << std::endl;
  std::cout << "secondOperandCaps=" << secondOperandCaps.size() << std::endl;
  std::cout << "sumCaps=" << sumCaps.size() << std::endl;
#endif

  // Now remove all cap half-spaces from the current polytope and create another one without them.
  boost::shared_ptr<Polytope_Rn> newSum;
  newSum.reset(new Polytope_Rn());
  newSum->createBoundingBox(bb_size);
  std::vector< boost::shared_ptr<HalfSpace_Rn> > tryCapsForNewSum;
  // Create the new cap half-spaces.
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS1(newSum->getListOfHalfSpaces());
  for (iteHS1.begin(); iteHS1.end()!=true; iteHS1.next()) {
    tryCapsForNewSum.push_back(iteHS1.current());
  }
  // Copy the non cap half-spaces.
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS2(_sum->getListOfHalfSpaces());
  for (iteHS2.begin(); iteHS2.end()!=true; iteHS2.next()) {
    if (sumCaps.find(iteHS2.currentIteratorNumber()) == sumCaps.end())
      newSum->addHalfSpace(iteHS2.current());
  }

  // Now the truncation.
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > lexmin_ite(newSum->getListOfHalfSpaces());
  StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
  // The number of facets for a cube.
  unsigned int truncationStep = 2*Rn::getDimension();
  DoubleDescription<
    boost::shared_ptr<PolyhedralCone_Rn>,
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
    StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
    DD(newSum, lexmin_ite, NRP, truncationStep);

  // Build the new list of cap half-spaces for _sum
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS3(newSum->getListOfHalfSpaces());
  for (iteHS3.begin(); iteHS3.end()!=true; iteHS3.next()) {
    std::vector< boost::shared_ptr<HalfSpace_Rn> >::const_iterator iteCaps = tryCapsForNewSum.begin();
    for (; iteCaps!=tryCapsForNewSum.end(); ++iteCaps) {
      if (iteHS3.current() == *iteCaps) {
        newCaps.insert(iteHS3.currentIteratorNumber());
        break;
      }
    }
  }
  _sum = newSum;

  return newSum;
}

//========================================================================================

PseudoIntersectionWithoutCaps::PseudoIntersectionWithoutCaps(
   const boost::shared_ptr<Polytope_Rn>& A,
   const boost::shared_ptr<Polytope_Rn>& B,
   boost::shared_ptr<Polytope_Rn>& C,
   const std::set< unsigned int >& firstOperandCaps,
   const std::set< unsigned int >& secondOperandCaps,
   std::set< unsigned int >& newCaps,
   double bb_size) {

  // Prepare the intersection polytope.
  C.reset(new Polytope_Rn());
  C->createBoundingBox(bb_size);
  std::vector< boost::shared_ptr<HalfSpace_Rn> > RnCaps;
  // Get the cube cap half-spaces
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS0(C->getListOfHalfSpaces());
  for (iteHS0.begin(); iteHS0.end()!=true; iteHS0.next()) {
    RnCaps.push_back(iteHS0.current());
  }
  // Copy the non cap half-spaces of A.
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS1(A->getListOfHalfSpaces());
  for (iteHS1.begin(); iteHS1.end()!=true; iteHS1.next()) {
    if (firstOperandCaps.find(iteHS1.currentIteratorNumber()) == firstOperandCaps.end())
      C->addHalfSpace(iteHS1.current());
  }
  // Copy the non cap half-spaces of B.
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS2(B->getListOfHalfSpaces());
  for (iteHS2.begin(); iteHS2.end()!=true; iteHS2.next()) {
    if (secondOperandCaps.find(iteHS2.currentIteratorNumber()) == secondOperandCaps.end())
      C->addHalfSpace(iteHS2.current());
  }

  // Now the truncation.
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > lexmin_ite(C->getListOfHalfSpaces());
  StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
  // The number of facets for a cube.
  unsigned int truncationStep = 2*Rn::getDimension();
  DoubleDescription<
    boost::shared_ptr<PolyhedralCone_Rn>,
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
    StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
    DD(C, lexmin_ite, NRP, truncationStep);

  // Build the new list of cap half-spaces for C
  newCaps.clear();
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS3(C->getListOfHalfSpaces());
  for (iteHS3.begin(); iteHS3.end()!=true; iteHS3.next()) {
    std::vector< boost::shared_ptr<HalfSpace_Rn> >::const_iterator iteCaps = RnCaps.begin();
    for (; iteCaps!=RnCaps.end(); ++iteCaps) {
      if (iteHS3.current() == *iteCaps) {
        newCaps.insert(iteHS3.currentIteratorNumber());
        break;
      }
    }
  }
}

//========================================================================================

int TopGeomTools::Translate(boost::shared_ptr<Polytope_Rn>& pol, const boost::numeric::ublas::vector<double>& v2t) {
  if (pol->numberOfGenerators() != 0) {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(pol->getListOfGenerators());
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
      const boost::shared_ptr<Generator_Rn>& v2_A = iteGN.current();
      boost::numeric::ublas::vector<double> coord = v2_A->vect() + v2t;
      v2_A->setCoordinates(coord);
    }}
  }
  if (pol->numberOfHalfSpaces() != 0) {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS(pol->getListOfHalfSpaces());
    {for (iteHS.begin(); iteHS.end()!=true; iteHS.next()) {
      const boost::shared_ptr<HalfSpace_Rn>& hs = iteHS.current();
      double prod = std::inner_product(hs->begin(), hs->end(), v2t.begin(), 0.);
      hs->setConstant( hs->getConstant()-prod );
    }}
  }
  return 0;
}

int TopGeomTools::GravityCenter(boost::shared_ptr<Polytope_Rn>& pol, boost::numeric::ublas::vector<double>& gravity_center) {
  if (pol->numberOfGenerators() != 0) {
    for (unsigned int i=0; i<gravity_center.size (); ++i)
      gravity_center(i) = 0.;
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(pol->getListOfGenerators());
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
      gravity_center += iteGN.current()->vect();
    }}
    gravity_center /= pol->numberOfGenerators();
  }
  else
    throw std::domain_error("DoubleDescriptionFromGenerators::Compute() the polytope already does not have generators.");
  return 0;
}

int TopGeomTools::PolarPolytope(
    const boost::shared_ptr<Polytope_Rn>& original_pol,
    boost::shared_ptr<Polytope_Rn>& polar_pol,
    bool forceComputation, double bb_size) throw (invalid_argument) {
  polar_pol.reset(new Polytope_Rn());
  unsigned int dim = original_pol->dimension();
  double TOL = Rn::getTolerance();
  if (original_pol->numberOfGenerators() != 0) {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(original_pol->getListOfGenerators());
    {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
      boost::shared_ptr<HalfSpace_Rn> HS;
      HS.reset(new HalfSpace_Rn(dim));
      HS->setConstant(1.);
      boost::numeric::ublas::vector< double >::const_iterator iteCoord = iteGN.current()->begin();
      unsigned int coord_count = 0;
      {for ( ; iteCoord != iteGN.current()->end(); ++iteCoord ) {
        //std::cout << *iteCoord << " ";
        HS->setCoefficient(coord_count, *iteCoord);
        ++coord_count;
      }}
      //std::cout << std::endl;
      polar_pol->addHalfSpace(HS);
    }}
    if (forceComputation == true) {
      try {
        polar_pol->createBoundingBox(bb_size);
        constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > lexmin_ite(polar_pol->getListOfHalfSpaces());
        StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
        // The number of facets for a cube.
        unsigned int truncationStep = 2*dim;
        DoubleDescription<
          boost::shared_ptr<PolyhedralCone_Rn>,
          constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
          StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
            DD(polar_pol, lexmin_ite, NRP, truncationStep);
      }
      catch(invalid_argument& except) {
        cerr << "Invalid argument exception in TopGeomTools::PolarPolytope() " << except.what() << endl;
        return -1;
      }
      catch(out_of_range& except) {
        cerr << "Out of range exception in TopGeomTools::PolarPolytope() " << except.what() << endl;
        return -1;
      }
      catch(ios_base::failure& except) {
        cerr << "In/out exception in TopGeomTools::PolarPolytope() " << except.what() << endl;
        return -1;
      }
      catch(logic_error& except) {
        cerr << "Logic error exception in TopGeomTools::PolarPolytope() " << except.what() << endl;
        return -1;
      }
      catch(...) {
        cerr << "Unexpected exception caught in TopGeomTools::PolarPolytope() !" << endl;
        return -1;
      }
    }
    // Whether there were or not half-spaces, we leave now.
    return 0;
  }
  // At this step we know that generators were not provided.
  if (original_pol->numberOfHalfSpaces() != 0) {
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS(original_pol->getListOfHalfSpaces());
    {for (iteHS.begin(); iteHS.end()!=true; iteHS.next()) {
      if (iteHS.current()->getConstant() < TOL)
        throw std::invalid_argument("TopGeomTools::PolarPolytope() the input polytope should contain the origin.");
      boost::numeric::ublas::vector<double> copyOfHSCoef(iteHS.current()->vect());
      // Scale the hyperplane such as the constant is equal to 1.
      copyOfHSCoef = copyOfHSCoef / iteHS.current()->getConstant();
      boost::shared_ptr<Generator_Rn> GN;
      GN.reset(new Generator_Rn(dim));
      GN->setCoordinates(copyOfHSCoef);
      //if (original_pol->numberOfGenerators() != 0) {
        //constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > polar_iteHS(polar_pol->getListOfHalfSpaces());
        //{for (polar_iteHS.begin(); polar_iteHS.end()!=true; polar_iteHS.next()) {
          //double halfSpaceNorm = std::inner_product(polar_iteHS.current()->begin(), polar_iteHS.current()->end(), polar_iteHS.current()->begin(), 0.);
          //halfSpaceNorm = sqrt(halfSpaceNorm);
          //double scalarProduct = std::inner_product(GN->begin(), GN->end(), polar_iteHS.current()->begin(), 0.);
          //double distanceToHyperplane = (scalarProduct+polar_iteHS.current()->getConstant()) / halfSpaceNorm;
          //if (distanceToHyperplane<TOL && distanceToHyperplane>-TOL)
            //GN->setFacet(polar_iteHS.current());
        //}}
      //}
      polar_pol->addGenerator(GN);
    }}
  }

  return 0;
}

int TopGeomTools::projectPolytopeOnCanonicalHyperplanes(
    const std::set< unsigned int >& listOfHyperplanes,
    const boost::shared_ptr<Polytope_Rn>& original_pol,
    boost::shared_ptr<Polytope_Rn>& proj_pol) throw (invalid_argument)
{
  // Here we really need to work with Rn::getDimension()
  unsigned int currentDimension = Rn::getDimension();
  if (listOfHyperplanes.empty()==true || listOfHyperplanes.size()>Rn::getDimension())
    throw std::invalid_argument("TopGeomTools::projectPolytopeOnCanonicalHyperplanes() wrong list of hyperplanes to project");
  unsigned int projectedDimension = listOfHyperplanes.size();

  Rn::setDimension(projectedDimension);
  proj_pol.reset(new Polytope_Rn());

  std::vector< boost::shared_ptr<Generator_Rn> > newArrayOfGN;
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(original_pol->getListOfGenerators());
  {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
    unsigned int j = 0;
    boost::shared_ptr<Generator_Rn> VX;
    VX.reset(new Generator_Rn(projectedDimension));
    for (unsigned int i=1; i<=currentDimension; ++i) {
      if (listOfHyperplanes.find(i) != listOfHyperplanes.end()) {
        // The direction i is significant so store the corresponding coordinate.
        VX->setCoordinate(j, iteGN.current()->getCoordinate(i-1));
        ++j;
      }
    }
    bool foundEqual = false;
    std::vector< boost::shared_ptr<Generator_Rn> >::const_iterator iteNew = newArrayOfGN.begin();
    for ( ; iteNew != newArrayOfGN.end() && foundEqual == false; ++iteNew) {
      if (VX->distanceFrom(**iteNew) < Rn::getTolerance())
        foundEqual = true;
    }
    if (foundEqual == false)
      newArrayOfGN.push_back( VX );
  }}

  std::vector< boost::shared_ptr<Generator_Rn> >::const_iterator iteNew2 = newArrayOfGN.begin();
  for ( ; iteNew2 != newArrayOfGN.end(); ++iteNew2)
    proj_pol->addGenerator(*iteNew2);
  ///@
  //std::cout << "DUMP" << std::endl;
  //proj_pol->dump(std::cout);

  DoubleDescriptionFromGenerators::Compute(proj_pol, 10000.);
  proj_pol->checkTopologyAndGeometry();

  Rn::setDimension(currentDimension);
  return 0;
}


//========================================================================================

int DoubleDescriptionFromGenerators::Compute(boost::shared_ptr<Polytope_Rn>& pol, double bb_size)
  throw (invalid_argument, out_of_range, ios_base::failure, logic_error) {
  unsigned int dim = pol->dimension();
  //double TOL = Rn::getTolerance();
  if (pol->numberOfHalfSpaces() != 0) 
    throw std::domain_error("DoubleDescriptionFromGenerators::Compute() the polytope already has half-spaces.");
  if (pol->numberOfGenerators() == 0) 
    throw std::domain_error("DoubleDescriptionFromGenerators::Compute() the polytope does not have generators.");
  // Step 1: computes the gravity center
  // Step 2: translate the polytope to have  gravity centered on the origin
  // Step 3: turn the V-description polytope into its H-description polar
  // Step 4: get the  V-description from the H-description polar
  // Step 5: get the polar of the polar
  // Step 6: translate it back
  
  // Step 1: compute the gravity center of the cloud of points.
  boost::numeric::ublas::vector<double> gravity_center(dim);
  TopGeomTools::GravityCenter(pol, gravity_center);
  
  // Step 2: translate the cloud to have in centered on the origin.
  TopGeomTools::Translate(pol, -gravity_center);
  
  // Step 3: compute the polar H-description and truncate it.
  boost::shared_ptr<Polytope_Rn> polar_pol(new Polytope_Rn());
  TopGeomTools::PolarPolytope(pol, polar_pol, true, bb_size);

  // Step 4: compute the polar of the polar.
  pol.reset(new Polytope_Rn());
  TopGeomTools::PolarPolytope(polar_pol, pol);
  
  // Step 5: translate it back.
  TopGeomTools::Translate(pol, gravity_center);
  
  return 0;
}

//========================================================================================

// set palette defined ( 0 "#A90D0D", 0.25 "#22A90D", 0.50 "#0DA9A9", 0.75 "#FFFFFF", 1.0 "#C0C0C0" )
// .. # The background color is set first, then the border colors, then the X & Y axis colors, then the plotting colors
// set terminal png size 1280,960 #ffffff #ffffff #ffffff #000000
// set output "biarn.png"
// plot "biarn.dat" using 3:4:(($5)/1000.):6 with circles linecolor palette z fill solid border lt 1 notitle, 'biarn.dat' using ($3):($4):1 with labels notitle;void Visualization::gnuplot2D(const boost::shared_ptr<Polytope_Rn>& polygon, const std::string& name, double col, std::ostream& out) throw (std::domain_error) {
void Visualization::gnuplot2D(const boost::shared_ptr<Polytope_Rn>& polygon, const std::string& name, double col, std::ostream& out) throw (std::domain_error) {
  if (Rn::getDimension() != 2)
    throw std::domain_error("Visualization::gnuplot(std::ostream& out) dimension is not 2 ");
  out << "set object " << name << " polygon from ";
  boost::shared_ptr<HalfSpace_Rn> oldHS;
  boost::shared_ptr<Generator_Rn> startVTX,referenceVTX;
  constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > LVX_A(polygon->getListOfGenerators());
  startVTX = LVX_A.current();
  std::set< boost::shared_ptr<Generator_Rn> > alreadyProcessedVTX;
  {for (LVX_A.begin(); LVX_A.end()!=true; LVX_A.next()) {
    if (alreadyProcessedVTX.find(LVX_A.current()) == alreadyProcessedVTX.end()) {
      referenceVTX = LVX_A.current();
      alreadyProcessedVTX.insert(referenceVTX);
      boost::shared_ptr<HalfSpace_Rn> referenceHS = referenceVTX->getFacet(0);
      if (referenceHS == oldHS)
        referenceHS = referenceVTX->getFacet(1);
      // Write the point
      out << LVX_A.current()->getCoordinate(0) << ","  << LVX_A.current()->getCoordinate(1) << " to ";
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > LVX_B(polygon->getListOfGenerators());
      {for (LVX_B.begin(); LVX_B.end()!=true; LVX_B.next()) {
        if (alreadyProcessedVTX.find(LVX_B.current()) == alreadyProcessedVTX.end()) {
          if (referenceHS == LVX_B.current()->getFacet(0)) {
            referenceHS = referenceVTX->getFacet(1);
            // Write the point
            out << LVX_B.current()->getCoordinate(0) << ","  << LVX_B.current()->getCoordinate(1) << " to ";
            referenceVTX = LVX_B.current();
            alreadyProcessedVTX.insert(referenceVTX);
          }
          else if (referenceHS == LVX_B.current()->getFacet(1)) {
            referenceHS = referenceVTX->getFacet(0);
            // Write the point
            out << LVX_B.current()->getCoordinate(0) << ","  << LVX_B.current()->getCoordinate(1) << " to ";
            referenceVTX = LVX_B.current();
            alreadyProcessedVTX.insert(referenceVTX);
          }
        }
      }}
    }
  }}
  out << startVTX->getCoordinate(0) << ","  << startVTX->getCoordinate(1);
  out << " lc palette " << col << std::endl;
  //set object 1 fc rgb '#000000' fillstyle solid lw 0
}


