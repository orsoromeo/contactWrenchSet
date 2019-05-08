// politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2012-2016 : Delos Vincent
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
/// \file DoubleDescription_Rn.h
/// \author Delos Vincent (vincent.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef DOUBLEDESCRIPTION_Rn
#define DOUBLEDESCRIPTION_Rn

#include <math.h>
#include <numeric>
#include <set>
#include <map>
#include "Tracking.h"
#include "Neighbours_Rn.h"



/// \brief The algorithm implemented here is an incremental algorithm as mentioned in
/// <i> How Good are Convex Hull Algorithms? </i> (1997) by <b>David Avis</b> and <b>David Bremner</b>.
/// Specific and efficient implementations can be found in <i>The double description method revisited</i>
/// (1996) written by <b> Komei Fukuda </b> and <b> Alain Prodon </b>. <br>
/// Incremental algorithms for the vertex enumeration problem compute the vertex description by
/// intersecting the defining half-spaces sequentially. An initial simplex is constructed from a subset of
/// <i>n+1</i> half-spaces and its vertices and 1-skeleton are computed. Additional half-spaces are introduced
/// sequentially and the vertex description and 1-skeleton are updated at each stage. Essentially such an
/// update amounts to identifying and removing all vertices that are not contained in the new half-space,
/// introducing new vertices for all intersections between edges and the bounding hyperplane of the
/// new half-space, and generating the new edges between these new vertices. <br>
/// This algorithm can be instantiated by polytopes or polyhedral cones, and as a second argument
/// can be instantiated by iterators such as minindex, lexmin, lexmax.
/// <ul> <li> minindex : Insert the half-spaces in the order given by the input. </li>
/// <li> lexmin   : Insert the half-spaces in the the lexicographic increasing order of coefficient vectors. </li>
/// <li> lexmax   : Insert the half-spaces in the the lexicographic decreasing order of coefficient vectors. </li> </ul>
template< class POLYHEDRON, class ITERATOR, class REDUNDANCY_PROCESSING > class DoubleDescription {

 public:
  DoubleDescription(POLYHEDRON poly, ITERATOR ite, REDUNDANCY_PROCESSING redproc, int truncationStep):_redundancyProcessing(redproc) {
    std::vector< boost::shared_ptr<Generator_Rn_SD> > listOfGenSD;
    // Generator_Rn => Generator_Rn_SD
    // From the current generators make a list of generators dedicated to work in the data structure algorithm.
    // Make sure the list of half-space numbers is sorted as we want to use algorithms further.
    poly->getListOfGeneratorsSD(listOfGenSD);
    _redundancyProcessing.initNumberOfVerticesPerHalfSpace(listOfGenSD, truncationStep); // Deal with redundancy.
    //_redundancyProcessing.dumpListOfRedundantHS( std::cout );

    // Do the hard job now.
    _isEmpty = !compute(poly, ite, truncationStep, listOfGenSD);

    if (_isEmpty == false) {
      // Set the pointers from the generators to the half-spaces before dealing with redundancy - as the redundancy
      // processing algorithm will change the numbers of the half-spaces. Here recycle the old pointers stored in
      // the polyhedron to build the new list with the created generators.
      poly->setListOfGeneratorsSD(listOfGenSD);
      // Include all the generators now.
      _redundancyProcessing.unhookRedundantHalfSpaces(poly);
      //poly->dump(std::cout);
    }
  }

  /// Compute the double description tracking all entities and considering the
  /// operator1 as the V-description and operator2 as the H-description.
  DoubleDescription(
      POLYHEDRON poly,
      ITERATOR ite,
      REDUNDANCY_PROCESSING redproc,
      int truncationStep,
      TrackingOperatorToResult& trackerVdesc,
      TrackingOperatorToResult& trackerHdesc):_redundancyProcessing(redproc) {
    std::vector< boost::shared_ptr<Generator_Rn_SD> > listOfGenSD;
    // Generator_Rn => Generator_Rn_SD
    // From the current generators make a list of generators dedicated to work in the data structure algorithm.
    // Make sure the list of half-space numbers is sorted as we want to use algorithms further.
    poly->getListOfGeneratorsSD(listOfGenSD);
    _redundancyProcessing.initNumberOfVerticesPerHalfSpace(listOfGenSD, truncationStep); // Deal with redundancy.
    //_redundancyProcessing.dumpListOfRedundantHS( std::cout );

    // Do the hard job now.
    _isEmpty = !compute(poly, ite, truncationStep, listOfGenSD);

    if (_isEmpty == false) {
      // Store the modifications in the binary operation tracker.
      // Operator1_Entity[i] = (Generator_Rn_SD::Status, Result_Entity[j])
      // Operator2_Entity[k] = (Generator_Rn_SD::Status, Result_Entity[l])
      //    Result_Entity[m] = (Generator_Rn_SD::Status, (Operator1_Entity[p], Operator2_Entity[q]))
      unsigned int nbRes=0;
      std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGenSD;
      {for (iteGenSD=listOfGenSD.begin(); iteGenSD!=listOfGenSD.end(); ++iteGenSD) {
        // In the tracker, the operand entities have been set to Operator_DELETED by default and the
        // result entities have been set by default to Result_UNCHANGED.
        Generator_Rn_SD::Status thisStatus = (*iteGenSD)->getStatus();
        unsigned int nbOp1 = (*iteGenSD)->getGeneratorNumber();
        if (thisStatus == Generator_Rn_SD::CREATED || thisStatus == Generator_Rn_SD::CREATED_AND_MODIFIED) {
          trackerVdesc.setResultEntityAsCreated(nbRes);
          trackerVdesc.setResult_Operator(nbRes,-1);
        }
        else if (thisStatus == Generator_Rn_SD::MODIFIED) {
          trackerVdesc.setResultEntityAsModified(nbRes);
          trackerVdesc.setResult_Operator(nbRes,nbOp1);
          trackerVdesc.setOperatorEntityAsModified(nbOp1);
          trackerVdesc.setOperator_Result(nbOp1,nbRes);
        }
        else if (thisStatus == Generator_Rn_SD::UNCHANGED) {
          trackerVdesc.setResultEntityAsUnchanged(nbRes);
          trackerVdesc.setResult_Operator(nbRes,nbOp1);
          trackerVdesc.setOperatorEntityAsUnchanged(nbOp1);
          trackerVdesc.setOperator_Result(nbOp1,nbRes);
        }
        ++nbRes;
      }}
      /* Set the pointers from the generators to the half-spaces before dealing with redundancy - as the
       *
       * redundancy processing algorithm will change the numbers of the half-spaces. */
      poly->setListOfGeneratorsSD(listOfGenSD);
      // Process the half-spaces now.
      _redundancyProcessing.markHdescription(trackerHdesc, truncationStep);
      _redundancyProcessing.unhookRedundantHalfSpaces(poly);
      //poly->dump(std::cout);
    }
  }

  bool getIsEmpty() const {return _isEmpty;}

  /// \brief For each generator compute its state according to the current half-space.
  void computeVertexStates(
    std::vector< boost::shared_ptr<Generator_Rn_SD> >& GN_list,
    const boost::shared_ptr<HalfSpace_Rn>& currentHalfSpace,
    std::vector<double>& GN_IN_sp,
    std::vector<double>& GN_OUT_sp,
    std::vector< boost::shared_ptr<Generator_Rn_SD> >& GN_IN,
    std::vector< boost::shared_ptr<Generator_Rn_SD> >& GN_OUT,
    std::vector< boost::shared_ptr<Generator_Rn_SD> >& GN_ON)
  {
    int GeneratorNumber=0;
    double TOL = Rn::getTolerance();
    double halfSpaceNorm =
      std::inner_product(currentHalfSpace->begin(), currentHalfSpace->end(), currentHalfSpace->begin(), 0.);
    halfSpaceNorm = sqrt(halfSpaceNorm);
    std::vector< boost::shared_ptr<Generator_Rn_SD> >::iterator iteGN;
    for (iteGN=GN_list.begin(); iteGN!=GN_list.end(); ++iteGN) {
      //const boost::shared_ptr<Generator_Rn_SD>& currentGenerator = iteGN.current();
      double scalarProduct =
          std::inner_product((*iteGN)->begin(), (*iteGN)->end(), currentHalfSpace->begin(), 0.);
          //boost::numeric::ublas::inner_prod(currentGenerator->vect(), currentHalfSpace->vect());
#ifdef DEBUG
      unsigned int RnDIM=currentHalfSpace->dimension();
      std::cout.precision(15);
      std::cout << "# V" << GeneratorNumber << " = [";
      {for (unsigned int ii=0; ii<RnDIM; ii++) {
        std::cout << (*iteGN)->getCoordinate(ii);
        if (ii != RnDIM-1)
          std::cout << ", ";
      }}
      std::cout << "]" << std::endl << "{ ";
      {for (unsigned int ii=0; ii<(*iteGN)->numberOfFacets(); ii++) {
        std::cout << (*iteGN)->getFacet(ii) << " ";
      }}
      std::cout << "}" << std::endl;
      std::cout << "dotP = " << scalarProduct;
      std::cout << ", dist = " << (scalarProduct+currentHalfSpace->getConstant()) / halfSpaceNorm;
#endif
      //std::cout << "sp = " << scalarProduct << " ";
      double distanceToHyperplane = (scalarProduct+currentHalfSpace->getConstant()) / halfSpaceNorm;
      if (distanceToHyperplane > TOL) {
        GN_IN.push_back((*iteGN));
        GN_IN_sp.push_back(scalarProduct);
      }
      else if (distanceToHyperplane < -TOL) {
        GN_OUT.push_back((*iteGN));
        GN_OUT_sp.push_back(scalarProduct);
      }
      else {
        GN_ON.push_back((*iteGN));
      }
#ifdef DEBUG
      HalfSpace_Rn::State currentState = HalfSpace_Rn::hs_UNKNOWN;
      if (distanceToHyperplane > TOL)
        currentState = HalfSpace_Rn::hs_IN;
      else if (distanceToHyperplane < -TOL)
        currentState = HalfSpace_Rn::hs_OUT;
      else
        currentState = HalfSpace_Rn::hs_ON;
      std::cout << ", state = " << HalfSpace_Rn::getStateAsText(currentState) << std::endl;
#endif
      GeneratorNumber++;
    }
  }

  /// \brief The main function splitting the polyhedron cone or polytope 1-skeleton with a list of half-spaces.
  bool compute(POLYHEDRON poly, ITERATOR iteHS, int truncationStep, std::vector< boost::shared_ptr<Generator_Rn_SD> >& listOfGenSD) {
    // Store the scalar product result.
    std::vector<double> GN_IN_sp;
    std::vector<double> GN_OUT_sp;
    // Store generators IN, ON, OUT or newly created.
    std::vector< boost::shared_ptr<Generator_Rn_SD> > GN_IN;
    std::vector< boost::shared_ptr<Generator_Rn_SD> > GN_OUT;
    std::vector< boost::shared_ptr<Generator_Rn_SD> > GN_ON;
    std::vector< boost::shared_ptr<Generator_Rn_SD> > GN_new;
    // Remove redundant half spaces.
    std::vector< boost::shared_ptr<HalfSpace_Rn> > redundantHS;
    int halfspaceNumber=truncationStep, generatorNumber=0;
    unsigned int RnDIM=poly->dimension();
    double TOL2=Rn::getTolerance()*Rn::getTolerance();

#ifdef DEBUG
    std::cout << std::endl << "DoubleDescription::compute(" << truncationStep << ")" << std::endl;
#endif
    // We use this algorithm to compute vertices for a single polyhedron or to intersect two polyhedra.
    // In the last case we truncate the polyhedron A generators net by the the polyhedron B half-spaces
    // so we need to step forward in the half-space list where B constraints begin.
    //constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS(poly->getListOfHalfSpaces());
    generatorNumber = poly->numberOfGenerators();
    iteHS.setStep(truncationStep);
    {for (iteHS.begin(); iteHS.end()!=true; iteHS.next()) {
      const boost::shared_ptr<HalfSpace_Rn>& currentHalfSpace = iteHS.current();
      _redundancyProcessing.addHalfSpace();
#ifdef DEBUG
      double halfSpaceNorm = std::inner_product(currentHalfSpace->begin(), currentHalfSpace->end(), currentHalfSpace->begin(), 0.);
      halfSpaceNorm = sqrt(halfSpaceNorm);
      std::cout << "Facet number " << halfspaceNumber << std::endl;
      std::cout << "H" << halfspaceNumber << " = ";
      currentHalfSpace->dump(std::cout);
      std::cout << ", norm=" << halfSpaceNorm;
      std::cout << std::endl;
#endif

      // Fill the arrays of generators according to their position with respect to the current half-space.
      computeVertexStates(listOfGenSD, currentHalfSpace, GN_IN_sp, GN_OUT_sp, GN_IN, GN_OUT, GN_ON);

      // Check whether the current constraint excludes all Generators.
      if (GN_IN.empty() == true) {
        // NO INTERSECTION
#ifdef DEBUG
        std::cout << "Truncation is empty." << std::endl;
#endif
        poly->reset();
        return false;
      }

      // Check whether the current constraint is redundant. Only available for polytopes.
      // For polyhedral cones, replace Rn::getDimension() by (Rn::getDimension()-1).
      if (GN_OUT.empty() == true) {// && VX_ON.size() < numberOfGeneratorsPerFacet()) {
        // REDUNDANT HALFSPACE (be careful when modifying a list under iterator)
        redundantHS.push_back(currentHalfSpace);
        //_listOfHalfSpaces.removeCurrentHalfSpace(iteHS.currentIteratorNumber());
        //_listOfHalfSpaces.erase(HSPos);
#ifdef DEBUG
        std::cout << "VX_OUT= " << GN_OUT.size() << std::endl;
        std::cout << "VX_IN = " <<  GN_IN.size() << std::endl;
        std::cout << "VX_ON = " <<  GN_ON.size() << std::endl;
        std::cout << "redundantHS = " <<  redundantHS.size() << std::endl;
        std::cout << "============================================" << std::endl;
#endif
      }
      else if (GN_IN.empty() != true) {
        // Build the list of all possible new generators.
        Neighbours_Rn newGeneratorsTool;
        {
          // OUT-IN & OUT-ON
          unsigned int count_OUT=0;
          std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_OUT;
          {for (iteGN_OUT=GN_OUT.begin(); iteGN_OUT!=GN_OUT.end(); ++iteGN_OUT, ++count_OUT) {
            unsigned int count_IN=0;
            std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_IN;
            {for (iteGN_IN=GN_IN.begin(); iteGN_IN!=GN_IN.end(); ++iteGN_IN, ++count_IN) {
              std::vector< unsigned int > commonPFacets;
              if (_redundancyProcessing.checkNeighbours(poly, *iteGN_IN, *iteGN_OUT, commonPFacets) == true)
                newGeneratorsTool.addGenerator(
                    commonPFacets,
                    count_IN,
                    count_OUT,
                    HalfSpace_Rn::hs_IN_OR_OUT);
            }}
            unsigned int count_ON=0;
            std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_ON;
            {for (iteGN_ON=GN_ON.begin(); iteGN_ON!=GN_ON.end(); ++iteGN_ON, ++count_ON) {
              std::vector< unsigned int > commonPFacets;
              if (_redundancyProcessing.checkNeighbours(poly, *iteGN_ON, *iteGN_OUT, commonPFacets) == true)
                newGeneratorsTool.addGenerator(
                    commonPFacets,
                    count_ON,
                    count_OUT,
                    HalfSpace_Rn::hs_ON);
            }}
          }}
          // ON-IN & ON-ON
          unsigned int count_ON1=0;
          std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_ON1;
          {for (iteGN_ON1=GN_ON.begin(); iteGN_ON1!=GN_ON.end(); ++iteGN_ON1, ++count_ON1) {
            unsigned int count_IN=0;
            std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_IN;
            {for (iteGN_IN=GN_IN.begin(); iteGN_IN!=GN_IN.end(); ++iteGN_IN, ++count_IN) {
              std::vector< unsigned int > commonPFacets;
              if (_redundancyProcessing.checkNeighbours(poly, *iteGN_IN, *iteGN_ON1, commonPFacets) == true)
                newGeneratorsTool.addGenerator(
                    commonPFacets,
                    count_IN,
                    count_ON1,
                    HalfSpace_Rn::hs_ON);
            }}
            unsigned int count_ON2=0;
            std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_ON2;
            {for (iteGN_ON2=GN_ON.begin(); iteGN_ON2!=GN_ON.end(); ++iteGN_ON2, ++count_ON2) {
              std::vector< unsigned int > commonPFacets;
              if (iteGN_ON1 != iteGN_ON2 &&
                  _redundancyProcessing.checkNeighbours(poly, *iteGN_ON1, *iteGN_ON2, commonPFacets) == true)
                newGeneratorsTool.addGenerator(
                    commonPFacets,
                    count_ON1,
                    count_ON2,
                    HalfSpace_Rn::hs_ON);
            }}
          }}
        }
        //newGeneratorsTool.dump(std::cout);
        for (newGeneratorsTool.begin(); newGeneratorsTool.end()!=true; newGeneratorsTool.next()) {
#ifdef DEBUG
          std::cout << "NGB1= " << newGeneratorsTool.currentGenInNumber() << ", NGB2= " << newGeneratorsTool.currentGenOutNumber() << std::endl;
#endif
          // Create a new Generator as a combination of currentGeneratorIn and currentGeneratorOut.
          double ay = GN_OUT_sp[newGeneratorsTool.currentGenOutNumber()];
          double az =  GN_IN_sp[newGeneratorsTool.currentGenInNumber()];
          // currentGeneratorIn and currentGeneratorOut are genuine neighbours.
          const boost::shared_ptr<Generator_Rn_SD>& currentGeneratorIn  =  GN_IN[newGeneratorsTool.currentGenInNumber()];
          const boost::shared_ptr<Generator_Rn_SD>& currentGeneratorOut = GN_OUT[newGeneratorsTool.currentGenOutNumber()];
          // Mark the new as CREATED.
          boost::shared_ptr<Generator_Rn_SD> newV(new Generator_Rn_SD(RnDIM, generatorNumber, Generator_Rn_SD::CREATED));
          ++generatorNumber;
          poly->createTruncatedGenerator(currentGeneratorOut, currentGeneratorIn, newV, ay, az, currentHalfSpace->getConstant());
          // Fill in the list of facets for the new Generator.
          std::vector< unsigned int > commonFacets;
          // Get facets.
          _redundancyProcessing.checkNeighbours(poly, currentGeneratorIn, currentGeneratorOut, commonFacets);
          newV->setAllFacets(commonFacets);
          // Do not add the current half-space at this step but wait a little bit more.

          // The new Generator is ready to be inserted in the data structure if not equal to
          // another one with state ON.
          bool notEq=true;
          std::set< unsigned int > setOfFacets;
          std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_ON_eq;
          {for (iteGN_ON_eq=GN_ON.begin(); iteGN_ON_eq!=GN_ON.end() && notEq==true; ++iteGN_ON_eq) {
            if ((*iteGN_ON_eq)->isEqual2(newV, RnDIM, TOL2)) {
              notEq = false;
              // newV has no future so just transfer its facets to the current equal generator.
              // Make sure not to include twice currentHalfSpace.
              newV->exportFacets(setOfFacets);
              (*iteGN_ON_eq)->exportFacets(setOfFacets);
              (*iteGN_ON_eq)->importFacets(setOfFacets);
              (*iteGN_ON_eq)->orderFacets();
            }
          }}
          std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_NEW_eq;
          {for (iteGN_NEW_eq=GN_new.begin(); iteGN_NEW_eq!=GN_new.end() && notEq==true; ++iteGN_NEW_eq) {
            if ((*iteGN_NEW_eq)->isEqual2(newV, RnDIM, TOL2)) {
              notEq = false;
              // newV has no future so just transfer its facets to the current equal generator.
              newV->setFacet(halfspaceNumber);
              newV->exportFacets(setOfFacets);
              (*iteGN_NEW_eq)->exportFacets(setOfFacets);
              (*iteGN_NEW_eq)->importFacets(setOfFacets);
              (*iteGN_NEW_eq)->orderFacets();
            }
          }}
          if (notEq == true) {
            newV->setFacet(halfspaceNumber);
            GN_new.push_back(newV);
#ifdef DEBUG
            std::cout << "New Generator = [";
            newV->save(std::cout);
            std::cout << "] GeneratorIN [";
            currentGeneratorIn->save(std::cout);
            std::cout << "], VX_IN_sp=" << az << "  GeneratorOUT [";
            currentGeneratorOut->save(std::cout);
            std::cout << "], VX_OUT_sp=" << ay << std::endl;
#endif
          }
          else {
#ifdef DEBUG
            std::cout << "No new Generator as [";
            newV->save(std::cout);
            std::cout << "] is equal to a previous one." << std::endl;
#endif
          }
        } // for (newGeneratorsTool.begin(); newGeneratorsTool.end()!=true; newGeneratorsTool.next()) {

#ifdef DEBUG
        std::cout << "VX_OUT= " << GN_OUT.size() << std::endl;
        std::cout << "VX_IN = " <<  GN_IN.size() << std::endl;
        std::cout << "VX_ON = " <<  GN_ON.size() << std::endl;
        std::cout << "redundantHS = " <<  redundantHS.size() << std::endl;
        std::cout << "============================================" << std::endl;
#endif

        // All the generators with state ON must have a new facet.
        std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_ON;
        for (iteGN_ON=GN_ON.begin(); iteGN_ON!=GN_ON.end(); ++iteGN_ON) {
          (*iteGN_ON)->setFacet(halfspaceNumber);
          if ((*iteGN_ON)->getStatus() == Generator_Rn_SD::CREATED)
            (*iteGN_ON)->setStatus(Generator_Rn_SD::CREATED_AND_MODIFIED);
          else if ((*iteGN_ON)->getStatus() == Generator_Rn_SD::UNCHANGED)
            (*iteGN_ON)->setStatus(Generator_Rn_SD::MODIFIED);
          // Insert the ON generators in the IN list as this list will be the official one soon.
          GN_IN.push_back(*iteGN_ON);
        }
        // The current face must mark all the vertices with state ON. 
        // No need to do more than marking as all vertices ON are now embedded in GN_IN.
        _redundancyProcessing.updateNumberOfVerticesPerHalfSpace(halfspaceNumber, GN_ON); // Deal with redundancy.

        std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_new;
        for (iteGN_new=GN_new.begin(); iteGN_new!=GN_new.end(); ++iteGN_new) {
          // Insert the new generators in the IN list as this list will be the official one soon.
          GN_IN.push_back(*iteGN_new);
          _redundancyProcessing.incrementNumberForVerticesForHalfSpace(*iteGN_new); // Deal with redundancy.
        }
        // For the redundancy processing :
        // . all the vertices with state NEW must have their face counter incremented++
        // . all the vertices with state ON  must have their face counter incremented++ (in a global way)
        // . all the vertices with state OUT must have their face counter decremented--

        std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_OUT;
        for (iteGN_OUT=GN_OUT.begin(); iteGN_OUT!=GN_OUT.end(); ++iteGN_OUT) {
          _redundancyProcessing.decrementNumberForVerticesForHalfSpace(*iteGN_OUT); // Deal with redundancy.
        }

        _redundancyProcessing.updateListOfRedundantHalfSpaces(poly->numberOfGeneratorsPerFacet());  // Deal with redundancy.
	
        // Remove all the generators OUT and replace them with the set IN-ON-new.
        listOfGenSD = GN_IN;
      } // else if (GN_IN.empty() != true) {

      GN_OUT_sp.clear();
      GN_OUT.clear();
      GN_IN_sp.clear();
      GN_IN.clear();
      GN_ON.clear();
      GN_new.clear();
      halfspaceNumber++;
      //_redundancyProcessing.dumpSD(std::cout);
    }} // for (iteHS.begin(); iteHS.end()!=true; iteHS.next()) {

    return true;
  }

 protected:
  /// Store the current state of the intersection.
  bool _isEmpty;
  /// This class is dedicated to dealing with redundant half-spaces with the desired policy.
  REDUNDANCY_PROCESSING _redundancyProcessing;

};

/// \brief Makes the assumption we do not need to process redundant half-spaces in a specific way.
template< class POLYHEDRON > class NoRedundancyProcessing {
 public:
  NoRedundancyProcessing() {}

  virtual ~NoRedundancyProcessing() {}

  /// Check whether two generators are neighbors in the context of not taking into account redundancy.
  virtual bool checkNeighbours(
    POLYHEDRON poly,
    const boost::shared_ptr<Generator_Rn>& genIn,
    const boost::shared_ptr<Generator_Rn>& genOut,
    std::vector< boost::shared_ptr<HalfSpace_Rn> >& commonFacets) {
    return poly->checkNeighbours(genIn, genOut, commonFacets);
  }

  /// Check whether two generators are neighbors in the context of not taking into account redundancy.
  virtual bool checkNeighbours(
    POLYHEDRON poly,
    const boost::shared_ptr<Generator_Rn>& genIn,
    const boost::shared_ptr<Generator_Rn>& genOut,
    std::vector< HalfSpace_Rn* >& commonFacets) {
    return poly->checkNeighbours(genIn, genOut, commonFacets);
  }

  /// Only useful in the case of dealing and processing redundancy.
  void initNumberOfVerticesPerHalfSpace(const std::vector< boost::shared_ptr<Generator_Rn> >&) {}

  /// Only useful in the case of dealing and processing redundancy.
  void updateNumberOfVerticesPerHalfSpace(const boost::shared_ptr<HalfSpace_Rn>& , unsigned int) {}

  /// Only useful in the case of dealing and processing redundancy.
  void incrementNumberForVerticesForHalfSpace(const boost::shared_ptr<Generator_Rn>&) {}

  /// Only useful in the case of dealing and processing redundancy.
  void decrementNumberForVerticesForHalfSpace(const boost::shared_ptr<Generator_Rn>&) {}

  /// Only useful in the case of dealing and processing redundancy.
  void updateListOfRedundantHalfSpaces(unsigned int) {}

  void unhookRedundantHalfSpaces(POLYHEDRON) {}

};


/// \brief This class can be more time-consuming than WeakRedundancyProcessing or NoRedundancyProcessing because it will perform extra
/// checks in the process of intersecting half-spaces. To determine if two vertices are neighbors it will make sure not to count
/// half-spaces marked as redundant.
template< class POLYHEDRON > class StrongRedundancyProcessing {

 public:
  StrongRedundancyProcessing():_numberOfHalfSpaces(0) {}

  virtual ~StrongRedundancyProcessing() {}

  void fillNumberOfVerticesPerHalfSpace(
    const POLYHEDRON,
    std::vector< unsigned int >& getNumberOfVerticesPerHalfSpace) {
    getNumberOfVerticesPerHalfSpace = _numberOfVerticesPerHalfSpace;
  }

  void fillListOfRedundantHS(
    const POLYHEDRON,
    std::vector< unsigned int >&,
    std::set< unsigned int >& getListOfRedundantHS) {
    getListOfRedundantHS = _listOfRedundantHS;
  }

  /// Make sure all back pointers from half-spaces to vertices are set.
  void initNumberOfVerticesPerHalfSpace(const std::vector< boost::shared_ptr<Generator_Rn_SD> >& LG, unsigned int nbHS) {
    _allGenPerHS.clear();
    _listOfRedundantHS.clear();
    _numberOfVerticesPerHalfSpace.clear();
    _numberOfHalfSpaces = nbHS;

    // First init to 0 to be able to increment
    {for (unsigned int j=0; j<nbHS; ++j) {
      _numberOfVerticesPerHalfSpace.push_back(0);
    }}
    std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator constITER_GN;
    for (constITER_GN=LG.begin(); constITER_GN!=LG.end(); ++constITER_GN) {
      {for (unsigned int j=0; j<(*constITER_GN)->numberOfFacets(); j++) {
        unsigned int nbF = (*constITER_GN)->getFacet(j);
        _numberOfVerticesPerHalfSpace[nbF]++;
      }}
    }
    //std::map< boost::shared_ptr<HalfSpace_Rn>, unsigned int >::const_iterator
    	//mit(numberOfVertexPerHalfSpace.begin()), mend(numberOfVertexPerHalfSpace.end());
    //for (; mit!=mend; ++mit)
      //std::cout << mit->first << '\t' << mit->second << std::endl;

    // Now deal with the set of vertices itself contained by the faces
    // to see if a given set can be included into another one. But first of
    // all set the list of back pointers.
    {for (unsigned int j=0; j<nbHS; ++j) {
      std::set< unsigned int > genSet;
      _allGenPerHS.push_back(genSet);
    }}
    for (constITER_GN=LG.begin(); constITER_GN!=LG.end(); ++constITER_GN) {
      {for (unsigned int j=0; j<(*constITER_GN)->numberOfFacets(); j++) {
        unsigned int F = (*constITER_GN)->getFacet(j);
        // GenPerFacet: number(HalfSpace_Rn) => {Generator_Rn_SD_0*, Generator_Rn_SD_1*, Generator_Rn_SD_2*, ...}
        // Fill the data structure storing all the generators for a given facet.
        _allGenPerHS[F].insert( (*constITER_GN)->getGeneratorNumber() );
      }}
    }
    //dumpSD(std::cout);
  }

  /// Make space for a new half-space.
  void addHalfSpace() {
    _numberOfVerticesPerHalfSpace.push_back(0);
    std::set< unsigned int > genSet;
    _allGenPerHS.push_back(genSet);
    _numberOfHalfSpaces++;
  }

  /// The current face must mark all the vertices with state ON.
  void updateNumberOfVerticesPerHalfSpace(
    unsigned int HS,
    const std::vector< boost::shared_ptr<Generator_Rn_SD> >& GN_ON) {
    _numberOfVerticesPerHalfSpace[ HS ] += GN_ON.size();
    std::vector< boost::shared_ptr<Generator_Rn_SD> >::const_iterator iteGN_ON;
    for (iteGN_ON=GN_ON.begin(); iteGN_ON!=GN_ON.end(); ++iteGN_ON)
      _allGenPerHS[ HS ].insert( (*iteGN_ON)->getGeneratorNumber() );
  }

  /// Make sure all the half-spaces belonging to a given generator have their vertices number incremented.
  void incrementNumberForVerticesForHalfSpace(const boost::shared_ptr<Generator_Rn_SD>& GEN) {
    for (unsigned int l=0; l<GEN->numberOfFacets(); ++l) {
      _numberOfVerticesPerHalfSpace[ GEN->getFacet(l) ]++;
      _allGenPerHS[ GEN->getFacet(l) ].insert( GEN->getGeneratorNumber() );
    }
  }

  /// Make sure all the half-spaces belonging to a given generator have their vertices number decremented.
  void decrementNumberForVerticesForHalfSpace(const boost::shared_ptr<Generator_Rn_SD>& GEN) {
    for (unsigned int l=0; l<GEN->numberOfFacets(); ++l) {
      _numberOfVerticesPerHalfSpace[ GEN->getFacet(l) ]--;
      _allGenPerHS[ GEN->getFacet(l) ].erase( GEN->getGeneratorNumber() );
    }
  }

  /// Call poly->checkNeighbours() with an extra argument not to count redundant
  /// half-spaces in the process of declaring two vertices as neighbors.
  bool checkNeighbours(
    POLYHEDRON poly,
    const boost::shared_ptr<Generator_Rn_SD>& genIn,
    const boost::shared_ptr<Generator_Rn_SD>& genOut,
    std::vector< unsigned int >& commonFacets) {
    // Compute the intersection on sorted lists of facets.
    std::set_intersection(genIn->facetsBegin(), genIn->facetsEnd(), genOut->facetsBegin(), genOut->facetsEnd(),
      std::inserter(commonFacets, commonFacets.end()));
    // For polyhedral cones neigbourhoodCondition() is Rn::getDimension()-1, for polytopes it is Rn::getDimension()-2.
    if (commonFacets.size() >= poly->neigbourhoodCondition())
      return true;

    return false;
  }

  /// Inspect all half-spaces, find their lists of generators and numbers of generators
  /// to know whether they are redundant or not.
  virtual void updateListOfRedundantHalfSpaces(unsigned int numberOfGeneratorsPerFacet) {
    {for(unsigned int i=0; i<_numberOfVerticesPerHalfSpace.size(); ++i) {
      // Stop when we reach the current half-space being processed, it was the last inserted.
      if (_numberOfVerticesPerHalfSpace[i] < numberOfGeneratorsPerFacet) {
        if (_listOfRedundantHS.insert( i ).second == true) {
#ifdef DEBUG
          std::cout << "Redundant " << i;
          std::cout << std::endl;
#endif
        }
      }
    }}

    // Now we've checked that half-spaces have the correct minimum number of generators we need to 
    // perform a more complicated task. Now run the map to see if a set of generators can be
    // included into another one. In this case it means that the current half-space is redundant.
    unsigned int i1,i2;
    std::vector< std::set< unsigned int > >::const_iterator it_1, it_2;
    {for (it_1=_allGenPerHS.begin(), i1=0; it_1!=_allGenPerHS.end(); ++it_1, ++i1) {
      {for (it_2=_allGenPerHS.begin(), i2=0; it_2!=_allGenPerHS.end(); ++it_2, ++i2) {
        if (it_1!=it_2 && _listOfRedundantHS.find(i1)==_listOfRedundantHS.end() && _listOfRedundantHS.find(i2)==_listOfRedundantHS.end()) {
          if (it_1->size() > it_2->size()) {
            if (std::includes(it_1->begin(), it_1->end(), it_2->begin(), it_2->end())) {
              _listOfRedundantHS.insert(i2);
#ifdef DEBUG
              std::cout << "New redundant facet: " << i2 << std::endl;
#endif
            }
          }
          else if (it_1->size() < it_2->size()) {
            if (std::includes(it_2->begin(), it_2->end(), it_1->begin(), it_1->end())) {
              _listOfRedundantHS.insert(i1);
#ifdef DEBUG
              std::cout << "New redundant facet: " << i1 << std::endl;
#endif
            }
          }
          else {
            // Unnecessary check for equality
          }
        }
      }}
    }}
  }

  void fillListOfRedundantHS(const POLYHEDRON poly) { updateListOfRedundantHalfSpaces(poly->numberOfGeneratorsPerFacet()); }

  void unhookRedundantHalfSpaces(POLYHEDRON poly) {
    //fillListOfRedundantHS( poly );
    // Unhook redundant half-spaces from the vertices.
    unsigned int hsNb;
    if (!_listOfRedundantHS.empty()) {
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS(poly->getListOfHalfSpaces());
      {for (iteHS.begin(), hsNb=0; iteHS.end()!=true; iteHS.next(), ++hsNb) {
        if (_listOfRedundantHS.find(hsNb) != _listOfRedundantHS.end()) {
          // We have a genuine redundant half-space.
          constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(poly->getListOfGenerators());
          {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
            boost::shared_ptr<Generator_Rn> gen = iteGN.current();
            {for (int j=gen->numberOfFacets()-1; j>=0; j--) {
              if (iteHS.current() == gen->getFacet(j))
                gen->removeFacet(j);
            }}
          }}
        }
      }}
    }
    // Add to the list of redundant half-spaces the ones which are not
    // referenced at all by vertices from the beginning of the process.
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS2(poly->getListOfHalfSpaces());
    {for (iteHS2.begin(), hsNb=0; iteHS2.end()!=true; iteHS2.next(), ++hsNb) {
      if (_listOfRedundantHS.find(hsNb) == _listOfRedundantHS.end()) {
        if (_numberOfVerticesPerHalfSpace[hsNb] == 0) {
          // The half-space has not been found so it's not in use at all and is redundant.
          _listOfRedundantHS.insert(hsNb);
        }
      }
    }}
    if (!_listOfRedundantHS.empty()) {
#ifdef DEBUG
      std::cout << "Remove " << _listOfRedundantHS.size() << " redundant facets: ";
#endif
      std::vector< unsigned int > HS2Remove;
      std::set< unsigned int >::const_iterator it;
      {for (it=_listOfRedundantHS.begin(); it!=_listOfRedundantHS.end(); ++it) {
#ifdef DEBUG
        std::cout << *it << " ";
#endif
        HS2Remove.push_back(*it);
      }}
      std::sort(HS2Remove.begin(), HS2Remove.end(), std::greater<unsigned int>());
      std::vector< unsigned int >::const_iterator it2;
      {for (it2=HS2Remove.begin(); it2!=HS2Remove.end(); ++it2) {
        // Do not forget to begin by the bigger to remove as we could change the ordering.
        poly->removeHalfSpace(*it2);
      }}
#ifdef DEBUG
      std::cout << std::endl;
#endif
    }
  }

  // Mark as UNCHANGED, CREATED or DELETED the half-spaces going through the double description process.
  void markHdescription(TrackingOperatorToResult& trackerHdesc, unsigned int truncationStep) {
    trackerHdesc.setNumbersOfEntities(_numberOfHalfSpaces, _numberOfHalfSpaces);
    // Before truncation step they are UNCHANGED by default, then they are CREATED by default.
    for (unsigned int i=0; i<truncationStep; ++i)
      trackerHdesc.setOperatorEntityAsUnchanged(i);
    for (unsigned int i=truncationStep; i<_numberOfHalfSpaces; ++i)
      trackerHdesc.setResultEntityAsCreated(i);
    std::set< unsigned int >::const_iterator iteRED=_listOfRedundantHS.begin();
    for ( ; iteRED!=_listOfRedundantHS.end(); ++iteRED)
      trackerHdesc.setOperatorEntityAsDeleted(*iteRED);
  }

  void dumpListOfRedundantHS(POLYHEDRON poly, std::ostream &this_stream) {
    unsigned int RnDIM=poly->dimension();
    this_stream << "List of redundant half-spaces :" << std::endl;
    unsigned int hsNb;
    if (!_listOfRedundantHS.empty()) {
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHS(poly->getListOfHalfSpaces());
      {for (iteHS.begin(), hsNb=0; iteHS.end()!=true; iteHS.next(), ++hsNb) {
        if (_listOfRedundantHS.find(hsNb) != _listOfRedundantHS.end()) {
          boost::shared_ptr<HalfSpace_Rn> currentHalfSpace = iteHS.current();
          this_stream << "H = (" << currentHalfSpace->getConstant() << ", ";
          {for (unsigned int ii=0; ii<RnDIM; ii++) {
            this_stream << currentHalfSpace->getCoefficient(ii);
            if (ii != RnDIM-1)
              this_stream << ", ";
          }}
          this_stream << ")" << std::endl;
        }
      }}
    }
  }

  void dumpSD(std::ostream &this_stream) {
    std::vector< std::set< unsigned int > >::const_iterator iteGHS;
    this_stream << "*** number(HalfSpace_Rn0) => {Generator_Rn0*, Generator_Rn1*, Generator_Rn2*, ...}" << std::endl;
    unsigned int counter=0;
    {for (iteGHS=_allGenPerHS.begin(); iteGHS!=_allGenPerHS.end(); ++iteGHS) {
      this_stream << counter <<  " => { ";
      std::copy(iteGHS->begin(), iteGHS->end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
      this_stream << "}" << std::endl;
      ++counter;
    }}
    this_stream << "*** number of generators per half-spaces : ";
    std::copy(_numberOfVerticesPerHalfSpace.begin(), _numberOfVerticesPerHalfSpace.end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
    this_stream << std::endl;
    this_stream << "*** list of redundant half-spaces : ";
    std::copy(_listOfRedundantHS.begin(), _listOfRedundantHS.end(), std::ostream_iterator<unsigned int>(std::cout, " ") );
    this_stream << std::endl;
  }

 protected:
  // Deal with the number of generators per half-space.
  /// Store all raw back pointers to know which vertices belong to a given half-space.
  // number(HalfSpace_Rn0) => {Generator_Rn0*, Generator_Rn1*, Generator_Rn2*, ...}
  std::vector< std::set< unsigned int > > _allGenPerHS;
  /// To know about how many vertices refer to a given half-space.
  std::vector< unsigned int > _numberOfVerticesPerHalfSpace;
  /// To know whether an half-space has been ticked redundant or not.
  std::set< unsigned int > _listOfRedundantHS;
  /// The total number of processed half-spaces.
  unsigned int _numberOfHalfSpaces;
};

#endif
