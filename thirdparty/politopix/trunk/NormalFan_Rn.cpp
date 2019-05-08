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
/// \file NormalFan_Rn.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include "GeometricObjectIterator_Rn.h"
#include "DoubleDescription_Rn.h"
#include "NormalFan_Rn.h"
#include "IO_Polytope.h"
#include <math.h>
#include <numeric>
#include <sstream>


// Constructor
NormalFan_Rn::NormalFan_Rn(const boost::shared_ptr<Polytope_Rn>& A) {
  boost::shared_ptr<PolyhedralCone_Rn> primalCone;
  for (unsigned int u=0; u<A->numberOfGenerators(); u++) {
    boost::shared_ptr<Generator_Rn> vx1 = A->getGenerator(u);
    //////////////////////////////////////////////////////////////////
    // STEP1 : split the polytope into its primal polyhedral cones //
    ////////////////////////////////////////////////////////////////
    primalCone.reset(new PolyhedralCone_Rn());
    // Insert all half-spaces connected to the current vertex into the primal cone.
    for (unsigned int i=0; i<vx1->numberOfFacets(); i++)
      primalCone->addHalfSpace(vx1->getFacet(i));
    //////////////////////////////////////////////
    // STEP 2 : build the dual polyhedral cone //
    ////////////////////////////////////////////
    boost::shared_ptr<PolyhedralCone_Rn> dualCone = primalCone->computeDualPolyhedralCone();
    // Insert the dual in its list.
    _listOfPolyhedralCones.push_back(dualCone);
    // Set the anchor for the corresponding cone.
    _listOfVertices.push_back(vx1);
    //std::cout << "$$$$$$$$$$$$ " << u << std::endl;
    //dualCone->dump(std::cout);
    //dualCone->checkTopologyAndGeometry();
    //for (unsigned int k=0; k<Rn::getDimension(); k++) {
      //vx1->setCoordinate(k, vx1->getCoordinate(k));
      //std::cout << vx1->getCoordinate(k) << " ";
    //}
    //std::cout << std::endl;
    //dualCone->checkTopologyAndGeometry();
    //std::cout << "###########################################"<< std::endl;
    //std::ostringstream stream_;
    //stream_ << "cone_p4ns";
    //stream_ << u;
    //std::string FileOut = stream_.str();
    //FileOut += ".pcon";
    //IO_Polytope::save(FileOut, dualCone);
  }
}

// CHECK POLYHEDRON
bool NormalFan_Rn::checkTopologyAndGeometry() const throw (std::domain_error) {
  return true;
}

void NormalFan_Rn::computeHyperplanesSeparationForProjection(
    const std::vector< boost::shared_ptr<HalfSpace_Rn> >& allHS,
    boost::shared_ptr<Polytope_Rn>& Pol) {
  unsigned int RnDIM = Rn::getDimension();
  double TOL = Rn::getTolerance();
  double TTOL = 2*TOL;
  Rn::setTolerance(TOL);
  std::vector< double > arrayOfNorms(allHS.size());
  std::vector< boost::shared_ptr<HalfSpace_Rn> >::const_iterator iteHS;
  unsigned int hs = 0;
  {for (iteHS = allHS.begin(), hs=0; iteHS!=allHS.end(); ++iteHS, ++hs) {
    double halfSpaceNorm = std::inner_product((*iteHS)->begin(), (*iteHS)->end(), (*iteHS)->begin(), 0.);
    arrayOfNorms[hs] = sqrt(halfSpaceNorm);
  }}
  unsigned int cone = 0;
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> >::iterator itePC;
  {for (itePC=_listOfPolyhedralCones.begin(), cone=0; itePC!=_listOfPolyhedralCones.end(); ++itePC, ++cone) {
    ///@
    //std::cout << std::endl << "### Cone" << cone << std::endl;
    const listOfGeometricObjects< boost::shared_ptr<Generator_Rn> >& listOfEdges = (*itePC)->getListOfGenerators();
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGN(listOfEdges);
    std::vector< bool > arrayOfEgdesOKforHS(allHS.size());
    {for (unsigned int i=0; i<arrayOfEgdesOKforHS.size(); ++i)
      arrayOfEgdesOKforHS[i] = false;
    }
    bool coneOKforHS = true;
    {for (iteHS = allHS.begin(), hs=0; iteHS!=allHS.end() && coneOKforHS==true; ++iteHS, ++hs) {
      unsigned int countON = 0;
      bool IN = false, OUT = false;
      {for (iteGN.begin(); iteGN.end()!=true; iteGN.next()) {
        double scalarProduct = std::inner_product(iteGN.current()->begin(), iteGN.current()->end(), (*iteHS)->begin(), 0.);
        double distanceToHyperplane = scalarProduct / arrayOfNorms[hs];

        ///@
        //std::cout.precision(15);
        //std::cout << "# E" << iteGN.currentIteratorNumber() << " = [";
        //{for (unsigned int ii=0; ii<RnDIM; ii++) {
          //std::cout << iteGN.current()->getCoordinate(ii);
          //if (ii != RnDIM-1)
            //std::cout << ", ";
        //}}
        //std::cout << "] => ";
        //(*iteHS)->dump(std::cout);
        //std::cout << ", dotP = " << scalarProduct;
        //std::cout << ", dist = " << distanceToHyperplane;
        ///

        if (distanceToHyperplane > TOL) {
          IN = true;

          ///@
          //std::cout << " (IN)" << std::endl;
          ///

        }
        else if (distanceToHyperplane < -TOL) {
          OUT = true;

          ///@
          //std::cout << " (OUT)" << std::endl;
          ///

        }
        else {
          ++countON;

          ///@
          //std::cout << " (ON)" << std::endl;
          ///

        }
      }}
      if ((IN==true && OUT==true) || (countON>=1))
        arrayOfEgdesOKforHS[hs] = true;
      else
        coneOKforHS = false;
    }}
    bool allHSOK = true;
    ///@
    //std::cout << "Array of boolean {";
    {for (unsigned int i=0; i<arrayOfEgdesOKforHS.size(); ++i) {
      ///@
      //std::cout << " " << (arrayOfEgdesOKforHS[i]==true? "T" : "F");
      if (arrayOfEgdesOKforHS[i] == false)
        allHSOK = false;
    }}
    ///@
    //std::cout << " }, result = " << (allHSOK==true? "T" : "F") << std::endl;
    //_listOfVertices[cone]->dump(std::cout);
    //std::cout<< std::endl;
    if (allHSOK == true)
      Pol->addGenerator(_listOfVertices[cone]);
  }}

  // Restore the previous configuration.
  Rn::setTolerance(TOL);
}


// ALGORITHMS
void NormalFan_Rn::computeCommonRefinement(const NormalFan_Rn& NF1, const NormalFan_Rn& NF2) {
#ifdef DEBUG
  std::cout << std::endl << "NormalFan_Rn::computeCommonRefinement()" << std::endl;
#endif
  unsigned int RnDIM=Rn::getDimension();
  unsigned int currentCone1Nb=0,currentCone2Nb=0;
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> >::const_iterator endLPC1=NF1.getListOfPolyhedralCones().end();
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> >::const_iterator endLPC2=NF2.getListOfPolyhedralCones().end();
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> >::const_iterator LPC1=NF1.getListOfPolyhedralCones().begin();
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> >::const_iterator LPC2=NF2.getListOfPolyhedralCones().begin();
  std::vector< boost::shared_ptr<Generator_Rn> >::const_iterator LVX1=NF1.getListOfGenerators().begin();
  std::vector< boost::shared_ptr<Generator_Rn> >::const_iterator LVX2=NF2.getListOfGenerators().begin();
  {for (; LPC1!=endLPC1; ++LPC1, ++LVX1) {
    {for (; LPC2!=endLPC2; ++LPC2, ++LVX2) {
      boost::shared_ptr<PolyhedralCone_Rn> intersectionCone(new PolyhedralCone_Rn());
      boost::shared_ptr<PolyhedralCone_Rn> Ca = *LPC1;
      boost::shared_ptr<PolyhedralCone_Rn> Cb = *LPC2;
      // Fill the data structures.
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<Generator_Rn> > iteGNA(Ca->getListOfGenerators());
      {for (iteGNA.begin(); iteGNA.end()!=true; iteGNA.next()) {
        // Make a deep copy of each generator
        boost::shared_ptr<Generator_Rn> gn(new Generator_Rn( *(iteGNA.current()) ));
        intersectionCone->addGenerator(gn);
      }}
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHSB(Cb->getListOfHalfSpaces());
      for (iteHSB.begin(); iteHSB.end()!=true; iteHSB.next())
        intersectionCone->addHalfSpace(iteHSB.current());
#ifdef DEBUG
      std::cout << "=== Cone from A number " << currentCone1Nb << std::endl;
      std::cout << "=== Cone from B number " << currentCone2Nb << std::endl;
#endif
      // Compute the intersection between the two polyhedral cones.
      TrackingOperatorToResult trackerVdesc,trackerHdesc;
      constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > lexmin_ite(intersectionCone->getListOfHalfSpaces());
      StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
      DoubleDescription< boost::shared_ptr<PolyhedralCone_Rn>, constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
        StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > > DD(intersectionCone, lexmin_ite, NRP, 0, trackerVdesc, trackerHdesc);
      bool notEmpty = !DD.getIsEmpty();
#ifdef DEBUG
      std::cout << "@@@ nbEDG=" << intersectionCone->numberOfGenerators() << "( ";
#endif
      // The second test is necessary when dealing with parallelism.
      if (notEmpty == true && intersectionCone->numberOfGenerators() >= RnDIM) {
      	// Set the anchor for the corresponding cone.
        boost::shared_ptr<Generator_Rn> VX(new Generator_Rn(RnDIM));
        VX->makeSum((*LVX1),(*LVX2));
#ifdef DEBUG
        for (unsigned int k=0; k<RnDIM; ++k)
      		std::cout << (*LVX1)->getCoordinate(k) << "+" << (*LVX2)->getCoordinate(k) << "=" << VX->getCoordinate(k) << "\t";
#endif
        _listOfVertices.push_back(VX);
        _listOfPolyhedralCones.push_back(intersectionCone);
        // Now deal with generators MODIFIED and CREATED.
        const std::vector< std::pair< StatusAfter, int > >& res = trackerVdesc.getResultToOperator();
        std::vector< std::pair< StatusAfter, int > >::const_iterator iteGN_INTER=res.begin();
        unsigned int counterV=0;
        // The generators created by the intersection have to be added just the way they are.
        // The generators modified by the intersection have to keep their old half-spaces and add the new ones.
        {for ( ; iteGN_INTER!=res.end(); ++iteGN_INTER) {
          if (iteGN_INTER->first == Result_CREATED) {
            Ca->addGenerator(intersectionCone->getGenerator(counterV));
          }
          else if (iteGN_INTER->first == Result_MODIFIED) {
            int nbModGEN = iteGN_INTER->second;
            // At this moment the new half-spaces are not marked as such, so just use the properties of the sets.
            std::set< boost::shared_ptr<HalfSpace_Rn> > setOfFacets1,setOfFacets2;
            Ca->getGenerator(nbModGEN)->exportFacets(setOfFacets1);
            intersectionCone->getGenerator(counterV)->exportFacets(setOfFacets2);
            setOfFacets1.insert(setOfFacets2.begin(), setOfFacets2.end());
            Ca->getGenerator(nbModGEN)->clearFacets();
            Ca->getGenerator(nbModGEN)->importFacets(setOfFacets1);
          }
          ++counterV;
        }}
        const std::vector< std::pair< StatusAfter, int > >& resH = trackerHdesc.getResultToOperator();
        std::vector< std::pair< StatusAfter, int > >::const_iterator iteHS_INTER=resH.begin();
        unsigned int counterH=0;
        // The new list of half-spaces contain the old list plus the ones of Cb which are not redundant.
        {for ( ; iteHS_INTER!=res.end(); ++iteHS_INTER) {
          if (iteHS_INTER->first == Result_CREATED) {
            Ca->addHalfSpace(intersectionCone->getHalfSpace(counterH));
          }
          ++counterH;
        }}
      }
#ifdef DEBUG
      std::cout << ")" << std::endl;
#endif
      ++currentCone2Nb;
    }}
    ++currentCone1Nb;
  }}
}

/// Dump the polyhedral structure on std::cout.
void NormalFan_Rn::dump(std::ostream& out) const {
  out << std::endl << "#NORMAL FAN" << std::endl;
  unsigned int nb=0;
  std::vector< boost::shared_ptr<Generator_Rn> >::const_iterator iteVX=_listOfVertices.begin();
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> >::const_iterator itePC=_listOfPolyhedralCones.begin();
  {for (; itePC!=_listOfPolyhedralCones.end(); ++itePC, ++iteVX) {
    out << "#" << nb << std::endl;
    out << "# anchor (";
    for (unsigned int k=0; k<Rn::getDimension(); k++)
      out << (*iteVX)->getCoordinate(k) << " ";
    out << ")" << std::endl;
    (*itePC)->dump(std::cout);
    ++nb;
  }}
}
