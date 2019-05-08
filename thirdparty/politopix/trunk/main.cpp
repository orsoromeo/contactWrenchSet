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
/// \file main.cpp
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#include <boost/shared_ptr.hpp>
#include <boost/timer.hpp>
#include <iostream>
#include <string.h>
#include <cmath>
#include "PolyhedralAlgorithms_Rn.h"
#include "Config.h"


/// This class stores static values to define the running context.
class RunOptions {

 public:
  RunOptions() {}

  /// Set the number of generators test after the algorithm.
  static void setCheckGeneratorsTest(bool b) {_checkGeneratorsTest = b;}

  /// Shall we test the number of generators after the algorithm ?
  static bool getCheckGeneratorsTest() {return _checkGeneratorsTest;}

  /// Set the number of generators to check after the algorithm.
  static void setCheckGeneratorsValue(unsigned int d) {_checkGeneratorsValue = d;}

  /// Get the number of generators to find after the algorithm.
  static unsigned int getCheckGeneratorsValue() {return _checkGeneratorsValue;}

 protected:
  static bool _checkGeneratorsTest;
  static unsigned int _checkGeneratorsValue;

};


bool RunOptions::_checkGeneratorsTest = false;
unsigned int RunOptions::_checkGeneratorsValue = 0;


int main(int argc, char* argv[]) {

  unsigned int dimension = 0;
  double tolerance= 0.;
  Rn::setDimension(0);
  std::string fileName1,fileName2;
  bool boundingBox = true;
  bool boundingVolume = false;
  double boundingSize = 1000;
  double volumeOfPolytope = -1.;
  unsigned int facetsToCheck = 0, generatorsToCheck = 0;
  int vertexToComputeDistances=-1, facetToComputeDistances=-1;
  bool computeDistances=false;
  bool p1=false, p2=false, c1=false, c2=false;
  bool COMPUTE_VOL=false, SUBSET=false;
  bool INTER=false, CHCK_EQ=false, MS=false, check_all=false, output=false;
  std::string outputFileName;

  std::ostringstream stream_;
  stream_ << "Version ";
  stream_ << politopix_VERSION_MAJOR;
  stream_ << ".";
  stream_ << politopix_VERSION_MINOR;
  stream_ << ".";
  stream_ << politopix_VERSION_PATCH;
  std::string version = stream_.str();

  // Parse the main arguments
  for(int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-p1") == 0 || strcmp(argv[i], "--polytope1") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid file " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      // Get file
      char* fileName_1 = argv[++i];
      std::string file1(fileName_1);
      fileName1 = file1;
      p1 = true;
    }
    else if (strcmp(argv[i], "-p2") == 0 || strcmp(argv[i], "--polytope2") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid file " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      char* fileName_2 = argv[++i];
      std::string file2(fileName_2);
      fileName2 = file2;
      p2 = true;
    }
    else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid file " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      char* fileName_3 = argv[++i];
      std::string file3(fileName_3);
      outputFileName = file3;
      output = true;
    }
    else if (strcmp(argv[i], "-c1") == 0 || strcmp(argv[i], "--polyhedralcone1") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid file " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      // Get file
      char* fileName_1 = argv[++i];
      std::string file1(fileName_1);
      fileName1 = file1;
      c1 = true;
    }
    else if (strcmp(argv[i], "-c2") == 0 || strcmp(argv[i], "--polyhedralcone2") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid file " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      char* fileName_2 = argv[++i];
      std::string file2(fileName_2);
      fileName2 = file2;
      c2 = true;
    }
    else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) {
      cout << version << std::endl;
      return EXIT_SUCCESS;
    }
    else if (strcmp(argv[i], "-ch") == 0 || strcmp(argv[i], "--check-all") == 0) {
      // We perform all checks, it can be slow.
      check_all = true;
    }
    else if (strcmp(argv[i], "-MS") == 0 || strcmp(argv[i], "--MinkowskiSum") == 0) {
      // We compute Minkowski sums.
      MS = true;
      INTER = false;
      CHCK_EQ = false;
    }
    else if (strcmp(argv[i], "-IN") == 0 || strcmp(argv[i], "--Intersection") == 0) {
      // We intersect polytopes or polyhedral cones, this is the default behaviour.
      MS = false;
      INTER = true;
      CHCK_EQ = false;
    }
    else if (strcmp(argv[i], "-SS") == 0 || strcmp(argv[i], "--Subset") == 0) {
      // We intersect polytopes or polyhedral cones, this is the default behaviour.
      SUBSET = true;
      MS = false;
      INTER = false;
      CHCK_EQ = false;
    }
    else if (strcmp(argv[i], "-EQ") == 0 || strcmp(argv[i], "--Equality") == 0) {
      // We intersect polytopes or polyhedral cones, this is the default behaviour.
      MS = false;
      INTER = false;
      CHCK_EQ = true;
    }
    else if (strcmp(argv[i], "-bb") == 0 || strcmp(argv[i], "--boundingbox") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid bounding box dimension " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      INTER = true;
      boundingBox = true;
      boundingVolume = true;
      boundingSize = atof(argv[++i]);
    }
    else if (strcmp(argv[i], "-bs") == 0 || strcmp(argv[i], "--boundingsimplex") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid bounding simplex dimension " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      // Here mark the fact we choose a simplex and not a cube.
      INTER = true;
      boundingBox = false;
      boundingVolume = true;
      boundingSize = atof(argv[++i]);
    }
    else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--dimension") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid cartesian space dimension " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      dimension = atoi(argv[++i]);
      Rn::setDimension(dimension);
    }
    else if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--tolerance") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid cartesian space tolerance " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      tolerance = atof(argv[++i]);
      Rn::setTolerance(tolerance);
    }
    else if (strcmp(argv[i], "-VO") == 0 || strcmp(argv[i], "--Volume") == 0) {
      COMPUTE_VOL = true;
      INTER = false;
    }
    else if (strcmp(argv[i], "-cg") == 0 || strcmp(argv[i], "--check-generators") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid number of generators to check " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      generatorsToCheck = atoi(argv[++i]);
      RunOptions::setCheckGeneratorsTest(true);
      RunOptions::setCheckGeneratorsValue(generatorsToCheck);
    }
    else if (strcmp(argv[i], "-cf") == 0 || strcmp(argv[i], "--check-facets") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid number of facets to check " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      facetsToCheck = atoi(argv[++i]);
      RunOptions::setCheckGeneratorsTest(true);
      RunOptions::setCheckGeneratorsValue(generatorsToCheck);
    }
    else if (strcmp(argv[i], "--generator") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid number of generator " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      vertexToComputeDistances = atoi(argv[++i]);
      computeDistances = true;
      MS = false;
      INTER = false;
      CHCK_EQ = false;
    }
    else if (strcmp(argv[i], "--facet") == 0) {
      if (i+1 == argc) {
        cerr << "Invalid number of facet " << argv[i] << std::endl;
        return EXIT_FAILURE;
      }
      facetToComputeDistances = atoi(argv[++i]);
      computeDistances = true;
      MS = false;
      INTER = false;
      CHCK_EQ = false;
    }
    else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      cout << version << std::endl;
      cout << "\t-d  [--dimension]        ARG\t:\tSet the cartesian space dimension" << std::endl;
      cout << "\t-t  [--tolerance]        ARG\t:\tSet the cartesian space tolerance [default: 1.e-06]" << std::endl;
      cout << "\t-o  [--output]           ARG\t:\tThe optional output file (ptop or pcon extension)" << std::endl;
      cout << "\t-p1 [--polytope1]        ARG\t:\tFirst  polytope input file (ptop extension)" << std::endl;
      cout << "\t-p2 [--polytope2]        ARG\t:\tSecond polytope input file (ptop extension)" << std::endl;
      cout << "\t-c1 [--polyhedralcone1]  ARG\t:\tFirst  polyhedral cone input file (pcon extension)" << std::endl;
      cout << "\t-c2 [--polyhedralcone2]  ARG\t:\tSecond polyhedral cone input file (pcon extension)" << std::endl;
      cout << "\t-cf [--check-facets]     ARG\t:\tUsed to test when we know the final number of facets" << std::endl;
      cout << "\t-cg [--check-generators] ARG\t:\tUsed to test when we know the final number of generators" << std::endl;
      cout << "\t-bs [--boundingsimplex]  ARG\t:\tBounding simplex size, containing the bounding box -bb (n+1 vertices)" << std::endl;
      cout << "\t-bb [--boundingbox]      ARG\t:\tBounding box size centered on the origin including the polytope (2^n vertices)" << std::endl;
      cout << "\t-ch [--check-all]\t\t:\tUsed to perform all tests (no arguments, can be slow)." << std::endl;
      cout << "\t-MS [--MinkowskiSum]        \t:\tSet the option to turn on Minkowski sums" << std::endl;
      cout << "\t-IN [--Intersection]        \t:\tSet the option to turn on intersections (default option)" << std::endl;
      cout << "\t-SS [--Subset]              \t:\tSet the option to test whether P1 is included in P2" << std::endl;
      cout << "\t-VO [--Volume]              \t:\tSet the option to turn on the volume computation for the polytope" << std::endl;
      cout << "\t-EQ [--Equality]            \t:\tSet the option to turn on the equality check between ptop or pcon" << std::endl;
      cout << "\t-v  [--version]\t\t\t:\tGive the current version" << std::endl;
      return 0;
    }
    else {
      cerr << "Unknown option " << argv[i] << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (Rn::getDimension() == 0) {
    cerr << "Dimension has not been set." << std::endl;
    return EXIT_FAILURE;
  }

  // Here we deal with the fact that the intersection is the default operation.
  // So if we have two operands and we don't perform neither a sum nor an equality/inclusion check then it's an intersection
  if (c2 == true || p2 == true)
    if (MS == false && CHCK_EQ == false && SUBSET == false)
      INTER = true;

  int truncationStep = 0;
  boost::shared_ptr<PolyhedralCone_Rn> A;
  boost::shared_ptr<PolyhedralCone_Rn> B;
  boost::shared_ptr<PolyhedralCone_Rn> C;
  try {
    if (p1==true) {
      A.reset(new Polytope_Rn());
      IO_Polytope::load(fileName1, A);
      // Is a bounding volume provided ?
      if (boundingVolume == true) {
        if (A->numberOfGenerators() == 0) {
          // By default the bounding volume is a cube (2n facets), if not a tetrahedron (n+1 facets).
          if (boundingBox == true) {
            truncationStep = 2*Rn::getDimension();
            A->createBoundingBox(boundingSize);
          }
          else {
            truncationStep = Rn::getDimension()+1;
            A->createBoundingSimplex(boundingSize);
          }
        }
        else {
          boost::shared_ptr<Polytope_Rn> PA = boost::static_pointer_cast<Polytope_Rn>(A);
          DoubleDescriptionFromGenerators::Compute( PA, boundingSize);
          if (check_all == true)
            PA->checkTopologyAndGeometry();
          std::string FileOut("output");
          FileOut += fileName1;
          if (output == true)
            FileOut = outputFileName;
          IO_Polytope::save(FileOut, PA);

          return EXIT_SUCCESS;
        }
      }
    }
    else if (c1==true) {
      A.reset(new PolyhedralCone_Rn());
      IO_Polytope::load(fileName1, A);
    }
    if (p2==true) {
      B.reset(new Polytope_Rn());
      IO_Polytope::load(fileName2, B);
    }
    else if (c2==true) {
      B.reset(new PolyhedralCone_Rn());
      IO_Polytope::load(fileName2, B);
    }
  }
  catch(std::ios_base::failure& e) {
    cerr << "In/out exception: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  if ((p1==true && p2==true) || (c1==true && c2==true)) {
    if (p1==true && p2==true) {
      C.reset(new Polytope_Rn());
    }
    else {
      C.reset(new PolyhedralCone_Rn());
    }
    {for (unsigned int i=0; i<A->numberOfGenerators(); i++) {
      C->addGenerator(A->getGenerator(i));
    }}
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHSA(A->getListOfHalfSpaces());
    for (iteHSA.begin(); iteHSA.end()!=true; iteHSA.next()) {
      C->addHalfSpace(iteHSA.current());
    }
    constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > iteHSB(B->getListOfHalfSpaces());
    for (iteHSB.begin(); iteHSB.end()!=true; iteHSB.next()) {
      C->addHalfSpace(iteHSB.current());
    }
    truncationStep = A->numberOfHalfSpaces();
  }
  else if ((p1==true && p2!=true) || (c1==true && c2!=true)) {
    // There is no polytope or polyhedral cone other than A.
    C = A;
  }
  else {
    cerr << "This program needs at least one polytope or one polyhedral cone and no mix between them." << std::endl;
    return EXIT_FAILURE;
  }

  boost::timer this_timer;
  try {
    // Do we have two operands and a sum or an intersection ?
    if ((MS == true) || (INTER == true) || (COMPUTE_VOL == true) || (SUBSET == true)) {
      if (MS == true) {
        // Minkowski sum calculation
        if (p1!=true || p2!=true) {
          cout << "ERROR 2 polytopes are needed to compute the minkowski sum." << std::endl;
          return -1;
        }
        //NormalFan_Rn NFA(boost::static_pointer_cast<Polytope_Rn>(A));
        //NFA.dump(cout);
        //NormalFan_Rn NFB(boost::static_pointer_cast<Polytope_Rn>(B));
        //NFB.dump(cout);
        //NormalFan_Rn NFC;
        //NFC.computeCommonRefinement(NFA,NFB);
        //NFC.dump(cout);
        C.reset(new Polytope_Rn());

        //cout << "VERTEX NUMBER = " << Polytope_Rn().computeMinkowskiVerticesNumber(A,B) << std::endl;

        boost::shared_ptr<Polytope_Rn> PA = boost::static_pointer_cast<Polytope_Rn>(A);
        boost::shared_ptr<Polytope_Rn> PB = boost::static_pointer_cast<Polytope_Rn>(B);
        boost::shared_ptr<Polytope_Rn> PC = boost::static_pointer_cast<Polytope_Rn>(C);
        MinkowskiSum Ope(PA,PB,PC);
        //return 0;
      }
      else if (COMPUTE_VOL == true) {
        boost::shared_ptr<Polytope_Rn> Pol = boost::static_pointer_cast<Polytope_Rn>(C);
        volumeOfPolytope = VolumeOfPolytopes_Rn::compute(Pol);
        cout << "VolumeOfPolytope=" << volumeOfPolytope << std::endl;
        return EXIT_SUCCESS;
      }
      else if (SUBSET == true) {
        bool isInside = A->isIncluded(B);
        std::string answer = (isInside == true) ? "" : "not";
        cout << "P1 is " << answer << " inside of P2" << std::endl;
        return EXIT_SUCCESS;
      }
      else if (INTER == false && MS == false && check_all == true && p2 == false && c2 == false) {
        // Here we just check the input body, we don't perform neither intersection nor sum.
        A->checkTopologyAndGeometry();
        return EXIT_SUCCESS;
      }
      else {
        //C->truncate(truncationStep);

        // Declare an iterator.
        //lexmaxIteratorOfListOfHalfSpaces lexmin_ite(C);
        //constIteratorOfListOfHalfSpaces lexmin_ite(C);
        constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> > lexmin_ite(C->getListOfHalfSpaces());
        //DoubleDescriptionWeakRedundancy< boost::shared_ptr<PolyhedralCone_Rn>, lexmaxIteratorOfListOfHalfSpaces > DD(C, lexmin_ite, truncationStep);
        //NoRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
        //WeakRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
        StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > NRP;
        DoubleDescription<
          boost::shared_ptr<PolyhedralCone_Rn>,
          constIteratorOfListOfGeometricObjects< boost::shared_ptr<HalfSpace_Rn> >,
          //constIteratorOfListOfHalfSpaces,
          //NoRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
          //WeakRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
          StrongRedundancyProcessing< boost::shared_ptr<PolyhedralCone_Rn> > >
            DD(C, lexmin_ite, NRP, truncationStep);
      }
      cout.precision(10 - (int)log10(Rn::getTolerance()) );
      cout << "TIME=" << this_timer.elapsed() << std::endl;

      if (C->numberOfGenerators()==0 && C->numberOfHalfSpaces()==0)
        cout << "The result is empty." << std::endl;
      else {
        if (generatorsToCheck != 0) {
          if (generatorsToCheck == C->numberOfGenerators())
            cout << "Found " << generatorsToCheck << " generators..... OK" << std::endl;
          else
            cout << "Found " << C->numberOfGenerators() << " generators instead of " << generatorsToCheck << " KO !!!!!" << std::endl;
        }
        if (facetsToCheck != 0) {
          if (facetsToCheck == C->numberOfHalfSpaces())
            cout << "Found " << facetsToCheck << "   facets  ..... OK" << std::endl;
          else
            cout << "Found " << C->numberOfHalfSpaces() << "   facets   instead of " << facetsToCheck << " KO !!!!!" << std::endl;
        }
        std::string FileOut("output");
        FileOut += fileName1;
        if (output == true)
          FileOut = outputFileName;
        IO_Polytope::save(FileOut, C);
      }
    }
  }
  catch(std::invalid_argument& excep) {
    cerr << "TIME=" << this_timer.elapsed() << std::endl;
    cerr << "Invalid argument exception " << excep.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(std::out_of_range& excep) {
    cerr << "TIME=" << this_timer.elapsed() << std::endl;
    cerr << "Out of range exception " << excep.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(std::ios_base::failure& excep) {
    cerr << "TIME=" << this_timer.elapsed() << std::endl;
    cerr << "In/out exception " << excep.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(std::domain_error& excep) {
    cerr << "TIME=" << this_timer.elapsed() << std::endl;
    cerr << "Domain error exception " << excep.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(std::logic_error& excep) {
    cerr << "TIME=" << this_timer.elapsed() << std::endl;
    cerr << "Logic error exception " << excep.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(...) {
    cerr << "TIME=" << this_timer.elapsed() << std::endl;
    cerr << "Unexpected exception caught !" << std::endl;
    return EXIT_FAILURE;
  }
  //C->dump();
  if (check_all == true)
    C->checkTopologyAndGeometry();
  else if (CHCK_EQ == true) {
    if ((p1==true && p2==true) || (c1==true && c2==true))
      A->checkEquality(B);
    else {
      cerr << "ERROR 2 polytopes or polyhedral cones are needed to check equality." << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (computeDistances == true) {
    if (vertexToComputeDistances >= 0) {
      C->checkGenerator(vertexToComputeDistances, std::cout);
    }
    else if (facetToComputeDistances >= 0) {
      C->checkFacet(facetToComputeDistances, std::cout);
    }
  }

  return EXIT_SUCCESS;
}
