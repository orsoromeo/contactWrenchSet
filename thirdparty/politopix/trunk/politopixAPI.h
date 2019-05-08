// politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2014-2016 : Delos Vincent
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
/// \file politopixAPI.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef INC_polAPI
#define INC_polAPI

#include "polito_Export.h"
#include "Rn.h"
#include "Voronoi_Rn.h"
#include "Polytope_Rn.h"
#include "IO_Polytope.h"
#include "Generator_Rn.h"
#include "HalfSpace_Rn.h"
#include "VolumeOfPolytopes_Rn.h"
#include "DoubleDescription_Rn.h"
#include "PolyhedralAlgorithms_Rn.h"
#include "GeometricObjectIterator_Rn.h"

using namespace std;

#define TEST_OK 0
#define TEST_KO -1

class politopixAPI {

 public:

  /// \brief Save a polytope in the corresponding file name.
  /// \param pathA The name of the current file
  /// \param A The corresponding polytope
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int savePolytope(const string& pathA,
					 boost::shared_ptr<Polytope_Rn>& A)  throw (ios_base::failure);

  /// \brief Load a polytope from the corresponding file name.
  /// \param pathA The name of the current file
  /// \param A The previously allocated polytope
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int loadPolytope(const string& pathA,
					 boost::shared_ptr<Polytope_Rn>& A)  throw (ios_base::failure);

  /// \brief Add a new half-space into the polytope data structure.
  /// \param A The current polytope
  /// \param HS The corresponding half-space
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int addHalfspace(boost::shared_ptr<Polytope_Rn>& A, const boost::shared_ptr<HalfSpace_Rn>& HS);

  /// \brief Add a new generator into the polytope data structure.
  /// \param A The current polytope
  /// \param GN The corresponding generator
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int addGenerator(boost::shared_ptr<Polytope_Rn>& A, const boost::shared_ptr<Generator_Rn>& GN);

  /// \brief Compute the HV-description for a given H-polytope or V-polytope with the Double Description algorithm.
  /// \param A The H-polytope or V-polytope
  /// \param bb_size The origin centered bounding box size providing the V-description
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int computeDoubleDescription(boost::shared_ptr<Polytope_Rn>& A, double bb_size=1000.) 
    throw (invalid_argument, out_of_range, ios_base::failure, logic_error);

  /// \brief Compute the HV-description for a given H-polytope or V-polytope with the Double Description algorithm.
  /// \param A The H-polytope or V-polytope
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int computeDoubleDescriptionWithoutCheck(boost::shared_ptr<Polytope_Rn>& A, double bb_size=1000.) 
    throw (invalid_argument, out_of_range, ios_base::failure, logic_error);

  /// \brief Compute the intersection between two HV-polytopes with the Double Description algorithm.
  /// \param A The 1st HV-polytope
  /// \param B The 2nd HV-polytope
  /// \param C The previously allocated result
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int computeIntersection(const boost::shared_ptr<Polytope_Rn>& A, const boost::shared_ptr<Polytope_Rn>& B,
      boost::shared_ptr<Polytope_Rn>& C)
    throw (invalid_argument, out_of_range, ios_base::failure, logic_error);

  /// \brief Compute the intersection between two HV-polytopes with the Double Description algorithm.
  /// \param A The 1st HV-polytope, then the result
  /// \param B The 2nd HV-polytope
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int computeIntersection(boost::shared_ptr<Polytope_Rn>& A, const boost::shared_ptr<Polytope_Rn>& B)
    throw (invalid_argument, out_of_range, ios_base::failure, logic_error);

  /// \brief Compute the intersection between two HV-polytopes with the Double Description algorithm.
  /// \param A The 1st HV-polytope, then the result
  /// \param B The 2nd HV-polytope
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int computeIntersectionWithoutCheck(boost::shared_ptr<Polytope_Rn>& A,
							   const boost::shared_ptr<Polytope_Rn>& B) 
    throw (invalid_argument, out_of_range, ios_base::failure, logic_error);

  /// \brief Test whether the polytope A V-description is inside the polytope B H-description.
  /// \param A The 1st V-polytope
  /// \param B The 2nd H-polytope
  /// \return true if \f[ A \subset B \f], false otherwise.
  static polito_EXPORT bool isIncluded(const boost::shared_ptr<Polytope_Rn>& A,
      const boost::shared_ptr<Polytope_Rn>& B)
    throw (invalid_argument, out_of_range, ios_base::failure, logic_error);

  /// \brief Compute the Minkowski sum between two HV-polytopes
  /// \param A The 1st HV-polytope
  /// \param B The 2nd HV-polytope
  /// \param C The previously allocated sum
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int computeMinkowskiSumOfPolytopes(const boost::shared_ptr<Polytope_Rn>& A,
      const boost::shared_ptr<Polytope_Rn>& B, boost::shared_ptr<Polytope_Rn>& C)
    throw (invalid_argument, out_of_range, ios_base::failure, logic_error);

  /// \brief Compute the Minkowski sum between two HV-polytopes tracing the generators
  /// \param A The 1st HV-polytope
  /// \param B The 2nd HV-polytope
  /// \param C The previously allocated sum
  /// \param traceGenerators Give for each generator of C, the list of numbers identifying its genitors in A and B
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int computeMinkowskiSumOfPolytopes(const boost::shared_ptr<Polytope_Rn>& A,
      const boost::shared_ptr<Polytope_Rn>& B,
      boost::shared_ptr<Polytope_Rn>& C,
      const std::vector< std::vector<int> >& genitorsOfGeneratorsA,
      const std::vector< std::vector<int> >& genitorsOfGeneratorsB,
      std::vector< std::vector<int> >& traceGenerators)
    throw (invalid_argument, out_of_range, ios_base::failure, logic_error);

  /// \brief Check whether two HV-polytopes are identical
  /// Check whether the vertices of A are inside B half-spaces and vice-versa. Perform also some topological verifications.
  /// \param A The 1st HV-polytope
  /// \param B The 2nd HV-polytope
  /// \param getFaceMapping If true, print the mapping between the generators and faces of both polytopes in case of equality.
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int checkEqualityOfPolytopes(const boost::shared_ptr<Polytope_Rn>& A,
      const boost::shared_ptr<Polytope_Rn>& B,
      bool getFaceMapping=false)
    throw (invalid_argument, out_of_range, ios_base::failure, logic_error);

  /// \brief Check whether two V-polytopes are identical
  /// Check whether the sets of vertices of A and B are equal.
  /// \param A The 1st V-polytope
  /// \param B The 2nd V-polytope
  /// \return true if \f[ \mathcal{V}_A = \mathcal{V}_B \f], false otherwise.
  static polito_EXPORT bool checkEqualityOfVertices(const boost::shared_ptr<Polytope_Rn>& A,
      const boost::shared_ptr<Polytope_Rn>& B)
    throw (invalid_argument, out_of_range, ios_base::failure, logic_error);

  /// \brief Check whether a HV-polytopes is correct
  /// \param A A HV-polytope
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int checkTopologyAndGeometry(const boost::shared_ptr<PolyhedralCone_Rn>& A);

  /// \brief Create a cube whose vertices will be (+-M, ..., +-M)
  /// \param A The cube variable
  /// \param M The cube half side length.
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int makeCube(boost::shared_ptr<Polytope_Rn>& A, double M);

  /// \brief Compute the polar polytope
  /// \param original_pol The input polytope
  /// \param polar_pol The polar polytope
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int PolarPolytope(const boost::shared_ptr<Polytope_Rn>& original_pol, boost::shared_ptr<Polytope_Rn>& polar_pol);

  /// \brief Translate a polytope or polyhedral cone by the given vector.
  /// \param pol The corresponding polytope
  /// \param v2t The translation vector
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int Translate(boost::shared_ptr<Polytope_Rn>& pol, const boost::numeric::ublas::vector<double>& v2t);

  /// \brief Return the volume of the given polytope P with its double description.
  /// The implemented algorithm can be found in <i> Volume Computation for Polytopes: Strategies and Performances </i>
  /// by <b>Andreas Enge</b> in <i> Encyclopedia of Optimization </i> 2nd edition, p 4032-4073.
  static polito_EXPORT double computeVolume(const boost::shared_ptr<Polytope_Rn> P)
    throw (invalid_argument, out_of_range, ios_base::failure);

  /// \brief Remove all cap half-spaces and then compute the intersection of two capped polytopes
  /// \param A The 1st HV-polytope
  /// \param B The 2nd HV-polytope
  /// \param C The previously allocated intersection
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int pseudoIntersection(
     const boost::shared_ptr<Polytope_Rn>& A,
     const boost::shared_ptr<Polytope_Rn>& B,
     boost::shared_ptr<Polytope_Rn>& C,
     const std::set< unsigned int >& firstOperandCaps,
     const std::set< unsigned int >& secondOperandCaps,
     std::set< unsigned int >& newCaps,
     double bb_size=1000.) throw (invalid_argument,out_of_range,ios_base::failure,logic_error);

  /// \brief Compute the Minkowski sum of two polytopes and then remove all cap half-spaces to truncate again.
  /// \param A The 1st HV-polytope
  /// \param B The 2nd HV-polytope
  /// \param C The previously allocated sum
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int pseudoSum(
     const boost::shared_ptr<Polytope_Rn>& A,
     const boost::shared_ptr<Polytope_Rn>& B,
     boost::shared_ptr<Polytope_Rn>& C,
     const std::set< unsigned int >& firstOperandCaps,
     const std::set< unsigned int >& secondOperandCaps,
     std::set< unsigned int >& newCaps,
     double bb_size=1000.) throw (invalid_argument,out_of_range,ios_base::failure,logic_error);

  /// \brief Compute the Voronoi Diagram in a n-dimensional space i.e. a partitioning of an input space into regions based on distance to points called seeds.
  /// \param inputSpace The input HV-polytope, most of the times a parallelepiped
  /// \param listOfSeeds The list of points to be considered as seeds
  /// \param VoronoiCells The list of the returned HV-polytopes partitioning the input space
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static polito_EXPORT int computeVoronoiDiagram(
     const boost::shared_ptr<Polytope_Rn>& inputSpace,
     const std::vector<Point_Rn>& listOfSeeds,
     std::vector< boost::shared_ptr<Polytope_Rn> >& VoronoiCells) throw (std::length_error);


};

#endif // INC_polAPI
