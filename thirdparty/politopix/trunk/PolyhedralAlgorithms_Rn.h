// politopix allows to make computations on polytopes such as finding vertices, intersecting, Minkowski sums, ...
//     Copyright (C) 2015-2016 : Delos Vincent
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
/// \file PolyhedralAlgorithms.h
/// \author Delos Vincent (v.delos@i2m.u-bordeaux1.fr)
// I2M (UMR CNRS 5295 / University of Bordeaux)

#ifndef INC_POLY_ALGO
#define INC_POLY_ALGO

#include <boost/numeric/ublas/vector.hpp>
#include "polito_Export.h"
#include "Rn.h"
#include "Polytope_Rn.h"
#include "IO_Polytope.h"
#include "Generator_Rn.h"
#include "NormalFan_Rn.h"
#include "HalfSpace_Rn.h"
#include "PolyhedralCone_Rn.h"
#include "VolumeOfPolytopes_Rn.h"
#include "DoubleDescription_Rn.h"
#include "PolyhedralAlgorithms_Rn.h"
#include "GeometricObjectIterator_Rn.h"

using namespace std;

typedef std::vector< unsigned int > ListOfFaces;

/// \class FaceEnumeration
/// \brief Combinatorial face enumeration for polytopes.
class polito_EXPORT FaceEnumeration {

 public:

  /// \brief General Face Enumeration Algorithm 
  FaceEnumeration(const boost::shared_ptr<Polytope_Rn>& A):_polytope(A) {}

  /// General Face Enumeration Algorithm from
  /// <i> Combinatorial face enumeration in convex polytopes </i> (1994) by <b>Komei Fukuda</b> and <b>Vera Rosta</b>. <br>
  ///  Input: the set \f$ \mathcal{P}_0 \f$ of vertices of a polytope <i>P</i>. <br>
  /// Output: the set \f$ \mathcal{P} \f$ of vertices of a polytope <i>P</i>. <br>
  /// <b>procedure</b> FaceEnumeration( \f$ \mathcal{P}_0 \f$: vertices) <br>
  /// Create a binary tree <i>T</i> with set of leaves \f$ \mathcal{P}_0 \f$ <br>
  /// k=0; <i>f<sub>0</sub></i> = \f$ | \mathcal{P}_0 | \f$; \f$ \mathcal{P}_{(0)} = \mathcal{P}_0 \f$ <br>
  /// <b> WHILE </b> <i> f<sub>k</sub> >= 2 </i> <b> DO </b> <br>
  /// &nbsp; &nbsp; <i>f<sub>k+1</sub> = 0 </i>; \f$ \mathcal{P}_{(k+1)} = \emptyset \f$ <br>
  /// &nbsp; &nbsp; <b> FOREACH </b> pair  <i> (F, F') </i> in \f$ \mathcal{P}_{(k)} \f$ <b> DO </b> <br>
  /// &nbsp; &nbsp; &nbsp; &nbsp; \f$ F'' = F \cap F' \f$ <br>
  /// &nbsp; &nbsp; &nbsp; &nbsp; <b> IF </b> \f$ F'' \notin T \f$ <b> THEN </b> <br>
  /// &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; <b> IF </b> \f$ F''== F \f$ or \f$ F''== F' \f$ <b> THEN </b> <br>
  /// &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Delete <i> F'' </i> from \f$ \mathcal{P}_{(k)} \f$; <br>
  /// &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; <i>f<sub>k</sub> = <i>f<sub>k</sub> - 1  <br>
  /// &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; <b> ENDIF </b> <br>
  /// &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Add <i> F'' </i> to <i>T</i> and to \f$ \mathcal{P}_{(k+1)} \f$ <br>
  /// &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; <i>f<sub>k+1</sub> = f<sub>k+1</sub> + 1</i> <br>
  /// &nbsp; &nbsp; &nbsp; &nbsp; <b> ENDIF </b> <br>
  /// &nbsp; &nbsp; <b> ENDFOR </b> <br>
  /// &nbsp; &nbsp; Output \f$ \mathcal{P}_{(k)} \f$; <i> k=k+1 </i> <br>
  /// <b> ENDWHILE </b> <br>
  static void Compute(const boost::shared_ptr<Polytope_Rn>& A);

  static void Compute(const boost::shared_ptr<Polytope_Rn>& A, FaceEnumeration& FE);

  const std::vector< std::vector< ListOfFaces > >& getFacesWithVertices() throw (std::domain_error) {
    if (_allFacesWithVertices.empty() == false)
      return _allFacesWithVertices;
    throw std::domain_error("FaceEnumeration::getFacesWithVertices() empty DS allFacesWithVertices.");
  }

  const std::vector< std::vector< ListOfFaces > >& getFacesWithFacets() throw (std::domain_error) {
    if (_allFacesWithFacets.empty() == false)
      return _allFacesWithFacets;
    throw std::domain_error("FaceEnumeration::getFacesWithFacets() empty DS allFacesWithFacets.");
  }

  void clear() {_allFacesWithFacets.clear(); _allFacesWithVertices.clear();}

  void printFacesWithVerticesToSage(std::ostream &this_ostream) const;

  void printFacesWithVertices(std::ostream &this_ostream) const;

  void printFacesWithFacets(std::ostream &this_ostream) const;

  /// \brief Save the polytope lattice
  static void save(const std::string& filename, const std::vector< std::vector< ListOfFaces > >& latt);

  /// \brief Load the polytope lattice
  static void load(const std::string& filename, std::vector< std::vector< ListOfFaces > >& latt)
    throw (std::ios_base::failure);

  /// \brief Load the polytope lattice
  /// 1st  line : comments = "# SpaceDimension NumberOfHalfspaces" <br>
  /// 2nd  line : SpaceDimension TotalNumberOfFaces <br>
  /// 3rd  line : FiDimension k V1 ... Vk <br>
  /// 4th  line : FjDimension l Vu ... V(u+l) <br>
  /// k-th line : ... <br>
  static void load(std::istream& this_stream, std::vector< std::vector< ListOfFaces > >& latt)
    throw (std::out_of_range);

  /// \brief Save the polytope lattice
  /// 1st  line : comments = "# SpaceDimension NumberOfHalfspaces" <br>
  /// 2nd  line : SpaceDimension TotalNumberOfFaces <br>
  /// 3rd  line : FiDimension k V1 ... Vk <br>
  /// 4th  line : FjDimension l Vu ... V(u+l) <br>
  /// k-th line : ... <br>
  static void save(std::ostream& this_stream, const std::vector< std::vector< ListOfFaces > >& latt);

  void save(std::ostream& this_stream) const { save(this_stream, _allFacesWithVertices); }


 protected:
  static void ComputeWithFacets(const boost::shared_ptr<Polytope_Rn>& A, FaceEnumeration& FaceEnum);

  static void ComputeWithVertices(const boost::shared_ptr<Polytope_Rn>& A, FaceEnumeration& FaceEnum);

  std::vector< std::vector< ListOfFaces > > _allFacesWithFacets;

  std::vector< std::vector< ListOfFaces > > _allFacesWithVertices;

  const boost::shared_ptr<Polytope_Rn>& _polytope;

};


/// \class MinkowskiSum
/// \brief Compute the Minkowski sum of two polytopes.
class polito_EXPORT MinkowskiSum {

 public:

  /// \brief Compute the Minkowski sum of two polytopes.
  MinkowskiSum(
      const boost::shared_ptr<Polytope_Rn>& A,
      const boost::shared_ptr<Polytope_Rn>& B,
      boost::shared_ptr<Polytope_Rn>& C):_firstOperand(A),_secondOperand(B),_sum(C)
  { C = compute(); }

  /// \brief Compute the Minkowski sum of two polytopes.
  MinkowskiSum(
      const boost::shared_ptr<Polytope_Rn>& A,
      const boost::shared_ptr<Polytope_Rn>& B,
      boost::shared_ptr<Polytope_Rn>& C,
      const std::vector< std::vector<int> >& genitorsOfGeneratorsA,
      const std::vector< std::vector<int> >& genitorsOfGeneratorsB,
      std::vector< std::vector<int> >& traceGenerators):_firstOperand(A),_secondOperand(B),_sum(C)
  {
    C = compute();
    traceGenerators.clear();
    //traceGenerators.resize(C->numberOfGenerators());
    unsigned int countOK=0;
    {for (unsigned int i=0; i<_MinkowskiDecomposition.size(); ++i) {
      if (_MinkowskiDecompositionOK[i] == true) {
	// vertex_i = vertex_a + vertex_b
	unsigned int a = _MinkowskiDecomposition[i].first;
	unsigned int b = _MinkowskiDecomposition[i].second;
	std::vector<int> genitorsAB;
	genitorsAB.insert(genitorsAB.end(), genitorsOfGeneratorsA[a].begin(), genitorsOfGeneratorsA[a].end());
	genitorsAB.insert(genitorsAB.end(), genitorsOfGeneratorsB[b].begin(), genitorsOfGeneratorsB[b].end());
	traceGenerators.push_back(genitorsAB);
	++countOK;
      }
    }}
  }

  /// \brief Remove the cap half-spaces stored in sets and then truncate again.
  /// \return The new sum
  boost::shared_ptr<Polytope_Rn> rebuildSum(
      const std::set< unsigned int >& firstOperandCaps,
      const std::set< unsigned int >& secondOperandCaps,
      std::set< unsigned int >& newCaps,
      double bb_size=1000.);

 protected:
  boost::shared_ptr<Polytope_Rn> compute() throw (std::domain_error);

  /// Do the final job after having intersected all dual cones. The reduction process simply compares all dual cones generators.
  void processNormalFan0();

  /// Do the final job after having intersected all dual cones. The reduction process uses neighbourhood properties to identify dual cones generators.
  void processNormalFan1();

  /// Do the final job after having intersected all dual cones. The reduction process builds half-spaces and identifies them with they generators lists.
  void processNormalFan2();

  /// \brief Return the cap half-spaces of the sum in function of the two operands cap half-spaces.
  void computeCapHalfSpaces(
      const std::set< unsigned int >& firstOperandCaps,
      const std::set< unsigned int >& secondOperandCaps,
      std::set< unsigned int >& sumCaps) const throw (std::domain_error);

 protected:
  const boost::shared_ptr<Polytope_Rn> _firstOperand;
  const boost::shared_ptr<Polytope_Rn> _secondOperand;
  boost::shared_ptr<Polytope_Rn> _sum;

  // Give the relation between a vertex of A and all of its associated Minkowski vertices i.e. its polyhedrical cap
  // A2C[i] = [ c_u, c_v, ... ]  }
  //                             } => c_u = a_i + b_j
  // B2C[j] = [ c_u, c_w, ... ]  }
  /// \brief Store the polyhedrical cap in C of each vertex of A
  std::vector< std::vector<unsigned int> > _A2C;
  /// \brief Store the polyhedrical cap in C of each vertex of B
  std::vector< std::vector<unsigned int> > _B2C;

  /// neighboursA(a_i,a_j)==1 <=> a_i R a_j
  /// \brief Store the neighbours of each vertex of A
  std::vector< std::vector<unsigned int> > _neighboursA;
  /// neighboursB(b_i,b_j)==1 <=> b_u R b_v
  /// \brief Store the neighbours of each vertex of B
  std::vector< std::vector<unsigned int> > _neighboursB;

  // minkowskiVertices(i,j)==1 <=> a_i + b_j is a Minkowski vertex.

  //  MinkowskiDecomposition[k] = (a_i,b_j) <=> c_k = a_i + b_j is a Minkowski vertex.
  // To store each Minkowski vertex decomposition
  /// \brief Store the genitors in A and B of each vertex of C
  std::vector< std::pair< unsigned int, unsigned int > > _MinkowskiDecomposition;
  /// Tell whether _MinkowskiDecomposition had to be considered or not.
  std::vector< bool > _MinkowskiDecompositionOK;

  /// The normal fan polyhedrical cones list
  std::vector< boost::shared_ptr<PolyhedralCone_Rn> > _NF_Cones;
  /// The list of C vertices
  std::vector< boost::shared_ptr<Generator_Rn> >      _NF_Vertices;

};


/// \class PseudoSumWithoutCaps
/// \brief Compute the Minkowski sum of two polytopes and then remove all cap half-spaces to truncate again.
class polito_EXPORT PseudoSumWithoutCaps : public MinkowskiSum {

public:

  /// \brief Compute the Minkowski sum of two polytopes and then remove all cap half-spaces to truncate again.
  PseudoSumWithoutCaps(
     const boost::shared_ptr<Polytope_Rn>& A,
     const boost::shared_ptr<Polytope_Rn>& B,
     boost::shared_ptr<Polytope_Rn>& C,
     const std::set< unsigned int >& firstOperandCaps,
     const std::set< unsigned int >& secondOperandCaps,
     std::set< unsigned int >& newCaps,
     double bb_size=1000.):MinkowskiSum(A,B,C)
   { C = rebuildSum(firstOperandCaps, secondOperandCaps, newCaps, bb_size); }


protected:

  /// \brief Remove the cap half-spaces stored in sets and then truncate again.
  /// \return The new sum
  boost::shared_ptr<Polytope_Rn> rebuildSum(
      const std::set< unsigned int >& firstOperandCaps,
      const std::set< unsigned int >& secondOperandCaps,
      std::set< unsigned int >& newCaps,
      double bb_size=1000.);

  /// \brief Return the cap half-spaces of the sum in function of the two operands cap half-spaces.
  void computeCapHalfSpaces(
      const std::set< unsigned int >& firstOperandCaps,
      const std::set< unsigned int >& secondOperandCaps,
      std::set< unsigned int >& sumCaps) throw (std::domain_error);

};


/// \class PseudoIntersectionWithoutCaps
/// \brief Remove all cap half-spaces and then compute the intersection of two capped polytopes
class polito_EXPORT PseudoIntersectionWithoutCaps {

public:

  /// \brief Remove all cap half-spaces and then compute the intersection of two capped polytopes
  PseudoIntersectionWithoutCaps(
     const boost::shared_ptr<Polytope_Rn>& A,
     const boost::shared_ptr<Polytope_Rn>& B,
     boost::shared_ptr<Polytope_Rn>& C,
     const std::set< unsigned int >& firstOperandCaps,
     const std::set< unsigned int >& secondOperandCaps,
     std::set< unsigned int >& newCaps,
     double bb_size=1000.);

};


/// \class TopGeomTools
/// \brief Basic tools for topology and geometry: translations, polarity, ...
class polito_EXPORT TopGeomTools {

public:

  /// \brief Translate a polytope or polyhedral cone by the given vector.
  /// \param pol The corresponding polytope
  /// \param v2t The translation vector
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static int Translate(boost::shared_ptr<Polytope_Rn>& pol, const boost::numeric::ublas::vector<double>& v2t);

  /// \brief Translate a polytope or polyhedral cone by the given vector.
  /// \param pol The corresponding polytope
  /// \param v2t The translation vector
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static int GravityCenter(boost::shared_ptr<Polytope_Rn>& pol, boost::numeric::ublas::vector<double>& gravity_center);

  /// \brief Compute the polar polytope
  /// \param original_pol The input polytope
  /// \param polar_pol The polar polytope
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static int PolarPolytope(
      const boost::shared_ptr<Polytope_Rn>& original_pol,
      boost::shared_ptr<Polytope_Rn>& polar_pol,
      bool forceComputation = true, double bb_size=1000.) throw (invalid_argument);

  /// \brief Compute the projection of a polytope on the intersection of canonical hyperplanes of the shape <i>x<sub>i</sub> = 0 </i>
  /// \param listOfHyperplanes The set of indices describing the canonical hyperplanes \f$ i \in listOfHyperplanes \Leftrightarrow x_i = 0 \f$
  /// \param original_pol The input polytope
  /// \param proj_pol The projected polytope
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static int projectPolytopeOnCanonicalHyperplanes(
      const std::set< unsigned int >& listOfHyperplanes,
      const boost::shared_ptr<Polytope_Rn>& original_pol,
      boost::shared_ptr<Polytope_Rn>& proj_pol) throw (invalid_argument);

};


/// \class DoubleDescriptionFromGenerators
/// \brief Compute the V-description from the H-description
class polito_EXPORT DoubleDescriptionFromGenerators {

public:

  /// \brief Use the polarity to get the facets from the generators
  /// \param pol The input polytope
  /// \return TEST_OK or 0 if the process was successful, TEST_KO or -1 if something went wrong.
  static int Compute(boost::shared_ptr<Polytope_Rn>& pol, double bb_size=1000.)
    throw (invalid_argument, out_of_range, ios_base::failure, logic_error);

};


/// \class Visualization
/// \brief 2D visualization tools
class polito_EXPORT Visualization {

public:

  /// \brief Provide the drwing of polygon under the gnuplot format
  /// \param pol The input 2D polygon
  /// \param name The polygon unique name
  /// \param col A number in [0,1] coding the color
  static void gnuplot2D(const boost::shared_ptr<Polytope_Rn>& polygon, const std::string& name, double col, std::ostream& out) throw (std::domain_error);

};

#endif // INC_POLY_ALGO
