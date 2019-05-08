/*  This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
 
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
 
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
 
#include <assert.h>
#include "CommandlineOptions.hh"
#include "Array.hh"
#include "VPolytope.hh"
#include "VPolytopeList.hh"
#include "IntegerMatrix.hh"
#include "LPinter.hh"
#include "RSinter.hh"

VPolytopeList polylist;

int main(const int argc, const char** argv) {
  CommandlineOptions::init(argc, argv);

  dd_set_global_constants();  // First, this must be called to use cdd

  // Read the polytopes

  if (!polylist.read(std::cin)) {
    std::cerr << "error while reading polytopes" << std::endl;
    return 1;
  }

  if (!CommandlineOptions::mute()) {
    std::cerr << "Sum of " << polylist.numPol() << " Polytopes" << std::endl;
  }

  polylist.incMinkSumSetup();
  
  if (CommandlineOptions::shortoutput()) {
    std::cout << "[";
    bool first = true;
    while (polylist.hasMoreVertices()) {
      if (!first)
	std::cout << "," << std::endl;
      std::cout << polylist.incExploreStepFast().coord();
      first = false;
    }
    std::cout << "]" << std::endl;
  } else {
    while (polylist.hasMoreVertices()) {
      if (CommandlineOptions::longoutput()) {
	std::cout << polylist.incExploreStep();
      } else {
	std::cout << polylist.incExploreStepFast();
      }
    }
  }

  dd_free_global_constants();  // Nice
  return 0;
}
