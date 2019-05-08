#!/bin/bash
PATH=$PATH:.
echo "Caracteristiques du jeu d'essai :"
echo "dimension : 6"
echo "Polytope_1, Polytope_2, Polytope_4, Polytope_5, Polytope_7 :"
echo "16 demi-espaces, dimension intrinseque : 4"
echo "Polytope_3, Polytope_6 :"
echo "8 demi-espaces, dimension intrinseque : 3"
echo ""
echo "Intersection des demi-espaces des operandes"
politopix -ch -p1 Polytope_1.ptop  -d 6 -bb 150 -t 1e-004 -cg 64
politopix -p1 outputPolytope_1.ptop -p2 pol1.ptop -d 6 -t 1e-004 -EQ

politopix -ch -p1 Polytope_2.ptop  -d 6 -bb 150 -t 1e-004 -cg 64
politopix -p1 outputPolytope_2.ptop -p2 pol2.ptop -d 6 -t 1e-004 -EQ

politopix -ch -p1 Polytope_3.ptop  -d 6 -bb 150 -t 1e-004 -cg 48
politopix -p1 outputPolytope_3.ptop -p2 pol3.ptop -d 6 -t 1e-004 -EQ

politopix -ch -p1 Polytope_4.ptop  -d 6 -bb 150 -t 1e-004 -cg 64
politopix -p1 outputPolytope_4.ptop -p2 pol4.ptop -d 6 -t 1e-004 -EQ

politopix -ch -p1 Polytope_5.ptop  -d 6 -bb 150 -t 1e-004 -cg 64
politopix -p1 outputPolytope_5.ptop -p2 pol5.ptop -d 6 -t 1e-004 -EQ

politopix -ch -p1 Polytope_6.ptop  -d 6 -bb 150 -t 1e-004 -cg 48
politopix -p1 outputPolytope_6.ptop -p2 pol6.ptop -d 6 -t 1e-004 -EQ

politopix -ch -p1 Polytope_7.ptop  -d 6 -bb 150 -t 1e-004 -cg 64
politopix -p1 outputPolytope_7.ptop -p2 pol7.ptop -d 6 -t 1e-004 -EQ

echo "Minkowski_8"
politopix -ch -p1 outputPolytope_2.ptop -p2 outputPolytope_7.ptop -o Minkowski_8.ptop  -MS -d 6 -t 0.0001
echo "Minkowski_9"
politopix -ch -p1 Minkowski_8.ptop -p2 outputPolytope_5.ptop -o Minkowski_9.ptop  -MS -d 6 -t 0.0001
echo "Minkowski_10"
politopix -ch -p1 outputPolytope_3.ptop -p2 outputPolytope_6.ptop -o Minkowski_10.ptop  -MS -d 6 -t 0.0001
echo "Intersection_11"
politopix -ch -p1 Minkowski_9.ptop -p2 Minkowski_10.ptop -o Intersection_11.ptop  -d 6 -t 0.0001
echo "Minkowski_12"
politopix -ch -p1   Intersection_11.ptop -p2 outputPolytope_1.ptop -o Minkowski_12.ptop  -MS -d 6 -t 0.0001
echo "Minkowski_13"
politopix -ch -p1 Minkowski_12.ptop -p2 outputPolytope_4.ptop -o Minkowski_13.ptop  -MS -d 6 -t 0.0001
