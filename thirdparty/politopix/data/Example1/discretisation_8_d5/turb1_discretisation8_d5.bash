#!/bin/bash
PATH=$PATH:.
politopix -ch -p1 Polytope_1.ptop  -d 5 -bb 150 -t 1e-006
politopix -ch -p1 Polytope_2.ptop  -d 5 -bb 150 -t 1e-006
politopix -ch -p1 Polytope_3.ptop  -d 5 -bb 150 -t 1e-006
politopix -ch -p1 Polytope_4.ptop  -d 5 -bb 150 -t 1e-006
politopix -ch -p1 Polytope_5.ptop  -d 5 -bb 150 -t 1e-006
politopix -ch -p1 Polytope_6.ptop  -d 5 -bb 150 -t 1e-006
politopix -ch -p1 Polytope_7.ptop  -d 5 -bb 150 -t 1e-006
politopix -ch -p1 outputPolytope_2.ptop -p2 outputPolytope_7.ptop -o Minkowski_8.ptop  -MS -d 5 -t 0.0001
politopix -ch -p1 Minkowski_8.ptop -p2 outputPolytope_5.ptop -o Minkowski_9.ptop  -MS -d 5 -t 0.0001
politopix -ch -p1 outputPolytope_3.ptop -p2 outputPolytope_6.ptop -o Minkowski_10.ptop  -MS -d 5 -t 0.0001
politopix -ch -p1 Minkowski_9.ptop -p2 Minkowski_10.ptop -o Intersection_11.ptop  -d 5 -t 0.0001
politopix -ch -p1   Intersection_11.ptop -p2 outputPolytope_1.ptop -o Minkowski_12.ptop  -MS -d 5 -t 0.0001
politopix -ch -p1 Minkowski_12.ptop -p2 outputPolytope_4.ptop -o Minkowski_13.ptop  -MS -d 5 -t 0.0001