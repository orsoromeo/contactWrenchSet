#!/bin/bash

echo "/// COMPUTE VERTICES FROM HALFSPACES ///" > MS_test_file.txt
file1=
file2=
boundingBox=
dimension=
tolerance=
fctnb=
vtxnb=

computeVertices() {
    echo "/// " $file1
    echo "///" >> MS_test_file.txt
    echo "/// " $file1 >> MS_test_file.txt
    echo "///" >> MS_test_file.txt
    if [ -e $file1 ]
    then
	./politopix -p1 $file1 -bb $boundingBox -d $dimension --check-generators $vtxnb -ch >> MS_test_file.txt
	#valgrind --leak-check=yes ./politopix -p1 $file1 -bb $boundingBox -d $dimension --check-generators $vtxnb -ch >> MS_test_file.txt
    else
	echo "### Missing file : " $file1 " ###" >> MS_test_file.txt
	echo "Missing file : " $file1
    fi
}

computePolytopeIntersection() {
    echo "////// " $file1 " " $file2
    echo "//////" >> MS_test_file.txt
    echo "////// " $file1 " " $file2 >> MS_test_file.txt
    echo "//////" >> MS_test_file.txt
    if [ -e $file1 ] && [ -e $file2 ]
    then
	./politopix -p1 $file1 -p2 $file2 -d $dimension -t $tolerance --check-generators $vtxnb -ch >> MS_test_file.txt
	#valgrind --leak-check=yes ./politopix -p1 $file1 -p2 $file2 -d $dimension -t $tolerance --check-generators $vtxnb -ch >> MS_test_file.txt
    else
	echo "### Missing file : " $file1 " or " $file2 " ###" >> MS_test_file.txt
	echo "Missing file : " $file1 " or " $file2
    fi
}

computePolyhedralConeIntersection() {
    echo "////// " $file1 " " $file2
    echo "//////" >> MS_test_file.txt
    echo "////// " $file1 " " $file2 >> MS_test_file.txt
    echo "//////" >> MS_test_file.txt
    if [ -e $file1 ] && [ -e $file2 ]
    then
	./politopix -c1 $file1 -c2 $file2 -d $dimension -t $tolerance --check-generators $vtxnb -ch >> MS_test_file.txt
	#valgrind --leak-check=yes ./politopix -c1 $file1 -c2 $file2 -d $dimension -t $tolerance --check-generators $vtxnb -ch >> MS_test_file.txt
    else
	echo "### Missing file : " $file1 " or " $file2 " ###" >> MS_test_file.txt
	echo "Missing file : " $file1 " or " $file2
    fi
}

computeMinkowskiSum() {
    echo "////// " $file1 " " $file2
    echo "//////" >> MS_test_file.txt
    echo "////// " $file1 " " $file2 >> MS_test_file.txt
    echo "//////" >> MS_test_file.txt
    if [ -e $file1 ] && [ -e $file2 ]
    then
	./politopix -p1 $file1 -p2 $file2 -d $dimension -t $tolerance --check-generators $vtxnb --check-facets $fctnb -MS -ch >> MS_test_file.txt
	#valgrind --leak-check=yes ./politopix -p1 $file1 -p2 $file2 -d $dimension -t $tolerance --check-generators $vtxnb --check-facets $fctnb -MS -ch >> MS_test_file.txt
    else
	echo "### Missing file : " $file1 " or " $file2 " ###" >> MS_test_file.txt
	echo "Missing file : " $file1 " or " $file2
    fi
}

testPolytopeEquality() {
    echo "////// " $file1 " " $file2
    echo "//////" >> MS_test_file.txt
    echo "////// " $file1 " " $file2 >> MS_test_file.txt
    echo "//////" >> MS_test_file.txt
    if [ -e $file1 ] && [ -e $file2 ]
    then
	./politopix -p1 $file1 -p2 $file2 -d $dimension -t $tolerance -EQ >> MS_test_file.txt
    else
	echo "### Missing file : " $file1 " or " $file2 " ###" >> MS_test_file.txt
	echo "Missing file : " $file1 " or " $file2
    fi
}

####################
# Compute vertices #
####################
#N_VERTICES 8
file1=cyclic.ptop
boundingBox=10000
dimension=4
vtxnb=8
computeVertices

#N_VERTICES 6
file1=lmp1.ptop
boundingBox=100
dimension=3
vtxnb=6
computeVertices

#N_VERTICES 6
file1=lmp2.ptop
boundingBox=100
dimension=3
vtxnb=6
computeVertices

#N_VERTICES 12
file1=lmp3.ptop
boundingBox=100
dimension=3
vtxnb=12
computeVertices

#N_VERTICES 8
file1=SharirCube.ptop
boundingBox=25
dimension=3
vtxnb=8
computeVertices

#N_VERTICES 8
file1=neighborly_4_8.ptop
boundingBox=10
dimension=4
vtxnb=8
computeVertices

#N_VERTICES 6
file1=BIR3-4-6.ptop
boundingBox=2
dimension=4
vtxnb=6
computeVertices

#N_VERTICES 4
file1=CUT3-3-4.ptop
boundingBox=2
dimension=3
vtxnb=4
computeVertices

#N_VERTICES 7
file1=p4_7_simplicial_1.ptop
boundingBox=2
dimension=4
vtxnb=7
computeVertices

#N_VERTICES 7
file1=p4_7_nonsimplicial_01.ptop
boundingBox=1
dimension=4
vtxnb=7
computeVertices

#N_VERTICES 6
file1=CR0-3-6.ptop
boundingBox=2
dimension=3
vtxnb=6
computeVertices

#N_VERTICES 600
file1=120-cell.ptop
boundingBox=10
dimension=4
vtxnb=600
computeVertices

#N_VERTICES 32
file1=not-cubical4.ptop
boundingBox=20
dimension=4
vtxnb=32
computeVertices

#N_VERTICES 11
file1=CF-10-11.ptop
boundingBox=2
dimension=10
#computeVertices

#N_VERTICES 8
file1=8D6.ptop
boundingBox=100.
dimension=6
vtxnb=8
computeVertices

#N_VERTICES 10
file1=10D6.ptop
boundingBox=100.
dimension=6
vtxnb=10
computeVertices

#N_VERTICES 20
file1=20D6.ptop
boundingBox=100.
dimension=6
vtxnb=20
computeVertices

###########################
# Test polytopes equality #
###########################
file1=output8D6.ptop
file2=8D6_v.ptop
dimension=6
tolerance=0.000001
testPolytopeEquality

file1=output10D6.ptop
file2=10D6_v.ptop
dimension=6
tolerance=0.000001
testPolytopeEquality

file1=output20D6.ptop
file2=20D6_v.ptop
dimension=6
tolerance=0.000001
testPolytopeEquality

###############################
# Intersect polytopes & cones #
###############################
file1=lmp1_v.ptop
file2=SharirCube_v.ptop
dimension=3
tolerance=0.0001
vtxnb=8
computePolytopeIntersection

file1=lmp1.pcon
file2=SharirCube.pcon
dimension=4
tolerance=0.0001
vtxnb=8
computePolyhedralConeIntersection

file1=SharirCube.pcon
file2=lmp1.pcon
dimension=4
tolerance=0.0001
vtxnb=8
computePolyhedralConeIntersection

#########################
# Compute Minkowski Sum #
#########################
file1=lmp1_v.ptop
file2=cube3D.ptop
dimension=3
tolerance=0.0001
fctnb=26
vtxnb=24
computeMinkowskiSum

file1=cube3D.ptop
file2=lmp1_v.ptop
dimension=3
tolerance=0.0001
fctnb=26
vtxnb=24
computeMinkowskiSum

file1=tetra3D_v_1.ptop
file2=tetra3D_v_2.ptop
dimension=3
tolerance=0.000001
fctnb=7
vtxnb=9
computeMinkowskiSum

file1=tetra3D_v_2.ptop
file2=tetra3D_v_1.ptop
dimension=3
tolerance=0.000001
fctnb=7
vtxnb=9
computeMinkowskiSum

file1=tetra3D_v_3.ptop
file2=tetra3D_v_4.ptop
dimension=3
tolerance=0.000001
fctnb=17
vtxnb=17
computeMinkowskiSum

file1=tetra3D_v_4.ptop
file2=tetra3D_v_3.ptop
dimension=3
tolerance=0.000001
fctnb=17
vtxnb=17
computeMinkowskiSum

file1=_p4s.ptop
file2=_p4ns.ptop
dimension=4
tolerance=0.0001
fctnb=33
vtxnb=25
computeMinkowskiSum

file1=_p4ns.ptop
file2=_p4s.ptop
dimension=4
tolerance=0.0001
fctnb=33
vtxnb=25
computeMinkowskiSum

file1=cyclic_v.ptop
file2=120-cell_v.ptop
dimension=4
tolerance=0.0001
fctnb=824
vtxnb=1467
#computeMinkowskiSum

file1=120-cell_v.ptop
file2=cyclic_v.ptop
dimension=4
tolerance=0.0001
fctnb=824
vtxnb=1467
#computeMinkowskiSum

file1=neighborly_4_8_v.ptop
file2=BIR3-4-6_v.ptop
dimension=4
tolerance=0.000001
fctnb=40
vtxnb=34
computeMinkowskiSum

file1=BIR3-4-6_v.ptop
file2=neighborly_4_8_v.ptop
dimension=4
tolerance=0.000001
fctnb=40
vtxnb=34
computeMinkowskiSum


# Do not leave here all the output files.
mv output* ./txt-output
# Print problems
grep KO MS_test_file.txt
grep exception MS_test_file.txt
