#!/bin/bash

#VolumeOfPolytope=2.55675e-05
#TIME=0.015
./politopix.exe -p1 teis1_3d.ptop -d 3 -VO

#VolumeOfPolytope=1.24114e-05
#TIME=0
./politopix.exe -p1 teis2_3d.ptop -d 3 -VO

#VolumeOfPolytope=1000
#TIME=0
./politopix.exe -p1 cube3D.ptop -d 3 -VO

#VolumeOfPolytope=10000
#TIME=0
./politopix.exe -p1 cube4D.ptop -d 4 -VO

#VolumeOfPolytope=100000
#TIME=0.015
./politopix.exe -p1 cube5D.ptop -d 5 -VO

#VolumeOfPolytope=1000000
#TIME=0.015
./politopix.exe -p1 cube6D.ptop -d 6 -VO

echo ""
echo "10 examples"
echo ""
#VolumeOfPolytope=665.111
./politopix.exe -p1 teis1_0_1-ctl_848_0_3d.ptop -d 3 -VO

#VolumeOfPolytope=1811.31
./politopix.exe -p1 teis3_100_2-ctl_573_0_3d.ptop -d 3 -VO

#VolumeOfPolytope=737.122
./politopix.exe -p1 teis4_104_0-ctl_973_3_3d.ptop -d 3 -VO

#VolumeOfPolytope=573.781
./politopix.exe -p1 teis5_107_2-ctl_186_0_3d.ptop -d 3 -VO

#VolumeOfPolytope=3087.92
./politopix.exe -p1 teis6_108_2-ctl_589_0_3d.ptop -d 3 -VO

#VolumeOfPolytope=1113.63
./politopix.exe -p1 teis7_109_1-ctl_872_3_3d.ptop -d 3 -VO

#VolumeOfPolytope=1544.25
./politopix.exe -p1 teis8_109_3-ctl_947_2_3d.ptop -d 3 -VO

#VolumeOfPolytope=607.856
./politopix.exe -p1 teis2_10_2-ctl_613_3_3d.ptop -d 3 -VO

#VolumeOfPolytope=1476.22
./politopix.exe -p1 teis9_110_0-ctl_95_0_3d.ptop -d 3 -VO

#VolumeOfPolytope=3216.24
./politopix.exe -p1 teis10_110_2-ctl_651_0_3d.ptop -d 3 -VO

#VolumeOfPolytope=0.0114731
# qhull Total volume:       0.011473107
#TIME=0.078
./politopix.exe -p1 20D6_v.ptop -d 6 -VO

#VolumeOfPolytope=0.125
# qhull Total volume:       0.125
#TIME=0
./politopix.exe -p1 BIR3-4-6_v.ptop -d 4 -VO

#VolumeOfPolytope=0.0159477
#TIME=0.421
./politopix.exe -p1 DG1.ptop -d 6 -VO
#$ qconvex s FA < qDG1.ptop
#  Approximate volume:       0.015947726

#VolumeOfPolytope=38.8411
#TIME=0.203
./politopix.exe -p1 DG2.ptop -d 6 -VO
#$ qhull s FA < qDG2.ptop
#  Approximate volume:       38.841074
