#!/bin/bash
#D4=../data/Turbine1/discretisation_4
#D5=../data/Turbine1/discretisation_5
#D6=../data/Turbine1/discretisation_6
#D8_d5=../data/Turbine1/discretisation_8_d5
#D8_d6=../data/Turbine1/discretisation_8_d6

#echo "Process discretisation_4"
#cd $D4
#./turb1_discretisation.bash > ./T4.txt
#echo "Process discretisation_5"
#cd ../discretisation_5
#./turb1_discretisation5.bash > ./T5.txt
#echo "Process discretisation_6"
#cd ../discretisation_6
#./turb1_discretisation6.bash > ./T6.txt
#echo "Process discretisation_8_d5"
#cd ../discretisation_8_d5
#./turb1_discretisation8_d5.bash > ./T8_d5.txt
#echo "Process discretisation_8_d6"
#cd ../discretisation_8_d6
#./turb1_discretisation8_d6.bash > ./T8_d6.txt

if [ ! -f politopix ]
then
echo "The executable politopix does not exist."
exit -1
fi

current_path=`pwd`

first_test=$current_path
first_test+=/../data
cp politopix $first_test
cd $first_test
echo "Work in" $first_test
./_test_MS.sh

INTERFD3R6=$current_path
INTERFD3R6+=/../../tolgeompub/TM/Interference_bride/DTL_BPTL/D3-R6/
cp politopix $INTERFD3R6
cd $INTERFD3R6
echo "Work in" $INTERFD3R6
# test whether a file beginning by output exists
if ls output* &> /dev/null
then
rm output*
fi
# test whether a file beginning by Intersection exists
if ls Intersection* &> /dev/null
then
rm Intersection*
fi
chmod +x ./brides.sh
./brides.sh > INTERFD3R6.txt
grep KO INTERFD3R6.txt
grep exception INTERFD3R6.txt

I5966=$current_path
I5966+=/../../tolgeompub/TM/Interference_bride/Bride_D_5.966/
cp politopix $I5966
cd $I5966
echo "Work in" $I5966
# test whether a file beginning by output exists
if ls output* &> /dev/null
then
rm output*
fi
# test whether a file beginning by Intersection exists
if ls Intersection* &> /dev/null
then
rm Intersection*
fi
chmod +x ./Bride.bash
./Bride.bash > I5966.txt
grep KO I5966.txt
grep exception I5966.txt

I5976=$current_path
I5976+=/../../tolgeompub/TM/Interference_bride/Bride_D_5.976/
cp politopix $I5976
cd $I5976
echo "Work in" $I5976
if ls output* &> /dev/null
then
rm output* 
fi
if ls Intersection* &> /dev/null
then
rm Intersection*
fi
chmod +x ./Bride.bash
./Bride.bash > I5976.txt
grep KO I5976.txt
grep exception I5976.txt

exit

I5986=$current_path
I5986+=/../../tolgeompub/TM/Interference_bride/Bride_D_5.986/
cp politopix $I5986
cd $I5986
echo "Work in" $I5986
if ls output* &> /dev/null
then
rm output* 
fi
if ls Intersection* &> /dev/null
then
rm Intersection*
fi
chmod +x ./Bride.bash
./Bride.bash > I5986.txt
grep KO I5986.txt
grep exception I5986.txt

I5996=$current_path
I5996+=/../../tolgeompub/TM/Interference_bride/Bride_D_5.996/
cp politopix $I5996
cd $I5996
echo "Work in" $I5996
if ls output* &> /dev/null
then
rm output* 
fi
if ls Intersection* &> /dev/null
then
rm Intersection*
fi
chmod +x ./Bride.bash
./Bride.bash > I5996.txt
grep KO I5996.txt
grep exception I5996.txt

I6006=$current_path
I6006+=/../../tolgeompub/TM/Interference_bride/Bride_D_6.006/
cp politopix $I6006
cd $I6006
echo "Work in" $I6006
if ls output* &> /dev/null
then
rm output* 
fi
if ls Intersection* &> /dev/null
then
rm Intersection*
fi
chmod +x ./Bride.bash
./Bride.bash > I6006.txt
grep KO I6006.txt
grep exception I6006.txt

I6016=$current_path
I6016+=/../../tolgeompub/TM/Interference_bride/Bride_D_6.016/
cp politopix $I6016
cd $I6016
echo "Work in" $I6016
if ls output* &> /dev/null
then
rm output* 
fi
if ls Intersection* &> /dev/null
then
rm Intersection*
fi
chmod +x ./Bride.bash
./Bride.bash > I6016.txt
grep KO I6016.txt
grep exception I6016.txt

I6026=$current_path
I6026+=/../../tolgeompub/TM/Interference_bride/Bride_D_6.026/
cp politopix $I6026
cd $I6026
echo "Work in" $I6026
if ls output* &> /dev/null
then
rm output* 
fi
if ls Intersection* &> /dev/null
then
rm Intersection*
fi
chmod +x ./Bride.bash
./Bride.bash > I6026.txt
grep KO I6026.txt
grep exception I6026.txt
