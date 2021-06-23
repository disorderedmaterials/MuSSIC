#!/bin/sh

## To compile 
make

## copy the executable to the test folder 

cp scattcg test_systems/water/.
cp scattcg test_systems/liquid_benzene
cp scattcg test_systems/watermeth
cp scattcg test_systems/C10TAB_water

## run the tests and check for the error message

cd test_systems/water
rm *.sq
rm *.gr
./scattcg
gfortran -o test_data test_data.F
var=$(./test_data)
if [ $var -eq 1 ]
then
echo "test for water failed"
exit 1
else
echo "test for water is successful"
fi
cd ..


cd liquid_benzene
rm *.sq
rm *.gr
./scattcg
gfortran -o test_data test_data.F
var=$(./test_data)
if [ $var -eq 1 ]
then
echo "test for liquid benzene failed"
exit 1
else
echo "test for liquid benzene is successful"
fi
cd ..

cd watermeth
rm *.sq
rm *.gr
./scattcg
gfortran -o test_data test_data.F
var=$(./test_data)
if [ $var -eq 1 ]
then
echo "test for water+methanol failed"
exit 1
else
echo "test for water+methanol is successful"
fi
cd ..

cd  C10TAB_water
rm *.sq
rm *.gr
./scattcg
gfortran -o test_data test_data.F
var=$(./test_data)
if [ $var -eq 1 ]
then
echo "test for C10TAB+ water failed"
exit 1
else
echo "test for C10TAB+water is successful"
fi

cd ..
cd ..

### Tests successful, continue running the code
echo "main code running for the current input"

./scattcg






