#!/bin/bash

No_Pro=$(awk -F "=" '/no_of_cores/ {print $2}' abe.ini)
No_Pro_cm=$(awk -F "=" '/no_of_cores_cm/ {print $2}' abe.ini)
OUTPUT="RunStatus"

if [ -d "$OUTPUT" ]; then
rm -r "$OUTPUT"
fi
mkdir "$OUTPUT"


echo "Running initial settings"
python initial.py &>> "$OUTPUT/initial.txt"
if [ $? != 0 ]; then
echo " Error while running initial.py; check Runstatus"
exit 1
fi

echo "Running CatalogMaker"
mpiexec -n $No_Pro_cm python CatalogMaker_MPI.py &>> "$OUTPUT/CatalogMaker.txt"
if [ $? != 0 ]; then
echo " Error while running CatalogMaker_MPI.py; check Runstatus"
exit 1
fi

echo "Running Data Combainer"
python Combainer.py &>> "$OUTPUT/Combainer.txt"
if [ $? != 0 ]; then
echo " Error while running Combianer.py; check Runstatus"
exit 1
fi

echo "Running CatalogAnalyser"
mpiexec -n $No_Pro python CatalogAnalyser_MPI.py &>> "$OUTPUT/CatalogAnalyser.txt"
if [ $? != 0 ]; then
echo " Error while running CatalogAnalyser_MPI.py; check Runstatus"
exit 1
fi

echo "Running Data Editor"
python DataEditor.py &>> "$OUTPUT/DataEditor.txt"
if [ $? != 0 ]; then
echo " Error while running DataEditor.py; check Runstatus"
exit 1
fi

echo "Running Data Analyser"
mpiexec -n $No_Pro python DataAnalysis_MPI.py &>> "$OUTPUT/DataAnalyser.txt"
if [ $? != 0 ]; then
echo " Error while running DataAnalysis_MPI.py; check Runstatus"
exit 1
fi 


echo "Running Luminosity dump combainer"
python Combainer.py lum &>> "$OUTPUT/luminosity.txt"
if [ $? != 0 ]; then
echo " Error while running Combianer.py(luminosity); check Runstatus"
exit 1
fi
