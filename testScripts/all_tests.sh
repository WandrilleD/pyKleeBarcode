#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`



sh $SCRIPT_DIR/test_computeSsum.sh > computeSsum.test_log.txt 2>&1
RV=$?
echo "test_computeSsum.sh"

if [ $RV -eq 0 ]
then
  echo "computeSsum ...................... OK"
else
  echo "Problem in computeSsum test. Consult computeSsum.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_computeIndicatorVector.sh > computeIndicatorVector.test_log.txt 2>&1
RV=$?
echo "test_computeIndicatorVector.sh"

if [ $RV -eq 0 ]
then
  echo "computeIndicatorVector ........... OK"
else
  echo "Problem in computeIndicatorVector test. Consult computeIndicatorVector.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_computeStructureMatrix.sh > computeStructureMatrix.test_log.txt 2>&1
RV=$?
echo "test_computeStructureMatrix.sh"

if [ $RV -eq 0 ]
then
  echo "computeStructureMatrix ........... OK"
else
  echo "Problem in computeStructureMatrix test. Consult computeStructureMatrix.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_convertStructureMat_CSV_to_bin.sh > convertStructureMat_CSV_to_bin.test_log.txt 2>&1
RV=$?
echo "test_convertStructureMat_CSV_to_bin.sh"

if [ $RV -eq 0 ]
then
  echo "convertStructureMat_CSV_to_bin ... OK"
else
  echo "Problem in convertStructureMat_CSV_to_bin test. Consult convertStructureMat_CSV_to_bin.test_log.txt"
  exit $RV
fi

sh $SCRIPT_DIR/test_pyKleeBarcode.sh > pyKleeBarcode.test_log.txt 2>&1
RV=$?
echo "test_pyKleeBarcode.sh"
if [ $RV -eq 0 ]
then
  echo "pyKleeBarcode .................... OK"
else
  echo "Problem in pyKleeBarcode test. Consult pyKleeBarcode.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_mergeSsum.sh > mergeSsum.test_log.txt 2>&1
RV=$?
echo "test_mergeSsum.sh"

if [ $RV -eq 0 ]
then
  echo "mergeSsum ........................ OK"
else
  echo "Problem in mergeSsum test. Consult mergeSsum.test_log.txt"
  exit $RV
fi

sh $SCRIPT_DIR/test_updateStructureMatrix.sh > updateStructureMatrix.test_log.txt 2>&1
RV=$?
echo "test_updateStructureMatrix.sh"

if [ $RV -eq 0 ]
then
  echo "updateStructureMatrix ............ OK"
else
  echo "Problem in updateStructureMatrix test. Consult updateStructureMatrix.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_heatmap_from_structureMatrix.sh > heatmap_from_structureMatrix.test_log.txt 2>&1
RV=$?
echo "test_heatmap_from_structureMatrix.sh"

if [ $RV -eq 0 ]
then
  echo "heatmap_from_structureMatrix ..... OK"
else
  echo "Problem in heatmap_from_structureMatrix test. Consult heatmap_from_structureMatrix.test_log.txt"
  exit $RV
fi

