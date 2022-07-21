#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`



sh $SCRIPT_DIR/test_computeRefMat.sh > computeRefMat.test_log.txt 2>&1
RV=$?
echo "test_computeRefMat.sh"

if [ $RV -eq 0 ]
then
  echo "computeRefMat .................................... OK"
else
  echo "Problem in computeRefMat test. Consult computeRefMat.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_computeIndicatorVector.sh > computeIndicatorVector.test_log.txt 2>&1
RV=$?
echo "test_computeIndicatorVector.sh"

if [ $RV -eq 0 ]
then
  echo "computeIndicatorVector ......................... OK"
else
  echo "Problem in computeIndicatorVector test. Consult computeIndicatorVector.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_computeStructureMatrix.sh > computeStructureMatrix.test_log.txt 2>&1
RV=$?
echo "test_computeStructureMatrix.sh"

if [ $RV -eq 0 ]
then
  echo "computeStructureMatrix ......................... OK"
else
  echo "Problem in computeStructureMatrix test. Consult computeStructureMatrix.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_convertStructureMat_CSV_to_bin.sh > convertStructureMat_CSV_to_bin.test_log.txt 2>&1
RV=$?
echo "test_convertStructureMat_CSV_to_bin.sh"

if [ $RV -eq 0 ]
then
  echo "convertStructureMat_CSV_to_bin ................. OK"
else
  echo "Problem in convertStructureMat_CSV_to_bin test. Consult convertStructureMat_CSV_to_bin.test_log.txt"
  exit $RV
fi

sh $SCRIPT_DIR/test_pyKleeBarcode.sh > pyKleeBarcode.test_log.txt 2>&1
RV=$?
echo "test_pyKleeBarcode.sh"
if [ $RV -eq 0 ]
then
  echo "pyKleeBarcode .................................. OK"
else
  echo "Problem in pyKleeBarcode test. Consult pyKleeBarcode.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_mergeRefMat.sh > mergeRefMat.test_log.txt 2>&1
RV=$?
echo "test_mergeRefMat.sh"

if [ $RV -eq 0 ]
then
  echo "mergeRefMat ...................................... OK"
else
  echo "Problem in mergeRefMat test. Consult mergeRefMat.test_log.txt"
  exit $RV
fi

sh $SCRIPT_DIR/test_updateStructureMatrix.sh > updateStructureMatrix.test_log.txt 2>&1
RV=$?
echo "test_updateStructureMatrix.sh"

if [ $RV -eq 0 ]
then
  echo "updateStructureMatrix .......................... OK"
else
  echo "Problem in updateStructureMatrix test. Consult updateStructureMatrix.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_heatmap_from_structureMatrix.sh > heatmap_from_structureMatrix.test_log.txt 2>&1
RV=$?
echo "test_heatmap_from_structureMatrix.sh"

if [ $RV -eq 0 ]
then
  echo "heatmap_from_structureMatrix ................... OK"
else
  echo "Problem in heatmap_from_structureMatrix test. Consult heatmap_from_structureMatrix.test_log.txt"
  exit $RV
fi



sh $SCRIPT_DIR/test_taxon_annotated_heatmap_from_structureMatrix.sh > taxon_annotated_heatmap_from_structureMatrix.test_log.txt 2>&1
RV=$?
echo "test_taxon_annotated_heatmap_from_structureMatrix.sh"

if [ $RV -eq 0 ]
then
  echo "taxon_annotated_heatmap_from_structureMatrix ... OK"
else
  echo "Problem in taxon_annotated_heatmap_from_structureMatrix test. Consult taxon_annotated_heatmap_from_structureMatrix.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_order_structureMatrix_with_list.sh > order_structureMatrix_with_list.test_log.txt 2>&1
RV=$?
echo "test_order_structureMatrix_with_list.sh"

if [ $RV -eq 0 ]
then
  echo "order_structureMatrix_with_list ................ OK"
else
  echo "Problem in order_structureMatrix_with_list test. Consult order_structureMatrix_with_list.test_log.txt"
  exit $RV
fi


sh $SCRIPT_DIR/test_order_structureMatrix_with_tree.sh > order_structureMatrix_with_tree.test_log.txt 2>&1
RV=$?
echo "test_order_structureMatrix_with_tree.sh"

if [ $RV -eq 0 ]
then
  echo "order_structureMatrix_with_tree ................ OK"
else
  echo "Problem in order_structureMatrix_with_tree test. Consult order_structureMatrix_with_tree.test_log.txt"
  exit $RV
fi



sh $SCRIPT_DIR/test_correspondenceFile_option.sh > correspondenceFile_option.test_log.txt 2>&1
RV=$?
echo "test_correspondenceFile_option.sh"

if [ $RV -eq 0 ]
then
  echo "correspondenceFile_option ...................... OK"
else
  echo "Problem in correspondenceFile_option test. Consult correspondenceFile_option.test_log.txt"
  exit $RV
fi

