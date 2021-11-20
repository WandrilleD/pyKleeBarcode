#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

REFSTM=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.stMat.csv
## generating structure matrix from the indicator vectors
echo "CSV to bin"
time python $srcDIR/convertStructureMat_CSV_to_bin.py -i $REFSTM -o stMat.tmp.bin
echo "bin to CSV"
time python $srcDIR/convertStructureMat_CSV_to_bin.py -i stMat.tmp.bin -o stMat.tmp.csv --bintocsv

## if everything is well, then the computed indicator vectors should be equivalent to the ones we stored
python $SCRIPT_DIR/compareMatrices.py -i1 $REFSTM -i2 stMat.tmp.csv --header --sep ,
RV=$?

## cleaning
rm stMat.tmp.bin
rm stMat.tmp.csv

exit $RV