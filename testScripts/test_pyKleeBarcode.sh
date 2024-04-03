#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

SOURCE=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas
REFSTM=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.stMat.bin
## generating structure matrix from the indicator vectors
python $srcDIR/pyKleeBarcode_MPI.py -i $SOURCE -o stMat.tmp -s 4

## if everything is well, then the computed indicator vectors should be equivalent to the ones we stored

## as they are stored in the binary format, we will convert to csv:
time python $srcDIR/convertStructureMat_CSV_to_bin.py -i stMat.tmp -o stMat.tmp.csv --bintocsv
time python $srcDIR/convertStructureMat_CSV_to_bin.py -i $REFSTM -o REFstMat.tmp.csv --bintocsv

python $SCRIPT_DIR/compareMatrices.py -i1 stMat.tmp.csv -i2 REFstMat.tmp.csv --header --sep ,
RV=$?

## cleaning
rm stMat.tmp
rm stMat.tmp.csv
rm REFstMat.tmp.csv

exit $RV