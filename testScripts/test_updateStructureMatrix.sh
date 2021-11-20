#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

SOURCE=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas
REFIVEC=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.indicatorVectors.csv

## generating 3 sets of indicator vectors:
head -n 10 $REFIVEC > iVec10.tmp.csv
head -n 20 $REFIVEC > iVec20.tmp.csv
tail -n 10 iVec20.tmp.csv > iVec1020.tmp.csv

## generating structure matrix from the different indicator vector sets
python $srcDIR/pyKleeBarcode_computeStructureMatrix.py -i iVec10.tmp.csv -o stMat10.tmp
python $srcDIR/pyKleeBarcode_computeStructureMatrix.py -i iVec20.tmp.csv -o stMat20.tmp

python $srcDIR/pyKleeBarcode_updateStructureMatrix.py -r iVec10.tmp.csv -n iVec1020.tmp.csv -m stMat10.tmp

## if everything is well, then the computed structure matrix vectors should be equivalent
## they will not be equal bit by bit, but up to (at least) a 10**-8 factor
python $srcDIR/convertStructureMat_CSV_to_bin.py -i stMat10.tmp -o stMat10.tmp.csv --bintocsv
python $srcDIR/convertStructureMat_CSV_to_bin.py -i stMat20.tmp -o stMat20.tmp.csv --bintocsv

python $SCRIPT_DIR/compareMatrices.py -i1 stMat10.tmp.csv -i2 stMat20.tmp.csv --header --sep ,
RV=$?

## cleaning
rm iVec10.tmp.csv
rm iVec20.tmp.csv
rm iVec1020.tmp.csv
rm stMat10.tmp.csv
rm stMat20.tmp.csv
rm stMat10.tmp
rm stMat20.tmp

exit $RV