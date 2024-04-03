#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

SOURCE=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas
REFIVEC=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.indicatorVectors.csv
REFSTM=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.stMat.bin
## generating structure matrix from the indicator vectors
python $srcDIR/pyKleeBarcode_computeStructureMatrix.py -i $REFIVEC -o stMat.tmp



## if everything is well, then the computed indicator vectors should be equivalent to the ones we stored
python $srcDIR/convertStructureMat_CSV_to_bin.py -i $REFSTM -o $REFSTM.csv --bintocsv
python $srcDIR/convertStructureMat_CSV_to_bin.py -i stMat.tmp -o stMat.tmp.csv --bintocsv

python $SCRIPT_DIR/compareMatrices.py -i1 $REFSTM.csv -i2 stMat.tmp.csv --header --sep ,
RV=$?

## cleaning
rm stMat.tmp
rm stMat.tmp.csv
rm $REFSTM.csv

exit $RV