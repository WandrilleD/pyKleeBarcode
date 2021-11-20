import numpy as np
from pyKleeBarcode_utils import loadSsum,saveSsum



if __name__ == "__main__":
	import sys
	import argparse

	parser = argparse.ArgumentParser(
			description="""Merges 2 Ssum matrices created from 2 different runs of pyKleeBarcode_computeSsum_MPI.py
			Warning : This script assumes that the 2 Ssum matrices where obtained from completely independent sets of sequences.
			""")
	parser.add_argument('-i1','--inputFile1', type=str, required=True,
		 help='input Ssum matrix 1 (expected: binary format as produced by the pyKleeBarcode_computeSsum_MPI.py script)')
	parser.add_argument('-i2','--inputFile2', type=str, required=True,
		 help='input Ssum matrix 2 (expected: binary format as produced by the pyKleeBarcode_computeSsum_MPI.py script)')
	parser.add_argument('-o','--outputFile', type=str, required=True,
		 help='output file name for the resulting Ssum matrix. NB: can be the same as -i1 or -i2')


	args = parser.parse_args()

	inFile1 = args.inputFile1
	inFile2 = args.inputFile2
	outFile= args.outputFile

	S_sum1 , N1 = loadSsum( inFile1 )
	S_sum2 , N2 = loadSsum( inFile2 )

	S_sum = S_sum1 + S_sum2
	N = N1 + N2
	saveSsum( S_sum , N , outFile )


