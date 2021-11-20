from collections import defaultdict
import numpy as np
from scipy.sparse import csr_matrix
import time

from pyKleeBarcode_utils import readIndicatorVectors, writeStructureMatrix_bin

from pyKleeBarcode_linearAlgebra import computeStructureMatrixFromIndicatorVectors, computeAndWriteStructureMatrixFromIndicatorVectors



if __name__ == "__main__":
	import sys

	
	import argparse
	parser = argparse.ArgumentParser(
			description="""Computes a structure matrix from a set of indicator vectors
						   according to the definitions of "A scalable method for analysis and display of DNA sequences" by Sirovitch et alii 
							""")
	parser.add_argument('-i','--inputFile', type=str, required=True,
		 help='input indicator vectors in csv format (one seq/species per line)')
	parser.add_argument('-o','--outputFile', type=str, required=True,
		 help='output file name for the structure matrix')
	args = parser.parse_args()
	inFile = args.inputFile
	outFile= args.outputFile






	timerDict = {}

	T = time.time()	
	global_indic_vectors , speciesOrder = readIndicatorVectors(args.inputFile)

	seqLen = global_indic_vectors[speciesOrder[0]].shape[0]
	nbSpecies = len(speciesOrder)

	timerDict["reading indicator vectors"] = time.time() -T
	print("indicator vectors read")

	
	## pushing all vectors in a single matrix and computing the structure matrix from them
	## at this point only the master process is working
	

	T = time.time()
	if nbSpecies <= 5000 : 

		structureMatrix = computeStructureMatrixFromIndicatorVectors( global_indic_vectors , speciesOrder , seqLen )

		timerDict["total structure matrix computation"] = time.time() - T

		print("Structure matrix computed")

		## reporting
		T = time.time()
		### old code writing the structure matrix as a csv file
		#OUT=open(outFile,'w')
		#CSVmat(structureMatrix,OUT,header=speciesOrder)
		#OUT.close()
		###
		writeStructureMatrix_bin( speciesOrder , structureMatrix , outFile )
		print("Structure matrix written in",outFile)
		timerDict["writing"] = time.time() - T
	else:
		
		OUT=open(outFile,'wb+')
		computeAndWriteStructureMatrixFromIndicatorVectors( global_indic_vectors , speciesOrder , seqLen , OUT,header=speciesOrder )
		OUT.close()
		print("Structure matrix written in",outFile)
		timerDict["structure matrix computation and writing"] = time.time() - T



	timeLines = ''
	for k in timerDict.keys():
		timeLines +=" ".join([ k , str(timerDict[k]) , 's']) + '\n'

	print(timeLines,end="")



        
