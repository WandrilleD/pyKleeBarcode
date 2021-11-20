from collections import defaultdict
import numpy as np
from scipy.sparse import csr_matrix
import time

from pyKleeBarcode_utils import yieldIndicatorVectors, writeStructureMatrix_line_bin

#from pyKleeBarcode_linearAlgebra import computeStructureMatrixFromIndicatorVectors, computeAndWriteStructureMatrixFromIndicatorVectors



def computeStructureMatrixSegment( indicVectorDict , speciesOrder , seqLen , outHandle,header=None ):
	"""
		Takes:
			- indicVectorDict (dict) : keys are species/group names, values are the corresponding indicator vector (numpy.array)
			- speciesOrder (list) : desired order of the indicVectorDict keys
			- seqLen (int) : length of the studied sequence (4* nucleotide size)
			- outHandle (file) : file to output to
			- header=None : optional, for a list of species name to wrtie as a header

		Returns:
			(None)
	"""

	nbSpecies = len(speciesOrder)

	V = np.zeros( ( seqLen , nbSpecies ) )

	for i,sp in enumerate(speciesOrder):
		V[:,i] =  indicVectorDict[sp]

	V=-V
	VT = V.T

	for i in range(nbSpecies):
		X =  VT[i,:] @ V[:,:i]

		writeStructureMatrix_line_bin( header[i] , np.abs(X[:i]) , outHandle )

	# if not header is None:
	# 	print(','.join(header),file=outHandle)

	# for i in range(nbSpecies):

	# 	X =  VT[i,:] @ V 
	# 	print(','.join([str(x) for x in abs(X) ]),file=outHandle)

	return


if __name__ == "__main__":
	import sys

	
	import argparse
	parser = argparse.ArgumentParser(
			description="""Updates a structure matrix from a set of reference indicator vectors (the sequences already in the structure matrix) amd a set of new structure vectors (to add to the structure matrix)
						   according to the definitions of "A scalable method for analysis and display of DNA sequences" by Sirovitch et alii .

						   IMPORTANT: this script assumes sequences are ordered similarly between the reference indicator vectors and structure matrix.
							""")
	parser.add_argument('-r','--reference-vectors', type=str, required=True,
		 help='input reference indicator vectors in csv format (one seq/species per line)')
	parser.add_argument('-n','--new-vectors', type=str, required=True,
		 help='input new indicator vectors in csv format (one seq/species per line)')
	parser.add_argument('-m','--structure-matrix', type=str, required=True,
		 help='structure matrix file to update')
	args = parser.parse_args()
	
	newFile = args.new_vectors
	refFile = args.reference_vectors
	outFile = args.structure_matrix




	new_indic_vectors , new_speciesOrder = None , []
	seqLen = 0

	T = time.time()	
	for s , arr in yieldIndicatorVectors( newFile ):
		new_speciesOrder.append(s)

		if new_indic_vectors is None:
			seqLen = len(arr)
			new_indic_vectors = arr.reshape( (seqLen,1) )
		else:
			new_indic_vectors = np.append(new_indic_vectors ,  arr.reshape( (seqLen,1) ) , axis = 1)
	
	new_indic_vectors = - new_indic_vectors

	print("new indicator vectors read ({:.2f}s)".format(time.time() -T))
	

	nbSpecies_new = len(new_speciesOrder)


	strMat_update = None

	T = time.time()	
	# going through reference vectors 1 by 1
	with open(refFile,'r') as IN:

		for l in IN:
			sl = l.strip().split(',')
			
			ref_arr = - np.array( [float(x) for x in sl[1:]] ).reshape(1,seqLen)

			## compute the correlation between this reference vector and the new ones 
			corr_arr = np.abs( (ref_arr @ new_indic_vectors).T )

			if strMat_update is None:
				strMat_update = corr_arr
			else:
				strMat_update = np.append( strMat_update , corr_arr , axis = 1 )
	
	nb_refSeq = strMat_update.shape[1]
	# now we want to compute the correlation of the new vectors between themselves:
	InnerV = np.abs( new_indic_vectors.T @ new_indic_vectors )
	strMat_update = np.append( strMat_update , InnerV , axis = 1 )

	print("correlation with and between new indicator vectors computed ({:.2f}s)".format(time.time() -T))

	T = time.time()
		
	OUT=open(outFile,'ab')

	for i,s in enumerate(new_speciesOrder):
		writeStructureMatrix_line_bin( s , strMat_update[i,:(nb_refSeq+i)] , OUT )
	
	OUT.close()

	print("Structure matrix updated in {} ({:.2f}s)".format(outFile , time.time()-T))
	