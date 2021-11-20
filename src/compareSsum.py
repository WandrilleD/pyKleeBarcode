import numpy as np
from pyKleeBarcode_utils import loadSsum,saveSsum



if __name__ == "__main__":
	import sys
	import argparse

	parser = argparse.ArgumentParser(
			description="""compares 2 Ssum matrices
			""")
	parser.add_argument('-i1','--inputFile1', type=str, required=True,
		 help='input Ssum matrix 1 (expected: binary format as produced by the pyKleeBarcode_computeSsum_MPI.py script)')
	parser.add_argument('-i2','--inputFile2', type=str, required=True,
		 help='input Ssum matrix 2 (expected: binary format as produced by the pyKleeBarcode_computeSsum_MPI.py script)')


	args = parser.parse_args()

	inFile1 = args.inputFile1
	inFile2 = args.inputFile2
	

	S_sum1 , N1 = loadSsum( inFile1 )
	S_sum2 , N2 = loadSsum( inFile2 )

	sameN = N1 == N2

	if not sameN :
		print("N differs between the two Ssum : {} vs. {}".format(N1,N2))


	closedness = np.isclose(S_sum1 , S_sum2 )
	sameMat = np.all( closedness )
	if not sameMat:
		s = closedness.size
		d = s - np.sum(closedness)
		print( "the two matrices differ on {} out of {} positions ({:.3f})".format( d, s , d / s ) )
	
	if sameMat and sameN:
		print( "the two Ssum matrices are equivalent" )
		exit(0)

	exit(1)



