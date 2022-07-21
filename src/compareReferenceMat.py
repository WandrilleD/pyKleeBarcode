import numpy as np
from pyKleeBarcode_utils import loadRefM,saveRefM



if __name__ == "__main__":
	import sys
	import argparse

	parser = argparse.ArgumentParser(
			description="""compares 2 RefM matrices
			""")
	parser.add_argument('-i1','--inputFile1', type=str, required=True,
		 help='input reference matrix 1 (expected: binary format as produced by the pyKleeBarcode_computeRefMat_MPI.py script)')
	parser.add_argument('-i2','--inputFile2', type=str, required=True,
		 help='input reference matrix 2 (expected: binary format as produced by the pyKleeBarcode_computeRefMat_MPI.py script)')


	args = parser.parse_args()

	inFile1 = args.inputFile1
	inFile2 = args.inputFile2
	

	S_sum1 , N1 = loadRefM( inFile1 )
	S_sum2 , N2 = loadRefM( inFile2 )

	sameN = N1 == N2

	if not sameN :
		print("N differs between the two reference matrix : {} vs. {}".format(N1,N2))


	closedness = np.isclose(S_sum1 , S_sum2 )
	sameMat = np.all( closedness )
	if not sameMat:
		s = closedness.size
		d = s - np.sum(closedness)
		print( "the two matrices differ on {} out of {} positions ({:.3f})".format( d, s , d / s ) )
	
	if sameMat and sameN:
		print( "the two reference matrices are equivalent" )
		exit(0)

	exit(1)



