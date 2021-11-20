import numpy as np



if __name__ == "__main__":
	import sys
	import argparse

	parser = argparse.ArgumentParser(
			description="""compares 2 csv matrices
			""")
	parser.add_argument('-i1','--inputFile1', type=str, required=True,
		 help='input csv matrix 1')
	parser.add_argument('-i2','--inputFile2', type=str, required=True,
		 help='input csv matrix 2')
	parser.add_argument('--header', action='store_true', required=False,
		help='assume the matrices have a header line and compare this as well')
	parser.add_argument('-s','--sep', type=str, default = ' ',
		help='separator between fields (default is whitespace)')

	args = parser.parse_args()

	inFile1 = args.inputFile1
	inFile2 = args.inputFile2
	

	hline1 = ''
	hline2 = ''

	mat1 = None
	mat1 = None

	with open(inFile1, 'r') as IN:
		if args.header:
			hline1 = IN.readline()
		mat1 = np.loadtxt( IN , delimiter = args.sep)

	with open(inFile2, 'r') as IN:
		if args.header:
			hline2 = IN.readline()
		mat2 = np.loadtxt( IN , delimiter = args.sep)
		


	sameH = hline1 == hline1

	if not sameH :
		print("the header differs between the two matrices")


	closedness = np.isclose(mat1 , mat2 )
	sameMat = np.all( closedness )
	if not sameMat:
		s = closedness.size
		d = s - np.sum(closedness)
		print( "the two matrices differ on {} out of {} positions ({:.3f})".format( d, s , d / s ) )
	
	if sameMat and sameH:
		print( "the two matrices are equivalent" )
		exit(0)

	exit(1)



