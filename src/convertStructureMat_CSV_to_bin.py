import numpy as np
from pyKleeBarcode_utils import writeStructureMatrix_line_bin,yieldStructureMatrix_bin




### testing code
# H = ['a','b','c']
# M = np.array([[1.0, 0., 0.],
#        [0.3945758 , 1.0, 0. ],
#        [0.89605101, 0.74844958, 1.0]])
# writeStructureMatrix_bin( H , M , 'test.bin' )
# i = 0
# for h,x in yieldStructureMatrix_bin( 'test.bin' ):
# 	assert H[i] == h
# 	assert np.allclose( M[i,:i] , x )
# 	print(h,x)
# 	i+=1

## reverse conversion.
def bin_to_csv_structure_matrix( fileBin , fileCSV ):
	"""
	This function reads a file containing a structure matrix in binary format and write it in csv format.
	
	Implementation note:
		structure matrices can get quite large in memory, (e.g. 3.1Gb for a 20 000 by 20 000 matrix),
		so we make the choice to have a "slow" algorithm which reads the binary file multiple time, but have some control over the used RAM.

	Takes:
		- fileBin : str . path to the structure matrix in binary format
		- fileCSV : str . path to the structure matrix in CSV format
	"""


	with open( fileCSV , 'w' ) as OUT:

		nb_species = 0 

		with open( fileBin , 'rb' ) as IN:
			for h,x in yieldStructureMatrix_bin( IN ):
				if nb_species != 0:
					print( ',' , end='', file = OUT )
				print(h , end='', file = OUT)

				nb_species+=1
			print('', file = OUT)

		for i in range(nb_species): # we read the entire bin file each time
			j = 0
			with open( fileBin , 'rb' ) as IN:
				for h,x in yieldStructureMatrix_bin( IN ): 
					#print(i,j,h,len(x))
					if j<i: # nothing : we a before the sequence of interest
						j+=1
						continue
					elif i == j: # sequence of interest : the entire line is interesting to us
						print( ','.join([str(v) for v in x]) , end='', file = OUT)
						if j != 0:
							print( ',' , end='', file = OUT )
						print( '1.0' , end='', file = OUT)
					else: # after the sequence of interest, we want the element of the column corresponding to the sequence of interest
						print( ',' + str(x[i]) , end='', file = OUT )

					j+=1
			
			print('', file = OUT)
			

	return		






if __name__ == "__main__":
	import sys
	import argparse

	parser = argparse.ArgumentParser(
			description="""converts a csv structure matrix to a binary one""")
	parser.add_argument('-i','--inputFile', type=str, required=True,
		 help='input csv structure matrix')
	parser.add_argument('-o','--outputFile', type=str, required=True,
		 help='output file name for the binary structure matrix')
	parser.add_argument('--bintocsv', action='store_true', required=False, 
		help='converts from binary to csv format')


	args = parser.parse_args()

	inFile = args.inputFile
	outFile = args.outputFile
	
	if args.bintocsv:
		bin_to_csv_structure_matrix( inFile , outFile )
	else:
		with open( inFile , 'r' ) as IN , open(outFile, 'wb+') as OUT:

			##getting the sequence names
			header = IN.readline().strip().split(',')

			for i,l in enumerate( IN ):

				sl = l.strip().split(',')
				arr = np.array( [float(x) for x in sl[:i]] )
				writeStructureMatrix_line_bin( header[i] , arr, OUT )



	exit(0)





