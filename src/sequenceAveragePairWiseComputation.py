from collections import defaultdict

from pyKleeBarcode_utils import iterFasta



def getSequenceSpecies( fastaId , sep , index ):
	""" 
		Takes :
			- fastaId (str) : id of the sequence
			- sep (str) : separator between fields in the fasta id
			- index (str) : index of the species name in the fasta id (starts at 0)

		Returns:
			(str) : species name
	"""
	return fastaId.split(sep)[index] # brutal

def getPairWiseDiff( seq1 , seq2 , countGapAsDiff=True):
	'''
		Takes:
			- seq1 (str) : a sequence
			- seq2 (str) : a sequence
			- countGapAsDiff (bool) (default=True) : if False, gaps do not add to the total differences
		
		Return:
			(float) : number of differences 
			OR
			 -1 : if the sequences do not have the same length


	>>> getPairWiseDiff( "a" , "at" , countGapAsDiff=True)
	-1
	>>> getPairWiseDiff( "ttg-" , "atgt" , countGapAsDiff=True)
	2
	>>> getPairWiseDiff( "ttg-" , "atgt" , countGapAsDiff=False)
	1

	'''
	if len(seq1) != len(seq2):
		return -1

	l = len(seq1)
	i=0
	diffs = 0
	while i < l:

		if seq1[i] != seq2[i]:
			if ( not countGapAsDiff ) and ( seq1[i]=='-' or seq2[i]=='-' ) :
				##we avoid counting this as a difference
				pass
			else:
				## difference
				diffs +=1
		i+=1

	return diffs


def getAveragePairWiseDiff( seqList1 , seqList2 , sameList = False , countGapAsDiff=True):
	'''
		Takes:
			- seqList1 (list) : list of sequences
			- seqList2 (list) : list of sequences
			- sameList (bool) (default=False) : True if the two lists are actually the same 
												in order to avoid counting the comparison of a sequence to itself
												or counting the differences twice
			- countGapAsDiff (bool) (default=True) : if False, gaps do not add to the total differences
		
		Return:
			(float) : average pairwise differences 
			OR
			 -1 : if some sequences do not have the same length
			OR
			 "NA" : if sameList is True and the lists contain a single sequence

	>>> getAveragePairWiseDiff( ["at"] , ["at"] , sameList = False , countGapAsDiff=True)
	0.0
	>>> getAveragePairWiseDiff( ["at","ag"] , ["at","t-"] , sameList = False , countGapAsDiff=True)
	1.25
	>>> getAveragePairWiseDiff( ["at","ag"] , ["at","t-"] , sameList = False , countGapAsDiff=False)
	0.75
	>>> getAveragePairWiseDiff( ["at","at","ag","ag","ag"] , ["at","at","ag","ag","ag"] , sameList = True , countGapAsDiff=True)
	0.6

	'''

	diffSum = 0
	nbComparisons = 0

	for i,s1 in enumerate(seqList1):
		for j,s2 in enumerate(seqList2):
			if sameList and i<=j:
				continue ## same sequence : ignore
			else:
				PD = getPairWiseDiff( s1 , s2 , countGapAsDiff )
				if PD == -1:
					## there was a length difference
					return -1
				else:
					diffSum += PD
					nbComparisons += 1
	#print(diffSum,nbComparisons)
	if nbComparisons == 0:
		return "NA"
	return diffSum/nbComparisons



if __name__ == "__main__":

	## testing
	import doctest
	doctest.testmod()
	## end testing


	import sys
	import argparse

	parser = argparse.ArgumentParser(
				description="""Computes an intra-species and inter-species sequence similarity matrix
								using average pairwise difference among individuals (often termed Pi in population genetics).
								""")
	parser.add_argument('-i','--inputFile', type=str, required=True,
			 help='input multiple sequence alignment in fasta format (preferably, trimmed)')
	parser.add_argument('-o','--outputFile', type=str, required=True,
			 help='output file name for the structural matrix')

	#parser.add_argument('-m','--max-seq-per-species', type=int, default=3,
	#		 help='maximum number of sequences kept per species. Default : 3. Set to 0 to have no limit ; be careful thought, as this parameter has a direct impact on speed.')

	parser.add_argument('-f','--field-delimitor', type=str, default="|",
			 help='field delimitor of the fasta sequence id lines')

	parser.add_argument('-s','--species-field-index', type=int, default=1,
			 help='index (starting at 0) of the species name in the fasta sequence id lines')


	args = parser.parse_args()

	inFile = args.inputFile
	outFile= args.outputFile


	#to be replaced by something better in time
	FASTA_ID_SEPARATOR= args.field_delimitor # separator between fields in the fasta id
	SPECIES_NAME_INDEX= args.species_field_index  # index of the species name in the fasta id (starts at 0)
	
	
	
	
	speciesSequences = defaultdict(list) # store ids of sequences, organized by species
	sequences = dict() # keys : fastaID , value : sequence
	speciesOrder = []
	seqLen = 0
	
	## reading data
	
	IN = open(inFile, 'r')
	for k,s in iterFasta(IN):
		sp = getSequenceSpecies(k, FASTA_ID_SEPARATOR , SPECIES_NAME_INDEX )
		
		speciesSequences[ sp ].append( s )
		seqLen = len(s)
		if len(speciesSequences[ sp ])==1:
			speciesOrder.append(sp) # keeping the order of species somewhere
	
	IN.close()
	

	sumWithin = 0
	nbWithin = 0
	sumBetween = 0
	nbBetween = 0

	APDmatrix = defaultdict( dict )

	for i , sp1 in enumerate( speciesOrder  ):
		for j , sp2 in enumerate( speciesOrder  ):
			if i<=j:
				APD = getAveragePairWiseDiff( speciesSequences[sp1] , speciesSequences[sp2] , sameList = (i==j)  , countGapAsDiff = True )
				if APD == -1:
					print( "Sequences length differ between species",sp1,"and",sp2,".","Make sure the sequences are aligned!" )
					exit(1)
				elif APD == "NA":
					APDmatrix[ sp1 ][ sp2 ] = APD ## same species and only 1 sequence ...
				else:
					APDmatrix[ sp1 ][ sp2 ] = APD / seqLen # normalizing by sequence length

					if i==j:
						sumWithin+=APDmatrix[ sp1 ][ sp2 ]
						nbWithin +=1
					else:
						sumBetween+=APDmatrix[ sp1 ][ sp2 ]
						nbBetween +=1

	## writing results
	OUT = open( outFile , "w" )
	print("sp1;sp2;n1;n2;APD" , file=OUT )
	for sp1 in APDmatrix.keys():
		for sp2 , apd in APDmatrix[sp1].items():
			print(  sp1 , sp2 , len(speciesSequences[sp1]) , len(speciesSequences[sp2]) , apd , file=OUT , sep = ";" )

	OUT.close()

	print( "average pairwise difference:" )
	print( "within species " , sumWithin / nbWithin )
	print( "between species" , sumBetween / nbBetween )

