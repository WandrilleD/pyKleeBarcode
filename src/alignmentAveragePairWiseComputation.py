from collections import defaultdict

from pyKleeBarcode_utils import iterFasta

IUPAC = {'a': ['a'],
		't': ['t'],
		'g': ['g'],
		'c': ['c'],
		'r': ['a','g'],
		'y': ['c','t'],
		's': ['g','c'],
		'w': ['a','t'],
		'k': ['g','t'],
		'm': ['a','c'],
		'b': ['t','g','c'],
		'd': ['a','t','g'],
		'h': ['a','t','c'],
		'v': ['a','g','c'],
		'n': ['a','t','g','c'],
		'-': ['-']}

DiffDict = {}

for l1 in IUPAC.keys():
	for l2 in IUPAC.keys():
		different = True
		for n in IUPAC[l1]:
			if n in IUPAC[l2]:
				different=False
				break
		DiffDict[l1+l2] = different



def getPairWiseDiff( seq1 , seq2 ):
	'''
		Takes:
			- seq1 (str) : a sequence
			- seq2 (str) : a sequence
			- countGapAsDiff (bool) (default=True) : if False, gaps do not add to the total differences
			- countNAsDiff=False (bool) (default=False) : if True, Ns do not add to the total differences
		
		Return:
			(float) : number of differences 
			OR
			 -1 : if the sequences do not have the same length

	'''
	if len(seq1) != len(seq2):
		return -1

	l = len(seq1)
	i=0
	diffs = 0
	while i < l:

		if isDiff(seq1[i],seq2[i]):
			diffs +=1
		i+=1

	return diffs

def countNumberComparablePositions( seqs ):
	"""
	Counts the number of positions in the list of aligned sequences that contains something else than only N or only gaps

		Takes:
			- seqs (list): list of sequences
	"""

	L = len(seqs[0])
	comparable=0
	for i in range(L):
		ref = seqs[0][i]
		if not ref in ['n','-']:
			comparable += 1
		else: #reference is either n or -
			for s in seqs:
				if s[i] != ref:
					comparable += 1
					break
	return comparable

def countNumberVariablePositions( seqs ):
	"""
	Counts the number of positions in the list of aligned sequences that contains something else than only N or only gaps

		Takes:
			- seqs (list): list of sequences
	"""

	L = len(seqs[0])
	variable=0
	for i in range(L):
		ref = seqs[0][i]
		if ref=='n': # reference is an n : this is undesirable
			j = 1
			while j < len(seqs):
				ref = seqs[j][i]
				if ref !='n':
					break
				j+=1
			if j == len(seqs):
				continue #only Ns


		for s in seqs:
			if isDiff(s[i],ref):
				variable += 1
				break

	return variable

def isDiff(c1,c2):
	return DiffDict[c1+c2]

if __name__ == "__main__":

	## testing
	import doctest
	doctest.testmod()
	## end testing


	import sys
	import argparse

	parser = argparse.ArgumentParser(
				description=""" computes average pairwise difference among sequences in an alignment (often termed Pi in population genetics) as well as other diversity metrics. """)
	parser.add_argument('-i','--inputFile', type=str, required=True,
			 help='input multiple sequence alignment in fasta format (preferably, trimmed)')
	parser.add_argument('-o','--outputFile', type=str, required=True,
			 help='output file name for all pairwise differences')



	args = parser.parse_args()

	inFile = args.inputFile
	outFile= args.outputFile

	sequences = [] 
	sequenceName = []
	seqLen = 0
	
	## reading data
	
	IN = open(inFile, 'r')
	for k,s in iterFasta(IN):		
		sequences.append(s.lower())
		sequenceName.append(k)
	IN.close()
	seqLen = len(sequences[0])	
	

	OUT = open(outFile,'w')

	nbComparable = countNumberComparablePositions( sequences )
	nbVariable = countNumberVariablePositions( sequences )
	nbSeqs = len(sequences)

	sumDiff = 0
	maxDiff = 0

	for i in range(nbSeqs-1):
		for j in range(i+1,nbSeqs):
			diff = getPairWiseDiff( sequences[i] , sequences[j] )

			print( sequenceName[i] , sequenceName[j] , diff , diff/nbComparable,  file=OUT , sep='\t')
			sumDiff+=diff
			if diff>maxDiff:
				maxDiff=diff

	Pi = 2*sumDiff / (nbSeqs**2)

	print( inFile , nbSeqs , seqLen , nbComparable , nbVariable , Pi , maxDiff , nbVariable/nbComparable , Pi/nbComparable , maxDiff/nbComparable )

# c('inFile' , 'nbSeqs' , 'seqLen' , 'nbComparable' , 'nbVariable' , 'Pi' , 'maxDiff' , 'frcVariable' , 'Pinorm' , 'maxDiffnorm' )