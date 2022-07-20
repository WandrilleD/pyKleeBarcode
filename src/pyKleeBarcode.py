
from collections import defaultdict
import numpy as np
from scipy.sparse import csr_matrix
import time
import warnings
warnings.filterwarnings('ignore')

from pyKleeBarcode_utils import iterFasta, writeStructureMatrix_bin, getSequenceSpecies , getModalValueIncludingSpecialRules , sequenceToBarcode, readSequenceSpeciesCorrespondenceFile

from pyKleeBarcode_linearAlgebra import computeSsum_sparse , computeSingleIndicatorVector_sparse , computeStructureMatrixFromIndicatorVectors, computeAndWriteStructureMatrixFromIndicatorVectors, computeIndicatorVectors




if __name__ == "__main__":

	import sys
	import argparse

	parser = argparse.ArgumentParser(
				description="""Computes a structure matrix between sets of DNA sequences (typically grouped by species)
							   according to the definitions of "A scalable method for analysis and display of DNA sequences" by Sirovitch et alii 
								""")
	parser.add_argument('-i','--inputFile', type=str, required=True,
			 help='input multiple sequence alignment in fasta format (preferably, trimmed)')
	parser.add_argument('-o','--outputFile', type=str, required=True,
			 help='output file name for the structural matrix')

	parser.add_argument('-m','--max-seq-per-species', type=int, default=5,
			 help='maximum number of sequences kept per species. Default : 3. Set to 0 to have no limit ; be careful thought, as this parameter has a direct impact on speed.')

	parser.add_argument('-f','--field-delimitor', type=str, default="|",
			 help='field delimitor of the fasta sequence id lines. default: "|"')

	parser.add_argument('-s','--species-field-index', type=int, default=1,
			 help='index (starting at 0) of the species name in the fasta sequence id lines. default: 1')

	parser.add_argument('-C','---correspondence-file', type=str, default="",
		 help='name of a comma-delimited file containing correspondence between sequence and species. Overrides option -s when used.\nOne sequence ID per line, sequence id and species name separated by a semicolumn (;).')


	parser.add_argument('--seed',type=int, default=None,
			 help='random seed used when selecting a species sequences if there is more than --max-seq-per-species. By default is it created using time.')


	args = parser.parse_args()

	inFile = args.inputFile
	outFile= args.outputFile

	# maximum of sequence per species/group
	MAX_SEQUENCE_PER_SPECIES = args.max_seq_per_species
	randomSeed = args.seed
	if MAX_SEQUENCE_PER_SPECIES > 0 and randomSeed is None:
		randomSeed = int(time.time())
	if MAX_SEQUENCE_PER_SPECIES > 0:
		print("random seed for species sequence selection :" , randomSeed)

	np.random.seed(randomSeed)

	#to be replaced by something better in time
	FASTA_ID_SEPARATOR= args.field_delimitor # separator between fields in the fasta id
	SPECIES_NAME_INDEX= args.species_field_index  # index of the species name in the fasta id (starts at 0)
	correspondence_file = args.correspondence_file	
	
	
	
	speciesIDs = defaultdict(list) # store ids of sequences, organized by species
	sequences = dict() # keys : fastaID , value : sequence
	speciesOrder = []
	seqLen = 0
	
	## reading data

	sid2sp = {}
	if len(correspondence_file) > 0 : ## we use the correspondence file instead of the Nth field of the fasta description line
		sid2sp = readSequenceSpeciesCorrespondenceFile( correspondence_file , sep = ';' )

	
	IN = open(inFile, 'r')
	for k,s in iterFasta(IN):
		sp=''
		if len(correspondence_file) > 0 : ## we use the correspondence file instead of the Nth field of the fasta description line
			sid = k.partition(FASTA_ID_SEPARATOR)[0]

			if not sid in sid2sp :
				print("Error : sequence {} is absent from the correspondence file {}".format( sid , correspondence_file ))
				exit(1)
			sp = sid2sp[sid]
		else:
			sp = getSequenceSpecies(k, FASTA_ID_SEPARATOR , SPECIES_NAME_INDEX )

		speciesIDs[ sp ].append( k )
		sequences[k] = s
		seqLen = len(s)
	
	IN.close()
	
	speciesOrder = [sp for sp in speciesIDs.keys()] # keeping the order of species somewhere

	## replacing unknown nucleotides (N) by either the modal value in that species at that position (provided there is 1 clear modal value which is neither - nor N)
	for sp,ids in speciesIDs.items():

		i = 0
		while i < seqLen :
			characterVector = []
			hasN=[]
			for ID in ids:
				if sequences[ID][i] == 'N' :
					hasN.append( ID ) # refer4encing sequences with a N
				characterVector.append(sequences[ID][i])
			if len(hasN) > 0 :
				replaceNcharacter = getModalValueIncludingSpecialRules(characterVector) # used to find the most frequent character which will replace the N characters in sequences
				for ID in hasN: # replacing the N character
					sequences[ID] = sequences[ID][:i] + replaceNcharacter + sequences[ID][i+1:]
					
			i+=1 

	## transforming nucleotide sequences into numeric vector and grouping them into species-level matrices
	speciesMatrices = {}

	for sp,ids in speciesIDs.items():
		
		#limiting number of sequences per species
		nbSeq = len(ids)
		if MAX_SEQUENCE_PER_SPECIES > 0 :

			np.random.seed( randomSeed )
			if nbSeq > MAX_SEQUENCE_PER_SPECIES :
				ids = np.random.choice(ids , MAX_SEQUENCE_PER_SPECIES , replace=False)

				nbSeq = MAX_SEQUENCE_PER_SPECIES


		#if nbSeq != MAX_SEQUENCE_PER_SPECIES:
		#	continue

		m = np.zeros(( nbSeq , 4 * (seqLen) ))
		
		for i,ID in enumerate(ids):
			m[i,] = sequenceToBarcode(sequences[ID])

		speciesMatrices[sp] = m # csr_matrix(m) < we could stock these as sparse already, but at the moment this is not too big a mem load so I prefer to keep them as ndarray for the sake of generality
	print("detected",len(speciesMatrices),"species")
	print("sequence of ",seqLen,"residues")



	if len(speciesOrder)<=50000:

		## computing indicator vectors
		structureMatrix = computeIndicatorVectors(speciesMatrices,speciesOrder)

		## reporting
		OUT=open(outFile,'w')
		writeStructureMatrix_bin( speciesOrder , structureMatrix , outFile )
		OUT.close()


	else: # alternative procedure where the structure matrix is not explicitely declared in memory, when it becomes too large
		
		S_sum =  computeSsum_sparse( speciesMatrices , speciesOrder , seqLen*4 )
		indicVec = {}

		for sp in speciesMatrices.keys():
			indicVec[sp] = computeSingleIndicatorVector_sparse( speciesMatrices[sp] , S_sum , len(speciesMatrices) )


		OUT=open(outFile,'w')
		computeAndWriteStructureMatrixFromIndicatorVectors( indicVec , speciesOrder , seqLen*4 , OUT,header=speciesOrder )
		OUT.close()
		
	print("Structure matrix written in",outFile)
		

	






