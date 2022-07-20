from mpi4py import MPI
from collections import defaultdict
import numpy as np
from scipy.sparse import csr_matrix
import time

from pyKleeBarcode_utils import iterFasta, CSVmat, getSequenceSpecies , getModalValueIncludingSpecialRules , sequenceToBarcode , saveSsum , readSequenceSpeciesCorrespondenceFile

from pyKleeBarcode_linearAlgebra import computeSsum_sparse


if __name__ == "__main__":
	import sys
        
	## init the world
	comm = MPI.COMM_WORLD
	MPI_size = comm.Get_size()
	MPI_rank = comm.Get_rank()

	## init some variables
	inFile = None
	outFile= None
	MAX_SEQUENCE_PER_SPECIES = None
	randomSeed = None

	#to be replaced by something better in time
	FASTA_ID_SEPARATOR= '|' # separator between fields in the fasta id
	SPECIES_NAME_INDEX= 1  # index of the species name in the fasta id (starts at 0)

	speciesPerProcess = [ [] for i in range(MPI_size)] # to keep track of which processes has what

	speciesSequences = {}
	# keys : sequence id
	# values: list of sequences in str format

	seqLen = 0 # length of the aligned sequence
	nbSpecies = 0 # total number of species
        
	print("process",MPI_rank,'/',MPI_size,"initalized.")
	## rank 0 is the master process


	timerDict = {}

	## master reads argument
	if MPI_rank == 0 :

	
		import argparse

		parser = argparse.ArgumentParser(
				description="""Computes Ssum, a matrix representation of the diversity of a DNA sequence across a number of individuals or groups or individuals (typically, species)
							   according to the definitions of "A scalable method for analysis and display of DNA sequences" by Sirovitch et alii 
								""")
		parser.add_argument('-i','--inputFile', type=str, required=True,
			 help='input multiple sequence alignment in fasta format (preferably, trimmed)')
		parser.add_argument('-o','--outputFile', type=str, required=True,
			 help='output file name for the Ssum matrix')

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

		#to be replaced by something better in time
		FASTA_ID_SEPARATOR= args.field_delimitor # separator between fields in the fasta id
		SPECIES_NAME_INDEX= args.species_field_index  # index of the species name in the fasta id (starts at 0)
	
		correspondence_file = args.correspondence_file

		# maximum of sequence per species/group
		MAX_SEQUENCE_PER_SPECIES = args.max_seq_per_species

		randomSeed = args.seed
		if MAX_SEQUENCE_PER_SPECIES > 0 and randomSeed is None:
			randomSeed = int(time.time())
		if MAX_SEQUENCE_PER_SPECIES > 0:
			print("random seed for species sequence selection :" , randomSeed)


	
	## master process reads the data
	if MPI_rank == 0:	


		T = time.time()

		sid2sp = {}
		if len(correspondence_file) > 0 : ## we use the correspondence file instead of the Nth field of the fasta description line
			sid2sp = readSequenceSpeciesCorrespondenceFile( correspondence_file , sep = ';' )


		speciesSequences = defaultdict(list) # stores sequences, organized by species
		seqLen = 0
		
		## reading data
		IN = open(inFile, 'r')
		for k,s in iterFasta(IN):
			sp=''
			if len(correspondence_file) > 0:
				sid = k.partition(FASTA_ID_SEPARATOR)[0]

				if not sid in sid2sp :
					print("Error : sequence {} is absent from the correspondence file {}".format( sid , correspondence_file ))
					exit(1)
				sp = sid2sp[sid]

			else:
				sp = getSequenceSpecies(k, FASTA_ID_SEPARATOR , SPECIES_NAME_INDEX )

			#if len ( speciesSequences[ sp ] ) == MAX_SEQUENCE_PER_SPECIES and MAX_SEQUENCE_PER_SPECIES >0 :
			#	continue ## this species already have the maximum number of sequences
			speciesSequences[ sp ].append( s )

			if seqLen == 0:
				seqLen = len(s)
			elif len(s) != seqLen:
				print( "Error : found sequence of differing length (",k,").","Please make sur the fasta file you provided contains aligned sequences only." , file = sys.stderr )
				exit(1)

		
		IN.close()
	
		speciesOrder = [sp for sp in speciesSequences.keys()] # keeping the order of species somewhere
		nbSpecies = len(speciesOrder)

		print("finished reading data.",len(speciesOrder),"species detected.")

		timerDict["reading"] = time.time() - T

	## broadcasting the sequence length
	seqLen = comm.bcast(seqLen, root=0)
	nbSpecies = comm.bcast(nbSpecies, root=0)

	## master process splits the species sequences between the different processes
	T = time.time()

	if MPI_rank == 0 :


		for i, sp in enumerate( speciesOrder ):
			speciesPerProcess[ i%MPI_size ].append( sp )

		## sending info to the other processes
		request_list = []

		for i in range(1,MPI_size):
			## first, build the dict we want to send
			localDict = { sp : speciesSequences.pop(sp) for sp in speciesPerProcess[i]}
			#print({ k:len(s) for k,s in localDict.items()})
			#localDict = { '>'*(x+1) : ['A'*(i+1) for j in range(3*x) ] for x in range(i**2)}
			#print(localDict)
			## then send the actual data
			request_list.append( comm.isend( localDict , dest= i , tag= i * 10 ) )
			#print('sending',len(localDict),'species to ',i)

		for r in request_list:
			r.wait()


		## the local data to process should be what is left in speciesSequences
		##print( set( speciesSequences.keys() ) == set( speciesPerProcess[0] )  )

	else:
		## receiving the species sequences
		print('process',MPI_rank,'waiting to receive')

		#req = comm.irecv(source=0, tag= MPI_rank * 10 )
		#print('...')
		#speciesSequences = req.wait()
		speciesSequences = comm.recv(source=0 , tag = MPI_rank*10)
		#print('...success',MPI_rank)
                
	local_speciesOrder = [sp for sp in speciesSequences.keys()]

	print("process",MPI_rank ,":", len(speciesSequences) ,"species")
	
	timerDict["species dispersion"] = time.time() - T

	## replacing unknown nucleotides (N) by either the modal value in that species at that position (provided there is 1 clear modal value which is neither - nor N)
	T = time.time()
	for sp,sequences in speciesSequences.items():

		## we scan the species alignment position by position
		i = 0
		while i < seqLen :
			characterVector = []
			hasN=[]
			for j,seq in enumerate(sequences):
				if seq[i] == 'N' :
					hasN.append( j ) # referencing sequences with a N
				characterVector.append(seq[i])
			if len(hasN) > 0 :
				replaceNcharacter = getModalValueIncludingSpecialRules(characterVector) # used to find the most frequent character which will replace the N characters in sequences
				for j in hasN: # replacing the N character
					sequences[j] = sequences[j][:i] + replaceNcharacter + sequences[j][i+1:]
			i+=1 
	timerDict["replacing Ns"] = time.time() - T

	## limiting the number of sequences per species
	MAX_SEQUENCE_PER_SPECIES = comm.bcast(MAX_SEQUENCE_PER_SPECIES, root=0)
	randomSeed = comm.bcast(randomSeed, root=0)
	if MAX_SEQUENCE_PER_SPECIES > 0 : #otherwise it is presumed that all sequences should be kept
		np.random.seed( randomSeed )
		for sp in local_speciesOrder:
			if len(speciesSequences[sp]) > MAX_SEQUENCE_PER_SPECIES :
				speciesSequences[sp] = np.random.choice(speciesSequences[sp] , MAX_SEQUENCE_PER_SPECIES , replace=False)



	## transforming nucleotide sequences into numeric vector and grouping them into species-level matrices
	T = time.time()

	speciesMatrices = {}

	for sp,sequences in speciesSequences.items():
		
		nbSeq = len(sequences)
		
		m = np.zeros(( nbSeq , 4 * (seqLen) ))
		
		for i,seq in enumerate(sequences):

			m[i,] = sequenceToBarcode(seq)

		speciesMatrices[sp] = csr_matrix(m) # we stock these as a sparse matrix

	timerDict["creating sparse matrices"] = time.time() - T

	print(MPI_rank,"created matrices representation of alignment")

	## computation of the matrix S_sum , which is the
	## sum of the product of the transpose of the species matrix with themselves	
	T = time.time()

	global_S_sum = None

	# first, local computation
	local_S_sum =  computeSsum_sparse( speciesMatrices , local_speciesOrder , seqLen*4 )

	

	timerDict["local S_sum computation"] = time.time() - T
	T2 = time.time()
	# then, repatriation to process 0 for the overall sum

	if MPI_rank != 0 :
		comm.send( local_S_sum , dest = 0 , tag= MPI_rank * 10 +1 )
	else:
		global_S_sum = local_S_sum
		for i in range(1,MPI_size):
			x = comm.recv(source = i , tag = i*10+1)
			#print(i,x.shape)
			global_S_sum += x

	## now, broadcasting the global S sum matrix
	global_S_sum = comm.bcast(global_S_sum, root=0)

	timerDict["local to global S_sum"] = time.time() - T2
	timerDict["total S_sum computation"] = time.time() - T

	if MPI_rank == 0 :
		print("S_sum matrix computed")


	print(timerDict)

	if MPI_rank == 0 :
		print( 'non-zeros', global_S_sum.count_nonzero() )
		print( 'shape' , global_S_sum.shape )
		print( 'shape' , global_S_sum.dtype )
		global_S_sum = global_S_sum.toarray()
		saveSsum( global_S_sum , nbSpecies , args.outputFile)


