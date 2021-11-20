from mpi4py import MPI
from collections import defaultdict
import numpy as np
from scipy.sparse import csr_matrix
import time
import warnings
warnings.filterwarnings('ignore')

from pyKleeBarcode_utils import iterFasta, writeStructureMatrix_bin, getSequenceSpecies , getModalValueIncludingSpecialRules , sequenceToBarcode

from pyKleeBarcode_linearAlgebra import computeSsum_sparse , computeSingleIndicatorVector_sparse , computeStructureMatrixFromIndicatorVectors, computeAndWriteStructureMatrixFromIndicatorVectors



def sendBigOneByOne( data , comm , dest , tag ):
	'''
		data (dict): dict of  
		comm : mpi communicator
		dest : destination of the transmission to send to
		tag (int) : base tag. message send will go from tag+1 to tag+len(data)+1
	'''
	nbChunks = len(data)

	#sending the number of chunks we will send
	comm.send( nbChunks  , dest = dest , tag = tag+1 )
	chunkSize = int(np.ceil(len(data)/nbChunks))
	print('sending',nbChunks , chunkSize)

	fixedKeys = list(data.keys())

	for i in range(nbChunks):

		minIndex = i*chunkSize
		maxIndex =  min( (i+1)*chunkSize , len(data) ) 

		#subData = {}
		#for j in range( minIndex , maxIndex ):
		#	k = fixedKeys[j]
		#	subData[ k ] = data.pop(k)
		#	
		#	if total_size(subData)/(2**31 - 1.) > 1.0 :
		#		print("error! even in chunks, the message size is above 2**31 - 1. This is outside the scope of this function.")
		k = fixedKeys[i]
		#print('sending sub',i,k  )
		comm.send( k , dest = dest , tag = tag+2*i+2 )
		comm.send( data[k]  , dest = dest , tag = tag+2*i+3 )
		#print('sent sub',i )

	return 

def recvBigOneByOne( data , comm , source , tag ):
	'''
		data (dict): dict of  
		comm : mpi communicator
		source : source of the transmission to receive from
		tag (int) : base tag. message send will go from tag+1 to tag+len(data)+1
	'''
	#receiving the number of chunks we will send
	nbChunks = comm.recv( source = source , tag = tag+1 )

	print('receiving',nbChunks)

	for i in range(nbChunks):
		#print('receiving sub',i , tag+i+2)
		k = comm.recv( source = source , tag = tag+2*i+2 )
		subData = comm.recv( source = source , tag = tag+2*i+3 )
		#print('received sub',i)
		data[k]=subData
	return 



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

		parser.add_argument('--seed',type=int, default=None,
				 help='random seed used when selecting a species sequences if there is more than --max-seq-per-species. By default is it created using time.')


		args = parser.parse_args()

		inFile = args.inputFile
		outFile= args.outputFile

		#to be replaced by something better in time
		FASTA_ID_SEPARATOR= args.field_delimitor # separator between fields in the fasta id
		SPECIES_NAME_INDEX= args.species_field_index  # index of the species name in the fasta id (starts at 0)
	
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

		speciesSequences = defaultdict(list) # stores sequences, organized by species
		seqLen = 0
		
		## reading data
		IN = open(inFile, 'r')
		for k,s in iterFasta(IN):
			sp = getSequenceSpecies(k, FASTA_ID_SEPARATOR , SPECIES_NAME_INDEX )

			#if len ( speciesSequences[ sp ] ) == MAX_SEQUENCE_PER_SPECIES and MAX_SEQUENCE_PER_SPECIES >0 :
			#	continue ## this species already have the maximum number of sequences
			speciesSequences[ sp ].append( s )
			if seqLen == 0:
				seqLen = len(s)
			elif len(s) != seqLen:
				print( "Error : found sequence of differing length (",k,").","Please make sur the fasta file you provided contians aligned sequences only." , file = sys.stderr )
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


	## Computing indicator vectors for each species

	T = time.time()
	global_indic_vectors = None
	local_indic_vectors = {}

	for sp in speciesMatrices.keys():
		local_indic_vectors[sp] = computeSingleIndicatorVector_sparse( speciesMatrices[sp] , global_S_sum , nbSpecies )

	timerDict["local indicator vectors computation"] = time.time() - T
	# then, repatriation to process 0
	if MPI_rank != 0 :
		#comm.send( local_indic_vectors , dest = 0 , tag= MPI_rank * 10 +2 )
		sendBigOneByOne( local_indic_vectors , comm , dest=0 , tag=MPI_rank * 10 +2 )
	else:
		T2 = time.time()
		
		global_indic_vectors = local_indic_vectors
		for i in range(1,MPI_size):
			#x = comm.recv(source = i , tag = i*10+2)
			x={}
			recvBigOneByOne( x ,  comm , source=i , tag=i*10+2 )	
			#print(i,x.shape)
			global_indic_vectors.update( x )

	
		timerDict["indicator vectors repatriation"] = time.time() - T2

	if MPI_rank == 0 :		
		print("Indicator vectors computed.")

	## pushing all vectors in a single matrix and computing the structure matrix from them
	## at this point only the master process is working
	if MPI_rank == 0 :

		T = time.time()
		if nbSpecies <= 5000 : 

			structureMatrix = computeStructureMatrixFromIndicatorVectors( global_indic_vectors ,
	                                                                              speciesOrder , seqLen*4 )
	
			timerDict["total structure matrix computation"] = time.time() - T
	
			print("Structure matrix computed")
	
			## reporting
			T = time.time()
			writeStructureMatrix_bin( speciesOrder , structureMatrix , outFile )
			print("Structure matrix written in",outFile)
			timerDict["writing"] = time.time() - T
		else:
			
			OUT=open(outFile,'wb+')
			computeAndWriteStructureMatrixFromIndicatorVectors( global_indic_vectors , speciesOrder , seqLen*4 , OUT,header=speciesOrder )
			OUT.close()
			print("Structure matrix written in",outFile)
			timerDict["structure matrix computation and writing"] = time.time() - T



	timeLines = ''
	for k in timerDict.keys():
		timeLines +=" ".join([str(MPI_rank) , k , str(timerDict[k]) , 's']) + '\n'

	print(timeLines,end="")

	print("process",MPI_rank,"exited normally")

        
