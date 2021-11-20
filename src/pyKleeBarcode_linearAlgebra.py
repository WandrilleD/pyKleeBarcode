import numpy as np
from numpy import linalg as LA
from scipy.sparse.linalg import eigs as eigSparse, ArpackNoConvergence
from scipy.sparse import csr_matrix
from pyKleeBarcode_utils import writeStructureMatrix_line_bin

#from memory_profiler import profile

def computeSsum( speciesMatrices , spList , seqLen ):
	S_Sum = np.zeros( ( seqLen , seqLen ) )
	for sp in spList :
		S_Sum += selfTransposeComput( speciesMatrices[sp] )
	return S_Sum


def selfTransposeComput( M ):
	return M.T @ M

def computeM( S_sum , S , nbSpecies ):
	return S - 1/(nbSpecies-1) * (S_sum-S) 


def selfTransposeComput_sparse( M ):
	"""
		Takes:
			- M (scipy.sparse.csr_matrix): a (rectangular) matrix

		Return:
			(scipy.sparse.csr_matrix) : product of the transpose of M with M
	"""
	return M.transpose().dot(M)

def computeSsum_sparse( speciesMatrices , spList , seqLen ):
	"""
		Takes:
			- speciesMatrices (dict) : keys are species/group names, values are the sparse matrix describing their associated sequences
			- spList (list) : ordered list of species/group names
			- seqLen (int) : length of the sequence in bool (4 times the nucleotide length)

		Return:
			(scipy.sparse.csr_matrix) : sum of the product of the transpose of the species matrix with themselves
	"""

	S_Sum = csr_matrix( ( seqLen , seqLen ) , dtype=float )
	for sp in spList :
		S_Sum += selfTransposeComput_sparse( speciesMatrices[sp] ) / speciesMatrices[sp].shape[0] 
	return S_Sum


def computeM_sparse( S_sum_sparse , S , nbSpecies ):
	return S - 1/(nbSpecies-1) * (S_sum_sparse-S) 


def computeSingleIndicatorVector_sparse( speciesM , S_sum_sparse , nbSpecies ):
	"""
		Takes:
			- speciesM (scipy.sparse.csr_matrix) : a (rectangular) matrix representing the species/group sequences
			- S_sum_sparse  (scipy.sparse.csr_matrix) : sum of the product of the transpose of the species matrix with themselves
			- nbSpecies (int) : number of species/groups in the computations

		Returns:
			(numpy.ndarray) : eigenVector corresponding to the species

	"""

	S = selfTransposeComput_sparse( speciesM ) / speciesM.shape[0] 
	
	M = computeM_sparse( S_sum_sparse , S , nbSpecies )
		
	# eigen value/vector computation
	# if the matrix is sparse, then using the dedicated scipy function is ~50x faster

	w , v = None , None
		
	try:
		w,v = eigSparse( M , k=1 , which='LR' , maxiter=20)
	except ArpackNoConvergence as e:
		#print("species {} , index {}".format(sp,i))
		#print(e)
		#w = e.eigenvalues
		#v = e.eigenvectors
		#print(w.shape, v.shape)
		
		print("species {} , index {} : ARPACK convergence problem".format(sp,i+1))
		print("no eigenvector : relying on non-sparse computations." )
		w, v = LA.eig( M.toarray() )
	except: 
		# I sometimes got a failure while testing when the matrix was not sparse enough, 
		# so I put this failsafe here
		# so far, I have not had any problems with data derived from actual sequence alignments
		w, v = LA.eig( M.toarray() )

	#We keep the maximum eigenvalue and associated vector
	maxIndex = np.argmax(w)
	return  v[:,maxIndex] 




def computeStructureMatrixFromIndicatorVectors( indicVectorDict , speciesOrder , seqLen ):
	"""
		Takes:
			- indicVectorDict (dict) : keys are species/group names, values are the corresponding indicator vector (numpy.array)
			- speciesOrder (list) : desired order of the indicVectorDict keys
			- seqLen (int) : length of the studied sequence (4* nucleotide size)

		Returns:
			(numpy.ndarray) : structure matrix
	"""
	nbSpecies = len(speciesOrder)
	V = np.zeros( ( seqLen , nbSpecies ) )
	for i,sp in enumerate(speciesOrder):
		V[:,i] =  indicVectorDict[sp]

	V=-V
	# InnerV is the structure matrix
	InnerV=V.T @ V

	return abs(InnerV)

def computeAndWriteStructureMatrixFromIndicatorVectors( indicVectorDict , speciesOrder , seqLen , outHandle,header=None ):
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

	return


#@profile
def computeIndicatorVectors( speciesMatrices , speciesOrder , verbose = True):
	"""
		Takes:
			- speciesMatrices (dict) : keys are species (or groups) names
									   values are a np.matrix containing the sequences associated to the species in [0,0,0,1]-like vector form
			- speciesOrder (list) : list of species name
			- verbose = True (bool)

		Returns:
			(np.array) : structure matrix, giving the correlations between the indicative vectors of species/groups
	"""


	nbSpecies = len(speciesMatrices)
	spList = speciesOrder
	seqLen = speciesMatrices[ spList[0] ].shape[1]


	#S_sum <- for each group/species : multiply with transposed and add 
	speciesMatrices_sparse = { sp : csr_matrix(m) for sp , m in speciesMatrices.items()}
	
	S_sum_sparse = computeSsum_sparse( speciesMatrices_sparse , speciesOrder , seqLen )

	V = np.zeros( ( seqLen , nbSpecies ) )
	D = np.zeros( ( nbSpecies , 1 ) )


	# core computationnal load : compute eigen vectors for each species
	#	the matrix we compute the eigen vector for is the matrix specific of that group minus the mean of the other species matrices
	# How would we handle species with a different number of individuals? I presume we divide each species matrix by the number of individuals...  

	for i,sp in enumerate( spList ) :
		
		

		S = selfTransposeComput_sparse( speciesMatrices_sparse[sp] ) / speciesMatrices_sparse[sp].shape[0] 
		# so this line repeated ... 
		# but we have to compare to cost of computing this twice with 
		# the cost of stocking a 2000 by 2000 matrix for each species, which is a lot

		M = computeM_sparse( S_sum_sparse , S , nbSpecies )
		
		# eigen value/vector computation
		# if the matrix is sparse, then using the dedicated scipy function is ~50x faster

		w , v = None , None
		
		try:
			w,v = eigSparse( M , k=1 , which='LR' , maxiter=20)
		except ArpackNoConvergence as e:
			#print("species {} , index {}".format(sp,i))
			#print(e)
			#w = e.eigenvalues
			#v = e.eigenvectors
			#print(w.shape, v.shape)
			
			print("species {} , index {} : ARPACK convergence problem".format(sp,i+1))
			print("no eigenvector : relying on non-sparse computations." )
			w, v = LA.eig( M.toarray() )

		except: 
			# I sometimes got a failure while testing when the matrix was not sparse enough, 
			# so I put this failsafe here
			# so far, I have not had any problems with data derived from actual sequence alignments
			if verbose:
				print( "species",sp,": needs to rely on non-sparse computations." )
			w, v = LA.eig( M.toarray() )

		#We keep the maximum eigenvalue and associated vector
		maxIndex = np.argmax(w)

		V[:,i] =  v[:,maxIndex] 
		D[i] = w[maxIndex]
		# now we have to extract the eigen vector associated to the highest value.

		if verbose:
			print("treated species : {:5}/{}".format(i+1,len(spList)) , end='\r')


	V=-V
	# InnerV is the structure matrix
	InnerV=V.T @ V

	return abs(InnerV)
