from collections import Counter
import sys
import numpy as np

def iterFasta(fileHandle):
	""" Takes :
			- fileHandle (file) : handle to input fasta file

		Yields :
			(str,str) : * sequence ID
						* sequence
	"""
	key=None
	seq=None
	l = fileHandle.readline()
	while l != "":
		
		if l.startswith(">"):
			if not key is None:
				yield key , seq
			key = l[1:].strip()
			seq = ""
		else:
			seq += l.strip().upper()
		l = fileHandle.readline()
	if not key is None:
		yield key,seq

def writeFasta( sequences , keyOrder , handleOUT ):
	"""
	Takes:
		- sequences (dict): keys are id ; values are sequences
		- keyOrder (list): list of sequence ids, giving the desired order of keys in the sequences dictionnary
		- handleOUT (file): handle to write to
	"""

	for k in keyOrder :
		print( ">" , k , file=handleOUT )
		print( sequences[k] , file=handleOUT )


def CSVmat(mat,outHandle,header=None):

	if not header is None:
		print(','.join(header),file=outHandle)
	for i in range(mat.shape[0]):
		print(','.join([str(x) for x in mat[i,:] ]),file=outHandle)
	


def readCSV(fname):
	""" reads a csv file as a numpy 2 dimensional array """	
	tabs = []

	IN=open(fname,'r')
	l = IN.readline()
	while l != "":
		tabs.append([ float(x) for x in l.strip().split(',') ])
		l = IN.readline()
	return np.array(tabs)

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


def getModalValueIncludingSpecialRules(V):
	"""
		Takes :
			- V (list) : list of nucleotides in a grouping at a given position

		Returns :
			(str) : the modal value in that species at that position 
					(provided there is 1 clear modal value which is neither - nor N, in which case N is returned)
	"""
	
	if len(V) ==1:
		if V[0] == '-':
			return 'N'
		return V[0]

	C = Counter()
	for c in V:
		C[c]+=1

	maxValue = []
	maxNumber=0

	for c,n in C.items():
		if n > maxNumber:
			maxNumber=n
			maxValue = [c]
		elif n == maxNumber:
			maxValue.append(c)


	if len(maxValue) == 1: # "normal case" : 1 mode 
		if maxValue[0] == '-': # if - -> do not change the Ns
			return 'N'
		return maxValue[0]
	elif 'N' in maxValue: 
		if len(maxValue) == 2 : # if there is 2 modes, but one of them is N
			other = maxValue[0] # finding the other non-N mode
			if other == 'N':
				other = maxValue[1]
			if other == '-': # if this mode is -, do not change N
				return 'N'
			return other #return the non-N mode

	# ambiguous mode 
	return 'N'

def sequenceToBarcode(seq):
	""" Takes :
			- seq (str) : DNA sequence ; authorized character are A,T,G,C,-,N

		Yields :
			(list) : list where each character of the sequence has been replaced by :
						A=[1,0,0,0];
						C=[0,1,0,0];
						G=[0,0,1,0];
						T=[0,0,0,1];
						-=[0,0,0,0];
						N=[0.25,0.25,0.25,0.25];

	"""
	encoding={'A':[1,0,0,0],
			  'C':[0,1,0,0],
			  'G':[0,0,1,0],
			  'T':[0,0,0,1],
			  'N':[0.25,0.25,0.25,0.25], 
			  '-':[0,0,0,0],
			  'Y':[0  ,0.5,0  ,0.5], #Pyrimidine (C or T)
			  'R':[0.5,0  ,0.5,0  ], #Purine (A or G)
			  'W':[0.5,0  ,0  ,0.5], #weak (A or T)
			  'S':[0  ,0.5,0.5,0  ], #strong (G or C)
			  'K':[0  ,0  ,0.5,0.5], #keto (T or G)
			  'M':[0.5,0.5,0  ,0  ], #amino (C or A)
			  'D':[1/3,0  ,1/3,1/3], #A, G, T (not C)
			  'V':[1/3,1/3,1/3,0  ], #A, C, G (not T)
			  'H':[1/3,1/3,0  ,1/3], #A, C, T (not G)
			  'B':[0  ,1/3,1/3,1/3]} #C, G, T (not A)





	vec = []
	for c in seq:
		if not c in encoding:
			print("!Warning! invalid character",c," will be replaced by N (only A,C,T,G,N,- and IUPAC authorized).", file=sys.stderr)
			c='N'
			#raise ValueError
		vec += encoding[c]
	return vec




def saveBinaryMatrix( arr, pth):
	"""function adapted from https://github.com/mverleg/array_storage_benchmark

	Takes:
		- arr : numpy array
		- pth : output path
	"""
	with open(pth, 'wb+') as fh:
		fh.write('{0:} {1:} {2:}\n'.format(arr.dtype, arr.shape[0], arr.shape[1]).encode('ascii'))
		fh.write(arr.data)
		

def loadBinaryMatrix( pth):
	"""function adapted from https://github.com/mverleg/array_storage_benchmark

	Takes:
		- pth : input path

	Returns:
		numpy array
	"""
	with open(pth, 'rb') as fh:
		header = fh.readline()
		data = fh.read()
	dtype, w, h = header.decode('ascii').strip().split()
	return np.frombuffer(data, dtype=dtype).reshape((int(w), int(h)))


def saveSsum( arr , N , pth):
	"""function adapted from https://github.com/mverleg/array_storage_benchmark

	Takes:
		- arr : numpy array . Ssum, a matrix representation of the diversity of a DNA sequence across a number of individuals or groups or individuals (typically, species)
							   according to the definitions of "A scalable method for analysis and display of DNA sequences" by Sirovitch et alii 
		- N : int . number of sequences (or group of sequences) used to create the array
		- pth : output path
	"""
	with open(pth, 'wb+') as fh:
		fh.write('{0:}\n'.format(N).encode('ascii'))
		fh.write('{0:} {1:} {2:}\n'.format(arr.dtype, arr.shape[0], arr.shape[1]).encode('ascii'))
		fh.write(arr.data)
		

def loadSsum( pth):
	"""function adapted from https://github.com/mverleg/array_storage_benchmark

	Takes:
		- pth : input path

	Returns:
		numpy array : Ssum, a matrix representation of the diversity of a DNA sequence across a number of individuals or groups or individuals (typically, species)
							   according to the definitions of "A scalable method for analysis and display of DNA sequences" by Sirovitch et alii 
		int : number of sequences (or group of sequences) used to create the array
	"""
	with open(pth, 'rb') as fh:
		header1 = fh.readline()
		header2 = fh.readline()
		data = fh.read()
	N = int( header1.decode('ascii').strip() )
	dtype, w, h = header2.decode('ascii').strip().split()
	return np.frombuffer(data, dtype=dtype).reshape((int(w), int(h))) , N


def writeIndicatorVectors( IVdict, speciesOrder , fileOut ):
	### simple csv writer for now
	with open(fileOut,'w') as OUT:
		for s in speciesOrder:
			#print( IVdict[s].shape )
			print( s , ','.join([str(float(x)) for x in IVdict[s] ]) , sep=',' , file=OUT )

def readIndicatorVectors( fileIn ):
	### simple csv reader for now
	IVdict = {}
	speciesOrder = []
	with open(fileIn,'r') as IN:

		for l in IN:
			sl = l.strip().split(',')
			s = sl[0]
			speciesOrder.append(s)
			IVdict[s] = np.array( [float(x) for x in sl[1:]] )
	return IVdict,speciesOrder

def yieldIndicatorVectors( fileIn ):
	### simple csv reader for now
	with open(fileIn,'r') as IN:

		for l in IN:
			sl = l.strip().split(',')
			s = sl[0]

			arr = np.array( [float(x) for x in sl[1:]] )
			yield s , arr
	



def writeStructureMatrix_line_bin( seqName , corVec , outHandle ):
	"""
	Takes:
		- seqName : str . name of the sequence
		- corVec : np.array . correlation between sequences and the current one
		- outHandle : writing handle 

	"""
	
	outHandle.write('{0:}\n'.format(seqName).encode('ascii'))
	outHandle.write( corVec.data )
		

def writeStructureMatrix_bin( seqNames , mat , path ):
	"""
	Takes:
		- seqNames : list . list of the name for each column/rows in the matrix
		- mat : numpy array . structure matrix : NxN matrix containing the correlation between sequences
		- path : output path
	"""
	with open(path, 'wb+') as fh:
		for i,n in enumerate( seqNames ):
			writeStructureMatrix_line_bin( n , mat[i,:i] , fh )


def readStructureMatrix_CSV( fileName ):
    """ reads a csv file as a numpy 2 dimensional array, with a header """ 
    tabs = []

    IN=open(fileName,'r')
    l = IN.readline()
    header = l.strip().split(',')
    l = IN.readline()
    while l != "":
        tabs.append([ float(x) for x in l.strip().split(',') ])
        l = IN.readline()
    return header , np.array(tabs)


def yieldStructureMatrix_bin( inHandle ):
	"""

	Takes:
		- inHandle : input binary handle

	Yield:
		str : sequence name
		numpy array : correlation of the sequence with some others 
					(corresponding to a line in the lower triangular matrix)
	"""

	i = 0
	name = inHandle.readline().decode('ascii').strip()
	while name:
		arr = np.fromfile( inHandle , dtype=np.float64 , count = i )
		yield name, arr
		i+=1
		name = inHandle.readline().decode('ascii').strip()

def readStructureMatrix_bin( pth ):
	"""

	Takes:
		- pth : str. name of the file to read

	Yield:
		list : list contnaing sequence names
		numpy array :  2D structure matrix

	"""

	seqList = []
	lines = []

	with open( pth , 'rb' ) as IN:
		for h,x in yieldStructureMatrix_bin( IN ):
			
			seqList.append(h)
			lines.append(x)

	M = np.ones( ( len(lines) , len(lines) ) )

	for i,l in enumerate( lines ):
		M[i,:i] = l
		M[:i,i] = l

	return seqList,M
