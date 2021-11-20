
from pyKleeBarcode_utils import iterFasta , writeFasta


if __name__ == "__main__" :

	import sys
	import argparse

	parser = argparse.ArgumentParser(description='Filter out invariant positions from a multiple sequence alignment.')
	parser.add_argument('-i','--inputFile', type=str, required=True,
			 help='input multiple sequence alignment in fasta format')
	parser.add_argument('-o','--outputFile', type=str, default="",
			 help='output file. Prints to stdout if omitted')
	parser.add_argument('-A','--start',type=int , required = True,
			 help='Start reading from the base pair in this position, inclusive')
	parser.add_argument('-O','--stop',type=int , required = True,
			 help='Stop reading from the base pair in this position, exclusive')
	
	args = parser.parse_args()

	lenSeq = 0

	## reading data
	sequences = dict() # keys : fastaID , value : sequence
	seqOrder=[]
	IN = open(args.inputFile)
	for k,s in iterFasta(IN):
		sequences[k] = s
		seqOrder.append(k)
		lenSeq = len(s)
		#print(k,len(s))
	IN.close()

	## checking parameters
	if lenSeq < args.start:
		print( "ERROR: the given start position (",args.start,") is above the sequence length (",lenSeq,")." , file=sys.stderr )
		exit(1)
	if args.start >= args.stop:
		print( "ERROR: the given start position (",args.start,") is above or equal to the stop position (",args.stop,")." , file=sys.stderr )
		exit(1)

	lenSeqRemaining = 0
	for k in seqOrder:
		sequences[k] = sequences[k][args.start : args.stop]
		lenSeqRemaining = len(sequences[k])


	## output

	print( "number of remaining position after trimming:" , lenSeqRemaining )
	OUT = sys.stdout 
	if args.outputFile != "":
		OUT = open( args.outputFile , "w" )
	writeFasta( sequences , seqOrder , OUT )

