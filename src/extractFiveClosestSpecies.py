import sys

def updateFiveBest( fiveBest , fiveBestVal , new , val ):
	
	if fiveBestVal[-1]>= val:
		return 

	if fiveBestVal[2]>= val:
		if fiveBestVal[3]>= val:
			fiveBest[-1] = new
			fiveBestVal[-1] = val
		else: # to index 3
			fiveBest.insert( 3 , new )
			fiveBestVal.insert( 3 , val )
			fiveBest.pop(-1)
			fiveBestVal.pop(-1)
	else: 
		if fiveBestVal[1]>= val :
				fiveBest.insert( 2 , new )
				fiveBestVal.insert( 2 , val )
				fiveBest.pop(-1)
				fiveBestVal.pop(-1)
		else :
			if fiveBestVal[0]>= val:
				fiveBest.insert( 1 , new )
				fiveBestVal.insert( 1 , val )
				fiveBest.pop(-1)
				fiveBestVal.pop(-1)
			else:
				fiveBest.insert( 0 , new )
				fiveBestVal.insert( 0 , val )
				fiveBest.pop(-1)
				fiveBestVal.pop(-1)


spList = []
fiveBest = []
fiveBestVal = []


with open(sys.argv[1],'r') as IN:
	l = IN.readline()
	spList = l.strip().split(',')

	nbSpecies = len(spList)
	for i in range(nbSpecies):
		fiveBest.append( [None,None,None,None,None] )
		fiveBestVal.append( [0,0,0,0,0] )

	lineIndex = 0

	for l in IN:

		sl = l.strip().split(',')

		for spIndex, val in enumerate(sl):

			if spIndex == lineIndex:
				continue
			val = float(val)

			updateFiveBest( fiveBest[spIndex] , fiveBestVal[spIndex] , lineIndex , val )

		lineIndex += 1

print(fiveBest)


OUT = open(sys.argv[2],'w')
for spIndex , sp in enumerate( spList ):
	print( sp , "\t" , ','.join([ spList[i] for i in fiveBest[spIndex] ]) , "\t" , ','.join([ str(val) for i in fiveBestVal[spIndex] ]) ,file=OUT)
