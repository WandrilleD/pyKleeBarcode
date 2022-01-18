import numpy as np
from pyKleeBarcode_utils import readStructureMatrix_bin, writeStructureMatrix_bin


def readStructureMatrix( fileName ):
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



def applyOrderToStructureMatrix( structureMatrix , currentOrder , newOrder ):
    """
        Very naive implementation, re-built the whole matrix

        Takes:
            structureMatrix (np.array) : 2D structure matrix
            currentOrder (list) : species name in the order that they have in the provided structure matrix
            newOrder (dict) : keys : species name , values : desired index

        Returns:
            (np.array) : 2D structure matrix, reordered
    """

    ## adding at the end sequences that are absent from the order 
    for i,sp in enumerate(currentOrder):
        if not sp in newOrder:
            print("WARNING :",sp,"is absent from the list : ignored in the new matrix")
            

    newMatrix = np.zeros( (len(newOrder),len(newOrder)) )

    for i,sp1 in enumerate(currentOrder):


        new1 = newOrder.get( sp1 , None)
        if new1 is None:
            continue

        for j,sp2 in enumerate(currentOrder):
            
            new2 = newOrder.get( sp2 )
            if new2 is None:
                continue

            newMatrix[new1][new2] = structureMatrix[i][j]

    return newMatrix


def CSVmat(mat,outHandle,header=None):

    if not header is None:
        print(','.join(header),file=outHandle)
    for i in range(mat.shape[0]):
        print(','.join([str(x) for x in mat[i,:] ]),file=outHandle)
    




if __name__ == "__main__" :

    import sys
    import argparse

    parser = argparse.ArgumentParser(
                description="""Reorders a structure matrix from a tree.
                                The tree leaves must correspond to names in the matrix header.
                                The new order corresponds to the Depth First traversal of the tree.
                            """)
    parser.add_argument('-m','--inputMatrix', type=str, required=True,
             help='input structure matrix obtained from an alignment')
    parser.add_argument('-i','--inputOrder', type=str, required=True,
             help='a file the desired order (one id per line)')

    parser.add_argument('-o','--outputFile', type=str, required=True,
             help='output file name for the reordered matrix.')
    parser.add_argument('--csv', action='store_true',
             help='specify that the structure matrix in in csv format (default expects a binary format)')


    args = parser.parse_args()

    STRMAT_FILE = args.inputMatrix
    ORDER_FILE = args.inputOrder

    outFile= args.outputFile
    
    # reading data
    orderDict = {}
    orderL = []
    with open(ORDER_FILE,'r') as IN:
        for l in IN:
            orderDict[ l.strip() ] = len(orderDict)
            orderL.append( l.strip() )
    print("order read")

    
    smreading = readStructureMatrix_bin
    if args.csv :
        smreading = readStructureMatrix

    matrixheader , M = smreading( STRMAT_FILE )

    print("matrix read")

    
    print( "The list and the matrix have {intSize} leaves in common".format( intSize= len(orderDict) ) )


    
    M = applyOrderToStructureMatrix( M , matrixheader , orderDict )

    print("matrix reordered")

    ## reporting
    if not args.csv:
        writeStructureMatrix_bin( orderL , M , outFile  )
    else:
        OUT=open(outFile,'w')
        CSVmat(M,OUT,header=orderL)
        OUT.close()

    print("ordered matrix written in",outFile)    
