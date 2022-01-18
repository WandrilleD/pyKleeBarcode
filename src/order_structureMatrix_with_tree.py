import numpy as np
from pyKleeBarcode_utils import readStructureMatrix_bin, writeStructureMatrix_bin

def readTree( treeLine ):
    c,t = readTreeAux( treeLine , current=0)
    return t

def readTreeAux( treeLine , current=0):
    """
    ## # testing 
    ## s = "((abc:0.3,def:1):12,hij,(klm,(nop,qrs:0.6))test:plop);"
    ## c,t = readTree( s , current=0)
    ## print( t==[['abc', 'def'], 'hij', ['klm', ['nop', 'qrs']]] )
    """
    #print(treeLine)
    #print(" "*current + "^")

    children = []

    if treeLine[current] == '(':
        current, ch = readTreeAux( treeLine , current= current + 1)
        children.append(ch)



    while treeLine[current] == ',':
        current , ch = readTreeAux( treeLine , current= current + 1)
        children.append(ch)

    #print("after ch", treeLine)
    #print("after ch", " "*current + "^")

    ## at this point, current is ')' or the beginning of a leaf 
    if len(children)>0: # expect a ')' which signals the end of children 
        current +=1

    i=current
    while not treeLine[i] in [',',')',';']:
        i+=1

    name,junk,blen = treeLine[current:i].partition(':')
    
    current = i
    

    if len(children) == 0 :
        #print(name)
        return current , name

    #print(children)
    return current , children


def treeTraversal( tree ):
    '''DFS traversal
    ## # testing
    ## s = "((abc:0.3,def:1):12,hij,(klm,(nop,qrs:0.6))test:plop);"
    ## c,t = readTree( s , current=0)
    ## print(  treeTraversal( t ) == ['abc', 'def', 'hij', 'klm', 'nop', 'qrs'] )

    '''
    if isinstance(tree,str):
        ## leaf
        return [tree]

    orderedLeaves = []
    for c in tree:
        orderedLeaves += treeTraversal( c )

    return orderedLeaves


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

    ## ignoring sequences that are absent from the order 
    for i,sp in enumerate(currentOrder):
        if not sp in newOrder:
            print("WARNING :",sp,"is absent from the tree : ignored in the new matrix")
            

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
    parser.add_argument('-t','--inputTree', type=str, required=True,
             help='a file containing a tree in newick format')

    parser.add_argument('-o','--outputFile', type=str, required=True,
             help='output file name for the reordered matrix.')
    parser.add_argument('--csv', action='store_true',
             help='specify that the structure matrix in in csv format (default expects a binary format)')


    args = parser.parse_args()

    STRMAT_FILE = args.inputMatrix
    TREE_FILE = args.inputTree

    outFile= args.outputFile
    
    # reading data
    IN = open(TREE_FILE,'r')
    treeLine = IN.readline()
    IN.close()
    tree = readTree( treeLine )
    #print(tree)

    print("tree read")
    smreading = readStructureMatrix_bin
    if args.csv :
        smreading = readStructureMatrix

    matrixheader , M = smreading( STRMAT_FILE )

    print("matrix read")

    #ordering
    DForder = treeTraversal(tree)

    DForderDict = { }
    orderL_filtered=[]
    for sp in DForder:
        if sp in matrixheader :
            DForderDict[sp]=len(DForderDict)
            orderL_filtered.append(sp)
        else:
            print("Warning: {} present in tree, but absent from structure matrix. ignored.".format(sp))


    print( "The tree and the matrix have {intSize} leaves in common".format( intSize= len(DForderDict) ) )


    
    M = applyOrderToStructureMatrix( M , matrixheader , DForderDict )

    print("matrix reordered")

    ## reporting
    if not args.csv:
        writeStructureMatrix_bin( orderL_filtered , M , outFile  )
    else:
        OUT=open(outFile,'w')
        CSVmat(M,OUT,header=DForder)
        OUT.close()

    print("ordered matrix written in",outFile)    
