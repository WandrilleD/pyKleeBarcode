import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict,Counter
from pyKleeBarcode_utils import readStructureMatrix_bin

class taxTree:
    def __init__( self ):
        self.root = None
        self.nodes = defaultdict(set)
        self.nodeParent = dict()
        self.genuses = set()

        self.splittingNodes = set()
        self.nodeSize = {}
        self.nodeHeight = {}
        #self.effectiveRoot = None


    def addPhylum(self , phylumLineage , isGenus):
        """ phylumLineage is an ordered list of taxon name, starting at the root and finishing at the phylum name,
        isGenus is a boolean. """

        for i,g in enumerate(phylumLineage):
            if i == 0:
                if self.root is None:
                    self.root = g 
                elif self.root != g:
                    print("problem: new phylum {} has a root ({}) different from that of the current tree ({})".format(g,phylumLineage[-1],self.root))
                    exit(1)
            else:
                self.nodes[phylumLineage[i-1]].add( g )
                if len(self.nodes[phylumLineage[i-1]])>1:
                    self.splittingNodes.add( phylumLineage[i-1] )

                ## adding a check for consistency in the tree structure. This can happen if the phyla names are not unique identifiers ...
                if not g in self.nodeParent:
                    self.nodeParent[g] = phylumLineage[i-1]
                elif self.nodeParent[g] != phylumLineage[i-1]:
                    print("warning : taxon",g,"appears with different parents")

                #if g == "Acanthocephala":
                #    print( g , self.nodeParent[g] , phylumLineage[i-1] )

        if isGenus:
            self.genuses.add( phylumLineage[-1] )
        return

    def sortNodeOrder(self):
        ''' alphabetically sorts the children of all nodes in order to stabililize plotting across datasets '''
        for n in self.nodes.keys():
            self.nodes[n] = list( self.nodes[n] )
            self.nodes[n].sort()



    def setNodeSizes(self , current = None,  genusSizes = {} ):

        if current in self.genuses:
            self.nodeSize[current] = genusSizes.get(current , 1 )## default size is 1
            return 

        if current is None:
            current = self.root

        size = 0
        for child in self.nodes[current]:
            if not child in self.nodeSize:
                self.setNodeSizes( child , genusSizes )
            size += self.nodeSize[child] 

        self.nodeSize[current]=size

        return

    def setNodeHeights(self , current = None,  ignoreNonSplitting = False ):
        """ if ignoreNonSplitting is True, taxons with a single child won't increase height 
            updates self.nodeHeight
        """

        if current in self.genuses:
            self.nodeHeight[current] = 0
            return 

        if current is None:
            current = self.root

        height = 0
        for child in self.nodes[current]:
            if not child in self.nodeHeight:
                self.setNodeHeights( child , ignoreNonSplitting )
            height = max( height , self.nodeHeight[child] )

        if ( not ignoreNonSplitting  ) or len( self.nodes[current] )>1:
            height +=1      
        self.nodeHeight[current]=height

        return 

    def getPhylumBrackets( self , current=None , currentMin=0 ):
        """ returns a dictionnary associating each phylum to its "bracket",
            which is tuple of integer representing the space taken by the phylum 
            in a DFT 2D representation 
         """

        if current in self.genuses:
            return {current:(currentMin, currentMin + self.nodeSize.get(current,1)-1 )}


        if current is None:
            current = self.root

        #print(current, self.nodes[current])

        phylumBrackets = {}

        localMin = currentMin
        for child in self.nodes[current]:
            phylumBrackets.update( self.getPhylumBrackets( child , localMin ) )
            localMin = phylumBrackets[child][1]+1 #incrementing the min


        phylumBrackets[current]=( currentMin , localMin-1  ) # -1 because there is no padding between the end of the parent bracket and the end of its last child's bracket


        return phylumBrackets
    





def readTaxonomy( inFile ):
    """ 
        Reads a file describing taxons such as the ones output by ete3 ncbiquery .
        
        Returns:
            (tuple):
                (dict) : keys are taxon ids, values are lists of taxon ids describing the lineage of the taxon (which is the last of the list) 
                (dict) : keys are taxon ids, values are the name of the taxon id
                (dict) : keys are taxon ids, values are the rank of the taxon id

    """
    phylaLineageDict = {}
    phylaNameDict = {}
    phylaRankDict = {}

    IN = open(inFile , 'r')
    l = IN.readline()
    if l.startswith("#  Taxid"):#header
        ## Taxid Sci.Name   Rank    Named Lineage   Taxid Lineage
        l = IN.readline()
    while l != "":
        sl = l.strip().split('\t')
        #print(sl)
        taxid = sl[0]
        name = sl[1]
        rank = sl[2]
        lineage = sl[4].split(',')
        lineageNames = sl[3].split(',')
        
        phylaLineageDict[taxid] = lineage
        phylaRankDict[taxid] = rank

        for i in range(len(lineage)):
            phylaNameDict[ lineage[i] ] = lineageNames[i]


        l = IN.readline()

    IN.close()

    return phylaLineageDict , phylaNameDict , phylaRankDict



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

def species2genus(spName , sep=" "):
    ''' dirt simple, just looks at the 1st work of the species name '''
    return spName.strip(sep).partition(sep)[0]

def getSpeciesOrderFromBrackets( representedGenus , phylaBrackets ):
    """
        Takes:
            - representedGenus (dict) : keys are taxid, values are list of species
            - phylaBrackets (dict) : keys are taxid, values are tuple of integer representing a 1D mapping of the genera
                                     NB : brackets start at 0 and move forward in increments of 1

        Returns:
            (dict): keys are species name, values are the index of the species in the order proposed by the brackets
    """
    speciesOrder = {}
    for genus,spList in representedGenus.items():

        currentPosition = phylaBrackets[ genus ][0]
        for sp in spList:
            speciesOrder[sp] = currentPosition
            currentPosition+=1
    return speciesOrder


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
            newOrder[sp] = len(newOrder)
    
    newMatrix = np.zeros( structureMatrix.shape )

    for i,sp1 in enumerate(currentOrder):


        new1 = newOrder.get( sp1 )

        for j,sp2 in enumerate(currentOrder):
            
            new2 = newOrder[ sp2 ]

            newMatrix[new1][new2] = structureMatrix[i][j]

    return newMatrix



## heatmap functions heavily inspired from
## https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/image_annotated_heatmap.html

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    if len(col_labels)>0:
        ax.set_xticks(np.arange(data.shape[1]))
        ax.set_yticks(np.arange(data.shape[0]))
        # ... and label them with the respective list entries.
        ax.set_xticklabels(col_labels)
        ax.set_yticklabels(row_labels)


        # Let the horizontal axes labeling appear on top.
        ax.tick_params(top=True, bottom=False,
                       labeltop=True, labelbottom=False)
    
        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
                 rotation_mode="anchor")
        ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
        ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
        ax.tick_params(which="minor", bottom=False, left=False)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    # Turn spines off and create white grid.
    #for edge, spine in ax.spines.items():
    #    spine.set_visible(False)

    #ax.grid(which="minor", color="w", linestyle='-', linewidth=3)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if data[i, j] is None or np.isnan(data[i, j]) :
                continue
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def plotBrackets( phylaBrackets , taxaTree ):

    for n in phylaBrackets.keys():
    
        x = - taxaTree.nodeHeight[n]
        Ys = phylaBrackets[n]
    
        if n in taxaTree.splittingNodes :
            plt.plot( [x,x] , Ys ,color="black")
            plt.plot( [x,0] , [Ys[0],Ys[0]] ,color="black")
            plt.plot( [x,0] , [Ys[1],Ys[1]] ,color="black")
    
        elif n in taxaTree.genuses:
            plt.plot( [x,x+0.5] , [Ys[0],Ys[0]] )
            plt.plot( [x,x+0.5] , [Ys[1],Ys[1]] )
            plt.plot( [x,x] , Ys )

    
    print("rendering ... ")
    plt.show()


def makeBracketBox( im , bracket , bracketLabel = "" , boxArgs={ "color":"black", "linewidth":1 } , textArgs = {"fontsize":6,"color":"white"} ):
    """
        Takes:
            - im (AxesImage) : The AxesImage to be labeled.
            - bracket (tuple) : couple of integer determining the bounding box of the bracket
            - bracketLabel="" (str)  : label to plot, leave empty to avoid plotting
            - boxArgs={ "color":"black", "linewidth"=1 } (dict) : arguments to give to the function plotting the bracket box
            - textArgs = {"fontsize":6,"color":"white"} : arguments to give to the function plotting the bracket label

    """

    boxes = []
    labels = []

    i = bracket[0] -0.5
    j = bracket[1] +0.5

    boxes.append( im.axes.plot( [i,j,j,i,i] , [i,i,j,j,i] , **boxArgs ) )
    
    if bracketLabel != "":
        if not "horizontalalignment" in textArgs:
            textArgs["horizontalalignment"]="center"
        if not "verticalalignment" in textArgs:
            textArgs["verticalalignment"]="top"
    
        labels.append( im.axes.text(i, j, bracketLabel , **textArgs) )


    return boxes,labels



if __name__ == "__main__" :

    import sys
    import argparse

    parser = argparse.ArgumentParser(
                description="""Produces the klee diagram from a structure matrix
                (heatmap-like visualisation of a sequence correlation measure), organized according to a provided taxonomic structure.
                Note that for readability purpose the names of groups/species will only be written if there are less than 100.""")
    parser.add_argument('-m','--inputMatrix', type=str, required=True,
             help='input structure matrix obtained from an alignment')
    parser.add_argument('-t','--inputTaxons', type=str, required=True,
             help='a file containing taxons information')

    parser.add_argument('-o','--outputFile', type=str, default=None,
             help='output file name for the image. Will be printed to screen if omitted.')

    parser.add_argument('--sep', type=str, default=" ",
             help='separator between genus and species name in the structure matrix (default:space).')
    
    parser.add_argument('-s','--size', type=float, default=640,
             help='height of the generated image (in pixels).')

    parser.add_argument('-S','--seq2species', type=str,default='',
             help='a file containing a correspondance between sequence ids (as they appear in the structure matrix), and species or genus name.\nOne sequence ID per line, sequence id and species name separated by a semicolumn (;).\nIf this option is omitted the sequence ids are presumed to reflect a species or genus name.')

    parser.add_argument('--label-with-species', action='store_true',
             help='label the image with the species name (given with the -S option)')


    parser.add_argument('--csv', action='store_true',
             help='specify that the structure matrix in in csv format (default expects a binary format)')

    args = parser.parse_args()

    STRMAT_FILE = args.inputMatrix
    TAXON_FILE = args.inputTaxons
    outFile= args.outputFile

    genusSep = args.sep

    ANNOTATION_THRESHOLD=100

    # reading the matrix
    smreading = readStructureMatrix_bin
    if args.csv :
        smreading = readStructureMatrix

    sqList , structureMatrix = smreading( STRMAT_FILE )

    print( "read the structure matrix containing",len(sqList),"elements." )
    # reading the taxons
    phylaLineageDict , phylaNameDict , phylaRankDict = readTaxonomy( TAXON_FILE )

    # tring to get a unique mapping between the names of the taxons   and taxids.
    # if this cannot be achieved, the script will fail (that could be circumvented by asking the user to explicitely provide genus taxids for each species)
    genusName2taxid = {}
    for taxid,rank in phylaRankDict.items():
        if rank == "genus":
            # this restricts us to the genus level only, which is necessary to deduce taxid from species name.
            name = phylaNameDict[taxid]
            if name in genusName2taxid:
                print("ERROR : the genus",name,"appear with different taxids :",taxid,genusName2taxid[name])
                exit(1)
            genusName2taxid[name] = taxid

    # handling case where a file giving correspondance between seqID and species name was given
    seq2species = {s:s for s in sqList}
    if not args.seq2species == '':
        try:
            with open(args.seq2species,'r') as IN:
                for l in IN:
                    sl = l.strip().split(';')
                    seq2species[sl[0]] = sl[1]
        except :
            print("ERROR : error while reading the sequence to species name correspondance file :",args.seq2species)
            exit(1)

    




    # creating a taxon tree from the genera that are present in the structure matrix
    representedGenus = defaultdict(list)
    unknownGenera = defaultdict(list)
    for sq in sqList:
        genusName = species2genus( seq2species[sq] , sep = genusSep)
        taxid = genusName2taxid.get(genusName, None)

        if taxid is None:
            unknownGenera[genusName].append( sq )
        else:
            representedGenus[ taxid ].append( sq )
    if len(unknownGenera)>0:
        print( "WARNING : the structure matrix contains unknown {} unknown genera.".format(len(unknownGenera)) )
        print('"'+'" "'.join([x for x in unknownGenera.keys()])+'"')
        if args.seq2species != '':
            print("resp. associated to species:", '"'+'" "'.join([ seq2species[ x ] for x in unknownGenera.keys()])+'"' )
        print('"'+'" "'.join([x for x in unknownGenera.keys()])+'"')
        print('The associated sequences will be pushed at the bottom/right of the plot')
        #exit(1)

    taxaTree = taxTree()
    for taxid in representedGenus.keys():
        
        lineage = phylaLineageDict[taxid]

        taxaTree.addPhylum( lineage , isGenus=True )

    print("taxa tree built")
    print("number of splits {}".format(len(taxaTree.splittingNodes)))
    print("number of genera {}".format(len(taxaTree.genuses)))

    taxaTree.sortNodeOrder()
    taxaTree.setNodeHeights(ignoreNonSplitting=True)
    taxaTree.setNodeSizes( genusSizes = {taxid:len(species) for taxid,species in representedGenus.items()} )
    #print(representedGenus)

    phylaBrackets = taxaTree.getPhylumBrackets(  )
    print("phyla brackets computed")



    ##simple plot of the brackets to check their validity
    #plotBrackets( phylaBrackets , taxaTree )

    # get the order of species from the brackets
    speciesOrderDict = getSpeciesOrderFromBrackets( representedGenus , phylaBrackets )

    #apply the order
    structureMatrix = applyOrderToStructureMatrix( structureMatrix , sqList , speciesOrderDict )

    for sp,i in speciesOrderDict.items():
        sqList[i]=sp

    print("structure matrix ordered, now plotting")
    ## in plt.subplot the figsize is in inches, and default dpi is 100
    fsize =  args.size / 100

    fig, ax = plt.subplots(figsize=(fsize*(5/4), fsize ) ) 
    #fig.patch.set_facecolor('grey')
    ## adding 1/4th additional width to account for the scale


    ## "masking" the lower triangle of the matrix to make room for the taxon annotation
    i=0
    while i < structureMatrix.shape[0]:
        j=0
        while j <= i:
            structureMatrix[i,j]=np.nan
            j+=1
        i+=1


    labels = sqList
    if args.label_with_species :
        labels = [seq2species[s] for s in sqList]



    if len(sqList)>ANNOTATION_THRESHOLD:
        labels = []

    colorMap = 'nipy_spectral' # lots of contrast, in particular close to 0
    #colorMap = 'gist_ncar' # lots of contrast, in particular close to 0
    #colorMap = 'gist_rainbow' # like jet, with an addotionnal color
    #colorMap = 'jet'
    #colorMap = 'cool'
    #colorMap = 'Wistia'
    #colorMap = 'viridis'
    #colorMap = 'PuOr'
    #colorMap = 'PiYG'
    #colorMap = 'YlGn'

    im, cbar = heatmap( structureMatrix , labels, labels, ax=ax,
                   cmap=colorMap, cbarlabel="indicator vector correlation")
    
    #print ing values, useful for debugging, but cluttering visualization otherwise
    #if len(spList)<=ANNOTATION_THRESHOLD:
    #    texts = annotate_heatmap(im,structureMatrix, valfmt="{x:.2f}" , fontsize=6)

    ## adding the bounding boxes of taxa
    for n in phylaBrackets.keys():
        if n in taxaTree.splittingNodes  or n in taxaTree.genuses : # restricting to nodes that split
            textArgs = {"fontsize":8,
                        "color":"black",
                        'weight':'bold',
                        'verticalalignment':"center",
                        'horizontalalignment':"center",
                        "rotation":-45
                        }
            if phylaRankDict.get(n) == "genus":
                textArgs['fontsize']=6
    
            makeBracketBox( im , phylaBrackets[n] , bracketLabel = phylaNameDict[ n ] , boxArgs={ "color":"grey", "linewidth":1 } , textArgs = textArgs )

    fig.tight_layout()

    if outFile is None:
        plt.show()
    else:
        plt.savefig(outFile)