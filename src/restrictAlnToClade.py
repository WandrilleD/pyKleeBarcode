import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from collections import defaultdict,Counter

from pyKleeBarcode_utils import iterFasta, getSequenceSpecies

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

    def extractGeneraFromClade( self , clade ):


        if clade in self.genuses:
            return [clade]

        if not clade in self.nodes :
            return None

        genera = []

        for child in self.nodes[clade]:
            genera += self.extractGeneraFromClade( child )

        return genera




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

def species2genus(spName , sep=" "):
    ''' dirt simple, just looks at the 1st work of the species name '''
    return spName.strip(sep).partition(sep)[0]


if __name__ == "__main__" :

    import sys
    import argparse

    parser = argparse.ArgumentParser(
                description="""Extract the sequences of specific clades from an alignment.""")
    parser.add_argument('-i','--inputAlignment', type=str, required=True,
             help='input structure alignment')
    parser.add_argument('-t','--inputTaxons', type=str, required=True,
             help='a file containing taxons information')

    parser.add_argument('-C','--cladeList', type=str, required=True,
             help='comma separated list of clades (no spaces)')

    parser.add_argument('-o','--outputFile', type=str, default=None,
             help='output aignment file name.')

    parser.add_argument('--sep', type=str, default=" ",
             help='separator between genus and species name in the structure matrix (default:space).')
    
    parser.add_argument('-f','--field-delimitor', type=str, default="|",
             help='field delimitor of the fasta sequence id lines')

    parser.add_argument('-s','--species-field-index', type=int, default=1,
             help='index (starting at 0) of the species name in the fasta sequence id lines')


    args = parser.parse_args()

    INPUT_FILE = args.inputAlignment
    TAXON_FILE = args.inputTaxons
    outFile= args.outputFile

    genusSep = args.sep
    FASTA_ID_SEPARATOR = args.field_delimitor
    SPECIES_NAME_INDEX = args.species_field_index

    selectedCladeList = args.cladeList.split(',')

    # reading the taxons
    phylaLineageDict , phylaNameDict , phylaRankDict = readTaxonomy( TAXON_FILE )

    reversePhylaNameDict=defaultdict(list)
    for taxid,name in phylaNameDict.items():
        reversePhylaNameDict[name].append( taxid )


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

    # creating the taxon tree from the genera that are present in the structure matrix
    taxaTree = taxTree()
    for taxid in phylaLineageDict.keys():
        
        lineage = phylaLineageDict[taxid]

        taxaTree.addPhylum( lineage , isGenus=True )

    print("taxa tree built")


    selectedGenera = set()
    for clade in selectedCladeList :
        taxidList = reversePhylaNameDict[clade]

        if len(taxidList)>1 :
            print( "Warning : clade",clade,"corresponds to",len(taxidList),'taxa. All will be included.' )
        for taxid in taxidList:

            childGenera = taxaTree.extractGeneraFromClade(  taxid )
            if childGenera is None:
                print("!Error! clade",clade ,'not found among the taxa tree.')
                exit(1)
            selectedGenera.update( set( [ phylaNameDict[genId] for genId in childGenera] ) )

    print("number of selected genera",len(selectedGenera))

    IN = open(INPUT_FILE, 'r')
    OUT = open( outFile , 'w' )
    for k,s in iterFasta(IN):
        sp = getSequenceSpecies(k, FASTA_ID_SEPARATOR , SPECIES_NAME_INDEX )
        genus = species2genus(sp , sep=genusSep)

        if genus in selectedGenera:
            print( '>' + k  , file = OUT)
            print( s , file = OUT)
    
    IN.close()
    OUT.close()


