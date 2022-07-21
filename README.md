
This code implements methods and scripts to compute a structure matrix from a multiple sequence alignment.

The structure matrix describe degrees of closedness between sequence, and is computed using the methods from:

 * Structural analysis of biodiversity , Sirovitch 2010
 * A scalable method for analysis and display of DNA sequences 2009


> Note on the usage of the term "species" in this notice: we use here the term species in the a sense, as a modality to group different sequences together. In practice, one can group sequences by any arbitrary unit they define, or not group them at all.

## Requirements : 

You need to have the scipy librairies installed. Look here for instructions : https://www.scipy.org/install.html (like them, I suggest using Anaconda).
For the MPI script you need to have MPI enable on the machine and the python library mpi4py installed.

## Test :

The `testscripts/` folder contains a number of tests for the different scripts. They provide both a way to verify that your environment is compatible with the software and provide usage examples.

To test all scripts at once, use `sh testScripts/all_tests.sh` (NB: the test scripts need to have write access to the current folder in order to write logs and temporary files).

## Usage 

### Structure Matrix computation

Structure matrix computation is performed in 3 steps:
 1. computing the reference matrix from a (trimmed) alignment
 2. computing indicator vectors from an Ssum matrix and an alignment
 3. computing the structure matrix from the indicator vectors


If the dataset you manipulate is of a reasonnable size (say, less than 1000 sequences), then you may use the `src/pyKleeBarcode_MPI.py` script which wraps all three steps together.

However if your dataset is larger we recommend you use the scripts corresponding to each step:
 1. `src/pyKleeBarcode_computeRefMat_MPI.py`
 2. `src/pyKleeBarcode_computeIndicatorVector_MPI.py`
 3. `src/pyKleeBarcode_computeStructureMatrix.py`

The reason why are that each step require different amount of resources and benefit differently from parallelization.
Also, the results from the first two steps (Ssum matrix and Indicator vector computation) can be aggregated with results from previous run (see [later section](...)).


#### step 1: pyKleeBarcode_computeRefMat_MPI.py

Computes the reference matrix, a matrix representation of the diversity of a DNA sequence across a number of individuals or groups or individuals (typically, species)
given in an input alignment.


```sh
python pyKleeBarcode_computeSsum_MPI.py --help
```
```
usage: pyKleeBarcode_computeRefMat_MPI.py [-h] -i INPUTFILE -o OUTPUTFILE
                                          [-m MAX_SEQ_PER_SPECIES]
                                          [-f FIELD_DELIMITOR]
                                          [-s SPECIES_FIELD_INDEX]
                                          [-C CORRESPONDENCE_FILE]
                                          [--seed SEED]

Computes a reference matrix: a matrix representation of the diversity of a DNA
sequence across a number of individuals or groups or individuals (typically,
species) according to the definitions of "A scalable method for analysis and
display of DNA sequences" by Sirovitch et alii

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputFile INPUTFILE
                        input multiple sequence alignment in fasta format
                        (preferably, trimmed)
  -o OUTPUTFILE, --outputFile OUTPUTFILE
                        output file name for the reference matrix
  -m MAX_SEQ_PER_SPECIES, --max-seq-per-species MAX_SEQ_PER_SPECIES
                        maximum number of sequences kept per species. Default
                        : 3. Set to 0 to have no limit ; be careful thought,
                        as this parameter has a direct impact on speed.
  -f FIELD_DELIMITOR, --field-delimitor FIELD_DELIMITOR
                        field delimitor of the fasta sequence id lines.
                        default: "|"
  -s SPECIES_FIELD_INDEX, --species-field-index SPECIES_FIELD_INDEX
                        index (starting at 0) of the species name in the fasta
                        sequence id lines. default: 1
  -C CORRESPONDENCE_FILE, ---correspondence-file CORRESPONDENCE_FILE
                        name of a comma-delimited file containing
                        correspondence between sequence and species. Overrides
                        option -s when used. One sequence ID per line,
                        sequence id and species name separated by a semicolumn
                        (;).
  --seed SEED           random seed used when selecting a species sequences if
                        there is more than --max-seq-per-species. By default
                        is it created using time.
```

The `-m`, `-f`, `-s`, `-C`, and `--seed` options are designed for the case where a number of sequence should be aggregated prior to the computation of the matrix.
This can serve as an elegant solution to the overrepresentation of certain species (or whichever grouping you choose) where these are prevented from taking too much space in the created reference space (ie. the reference matrix), while still conserving all the information about the diversity of the sequence in these species.


**Aparte : grouping of sequences.**
By default, the grouping of sequences is governed by the field delimitor (`-f`) and field index (`-s`) options, whose default are, respectively `'|'` and `1`.
This means that the sequences will be grouped according to the second field (indexing begins at 0) in a | delimited fasta id line. For instance:

`>GBSP13299-19|Pseudocorynosoma anatarium|COI-5P|KX688147`

With the default values the group will be : `Pseudocorynosoma anatarium`.

With this form of sequence ids, one can avoid any kind of grouping (ie. keep all sequences as separate), one could set the `-s` option to `0` : the grouping would be the unique sequence id.

Alternatively, you can provide the associations between sequence and group in an external file using the `-C` option.
The file is expected to contain one association (ie, 1 sequence and 1 group) per line, separated by semicolons (`;`).

The sequence id used should correspond to the first field of the corresponding the fasta id line (field delimitor:`|`, can be changed with `-f`).

For instance if the correspondence file contains the line:

`GBSP13299-19;specieA`

Then a sequence with the fasta id line: `>GBSP13299-19|Pseudocorynosoma anatarium|COI-5P|KX688147` will be associated to `speciesA`.

---


As the title entails, this script can be parallelized using MPI if you have installed the proper library on your system.

**example usage:**
```
mpirun -np 4 python src/pyKleeBarcode_computeRefMat_MPI.py -i testData/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas -o 1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.Ssum 
```

> See the various scripts in `testScripts/` for more example usages.

#### step 2: pyKleeBarcode_computeIndicatorVector_MPI.py


```sh
python pyKleeBarcode_computeIndicatorVector_MPI.py --help
```
```
usage: pyKleeBarcode_computeIndicatorVector_MPI.py [-h] -i INPUTFILE -S
                                                   INPUTSSUMFILE -o OUTPUTFILE
                                                   [-m MAX_SEQ_PER_SPECIES]
                                                   [-f FIELD_DELIMITOR]
                                                   [-s SPECIES_FIELD_INDEX]
                                                   [-C CORRESPONDENCE_FILE]
                                                   [--seed SEED]

Computes the indicator vectors of the DNA sequences of a number of individuals
or groups or individuals (typically, species), from the diversity represented
in a given Ssum matrix according to the definitions of "A scalable method for
analysis and display of DNA sequences" by Sirovitch et alii

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputFile INPUTFILE
                        input multiple sequence alignment in fasta format
                        (preferably, trimmed)
  -S INPUTSSUMFILE, --inputSsumFile INPUTSSUMFILE
                        input Ssum matrix (expected: binary format as produced
                        by the pyKleeBarcode_computeSsum_MPI.py script)
  -o OUTPUTFILE, --outputFile OUTPUTFILE
                        output file name for the Ssum matrix
  -m MAX_SEQ_PER_SPECIES, --max-seq-per-species MAX_SEQ_PER_SPECIES
                        maximum number of sequences kept per species. Default
                        : 3. Set to 0 to have no limit ; be careful thought,
                        as this parameter has a direct impact on speed.
  -f FIELD_DELIMITOR, --field-delimitor FIELD_DELIMITOR
                        field delimitor of the fasta sequence id lines. default: "|"
  -s SPECIES_FIELD_INDEX, --species-field-index SPECIES_FIELD_INDEX
                        index (starting at 0) of the species name in the fasta
                        sequence id lines. default: 1
  -C CORRESPONDENCE_FILE, ---correspondence-file CORRESPONDENCE_FILE
                        name of a comma-delimited file containing
                        correspondence between sequence and species. Overrides
                        option -s when used. One sequence ID per line,
                        sequence id and species name separated by a semicolumn
                        (;).
  --seed SEED           random seed used when selecting a species sequences if
                        there is more than --max-seq-per-species. By default
                        is it created using time.
```

The `-m`, `-f`, `-s`, `-C` and `--seed` options are designed for the case where a number of sequence should be aggregated prior to the computation of the indicator vector.
In general they should be set to(although this is not mandatory), the same values as the one given during the computation of the reference matrix in the previous step.

**example usage:**
```
mpirun -np 4 python src/pyKleeBarcode_computeIndicatorVector_MPI.py -i testData/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas -S 1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.Ssum -o 1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.indicatorVectors.csv
```

#### step 3: pyKleeBarcode_computeStructureMatrix.py

```sh
python src/pyKleeBarcode_computeStructureMatrix.py --help
```
```
usage: pyKleeBarcode_computeStructureMatrix.py [-h] -i INPUTFILE -o OUTPUTFILE

Computes a structure matrix from a set of indicator vectors according to the
definitions of "A scalable method for analysis and display of DNA sequences"
by Sirovitch et alii

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputFile INPUTFILE
                        input indicator vectors in csv format (one seq/species
                        per line)
  -o OUTPUTFILE, --outputFile OUTPUTFILE
                        output file name for the structure matrix
```

**example usage:**
```
python src/pyKleeBarcode_computeStructureMatrix.py -i testData/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.indicatorVectors.csv -o testData/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.stMat.bin
```

#### all steps in one go (recommended for small dataset only)

```sh
python src/pyKleeBarcode_MPI.py --help
```
```
usage: pyKleeBarcode_MPI.py [-h] -i INPUTFILE -o OUTPUTFILE
                            [-m MAX_SEQ_PER_SPECIES] [-f FIELD_DELIMITOR]
                            [-s SPECIES_FIELD_INDEX] [-C CORRESPONDENCE_FILE]
                            [--seed SEED]

Computes a structure matrix between sets of DNA sequences (typically grouped
by species) according to the definitions of "A scalable method for analysis
and display of DNA sequences" by Sirovitch et alii

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputFile INPUTFILE
                        input multiple sequence alignment in fasta format
                        (preferably, trimmed)
  -o OUTPUTFILE, --outputFile OUTPUTFILE
                        output file name for the structural matrix
  -m MAX_SEQ_PER_SPECIES, --max-seq-per-species MAX_SEQ_PER_SPECIES
                        maximum number of sequences kept per species. Default
                        : 3. Set to 0 to have no limit ; be careful thought,
                        as this parameter has a direct impact on speed.
  -f FIELD_DELIMITOR, --field-delimitor FIELD_DELIMITOR
                        field delimitor of the fasta sequence id lines.
                        default: "|"
  -s SPECIES_FIELD_INDEX, --species-field-index SPECIES_FIELD_INDEX
                        index (starting at 0) of the species name in the fasta
                        sequence id lines. default: 1
  -C CORRESPONDENCE_FILE, ---correspondence-file CORRESPONDENCE_FILE
                        name of a comma-delimited file containing
                        correspondence between sequence and species. Overrides
                        option -s when used. One sequence ID per line,
                        sequence id and species name separated by a semicolumn
                        (;).
  --seed SEED           random seed used when selecting a species sequences if
                        there is more than --max-seq-per-species. By default
                        is it created using time.

```

The different options follow the same logic as the separate steps.

**example usage:**
```
python src/pyKleeBarcode_MPI.py --inputFile testData/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas -o 1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.stMat.bin
```

### Handling large datasets & updating an existing system

`pyKleeBarcode` offers utilities to allows massive distribution of the computations needed to compute indicator vectors and a structure matrix,
as well as the possibility to update an existing structure matrix with new data without having to re-compute everything from scratch.

#### Ssum matrix merging : `src/mergeReferenceMat.py`

Provided that they were computed on an independent set of sequences AND that no grouping/species are shared, then two reference matrices can be merged (at a very reasonnable computational cost, because the only required operation is a simple matrix addition).

This allows the computation of the reference matrix to be subdivided in any number of subtasks (as long as the species/groups are part of the same subtask), by simply splitting the input alignment, and subsequently merging the resulting matrices.

This also allows one to update a previously obtained matrix with a new set of sequences, as long as these new sequences only belong to species/groups hitherto absent from the reference matrix.

#### Indicator vector merging

The computation of indicator vectors can be made completely independently between each sequence, provided the script are given the same Ssum matrix.
Thus, after the Ssum matrix computation one can split the input alignment however they see fit in order to distribute the computation of the indicator vectors as needed.

Indicator vector files are simple `.csv` files which can be merged by concatenation, for example with a simple `cat` command call.

Given this, updating a set of indicator vectors with new sequences is fairly simple. One would have to compute the indicator vectors of the new sequences (using the same Ssum matrix as the one used for the already existing indicator vectors), and then concatenate the created indicator vector file to the existing one.


#### Structure matrix updating

While the previous steps could be split fairly arbitrarily with at little computational cost, the structure matrix computation matrix computation is a more complicated affair, because it must look at all pairs of indicator vectors.

Nevertheless, we have devised an algorithm that allows the update of an existing structure matrix with new indicator vector (thus, new sequences). 
Because structure matrix files quickly become large (nb: it grows quadratically with the number of sequences) the structure matrix file format has been thought to allow and updating operation which does not have to re-write the whole file but only append the new information to it (long story short : we represent the lower triangular portion of the matrix only).



```sh
python src/pyKleeBarcode_updateStructureMatrix.py --help
```
```
usage: pyKleeBarcode_updateStructureMatrix.py [-h] -r REFERENCE_VECTORS -n
                                              NEW_VECTORS -m STRUCTURE_MATRIX

Updates a structure matrix from a set of reference indicator vectors (the
sequences already in the structure matrix) amd a set of new structure vectors
(to add to the structure matrix) according to the definitions of "A scalable
method for analysis and display of DNA sequences" by Sirovitch et alii .
IMPORTANT: this script assumes sequences are ordered similarly between the
reference indicator vectors and structure matrix.

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE_VECTORS, --reference-vectors REFERENCE_VECTORS
                        input reference indicator vectors in csv format (one
                        seq/species per line)
  -n NEW_VECTORS, --new-vectors NEW_VECTORS
                        input new indicator vectors in csv format (one
                        seq/species per line)
  -m STRUCTURE_MATRIX, --structure-matrix STRUCTURE_MATRIX
                        structure matrix file to update
```



### Other scripts

#### Alignment utilities


 * `trimAlignment.py` : to trim an alignment
 * `sequenceAveragePairWiseComputation.py` : to compute average pairwise distances among and between species in an alignment
 * `alignmentAveragePairWiseComputation.py` : computes average pairwise difference among sequences in an alignment (often termed Pi in population genetics) as well as other diversity metrics.
 * `restrictAlnToClade.py` : Extract the sequences of specific clades from an alignment.
 * `convertStructureMat_CSV_to_bin.py` : utility to convert a structure matrix between the binary and a csv format useful for debugging


#### Representation of structure matrix

* `order_structureMatrix_with_tree.py` : reorders the rows inside a structure matrix according to a tree in newick format
* `order_structureMatrix_with_list.py` : reorders the rows inside a structure matrix according to a list (file with one id per line)
* `heatmap_from_structureMatrix.py` : produces a klee diagram image from a structure matrix
* `taxon_annotated_heatmap_from_structureMatrix.py` : produces a klee diagram image from a structure matrix, with taxonomic annotations from a file.
    see an example of annotation file in `testData/Animal.taxons.ok_and_resolved.txt`, which should work for most animal taxons following the ncbi taxonomy.
