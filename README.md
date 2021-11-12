# PhyloAnalyser
PhyloAnalyser by erpeg
### Tool to perform phylogenetic analysis of chosen organisms.
#### Description
Main objective of this tool was to create software, that incorporates multiple fast computing softwares into one, enabling fast analysis.

PhyloAnalyser workflow consists of a few steps:
1. Retrieve data from Uniprot based on taxon identifiers
2. Process .fasta files so they are suitable for clustering process (remove too short/long sequences)
3. Perform clustering process using MMSeq2 easy-cluster option
4. Process clustering output (remove too short/long clusters, obtain one-to-one clusters)
5. Align sequences using MAFFT (through MuscleCommandLine)
6. Infer trees using FastTree and save output in two distinct directories (all for all trees and filtered for one-to-one clusters)
7. Based on bootstrap values obtained from FastTree select trees with bootstrap values over specified threshold
8. Compute consensus trees based on inferred trees using Bio.Phylo.Consensus majority_consensus()
 function, with strictness 0.3 - can be changed
9. Compute supertree using fxxxxxxxx2 software
10. Calcualte Robinson-Foulds for generated species trees and compare with known_topology.nwk and generate RF matrix
11. Generate pictures of each tree in .png and .txt (ascii) format.

For 9k clusters with fungi as input, whole analysis takes around 20 minutes.
 
 #### Prerequisities
In order to function properly this application requires a few packages and tools preinstaled (mainly inside conda environment).

* Conda environment with >= Python 3.7
* WKHTLMTOPDF tool - `conda install -c bioconda wkhtmltopdf`
* Biopython - `conda install -c conda-forge biopython`
* SeqKit - `conda install -c bioconda seqkit`
* MmSeq2 - `conda install -c conda-forge -c bioconda mmseqs2`
* MAFFT - `conda install -c bioconda mafft`
* FastTree - `conda install -c bioconda fasttree`
* fxxxxxxxx2 in folder where analysis is going to be performed (added to PATH fxxxxxxx2 directory) - title of software is  filled with "x" due to the fact this tool is not officialy released
* other Python libraries requirements can be found in `requirements.txt`

Furthermore in order to get proper analysis, in place where script is run newick file with known toplogoy named as "known_topology.nwk" has to be present
 
 #### Input file format
 As input file .txt file is accepted with NCBI taxon id in each line
 ```shell script
765915
578462
1314771
747725
```
File is passed into parser by using -i parameter or --input_filename

#### Output
As output many directories are received in below structure.
Names of each directories are quite self-explenatory.

"ft" directory has been created in order to enable future implementation of other methods used to inferr trees.

```shell script
.
├── alignment
│   ├── all
│   └── filtered
├── archive
├── cluster
│   ├── all
│   └── filtered
├── seq
│   └── organisms
│       ├── concatenated
│       ├── renamed
│       └── stock_ids
├── tree
│   ├── calculated
│   ├── concatenated
│   └── ft
│       ├── all
│       │   ├── boot_done
│       │   │   └── supertree
│       │   ├── boot_nondone
│       │   │   └── supertree
│       │   ├── bootstrapped
│       │   └── inferred
│       └── filtered
│           ├── boot_done
│           │   ├── consensus
│           │   └── supertree
│           ├── boot_nondone
│           │   ├── consensus
│           │   └── supertree
│           ├── bootstrapped
│           └── inferred
└── tree_figures

```
In output directory file all_trees_rf.png represents matrix of Robinson-Foulds calculated between each of trees

In tree_figures are pictures of computed trees and computed trees are in tree/calculated directory.

#### Example usage

```shell script
./app.py -i taxid_to_analyse.txt -o example_directory 
```
or
```shell script
./app.py --input_filename taxid_to_analyse.txt --output_dir last_test 
```
