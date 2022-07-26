# LiftoffTools

LiftoffTools is a toolkit to compare genes lifted between genome assemblies. Specifically it is designed to compare genes lifted over using [Liftoff](https://github.com/agshumate/Liftoff) although it is also compatible with other lift-over tools such as UCSC liftOver as long as the feature IDs are the same. LiftoffTools provides 3 different modules. The first identifies variants in protein-coding genes and their effects on the gene. The second compares the gene synteny, and the third clusters genes into groups of paralogs to evaluate gene copy number gain and loss. The input for all modules is the reference genome assembly (FASTA), target genome assembly (FASTA), reference annotation (GFF/GTF), and target annotation (GFF/GTF). 


## Installation

```
git clone https://github.com/agshumate/LiftoffTools liftofftools 
cd liftofftools
python setup.py install
```
For the synteny module, you must also have [MMSeqs2](https://github.com/soedinglab/MMseqs2) installed and in your path 

## Usage

```
usage: liftofftools [-h] -r R -t T -rg GFF/GTF or DB -tg GFF/GTF or DB [-c]
                    [-f F] [-infer-genes] [-dir DIR] [-force]
                    [-mmseqs_path MMSEQS_PATH] [-mmseqs_params =STR]
                    [-edit-distance] [-V]
                    {clusters,variants,synteny,all}

Compare gene annotations across genome assemblies

Subcommands:
  {clusters,variants,synteny,all}

optional arguments:
  -h, --help            show this help message and exit
  -r R                  reference fasta
  -t T                  target fasta
  -rg GFF/GTF or DB     reference annotation file to lift over in GFF or GTF
                        format or gffutils database created in previous
                        liftoff or liftofftools run
  -tg GFF/GTF or DB     target annotation file to lift over in GFF or GTF
                        format or gffutils databased created in previous
                        liftoff or liftofftools run
  -c                    analyze protein coding gene clusters only
  -f F                  text file with additional feature types besides genes
                        to analyze
  -infer-genes
  -dir DIR              output directory
  -force                force overwrite of output/intermediate files in -dir
  -V, --version         show program version

clusters arguments:
  -mmseqs_path MMSEQS_PATH
                        mmseqs path if not in working directory or PATH
  -mmseqs_params =STR   space delimited list of additional mmseqs parameters.
                        Default="--min-seq-id 0.9 -c 0.9"

synteny arguments:
  -edit-distance        calculate edit distance between reference gene order
                        and target gene order
   -r-sort R_SORT        txt file with the order of the reference chromosomes
                        to be plotted on the x-axis
  -t-sort T_SORT        txt file with the order of the target chromosomes to
                        be plotted on the y-axis
```
### Output Directory
By default, all output files will be written to a directory within the current working directory called liftofftools_output. This can be changed with the -dir parameter. By default LiftoffTools will not overwrite files in the output directory. To enable overwrite, use the -force option. 


## Modules
To run all three modules, use the following command:

```
liftofftools all -r <reference.fa> -t <target.fa> -rg <reference.gff3> -tg <target.gff3>
```
Each module can also be run separately. See the information below for a description of each module and the command to run it. 

### Variants Module
The variants module can be run with the following command:

```
liftofftools variants -r <reference.fa> -t <target.fa> -rg <reference.gff3> -tg <target.gff3>
```

The variants module calculates the sequence identity between transcripts in the reference genome and the corresponding transcript in the target genome and for protein-coding genes, identifies variants and their effect on the gene. The effects we look for are defined below.

```
synonymous - A point mutation in the target transcript that does not change the amino acid sequence.
nonsynonymous - A point mutation in the target transcript that changes the amino acid sequence.
inframe deletion - A deletion in the target transcript sequence that is of some length divisible by 3.
inframe insertion - An insertion in the target transcript sequence that is of some length divisible by 3.
start codon loss - A point mutation in the start codon of the target transcript sequence.
5' truncation - A deletion of the 5' end of the target transcript. 
3' truncation - A deletion of the 3' end of the target transcript. 
frameshift - A deletion or insertion in the target transcript that is of some length not divisible by 3. 
stop codon gained - A point mutation in the target transcript that introduces a premature stop codon. 
```

If there is more than one variant, we report only the most severe. For example, if a transcript has a synonymous mutation and a frameshift mutation, we output ‘frameshift’ for that transcript as this would be more disruptive to gene function.

#### Variants Output
The output of the variants module is a file called 'variant_effects' located in the output directory (liftofftools_output by default). This is a tab separated file with the following fields.

|Col|Description                                                                       |
|--:|:---------------------------------------------------------------------------------|
|1  |Reference transcript ID                                                           |
|2  |Target transcript ID                                                              |
|3  |Sequence identity at the DNA level (0-1.0)                                        |
|4  |Sequence identity at the amino acid level (0-1.0). 'NA' for non-coding transcripts|
|5  |Variant effect (mose severe variant). 'NA' for non-coding transcripts             |



### Synteny Module

The synteny module can be run with the following command:

```
liftofftools synteny -r <reference.fa> -t <target.fa> -rg <reference.gff3> -tg <target.gff3>
```

This module compares the gene order in the reference annotation to the order in the target annotation. The genes present in both annotations are sorted first by chromosome and then by start coordinate in each annotation. 

#### Synteny Output

##### Output File 1
The first output of the synteny module is a dot plot called gene_order_plot.pdf located in the output directory (liftofftools_output by default).
In this file each gene is a point on a 2D plot where the x-coordinate is the ordinal position (e.g., 1st, 2nd, 3rd, etc.) in the reference genome and the y-coordinate is the ordinal position in the target genome. The color of the point corresponds to the sequence identity between reference gene and the target gene where green indicates higher identity and red indicates lower identity. Note this color feature is only available for target annotations created by Liftoff which have the sequence identity information in the GTF/GFF3.

##### Output File 2 
The second output of the synteny module is a tab separated file called 'gene_order' located in the output directory. If the -edit_distance option is used, the first line is the edit distance between the gene order in the reference and the gene order in the target which provides an estimate of how many genes in the target genome are in a different order compared to the reference. This file has the following fields:

#### Chromosome Ordering
the -r-sort and -t-sort options can be provided to establish the order in which the chromosomes are plotted for the reference genes and target genes respectively. These options must be provided together. The input format for these options is a simple .txt file with the order of the chromosomes like the example here. 
'''
chr1
chr2
chr3
chr4
''''
If -r-sort and t-sort are provided, only genes on the chromosomes in the .txt files will be plotted. Thus, if you do not wish to plot genes on all contigs/chromosomes, simply omit them from the .txt files. If these options are not provided, liftofftools will try to infer the matching order of the reference chromosomes and the target chromosomes by looking at where the genes from each reference chromosome mapped. 



|Col|Description                                                                                                                                  |
|--:|:--------------------------------------------------------------------------------------------------------------------------------------------|
|1  |Gene ID                                                                                                                                      |
|2  |Reference ordinal position ('NA' if gene is not present in the reference)                                                                    |
|3  |Target ordinal position ('NA' if gene is not present in the target)                                                                         |
|4  |Reference chromosome ('NA' if gene is not present in the reference)                                                                          |
|5  |Target chromosome ('NA' if gene is not present in the target)                                                                                |
|6  |Sequence identity ('NA' if gene is not present in both reference and target or if sequence identity information is not present in annotation)|



### Clusters Module
The clusters module can be run with the following command:

```
liftofftools clusters -r <reference.fa> -t <target.fa> -rg <reference.gff3> -tg <target.gff3>
```

This module clusters the genes into paralogous groups using [MMSeqs2](https://github.com/soedinglab/MMseqs2) to evaluate gene copy number gain and loss. MMSeqs2 clusters the amino acid sequences of the protein-coding genes, and the nucleotide sequences of noncoding genes. For each gene we select only the longest isoform to be included in the clustering. By default, genes are clustered together if they are at least 90% identical across at least 90% of both of their lengths. These parameters can be adjusted with the -mmseqs_params option. For example, ```-mmseqs_params "--min-seq-id 0.7 -c 0.8" ``` would cluster genes together that are at least 70% identical across 80% of their lengths. The -mmseqs_params option accepts any parameters accepted by the linclust subcommand of MMseqs2. This module also finds the closest paralog for every gene that failed to map to the target genome. 

#### Clusters Output

##### Output File 1
The first file output by the clusters module is a tab separated file called 'clusters' in the output directory (liftofftools_output by default). This file has the following fields for each cluster:

|Col|Description                                                                               |
|--:|:-----------------------------------------------------------------------------------------|
|1  |Cluster ID (prefixed by 'protein' for protein coding genes and 'DNA' for non-coding genes)|                                                         
|2  |The number of reference genes in the cluster                                              |
|3  |The IDs of the reference genes in the cluster (comma separated)                           |
|4  |The number of target genes in the cluster                                                 |
|5  |The IDs of the target genes in the cluster  (comma separated)                             |

The number of gene IDs in columns 3 and 5 will be the same as the integers in columns 2 and 4 respectively. Columns 2 and 4 are provided to easily identify clusters which have gained or lost copies in the target genome. 


##### Output File 2 
The second file output by the clusters module is a tab separated file called 'unmapped_closest_paralogs'. This file contains a line for every gene in the reference that failed to map to the target with the following fields:


|Col|Description                                                 |
|--:|:-----------------------------------------------------------|
|1  |ID of gene that failed to map to the target genome          |                                                         
|2  |ID of mapped paralog with the highest sequence identity     |
|3  |DNA sequence identity between unmapped gene and closest mapped paralog (0-1.0)     |
|4  |Protein sequence identity between unmapped gene and closest mapped paralog (0 -1.0)|



