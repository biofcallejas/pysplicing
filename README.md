# TACOS
## Trichomonad Assembler for COmplex Splicing

This semi-automated pipeline was designed to filter the Trichomonas vagialis transcriptome for chimeric alignments. 
The Trichomonas genome, full of complex repeats, can produce hundreds of thousands of chimeric mappings from RNAseq data, making it difficult to produce accurate transcript quantifications and intron detection.
The main part of this pipeline includes a filtering process for those chimeric mappings by leveraging the conserved intron motifs, the genome sequence, and the transcript evidence.

***The introns in these species have conserved intron motifs at 5' and 3' ends beyond the universal AG/GT dinucleotides.***

Although it was designed for Trichomonas vaginalis, the main python script can be used for other species and custom intron motifs, up to seven nucleotides at the 5' and up to twelve at the 3' end. 

## Installation

Before running TACOS, please make sure you meet the following requirements:

TACOS has been tested on Python 3.8.12 and the following libraries that might require manual installation:

[pysam 0.20.0](https://github.com/pysam-developers/pysam)
[six 1.16.0](https://pypi.org/project/six/)

```
STAR
StringTie
Samtools
```

## TACOS pipeline is divided into three main steps:

```
1.- RNAseq mapping to the reference genome.
2.- Processing of the SJs.
3.- Transcript assembly
```

![tacos_summary](https://user-images.githubusercontent.com/45425927/219090905-d3c7e9dd-7d35-4b2a-929c-e44cd968ffc0.jpg)

The complete pseudo-code can be summaryzed as follows:

![pseudo_code](https://user-images.githubusercontent.com/45425927/219105430-0dc3a9c5-fd9c-44e5-916b-74207004f82d.jpg)


## Running TACOS

###### 1. Mapping 
Mapping parameters. TACOS has been tested with STAR 2.7.6a.
Custom parameters for mapping are not needed (nor reccomended), however the use of the following parameters is mandatory:

```
--outSAMtype BAM SortedByCoordinate 
--outSAMstrandField intronMotif 
--outSAMattributes All
```

BAM file must be indexed, and the index must be present in the same path.
It can be obtained with samtools by running:

samtools index file.bam

###### 2. Running the main script

```
python tacos_v2.py -h 

usage: tacos_v2.py [-h] -f INPUT -sj SJ_STAR -b BAM_F -o OUTPUT -5m MOTIF5P -3m MOTIF3P

Trichomonad assembler for complex splicing
---------------------
Tested on python 3.8.12

optional arguments:
  -h, --help            show this help message and exit

Mandatory arguments:
  -f INPUT, --input, Input fasta file (Genome reference)
                        
  -sj SJ_STAR, --input_sj, SJ.out.tab from STAR mapping
                        
  -b BAM_F, --input_bam, BAM file from STAR mapping
                        
  -o OUTPUT, --output, Prefix: Prefix for output files
                        
  -5m MOTIF5P, --5p-motif, String: Splicing motif at 5p (upper case nucleotides only, no spaces)
                        
  -3m MOTIF3P, --3p-motif, String: Splicing motif at 3p (upper case nucleotides only, no spaces)
```

###### 3. Transcript assembly 

*In process ...*

## Running TACOS with test files

*In process ...*

