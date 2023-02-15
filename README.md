# TACOS
## Trichomonad Assembler for COmplex Splicing

Before running TACOS (Trichomonad Assembler for COmplex Splicing), please make sure you meet the following requirements:

TACOS has been tested on python 3.8.12, and needs the following libraries/spftware (some might require manual installation):

```
import tabulate
import plotext
import pysam
import six

STAR
StringTie
Samtools
```

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

TACOS pipeline is divided into three main steps:
```
1.- RNA seq mapping to the reference genome.
2.- Processing of the SJs.
3.- Transcript assembly
```

![tacos_summary](https://user-images.githubusercontent.com/45425927/219090905-d3c7e9dd-7d35-4b2a-929c-e44cd968ffc0.jpg)


