# TACOS
## Trichomonad Assembler for COmplex Splicing

Before running TACOS (Trichomonad Assembler for COmplex Splicing), please make sure you meet the following requirements:

TACOS has been tested on python 3.8.12, and needs the following libraries (some might require manual installation):

```
import tabulate
import plotext
import pysam
import six 
```

Mapping parameters. TACOS has been tested with STAR 2.7.6a.
Custom mapping parameters are not needed (nor reccomended), however the following parameters are mandatory:

```
--outSAMtype BAM SortedByCoordinate 
--outSAMstrandField intronMotif 
--outSAMattributes All
```

BAM file must be indexed, and the index must be present in the same path.
It cab be obtained with samtools by running:

samtools index file.bam