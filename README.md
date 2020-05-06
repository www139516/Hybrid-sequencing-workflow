# Hybrid-sequencing-workflow
These scripts are used for processing the hybrid sequencing data from maize, utilizing the RNA-seq data to correct and checking the whole chain of splicing junction of reads identified in PacBio-sequencing.

The workflow are based on the STAIR pipeline. We have edited and optimized it to ensure it fits for the IsoSeq3 and sequencing data from maize.

The workflow is depended on the following softwares/packages:
1. IsoSeq3
2. python 3.0 or higher
3. pandas
4. click
5. gmap
6. proovread
### Design
* /main.py
* /funcs/
* --funcs.py some common functions used by different classes
* /modules/
* --sub2flnc/: generates FLNC based on the subreads
* --sj_validation/: validate the whole chain of splicing junctions in each isoform
* --sort_sam/: sort the .sam file
* --improve_long_read_seq/: using RNA-seq data to improve the quality of long-read sequencing, using "proovread"
* --mapping/: map the long reads to the genome to generate .sam file, using "gmap"

