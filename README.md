# proseq_2.0
Preprocesses and Aligns Run-On Sequencing (PRO/GRO/ChRO-seq) data from Single-Read or Paired-End Illumina Sequencing

```
Preprocesses and aligns PRO-seq data.

Takes *.fastq.gz in the current working directory as input and writes
BAM and bigWig files as output to the user-assigned output-dir.

Requirements in current working directory:
cutadapt 1.8.3, prinseq-lite.pl 0.20.2, bwa, samtools, bedtools, bedGraphToBigWig

bash proseq.bsh [options]

options:

To get help:
-h, --help             Show this brief help menu.

Required options:
-SE, --SEQ=SE          Single-end sequencing
-PE, --SEQ=PE          Paired-end sequencing
-i, --bwa-index=PATH   Path to the BWA index of the target genome (i.e., bwa index).
-c, --chrom-info=PATH  Location of the chromInfo table.

I/O options:
-I, --fastq=PREFIX     Prefix for input files.
                       Paired-end files require identical prefix and end with _R1.fastq.gz and _R2.fastq.gz
                       eg: PREFIX_R1.fastq.gz, PREFIX_R2.fastq.gz
-T, --tmp=PATH         Path to a temporary storage directory.
-O, --output-dir=DIR   Specify a directory to store output in.

Required options for SE
-G, --SE_READ=RNA_5prime Single-end sequencing from 5' end of nascent RNA, like GRO-seq
-P, --SE_READ=RNA_3prime Single-end sequencing from 3' end of nascent RNA, like PRO-seq

Options for PE
--RNA5=R1_5prime    Specify the location of the 5' end of RNA [default: R1_5prime]
--RNA3=R2_5prime    Specify the location of the 3' end of RNA [default: R2_5prime]
                    Available options: R1_5prime: the 5' end of R1 reads
                                       R2_5prime: the 5' end of R2 reads
-5, --map5=TRUE     Report the 5' end of RNA [default on, --map5=TRUE]
-3, --map5=FALSE    Report the 3' end of RNA, only available for PE [default off, --map5=TRUE]
-s, --opposite-strand=TRUE
                    Enable this option if the RNA are at the different strand
                    as the reads set at RNA5 [default: disable]

Optional operations:
--ADAPT1=GATCGTCGGACTGTAGAACTCTGAACG
                    3' adapter to be removed from the 3' end of R2.
--ADAPT2=AGATCGGAAGAGCACACGTCTGAACTC
                    3' adapter to be removed from the 3' end of R1.

--UMI1=0            The length of UMI barcode on the 5' of R1 read. [default: 0]
--UMI2=0            The length of UMI barcode on the 5' of R2 read. [default: 0]
--add_barcode1=0    The length of additional barcode that will be trimmed
                    on the 5' of R1 read. [default: 0]
--add_barcode2=0    The length of additional barcode that will be trimmed
                    on the 5' of R2 read. [default: 0]
--thread=1          Number of threads can be used [default: 1]
```
