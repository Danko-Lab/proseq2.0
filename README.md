# Proseq2.0
Preprocesses and Aligns Run-On Sequencing (PRO/GRO/ChRO-seq) data from Single-Read or Paired-End Illumina Sequencing. It is a new version of ProseqMapper (https://github.com/Danko-Lab/utils/tree/master/proseq).

Currently we provide two commands: proseq mapper and bigWig merge.

# NOVOGENE MODIFICATION:

**Novogene paired end files are named PREFIX_1.fq.gz and PREFIX_2.fq.gz. This fork changes the references in proseq2.0.bsh to use this standard.**

## Overview
Our proseq2.0 pipeline will take single-end or paired-end sequencing reads in fq.gz format as input. The pipeline will automate three routine pre-processing and alignment options, including
+ pre-processing reads: remove the adapter sequence and quality trim the reads (cutadapt), deduplicate the reads if UMI barcodes are used (prinseq-lite.pl)
+ mapping reads to a reference genome (BWA)
+ converting BAM files into bedGraph and BigWig formats (kentsource). When converting to bedGraph and BigWig, the pipeline only report the 5’ end position of the reads after UMI/adapter removal. For pair-end sequencing, user can choose to report the 5’ end of R1 or R2 reads. The output BigWig files ending with _minus.bw or _plus.bw are raw read counts without normalization. The RPM normalized outputs end with a suffix of .rpm.bw.

To run our pipeline users must first download the pipeline files and install dependencies indicated in this README.md. In addition, user need to provide a path to a BWA index file and the path to the chromInfo file for the genome of choice. After running this pipeline, users should have processed data files in the specified output directory.

Please note that proseq2.0 does __NOT__ include steps to ensure chrM reads are really from mitochondrial RNA polymerase. All reads mapped to rRNA or chrM will be removed from bigWig files.

## Citation
Chu, T., Wang, Z., Chou, S. P., & Danko, C. G. (2018). Discovering Transcriptional Regulatory Elements From Run‐On and Sequencing Data Using the Web‐Based dREG Gateway. Current protocols in bioinformatics, e70.

## Dependencies

The pipelines depend on several common bioinformatics tools:
- [ ] cutadapt (https://cutadapt.readthedocs.io/en/stable/installation.html)
- [ ] fastx_trimmer (http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)
- [ ] seqtk (https://github.com/lh3/seqtk)
- [ ] prinseq-lite.pl (https://sourceforge.net/projects/prinseq/files/standalone/)
- [ ] bwa (https://sourceforge.net/projects/bio-bwa/files/)
- [ ] samtools version: 1.9 (http://www.htslib.org/download/)
- [ ] bedtools v2.28.0 (http://bedtools.readthedocs.org/en/latest/)
- [ ] bedGraphToBigWig (from the Kent source utilities http://hgdownload.cse.ucsc.edu/admin/exe/)

Please make sure you can call the bioinformatics tools from your current working directory.

## Usage
```
Preprocesses and aligns PRO-seq data.

Takes PREFIX.fq.gz (SE),  PREFIX_1.fq.gz, PREFIX_2.fq.gz (PE)
or *.fq.gz in the current working directory as input and writes
BAM and bigWig files as output to the user-assigned output-dir.
The output bigWig files ending with _minus.bw or _plus.bw are raw read counts without normalization.
The RPM normalized outputs end with a suffix of .rpm.bw.


Requirements in current working directory:
cutadapt 1.8.3, fastx_trimmer, seqtk, prinseq-lite.pl 0.20.2, bwa, samtools, bedtools, and bedGraphToBigWig.

bash proseq2.0.bsh [options]

options:

To get help:
-h, --help             Show this brief help menu.

Required options:
-SE, --SEQ=SE          Single-end sequencing.
-PE, --SEQ=PE          Paired-end sequencing.
-i, --bwa-index=PATH   Path to the BWA index of the target genome
                       (i.e., bwa index).
-c, --chrom-info=PATH  Location of the chromInfo table.

I/O options:
-I, --fastq=PREFIX     Prefix for input files.
                       Paired-end files require identical prefix
                       and end with _1.fq.gz and _2.fq.gz
                       eg: PREFIX_1.fq.gz, PREFIX_2.fq.gz.
-T, --tmp=PATH         Path to a temporary storage directory.
-O, --output-dir=DIR   Specify a directory to store output in.

Required options for SE
-G, --SE_READ=RNA_5prime Single-end sequencing from 5' end of
                         nascent RNA, like GRO-seq.
-P, --SE_READ=RNA_3prime Single-end sequencing from 3' end of
                         nascent RNA, like PRO-seq.

Options for PE
--RNA5=R1_5prime    Specify the location of the 5' end of RNA
                    [default: R1_5prime].
--RNA3=R2_5prime    Specify the location of the 3' end of RNA
                    [default: R2_5prime].
                    Available options: R1_5prime: the 5' end of R1 reads
                                       R2_5prime: the 5' end of R2 reads
-5, --map5=TRUE     Report the 5' end of RNA [default on, --map5=TRUE].
-3, --map5=FALSE    Report the 3' end of RNA,
                    only available for PE [default off, --map5=TRUE].
-s, --opposite-strand=TRUE
                    Enable this option if the RNA are at the different strand
                    as the reads set at RNA5 [default: disable].

Optional operations:
--ADAPT_SE=TGGAATTCTCGGGTGCCAAGG
                    3' adapter to be removed from the 3' end of SE reads.
                   [default:TGGAATTCTCGGGTGCCAAGG]
--ADAPT1=GATCGTCGGACTGTAGAACTCTGAACG
                    3' adapter to be removed from the 3' end of R2.
                   [default:GATCGTCGGACTGTAGAACTCTGAACG]
--ADAPT2=AGATCGGAAGAGCACACGTCTGAACTC
                    3' adapter to be removed from the 3' end of R1.
                   [default:AGATCGGAAGAGCACACGTCTGAACTC]

--UMI1=0           The length of UMI barcode on the 5' of R1 read.
                   [default: 0]
--UMI2=0           The length of UMI barcode on the 5' of R2 read.
                   [default: 0]
When UMI1 or UMI2 are set > 0, the pipeline will perform PCR deduplicate.

--Force_deduplicate=FALSE
                   When --Force_deduplicate=TRUE, it will force the pipeline to
                   perform PCR deduplicate even there is no UMI barcode
                   (i.e. UMI1=0 and UMI2=0). [default: FALSE]
--ADD_B1=0         The length of additional barcode that will be trimmed
                   on the 5' of R1 read. [default: 0]
--ADD_B2=0         The length of additional barcode that will be trimmed
                   on the 5' of R2 read. [default: 0]
--thread=1         Number of threads can be used [default: 1]

-4DREG             Using the pre-defined parameters to get the most reads
                   for dREG package. Please use this flag to make the bigWig
                   files compatible with dREG algorithm. Only available for
                   Single-end sequencing.[default: off]
-aln               Use BWA-backtrack [default: SE uses BWA-backtrack (aln), PE uses BWA-MEM (mem)]
-mem               Use BWA-MEM [default: SE uses BWA-backtrack (aln), PE uses BWA-MEM (mem)]
--MAP_LENGTH       Set a data-set wide length cutoff for mapping (e.g. --MAP_LENGTH=36)
                   or map the whole reads after QC (off, --MAP_LENGTH=0). [default: off]

```
<img src="images/lib.png">


## Examples
The pipeline requires two parameters for genome information, including BWA index (--bwa-index) and chrom info (--chrom-info).

__BWA index__ should be generated using the __bwa index__ command according to BWA manual at http://bio-bwa.sourceforge.net/bwa.shtml . Please note that the program only take in the prefix when you assign the index, no ".bwt" in the end. See the BWA manual for more details.

__Chrom info__ is a tab-delimited file with two columns. The first column is the chromosome name and the second is the size of the chromosome. Please see see /input_file_exmaples/mm10.chromInfo for example

```
bwa index Genome.fa.gz                # Don't do this if you already have bwa index
export bwaIndex=Genome.fa.gz          # Or path to your bwa index
export chromInfo=PathToChromInfo
```

### Example 1

PREFIX.fq.gz were made according to GRO-seq protocol as in  https://www.ncbi.nlm.nih.gov/pubmed/19056941
Give UMI1=6, the pipeline will remove PCR duplicates and trim the 6bp UMI barcode.
```
bash proseq2.0.bsh -i $bwaIndex -c $chromInfo -SE -G -T myOutput1 -O myOutput1 --UMI1=6 -I PREFIX
```
### Example 2

PREFIX.fq.gz were made according to PRO-seq protocol as in  https://www.ncbi.nlm.nih.gov/pubmed/23430654
UMI1 was set to 0 by default. The pipeline will NOT remove PCR duplicates.
```
bash proseq2.0.bsh -i $bwaIndex -c $chromInfo -SE -P -T myOutput2 -O myOutput2 -I PREFIX
```
### Example 3

__PREFIX_1.fq.gz__ and __PREFIX_2.fq.gz__ were Paired-End sequenced as in chromatin run-on and sequencing (ChRO-seq) in https://www.biorxiv.org/content/early/2017/09/07/185991
* Please note that Paired-end files require identical PREFIX and end with _1.fq.gz and _2.fq.gz.

  Assign the file use __-I PREFIX__. No _1.fq.gz, _2.fq.gz, nor *fq.gz is in the end.
* There is a 6N UMI barcode on R1. Pipeline will perform PCR deduplicat.
```
bash proseq2.0.bsh -i $bwaIndex -c $chromInfo -PE --RNA3=R1_5prime -T myOutput3 -O myOutput3 -I PREFIX --UMI1=6 --ADAPT1=GATCGTCGGACTGTAGAACTCTGAAC --ADAPT2=TGGAATTCTCGGGTGCCAAGG
```
### Example 4
Same as in Example 3 but without UMI barcode.
* UMI1 and UMI2 were set to 0 by default. The pipeline will NOT remove PCR duplicates.
```
bash proseq2.0.bsh -i $bwaIndex -c $chromInfo -PE --RNA3=R1_5prime -T myOutput4 -O myOutput4 -I PREFIX --ADAPT1=GATCGTCGGACTGTAGAACTCTGAAC --ADAPT2=TGGAATTCTCGGGTGCCAAGG
```
### Example 5
If you made ChRO-seq library with High-throughput Adapters from Danko lab __AND all the samples have distinct i7 barcode.__
* Please note that proseq2.0.bsh does __NOT__ work if there your samples share i7 barcode.
```
bash proseq2.0.bsh -i $bwaIndex -c $chromInfo -PE --UMI1=4 --UMI2=4 --ADD_B1=6 -T myOutput5 -O myOutput5 -I PREFIX
```

## Useful references

* GRO-seq: http://www.sciencemag.org/content/322/5909/1845.long
* PRO-seq: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3974810/
* dREG: http://www.nature.com/nmeth/journal/v12/n5/full/nmeth.3329.html

## Notes for **CBSUdanko** users:

1. Setup your environment to use the bioinformatics tools (e.g. prinseq-lite.pl,bedGraphToBigWig,samtools...)
```
export PATH=$PATH:/programs/prinseq-lite-0.20.2:/programs:/home/zw355/lib/bin:/home/zw355/lib/ucsc
```

2. Find the BWA index and chromosome table in the server:
```
export human_genome=/local/storage/data/short_read_index/hg19/bwa.rRNA-0.7.5a-r405/hg19.rRNA
export human_chinfo=/local/storage/data/hg19/hg19.chromInfo

export mouse_genome=/local/storage/data/short_read_index/mm10/bwa.rRNA-0.7.8-r455/mm10.rRNA.fa.gz
export mouse_chinfo=/local/storage/data/mm10/mm10.chromInfo

export dog_genome=/local/storage/data/short_read_index/canFam3/bwa.rRNA-0.7.8-r455/canFam3.rRNA.fa
export dog_chinfo=/local/storage/data/canFam3/canFam3.chromInfo
```

3. Using --UMI1=6 to replace -b6 if you have used it in the old version (proseqMapper.bsh).

## Notes for **dREG** users:

In order to make the most compatible with dREG algorithm, please use **-4DREG** flag when you process the PRO-seq and GRO-seq reads. The dREG package needs enriched reads to
detect the transcriptional peaks, we use the "bwa aln" to do mappping and set lower filtering score (0) to get the most reads in this pipeline. Only available for Single-end sequencing.

Here is an examples to generate the bigWig for dREG.

```
bash proseq2.0.bsh -SE -G -4DREG -i $dog_genome -c $dog_chinfo -I ./example1_R1  -T ./tmpdir -O ./outputdir
```

