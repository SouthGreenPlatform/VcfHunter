Tutorial for VcfHunter
======================

This tutorial go through all steps from variant calling from DNAseq
to marker generation for genetic mapping.
***VcfHunter tools***

Go to the VcfHunter folder (Scripts can be run from any folder but the
command lines in this tutorial assume you are in this folder and that
you have python 2.7 and 3 versions).

### Available data:

*data/reference/* is a folder containing a fasta file with the reference
sequence and a gff file containing the reference sequence annotation
*data/reads/* is a folder containing DNAseq paired reads simulated from
3 gene pools of diploid ancestral accessions (11, 10 and 9 accessions
from each gene pool respectively). *data/config/* is a folder containing
one configuration files needed for the DNAseq variant calling.

A specific folder as been created to test the tools (*TestTools*). Go to
this folder to run all command lines of this toolbox. But these tools
can also be launched from everywhere you have permission to write, you
will only have to change path to the configuration files.

For the following example, we will assume that you have a computer with
8 processors. If it is not the case or is you want to use less
processors than available you should/can change the number of processors
you allow the program to use.

In this tutorial we assume that you have performed the preliminary steps
of checking your reads and filtered if necessary. After that, it is time
to begin the variant calling.

A - Variant calling
-------------------

Because different mapping tools are needed to perfom the mapping
depending on the type of data (DNA), two different pipeline
shoold be use for read mapping and post mapping processing.

### Variant calling

a - Mapping of the DNA reads: In this step, reads are aligned against
the reference sequence using BWA mem algorithm

    python2 ../bin/process_reseq_1.0.py -c ../data/config/DNAseq.conf -t 8 -p DNAseq -s a

40 folders have been created containing each:

-   \*.sam.stat files with mapping statistics
-   \*\_merged.bam file containing aligned read merged in one bam file
-   \*\_merged.bai file (index of the bam file)
-   STAT folder containing a .html file summarising statistics on
    libraries.

b - Removing duplicates reads: In this step, read duplicates are removed

    python2 ../bin/process_reseq_1.0.py -c ../data/config/DNAseq.conf -t 8 -p DNAseq -s b

3 new files have been added in each folders:

-   \*\_rmdup.bam file containing non redundant reads resulting from PCR
    and optical duplicates
-   \*\_rmdup.bai file (index of the bam file)
-   \*\_duplicate file summarising read duplication statistics

c - Removing duplicates reads: In this step, read are realigned around
indels

    python2 ../bin/process_reseq_1.0.py -c ../data/config/DNAseq.conf -t 8 -p DNAseq -s c

2 new files have been added in each folders:

-   \*\_realigned.bam file reads realigned around indels
-   \*\_realigned.bai file (index of the bam file)

d - Base recalibration: In this step, reads base sequencing quality are
recalculated. **This step is recommended by GATK best practice but we do
not recommend to use it if you use our pipeline.**

    python2 ../bin/process_reseq_1.0.py -c ../data/config/DNAseq.conf -t 8 -p DNAseq -s d

5 new files have been added in each folders:

-   \*\_real\_recal.bam file reads bases recalibrated
-   \*\_real\_recal.bai file (index of the bam file)
-   \*\_recalibration\_plot.pdf file
-   \*\_selected.vcf file containing variant sites used for base
    recalibration (if no vcf have been provided)
-   \*\_selected.vcf.idx file (index of the vcf file)

e - Allele counting: In this step, reads are used to count for each
covered sites the number of reads supporting each bases
(A,T,G,C,N,\*=deletion).

    python2 ../bin/process_reseq_1.0.py -c ../data/config/DNAseq.conf -t 8 -p DNAseq -s e

Several new files have been added in each folders:

-   \*.gz files recording, for each chromosomes/sequences in the fasta
    provided as reference, the count of the number of reads supporting
    each bases (A,T,G,C,N,\*=deletion) at each covered site.

f - Vcf generation: generate the vcf for all accessions in the
configuration file

    python2 ../bin/process_reseq_1.0.py -c ../data/config/DNAseq.conf -t 8 -p DNAseq -s f

Several new files in the current directory:

-   prefix\_\*\_all\_allele\_count.vcf file for each
    chromosomes/sequences in the fasta provided as reference.

g - Vcf merging: generate a single vcf from all chromosome/sequences in
the fasta provided as reference. This step is not needed here has only
one sequence is passed but you can try the command line anyway ;-)

    python2 ../bin/process_reseq_1.0.py -c ../data/config/DNAseq.conf -t 8 -p DNAseq -s g

1 new file is created in the current directory:

-   prefix\_all\_allele\_count.vcf file

h - Statistics: Compute mapping statistics

    python2 ../bin/process_reseq_1.0.py -c ../data/config/DNAseq.conf -t 8 -p DNAseq -s h

Several new files in the current directory:

-   \*\_lib.stats file containing mapping statistics per libraries
-   \*\_acc.stats file containing mapping statistics per accessions

***All the steps can be lauched in one command line. First remove the
generated files:***

    rm -rf *

Run the pipeline again (without step d)

    python2 ../bin/process_reseq_1.0.py -c ../data/config/DNAseq.conf -t 8 -p DNAseq -s abcefgh


B - VCF prefiltering
--------------------

Once the variant calling is performed, it is now time to filter the
calling. This is a necessary step with this pipeline because all variant
sites are reported (one read, on one accession supporting a variant).
The phylosophie of ***process\_RNAseq*** and ***process\_reseq*** is to
report all variant sites and then decide to filter the vcf (what is
error and what is not) based on understandable and simple parametres
using ***VcfPreFilter*** tool. For the following step we will work on
the vcf file generated on both DNA and RNA seq data. VcfPreFilter can be
launched as followed:

    python2 ../bin/VcfPreFilter.1.0.py -v DNAseq_all_allele_count.vcf -m 10 -M 10000 -f 0.05 -c 3 -o DNAseq_prefiltered.vcf

The outpout is a vcf file (named as filled in -o option) in which the
variant line are filtered as followed:

-   1 - an accession of a variant line is considered only if its site
    coverage is comprised in [-m : -M] (in this example [10:10000]),
-   2 - for each considered accession, an allele is considered only if
    its frequency at in the accession if equal or greater than -f (in
    this example: 0.05) and its coverage is equal or greater than -c (in
    this example: 3),
-   3 - the line is kept only if their is at least one allele different
    from the reference kept in at least one accession at the end of the
    process,

The variant calling is recalculated **for each accessions** based on
selected alleles (even accessions which are not enough covered according
our parameters will have a genotype). An additional tag (GC) is added in
the FORMAT column of the VCF. This tag is calculated as followed :
-log10(best genotype probability) + log10(second best genotype
probability) and give an idea of the quality of the calling at the
datapoint.
