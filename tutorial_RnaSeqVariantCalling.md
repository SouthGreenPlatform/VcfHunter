Tutorial for VcfHunter
======================

This tutorial go through all steps of RNA seq variant calling from read
mapping to vcf filtration.

***VcfHunter tools***

Go to the VcfHunter folder (Scripts can be run from any folder but the
command lines in this tutorial assume you are in this folder and that
you have python 3 version).

### Available data:

- *data/reference/* is a folder containing a fasta file with the reference
sequence and a gff file containing the reference sequence annotation
- *data/reads/* is a folder containing RNAseq paired reads simulated from
10 hybrids obtained from admixture events between 3 gene pools of diploid
ancestral accessions. 
- *data/config/* is a folder containing two configuration files needed for
the RNAseq and DNAseq variant calling.

A specific folder as been created to test the tools (*TestTools*). **Go to
this folder to run all command lines of this toolbox**. But these tools
can also be launched from everywhere you have permission to write, you
will only have to change path to the configuration files.

For the following example, we will assume that you have a computer with
8 processors. If it is not the case or is you want to use less
processors than available you should/can change the number of processors
you allow the program to use.

In this tutorial we assume that you have performed the preliminary steps
of checking your reads and filtered if necessary. After this small introduction,
it is time to begin the variant calling.

# A - Variant calling
----------------------

## a - Reference indexation
----------------------------

Create an index for STAR RNAseq mapping

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s a

1 folder is created in the current directory (name =
prefix+_ref_star\_1). This folder contains index constructed and used by
STAR to align RNAseq reads

## b - Identification of splicing sites
---------------------------------------

During this step, all reads from
all libraries are merged and aligned against the reference sequence
which allowed STAR to identify splicing sites which will be used to
build a new index (step c) taking in account the splicing site.

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s b

A file (prefix+"_JUNC_ESTIMATION_SJ.out.tab") containing the splicing
sites is generated

## c - Reference indexation
---------------------------

This time the reference is indexed but with the splicing site information identified with previous step

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s c

1 folder is created in the current directory (name =
prefix+_ref_star_2). This folder contains index constructed and used by
STAR to align RNAseq reads

## d - Read mapping
-------------------

Reads are aligned against the reference sequence using STAR software. Read are indexed one library after the other.

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s d

10 folders have been created containing each:

-   \*Aligned.out.sam files containing aligned reads (1 per libraries)
-   \*Log.final.out file containing mapping statistics of STAR (1 per
    libraries)
-   \*.SJ.out.tab file containing splicing sites (1 per libraries)

1 additional folder was created and named as passed to -p option
containing a file named prefix + "\_mapping.tab" which summarise mapping
statistics for all accessions.

## e - Read merging and ordering
--------------------------------

Reads from different libraries but the same accessions are merged and sorted on coordinates.

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s e

Several new files have been added in each folders:

-   \*_merged.bam files containing all aligned reads from all libraries
    of one accession.
-   \*_merged.bai files (index of the bam file)

## f - Removing duplicates reads
--------------------------------

In this step, read duplicates (resulting from PCR before sequencing) are removed.

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s f

2 new files have been added in each folders:

-   \*_rmdup.bam file containing non redundant reads resulting from PCR
    and optical duplicates
-   \*_duplicate file summarising read duplication statistics

An additional file was created in the folder named as passed to -p
option. This file named prefix + "*rmdup*stat.tab" contained duplicates
statistics for all accessions.

## g - Bam reordering
---------------------

In this step reads are reordered in bam files. This is necessary for subsequent steps.

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s g

2 new files have been added in each folders:

-   \*_reorder.bam file containing reordered reads based on the order
    sequences are provided in the multifasta
-   \*_reorder.bai file (index of the bam file)

## h - Reads splitting
----------------------

In this step, reads are split on splicing sites

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s h

2 new files have been added in each folders:

-   \*_trim.bam file containing splitted reads
-   \*_trim.bai file (index of the bam file)

## i - Indel realignment
------------------------

In this step, reads realigned around indels. Countrary to **process\_reseq** a vcf on known
indel could not be provided. It will be implemented in future...

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s i

2 new files have been added in each folders:

-   \*_realigned.bam file containing realigned reads
-   \*_realigned.bai file (index of the bam file)

## j - Allele counting
----------------------

In this step, reads are used to count for each covered sites the number of reads supporting
each bases (A,T,G,C,N,\*=deletion).

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s j

Several new files have been added in each folders:

-   \*.gz files recording, for each chromosomes/sequences in the fasta
    provided as reference, the count of the number of reads supporting
    each bases (A,T,G,C,N,\*=deletion) at each covered site.

## k - Vcf generation
----------------------

Generate the vcf for all accessions in the configuration file

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s k

Several new files in the current directory:

-   prefix_\*_all_allele_count.vcf file for each
    chromosomes/sequences in the fasta provided as reference.

## l - Vcf merging
------------------

Generate a single vcf from all chromosome/sequences in
the fasta provided as reference. This step is not needed here has only
one sequence is passed but you can try the command line anyway ;-)

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s l

1 new file is created in the current directory:

-   prefix_all_allele_count.vcf file

## m - Exon coverage
--------------------

Calculate for each gene filled in the gff file passed
in the configuration file and for each accessions, the proportion of the
exon covered by the library.

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s m

1 new file is created in the directory named RNAseq passed with -p options:

-   RNAseq_exon_coverage.cov file which contained information on exon coverage

***All the steps can be launched in one command line. First remove the
generated files:***

    rm -rf sample71 sample72 sample73 sample74 sample75 sample76 sample77 sample78 sample79 sample80

Run the pipeline again

    python3 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s abcdefghijklm


# B - VCF prefiltering
----------------------

Once the variant calling is performed, it is now time to filter the
calling. This is a necessary step with this pipeline because all variant
sites are reported (one read, on one accession supporting a variant).
The phylosophie of ***process\_RNAseq*** is to
report all variant sites and then decide to filter the vcf (what is
error and what is not) based on understandable and simple parameters
using ***VcfPreFilter*** tool which also performs a better calling
with these selected sites. For the following step we will work on
the vcf file generated on **Variant calling** section. VcfPreFilter can be
launched as followed:

    python3 ../bin/VcfPreFilter.1.0.py -v RNAseq_all_allele_count.vcf -m 10 -M 10000 -f 0.05 -c 3 -o RNAseq_prefiltered.vcf -d y

The output is a vcf file (named as filled in -o option) in which the
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
(best genotype probability)/(second best genotype probability) and give an
idea of the quality of the calling at the data point.
Two "types" of genotypes can be called: if the ***-d y*** option is passed (default),
only two alleles are authorized in a genotype (*ie* a triploid accession cannot
be 0/1/2). This gain computation time but genotype can be erroneous if divergent
time between haplotypes is great. If ***-d n*** option is passed, all possible
combinations of each alleles are tested and not only bi-allelic combination. In
this case the computation time can be a lot more important if the ploidy level is
high.

# C - VCF Filtering
-------------------

Sometime we want to filter the VCF to keep only calling for some accession,
only biallelic sites, convert to unknown genotype data-point which are not
sufficiently covered, remove variant sites with too much missing data, ...
The purpose of the following script is to perform those filter. As there are
several filtering options possible, we will show an example of filtering that
can be applied but you can have a look at filtering options available in the
README file or directly in with the program with the following command line:

    python3 ../bin/vcfFilter.1.0.py -h

In this tutorial we will use vcfFilter.1.0.py program to apply some filter to
all accessions of the vcf. To do so, we first need to create a file in which
we list accessions to perform filter on. This can be done with the following
command line:

    head -n 1000 RNAseq_prefiltered.vcf | grep "#CHROM" | sed 's/\t/\n/g' | tail -n +10 > all_names.tab

When performing this filter, we only want to keep bi-allelique sites. In other
word we want to remove mono-allelic, tri-allelic, tetra-allelic sites. However
as there are two other type of variant state possible with our pipeline: unknown
base (N) or deletion relative to reference (*) penta and hexa allelic state can
also be possible we also need to remove these state. This can be done with the
following option: *--RmAlAlt 1:3:4:5:6.*
We also want to convert to missing data, all data-points which are too much covered
(probably resulting from repeat sequences and thus being multiloci) and those who
are not enough covered (in which variant calling may be approximative). This can be
done with the following options: *--MinCov 10 --MaxCov 300*.
In addition, we want that each alleles of a genotype is covered by at least 3 reads.
If not, this genotype is converted to missing data. This can be done with the
following option: *--MinAl 3*.
At last we do not want more than 1 missing data per variant line. This can be done
with the following option: *--nMiss 1*.
We want the new vcf to be generated in a file prefixed DNAseq_Filtered. This can be
done with the following option: *--prefix DNAseq_Filtered*.
At last we want to compress the final vcf to gain disc space. This can be done with
the following command line: -g y.
Now we run the analysis:
 
    python3 ../bin/vcfFilter.1.0.py --vcf RNAseq_prefiltered.vcf --names all_names.tab --MinCov 10 --MaxCov 300 --MinAl 3 --nMiss 1 --RmAlAlt 1:3:4:5:6 --prefix RNAseq_Filtered -g y

The output is a vcf file (named ***RNAseq_Filtered_filt.vcf.gz***) in which the
variant line have been filtered.


