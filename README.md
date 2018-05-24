Purpose of VcfHunter
====================

VcfExplorer regroups several programs which principale aims are to map
RNAseq data onto reference genome sequence, perform variant calling,
manipulate vcf files, perform chromosome painting of accessions based
on the contribution of ancestral groups and select marker for genetic
map analysis.
<br><br><br>

Installation
------------

All proposed tools described here are written in python and work on
linux system To install the tools:

1.  open the loca\_programs.conf file
2.  set the path to each programs required
<br><br><br>

Dependencies
------------

1.  PicarTools, https://broadinstitute.github.io/picard/
2.  GATK, https://software.broadinstitute.org/gatk/
3.  Samtools, https://github.com/samtools/samtools
4.  Bamtools, https://github.com/pezmaster31/bamtools
5.  bam-readcount, https://github.com/genome/bam-readcount
6.  gnuplot, http://www.gnuplot.info/

Python2, Python3, Java and Biopython are also required.
<br><br><br>

Description
-----------

The package provided comprised X programs listed here:

-   process\_reseq.1.0.py (python2)
-   VcfPreFilter.1.0.py (python2)
-   vcfFilter.1.0.py (python3)
-   vcf2pop.1.0.py (python3)

All X programs run using the following command: python program-name <--options-name value>
<br><br><br>

Programs
--------
<br><br>

### process\_reseq.1.0.py (python2)

This program takes a reference DNA sequence multifasta file and several
fastq files and returns a bam file for each accessions and a final VCF
file containing alleles count at each variant site having at least one
variant allele supported by at least one read.

![](http://banana-genome-http.cirad.fr/image/Process_ReSeq_Fig1.png)

#### Options:

    --conf: A configuration file containing path to references sequence (multifasta file) and RNAseq reads (fastq files).
    --thread: Max number of accessions treated at the same time. Do not exceed the number of processors available! [default: 1] 
    --queue: If you are using SGE sheduler: the queue name for the job to perform parallelisation. If not do not fill.
    --prefix: Prefix for vcf file and statistics folders.
	--chrom: Chromosomes to work with (only for step f). If "all" : all chromosomes will be used for calling. Otherwise : a list of chromosome names separated by ":" [Default:all]
    --steps: A string containing steps to perform:
        a: Aligning libraries
        b: Removing duplicates
        c: Indel realignment
        d: Bases recalibration
        e: Allele counting
        f: Genotype calling
        g: Merging genotype calling
        h: Mapping statistics calculation

#### Configuration file description:

The configuration file should contain 4 sections and must be formated as
followed:

    [Libraries]
    lib1 = genome_name path_to_mate1 path_to_mate2 ploidy
    lib2 = genome_name path_to_single ploidy
    ...
    [Reference]
    genome = path_to_the_reference_sequence

#### Output:

**Warning:** This program need to create a .dict and a .fai file un the
folder were the reference sequence is stored if they do not already
exist. Make sure that you have right to write in this folder!

Outputs are dependent of the steps you are running and each steps use
the output of the preceding one.

-   **step a:** generates a folder for each accession, (names filled in
    column 3 "genome\_name") filled in the configuration file, which
    contained the bam (\*\_merged.bam) and bai (\*\_merged.bai) files of
    aligned reads, a .stat file generated for each libraries with
    **samtools stat** program and a STAT folder containing a html files
    summarising mapping statistics,
-   **step b:** generates a bam (\*\_rmdup.bam) and bai (\*\_rmdup.bai)
    files for each accessions with duplicated reads removed. In addition
    in each folder duplicate statistics wer recorder in a file named
    \*\_duplicate,
-   **step c:** generates a bam (\*\_realigned.bam) and bai
    (\*\_realigned.bai) files realigned around indel for each
    accessions,
-   **step d:** generates a ,
-   **step e:** generates for each accessions and each chromosomes a
    file (Accession + "*allele*count\_"+ chromosome + ".gz") counting
    variant at each covered bases,
-   **step f:** generates a generates several vcf (one for each
    reference sequences) (--prefix + reference sequence +
    "*allele*count.vcf") file counting for each variant sites (at least
    on reads supporting a variant) and each accession the number of
    reads supporting each allele. For each accession in the vcf a
    genotype (GT tag) was called based on a binomial test, allelic depth
    was counted (AD tag) and total depth was rapported (DP tag),
-   **step g:** generates a generates a uniq vcf file resulting in the
    concatenation of all vcf files generated at step g,
-   **step h:** generates two files (\*\_acc.stats and \*\_lib.stats)
    collecting mapping statistics on each libraries and accessions
    respectively
<br><br>

### VcfPreFilter.1.0.py (python2)

This script filter VCF file generated by Process\_RNAseq.1.0.py by
removing homozygous sites for all accessions based on accession minimal
and maximal coverage, accession minimal allele coverage and frequency
parameters passed. The idea of this tool is to filter out variant lines
resulting from sequencing errors. Filter are applicated as followed:

-   only data points covered by at least "minimal coverage" reads were
    considered,
-   only data points covered by at most "maximal coverage" reads were
    considered,
-   only variant alleles supported by at least "minimal allele coverage"
    reads and having a frequency equal or greater to "minimal frequency"
    in one accession were kept as variant,
-   sites showing at least one variant allele were kept for variant
    calling. For each accession and at each variant site, a genotype was
    called based on the maximum likelihood of all possible genotype
    calculated based on a binomial distribution assuming a sequencing
    error rate of 0.005. The variant calling file was formated to vcf
    format.

#### Options:

    --vcf: The VCF file
    --MinCov: Minimal read coverage for site. [Default: 10]
    --MaxCov: Maximal read coverage for site. [Default: 1000]
    --minFreq: Minimal allele frequency in an accession to keep the allele for calling in the row
    --MinAlCov: Minimal read number of minor allele to call variant heterozygous (between 1 and infinity). [Default: 3]
    --out: Prefix for output files. [Default: Pop]
<br><br>

### vcfFilter.1.0.py (python3)

This program filter a vcf file based on several criterias. This
is a improved version of FILTER option of vcf2struct.1.0.py

*Options:*

    --vcf: A standard vcf file
    --names: A one column file containing accessions to treat.
    --outgroup: (optional) A one column file containing accession names that will not be used for filtering but will remain in the
	 output file.
    --RmType: (optional) Variant status to filter out (several values can be passed in this case they should be separated by ","). 
        Possible values: 
            *Values which can be found in the FILTER column: PASS, DP_FILTER, QD_FILTER, SnpCluster, 
            *Other values: INDELS, SNP, AUTAPO (accession specific variant site).
    --RmAlAlt: (optional) Number of alleles at the site in filtered accessions to remove the variant site (several values can be
	 passed and should be sepatated by ":"). Values: 1,2,3,...,n
    --MinCov: Minimal coverage by accession to keep genotype calling (integer). If the value is lower, genotype will be converted to
	 unknown for the concerned accession. [Default: 10]
    --MaxCov: Maximal coverage by accession to keep genotype calling (integer). If the value is greater, genotype will be converted to
	 unknown for the concerned accession. [Default: 1000]
    --MinFreq: Minimal allele frequency to keep genotype calling (float). If the value is lower, genotype will be converted to unknown
	 for the concerned accession. [Default: 0.05]
    --MinAl: Minimal allele coverage by accession to keep genotype calling (integer). If the value is lower for at least one allele,
	 genotype will be converted to unknown for the concerned accession. [Default: 3]
    --nMiss: Maximal number of missing genotype in a line to keep the line (integer). [Default: 0]
    --prefix: The prefix for output files. [Default: WorkOnVcf]

*Output:*\
 **\*\_filt.vcf:** a filtered vcf file based on passed options.
<br><br>
 
### vcf2pop.1.0.py (python3)

This program will select markers for genetical mapping
analysis from a vcf file based on several criterias. It will
outpout coded markers for two genetical mapping software
(onemap and joinmap).

*Options:*

    --vcf: The vcf file
	--MinCov: Minimal read coverage for a marker in an accession (interger). If a lower value is
	 found data point is converted to missing. [Default: 10]
	--MaxCov: Maximal read coverage for a marker in an accession (interger). If a greater value is
	 found data point is converted to missing. [Default: 1000]
	--WinFreq: Window for minority allele coverage frequency to be insufficient to call a
	 heterozygous but to high to call an homozygous (example: "0.05:0.1"). With the example if
	 minority allele is in ]0.05:0.1] calling will become missing for this data point.
	--MinAlCov: Minimal read number of minor allele to call variant heterozygous (between 1 and
	 infinity). [Default: 1]
	--miss: Maximal missing data proportion in the progeny (Excluding parents) (between 0 and 1).
	 greater missing proportion will result in removing the marker. [Default: 0.2]
	--pValue: P-value threshold to keep marker (between 0 and 1). This p-value is calculated to
	 based on a Khi2 test comparing the marker segregation to expected segregation. [Default: 0.0001]
	--pop: Population type (Possible values: SELFPOL, SELF, BiP). [Default: BiP]
		Possible values:
			BiP: bi-parental cross. Expected segregation tested: 0.5/0.5 (parental markers) and 0.25/0.5/0.25 (bridge markers).
			SELF: selfing population. Expected segregation tested: 0.25/0.5/0.25 (bridge markers).
	--prefix: Prefix for output files. [Default: Pop]
	--addcov: A tabulated file containing genotype of all markers passing filter is outputed.
	 If this option is passed, in addition to genotypes, alleles coverage information is also filled.
		Possible values:
			y: add this information
			n: do not add this information
		[default: n]
	--drawplot: Draw statistic plot (y or n).
		Possible values:
			y: add this information
			n: do not add this information
		[Default: n]
	--parent: (optional) Names of the parents of the population (separated by ":"). If passed,
	 these names will be used to parse marker depending of there segregation and heterozygosity
	 in the parents.
	--NoUsed: (optional) A tabultated file containing in one column, names of accessions to exclude from the
	 filtration (based on missing data and p-value) but which will be kept in final files.
	--exclude: (optional) A tabultated file containing in one column, names of accessions to exclude from the
	 analysis and the files.
	--ref: (optional) The reference fasta file. If passed, a tag associated to the marker will be
	 outpouted in a fasta file. This tag will contained 125 bases before the marker and 125 bases
	 after.
	--remove: (optional) For some programs, marker name length is limited. This option helps you to reduce marker
	 names. By default marker name is "chromosome name"+"M"+"site position". A string can be passed
	 that will be searched and removed from all marker name. This is not neccessary if your chromosome
	 name is not to long.

*Outputs:*\
 **\*\_JM\_Bridge.loc:** A .loc file that can be passed to joinmap that contained bridge markers.\
 **\*\_JM\_*Parent*.loc:** Two .loc file that can be passed to joinmap that contained parent1 and
parent2 markers respectively. Only if parent option is filled.\
 **\*\_JM\_unknown.loc:** A .loc file that can be passed to joinmap that contained unknown parent
markers (missing data for both parents). Only if parent option is filled.\
 **\*\_onemap\_Bridge.tab:** A .tab file that can be passed to onemap that contained bridge markers.\
 **\*\_onemap\_*Parent*.tab:** Two .tab file that can be passed to onemap that contained parent1 and
parent2 markers respectively. Only if parent option is filled.\
 **\*\_onemap\_unknown.tab:** A .tab file that can be passed to onemap that contained unknown parent
markers (missing data for both parents). Only if parent option is filled.\
 **\*\_tab\_Bridge.tab:** A .tab file correponding to a simplified joinmap format that contained bridge markers.\
 **\*\_tab\_*Parent*.tab:** Two .tab filecorreponding to a simplified joinmap format that contained parent1 and
parent2 markers respectively. Only if parent option is filled.\
 **\*\_tab\_unknown.tab:** A .tab file correponding to a simplified joinmap format that contained unknown parent
markers (missing data for both parents). Only if parent option is filled.\
 **\*\_report.tab:** A file report.\
 **\*\_sub.vcf:** A sub vcf corresponding to the original vcf with only lines corresponding toconserved markers
(no filtering applied in this vcf).\
 **\*\.tab:** A file containing for aech selected marker, the genotype of each accessions based on filter applied.
Two additional values are added at the end of the file: the Khi-Square value and the P-value of the test.\
 **\*\.pdf:** A pdf file containing various statistics on the vcf filtrations. Only if --drawplot=y.\
 **\*\_tags.fasta:** A fasta file containing marker tags (to align against another reference genome for example).
Only if --ref option is filled.
<br><br>
