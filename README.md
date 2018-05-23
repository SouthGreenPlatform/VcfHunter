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

1.  unzip the tar.gz file with the following command line : tar -xvzf
2.  open the loca\_programs.conf file
3.  set the path to each programs required
<br><br><br>

Dependencies
------------

1.  STAR, https://github.com/alexdobin/STAR
2.  PicarTools, https://broadinstitute.github.io/picard/
3.  GATK, https://software.broadinstitute.org/gatk/
4.  Samtools, https://github.com/samtools/samtools
5.  Bamtools, https://github.com/pezmaster31/bamtools
6.  bam-readcount, https://github.com/genome/bam-readcount
7.  gnuplot, http://www.gnuplot.info/
8.  circos-0.67 or greater, http://circos.ca/software/download/circos/

Python2, Python3, Java and Biopython are also required.
<br><br><br>

Description
-----------

The package provided comprised X programs listed here:

-   process\_RNAseq.1.0.py (python2)
-   process\_reseq.1.0.py (python2)
-   VcfPreFilter.1.0.py (python2)
-   vcf2struct.1.0.py (python3)
-   vcf2linear.1.0.py (python3)
-   haplo2kar.1.0.py (python3)
-   haplo2karByChr.1.0.py (python3)
-   haplo2Circos.1.0.py (python3)
-   vcfFilter.1.0.py (python3)
-   vcf2allPropAndCov.py (python3)
-   vcf2pop.1.0.py (python3)
-   vcf2popNew.1.0.py (python3)
-   RecombCalculatorDDose.py (python3)
-   Draw_dot_plot.py (python3)

All X programs run using the following command: ~\~~ python program-name
\<--options-name value\> ~\~~
<br><br><br>

Programs
--------
<br><br>

### process\_RNAseq.1.0.py (python2)

This program takes a reference DNA sequence multifasta file and several
fastq files and returns a bam file for each accessions and a final VCF
file containing alleles count at each variant site having at least one
variant allele supported by at least one read.

![](http://banana-genome-http.cirad.fr/image/Process_RNAseq.png)

#### Options:

    --conf: A configuration file containing path to references sequence (multifasta file) and RNAseq reads (fastq files).
    --thread: Number of processor to use (i.e. Max number of accessions treated at the same time). Do not exceed the number of processors available! [default: 1] 
    --queue: If you are using SGE sheduler: the queue name for the job to perform parallelisation. If not do not fill.
    --prefix: Prefix for vcf file and statistics folders.
    --steps: A string containing steps to perform:
        a: STAR indexing reference,
        b: Merging fastq from first mapping step to identify splicing sites,
        c: STAR Indexing reference with identified splicing sites,
        d: Second mapping step,
        e: Merging bam libraries with identical identifier,
        f: Removing duplicates,
        g: Reordering reads,
        h: Splitting and trimming reads,
        i: Indel realignment,
        j: Allele counting,
        k: Genotype calling,
        l: Merging genotype calling,
        m: Gene exon coverage statistics calculation.

#### Configuration file description:

The configuration file should contain 4 sections and must be formated as
followed:

    [Libraries]
    lib1 = genome_name path_to_mate1 path_to_mate2 ploidy
    lib2 = genome_name path_to_single ploidy
    ...
    [Reference]
    genome = path_to_the_reference_sequence
    [star]
    options = additional_options_to_pass
    [General]
    max_size = max_read_length
    gff3 (optional) = path to a gff3 file used to calculate statistics on genes coverage in step "m"

#### Output:

**Warning:** This program need to create a .dict and a .fai file un the
folder were the reference sequence is stored if they do not already
exist. Make sure that you have right to write in this folder!

Outputs are dependent of the steps you are running and each steps use
the output of the preceding one.

-   **step a:** generates a folder (--prefix + *ref*star\_1) in which
    the reference sequence is indexed,
-   **step b:** generates a file (--prefix +
    *JUNC*ESTIMATION\_SJ.out.tab) containing splicing sites detected by
    STAR on the complete dataset,
-   **step c:** generates a folder (--prefix + *ref*star\_2) in which
    the reference sequence is indexed with splicing sites detected,
-   **step d:** generates a folder for each accession, (names filled in
    column 3 "genome\_name") filled in the configuration file, which
    contained the sam files of aligned reads and a .final.out file of
    mapping statistice for each libraries. In addition a (--prefix)
    folder containing a mapping statistics file (--prefix + mapping.tab)
    for all accession is generated.
-   **step e:** generates a merged bam (\*\_merged.bam) and bai
    (\*\_merged.bai) files containing all reads of all libraries of the
    same accession,
-   **step f:** generates a bam (\**rmdup.bam) files for each accessions
    with duplicated reads removed. In addition, a file named (--prefix +
    rmdup*stat.tab) file containing duplicate statistics for each
    accessions was generated in the (--prefix) folder.
-   **step g:** generates a reodered (\*\_reorder.bam) and bai
    (\*\_reorder.bai) files for each accessions,
-   **step h:** generates a splitted and trimmed (on splicing sites) bam
    (\*\_trim.bam) and bai (\*\_trim.bai) files for each accessions,
-   **step i:** generates a bam (\*\_realigned.bam) and bai
    (\*\_realigned.bai) files realigned around indel for each
    accessions,
-   **step j:** generates for each accessions and each chromosomes a
    file (Accession + "*allele*count\_"+ chromosome + ".gz") counting
    variant at each covered bases,
-   **step k:** generates several vcf (one for each reference sequences)
    (--prefix + reference sequence + "*allele*count.vcf") file counting
    for each variant sites (at least on reads supporting a variant) and
    each accession the number of reads supporting each allele. For each
    accession in the vcf a genotype (GT tag) was called based on a
    binomial test, allelic depth was counted (AD tag) and total depth
    was rapported (DP tag),
-   **step l:** generates a uniq vcf file (--prefix +
    "*all*allele\_count.vcf") resulting in the concatenation of all vcf
    files generated at step k,
-   **step m:** Calculate exon coverage proportion of each accession for
    genes provided in the gff file
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

### vcf2struct.1.0.py (python3)

This program has been designed to perform statistics, filter, manipulate
vcf files as well as analysing the mosaique structure of genomes.

#### Mandatory Options:

    --type: A string corresponding to the type of analysis performed. 
    Possible values are: 
        RANDOM_SUB_SET 
        STAT
        FILTER
        COMPARE
        ADD_REF
        AL_IDENTITY
        FACTORIAL
        VISUALIZE_VAR_3D
        VISUALIZE_VAR_2D
        SNP_CLUST-Kmean
        SNP_CLUST-MeanShift
        FILTER_ON_MAX_GP_PROP
        MERGE_VCF
        ALL_PROP
        GET_GENOTYPE
        GET_GENOTYPE_AND_GROUP

#### Options and outputs depend on analysis type:

-   **RANDOM\_SUB\_SET:** Generate a vcf subset from the original vcf
    file in which variant line are sampled.\

*Options:*

    --vcf: A vcf file.
    --nRand: Number of variant site to get randomly from the vcf file (integer). [Default: 1000]
    --prefix: The prefix for output files. [Default: WorkOnVcf]

*Outputs:*\
 **\*\_subset.vcf:** A vcf file corresponding to a random line subset
of the first vcf

-   **STAT:** Calculate statistics on the vcf file.\

*Options:*

    --vcf: A vcf file.
    --names: A one column file containing accession names to treat.
    --prefix: The prefix for output files. [Default: WorkOnVcf]
    --gff3: (optional) A gff3 file containing gene annotation. If not filled, the statistics will return 0 for these coding sequence values.

*Outputs:*\
 **\*\_general.stat:** A file containing vcf global statistics such as
number of indel and SNP sites, number of different tags in the FILTER
columns, number of sites with 1,2,3,4, ... variants, number of
transitions, transversions and variant in annotated regions.\
 **\*\_accession.stat:** A file containing for each accessions missing
data number, number of alleles specific to this accession in the vcf,
number of homozygous sites identical to the reference, number of
homozygous sites different from the reference and number of heterozygous
sites.

-   **FILTER:** Filter a vcf file based on several parameters such as
    datapoint coverage and allele coverage, number of variant and
    variant type.\

*Options:*

    --vcf: A vcf file.
    --names: A one column file containing accession names to treat.
    --outgroup: (optional) A one column file containing accession names that will not be used for filtering but will remain in the
	 output file.
    --prefix: The prefix for output files. [Default: WorkOnVcf]
    --MinCov: Minimal coverage by accession to keep genotype calling (integer). If the value is lower, genotype will be converted to
	 unknown for the concerned accession. [Default: 10]
    --MinAl: Minimal allele coverage by accession to keep genotype calling (integer). If the value is lower for at least one allele,
	 genotype will be converted to unknown for the concerned accession. [Default: 3]
    --nMiss: Maximal number of missing genotype in a line to keep the line (integer). [Default: 0]
    --RmAlAlt: (optional) Number of alleles at the site in filtered accessions to remove the variant site (several values can be
	 passed and should be sepatated by ":"). Values: 1,2,3,...,n
    --RmType: (optional) Variant status to filter out (several values can be passed in this case they should be separated by ":"). 
        Possible values: 
            *Values which can be found in the FILTER column: PASS, DP_FILTER, QD_FILTER, SnpCluster, 
            *Other values: INDELS, SNP, AUTAPO (accession specific variant site).

*Outputs:*\
 **\*\_filt.vcf:** a filtered vcf file based on passed options.

-   **COMPARE:** Compare two variant accessions (from two different vcf
    files, in the same vcf file). This will output in standard output
    specific variant lines for accession1 and specific variant lines for
    accession2 as well as shared variant sites. In addition, identical
    variant calling will be counted and proportion among shared variant
    sites will be calculated.\

*Options:*

    --vcf: A vcf file.
    --comp1: Accession name to compare. If 2 vcf files are provided and only one accession name is passed, this name will be
	 searched in both vcf files.
    --vcf2: (optional if --comp2 filled) A second vcf file.
    --comp2:(optional if --vcf2 filled) Second accession to compare.
    If 2 vcf and 2 names are passed, comp1 will be searched in vcf and comp2 will be searched in vcf2.\

-   **ADD\_REF:** Add a haploid accession corresponding to the reference
    to the vcf called ref\_silico\_call.\

*Options:*

    --vcf: A vcf file.
    --prefix: The prefix for output files. [Default: WorkOnVcf]
    --ref_cov: Value to put to AD and DP flags (integer). [Default: 1000]

*Outputs:*\
 **\*\_add\_ref.vcf:** A new vcf file.

-   **AL\_IDENTITY:** Calculate genotype identity.\

*Options:*

    --vcf: A vcf file.
    --prefix: The prefix for output files. [Default: WorkOnVcf]

*Outputs:*\
 **\*\_ident.mat:** A matrix containing genotype pairwise identity.

-   **FACTORIAL:** Perform factorial analysis on the vcf file. This
    analysis is performed in several steaps.\

1- First the vcf file is recoded as followed: For each allele at
 each variants site two markers were generated; One marker for the
 presence of the allele (0/1 coded) and one for the absence of the
 allele (0/1 coded).\
 ![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig1.png)\
 Only alleles present or absent in **part** (not all) of selected
 accessions were included in the final matrix file.\
 If groups information was passed to the script, alleles groups were
 attributed based on the following rule: the allele is attributed to
 a group if it is only present in this group but not in other defined
 groups. If no grouping information the GROUP column is filled with
 UN value. This grouping value **doesn't have any influence on the
 analysis**, it only allows to add colors graphs drawn. It will also
 help to validate if the structure of your data correspond to the one
 you suspect.\
 2- The factorial analysis was performed on the transposed matrix
 using R. Graphical outputs of the analysis were draw and for example
 accessions and alleles can be projected along axis in the following
 picture.\
 ![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig2.png)\
 In this example allele projected along synthetic axis were colorated
 if group informations were passed to the program.\

*Options:*

    --vcf: A vcf file.
    --names: A one column file containing accession names to work with.
    --prefix: The prefix for output files. [Default: WorkOnVcf]
    --nAxes: Axis number to keep for the factorial analysis (integer). [Default: 4]
    --mulType: Multivariate analysis type. Possible values: coa, pca and pca_normed [Default: coa]
    --group: (optional) A file containing two sections: A section[group] with in col 1 accession name ; col 2 group (UN for
	 unknown group). All group should be in capital letters. A section [color], that define for each group a color for pca
	 drawing (in RGB+alpha percentage, ex: red=1:green=0:blue=0:alpha=0.1)
	Optional and dependant parameters:
	--dGroup: If passed, all alleles belonging to groups passed to this option will be removed.
	--mat: Matrix of grouped alleles (with either a GROUP or a K-mean_GROUP column). If a K-mean_GROUP column is found, the
	 filter will be performed on this column, else it will be performed on the GROUP one.
	
*Outputs:*\
 **\*\_axis*x*\_vs\_*y*\_accessions.pdf:** Several pdf files showing
accessions projected along x and y axis.\
 **\*\_axis*x*\_vs\_*y*.pdf:** Several pdf files showing accessions projected
along x and y synthetic axis in a first graphe and allele projection
along x and y synthetic axis.\
 **\*\_inertia.pdf:** A pdf files showing axis inertia.\
 **\*\_matri\_4\_PCA.tab:** A tabulated file of the recoded vcf file passed
to R for the factorial analysis.\
 **\*\_multivariate.R:** The R script file passed to R to do the
abalysis.\
 **\*\_multivariate.Rout:** The R script log file.\
 **\*\_individuals\_coordinates.tab:** A tabulated file of individuals
coordinates along synthetic axis.\
 **\*\_variables\_coordinates.tab:** A tabulated file of allele coordinates
along synthetic axis.\
 **\*\_variables\_coordinates\_scaled.tab:** A tabulated file of allele scaled
coordinates (colums centered and reduced) along synthetic axis.

-   **SNP\_CLUST-Kmean:** Perform a k-mean clustering of allele based on
    their coordinates on synthetic axis. It uses the k-mean algorithm of
    scikit-learn.\

*Options:*

    --VarCoord: The *_variables_coordinates.tab file generated by FACTORIAL.
    --dAxes: Axes to use in kmean clustering. Axis should be separated by ":".
    --nGroup: Group number for the k-mean algorithm that will cluster variables (alleles) based on their coordinates. [Default: 2]
    --mat: The *_matrix_4_PCA.tab tabulated file containing variant allele encoded generated by FACTORIAL.
    --prefix: The prefix for output files. [Default: WorkOnVcf]
    --thread: Number of processor available. [Default: 1]
    --iter: Parallele k-mean clustering different startpoints performed. [Default: 100]
	--AP: Cluster absent (A) and present (P) lines. Possible values, "y" or "n" [Default: y]

*Outputs:*\
 **\*\_centroid\_coordinates.tab:** Coordinates of the distincts final
centroides calculated for the X k-mean random starpoints.\
 **\*\_centroid\_iteration\_grouping.tab:** A tabulated file indicating
centroid group attribution.\
 **\*\_kMean\_allele.tab:** A tabulated file equivalent to
 **\*\_matrix\_4\_PCA.tab** in which an column (named K-mean\_GROUP) containing
k-mean allele grouping has been added.\
 **\*\_group\_color.tab:** Color file in which a RGB color as been
attributed to each k-mean cluster group. This file can be used for
VISUALIZE\_VAR\_2D and VISUALIZE\_VAR\_3D tools.\
 **\*\_kMean\_gp\_prop.tab:** A tabulated file reporting for each allele the
probability to be in each groups. This is not a "real" probability, the
idea was to have a statistics in case you want to filter alleles. This
value was calculated as the inverse of the euclidian distance of one
point and each centroids and these values were normalized so that the
sum is equal to 1.

-   **SNP\_CLUST-MeanShift:** Perform a clustering using the MeanShift
    algorithm of scikit-learn of allele based on their coordinates on
    synthetic axis.\

*Options:*

    --VarCoord: The *_variables_coordinates.tab file generated by FACTORIAL.
    --dAxes: Axes to use in mean shift clustering. Axis should be separated by ":".
    --quantile: The quantile value to estimate de bandwidth parameters used in the MeanShift. Value should be in [0:1]. [Default: 0.2]
    --mat: The *_matrix_4_PCA.tab tabulated file containing variant allele encoded generated by FACTORIAL.
    --prefix: The prefix for output files. [Default: WorkOnVcf]
    --thread: Number of processor available. [Default: 1]
	--MeanShiftAll: Cluster all point in the MeanShift. Possible values, "y" or "n" [Default: y]
	--bandwidth: Bandwidth value used for mean shift. If filled, the --quantile parameter is ignored. [Default: None]
	--AP: Cluster absent (A) and present (P) lines. Possible values, "y" or "n" [Default: y]

*Outputs:*\
 **\*\_centroid\_coordinates.tab:** Coordinates of centroides calculated.\
 **\*\_centroid\_iteration\_grouping.tab:** A tabulated file indicating
centroid group attribution.\
 **\*\_kMean\_allele.tab:** A tabulated file equivalent to
**\*\_matrix\_4\_PCA.tab** in which an column (named K-mean\_GROUP) containing
k-mean allele grouping has been added.\
 **\*\_group\_color.tab:** Color file in which a RGB color as been
attributed to each k-mean cluster group. This file can be used for
VISUALIZE\_VAR\_2D and VISUALIZE\_VAR\_3D tools.\
 **\*\_kMean\_gp\_prop.tab:** A tabulated file reporting for each allele the
probability to be in each groups. This is not a "real" probability, the
idea was to have a statistics in case you want to filter alleles. This
value was calculated as the inverse of the euclidian distance of one
point and each centroids and these values were normalized so that the
sum is equal to 1.

-   **VISUALIZE\_VAR\_2D:** Perform plots of alleles projected along
    synthetic axes.\

*Options:*

    --VarCoord: The *_variables_coordinates.tab file generated by FACTORIAL.
    --dAxes: Axes to plot. Axis should be separated by ":". All conbination between axis will be drawn.
    --mat: The *_matrix_4_PCA.tab or *_kMean_allele.tab.
    --group: (optional) A file containing at least the a section [color], that define for each group a color
	 (in RGB+alpha percentage, ex: red=1:green=0:blue=0:alpha=0.1). This file can be the *_group_color.tab
	 generated by SNP_CLUST.
    --dGroup: Groups IDs to draw. Groups names should be separated by ":". Groups names will be searched in
	 the file provided in --mat argument in the "K-mean_GROUP" column, and if such column is not present, they
	 will be searched in "GROUP" columns. If none of these columns are found the program will exit without finishing.
    --prefix: The prefix for output files. [Default: WorkOnVcf]

*Outputs:*\
 \***axis*x*\_vs\_axis*y*.png:** Several png files showing
grouped alleles projected along x and y synthetic axes.

-   **VISUALIZE\_VAR\_3D:** Perform an interactive 3d plots of alleles
    projected along 3 synthetic axes.\

*Options:*

    --VarCoord: The *_variables_coordinates.tab file generated by FACTORIAL.
    --dAxes: 3 axes to plot. Axis should be separated by ":". All conbination between axis will be drawn.
    --mat: The *_matrix_4_PCA.tab or *_kMean_allele.tab.
    --dGroup: Groups IDs to draw. Groups names should be separated by ":". Groups names will be searched in the
	 file provided in --mat argument in the "K-mean_GROUP" column, and if such column is not present, they will
	 be searched in "GROUP" columns. If none of these columns are found the program will exit without finishing.
    --group: (optional) A file containing at least the a section [color], that define for each group a color
	 (in RGB+alpha percentage, ex: red=1:green=0:blue=0:alpha=0.1).

-   **FILTER\_ON\_MAX\_GP\_PROP:** Filter the matrix file by removing
    ambiguous grouped alleles. *i.e.* alleles which changed of groups
    during the k-mean distinct attempts.\

*Options:*

    --gpPropFile: The group file proportion, *_kMean_gp_prop.tab file generated by SNP_CLUST.
    --mat: The *_kMean_allele.tab file generated by SNP_CLUST.
    --gpPropValue: Minimal value to keep the allele grouped. This value should be comprised between 0 and 1 [Default: 0.95]
    --prefix: The prefix for output files. [Default: WorkOnVcf]

*Outputs:*\
 **\*\_kMean\_allele\_filtered\_with\_*x*\_value.tab:** A filtered
file of \***\_kMean_allele.tab** file in which ambiguous alleles were
removed.

-   **MERGE\_VCF:** Add an accession from a vcf to a second one vcf
    file. If the variant line is absent for the added accession, missing
    genotype if filled. If a variant line is present for the accession
    but absent in the vcf, this line will not be added.\

*Options:*

    --vcf: A vcf file in which you want to add an accession.
    --vcf2: A vcf file containing accession you want to add to the fist vcf file.
    --comp1: Accession name to add to the first vcf.
    --prefix: The prefix for output files. [Default: WorkOnVcf]

*Outputs:*\
 **\*\_merged.vcf:** A new vcf file with the accessions variant
calling added.

-   **ALL\_PROP:** Calculate for each accessions the number of grouped
    allele for each groups.\

*Options:*

    --vcf: A vcf file.
    --mat: A matrix file containing allele grouping with either "K-mean_GROUP" or "GROUP" columns, if "K-mean_GROUP" does not
	 exists the "GROUP" column will be used. These files are those generated by FACTORIAL (*_matrix_4_PCA.tab), SNP_CLUST (*_kMean_allele.tab)
	 or  FILTER_ON_MAX_GP_PROP (*_kMean_allele_filtered_with_x_value.tab). 
    --prefix: The prefix for output files. [Default: WorkOnVcf]
    --names: (optional) A one column file containing accession names to calculate proportion on (it can be accessions that were
	 not used for PCA analysis but originating from the same original vcf).
    --dGroup: (optional) Groups IDs to calculate proportion on. Groups names should be separated by ":". Groups names will be searched
	 in the file provided in --mat argument in the "K-mean_GROUP" column, and if such column is not present, they will be searched in
	 "GROUP" columns. If none of these columns are found the program will exit without finishing.
    --ExclChr: (optional) A list of chromosome to exclude from the analysis separated by ":".

*Outputs:*\
 **\*\_gp\_prop.tab:** A tabulated file containing for each accessions the
number of grouped allele for each groups.

-   **GET\_GENOTYPE:** Sometimes as a biologist you want to see the data
    and more precisely have a look at the genotypes (in term of A,T,G,C)
    along chromosomes. **This tool if for you!**\

*Options:*

    --vcf: A vcf file.
    --names: A one column file containing accession names to plot genotypes
    --prefix: The prefix for output files. [Default: WorkOnVcf]

*Outputs:*\
 **\*\_genotype.gen:** A tabulated file containing, for each
accessions (columns) and at each positions (rows) along chromosomes, the
genotype.

-   **GET\_GENOTYPE\_AND\_GROUP:** Sometimes as a biologist you want to
    see the data and more precisely have a look at the genotypes (in
    term of A,T,G,C)and allele grouping along chromosomes. **This tool
    if for you!**\

*Options:*

    --vcf: A vcf file.
    --names: A one column file containing accession names to plot genotypes
    --mat: A matrix file containing allele grouping with either "K-mean_GROUP" or "GROUP" columns, if "K-mean_GROUP" does not
	 exists the "GROUP" column will be used. These files are those generated by FACTORIAL (*_matrix_4_PCA.tab), SNP_CLUST (*_kMean_allele.tab)
	 or  FILTER_ON_MAX_GP_PROP (*_kMean_allele_filtered_with_x_value.tab).
    --dGroup: Groups IDs to calculate proportion on. Groups names should be separated by ":". Groups names will be searched
	 in the file provided in --mat argument in the "K-mean_GROUP" column, and if such column is not present, they will be searched
	 in "GROUP" columns. If none of these columns are found the program will exit without finishing.
    --prefix: The prefix for output files. [Default: WorkOnVcf]

*Outputs:*\
 **\*\_genotype.gen:** A tabulated file containing, for each
accessions (2 columns per accessions) and at each positions (rows) along
chromosomes, the genotype (first of the two columns) and allele grouping
(second of the two columns).
<br><br>

### vcf2linear.1.1.py (python3)

This programs aims at performing a chromosome painting of accessions
along chromosome based on the allele grouping correponding to the
ancestral groups. The idea is to search, for distinct regions (along the
chromosomes) of the accession we want to study, the probability of an
ancestral origin (homozygous or heterozygous). Because, distinct
ancestral groups can contain different specific alleles, show distinct
heterozygosity levels and allele fixation (du to distinct reproductive
mode for example) such probability is not straightforward. An idea can
be to calculate for each regions, the expected ancestral allele grouping
proportion (based on ancestral accessions) and compare this value to the
observed one. Because dataset can be small and/or representatives of
ancestral group are very few, the expected ancestral allele grouping
proportion can be difficult to estimate.\
 This programs solves the problem by simulating ancestral population
from accessions identified as representatives of the ancestral groups
under panmixie hypothesis. For each ancestral populations a total of 100
individuals were simulated. The same was performed for hybrids between
two populations. At the and of this simulation process mean and standard
deviation of allele grouping proportions were calculated for each
ancestral group on sliding windows of size n overlapping of n-1. Because
the simulation process and access to the simulation is a long process,
we also implemented a function that estimated grouping proportions based
on sums of binomial densities. The mean values estimated are exact but
standard deviations are approximated.\
 These grouping proportions were then calculated on tested accessions
and estimated proportion were used to calculate the ancestry probability
as followed.\

With:\
 ***Hmu***: Mean gX expected number for the group gX based on
simulations on homozygous groups,\
 ***Hsd***: Maximal standard deviation value calculated for all groups
based on simulations on homozygous groups,\
 ***Emu***: Mean gX expected number for the group gX based on
simulations on heterozygous groups,\
 ***Esd***: Maximal standard deviation value calculated for all groups
based on simulations on heterozygous groups,\
 ***Nmu***: Mean noisy allele grouped based on simulations (calculated
as the number of alleles not from the right group according to the
simulated population),\
 ***Nsd***: Noise standard deviation,\
 ***obs***: gX number observed for the accession on the window,\
 ***PgX***: Probability to be at least gX,\
 ***PgN***: Probability to be in the noise (unknown/additionnal
ancestor),\
 ***N(x, mu, sd)***: Probability density function of the normal
distribution for mean = *mu* and standard deviation = *sd*.

-   probability to be homozygous for group gX:\
    -   if *obs* ***superior or equal*** (*Hmu* - *Hsd*): *PgX* = 1\
    -   if *obs* ***inferior*** (*Hmu* - *Hsd*) : *PgX = N(obs, Hmu-Hsd,
        Hsd)/N(Hmu-Hsd, Hmu-Hsd, Hsd)*\
-   probability to be heterozygous for group gX:\
    -   if *obs* ***superior or equal*** (*Emu* - *Esd*): *PgX* = 1\
    -   if *obs* ***inferior*** (*Emu* - *Esd*) : *PgX = N(obs, Emu-Esd,
        Esd)/N(Emu-Esd, Emu-Esd, Esd)*\
-   probability to be in the noise:\
    -   if *obs* ***inferior or equal*** (*Nmu* + *Nsd*): *PgN* = 1\
    -   if *obs* ***superior*** (*Nmu* + *Nsd*): *PgN = N(obs, Emu+Esd,
        Esd)/N(Emu+Esd, Emu+Esd, Esd)*\

Genotype grouping at the window was then attributed based on the maximal
probability and haplotypes representation was performed trying to
minimize the recombination events.

*Options:*

    --vcf: A vcf file.
    --mat: A matrix file containing allele grouping with either "K-mean_GROUP" or "GROUP" columns, if "K-mean_GROUP" does not
	 exists the "GROUP" column will be used. These files are those generated by FACTORIAL (*_matrix_4_PCA.tab), SNP_CLUST (*_kMean_allele.tab)
	 or  FILTER_ON_MAX_GP_PROP (*_kMean_allele_filtered_with_x_value.tab).
    --names: A one column file containing accession names to treat (it is not necessarilly thoses used in the Factorial analysis
	but they must share variant sites).
    --namesH: A two column file containing accession names used to simulate populations. 
    --chr: Chromosomes names to work with (each chromosome names should be separated by ":"). If not filled, all chromosomes will be used.
    --win: Half window size (allele number) around a variant site to evaluate the structure at the site (integer). [Default: 25]
    --ploidy: Ploidy level (integer). [Default: 2]
    --thread: Number of processors to use (integer), [default: 1]
    --type: Type of estimation performed: "Simul", "Binom". If "Simul", a total of 100 individuals are simulated for
     each combinations of haplotype and mean values and sd values are estimated based on these simulation. If
     "Binom", mean value is calculated as the sum of binomial mean at each point (exact estimator) and sd
     value is estimated as sqrt(sum variance at each point). This is not the exact sd but the analysis is a
     lot more faster! If you do not trust this sd estimation, you can choose to change this estimator
     by filling a value between ]0,1] to --prop parameter. In this case, the program will use the
     maximal_expected_value*prop_argument instead of using the maximal sd observed for all groups for probability
     calculation [default: Binom]
    --prop: Estimator different from sd calculated as mean_value*--prop. Value should be comprised in ]0,1].
     A value of 0, means that this parameter is not used. [default: 0]
    --gcol: A file containing at least the a section [color], that define for each group a color (in RGB+alpha percentage,
	ex: red=1:green=0:blue=0:alpha=0.1). This file can be the *_group_color.tab generated by SNP_CLUST.
	--Ambiguous: A 2 column file containing accession names used to simulate populations and their group.
	 In this file we only pass admixed representative of a group not represented in nameH. These accessions will
	 be used to infer expected allele of a group at the observed position. This inference is calculated as followed:
	 a probability to have an allele of this group in an accession of this group is of 1 an allele of the group is
	 found in the passed accession. Accessions passed to this option are not used for noise inference.
	--FormerFolder: Folder containing previous chromosome painting. This chromosome painting will be used to better
	estimate the expected number of alleles based on admixed accessions representative of a group. If only the ancestral
	 group is found in the haplotypes of the admixed accessions (passed to --Ambiguous), the probability is calculated has
	 the number of alleles of the groupe observed at this position. In others cases, the probability is still 1 if an alleles
	 of the group is observed.
	--sdMult: Multiplicator of standard deviation for probability calculation of a segment to be of a group. [Default: 1]
    --prefix: The prefix for output files. [Default: WorkOnVcf]

*Output*\

A folder with the name passed in --prefix options which contained
several files for each accessions and each chromosomes:\

**\*_Accession\_chromosome.tab:** A tabulated file with for each window
around a given position:

-   1- the count of each grouped alleles *(as much columns as ancestral
    groups)*,
-   2- the expected number of grouped alleles, for a given group, in
    case of single ancestral origin (homozygous) *(as much column as
    ancestral groups)*,
-   3- the expected number of grouped alleles in case of single
    ancestral origin (homozygous) *(as much columns as ancestral
    groups)*,
-   4- the corresponding expected standard deviation *(as much column as
    ancestral groups)*,
-   5- the retained maximal standart deviation for only one ancestry,
-   6- the expected number of grouped alleles, for a given group, in
    case of two ancestral origin (heterozygous) *(as much columns as
    ancestral groups)*,
-   7- the expected number of grouped alleles in case of two ancestral
    origin (heterozygous) *(as much columns as ancestral groups)*,
-   8- the corresponding expected standard deviation *(as much columns
    as ancestral groups)*,
-   9- the retained maximal standart deviation for two distinct ancestry
    ,
-   10- the probability to have a single origin for the group X *(as
    much columns as ancestral groups)*,
-   11- the probability to have a group X origin and a second origin
    *(as much columns as ancestral groups)*,
-   12- the heterozygosity level (in the studied window),
-   13- expected noise count for each groups *(as much columns as
    ancestral groups)*,
-   14- expected noise standard deviation *(as much columns as ancestral
    groups)*,
-   15- the retained maximal noise standart deviation,
-   16- accession noise probabilities *(as much columns as ancestral
    groups)*.\

**Accession\_chromosome\_haplo1.tab:** A tabulated file listing grouped
region bloacs for haplotype1\
 **Accession\_chromosome\_haplo2.tab:** A tabulated file listing grouped
region bloacs for haplotype2\
 **Accession\_chromosome\_density.pdf:** A pdf file summarysing all
statistics (heterozygosity, expected and observed ancestries along
chromosomes, accessions ancestry probabilities and infered accession
haplotypes)
<br><br>

### haplo2kar.1.0.py

This program perform a synthesis of the haplotypes reconstructed by
vcf2linear.1.0.py by drawing for an accession the chromosomal painting
for all its chromosomes.

*Options:*

    --acc: Accession name.
    --chr: Chromosomes list to draw (separated by ":")
    --gcol: A file containing at least the a section [color], that define for each group a color (in RGB+alpha percentage,
	 ex: red=1:green=0:blue=0:alpha=0.1). This file can be the *_group_color.tab generated by SNP_CLUST.
    --dg: Groups to draw, separated by ":". Groups not passed here will be replaced by grey.
    --centro: A tabulated file locating pericentromeric regions. In column 1: chromosome name, column 2: start, column 3: end.

*Output:*\
 **Accession.pdf:** A pdf file with chromosome painting for all
chromosomes of this accession.
<br><br>

### haplo2karByChr.1.0.py

This program perform a synthesis of the haplotypes reconstructed by
vcf2linear.1.0.py by drawing for a chromosome the chromosomal painting
for all selected accessions.

*Options:*

    --acc: Accession names configuration file with Column1: accession name and column2: ploidy level.
    --chr: Chromosome to draw
    --gcol: A file containing at least the a section [color], that define for each group a color (in RGB+alpha percentage,
	 ex: red=1:green=0:blue=0:alpha=0.1). This file can be the *_group_color.tab generated by SNP_CLUST.
    --dg: Groups to draw, separated by ":". Groups not passed here will be replaced by grey.
    --centro: A tabulated file locating pericentromeric regions. In column 1: chromosome name, column 2: start, column 3: end.
    --prefix: Prefix for output files. [Default: All_acc]
    --maxChr: Maximal haplotype number to draw in a pdf. [Default: 53]

*Output:*\
 **\*\_chromosme\_X.pdf:** As much pdf file as necessary with
chromosome painting for all chromosomes of this accession.
<br><br>

### haplo2Circos.1.0.py

This program perform a synthesis of the haplotypes reconstructed by
vcf2linear.1.0.py by drawing a circos representation of the chromosomal
painting for all selected accessions and chromosomes.

*Options:*

    --acc: Accession names configuration file with Column1: accession name and column2: ploidy level.
    --chr: Chromosome(s) to draw, separated by ":"
    --gcol: A file containing at least the a section [color], that define for each group a color (in RGB+alpha percentage,
	 ex: red=1:green=0:blue=0:alpha=0.1). This file can be the *_group_color.tab generated by SNP_CLUST.
    --dg: Groups to draw, separated by ":". Groups not passed here will be replaced by grey.
    --centro: A tabulated file locating pericentromeric regions. In column 1: chromosome name, column 2: start, column 3: end.
    --prefix: Prefix for output files. [Default: All_acc]

*Output:*\
 **\*.conf**: the configuration file used by circos.
I choose to keep this file so that it can be edited to change aspect of
the Figure if you want. After editing this file you will only have to
run the following command line to have your new Figure : *circos -conf Circos\_All\_Admix.conf -noparanoid*\
**\*\_housekeeping.conf**: a second file used by circos\
**\*\.kar**: a third file also used by circos\
**\*\.png**: the circos Figure
<br><br>

### vcfFilter.1.0.py

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
 
### vcf2allPropAndCov.py

This program perform two things based on a vcf. 1) It plots for an
accession, the allele coverage alongs its chromosomes. 2) It identify,
based on known ancestral accessions in the vcf, the alleles specific to
each groups and plot the alleles proportion at a site in the accession
along chromosomes.

*Options:*

    --conf: Conf file containing vcf location (one per chromosome or a single vcf for all chromosomes),
	--origin: A 2 column file containing accession name (col1), origin/group (Col2),
	--acc: Accession to work with,
	--ploidy: Accession ploidy (integer),
	--all: Allele should be present in all accessions to attribute to a group,
	--prefix: Prefix for output files. [Default: AlleleOrigin]

*Output:*\
 **\*\Cov.png:** a png file presenting SNP coverage along the chromosomes.\
 **\*\Ratio.png:** a png file presenting ancestral allele proportion at a site along the chromosomes.\
 **\*\_AlleleOriginAndRatio.tab:** a tabulated file repporting for each sites were an ancestral allele
has been attributed, its origin and the proportion of reads supporting this allele. This files contains
chromosome (col1), position (col2), allele (col3), ancestral origin (col4) and allele ratio (col5).\
 **\*\_stats.tab:** a tabulated file repporting various statistics on the alleles of the accession.
<br><br>

### vcf2allPropAndCovByChr.py

This program perform two things based on a vcf. 1) It plots for a chromosome of all accessions in a vcf,
the allele coverage alongs its chromosomes. 2) It identify,based on known ancestral accessions in the vcf,
the alleles specific to each groups and plot the alleles proportion at a site along chromosomes for all
accessions.

*Options:*

    --conf: Conf file containing vcf location (one per chromosome or a single vcf for all chromosomes),
	--origin: A 2 column file containing accession name (col1), origin/group (Col2),
	--ploidy: Accession ploidy (integer). If not all accessions have the same ploidy, this is not a problem. This ploidy information
	 is only used to draw vertical lines in the coverage plot that help to identify ploidy change,
	--all: Allele should be present in all accessions to attribute to a group,
	--prefix: Prefix for output files. [Default: RatioAndCov]

*Output:*\
 **prefix\_chromosomeN\_X\_Cov.png:** X png files presenting SNP coverage along chromosomeN.\
 **prefix\_chromosomeN\_X\_Ratio.png:** X png files presenting ancestral allele proportion at a site along chromosomeN.
<br><br>


### vcf2pop.1.0.py

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


### vcf2popNew.1.0.py

This program will select markers for genetical mapping
analysis from a vcf file based on several criterias including segregation
 ratio. It will outpout coded markers as requested by the user.

*Options:*

    --vcf: The vcf file
	--seg: Segregation tested. Several segregations can be passed and should be separated by "/".
	 A segregation should look like as follows: Name:Parents:MarkerCoding:MarkerSegregation:PvalueForTest.
	 With a real example: SimpleDose:P1,P2:Ho,He@nn,np:0.5,0.5:1e-5/Bridge:P1,P2:Ho,He,Ho@hh,hk,kk:0.25,0.5,0.25:1e-5
	 (Ho for homozygous, He for heterozygous)
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
	--prefix: Prefix for output files. [Default: Pop]
	--addcov: A tabulated file containing genotype of all markers passing filter is outputed.
	 If this option is passed, in addition to genotypes, alleles coverage information is also filled.
		Possible values:
			y: add this information
			n: do not add this information
		[default: n]
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
 **\*\_tab\_Bridge.tab:** A .tab file correponding to a simplified joinmap format that contained bridge markers.\
 **\*\_tab\_*SegregationName*\_*Parent*.tab:** Two .tab filecorreponding to a simplified joinmap format that contained parent1 and
parent2 markers respectively. Only if parent option is filled.\
 **\*\_tab\_unknown.tab:** A .tab file correponding to a simplified joinmap format that contained unknown parent
markers (missing data for both parents). Only if parent option is filled.\
 **\*\_report.tab:** A file report.\
 **\*\_sub.vcf:** A sub vcf corresponding to the original vcf with only lines corresponding toconserved markers
(no filtering applied in this vcf).\
 **\*\.tab:** A file containing for aech selected marker, the genotype of each accessions based on filter applied.
Two additional values are added at the end of the file: the Khi-Square value and the P-value of the test.\
 **\*\_tags.fasta:** A fasta file containing marker tags (to align against another reference genome for example).
Only if --ref option is filled.
<br><br>


### RecombCalculatorDDose.py

This program perform is designed to calculate frequencies of recombination observed between two pairs of markers.
It can also calculate marker segregation distortion.

*Options:*

    --matrix: The marker file matrix (output of vcf2pop.1.0.py or vcf2popNew.1.0.py),
	--output: Prefix for output files,
	--phased: Are marker phased,
		Possible values:
			y: markers are phased
			n: markers are not phased
		[default: n]
	--steps: Analysis to perform,
		Possible values:
			R: Calculate recombination rate
			S: Calculate segregation distortions

*Output:*\
 **\_REC.tab:** A tabulated file of pairwise marker recombination rate (--setp R).\
 **\_SegDist.tab:** A tabulated file of marker segragtion distortions (--setp S).
<br><br>


### Draw_dot_plot.py

This program draw a dotplot based on marker pairwise recombination file obtained from ***RecombCalculatorDDose***.

*Options:*

    --matrix: The pairwise matrix marker file (generated by RecombCalculatorDDose.py),
	--loc: Loci to plot with their locations (col1: marker name, col2: chromosome, col3: position),
	--chr: (Optional) List of chromosomes to draw in this order (separated by ":"),
	--agp: (Optional) Agp file locating scaffolds in the reference sequence,
	--stat: (Optional) A two column file with column 1: marker name, column2: statistics,
	--phys: A value specifying if the marker position should be defined based on physical position or not.
		Possible
			y: yes
			n: no: markers are ordered along chromosomes with positions incremented by 1 
		[default: y]
	--output: Name of output file. Extension (.png, .svg, .pdf) of the outpout is autaomatically identified by the script.

*Output:*\
 A dotplot file representing pairwise marker linkage named according to --outpout option.\
 A heatmap file representing pairwise statistic coding named according to --outpout option.
<br><br>








