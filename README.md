# Purpose of VcfExplorer

VcfExplorer regroups several programs which principale aims are to map RNAseq data onto reference genome sequence, perform variant calling, manipulate vcf files and perform chromosome painting of accessions based on the contribution of ancestral groups.


## Installation

All proposed tools described here are written in python and work on linux system To install the tools:

1. unzip the tar.gz file with the following command line : tar -xvzf
2. open the loca_programs.conf file
3. set the path to each programs required

## Dependencies

1. STAR, https://github.com/alexdobin/STAR
2. PicarTools, https://broadinstitute.github.io/picard/
3. GATK, https://software.broadinstitute.org/gatk/
4. Samtools, https://github.com/samtools/samtools
5. Bamtools, https://github.com/pezmaster31/bamtools
6. bam-readcount, https://github.com/genome/bam-readcount
7. gnuplot, http://www.gnuplot.info/
8. circos-0.67 or greater, http://circos.ca/software/download/circos/

Python2, Python3, Java and Biopython are also required.

## Description

The package provided comprised X programs listed here:

* process_RNAseq.1.0.py (python2)
* process_reseq.1.0.py (python2)
* VcfPreFilter.1.0.py (python2)
* vcf2struct.1.0.py (python3)
* vcf2linear.1.0.py (python3)
* haplo2kar.1.0.py (python3)
* haplo2karByChr.1.0.py (python3)
* haplo2Circos.1.0.py (python3)

All X programs run using the following command:
~~~
python program-name <--options-name value>
~~~


## Programs


### process_RNAseq.1.0.py (python2)

This program takes a reference DNA sequence multifasta file and several fastq files and returns a bam file for each accessions and a final VCF file containing alleles count at each variant site having at least one variant allele supported by at least one read. 


<img src="http://banana-genome-http.cirad.fr/image/Process_RNAseq.png" height="600" width="450">


#### Options:
```
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
```

#### Configuration file description:

The configuration file should contain 4 sections and must be formated as followed:
```
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
```

#### Output:

**Warning:** This program need to create a .dict and a .fai file un the folder were the reference sequence is stored if they do not already exist. Make sure that you have right to write in this folder!

Outputs are dependent of the steps you are running and each steps use the output of the preceding one.

* **step a:** generates a folder (--prefix + _ref_star_1) in which the reference sequence is indexed,
* **step b:** generates a file (--prefix + _JUNC_ESTIMATION_SJ.out.tab) containing splicing sites detected by STAR on the complete dataset,
* **step c:** generates a folder (--prefix + _ref_star_2) in which the reference sequence is indexed with splicing sites detected,
* **step d:** generates a folder for each accession, (names filled in column 3 "genome_name") filled in the configuration file, which contained the sam files of aligned reads and a .final.out file of mapping statistice for each libraries. In addition a (--prefix) folder containing a mapping statistics file (--prefix + mapping.tab) for all accession is generated.
* **step e:** generates a merged bam (*_merged.bam) and bai (*_merged.bai) files containing all reads of all libraries of the same accession,
* **step f:** generates a bam (*_rmdup.bam) files for each accessions with duplicated reads removed. In addition, a file named (--prefix + rmdup_stat.tab) file containing duplicate statistics for each accessions was generated in the (--prefix) folder.
* **step g:** generates a reodered (*_reorder.bam) and bai (*_reorder.bai) files for each accessions,
* **step h:** generates a splitted and trimmed (on splicing sites) bam (*_trim.bam) and bai (*_trim.bai) files for each accessions,
* **step i:** generates a bam (*_realigned.bam) and bai (*_realigned.bai) files realigned around indel for each accessions,
* **step j:** generates for each accessions and each chromosomes a file (Accession + "_allele_count_"+ chromosome + ".gz") counting variant at each covered bases,
* **step k:** generates several vcf (one for each reference sequences) (--prefix + reference sequence + "_allele_count.vcf") file counting for each variant sites (at least on reads supporting a variant) and each accession the number of reads supporting each allele. For each accession in the vcf a genotype (GT tag) was called based on a binomial test, allelic depth was counted (AD tag) and total depth was rapported (DP tag),
* **step l:** generates a uniq vcf file (--prefix + "_all_allele_count.vcf") resulting in the concatenation of all vcf files generated at step k,
* **step m:** Calculate exon coverage proportion of each accession for genes provided in the gff file


### process_reseq.1.0.py (python2)

This program takes a reference DNA sequence multifasta file and several fastq files and returns a bam file for each accessions and a final VCF file containing alleles count at each variant site having at least one variant allele supported by at least one read. 


<img src="http://banana-genome-http.cirad.fr/image/Process_ReSeq_Fig1.png" height="600" width="450">


#### Options:

```
--conf: A configuration file containing path to references sequence (multifasta file) and RNAseq reads (fastq files).
--thread: Max number of accessions treated at the same time. Do not exceed the number of processors available! [default: 1] 
--queue: If you are using SGE sheduler: the queue name for the job to perform parallelisation. If not do not fill.
--prefix: Prefix for vcf file and statistics folders.
--steps: A string containing steps to perform:
	a: Aligning libraries
    b: Removing duplicates
    c: Indel realignment
    d: Bases recalibration
    e: Allele counting
    f: Genotype calling
    g: Merging genotype calling
    h: Mapping statistics calculation
```

#### Configuration file description:

The configuration file should contain 4 sections and must be formated as followed:

```
[Libraries]
lib1 = genome_name path_to_mate1 path_to_mate2 ploidy
lib2 = genome_name path_to_single ploidy
...
[Reference]
genome = path_to_the_reference_sequence
```

#### Output:

**Warning:** This program need to create a .dict and a .fai file un the folder were the reference sequence is stored if they do not already exist. Make sure that you have right to write in this folder!

Outputs are dependent of the steps you are running and each steps use the output of the preceding one.

* **step a:** generates a folder for each accession, (names filled in column 3 "genome_name") filled in the configuration file, which contained the bam (*_merged.bam) and bai (*_merged.bai) files of aligned reads, a .stat file generated for each libraries with **samtools stat** program and a STAT folder containing a html files summarising mapping statistics,
* **step b:** generates a bam (*_rmdup.bam) and bai (*_rmdup.bai) files for each accessions with duplicated reads removed. In addition in each folder duplicate statistics wer recorder in a file named *_duplicate,
* **step c:** generates a bam (*_realigned.bam) and bai (*_realigned.bai) files realigned around indel for each accessions,
* **step d:** generates a ,
* **step e:** generates for each accessions and each chromosomes a file (Accession + "_allele_count_"+ chromosome + ".gz") counting variant at each covered bases,
* **step f:** generates a generates several vcf (one for each reference sequences) (--prefix + reference sequence + "_allele_count.vcf") file counting for each variant sites (at least on reads supporting a variant) and each accession the number of reads supporting each allele. For each accession in the vcf a genotype (GT tag) was called based on a binomial test, allelic depth was counted (AD tag) and total depth was rapported (DP tag),
* **step g:** generates a generates a uniq vcf file resulting in the concatenation of all vcf files generated at step g,
* **step h:** generates two files (*_acc.stats and *_lib.stats) collecting mapping statistics on each libraries and accessions respectively,


### VcfPreFilter.1.0.py (python2)

This script filter VCF file generated by Process_RNAseq.1.0.py by removing homozygous sites for all accessions based on accession minimal and maximal coverage, accession minimal allele coverage and frequency parameters passed. The idea of this tool is to filter out variant lines resulting from sequencing errors. Filter are applicated as followed:

* only data points covered by at least "minimal coverage" reads were considered,
* only data points covered by at most "maximal coverage" reads were considered,
* only variant alleles supported by at least "minimal allele coverage" reads and having a frequency equal or greater to "minimal frequency" in one accession were kept as variant,
* sites showing at least one variant allele were kept for variant calling.
For each accession and at each variant site, a genotype was called based on the maximum likelihood of all possible genotype calculated based on a binomial distribution assuming a sequencing error rate of 0.005. The variant calling file was formated to vcf format.


#### Options:

```
--vcf: The VCF file
--MinCov: Minimal read coverage for site. [Default: 10]
--MaxCov: Maximal read coverage for site. [Default: 1000]
--minFreq: Minimal allele frequency in an accession to keep the allele for calling in the row
--MinAlCov: Minimal read number of minor allele to call variant heterozygous (between 1 and infinity). [Default: 3]
--out: Prefix for output files. [Default: Pop]
```


### vcf2struct.1.0.py (python3)


This program has been designed to perform statistics, filter, manipulate vcf files as well as analysing the mosaique structure of genomes.

#### Mandatory Options:

```
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
```

#### Options and outputs depend on analysis type:

* **RANDOM_SUB_SET:** Generate a vcf subset from the original vcf file in which variant line are sampled.<br/>

*Options:*

```
--vcf: A vcf file.
--nRand: Number of variant site to get randomly from the vcf file (integer). [Default: 1000]
--prefix: The prefix for output files. [Default: WorkOnVcf]
```

*Outputs:*<br/>
**_subset.vcf:** A vcf file corresponding to a random line subset of the first vcf<br/>
<br/>


* **STAT:** Calculate statistics on the vcf file.<br/>

*Options:*

```
--vcf: A vcf file.
--names: A one column file containing accession names to treat.
--prefix: The prefix for output files. [Default: WorkOnVcf]
--gff3: (optional) A gff3 file containing gene annotation. If not filled, the statistics will return 0 for these coding sequence values.
```

*Outputs:*<br/>
**_general.stat:** A file containing vcf global statistics such as number of indel and SNP sites, number of different tags in the FILTER columns, number of sites with 1,2,3,4, ... variants, number of transitions, transversions and variant in annotated regions.<br/>
**_accession.stat:** A file containing for each accessions missing data number, number of alleles specific to this accession in the vcf, number of homozygous sites identical to the reference, number of homozygous sites different from the reference and number of heterozygous sites.<br/>
<br/>


* **FILTER:** Filter a vcf file based on several parameters such as datapoint coverage and allele coverage, number of variant and variant type.<br/>

*Options:*

```
--vcf: A vcf file.
--names: A one column file containing accession names to treat.
--outgroup: (optional) A one column file containing accession names that will not be used for filtering but will remain in the output file.
--prefix: The prefix for output files. [Default: WorkOnVcf]
--MinCov: Minimal coverage by accession to keep genotype calling (integer). If the value is lower, genotype will be converted to unknown for the concerned accession. [Default: 10]
--MinAl: Minimal allele coverage by accession to keep genotype calling (integer). If the value is lower for at least one allele, genotype will be converted to unknown for the concerned accession. [Default: 3]
--nMiss: Maximal number of missing genotype in a line to keep the line (integer). [Default: 0]
--RmAlAlt: (optional) Number of alleles at the site in filtered accessions to remove the variant site (several values can be passed and should be sepatated by ":"). Values: 1,2,3,...,n
--RmType: (optional) Variant status to filter out (several values can be passed in this case they should be separated by ":"). 
	Possible values: 
		*Values which can be found in the FILTER column: PASS, DP_FILTER, QD_FILTER, SnpCluster, 
    	*Other values: INDELS, SNP, AUTAPO (accession specific variant site).
```

*Outputs:*<br/>
**_filt.vcf:** a filtered vcf file based on passed options.<br/>
<br/>


* **COMPARE:** Compare two variant accessions (from two different vcf files, in the same vcf file). This will output in standard output specific variant lines for accession1 and specific variant lines for accession2 as well as shared variant sites. In addition, identical variant calling will be counted and proportion among shared variant sites will be calculated.<br/>

*Options:*

```
--vcf: A vcf file.
--comp1: Accession name to compare. If 2 vcf files are provided and only one accession name is passed, this name will be searched in both vcf files.
--vcf2: (optional if --comp2 filled) A second vcf file.
--comp2:(optional if --vcf2 filled) Second accession to compare.
If 2 vcf and 2 names are passed, comp1 will be searched in vcf and comp2 will be searched in vcf2.
```
<br/>
<br/>


* **ADD_REF:** Add a haploid accession corresponding to the reference to the vcf called ref_silico_call.<br/>

*Options:*

```
--vcf: A vcf file.
--prefix: The prefix for output files. [Default: WorkOnVcf]
--ref_cov: Value to put to AD and DP flags (integer). [Default: 1000]
```

*Outputs:*<br/>
**_add_ref.vcf:** A new vcf file.<br/>
<br/>


* **AL_IDENTITY:** Calculate genotype identity.<br/>

*Options:*

```
--vcf: A vcf file.
--prefix: The prefix for output files. [Default: WorkOnVcf]
```

*Outputs:*<br/>
**_ident.mat:** A matrix containing genotype pairwise identity.<br/>
<br/>


* **FACTORIAL:** Perform factorial analysis on the vcf file. This analysis is performed in several steaps.<br/>

	1- First the vcf file is recoded as followed: For each allele at each variants site two markers were generated; One marker for the presence of the allele (0/1 coded) and one for the absence of the allele (0/1 coded). <br/> <img src="http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig1.png"  width="700"> <br/>Only alleles present or absent in **part** (not all) of selected accessions were included in the final matrix file.<br/> If groups information was passed to the script, alleles groups were attributed based on the following rule: the allele is attributed to a group if it is only present in this group but not in other defined groups. If no grouping information the GROUP column is filled with UN value. This grouping value **doesn't have any influence on the analysis**, it only allows to add colors graphs drawn. It will also help to validate if the structure of your data correspond to the one you suspect.<br/>
    2- The factorial analysis was performed on the transposed matrix using R. Graphical outputs of the analysis were draw and for example accessions and alleles can be projected along axis in the following picture. <br/> <img src="http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig2.png"  width="700"> <br/>In this example allele projected along synthetic axis were colorated if group informations were passed to the program.<br/>

*Options:*

```
--vcf: A vcf file.
--names: A one column file containing accession names to work with.
--prefix: The prefix for output files. [Default: WorkOnVcf]
--nAxes: Axis number to keep for the factorial analysis (integer). [Default: 4]
--mulType: Multivariate analysis type. Possible values: coa, pca and pca_normed [Default: coa]
--group: (optional) A file containing two sections: A section[group] with in col 1 accession name ; col 2 group (UN for unknown group). All group should be in capital letters. A section [color], that define for each group a color for pca drawing (in RGB+alpha percentage, ex: red=1:green=0:blue=0:alpha=0.1)
```

*Outputs:*<br/>
**_axis_x_vs_y_accessions.pdf:** Several pdf files showing accessions projected along x and y axis.<br/>
**_axis_x_vs_y.pdf:** Several pdf files showing accessions projected along x and y synthetic axis in a first graphe and allele projection along x and y synthetic axis.<br/>
**_inertia.pdf:** A pdf files showing axis inertia.<br/>
**_matrix_4_PCA.tab:** A tabulated file of the recoded vcf file passed to R for the factorial analysis.<br/>
**_multivariate.R:** The R script file passed to R to do the abalysis.<br/>
**_multivariate.Rout:** The R script log file.<br/>
**_individuals_coordinates.tab:** A tabulated file of individuals coordinates along synthetic axis.<br/>
**_variables_coordinates.tab:** A tabulated file of allele coordinates along synthetic axis.<br/>
<br/>


* **SNP_CLUST-Kmean:** Perform a k-mean clustering of allele based on their coordinates on synthetic axis. It uses the k-mean algorithm of scikit-learn.<br/>

*Options:*

```
--VarCoord: The *_variables_coordinates.tab file generated by FACTORIAL.
--dAxes: Axes to use in kmean clustering. Axis should be separated by ":".
--nGroup: Group number for the k-mean algorithm that will cluster variables (alleles) based on their coordinates. [Default: 2]
--mat: The *_matrix_4_PCA.tab tabulated file containing variant allele encoded generated by FACTORIAL.
--prefix: The prefix for output files. [Default: WorkOnVcf]
--thread: Number of processor available. [Default: 1]
--iter: Parallele k-mean clustering different startpoints performed. [Default: 100]
```

*Outputs:*<br/>
**_centroid_coordinates.tab:** Coordinates of the distincts final centroides calculated for the X k-mean random starpoints.<br/>
**_centroid_iteration_grouping.tab:** A tabulated file indicating centroid group attribution.<br/>
**_kMean_allele.tab:** A tabulated file equivalent to *_matrix_4_PCA.tab in which an column (named K-mean_GROUP) containing k-mean allele grouping has been added.<br/>
**_group_color.tab:** Color file in which a RGB color as been attributed to each k-mean cluster group. This file can be used for VISUALIZE_VAR_2D and VISUALIZE_VAR_3D tools.<br/>
**_kMean_gp_prop.tab:** A tabulated file reporting for each allele the probability to be in each groups. This is not a "real" probability, the idea was to have a statistics in case you want to filter alleles. This value was calculated as the inverse of the euclidian distance of one point and each centroids and these values were normalized so that the sum is equal to 1.<br/>
<br/>


* **SNP_CLUST-MeanShift:** Perform a clustering using the MeanShift algorithm of scikit-learn of allele based on their coordinates on synthetic axis.<br/>

*Options:*

```
--VarCoord: The *_variables_coordinates.tab file generated by FACTORIAL.
--dAxes: Axes to use in mean shift clustering. Axis should be separated by ":".
--quantile: The quantile value to estimate de bandwidth parameters used in the MeanShift. Value should be in [0:1]. [Default: 0.2]
--mat: The *_matrix_4_PCA.tab tabulated file containing variant allele encoded generated by FACTORIAL.
--prefix: The prefix for output files. [Default: WorkOnVcf]
--thread: Number of processor available. [Default: 1]
```

*Outputs:*<br/>
**_centroid_coordinates.tab:** Coordinates of centroides calculated.<br/>
**_centroid_iteration_grouping.tab:** A tabulated file indicating centroid group attribution.<br/>
**_kMean_allele.tab:** A tabulated file equivalent to *_matrix_4_PCA.tab in which an column (named K-mean_GROUP) containing k-mean allele grouping has been added.<br/>
**_group_color.tab:** Color file in which a RGB color as been attributed to each k-mean cluster group. This file can be used for VISUALIZE_VAR_2D and VISUALIZE_VAR_3D tools.<br/>
**_kMean_gp_prop.tab:** A tabulated file reporting for each allele the probability to be in each groups. This is not a "real" probability, the idea was to have a statistics in case you want to filter alleles. This value was calculated as the inverse of the euclidian distance of one point and each centroids and these values were normalized so that the sum is equal to 1.<br/>
<br/>


* **VISUALIZE_VAR_2D:** Perform plots of alleles projected along synthetic axes.<br/>

*Options:*

```
--VarCoord: The *_variables_coordinates.tab file generated by FACTORIAL.
--dAxes: Axes to plot. Axis should be separated by ":". All conbination between axis will be drawn.
--mat: The *_matrix_4_PCA.tab or *_kMean_allele.tab.
--group: (optional) A file containing at least the a section [color], that define for each group a color (in RGB+alpha percentage, ex: red=1:green=0:blue=0:alpha=0.1). This file can be the *_group_color.tab generated by SNP_CLUST.
--dGroup: Groups IDs to draw. Groups names should be separated by ":". Groups names will be searched in the file provided in --mat argument in the "K-mean_GROUP" column, and if such column is not present, they will be searched in "GROUP" columns. If none of these columns are found the program will exit without finishing.
--prefix: The prefix for output files. [Default: WorkOnVcf]
```

*Outputs:*<br/>
**_2d_axis***x***_vs_axis***y***.png:** Several png files showing grouped alleles projected along x and y synthetic axes.<br/>
<br/>


* **VISUALIZE_VAR_3D:** Perform an interactive 3d plots of alleles projected along 3 synthetic axes.<br/>

*Options:*

```
--VarCoord: The *_variables_coordinates.tab file generated by FACTORIAL.
--dAxes: 3 axes to plot. Axis should be separated by ":". All conbination between axis will be drawn.
--mat: The *_matrix_4_PCA.tab or *_kMean_allele.tab.
--dGroup: Groups IDs to draw. Groups names should be separated by ":". Groups names will be searched in the file provided in --mat argument in the "K-mean_GROUP" column, and if such column is not present, they will be searched in "GROUP" columns. If none of these columns are found the program will exit without finishing.
--group: (optional) A file containing at least the a section [color], that define for each group a color (in RGB+alpha percentage, ex: red=1:green=0:blue=0:alpha=0.1).
```


* **FILTER_ON_MAX_GP_PROP:** Filter the matrix file by removing ambiguous grouped alleles. *i.e.* alleles which changed of groups during the k-mean distinct attempts.<br/>

*Options:*

```
--gpPropFile: The group file proportion, *_kMean_gp_prop.tab file generated by SNP_CLUST.
--mat: The *_kMean_allele.tab file generated by SNP_CLUST.
--gpPropValue: Minimal value to keep the allele grouped. This value should be comprised between 0 and 1 [Default: 0.95]
--prefix: The prefix for output files. [Default: WorkOnVcf]
```

*Outputs:*<br/>
**_kMean_allele_filtered_with_** *x* **_value.tab:** A filtered file of *_kMean_allele.tab file in which ambiguous alleles were removed.<br/>
<br/>


* **MERGE_VCF:** Add an accession from a vcf to a second one vcf file. If the variant line is absent for the added accession, missing genotype if filled. If a variant line is present for the accession but absent in the vcf, this line will not be added.<br/>

*Options:*

```
--vcf: A vcf file in which you want to add an accession.
--vcf2: A vcf file containing accession you want to add to the fist vcf file.
--comp1: Accession name to add to the first vcf.
--prefix: The prefix for output files. [Default: WorkOnVcf]
```

*Outputs:*<br/>
**_merged.vcf:** A new vcf file with the accessions variant calling added.<br/>
<br/>


* **ALL_PROP:** Calculate for each accessions the number of grouped allele for each groups.<br/>

*Options:*

```
--vcf: A vcf file.
--mat: A matrix file containing allele grouping with either "K-mean_GROUP" or "GROUP" columns, if "K-mean_GROUP" does not exists the "GROUP" column will be used. These files are those generated by FACTORIAL (*_matrix_4_PCA.tab), SNP_CLUST (*_kMean_allele.tab) or  FILTER_ON_MAX_GP_PROP (*_kMean_allele_filtered_with_x_value.tab). 
--prefix: The prefix for output files. [Default: WorkOnVcf]
--names: (optional) A one column file containing accession names to calculate proportion on (it can be accessions that were not used for PCA analysis but originating from the same original vcf).
--dGroup: (optional) Groups IDs to calculate proportion on. Groups names should be separated by ":". Groups names will be searched in the file provided in --mat argument in the "K-mean_GROUP" column, and if such column is not present, they will be searched in "GROUP" columns. If none of these columns are found the program will exit without finishing.
--ExclChr: (optional) A list of chromosome to exclude from the analysis separated by ":".
```

*Outputs:*<br/>
**_gp_prop.tab:** A tabulated file containing for each accessions the number of grouped allele for each groups.<br/>
<br/>


* **GET_GENOTYPE:** Sometimes as a biologist you want to see the data and more precisely have a look at the genotypes (in term of A,T,G,C) along chromosomes. **This tool if for you!**<br/>

*Options:*

```
--vcf: A vcf file.
--names: A one column file containing accession names to plot genotypes
--prefix: The prefix for output files. [Default: WorkOnVcf]
```

*Outputs:*<br/>
**_genotype.gen:** A tabulated file containing, for each accessions (columns) and at each positions (rows) along chromosomes, the genotype.<br/>
<br/>


* **GET_GENOTYPE_AND_GROUP:** Sometimes as a biologist you want to see the data and more precisely have a look at the genotypes (in term of A,T,G,C)and allele grouping along chromosomes. **This tool if for you!**<br/>

*Options:*

```
--vcf: A vcf file.
--names: A one column file containing accession names to plot genotypes
--mat: A matrix file containing allele grouping with either "K-mean_GROUP" or "GROUP" columns, if "K-mean_GROUP" does not exists the "GROUP" column will be used. These files are those generated by FACTORIAL (*_matrix_4_PCA.tab), SNP_CLUST (*_kMean_allele.tab) or  FILTER_ON_MAX_GP_PROP (*_kMean_allele_filtered_with_x_value.tab).
--dGroup: Groups IDs to calculate proportion on. Groups names should be separated by ":". Groups names will be searched in the file provided in --mat argument in the "K-mean_GROUP" column, and if such column is not present, they will be searched in "GROUP" columns. If none of these columns are found the program will exit without finishing.
--prefix: The prefix for output files. [Default: WorkOnVcf]
```

*Outputs:*<br/>
**_genotype.gen:** A tabulated file containing, for each accessions (2 columns per accessions) and at each positions (rows) along chromosomes, the genotype (first of the two columns) and allele grouping (second of the two columns).<br/>
<br/>


### vcf2linear.1.0.py (python3)

This programs aims at performing a chromosome painting of accessions along chromosome based on the allele grouping correponding to the ancestral groups. The idea is to search, for distinct regions (along the chromosomes) of the accession we want to study, the probability of an ancestral origin (homozygous or heterozygous). Because, distinct ancestral groups can contain different specific alleles, show distinct heterozygosity levels and allele fixation (du to distinct reproductive mode for example) such probability is not straightforward. An idea can be to calculate for each regions, the expected ancestral allele grouping proportion (based on ancestral accessions) and compare this value to the observed one. Because dataset can be small and/or representatives of ancestral group are very few, the expected ancestral allele grouping proportion can be difficult to estimate. <br/>
This programs solves the problem by simulating ancestral population from accessions identified as representatives of the ancestral groups under panmixie hypothesis. For each ancestral populations a total of 100 individuals were simulated. The same was performed for hybrids between two populations. At the and of this simulation process mean and standard deviation of allele grouping proportions were calculated for each ancestral group on sliding windows of size n overlapping of n-1.<br/>
These grouping proportions were then calculated on tested accessions and estimated proportion were used to calculate the ancestry probability as followed.<br/>

With:<br/>
**_Hmu_**: Mean gX expected number for the group gX based on simulations on homozygous groups,<br/>
**_Hsd_**: Maximal standard deviation value calculated for all groups based on simulations on homozygous groups,<br/>
**_Emu_**: Mean gX expected number for the group gX based on simulations on heterozygous groups,<br/>
**_Esd_**: Maximal standard deviation value calculated for all groups based on simulations on heterozygous groups,<br/>
**_Nmu_**: Mean noisy allele grouped based on simulations (calculated as the number of alleles not from the right group according to the simulated population),<br/>
**_Nsd_**: Noise standard deviation,<br/>
**_obs_**: gX number observed for the accession on the window,<br/>
**_PgX_**: Probability to be at least gX,<br/>
**_PgN_**: Probability to be in the noise (unknown/additionnal ancestor),<br/>
**_N(x, mu, sd)_**: Probability density function of the normal distribution for mean = *mu* and standard deviation = *sd*.

* probability to be homozygous for group gX:<br/>
	* if *obs* **_superior or equal_** (*Hmu* - *Hsd*): *PgX* = 1<br/>
	* if *obs* **_inferior_** (*Hmu* - *Hsd*) : *PgX = N(obs, Hmu-Hsd, Hsd)/N(Hmu-Hsd, Hmu-Hsd, Hsd)*<br/>
* probability to be heterozygous for group gX:<br/>
	* if *obs* **_superior or equal_** (*Emu* - *Esd*): *PgX* = 1<br/>
	* if *obs* **_inferior_** (*Emu* - *Esd*) : *PgX = N(obs, Emu-Esd, Esd)/N(Emu-Esd, Emu-Esd, Esd)*<br/>
* probability to be in the noise:<br/>
	* if *obs* **_inferior or equal_** (*Nmu* + *Nsd*): *PgN* = 1<br/>
	* if *obs* **_superior_** (*Nmu* + *Nsd*): *PgN = N(obs, Emu+Esd, Esd)/N(Emu+Esd, Emu+Esd, Esd)*<br/>

Genotype grouping at the window was then attributed based on the maximal probability and haplotypes representation was performed trying to minimize the recombination events.


*Options:*

```
--vcf: A vcf file.
--mat: A matrix file containing allele grouping with either "K-mean_GROUP" or "GROUP" columns, if "K-mean_GROUP" does not exists the "GROUP" column will be used. These files are those generated by FACTORIAL (*_matrix_4_PCA.tab), SNP_CLUST (*_kMean_allele.tab) or  FILTER_ON_MAX_GP_PROP (*_kMean_allele_filtered_with_x_value.tab).
--names: A one column file containing accession names to treat (it is not necessarilly thoses used in the Factorial analysis but they must share variant sites).
--namesH: A two column file containing accession names used to simulate populations. 
--chr: Chromosomes names to work with (each chromosome names should be separated by ":"). If not filled, all chromosomes will be used.
--win: Half window size (allele number) around a variant site to evaluate the structure at the site (integer). [Default: 25]
--gcol: A file containing at least the a section [color], that define for each group a color (in RGB+alpha percentage, ex: red=1:green=0:blue=0:alpha=0.1). This file can be the *_group_color.tab generated by SNP_CLUST.
--prefix: The prefix for output files. [Default: WorkOnVcf]
```

*Output* <br/>

A folder with the name passed in --prefix options which contained several files for each accessions and each chromosomes:<br/>

**_Accession_chromosome.tab:** A tabulated file with for each window around a given position:

* 1- the count of each grouped alleles *(as much columns as ancestral groups)*, 
* 2- the expected number of grouped alleles, for a given group, in case of single ancestral origin (homozygous) *(as much column as ancestral groups)*,
* 3- the expected number of grouped alleles in case of single ancestral origin (homozygous) *(as much columns as ancestral groups)*,
* 4- the corresponding expected standard deviation *(as much column as ancestral groups)*,
* 5- the retained maximal standart deviation for only one ancestry,
* 6- the expected number of grouped alleles, for a given group, in case of two ancestral origin (heterozygous) *(as much columns as ancestral groups)*,
* 7- the expected number of grouped alleles in case of two ancestral origin (heterozygous) *(as much columns as ancestral groups)*,
* 8- the corresponding expected standard deviation *(as much columns as ancestral groups)*,
* 9- the retained maximal standart deviation for two distinct ancestry ,
* 10- the probability to have a single origin for the group X *(as much columns as ancestral groups)*,
* 11- the probability to have a group X origin and a second origin *(as much columns as ancestral groups)*,
* 12- the heterozygosity level (in the studied window),
* 13- expected noise count for each groups *(as much columns as ancestral groups)*,
* 14- expected noise standard deviation *(as much columns as ancestral groups)*,
* 15- the retained maximal noise standart deviation,
* 16- accession noise probabilities *(as much columns as ancestral groups)*.<br/>

**Accession_chromosome_haplo1.tab:** A tabulated file listing grouped region bloacs for haplotype1<br/>
**Accession_chromosome_haplo2.tab:** A tabulated file listing grouped region bloacs for haplotype2<br/>
**Accession_chromosome_density.pdf:** A pdf file summarysing all statistics (heterozygosity, expected and observed ancestries along chromosomes, accessions ancestry probabilities and infered accession haplotypes)<br/>
<br/>

### haplo2kar.1.0.py

This program perform a synthesis of the haplotypes reconstructed by vcf2linear.1.0.py by drawing for an accession the chromosomal painting for all its chromosomes.

*Options:*

```
--acc: Accession name.
--chr: Chromosomes list to draw (separated by ":")
--gcol: A file containing at least the a section [color], that define for each group a color (in RGB+alpha percentage, ex: red=1:green=0:blue=0:alpha=0.1). This file can be the *_group_color.tab generated by SNP_CLUST.
--dg: Groups to draw, separated by ":". Groups not passed here will be replaced by grey.
--centro: A tabulated file locating pericentromeric regions. In column 1: chromosome name, column 2: start, column 3: end.
```

*Output:* <br/>
**Accession.pdf:** A pdf file with chromosome painting for all chromosomes of this accession.<br/><br/>



### haplo2karByChr.1.0.py

This program perform a synthesis of the haplotypes reconstructed by vcf2linear.1.0.py by drawing for a chromosome the chromosomal painting for all selected accessions.

*Options:*

```
--acc: Accession names configuration file with Column1: accession name and column2: ploidy level.
--chr: Chromosome to draw
--gcol: A file containing at least the a section [color], that define for each group a color (in RGB+alpha percentage, ex: red=1:green=0:blue=0:alpha=0.1). This file can be the *_group_color.tab generated by SNP_CLUST.
--dg: Groups to draw, separated by ":". Groups not passed here will be replaced by grey.
--centro: A tabulated file locating pericentromeric regions. In column 1: chromosome name, column 2: start, column 3: end.
--prefix: Prefix for output files. [Default: All_acc]
--maxChr: Maximal haplotype number to draw in a pdf. [Default: 53]
```

*Output:* <br/>
**prefix_chromosme_X.pdf:** As much pdf file as necessary with chromosome painting for all chromosomes of this accession.<br/><br/>


### haplo2Circos.1.0.py

This program perform a synthesis of the haplotypes reconstructed by vcf2linear.1.0.py by drawing a circos representation of the chromosomal painting for all selected accessions and chromosomes.

*Options:*

```
--acc: Accession names configuration file with Column1: accession name and column2: ploidy level.
--chr: Chromosome(s) to draw, separated by ":"
--gcol: A file containing at least the a section [color], that define for each group a color (in RGB+alpha percentage, ex: red=1:green=0:blue=0:alpha=0.1). This file can be the *_group_color.tab generated by SNP_CLUST.
--dg: Groups to draw, separated by ":". Groups not passed here will be replaced by grey.
--centro: A tabulated file locating pericentromeric regions. In column 1: chromosome name, column 2: start, column 3: end.
--prefix: Prefix for output files. [Default: All_acc]
```

*Output:* <br/>
* **Circos\_All\_Admix.conf**: the configuration file used by circos. I choose to keep this file so that it can be edited to change aspect of the Figure if you want. After editing this file you will only have to run the following command line to have your new Figure : *circos -conf Circos_All_Admix.conf -noparanoid*
* **Circos\_All\_Admix\_housekeeping.conf**: a second file used by circos,
* **Circos\_All\_Admix.kar**: a third file also used by circos,
* **Circos\_All\_Admix.png**: the circos Figure.<br/><br/>
