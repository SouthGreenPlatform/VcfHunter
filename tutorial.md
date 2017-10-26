Tutorial for VcfHunter
======================

This tutorial go through all steps from variant calling from DNAseq
and/or RNAseq to caracterization of genome mosaic structure using
***VcfHunter tools***

Go to the VcfHunter folder (Scripts can be run from any folder but the
command lines in this tutorial assume you are in this folder and that
you have python 2.7 and 3 versions).

### Available data:

*data/reference/* is a folder containing a fasta file with the reference
sequence and a gff file containing the reference sequence annotation
*data/reads/* is a folder containing DNAseq paired reads simulated from
3 gene pools of diploid ancestral accessions (11, 10 and 9 accessions
from each gene pool respectively) and 10 DNAseq and 10 RNAseq paired
reads simulated from hybrids diploids accessions resulting from several
crosses of accessions from the ancestral pool. *data/config/* is a
folder containing two configuration files needed for the RNAseq and
DNAseq variant calling

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
depending on the type of data (RNA or DNA), two different pipeline
shoold be use for read mapping and post mapping processing.

### DNA seq variant calling

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

### RNA seq variant calling

a - Reference indexation: Create an index for STAR RNAseq mapping

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s a

1 folder is created in the current directory (name =
prefix+*ref*star\_1). This folder contains index constructed and used by
STAR to align RNAseq reads

b - Identification of splicing sites. During this step, all reads from
all libraries are merged and aligned against the reference sequence
which allowed STAR to identify splicing sites which will be used to
build a new index (step c) taking in account the splicing site.

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s b

A file (prefix+"*JUNC*ESTIMATION\_SJ.out.tab") containing the splicing
sites is generated

c - Reference indexation: with the splicing site information

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s c

1 folder is created in the current directory (name =
prefix+*ref*star\_2). This folder contains index constructed and used by
STAR to align RNAseq reads

d - Read mapping: Reads are aligned against the reference sequence using
STAR software

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s d

10 folders have been created containing each:

-   \*Aligned.out.sam files containing aligned reads (1 per libraries)
-   \*Log.final.out file containing mapping statistics of STAR (1 per
    libraries)
-   \*.SJ.out.tab file containing splicing sites (1 per libraries)

1 additional folder was created and named as passed to -p option
containing a file named prefix + "\_mapping.tab" which summarise mapping
statistics for all accessions.

e - Read merging and ordering: Reads from different libraries but the
same accessions are merged and sorted on coordinate.

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s e

Several new files have been added in each folders:

-   \*\_merged.bam files containing all aligned reads from all libraries
    of one accession.
-   \*\_merged.bai files (index of the bam file)

f - Removing duplicates reads: In this step, read duplicates are removed

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s e

2 new files have been added in each folders:

-   \*\_rmdup.bam file containing non redundant reads resulting from PCR
    and optical duplicates
-   \*\_duplicate file summarising read duplication statistics

An additional file was created in the folder named as passed to -p
option. This file named prefix + "*rmdup*stat.tab" contained duplicates
statistics for all accessions.

g - Bam reordering: In this step, is necessary for subsequent steps

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s g

2 new files have been added in each folders:

-   \*\_reorder.bam file containing reordered reads based on the order
    sequences are provided in the multifasta
-   \*\_reorder.bai file (index of the bam file)

h - Reads splitting: In this step, reads are splitted on splicing sites

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s h

2 new files have been added in each folders:

-   \*\_trim.bam file containing splitted reads
-   \*\_trim.bai file (index of the bam file)

i - Indel realigment: In this step, reads realigned around indels.
Countrary to **process\_reseq** a vcf on known indel could not be
provided. It will be implemented in futur...

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s i

2 new files have been added in each folders:

-   \*\_realigned.bam file containing realigned reads
-   \*\_realigned.bai file (index of the bam file)

j - Allele counting: In this step, reads are used to count for each
covered sites the number of reads supporting each bases
(A,T,G,C,N,\*=deletion).

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s j

Several new files have been added in each folders:

-   \*.gz files recording, for each chromosomes/sequences in the fasta
    provided as reference, the count of the number of reads supporting
    each bases (A,T,G,C,N,\*=deletion) at each covered site.

k - Vcf generation: generate the vcf for all accessions in the
configuration file

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s k

Several new files in the current directory:

-   prefix\_\*\_all\_allele\_count.vcf file for each
    chromosomes/sequences in the fasta provided as reference.

l - Vcf merging: generate a single vcf from all chromosome/sequences in
the fasta provided as reference. This step is not needed here has only
one sequence is passed but you can try the command line anyway ;-)

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s l

1 new file is created in the current directory:

-   prefix\_all\_allele\_count.vcf file

m - Exon coverage: Calculate for each gene filled in the gff file passed
in the configuration file and for each accessions, the proportion of the
exon covered by the library.

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p RNAseq -s m

1 new file is created in the current directory:

-   prefix\_all\_allele\_count.vcf file

***All the steps can be lauched in one command line. First remove the
generated files:***

    rm -rf sample71 sample72 sample73 sample74 sample75 sample76 sample77 sample78 sample79 sample80

Run the pipeline again

    python2 ../bin/process_RNAseq.1.0.py -c ../data/config/RNAseq.conf -t 8 -p DNAseq -s abcdefghijklm

***RNAseq and DNAseq bam can be used to generate a single vcf:*** At
this stage we assume that RNAseq and DNAseq have been processed until
step j and e for **process\_RNAseq** and **process\_reseq
respectively**. To generate a single vcf a new configuration file should
be created containing informations for RNAseq and DNAseq accession. This
file is available in *data/config/DNA\_RNAseq.conf*. To obtain the vcf
run the following command line:

    python2 ../bin/process_reseq_1.0.py -c ../data/config/DNA_RNAseq.conf -t 8 -p DNA_RNAseq -s fg

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

    python2 ../bin/VcfPreFilter.1.0.py -v DNA_RNAseq_all_allele_count.vcf -m 10 -M 10000 -f 0.05 -c 3 -o DNA_RNAseq_prefiltered.vcf

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

C - VCF statistics
------------------

Now that we have a vcf file in which we have only kept variant lines in
which we are relatively confident, it is time to perform some statistics
on the vcf. To do so, we first need to create a file in which we list
accessions to perform statistics on. This can be done with the following
command line:

    head -n 1000 DNA_RNAseq_prefiltered.vcf | grep "#CHROM" | sed 's/\t/\n/g' | tail -n 50 > all_names.tab

Statistics calculation can be done with the following command line:

    python3 ../bin/vcf2struct.1.0.py --vcf DNA_RNAseq_prefiltered.vcf --names all_names.tab --type STAT --prefix DNA_RNAseq_prefiltered

2 files are generated:

-   prefix+"\_general.stat" which regroup global statistics
-   prefix+"\_accession.stat" which regroup statistics per accessions
    such as missing data, accession specific alleles, homozygous sites
    and heterozygous sites.

When looking at the prefix+"\_general.stat" file we observe that there
is around 60000 sites in which we detected only one allele. This
suggests that either this is reference specific alleles or sequencing
errors or that it is monomorphous sites that have been identified
polymorphous based on the criteria passed in during vcf prefilter (but
binomial probabilities identified monorphous).

2 solution are available:

-   running again ***VcfPreFilter*** with more stringent parameters
-   filtering the vcf removing monoallelic variant. This is the solution
    we will choose in this tutorial.

In this same filtering step, we will also identify variant sites in
which we dont have enough confidence (--MinCov and MinAl options) and we
will remove also tri-allelic sites. To perform the filter, run the
following command line:

    python3 ../bin/vcf2struct.1.0.py --vcf DNA_RNAseq_prefiltered.vcf --names all_names.tab --type FILTER --MinCov 10 --MinAl 3 --nMiss 50 --RmAlAlt 1:3 --prefix DNA_RNAseq_prefiltered

The output is a filtered vcf file (prefix+"\_filt.vcf") in which all
mono- or tri- allelic lines are removed. In addition, all datapoints
which have a coverage lower than 10 reads or a minor allele frequency
lower than 3 reads are converted to missing data.

We will perform a final filter on the vcf file because the analysis that
perform the chromosome painting does not support the missing data.
Because the vcf is a mix of DNA and RNAseq, the accession genotyped with
the RNAseq will have more missing data. In this context we will filter
the DNAseq accession for no missing data. This is not clear why we do
that only on the DNAseq but it will become clearer and clearer along the
process why we do that.

So to filter the vcf on DNAseq data only we should generate 2 files. A
file which regroup accession names to use for filtering (DNAseq
accession) and accessions we do not want to use for filtering but which
we want to keep in the vcf anyway (RNAseq accessions). This can be done
with the following command line:

    grep Lib ../data/config/DNAseq.conf | cut -f 3 > DNAseq_names.tab
    grep Lib ../data/config/RNAseq.conf | cut -f 3 > RNAseq_names.tab

And to perform the filter, run the following command line:

    python3 ../bin/vcf2struct.1.0.py --vcf DNA_RNAseq_prefiltered_filt.vcf --names DNAseq_names.tab --outgroup RNAseq_names.tab --type FILTER --MinCov 10 --MinAl 3 --nMiss 0 --RmAlAlt 1:3 --prefix DNA_RNAseq_final

The output is a filtered vcf file (prefix+"\_filt.vcf") in which there
is no missing datas for the DNAseq accession but it remains missing data
for RNAseq data. You can calculate statistics if you want:

    python3 ../bin/vcf2struct.1.0.py --vcf DNA_RNAseq_final_filt.vcf --names all_names.tab --type STAT --prefix DNA_RNAseq_final_filt

D - PCA analysis
----------------

Now that we have a "good" vcf in which we are confident in the variant
line, it is time to analyse the dataset. In another way to say: to
perform the chromosome painting! The fisrt step of the chromosome
painting is to perform a COA analysis on the dataset to cluster the
alleles and the accession. Create a folder in which the analysis will be
performed and run the following command line:

    mkdir AllClust
    python3 ../bin/vcf2struct.1.0.py --vcf DNA_RNAseq_final_filt.vcf --names DNAseq_names.tab --type FACTORIAL --prefix AllClust/ClustAnalysis --nAxes 6 --mulType coa

The last command line run the factorial analysis (--type FACTORIAL
option). During this analysis the vcf file is recoded as followed: For
each allele at each variants site two markers were generated; One marker
for the presence of the allele (0/1 coded) and one for the absence of
the allele (0/1 coded).

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig3.png)

Only alleles present or absent in **part** (not all) of selected
accessions were included in the final matrix file named
***AllClust/ClustAnalysis\_matrix\_4\_PCA.tab*** in this example. An
additional column named "GROUP" can be identified. This column is filled
with "UN" value if no --group argument is passed. We will explain later
this arguement.

The factorial analysis (here a COA, --mulType option) was performed on
the transposed matrix using R (The R script is generated by the sript
and can be found here: ***AllClust/ClustAnalysis\_multivariate.R***). R
warning messages and command lines are recorded in the file named
***ClustAnalysis\_multivariate.Rout***. Graphical outputs of the
analysis were draw and for example accessions and alleles can be
projected along axis in the following picture.

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig4.png)

In this example the left graphe represent accessions projected along
axis 1 and 2 and the right represent the allele projected along
synthetic axis. A graphical representation is performed for each axis
combinations and each file is named according to the following
nomenclature ***prefix + \_axis\_X\_vs\_Y.pdf***. Several pdf for
accessions along axis only is also generated and are named according to
the following nomenclature ***prefix +
\_axis\_X\_vs\_Y\_accessions.pdf***.

The cumulated inertia is plotted in the graphe named
***AllClust/ClustAnalysis\_inertia.pdf***

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig5.png)

Individual and variables coordinates for the selected 6 fisrt axis
(--nAxes option) are recorded in files named
***AllClust/ClustAnalysis\_individuals\_coordinates.tab*** and
***AllClust/ClustAnalysis\_variables\_coordinates.tab*** respectively.
A third file named
***AllClust/ClustAnalysis\_variables\_coordinates\scaled.tab***
containing allele scaled coordinates (colums centered and reduced)
along synthetic axis is generated.

### The --group option

We assume that in some case you have additionnal informations on your
dataset such as which accessions are adimxed and which accessions are
likely to be the ancestral one. And maybe you want to verify/project
this information in your analysis. This can be done passing a
configuration file with two section to the --group option. This file can
be found in the data/config/ folder and is named ***AncestryInfo.tab***.
You can have a look at the file if you want but basically the two
sections are named [group] and [color] and contained respectively the
accession suspected grouping and a color (in RGB proportion) you want to
attribute to each group. Accessions with no group should filled with
"UN" value. **Warning: group name should be written in upper case (du to
R sorting)**.

**This additional information will not change the results of the
analysis.** It will only add color to the figures and fill the GROUP
column: in the GROUP column alleles groups were attributed based on the
following rule: the allele is attributed to a group if it is only
present in this group but not in other defined groups.

you can try this option running the following command line:

    mkdir AllClust_group
    python3 ../bin/vcf2struct.1.0.py --vcf DNA_RNAseq_final_filt.vcf --names DNAseq_names.tab --type FACTORIAL --prefix AllClust_group/ClustAnalysis --nAxes 6 --mulType coa --group ../data/config/AncestryInfo.tab

Outpout are strictely the same with two exceptions: pdf files are
colorated according to the groups passed in --group option. And the
"GROUP" column has been filled if the grouping rules are fulfilled.

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig6.png)

You can observe that if you compare this example with the preceding one
that only the orientation of the axis changed but the accessions and
allele coordiantes remained the same.

D - Allele clustering
---------------------

### Mean Shift clustering

Now that allele have been projected along synthetic axes, it is time to
clusterize these alleles. The idea is that the structure reflected by
the synthetic axis represent the ancestral structure. In this context,
the alleles at the extremities of the cloud of points will be the
ancestral ones. These alleles can be clusterized using several
approaches. In this turorial we will use a Mean Shift clustering
approach.

    python3 ../bin/vcf2struct.1.0.py --type SNP_CLUST-MeanShift --VarCoord AllClust/ClustAnalysis_variables_coordinates.tab --dAxes 1:2 --mat AllClust/ClustAnalysis_matrix_4_PCA.tab --thread 8 --prefix AllClust/ClustAnalysis

The Mean Shift clustering is performed with only the 2 first axes of the
COA (--dAxes 1:2) because the analysis showed that most of the inertia
is on these axes. With a mean shift approach, the number of group is
automatically detected.

During the process, several informations are returned to stdout, but at
the end of the process three main informations are returned: (i) the
number of alleles used for the analysis. Allele present or absent in all
accessions are removed. (ii) the number of estimated clusters which can
be found in the line:

    number of estimated clusters : 3

3.  the number of allele grouped within each group is returned and
    should look like as followed:

<!-- -->

    Group g0 contained 11439 dots
    Group g1 contained 2060 dots
    Group g2 contained 1657 dots

These results suggested a large desequilibrium in the grouping and we
will discuss latter the reason of such desequilibrium. We will focus at
the moment on the output and on the way we can interpret the analysis.
Five file are generated and can be found in the **AllClust** folder:

-   **ClustAnalysis\_kMean\_allele.tab** file which correspond to the
    ***ClustAnalysis\_matrix\_4\_PCA.tab*** in which the allele grouping
    has been recorded.
-   **ClustAnalysis\_centroid\_coordinates.tab** file which regroup the
    centroids coordinates.
-   **ClustAnalysis\_centroid\_iteration\_grouping.tab** file which
    records for each centroid its grouping.
-   **ClustAnalysis\_group\_color.tab** file that attribute a color to
    the groups.
-   **ClustAnalysis\_kMean\_gp\_prop.tab** file that report for each
    allele the probability to be in each groups. This is not a "real"
    probability, the idea was to have a statistics in case you want to
    filter alleles. This value was calculated as the inverse of the
    euclidian distance of one point and each centroids and these values
    were normalized so that the sum is equal to 1.

### Clustering visualization

The clustering visualization can be performed using 2 tools of
vcf2struct, a 2d visualization can be performed with --type
VISUALIZE\_VAR\_2D, and a 3d interactive visualization can be performed
using --type VISUALIZE\_VAR\_3D. The 2d of centroids visualization can
be launched with the following command line on centroids.

    python3 ../bin/vcf2struct.1.0.py --type VISUALIZE_VAR_2D --VarCoord AllClust/ClustAnalysis_centroid_coordinates.tab --dAxes 1:2 --mat AllClust/ClustAnalysis_centroid_iteration_grouping.tab --group AllClust/ClustAnalysis_group_color.tab --prefix AllClust/CentroidGrouping

Output is as much files as axis combinations named according to the
following nomenclature ***prefix + \_axisX\_vs\_axisY.png***. The output
should look like this:

![](http://banana-genome-http.cirad.fr/image/CentroidGrouping_axis1_vs_axis2.png)

This Figure represent the centroids location. Colors may not be the same
when you run the analysis because the color attribution is random.

Visualization of the allele grouping can be done as followed:

    python3 ../bin/vcf2struct.1.0.py --type VISUALIZE_VAR_2D --VarCoord AllClust/ClustAnalysis_variables_coordinates.tab --dAxes 1:2 --mat AllClust/ClustAnalysis_kMean_allele.tab --group AllClust/ClustAnalysis_group_color.tab --prefix AllClust/AlleleGrouping

The outpout named ***prefix + \_axis1\_vs\_axis2.png*** look like this
(color may be different):

![](http://banana-genome-http.cirad.fr/image/AlleleGrouping_axis1_vs_axis2.png)

In this picture you have the representation of the allele clustering
performed by the Mean Shift approach. We can clearly observe that the
red group (g0) is over-represented and may be split in two to have
alleles representing ancetral group g0 and en central cluster in which
allele not fixed/not representing ancestral groups are clustered.

As in this example we choose to work only with 2 axis, it is not
necessary to have a 3d visualization but we can try the command anyway:

    python3 ../bin/vcf2struct.1.0.py --type VISUALIZE_VAR_3D --VarCoord AllClust/ClustAnalysis_variables_coordinates.tab --dAxes 1:2:3 --mat AllClust/ClustAnalysis_kMean_allele.tab --group AllClust/ClustAnalysis_group_color.tab

A window which should look like this should open:

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig7.png)

This 3d visualisation can be rotated with the mouse.

### Important remark

In this example we observed that the allele of the g0 group are over
represented. And maybe this group should be splitted. In fact this group
should be splitted as there is 3 ancestral groups in the simulated data!
In real dataset you don't have this information but an admixture or SNMF
analysis can be carried first which should give you an idea of the group
number you expect. In addition the overrepresentation of a group (which
is not at the center of the datapoint) is a strong indicator of problem
in allele grouping! This problem in grouping may be due to the strong
contribution of this group to the hybrid accessions (cf accessions
projection along synthetic axis) and/or a structuration in hybrid
accessions that generate noise in wild accessions structure or problem
in calibration of the clustering (in this example here: estimation of
the bandwidth). This is a problem because during chromosome painting, g0
regions may be over represented and/or the chromosome painting may
represent a structure that is not related to the ancestral contribution
but rather a mix between ancestral contribution and hybrid structure!!!
In this context, it may be clever test several clustering parameters
until you reach the good number of cluster and even to re-run the COA
analysis and clustering clustering but only on accessions that are not
supposed hybrids to have a better allele grouping and ultimately a
better chromosome painting.

At this point, you may ask why not suggesting to do this directely as I
knew this will happen in our dataset. To this I respond that it is not
so long to run and this example provide a limit of the method explained
here that you should take in account! The idea is that you can try
several multivariate analysis and clustering methods and check which
method group the best your alleles.

We will do this in several steps: Firts we will vary the --quantile
parameter which will allow to change bandwidth estimation parameter used
in mean shift clustering. By default this value is put to 0.2. Lowering
it should allow to have more clusters. Try the following command line:

    python3 ../bin/vcf2struct.1.0.py --type SNP_CLUST-MeanShift --VarCoord AllClust/ClustAnalysis_variables_coordinates.tab --dAxes 1:2 --mat AllClust/ClustAnalysis_matrix_4_PCA.tab --thread 8 --prefix AllClust/ClustAnalysis --quantile 0.1

The output should look like this:

    loading modules
    modules loaded
    number of estimated clusters : 5
    Printing files
    Group g0 contained 9329 dots
    Group g1 contained 2124 dots
    Group g2 contained 1954 dots
    Group g3 contained 850 dots
    Group g4 contained 899 dots

At this point you have 5 clusters, while in theory you expected 4 (3
corresponding to ancestral groups and 1 corresponding to unassigned
alleles). This clustering can be visualized using the following command
line:

    python3 ../bin/vcf2struct.1.0.py --type VISUALIZE_VAR_2D --VarCoord AllClust/ClustAnalysis_variables_coordinates.tab --dAxes 1:2 --mat AllClust/ClustAnalysis_kMean_allele.tab --group AllClust/ClustAnalysis_group_color.tab --prefix AllClust/AlleleGrouping

![](http://banana-genome-http.cirad.fr/image/sub1_AlleleGrouping_axis1_vs_axis2.png)

You can observe that g3 and g4 are close and could be merged if we want
4 clusters. I don't say that it is the solution but for our example, we
will take this assertion. The idea is just to show you how to use the
program and the different command line. So, to reduce the cluster
number, we should increase the quantile value. We should then look a
value between 0.1 and 0.2. Let's take 0.15!

    python3 ../bin/vcf2struct.1.0.py --type SNP_CLUST-MeanShift --VarCoord AllClust/ClustAnalysis_variables_coordinates.tab --dAxes 1:2 --mat AllClust/ClustAnalysis_matrix_4_PCA.tab --thread 8 --prefix AllClust/ClustAnalysis --quantile 0.15

The output should look like this:

    loading modules
    modules loaded
    number of estimated clusters : 4
    Printing files
    Group g0 contained 9427 dots
    Group g1 contained 2196 dots
    Group g2 contained 1892 dots
    Group g3 contained 1641 dots

This time we have the expected group number. Data visualization should
be like this:

    python3 ../bin/vcf2struct.1.0.py --type VISUALIZE_VAR_2D --VarCoord AllClust/ClustAnalysis_variables_coordinates.tab --dAxes 1:2 --mat AllClust/ClustAnalysis_kMean_allele.tab --group AllClust/ClustAnalysis_group_color.tab --prefix AllClust/AlleleGrouping

![](http://banana-genome-http.cirad.fr/image/sub2_AlleleGrouping_axis1_vs_axis2.png)

We have an over represented central group g0, which correspond to
unassigned alleles. However, we can observe that the group g2 is steel a
little over-represented. This can be due to the presence of admixed
accessions with a strong contribution of this ancestral group, that
perturbate the multivatiate analysis. We can try to solves this problem
by running again the analysis but this time only on
homogeneous/ancestral accessions. This is what we will do in the next
point.

E - Runnig again Multivariate analysis and clustering
-----------------------------------------------------

First we need to create a new name file in which their will only be
"ancestral" accessions. Some accession may have introgression but they
should not be to much (this can be verified with several methods such as
the first COA, adixture or SNMF analysis). In this example, we the COA
analysis, we observed that acccession from sample61 to sample70 seemed
not "pure". I also know this from the simulations we have made to
generate the dataset ;-). These accessions should be removed from the
name file. To do this, run :

    head -n 31 DNAseq_names.tab > DNAseqFinalName.tab

Once this is done do again the Multivariate analysis:

    mkdir Final
    python3 ../bin/vcf2struct.1.0.py --vcf DNA_RNAseq_final_filt.vcf --names DNAseqFinalName.tab --type FACTORIAL --prefix Final/ClustAnalysis --nAxes 6 --mulType coa --group ../data/config/AncestryInfo.tab

When you look at the
***Final/ClustAnalysis\_axis\_1\_vs\_2\_accessions.pdf***, you can
observe that accessions are separated on the 2 axis with no intermediate
accessions that could perturbate the analysis.

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig8.png)

And the clustering:

    python3 ../bin/vcf2struct.1.0.py --type SNP_CLUST-MeanShift --VarCoord Final/ClustAnalysis_variables_coordinates.tab --dAxes 1:2 --mat Final/ClustAnalysis_matrix_4_PCA.tab --thread 8 --prefix Final/ClustAnalysis  --quantile 0.15

    loading modules
    modules loaded
    number of estimated clusters : 4
    Printing files
    Group g0 contained 9134 dots
    Group g1 contained 1837 dots
    Group g2 contained 1899 dots
    Group g3 contained 1957 dots

At this point you can observe that we have 4 clusters and clusterised
allele are more homogenous with one over-represented group (g0)
corresponding to unassigned alleles (in the center of the cloud of
point). You can observe centroids and allele clusering with the
previously described command lines:

    python3 ../bin/vcf2struct.1.0.py --type VISUALIZE_VAR_2D --VarCoord Final/ClustAnalysis_centroid_coordinates.tab --dAxes 1:2 --mat Final/ClustAnalysis_centroid_iteration_grouping.tab --group Final/ClustAnalysis_group_color.tab --prefix Final/CentroidGrouping
    python3 ../bin/vcf2struct.1.0.py --type VISUALIZE_VAR_2D --VarCoord Final/ClustAnalysis_variables_coordinates.tab --dAxes 1:2 --mat Final/ClustAnalysis_kMean_allele.tab --group Final/ClustAnalysis_group_color.tab --prefix Final/AlleleGrouping

The resulting centroid grouping:

![](http://banana-genome-http.cirad.fr/image/FinaleCentroidGrouping_axis1_vs_axis2.png)

The final allele grouping:

![](http://banana-genome-http.cirad.fr/image/FinaleAlleleGrouping_axis1_vs_axis2.png)

F - Performing the chromosome painting
--------------------------------------

Now that allele are clusterized, it is time to perform the chromosome
painting. The idea is to count, on a sliding window along chromosomes,
the number of grouped alleles of each groups for each hybrid accession
and then to attribute a chromosome region origin based on the maximal
grouped alleles. However it is not so simple, determining the expected
grouped allele number is not necessary an easy task: Sometime the
representant of the ancestral groups are very few, or the diversity of
these groups is not well represented or the degree of ancestral allele
fixation is not the same between accession depending if the accessions
are allogamous or autogamous. This program has been designed to solve
these problems (at least try) as explained in the following Figure:

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig9.png)

This figure explained the impact of an unbalance in the allele fixation
levels of ancestral population on the expectation of grouped allele
number at a position. We can observe that depending on the allele
fixation level in the population, for the same number of allele grouped,
the expected allele number of Gp2 is lower than Gp1 because all alleles
of Gp2 are not fixed! In the case you have a large sample of your
ancestral populations, estimating the expected grouped allele number on
the window is relatively easy, but this is not always the case. In this
context the program try to recreate the ancestral population from the
few representatives provided. This can be done by simulating ancestral
(our programs creates populations or 100 individuals) based on alleles
sampling from ancestral population representatives. Then in each
populations, the mean and standard deviation number of grouped alleles
is calculated on given window size which give the expectation of grouped
alleles on these windows to attribute 1, 2, 3, ... haplotypes. This
simulation step can be called with the *-T Simul* option in the program.
Because this simulation step is long and because sampling alleles at a
site followed a binomial test, the mean and sd expected groups at a
site, can be directely calculated from grouped alleles at the site. Then
the expected number of grouped alleles on a window is the sum of the
mean value obtained at a site. However, the standard deviation on the
window is not easy to calculate and then this value is approximated by
calcutating the square root of the sum of the variances at each sites of
the windows. This algorithm as also been implemented in our program and
can be called with *-T Binom* option. Because the estimation of standard
deviation is not completly satisfactory we also added a *--prop* option
which allowed to define an arbitrary *standard deviation* value
calculated as --prop \* Mean\_Expected\_value. All these calculations or
simulations are also done for hybrids population, to estimate the number
of grouped alleles in case of several ancestral origin at a site
(example: one haploype of Gp1 and one haplotype of Gp2).

An additional complexity level is that depending on the sampling of the
ancestral representative of the ancestral groups, the variance of the
expected number of alleles of a group can be greatly underestimated as
well as the mean value (cf Figure below). In this context, the --prop
value may be more appropriate.

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig10.png)

This figure showed the impact of the ancestral accession sampling on the
simulating of ancestral population. To partially solve this problem, we
chose to attribute the maximal variance found in a group and on a
window, to all groups in the same window.

After this short introduction, it is time to try the command line! First
you need to create a file which identify which accession belong to which
ancestral group. This can be donne by comparing the clustering and
accession grouping. In the following picture we showed the clustering of
alleles and accession projections allong synthetic axis.

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig10bis.png)

In this picture we can observe that sample1 to sample20 are accession
corresponding to the group **g1**, accessions from sample21 to sample40
belonged to group **g2** and accessions from sample41 to sample 49
corresponded to group **g3**. You can create this file manually, by
because I am lazy and I already have a preformated file, you can simply
run the following command line:

    grep sample ../data/config/AncestryInfo.tab | sed 's/X1/g1/' | sed 's/Y1/g2/' | sed 's/Z1/g3/' | grep -v 'UN' > ancestor.gp

**Be carreful, as the group attribution is random, you will have to
check to which group (X1, Y1 and Z1) correspond the allele group (g0,
g1, g2 and g3) to adjust the command line**

Once this file is ready, it is time to run the chromosome painting. In
this example, we will use the fast algorithm (*-T Binom* option). This
can be done with the \***vcf2linear.1.1.py** program by running the
following command line:

    mkdir Painting
    python3 ../bin/vcf2linear.1.1.py --vcf DNA_RNAseq_final_filt.vcf --names DNAseq_names.tab --namesH ancestor.gp --win 50 --mat Final/ClustAnalysis_kMean_allele.tab --prefix Painting --gcol Final/ClustAnalysis_group_color.tab --chr RefSeq --ploidy 2 -t 8 -T Binom

Several output are generated in the ***Painting*** folder. For each
accessions name in the ***DNAseq\_names.tab*** file, a file named
according to the following nomenclature
***prefix\_accessionName\_chr.tab*** is generated. This file record, for
each allele at the center of the window of 101 (--win 50 means a window
of 50 alleles before and after the considered allele), the number of
counted grouped alleles in the accession for each groups (g1, g2, g3),
the expected grouped alleles in case of only one haplotype (**H1**) of
the considered group is expected (Loc-mu-H1-g3, Loc-mu-H1-g1,
Loc-mu-H1-g2), the corresponding standard deviation value (Loc-sd-H1-g3,
Loc-sd-H1-g1, Loc-sd-H1-g2) and the maximal standard deviation retained
(Loc-max-sd-H1). Depending on the ploidy level of the studied accessions
(--ploidy) these values are calculated for 2 haplotypes of a same group
(**H2**), 3 haplotypes (**H3**), 4 haplotypes (**H4**), ... These values
are followed by a column named *hetero*, which calculated the
heterozygosity level in the window (proportion of heterozygous sites).
Then, probalility to have at least 1, 2, 3, ... (depending on the
--ploidy option) haplotypes for each group is rapported (Prob-H1-g3,
Prob-H1-g1, Prob-H1-g2) for at least one haplotype and (Prob-H2-g3,
Prob-H2-g1, Prob-H2-g2) for at least 2 haplotypes. The expected counted
alleles of a group resulting from noise is calculated in all simulated
accessions (in case of simulation) that have no contributors of these
groups or based on the sum of the binomial expected values calculated as
the number of time an allele of a group appear in individuals of other
groups. The mean values are reported in columns named :
*Loc-mu-noise-g3*, *Loc-mu-noise-g1*, *Loc-mu-noise-g2* and the
coresponding standard deviation can be found in columns named :
*Loc-sd-noise-g3*, *Loc-sd-noise-g1*, *Loc-sd-noise-g2*. The maximal
variance is also reported in column *Loc-max-sd-noise*. These values
allowed to calculate the probability that the counted alleles of a group
is from noise and these values are reported in column named:
*Prob-noise-g3*, *Prob-noise-g1*, *Prob-noise-g2*.

Haplotype probabilities are calculated as followed:

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig11.png)

If the observed grouped allele number is higher than the **expected
value - the maximal standard deviation** a probability of 1 is
attributed. This probability is the probability to have at least X
haplotypes in the analysed region. If the observed value is lower the
probability is calculated based on a density probability of a normal
distribution of mean value = **expected value - the maximal standard
deviation** and sd = **the maximal standard deviation**. This density
probability is normalized to reach a probability of 1 when the observed
value is equal to **expected value - the maximal standard deviation**.
The noise probability is calculated following the same phylosophie but
this time a probability of 1 is attributed if the observed grouped
allele number is lower than the **expected value + the maximal standard
deviation**. Else, the probability to be in the noise is calculated
based on a density probability of a normal distribution of mean value =
**expected value + the maximal standard deviation** and sd = **the
maximal standard deviation**. This density probability is normalized to
reach a probability of 1 when the observed value is equal to **expected
value + the maximal standard deviation**. Based on these probabilities,
an ancestral origin is attributed to each haplotypes of the studied
accession. The attribution is performed as followed, for each ancestral
group if the probability to be in the noise is higher than the
probability to be in the group, the probability to be in the group is
converted to 0. Then, the haplotype probability is calculated based on
the sorting of the maximal probability. If no consensus can be found
(*ex.* incompatible three best probabilities for a diploid) unknown
haplotypes are filled. Haplotypes are ordered trying to minimize
recombination events. Each haplotypes for or each accessions are
outputed in files found in the folder passed in --prefix option and are
named as followed: ***accessionName\_chr\_haploX.tab***. There are as
much file as there are haplotypes. In each files there are 5 columns:

-   1 - acessions name
-   2 - chromosome name
-   3 - start position of the bloc
-   4 - end position of the bloc
-   5 - ancestral origin of the bloc

These files will be used to draw figures in the following steps.

All these statistics are also plotted in files named according to the
following nomenclature ***accessionName\_chr\_density.pdf*** and can be
found in the folder passed in --prefix option. The following Figure is
obtained from the accession *sample68*

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig12.png)

This Figure regroup several informations such as the heterozygosity
level along the chromosome, for each ancestry the expected number of
grouped allele in case of 1 haplotype of this origin (green band), 2
haplotypes of this origin (red band) or the noise (grey band) and the
observed number of grouped allele in the accession (black curve). The
last graph represent the calculated probabilities and based on these
probabilities a chromosome painting is performed. *As the color
attribution is random, the colors of the chromosome painting may not be
the same when your run the program but the blocs should be the same.*

G - Chromosome painting visualization
-------------------------------------

The chromsome painting can be vizualized from distinct ways. For example
you can have several chromosomes from one accession and you may want to
plot all your chromosomes on a single file. This can be done with the
following command line:

    cd Painting
    python3 ../../bin/haplo2kar.1.0.py --acc sample68 --chr RefSeq --gcol ../Final/ClustAnalysis_group_color.tab --dg g1:g2:g3 --centro ../../data/reference/centro_pos.tab --ploidy 2

This command line outpout a pdf named ***sample68.pdf*** which should
look like this:

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig13.png)

This is not very different from the previous Figure but the idea was
just to try the command line. You can also observe that pericentromeric
regions (located by the ../../data/reference/centro\_pos.tab file) are
located in the figure with the grey lines while chromosome arms are
identified with black lines.

Just to have an idea of what it looks like when you have several
chromosomes you can do the following. Create new files from other
accessions that will simulate new chromosomes:

    sed 's/RefSeq/RefSeq1/' sample7_RefSeq_haplo1.tab > sample68_RefSeq1_haplo1.tab
    sed 's/RefSeq/RefSeq1/' sample7_RefSeq_haplo2.tab > sample68_RefSeq1_haplo2.tab
    sed 's/RefSeq/RefSeq2/' sample66_RefSeq_haplo1.tab > sample68_RefSeq2_haplo1.tab
    sed 's/RefSeq/RefSeq2/' sample66_RefSeq_haplo2.tab > sample68_RefSeq2_haplo2.tab

At this point we have generated 2 additionnals chromosomes (RefSeq1 and
RefSeq2) for sample68 which are in fact the chromosome painting of
sample7 and sample66. To do the chromosoime painting of this "false"
accession run the following command line:

    python3 ../../bin/haplo2kar.1.0.py --acc sample68 --chr RefSeq:RefSeq1:RefSeq2 --gcol ../Final/ClustAnalysis_group_color.tab --dg g1:g2:g3 --centro ../../data/reference/centro_pos.tab --ploidy 2

This command line overwrite the original Figure and the outpout should
look like this:

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig14.png)

You can observe our three chromosomes (RefSeq, RefSeq1 and RefSeq2) of
the same accession in one Figure.

You may also want compare the chromosome painting of several accessions.
*i.e.* to plot the chromosome painting of several accessions but for one
chromosome in one figure. For example you want to see the chromosome
painting of all admixed accessions.

Fisrt create the liste of admixed accession with in column 2 the ploidy
level:

    tail -n 10 ../DNAseq_names.tab | sed 's/$/\t2/' > ../DNAseq_names_admix.tab

Then create the Figure:

    python3 ../../bin/haplo2karByChr.1.0.py --acc ../DNAseq_names_admix.tab --chr RefSeq --gcol ../Final/ClustAnalysis_group_color.tab --dg g1:g2:g3 --centro ../../data/reference/centro_pos.tab --prefix All_Admix 

This command creates a file named ***All\_Admix\_RefSeq\_1.pdf*** which
contained haplotypes for accessions passed in --acc option. The outpout
should look like this (at the bottom of the page):

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig15.png)

This program is designed to print at most 53 (default parameters)
chromosome per pages, if you have more chromosomes to draw, if this
parameters is not changed the first 53 haplotypes will be drawn in a
fisrt pdf called ***All\_Admix\_RefSeq\_1.pdf*** and the remaining will
be drawn in a pdf called ***All\_Admix\_RefSeq\_2.pdf***. You can also
change the --maxChr parameters to draw more chromosomes per pdf.

And finally, you may also want to generate a circos representation of
your data with several accessions and chromosomes in the same Figure.
This can be done with the haplo2Circos.1.0.py. This can be done with the
following command line:

    python3 ../../bin/haplo2Circos.1.0.py --acc ../DNAseq_names_admix.tab --chr RefSeq --gcol ../Final/ClustAnalysis_group_color.tab --dg g1:g2:g3 --centro ../../data/reference/centro_pos.tab --prefix Circos_All_Admix

This programs outpouts 4 files:

-   **Circos\_All\_Admix.conf**: the configuration file used by circos.
    I choose to keep this file so that it can be edited to change aspect
    of the Figure if you want. After editing this file you will only
    have to run the following command line to have your new Figure :
    *circos -conf Circos\_All\_Admix.conf -noparanoid*
-   **Circos\_All\_Admix\_housekeeping.conf**: a second file used by
    circos,
-   **Circos\_All\_Admix.kar**: a third file also used by circos,
-   **Circos\_All\_Admix.png**: the circos Figure

The Figure obtained should look like this:

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig16.png)

In the circos picture, accessions are ordered from the outside to the
inside in following the order passed in --acc option.

Again in this example, because we do not have more than one chromosome,
we will do the same thing than previously for sample68: we will generate
2 additionnal chromosomes for 4 other accessions to have an idea of what
the circos looks like whit several chromosomes.

    sed 's/RefSeq/RefSeq1/' sample21_RefSeq_haplo1.tab > sample67_RefSeq1_haplo1.tab
    sed 's/RefSeq/RefSeq1/' sample21_RefSeq_haplo2.tab > sample67_RefSeq1_haplo2.tab
    sed 's/RefSeq/RefSeq2/' sample65_RefSeq_haplo1.tab > sample67_RefSeq2_haplo1.tab
    sed 's/RefSeq/RefSeq2/' sample65_RefSeq_haplo2.tab > sample67_RefSeq2_haplo2.tab
    sed 's/RefSeq/RefSeq1/' sample40_RefSeq_haplo1.tab > sample69_RefSeq1_haplo1.tab
    sed 's/RefSeq/RefSeq1/' sample40_RefSeq_haplo2.tab > sample69_RefSeq1_haplo2.tab
    sed 's/RefSeq/RefSeq2/' sample64_RefSeq_haplo1.tab > sample69_RefSeq2_haplo1.tab
    sed 's/RefSeq/RefSeq2/' sample64_RefSeq_haplo2.tab > sample69_RefSeq2_haplo2.tab
    sed 's/RefSeq/RefSeq1/' sample5_RefSeq_haplo1.tab > sample61_RefSeq1_haplo1.tab
    sed 's/RefSeq/RefSeq1/' sample5_RefSeq_haplo2.tab > sample61_RefSeq1_haplo2.tab
    sed 's/RefSeq/RefSeq2/' sample63_RefSeq_haplo1.tab > sample61_RefSeq2_haplo1.tab
    sed 's/RefSeq/RefSeq2/' sample63_RefSeq_haplo2.tab > sample61_RefSeq2_haplo2.tab
    sed 's/RefSeq/RefSeq1/' sample45_RefSeq_haplo1.tab > sample70_RefSeq1_haplo1.tab
    sed 's/RefSeq/RefSeq1/' sample45_RefSeq_haplo2.tab > sample70_RefSeq1_haplo2.tab
    sed 's/RefSeq/RefSeq2/' sample62_RefSeq_haplo1.tab > sample70_RefSeq2_haplo1.tab
    sed 's/RefSeq/RefSeq2/' sample62_RefSeq_haplo2.tab > sample70_RefSeq2_haplo2.tab

We also need to generate the configuration file to pass to the --acc
option:

    echo "sample67 2" > for_circos.tab
    echo "sample68 2" >> for_circos.tab
    echo "sample69 2" >> for_circos.tab
    echo "sample61 2" >> for_circos.tab
    echo "sample70 2" >> for_circos.tab

And then, it is time to draw the Circos:

    python3 ../../bin/haplo2Circos.1.0.py --acc for_circos.tab --chr RefSeq:RefSeq1:RefSeq2 --gcol ../Final/ClustAnalysis_group_color.tab --dg g1:g2:g3 --centro ../../data/reference/centro_pos.tab --prefix Circos_Acc

The output should look like this:

![](http://banana-genome-http.cirad.fr/image/Vcf2struct_Fig17.png)

H - Miscellaneous
-----------------

In this tutorial I have described main tools starting from fastq and
finishing with chromosome painting. However, the the Vcf2struct.1.0.py
allowed to perform several other tasks which are described in the
***README*** file. You can try them if you want.

In addition, once allele grouping as been performed, the analysis can be
also perfomed on the RNAseq data, even if they have not been used in the
analysis! Go to the TestTools folder and run the following command line:

    python3 ../bin/vcf2linear.1.1.py --vcf DNA_RNAseq_final_filt.vcf --names RNAseq_names.tab --namesH ancestor.gp --win 100 --mat Final/ClustAnalysis_kMean_allele.tab --prefix Painting --gcol Final/ClustAnalysis_group_color.tab --chr RefSeq --ploidy 2

The ouput are the same as for the DNAseq data. The only difference is
that the analysis will only be performed on site of grouped alleles also
present in the RNAseq dataset.

