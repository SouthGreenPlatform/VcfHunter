Tutorial for VcfHunter chromosome painting
==========================================

This tutorial aimed at studying the ancestral contribution along chromosomes
of one accession or a group of accessions. This program will work if a set
of ancestral accessions is available. These accessions should have no 
introgressions.

**Go to the VcfHunter/TestTools/ folder** (Scripts can be run from any
folder but the command lines in this tutorial assume you are in this
folder and that you have python 3 version).


Chromosome painting using non admixed ancestral accessions
----------------------------------------------------------------

### Available data:

*../data/config/Origin.tab* is a file which contained two column: a first
column containing ancestral accession names and a second column containing
their ancestral origin (this program can work until 8 distinct origins).
*../data/config/Vcf.conf* is a file which contained path to vcf files which
will be used for e-chromosome painting.
*../data/vcf/* is a folder containing the vcf for 5 chromosomes on 15 accessions
which will be used in this tutorial.

### Principle:
The programs described in this section performed two types of analysis: it
performs a chromosome painting along accessions according to ancestral
accessions defined and it also perform a plot of read coverage along the
chromosome to identify aneuploidy.

These programs worked as follows: based on the vcf file provided and a file providing the names of ancestral
accessions, the origin of each allele is attributed to an ancestral group according to the
following rule (two distinct rules are available):

rule1 - An allele is attributed to a group if it has been found only in this group.
or
rule2 - An allele is attributed to a group if it has been found in all members of
the group and absent from members of others groups.

Hence, for each allele of each group, and for the studied accession, the number of read
having this allele is calculated and divided by the number of reads at the position
in the studied accession. In this context, a diploid accession having two
chromosomes of the same origin (ex. red) should have a red allele ratio near 1
(see Figure below (A)). A diploid accessions with chromosomes of two distinct
origin (ex. one red and one green) should have green allele ratio near 0.5 and
red allele ratio near 0.5 too (see Figure below (B)). A same approach can be applied
to polyploid accessions: with a ratio of 0.33 and 0.66 for one and two ancestral
chromosomes of the same origin respectively for triploid (see Figure below (C)) and
a ratio of 0.25, 0.5, 0.75 for one, two and three ancestral chromosomes of the same
origin respectively for tetraploid accession (see Figure below (D)).

![](http://banana-genome-http.cirad.fr/image/AllelePropAndCov_Fig1.png)

In addition, for each position where an allele can be attributed to a group, the read
coverage for the accession is calculated. An average coverage is then calculated for
the accession and read coverage along the chromosome of the accession is then plotted
relative to the average coverage of the accession. This approach allow to identify missing
or supernumerary chromosomes (or chromosomal region) (See Figure below). 

![](http://banana-genome-http.cirad.fr/image/AllelePropAndCov_Fig2.png)

### Running analysis
Two programs are available and can be run depending on what figures you expect. The first
program allowed to perform a chromosome painting for all chromosome of one accessions.

Go to the TestTools folder and run the following command line:
```
../bin/vcf2allPropAndCov.py --conf ../data/config/Vcf.conf --origin ../data/config/Origin.tab --acc Kunnan --ploidy 2 --NoMiss n --all y
```

This command line analyse the "Kunnan" accession (--acc Kunnan) on chromosome 1,2,3,4
and 9 found in the vcf passed as a list to --conf ../data/config/Vcf.conf option.
Alleles were grouped according to accessions origin file passed to
--origin ../data/config/Origin.tab option. To attribute allele origin, ancestral
accession can have missing data (--NoMiss n) but all accessions without missing
data should have the allele (--all y).
**Remark:**  For any options passed, for a SNP position, if an ancestral group
as all its accessions with missing data, the SNP position is not used by the program.


This programs outpouts 4 files:

-   **Kunnan\_AlleleOriginAndRatio.tab** is a file describing for each grouped
    allele, its origin and the proportion of reads having this allele at the
    studied position in the accession.
-   **Kunnan\_stats.tab** is a file reporting statistics on SNP sites used,
    sites where an allele is attributed to each groups and alleles number
    attributed to each groups in the accession.
-   **KunnanCov.png** is a figure showing read coverage along chromosomes (see
    figure below for interpretation).
-   **KunnanRatio.png**: is a figure showing grouped allele read ratio along
    chromosomes (see figure below for interpretation).

The following figure describe the two outputs of the program and how to interpret
these outputs. In this tutorial, the two ancestral groups are named "AA" and "BB"
and as colors are arbitrarily attributed based on an alphanumeric sorting of
ancestor names and first color is green and second is red, then "AA" is green and
"BB" is red.

![](http://banana-genome-http.cirad.fr/image/AllelePropAndCov_Fig3.png)

Their is 3 additional admixed accessions, two triploids (GP1 and P025) and a
tetraploid one P1. These accessions can also be analysed  as Kunnan with the
following command lines:

```
../bin/vcf2allPropAndCov.py --conf ../data/config/Vcf.conf --origin ../data/config/Origin.tab --acc GP1 --ploidy 3 --NoMiss n --all y
../bin/vcf2allPropAndCov.py --conf ../data/config/Vcf.conf --origin ../data/config/Origin.tab --acc P025 --ploidy 3 --NoMiss n --all y
../bin/vcf2allPropAndCov.py --conf ../data/config/Vcf.conf --origin ../data/config/Origin.tab --acc P1 --ploidy 4 --NoMiss n --all y
```

The following figure is the picture of SNP read coverage you should obtain
and the interpretation that can be made from this coverage. We can observe
that for accession P025 a chromosome region is missing for the start of
chromosome 1 and the major part of chromosome 3.

![](http://banana-genome-http.cirad.fr/image/AllelePropAndCov_Fig4.png)

The following figure is the picture of allele ratio you should obtain and
green and red bars represent the interpretation that can be made from this
allele ratio and read coverage.

![](http://banana-genome-http.cirad.fr/image/AllelePropAndCov_Fig5.png)


One also want to perform this same analysis but by comparing one chromosome
from several accessions. This can be performed with the following command
line:

```
../bin/vcf2allPropAndCovByChr.py --conf ../data/config/Vcf.conf --origin ../data/config/Origin.tab --ploidy 3 --NoMiss n --all y --acc T04,T02,Kunnan,GP1,P025,P1
```

This command line analyse T04, T02, Kunnan, GP1, P025 and P1 accessions
(--acc T04,T02,Kunnan,GP1,P025,P1) on chromosome 1,2,3,4 and 9 found in the
vcf passed as a list to --conf ../data/config/Vcf.conf option. Alleles
were grouped according to accessions origin file passed to --origin
../data/config/Origin.tab option. To attribute allele origin, ancestral
accession can have missing data (--NoMiss n) but all accessions without
missing data should have the allele (--all y).
**Remark:** If the --acc option is omitted, the picture will be drawn for
all accessions found in the vcf. No more than 15 chromosomes are drawn per
figures. If there are more than 15 chromosomes, they are drawn in a second
file, a thrid, etc...

Two types of file are generated for each chromosomes. The coverage file
and the allele ratio file.
Here is an exemple of file generated for chromosome 3 on the 4 tested
accessions:

![](http://banana-genome-http.cirad.fr/image/AllelePropAndCov_Fig6.png)

In this figure, chromosome 3 of all tested accessions are represented in
the same figure. We added T04 and T02 which are accessions representative
of ancestors green and red, respectively. These accessions are not required
but I think that it might be important to show you that we can also analyze
accessions which are used as ancestor. It is also important to verify that
their is no introgression in accessions used as ancestor. by doing this, this
can be verified because according to our rules of ancestor specific allele
attribution, an introgression of a red region in a green ancestor would result
in the absence of red allele specific attribution in the introgressed region.
And this can be verified studying the allele ratio in ancestor.
One can also observe that an arbitrary ploidy level as been passed. This is
due to the fact that this program cannot manage distinct ploidy level at the
same time. However, this is not a problem, because this ploidy level is only
used to estimate the expected coverage in case of supernumerary or missing
regions in coverage analysis (horizontal lines in coverage figures). Thus
these line will not be exact if the ploidy level is not the good one but the
identification of supernumerary or missing regions will not be perturbed.
