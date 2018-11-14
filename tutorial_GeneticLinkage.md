Tutorial for VcfHunter Genetic linkage analysis
===============================================

This tutorial aimed at studying the genetic linkage of genotyped
individuals from a mapping population. Tools presented in this
tutorial goes from markers selections from vcf, marker coding,
linkage calculation and marker linkage representations along
chromosome/genetic map.

Go to the VcfHunter folder (Scripts can be run from any folder but the
command lines in this tutorial assume you are in this folder and that
you have python 3 version).


Available data:
---------------
*data/vcf/* is a folder containing several vcf including Carto.vcf.gz
that contained genotyping on 188 individuals from a triploid mapping
population (yes for the fun we will not use standard mapping pop)which
will be used in this tutorial.

A- Selecting segregating marker:
--------------------------------

### Context/Principle:
At the start of this section, we have 186 triploid genomes issued from
a cross between a diploid accession an a tetraploid accession. This
progeny (named from P001 to P204) and their parents (P1 and P2) have
been sequenced. Sequencing data were aligned against reference sequence
comprising 3 chromosomes and a variant calling was performed. The aim is
to use the variant calling file (vcf) generated to select segregating
markers.
Why selecting segregating marker and not using directly? Because, in the
real world their are missing data, not all polymorphous sites are true
polymorphism (i.e there are sequencing errors or multiple loci aligned at
the same position due to repeat sequences generated false polymorphism).
We thus, need to select the good markers and convert genotype matrix in
a way that programs that performed genetic map can understand.
To do so, two programs have been developed (***vcf2pop.1.0.py*** and
***vcf2popNew.1.0.py***) that do nearly the same thing but not exactly:

***vcf2pop.1.0.py*** has been designed to work specifically on diploid
data and as it is simpler to use, I will not describe this program in
this tutorial. (I will probably write this section in the future). This
program generated additional files relative to ***vcf2popNew.1.0.py***
that can be directly passed to Onemap or JoinMap genetic mapping algorithm
(see README for detailed informations).

***vcf2popNew.1.0.py*** can work on any segregating population which
implicated bi-allelic markers but required parent to be present in the
vcf file as well as the progeny and do not return file usable by Onemap
and JoinMap. But they can be easily eddited with excel or any open source
table manager to match their standards.

Both programs work on the same principle and are based on a step by step
selection: 

**(1)** Each data-point of the vcf is analyzed as follows: data-point/genotype
sufficiently covered and for which each allele are also sufficiently covered
are conserved. Otherwise the data-point is converted to missing data. 

**(2)** for each variant line, if to much missing data are found, the
variant site is removed.

**(3)** for each variant line passing the missing data threshold, the
marker segregation is compared to the ones proposed to the script using
a Chi-square test.

**(4)** marker having a p-value superior to the p-value threshold passed
to the script for a segregation are identified as markers corresponding
to this segregation. If several segregations are possible according to
p-values passed, marker is attributed to the segregation which have the
higher p-value.

**(5)** if a bi-parental population is studied, markers are attributed to
parental markers according to following rules: (i) if both parent are
heterozygous, markers are attributed to bridge markers if the segregation
corresponded, else it is attributed to a file containing marker of unknown
origin. (ii) if one parent is heterozygous and the other one homozygous,
marker is attributed to the heterozygous one. (iii) if one parent is
homozygous and the other one is missing data, marker is attributed to
the parent with missing data as if the marker segregate and one parent
is homozygous, the second one should be heterozygous. (iv) in a similar
way, if one parent is heterozygous and the other one is missing data,
marker is attributed to the heterozygous parent. (v) if all parent have
missing data, marker is attributed to a file containing marker of unknown
origin.


### Analysis:
To run the analysis, run the following command line:

    python3 ../bin/vcf2popNew.1.0.py -v ../data/vcf/Carto.vcf.gz -S SimpleDose:P1,P2:Ho,He@nn,np:0.5,0.5:1e-10/Bridge:P1,P2:Ho,He,Ho@hh,hk,kk:0.25,0.5,0.25:1e-10/DoubleDose:P1,P2:Ho,He@hh,k-:0.1667,0.8333:1e-7 -m 15 -M 300 -f 0.01:0.1 -c 3 -s 0.05 -o Pop

The ***-v*** option tells the program that the vcf it will be working on is here ***../data/vcf/Carto.vcf.gz***.

The ***-S*** option tells the program that 3 distinct segregations will be tested:

(1) ***SimpleDose:P1,P2:Ho,He@nn,np:0.5,0.5:1e-10*** the first segregation tested is names ***SimpleDose***, the first parent is ***P1*** and the second is ***P2***. Only one type of homozygous genotype could be found and homozygous (***Ho***) genotypes should be coded ***nn*** and their expected proportion is ***0.5***, heterozygous (***He***) genotypes should be coded ***nn*** and their expected proportion is ***0.5***, the p-value threshold to retain this marker to belong to this segregation is ***1e-10***

(2) ***Bridge:P1,P2:Ho,He,Ho@hh,hk,kk:0.25,0.5,0.25:1e-10*** the second segregation tested is names ***Bridge***, the first parent is ***P1*** and the second is ***P2***. Two type of homozygous genotype can be found: first homozygous type (***Ho***) genotypes should be coded ***hh*** and their expected proportion is ***0.25***, heterozygous (***He***) genotypes should be coded ***hk*** and their expected proportion is ***0.5***, second homozygous type (***Ho***) genotypes should be coded ***kk*** and their expected proportion is ***0.25***, the p-value threshold to retain this marker to belong to this segregation is ***1e-10***

(3) ***DoubleDose:P1,P2:Ho,He@hh,k-:0.1667,0.8333:1e-7*** the third segregation tested is names ***DoubleDose***, the first parent is ***P1*** and the second is ***P2***. Only one type of homozygous genotype could be found and homozygous (***Ho***) genotypes should be coded ***hh*** and their expected proportion is ***0.1667***, heterozygous (***He***) genotypes should be coded ***k-*** and their expected proportion is ***0.8333***, the p-value threshold to retain this marker to belong to this segregation is ***1e-7***

Why three segregation are tested here? Remember, we are working on a cross implication a diploid parent and a tetraploid one. So, for the diploid parent, only two segregations are possible, the one named ***SimpleDose*** (heterozygous marker in the diploid and homozygous in the polyploid) and the one named ***Bridge*** (heterozygous in both parent). However, for the tetraploid parent the number of haplotypes having each alleles when it is heterozygous has an impact on the marker segregation. I will pass over the polysomal versus disomal chromosomal pairing and directly say that we are in polysomal pairing of (home)omologs. If only one haplotype as the alternative allele, and in case the diploid parent is homozygous, the expected segregation is the ***SimpleDose*** one. But if two haplotypes have the alternate allele and the diploid parent is heterozygous, the expected segregation is the ***DoubleDose*** one.

The ***-m*** option tells the program that data-point covered by less than ***15*** reads were converted to missing.

The ***-M*** option tells the program that data-point covered by more than ***300*** reads were converted to missing.

The ***-f*** option tells the program accession having a minor allele frequency comprised between ***0.01*** and ***0.1*** were converted to missing. This interval represent the interval in which we consider that we cannot discriminate between homozygous and heterozygous state of the data-point.

The ***-c*** option tells the program that data-point having a minor allele supported by less than ***3*** reads were converted to missing data.

The ***-s*** option tells the program that variant site with more than ***0.05*** (i.e. 5%) missing data were not analyzed.

The ***-o*** option tells the program output files should be prefixed by ***Pop***.

This programs outpouts 4 files:

-   **Pop_report.tab** a file reporting statistics on marker selection.
-   **Pop_sub.vcf** is a sub vcf file of marker passing filters.
-   **Pop.tab** is file showing genotype for selected markers, best Chi-square value and p-value.
-   **Pop_tab_Bridge.tab**: is a marker coded file corresponding to selected marker identified as the ***Bridge*** segregation.
-   **Pop_tab_DoubleDose_P1.tab**: is a marker coded file corresponding to selected marker identified as the ***DoubleDose*** segregation and hetrozygous in ***P1***.
-   **Pop_tab_DoubleDose_P2.tab**: is a marker coded file corresponding to selected marker identified as the ***DoubleDose*** segregation and hetrozygous in ***P2***.
-   **Pop_tab_SimpleDose_P1.tab**: is a marker coded file corresponding to selected marker identified as the ***SimpleDose*** segregation and hetrozygous in ***P1***.
-   **Pop_tab_SimpleDose_P2.tab**: is a marker coded file corresponding to selected marker identified as the ***SimpleDose*** segregation and hetrozygous in ***P2***.
-   **Pop_tab_unknown.tab**: is a marker coded file corresponding to selected marker that cannot be attributed to one of the two parent.


B- Calculating recombination rate and segregation distortion:
-------------------------------------------------------------

At this stage we thus have several file containing segregating markers in distinct dose for both parent. P2 is the diploid parent while P1 is the tetraploid one. For simplicity and to make this tutorial available as soon as possible we will work on SimpleDose markers from the diploid parent (***Pop_tab_SimpleDose_P2.tab*** file).

The aims of this section is to calculate pairwise marker recombination rate and marker segregation distortion relative.
Pairwise recombination rate can be calculated with the following command line:

    python3 ../bin/RecombCalculatorDDose.py -m Pop_tab_SimpleDose_P2.tab -o Pop_tab_SimpleDose_P2 -p n -s R

The output is a file named ***Pop_tab_SimpleDose_P2_REC.tab*** containing the matrix of pairwise marker correlation.

Segregation distortion can be calculated with the following command line:

    python3 ../bin/RecombCalculatorDDose.py -m Pop_tab_SimpleDose_P2.tab -o Pop_tab_SimpleDose_P2 -p n -s S

The output is a file named ***Pop_tab_SimpleDose_P2_SegDist.tab*** containing for each marker its level of segregation distortion calculated as -log10(Chi-square test).


C- Calculating recombination rate and segregation distortion:
-------------------------------------------------------------

Once all these informations are calculated what we want to doo is to represent these data. in a comprehensive way. To do so, an additional file locating markers along chromosome is required. Luckily, this information is contained in marker name! The required file should contain in column1 marker name, in column2 chromosome name and in column2 position on chromosome. With our data this file can be generated with the following command line:

    cut -f 1 Pop_tab_SimpleDose_P2.tab | grep -v Marker > toto1
    cut -f 1 Pop_tab_SimpleDose_P2.tab | grep -v Marker | sed 's/M/\t'/ > toto2
    paste toto1 toto2 > MarkerOrder.tab
	rm toto1 toto2

To draw the picture run the following command line (This can be long when there are lots of point to draw, here around 4 mins):

    python3 ../bin/Draw_dot_plot.py -m Pop_tab_SimpleDose_P2_REC.tab -l MarkerOrder.tab -o DiploNoPhys.png -s Pop_tab_SimpleDose_P2_SegDist.tab -p n

This command line outputs two files:

-   **DiploNoPhys_heatmap.png** a png file showing color code relative to the linkage coding.
-   **DiploNoPhys.png** a png file marker linkage ordered along chromosomes as well as the marker segregation distortion calculated earlier.

The following picture is the picture you should obtain. In this picture you can observe the marker linkage of marker ordered along chromosome 1, 2 and 3. As the ***-p n*** option as been passed, marker spacing equal between two contiguous markers. Another way to say it is that marker spacing is not the physical one. Each dot represent the marker linkage intensity between two markers. A warm color represent a strong linkage, a cold one represent a weak linkage. The segregation distortion calculated earlier is represented on the graph at the right of the picture.
![](http://banana-genome-http.cirad.fr/image/DrawDotPlot1.png)









