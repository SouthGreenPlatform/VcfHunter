#!/usr/bin/env python3

#-----------
# Parameters

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c","--conf-file",
                    help="conf file, list of vcf path, one per line",
                    type=str,required=True)
parser.add_argument("-g","--group-file",
                    help="group file, with a pair accession - group per line",
                    type=str,required=True)
parser.add_argument("-i","--input-dir",
                    help="input directory, default is \"./step2\"",
                    type=str,default="./step2/")
parser.add_argument("-o","--output-dir",
                    help="output directory, default is \"./step3\"",
                    type=str,default="./step3/")
parser.add_argument("-a","--accession",
                    help="accession to scan from vcf",
                    type=str,required=True)
parser.add_argument("-d","--depth",
                    help="Minimal depth to consider a site. value comprised between 1 to infinite. If only considering sites with genotype, omit this argument.",
                    type=str,default="n")
args = parser.parse_args()

#-------
# Import

import os
import re
import sys
import gzip

#------
# Input

conf_file = args.conf_file
group_file = args.group_file
in_dir = args.input_dir
wk_dir = args.output_dir
individual = args.accession
depthValue = args.depth

#-------
# Output

os.makedirs(wk_dir,exist_ok=True)
out_file = os.path.join(wk_dir,individual + "_ratio.tab.gz")

#----------
# Conf file

vcf_list = []
print("Read",conf_file)
with open(conf_file,"r") as cf:
    for line in cf:
        line = line.rstrip()
        if line:
            vcf_list.append(line)

#-----------------
# Group ratio file

group_list = []
print("Read",group_file)
with open(group_file,"r") as gf:
    for line in gf:
        line = line.rstrip()
        group = line.split("\t")[1]
        if not group in group_list:
          group_list.append(group)

#-----------------
# Load group ratio

grp_ratio_d = {}
for group in group_list:
    gfile = os.path.join(in_dir,group+"_ratio.tab.gz")
    print("Read",gfile)
    with gzip.open(gfile,"rt") as g:
        for line in g:
            line = line.rstrip()
            (chrm,pos,allele,grp,ratio,ancnum) = line.split()
            if not (chrm,pos,allele) in grp_ratio_d:
                grp_ratio_d[(chrm,pos,allele)] =  {}
            grp_ratio_d[(chrm,pos,allele)][grp] = ratio

#check if an allele is in more than one group
to_remove = set()
for key in grp_ratio_d:
    if len(grp_ratio_d[key]) > 1:
        to_remove.add(key)

for key in to_remove:
    del grp_ratio_d[key]

len_tr = len(to_remove)
if len_tr:
    print("Warning: there is ",len_tr," alleles belonging",
          "to more than one group")

#----------
# Parse VCF

CHROM = 0
POS = 1
ID = 2
REF = 3
ALT = 4
QUAL = 5
FILTER = 6
INFO = 7
FORMAT = 8


#Write output file header
with gzip.open(out_file,"wt") as outf:
    header = "\t".join(["chr","pos","allele","obs_ratio","exp_ratio","grp"])
    header += "\n"
    outf.write(header)

#loop on vcf
for vcf_file in vcf_list:
    print("Read",vcf_file)
    #First we search for '#CHR ...' line
    if vcf_file[-2:] == 'gz':
        with gzip.open(vcf_file,"rt") as vf:
            line = vf.readline().rstrip()
            search_header = True
            while search_header:
                if re.match(r"^#CHROM",line):
                    search_header = False
                    #get individual we are looking for
                    ind_pos = line.split().index(individual)
                    #keep position to restart after header
                    vcf_data_pos = vf.tell()
                line = vf.readline().rstrip()
    else:
        with open(vcf_file,"r") as vf:
            line = vf.readline().rstrip()
            search_header = True
            while search_header:
                if re.match(r"^#CHROM",line):
                    search_header = False
                    #get individual we are looking for
                    ind_pos = line.split().index(individual)
                    #keep position to restart after header
                    vcf_data_pos = vf.tell()
                line = vf.readline().rstrip()
    # Use only genotyped sites
    if depthValue == 'n':
        #Then we parse the data
        if vcf_file[-2:] == 'gz':
            with gzip.open(vcf_file,"rt") as vf, \
                 gzip.open(out_file,"at") as outf:
                vf.seek(vcf_data_pos)
                for line in vf:
                    #split line
                    spl_line = line.split()
                    f_chr = spl_line[CHROM]
                    f_ref = spl_line[REF]
                    f_alt = spl_line[ALT]
                    f_pos = spl_line[POS]
                    f_format = spl_line[FORMAT]
                    f_data = spl_line[ind_pos]
                    alleles = tuple([f_ref] + f_alt.split(","))
                    #get AD and DP position and extract them
                    spl_format = f_format.split(":")
                    spl_data = f_data.split(":")
                    # work only on genotyped lines
                    if not('.' in spl_data[spl_format.index("GT")]):
                        #search if variant in table
                        for allele_pos in range(len(alleles)):
                            allele = alleles[allele_pos]
                            key = (f_chr,f_pos,allele)
                            if key in grp_ratio_d:
                                #get the good AD field
                                AD_tmp = spl_data[spl_format.index("AD")]
                                AD_val = AD_tmp.split(",")[allele_pos]
                                DP_val = spl_data[spl_format.index("DP")]
                                #get group and expected ratio
                                subkey = list(grp_ratio_d[key].keys())
                                group = subkey[0]
                                exp_ratio = grp_ratio_d[key][group]
                                #print
                                ratio = str(float(AD_val)/float(DP_val))
                                str_out = "\t".join([f_chr,f_pos,
                                                    allele,ratio,
                                                    exp_ratio,group])
                                outf.write(str_out + "\n")


            
        else:
            with open(vcf_file,"r") as vf, \
                 gzip.open(out_file,"at") as outf:
                vf.seek(vcf_data_pos)
                for line in vf:
                    #split line
                    spl_line = line.split()
                    f_chr = spl_line[CHROM]
                    f_ref = spl_line[REF]
                    f_alt = spl_line[ALT]
                    f_pos = spl_line[POS]
                    f_format = spl_line[FORMAT]
                    f_data = spl_line[ind_pos]
                    alleles = tuple([f_ref] + f_alt.split(","))
                    #get AD and DP position and extract them
                    spl_format = f_format.split(":")
                    spl_data = f_data.split(":")
                    # work only on genotyped lines
                    if not('.' in spl_data[spl_format.index("GT")]):
                        #search if variant in table
                        for allele_pos in range(len(alleles)):
                            allele = alleles[allele_pos]
                            key = (f_chr,f_pos,allele)
                            if key in grp_ratio_d:
                                #get the good AD field
                                AD_tmp = spl_data[spl_format.index("AD")]
                                AD_val = AD_tmp.split(",")[allele_pos]
                                DP_val = spl_data[spl_format.index("DP")]
                                #get group and expected ratio
                                subkey = list(grp_ratio_d[key].keys())
                                group = subkey[0]
                                exp_ratio = grp_ratio_d[key][group]
                                #print
                                ratio = str(float(AD_val)/float(DP_val))
                                str_out = "\t".join([f_chr,f_pos,
                                                    allele,ratio,
                                                    exp_ratio,group])
                                outf.write(str_out + "\n")
    else:
        MinDPValue = int(depthValue)
        #Then we parse the data
        if vcf_file[-2:] == 'gz':
            with gzip.open(vcf_file,"rt") as vf, \
                 gzip.open(out_file,"at") as outf:
                vf.seek(vcf_data_pos)
                for line in vf:
                    #split line
                    spl_line = line.split()
                    f_chr = spl_line[CHROM]
                    f_ref = spl_line[REF]
                    f_alt = spl_line[ALT]
                    f_pos = spl_line[POS]
                    f_format = spl_line[FORMAT]
                    f_data = spl_line[ind_pos]
                    alleles = tuple([f_ref] + f_alt.split(","))
                    #get AD and DP position and extract them
                    spl_format = f_format.split(":")
                    spl_data = f_data.split(":")
                    # work only on sufficiently covered sites
                    DP_val = int(spl_data[spl_format.index("DP")])
                    if DP_val >= MinDPValue:
                        #search if variant in table
                        for allele_pos in range(len(alleles)):
                            allele = alleles[allele_pos]
                            key = (f_chr,f_pos,allele)
                            if key in grp_ratio_d:
                                #get the good AD field
                                AD_tmp = spl_data[spl_format.index("AD")]
                                AD_val = AD_tmp.split(",")[allele_pos]
                                #get group and expected ratio
                                subkey = list(grp_ratio_d[key].keys())
                                group = subkey[0]
                                exp_ratio = grp_ratio_d[key][group]
                                #print
                                ratio = str(float(AD_val)/float(DP_val))
                                str_out = "\t".join([f_chr,f_pos,
                                                    allele,ratio,
                                                    exp_ratio,group])
                                outf.write(str_out + "\n")


            
        else:
            with open(vcf_file,"r") as vf, \
                 gzip.open(out_file,"at") as outf:
                vf.seek(vcf_data_pos)
                for line in vf:
                    #split line
                    spl_line = line.split()
                    f_chr = spl_line[CHROM]
                    f_ref = spl_line[REF]
                    f_alt = spl_line[ALT]
                    f_pos = spl_line[POS]
                    f_format = spl_line[FORMAT]
                    f_data = spl_line[ind_pos]
                    alleles = tuple([f_ref] + f_alt.split(","))
                    #get AD and DP position and extract them
                    spl_format = f_format.split(":")
                    spl_data = f_data.split(":")
                    # work only on sufficiently covered sites
                    DP_val = int(spl_data[spl_format.index("DP")])
                    if DP_val >= MinDPValue:
                        #search if variant in table
                        for allele_pos in range(len(alleles)):
                            allele = alleles[allele_pos]
                            key = (f_chr,f_pos,allele)
                            if key in grp_ratio_d:
                                #get the good AD field
                                AD_tmp = spl_data[spl_format.index("AD")]
                                AD_val = AD_tmp.split(",")[allele_pos]
                                #get group and expected ratio
                                subkey = list(grp_ratio_d[key].keys())
                                group = subkey[0]
                                exp_ratio = grp_ratio_d[key][group]
                                #print
                                ratio = str(float(AD_val)/float(DP_val))
                                str_out = "\t".join([f_chr,f_pos,
                                                    allele,ratio,
                                                    exp_ratio,group])
                                outf.write(str_out + "\n")