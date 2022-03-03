#!/usr/bin/env python3
#
#  Copyright 2014 CIRAD
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/> or
#  write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# -*- coding: utf-8 -*-

print("---------------------------------------------------------------")
print("WARNING:\n", 
	  "this script don't handle missing allele ratio from each group,\n",
	  "and set them to '0'. With this we assume that there are no \n",
	  "missing data in ancestral accession.")
print("---------------------------------------------------------------")

#-----------
# Parameters

#need one import
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-g","--group-file",
					help="group file, with a pair accession - group per line",
					type=str,required=True)
parser.add_argument("-p","--stat-file-pattern",
					help="pattern of the statistics file from vcf2allPropAndCov",
					type=str,default="_AlleleOriginAndRatio.tab.gz")
parser.add_argument("-i","--input-dir",
					help="input directory, default is \"./step1\"",
					type=str,default="./step1/")
parser.add_argument("-e","--excl",
					help="A file containing region to exclude for ancestry attribution in some introgressed accessions.",
					type=str,default=None)
parser.add_argument("-o","--output-dir",
					help="output directory, default is \"./step2\"",
					type=str,default="./step2/")
args = parser.parse_args()


#-------
# Import

import os
import sys
import gzip

#------
# Input

group_file = args.group_file
stat_file_tail = args.stat_file_pattern
in_dir = args.input_dir
wk_dir = args.output_dir

#-------------------------------------------------------
# Loading introgression information if a file was passed

DicoRegion = {}
if args.excl != None:
	file = open(args.excl)
	for line in file:
		data = line.split()
		if data:
			if len(data) != 8:
				sys.exit('The program exited without finishing because the following line does not have a correct number of columns:\n'+line)
			else:
				START = int(data[6])
				END = int(data[7])
				CHR = data[1]
				ACC = data[0]
				if not (ACC in DicoRegion):
					DicoRegion[ACC] = {}
				if not (CHR in DicoRegion[ACC]):
					DicoRegion[ACC][CHR] = set()
				DicoRegion[ACC][CHR].add((START, END))
	file.close()

#-----------
# Group file

#store accession name per group
group_d = {}
with open(group_file,"r") as gf:
	for line in gf:
		line.rstrip()
		(acc,group) = str.split(line)
		if(group in group_d):
			group_d[group].append(acc)
		else:
			group_d[group]=[acc]

#-----------
# Group file

os.makedirs(wk_dir,exist_ok=True)

#get all file from a group and write mean ratio
for group in group_d.keys():
	print("Compute ratio for group",group)
	ratio_d = {}
	grplen = len(group_d[group])
	
	#store all data from each file
	for acc in group_d[group]:
		print("Read file",acc)
		stat_file = (os.path.join(in_dir,acc,acc + stat_file_tail))
		with gzip.open(stat_file,"rt") as sf:
			for line in sf:
				line.rstrip()
				(chro,pos,allele,grp,ratio) = str.split(line)
				key=(chro,pos,allele,grp)
				
				# Verification that we are not in introgressed region for the accession
				notonintrog = True
				if acc in DicoRegion:
					if key[0] in DicoRegion[acc]:
						for reg in DicoRegion[acc][key[0]]:
							if reg[0] < int(key[1]) and int(key[1]) < reg[1]:
								notonintrog = False
				if notonintrog:
					if(key in ratio_d):
						ratio_d[key] += float(ratio)
					else:
						ratio_d[key] = float(ratio)
	
	#print result
	out_file=group+"_ratio.tab.gz"
	print("Write",out_file)
	with gzip.open(os.path.join(wk_dir,out_file),"wt") as outf:
		for key in ratio_d.keys():
			# Counting ancestor in introgressed region
			dicoIntro = set()
			for acc in group_d[group]:
				if acc in DicoRegion:
					if key[0] in DicoRegion[acc]:
						for reg in DicoRegion[acc][key[0]]:
							if reg[0] < int(key[1]) and int(key[1]) < reg[1]:
								dicoIntro.add(acc)
			outf.write("{}\t{}\t{}\t{}".format(*key)+"\t"+str(ratio_d[key]/(grplen-len(dicoIntro)))+"\t"+str(grplen-len(dicoIntro))+"\n")
	print("---")

