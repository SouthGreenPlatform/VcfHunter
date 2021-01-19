#!/usr/bin/env python
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

import sys
sys.stdout.write('loading modules\n')
import optparse
import os
import shutil
import subprocess
import tempfile
import fileinput
import time
import random
import math
import gzip
sys.stdout.write('modules loaded\n')

def RecordAccession(VCF, NAMES, MAT):
	
	"""
		Fill a list with accession to work with
		
		:param NAMES: Path to file containing accession names
		:type NAMES: str
		:param VCF: Path to vcf file
		:type VCF: str
		:param MAT: Path to matrix file containing grouping informations
		:type MAT: str
		:return: A list containing accession to work with
		:rtype: list
	"""
	
	ACC_TO_WORK = []
	if NAMES == None:
		sys.stdout.write('No file name was provided in --names argument. All accessions in vcf file (if provided) will be drawn. If no vcf too, all accessions in matrix will be drawn.\n')
		if VCF == None:
			file = open(MAT)
			ACC_TO_WORK = file.readline().replace('K-mean_GROUP','').replace('GROUP','').split()
			file.close()
		else:
			if VCF[-3:] == '.gz':
				file = gzip.open(VCF,'rt')
			else:
				file = open(VCF)
			for line in file:
				data = line.split()
				if data:
					if data[0] == '#CHROM':
						ACC_TO_WORK = data[data.index('FORMAT')+1:]
						break
			file.close()
	else:
		file = open(NAMES)
		header = file.readline()
		for line in file:
			data = line.split()
			if data:
				ACC_TO_WORK.append(data[0])
		file.close()
	
	return ACC_TO_WORK

def CaclNbId(GENO1, GENO2):
	
	# Warning : stats to calculate also with different ploidy level in the sens "the diploid can be parent of the polyploid"
	
	geno1 = list(GENO1)
	geno2 = list(GENO2)
	while '.' in geno1:
		geno1.remove('.')
	while '.' in geno2:
		geno2.remove('.')
	ident = 0
	if len(geno1) > len(geno2):
		comp = len(geno2)
		while geno2:
			value = geno2[0]
			if value in geno1:
				ident += 1
				geno1.remove(value)
			del geno2[0]
	else:
		comp = len(geno1)
		while geno1:
			value = geno1[0]
			if value in geno2:
				ident += 1
				geno2.remove(value)
			del geno1[0]
	return [ident,comp]

def calculDis(VCF, NAMES, PREFIX):
	
	# recording accession names to draw
	acc_to_work = RecordAccession(VCF, NAMES, None)
	
	# Creating the count matrix
	dico_ident = {}
	dico_total = {}
	dico_dis = {}
	liste = [0]*len(acc_to_work)
	for i in range(len(acc_to_work)):
		dico_ident[i] = list(liste)
		dico_total[i] = list(liste)
		dico_dis[i] = list(liste)
	
	# filling the count matrix
	
	if VCF[-3:] == '.gz':
		file = gzip.open(VCF,'rt')
	else:
		file = open(VCF)
	
	j = 0
	dicoAccPos = {}
	AccNumber = len(acc_to_work)
	for line in file:
		data = line.split()
		if data:
			if data[0] == "#CHROM":
				header = data
				for i in range(len(acc_to_work)):
					name = acc_to_work[i]
					dicoAccPos[name] = header.index(name)
				FORMATPOS = header.index("FORMAT")
			elif not(data[0][0] == "#"):
				j += 1
				sys.stdout.write(str(j)+'\n')
				sys.stdout.flush()
				FORMAT = data[FORMATPOS].split(":")
				GTPOS = FORMAT.index('GT')
				GENOTYPES = []
				for i in range(AccNumber):
					name1 = acc_to_work[i]
					CALL1 = data[dicoAccPos[name1]].split(":")
					GENO1 = list(CALL1[GTPOS].replace('|','/').split('/'))
					GENOTYPES.append(GENO1)
				for i in range(AccNumber):
					for k in range(AccNumber):
						if k >= i:
							to_add = CaclNbId(GENOTYPES[i], GENOTYPES[k])
							dico_ident[i][k] += to_add[0]
							dico_total[i][k] += to_add[1]
							dico_ident[k][i] += to_add[0]
							dico_total[k][i] += to_add[1]
	
	# Calculating final values
	for i in dico_dis:
		for k in range(len(dico_dis[i])):
			print(acc_to_work[i], acc_to_work[k] , dico_total[i][k])
			if float(dico_total[i][k]) == 0:
				dico_dis[i][k] = '999'
				print('No dissimilarity found for accessions', acc_to_work[i], acc_to_work[k])
			else:
				dico_dis[i][k] = 1-(dico_ident[i][k]/float(dico_total[i][k]))
	
	# Printing results
	outfile = open(PREFIX+'.dis','w')
	debut = ['N']+list(range(len(acc_to_work)))
	outfile.write('\t'.join(list(map(str,debut))))
	outfile.write('\n')
	
	for i in range(len(acc_to_work)):
		mot = [str(i)]
		for k in range(len(acc_to_work)):
			mot.append(str(dico_dis[i][k]))
		outfile.write('\t'.join(mot))
		outfile.write('\n')
	outfile.close()
	
	# Printing correspondance file
	outfile = open(PREFIX+'.cor','w')
	file = open(NAMES)
	header = file.readline().split()
	ToPrint = ["N"]+ header
	outfile.write('\t'.join(ToPrint)+'\n')
	for line in file:
		data = line.split()
		if data:
			outfile.write(str(acc_to_work.index(data[0]))+'\t'+line)
	outfile.close()
	

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n"
	"This program calculate a simple matching dissimilarity between accessions passed in a vcf. The dissimilarity is calculated "
	"as the sum of 1-(shared alleles between two individuals) divided by the lower ploidy (of the two individuals) multiplication "
	"of the number of compared sites.")
	# Wrapper options. 
	parser.add_option( '-v',	'--vcf',			dest='vcf',			default=None,			help='The vcf file. [Default: %default]')
	parser.add_option( '-n',	'--names',			dest='names',		default=None,			help='A one column file containing accession names to treat. [Default: %default]')
	parser.add_option( '-p',	'--prefix',			dest='prefix',		default='WorkOnVcf', 	help='The prefix for output files. [Default: %default]')
	(options, args) = parser.parse_args()
	
	
	# Identify genome blocs in accessions and draw a circos by accessions
	if options.vcf == None:
		sys.exit('Please provide a matrix file to --mat argument')
	if options.names == None:
		sys.exit('Please provide a matrix file to --names argument')
	calculDis(options.vcf, options.names, options.prefix)
		
if __name__ == "__main__": __main__()