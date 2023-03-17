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
import datetime
import traceback
import multiprocessing as mp
from inspect import currentframe, getframeinfo
from operator import itemgetter

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
from scipy.interpolate import spline
from textwrap import wrap


def draw_chr(DICO_INFO, CHR_INFO, DICO_GROUP, OUT, CHROM, PLOIDY, ACC_order):
	
	# Color definition
	color = ((0,0.8,0), (1,0,0),(0,0,1))

	NB=15

	# Calcul de la taille max
	MAXI = CHR_INFO[CHROM]

	# Getting group order
	group = sorted(list(DICO_GROUP.keys()))

	# Calculating ticks
	i = 5000000
	ticks_pos = [0]
	ticks_labels = [0]

	while i < MAXI:
		ticks_pos.append(i)
		ticks_labels.append(str(int(i/1000000)))
		i += 5000000
	i = 1000000
	
	minor_ticks = []
	while i < MAXI:
		minor_ticks.append(i)
		i += 1000000
	
	######Drawing allele ratio######
	
	# On fait la figure
	POSSPAN = 0
	fig = plt.figure(figsize=(10.5, 14.85))
	fig.subplots_adjust(left=0.015, right=0.985, top=0.98, bottom=0.05)

	
	for i, acc in enumerate(ACC_order):
		# print (i, acc)
		ax = plt.subplot2grid((NB,15),(POSSPAN,1), colspan=14, rowspan=1)
		ax.set_ylim(0, 1.05)
		ax.set_xlim(0, MAXI)
		ax.axhline(y=0.25, linewidth=0.1, color = (0,0,0.5))
		ax.axhline(y=0.75, linewidth=0.1, color = (0,0,0.5))
		ax.axhline(y=0.5, linewidth=0.1, color = (0,0,0.5))
		ax.axhline(y=0.33, linewidth=0.1, color = (0.5,0.5,0.5))
		ax.axhline(y=0.66, linewidth=0.1, color = (0.5,0.5,0.5))
		for gp  in group:
			# Getting position
			if not (acc in DICO_INFO):
				value = [0]
				final_pos = [0]
			elif not(CHROM in DICO_INFO[acc]):
				value = [0]
				final_pos = [0]
			else:
				position = sorted(list(DICO_INFO[acc][CHROM].keys()))
				#print position
				value = []
				final_pos = []
				for pos in position:
					#print "coucou for pos in position"
					if gp in DICO_INFO[acc][CHROM][pos]:
						value.append(DICO_INFO[acc][CHROM][pos][gp])
						final_pos.append(pos)
			ax.plot(final_pos, value, 'o', ms=1.5, mew=0, mfc=color[group.index(gp)])
			ax.axes.yaxis.set_ticklabels([])
			ax.axes.xaxis.set_ticklabels([])
		
		# if CHROM == "chr01":
			# ax.axvline(x=650000, ymin=0, ymax = 1.05, linewidth=1, color='k')
		# elif CHROM == "chr03":
			# ax.axvline(x=26675034, ymin=0, ymax = 1.05, linewidth=1, color='k')

		ax = plt.subplot2grid((NB,15),(POSSPAN,0), colspan=1, rowspan=1)
		ax.axis('off')	
		ax.axis([0, 1, 0, 1])
		ax.text(0, 0.5, acc, size=12, va='center', fontweight='bold')
		

		POSSPAN += 1	
		if (i+1)%15==0:   
			OUT2=OUT+'_'+CHROM+'_'+str(int(i/15)+1)+'_Ratio.png'
			fig.savefig(OUT2)
			plt.close(fig)

			POSSPAN = 0
			fig = plt.figure(figsize=(10.5, 14.85))
			fig.subplots_adjust(left=0.015, right=0.985, top=0.98, bottom=0.05)
		j=i	

	# print (j)
	if (i+1)%15 != 0:
		fig.savefig(OUT+'_'+CHROM+'_'+str(int(j/15)+1)+'_Ratio.png')
		plt.close(fig)	

	######Drawing coverage######
	
	# On fait la figure
	POSSPAN = 0
	fig = plt.figure(figsize=(10.5, 14.85))
	fig.subplots_adjust(left=0.015, right=0.985, top=0.98, bottom=0.05)

	
	for i, acc in enumerate(ACC_order):
		MEAN_COV = sum(DICO_INFO[acc]['total'])/float(len(DICO_INFO[acc]['total']))
		#print (i, acc)
		ax = plt.subplot2grid((NB,15),(POSSPAN,1), colspan=14, rowspan=1)
		ax.set_ylim(0, 2*MEAN_COV)
		ax.set_xlim(0, MAXI)
		# Getting position
		if not (acc in DICO_INFO):
			value = [0]
			final_pos = [0]
		elif not(CHROM in DICO_INFO[acc]):
			value = [0]
			final_pos = [0]
		else:
			position = sorted(list(DICO_INFO[acc][CHROM].keys()))
			#print position
			value = []
			final_pos = []
			fill_plus1 = [(PLOIDY+0.5)/PLOIDY*MEAN_COV]*(len(DICO_INFO[acc][CHROM])+2) ## For band
			fill_moins1 = [(PLOIDY-0.5)/PLOIDY*MEAN_COV]*(len(DICO_INFO[acc][CHROM])+2) ## For band
			for pos in position:
				value.append(DICO_INFO[acc][CHROM][pos]['cov'])
				final_pos.append(pos)
			ax.fill_between([0]+final_pos+[CHR_INFO[CHROM]], fill_plus1, fill_moins1, color=(1,0,0), alpha=0.20) ## For band
		ax.axhline(y=MEAN_COV, linewidth=0.5, color = (0,1,0))
		ax.axhline(y=(PLOIDY+1)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,1,0))
		ax.axhline(y=(PLOIDY-1)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,1,0))
		ax.axhline(y=(PLOIDY+2)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,0,1))
		ax.axhline(y=(PLOIDY-2)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,0,1))
		ax.plot(final_pos, value, 'o', ms=1.5, mew=0, mfc='black')
		ax.axes.yaxis.set_ticklabels([])
		ax.axes.xaxis.set_ticklabels([])
		
		# if CHROM == "chr01":
			# ax.axvline(x=650000, ymin=0, ymax = 2*MEAN_COV, linewidth=1, color='k')
		# elif CHROM == "chr03":
			# ax.axvline(x=26675034, ymin=0, ymax = 2*MEAN_COV, linewidth=1, color='k')

		ax = plt.subplot2grid((NB,15),(POSSPAN,0), colspan=1, rowspan=1)
		ax.axis('off')	
		ax.axis([0, 1, 0, 1])
		ax.text(0, 0.5, acc, size=12, va='center' ,fontweight='bold')

		POSSPAN += 1	
		if (i+1)%15==0:   
			OUT2=OUT+'_'+CHROM+'_'+str(int(i/15)+1)+'_Cov.png'
			fig.savefig(OUT2)
			plt.close(fig)

			POSSPAN = 0
			fig = plt.figure(figsize=(10.5, 14.85))
			fig.subplots_adjust(left=0.015, right=0.985, top=0.98, bottom=0.05)
		j=i	

	print (j)
	if (i+1)%15 != 0:
		fig.savefig(OUT+'_'+CHROM+'_'+str(int(j/15)+1)+'_Cov.png')
		plt.close(fig)

def get_chr_size(LINE):
	
	"""
		Get reference chromosome size from vcf file
		
		:param LINE: The ##contig line of the vcf file splitted on spaces
		:type LINE: list
		:return: A list with [0] --> chromosome name, [1] --> chromosome size
		:rtype: list
	"""
	
	split_on_chev = LINE[0].split('>')[0]
	
	split_on_eq = split_on_chev.split('=')
	
	return [split_on_eq[2].split(',')[0], int(split_on_eq[3])]

def __main__():

	#print "prout"

	#Parse Command Line
	parser = optparse.OptionParser(
		description="This program perform two things based on a vcf. 1) It plots for a chromosome of all "
		"accessions in a vcf, the allele coverage along its chromosomes. 2) It identify, based on known "
		"ancestral accessions in the vcf, the alleles specific to each groups and plot the alleles "
		"proportion at a site along chromosomes for all accessions.",
		epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr), Marion Dupouy and Franc-Christophe Baurens (baurens@cirad.fr)")

	# Wrapper options. 
	parser.add_option( '',	'--conf',			dest='conf',		default=None,			help='Conf file containing vcf(s) location(s)')
	parser.add_option( '',	'--vcf',			dest='vcf',			default=None,			help='Path to uniq vcf file. (--conf and --vcf are mutually exclusive). If --vcf option is passed, --conf will be omited')
	parser.add_option( '',	'--origin',			dest='origin',		default=None,			help='A 2 column file containing accession name (col1) origin/group (Col2). [Default: %default]')
	parser.add_option( '',	'--ploidy',			dest='ploidy',		default=None,			help='Accession ploidy')
	parser.add_option( '',	'--NoMiss',			dest='NoMiss',		default='n',			help='No missing data are allowed in accessions used to group alleles. [Default: %default]')
	parser.add_option( '',	'--all',			dest='all',			default='n',			help='Allele should be present in all accessions of the group. [Default: %default]')
	parser.add_option( '',	'--acc',			dest='acc',			default=None,			help='Accession to work with. If ignored, all accessions in the vcf will be used. Else accessions should be separated by ",". [Default: %default]')
	parser.add_option( '',	'--prefix',			dest='prefix',		default='RatioAndCov',	help='Prefix for output file name. [Default: %default]')
	(options, args) = parser.parse_args()

	if options.conf == None and options.vcf == None:
		sys.exit('Please provide either a conf file to --conf argument or a vcf path to --vcf argument')
	if options.origin == None:
		sys.exit('Please provide a origin file to --origin argument')
	if options.ploidy == None:
		sys.exit('Please provide a ploidy level to --ploidy argument')

	# recording accession group
	dico_group = {}
	file = open(options.origin)
	for line in file:
		data = line.split()
		if data:
			if not(data[1] in dico_group):
				dico_group[data[1]] = []
			dico_group[data[1]].append(data[0])
	file.close()

	# recording vcf file to work with
	dico_vcf = set()
	if options.conf == None:
		dico_vcf.add(options.vcf)
	else:
		file = open(options.conf)
		for line in file:
			data = line.split()
			if data:
				dico_vcf.add(data[0])
		file.close()


	dico_chr = {}	#information from ##contig lines
	dico_acc = {}	#contain for each accession (keys) drawing information for the wanted chr
	
	recorded_chr = set()
	
	#pour chaque ligne:

	for vcf in dico_vcf:
		if vcf[-3:] == '.gz':
			file = gzip.open(vcf,'rt')
		else:
			file = open(vcf)
		for line in file:
			data = line.split()
			if data:
				#si ligne ##contig : stocker taille chr
				if data[0][0:9] == '##contig=':
					chr_info = get_chr_size(data)
					dico_chr[chr_info[0]] = chr_info[1]
				#si ligne #CHROM : stocker header (liste des accessions)
				elif data[0] == "#CHROM":
					header = data
					if options.acc == None:
						acc_list=list(sorted(header[header.index("FORMAT")+1:]))
					else:
						acc_list = options.acc.split(',')
				elif data[0][0] != "#":
					FORMAT = data[header.index("FORMAT")].split(':')
					CHR = data[header.index("#CHROM")]
					POS = int(data[header.index("POS")])
					REF = [data[header.index("REF")]]
					ALT = data[header.index("ALT")].split(',')
					VARIANT = REF + ALT

					# recording group1 and group2 alleles
					dico_allele = {}
					for gp in dico_group:
						dico_allele[gp] = set()
						for ACC in dico_group[gp]:
							#print ACC
							ACCESSION = data[header.index(ACC)].split(':')
							GENOTYPE = set(ACCESSION[FORMAT.index("GT")].replace('|','/').split('/'))
							for geno in GENOTYPE:
								dico_allele[gp].add(geno)
					
					# to manage missing data in groups
					Good = 1
					for gp in dico_allele:
						if options.NoMiss == 'n':
							if '.' in dico_allele[gp] and len(dico_allele[gp]) == 1:
								Good = 0
						elif options.NoMiss == 'y':
							if '.' in dico_allele[gp]:
								Good = 0
						else:
							sys.exit('Oups, their is a bug... either the vcf has a problem or I made a mistake in the programing')
					if Good:
						########################################################################################################################### Two distinct philosophy of work
						if options.all == 'y': # to manage alleles specific to all accessions of a group
							# looking for alleles specific of groups
							dico_alleles_groups = {}
							for gp in dico_allele:
								for allele in dico_allele[gp]:
									if allele != '.':
										if not (allele in dico_alleles_groups):
											dico_alleles_groups[allele] = set()
										dico_alleles_groups[allele].add(gp)
							
							dico_allele_specific = {}
							for allele in dico_alleles_groups:
								if len(dico_alleles_groups[allele]) == 1: #It is an allele specific of a group because it has been found only in one group
									GROUP = list(dico_alleles_groups[allele])[0]
									if not(GROUP in dico_allele_specific):
										dico_allele_specific[GROUP] = set(allele)
							
							
							# Now we should verify that all accessions (excluding accessions with missing data) have this allele
							DICO_ACC = {} # for counting the number of accession of a group without missing data
							DICO_ALL_IN_ACC_COUNT = {} # for counting the number of accession of a group without missing data
							for gp in dico_group:
								DICO_ACC[gp] = 0
								for ACC in dico_group[gp]:
									ACCESSION = data[header.index(ACC)].split(':')
									GENOTYPE = set(ACCESSION[FORMAT.index("GT")].replace('|','/').split('/'))
									if not('.' in  GENOTYPE):
										DICO_ACC[gp] += 1
										for geno in GENOTYPE:
											if not(geno in DICO_ALL_IN_ACC_COUNT):
												DICO_ALL_IN_ACC_COUNT[geno] = 0
											DICO_ALL_IN_ACC_COUNT[geno] += 1
							
							# Selecting allele to work with
							dico_allele = {}
							for gp in dico_allele_specific:
								dico_allele[gp] = set()
								ExpectedGpNumber = DICO_ACC[gp]
								for allele in dico_allele_specific[gp]:
									if DICO_ALL_IN_ACC_COUNT[allele] == ExpectedGpNumber:
										dico_allele[gp].add(allele)

						###################################################################################################################
						
						elif options.all == 'n':
							pass
						else:
							sys.exit('Oups, the program exited without finishing: please provide either "y" ot "n" to --all options\n')
						
						#pour chaque accession
						for acc in acc_list :
							ACCESSION = data[header.index(acc)].split(':')
							GENOTYPE = set(ACCESSION[FORMAT.index("GT")].replace('|','/').split('/'))
							COVERAGE = list(map(int, ACCESSION[FORMAT.index("AD")].split(',')))
							add_cov = 0
							for geno in GENOTYPE:
								if geno != '.':
									group = []
									for gp in dico_allele:
										if geno in dico_allele[gp]:
											group.append(gp)
									RATIO = COVERAGE[int(geno)]/float(sum(COVERAGE))
									if len(group) == 1:
										add_cov = 1
										if not(acc in dico_acc):
											dico_acc[acc] = {}
											dico_acc[acc]['total'] = []
										if not (CHR in dico_acc[acc]):
											dico_acc[acc][CHR] = {}
											recorded_chr.add(CHR)
											# print (CHR)
										if not(POS in dico_acc[acc][CHR]):
											dico_acc[acc][CHR][POS]={}
										dico_acc[acc][CHR][POS][group[0]] = RATIO
										# outfile.write('\t'.join([CHR, acc, str(POS), VARIANT[int(geno)], group[0], str(RATIO)]))
										# outfile.write('\n')
							if add_cov:
								dico_acc[acc][CHR][POS]['cov'] = sum(COVERAGE)
								dico_acc[acc]['total'].append(sum(COVERAGE))
								# print (len(dico_acc[acc]['total']))
		file.close()
	
	for chr in recorded_chr:
		draw_chr(dico_acc, dico_chr, dico_group, options.prefix, chr, int(options.ploidy), acc_list)

if __name__ == "__main__": __main__()







