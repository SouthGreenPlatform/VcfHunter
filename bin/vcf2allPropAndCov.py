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


def draw_chr(DICO_INFO, CHR_INFO, DICO_GROUP, OUT, MEAN_COV, PLOIDY):
	
	
	
	# Color definition
	color = ((0,0.8,0),(1,0,0),(0,0,1),(0.780392156862745,0.0823529411764706,0.52156862745098),(0.541176470588235,0.211764705882353,0.0588235294117647),(1,0.498039215686275,0),(0.5,0.5,0.5),(1,1,0))
	
	# getting chromosomes list
	chr_list_temp = sorted(list(CHR_INFO.keys()))
	chr_list = []
	for chr in chr_list_temp:
		if chr in DICO_INFO:
			chr_list.append(chr)
	
	# On calcule le nombre de fenetres
	NB = len(chr_list)
	
	# Calcule de la taille max
	MAXI = 0
	for n in CHR_INFO:
		if MAXI < CHR_INFO[n]:
			MAXI = CHR_INFO[n]
	
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
	fig.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)
	
	# On dessine par chromosomes
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,1), colspan=14, rowspan=1)
		ax.set_ylim(0, 1.1)
		ax.set_xlim(0, MAXI)
		for gp  in group:
			# Getting position
			position = sorted(list(DICO_INFO[chr].keys()))
			value =  []
			final_pos = []
			ax.axhline(y=0.25, linewidth=0.1, color = (0,0,0.5))
			ax.axhline(y=0.75, linewidth=0.1, color = (0,0,0.5))
			ax.axhline(y=0.5, linewidth=0.1, color = (0,0,0.5))
			ax.axhline(y=0.33, linewidth=0.1, color = (0.5,0.5,0.5))
			ax.axhline(y=0.66, linewidth=0.1, color = (0.5,0.5,0.5))
			for pos in position:
				if gp in DICO_INFO[chr][pos]:
					value.append(DICO_INFO[chr][pos][gp])
					final_pos.append(pos)
			ax.plot(final_pos, value, 'o', ms=1.5, mew=0, mfc=color[group.index(gp)])
			ax.axes.yaxis.set_ticklabels([])
			ax.axes.xaxis.set_ticklabels([])
		
		# if chr == "chr01":
			# ax.axvline(x=650000, ymin=0, ymax = 1.05, linewidth=1, color='k')
		# elif chr == "chr03":
			# ax.axvline(x=26675034, ymin=0, ymax = 1.05, linewidth=1, color='k')
		
		POSSPAN += 1
	ax.set_xticks(ticks_pos)
	ax.set_xticks(minor_ticks, minor=True)
	ax.set_xticklabels(ticks_labels)
	
	# On mets le nom des chromosomes
	POSSPAN = 0
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,0), colspan=1, rowspan=1)
		ax.axis('off')	
		ax.axis([0, 1, 0, 1])
		ax.text(0, 0.5, chr, size=12, va='center', fontweight='bold')
		
		POSSPAN += 1
	
	fig.savefig(OUT+'Ratio.png')
	plt.close(fig)
	
	
	
	######Drawing coverage######
	
	# On fait la figure
	POSSPAN = 0
	fig = plt.figure(figsize=(10.5, 14.85))
	fig.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)
	
	# On dessine par chromosomes
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,1), colspan=14, rowspan=1)
		ax.set_ylim(0, 2*MEAN_COV)
		ax.set_xlim(0, MAXI)
		position = sorted(list(DICO_INFO[chr].keys()))
		value =  []
		final_pos = []
		fill_plus1 = [(PLOIDY+0.5)/PLOIDY*MEAN_COV]*(len(DICO_INFO[chr])+2) ## For band
		fill_moins1 = [(PLOIDY-0.5)/PLOIDY*MEAN_COV]*(len(DICO_INFO[chr])+2) ## For band
		fill_plus2 = [(PLOIDY+1.5)/PLOIDY*MEAN_COV]*(len(DICO_INFO[chr])+2) ## For band
		fill_moins2 = [(PLOIDY-1.5)/PLOIDY*MEAN_COV]*(len(DICO_INFO[chr])+2) ## For band
		for pos in position:
			value.append(DICO_INFO[chr][pos]['cov'])
			final_pos.append(pos)
		ax.fill_between([0]+final_pos+[CHR_INFO[chr]], fill_plus1, fill_moins1, color=(1,0,0), alpha=0.20) ## For band
		# ax.fill_between([0]+final_pos+[CHR_INFO[chr]], fill_plus2, fill_moins2, color=(0,0,1), alpha=0.10) ## for band
		ax.axhline(y=MEAN_COV, linewidth=0.5, color = (0,0,0))
		ax.axhline(y=(PLOIDY+1)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,1,0))
		ax.axhline(y=(PLOIDY-1)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,1,0))
		ax.axhline(y=(PLOIDY+2)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,0,1))
		ax.axhline(y=(PLOIDY-2)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,0,1))
		ax.plot(final_pos, value, 'o', ms=1.5, mew=0, mfc='black')
		ax.axes.yaxis.set_ticklabels([])
		ax.axes.xaxis.set_ticklabels([])
		
		# if chr == "chr01":
			# ax.axvline(x=650000, ymin=0, ymax = 2*MEAN_COV, linewidth=1, color='k')
		# elif chr == "chr03":
			# ax.axvline(x=26675034, ymin=0, ymax = 2*MEAN_COV, linewidth=1, color='k')
		
		POSSPAN += 1
	ax.set_xticks(ticks_pos)
	ax.set_xticks(minor_ticks, minor=True)
	ax.set_xticklabels(ticks_labels)
	
	# On mets le nom des chromosomes
	POSSPAN = 0
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,0), colspan=1, rowspan=1)
		ax.axis('off')	
		ax.axis([0, 1, 0, 1])
		ax.text(0, 0.5, chr, size=12, va='center', fontweight='bold')
		
		POSSPAN += 1
	
	fig.savefig(OUT+'Cov.png')
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
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n This program count allele ratio based on two origins")
	# Wrapper options. 
	parser.add_option( '',	'--conf',			dest='conf',		default=None,			help='Conf file containing vcf location (one per chromosome). [Default: %default]')
	parser.add_option( '',	'--origin',			dest='origin',		default=None,			help='A 2 column file containing accession name (col1) origin/group (Col2). [Default: %default]')
	parser.add_option( '',	'--acc',			dest='acc',			default=None,			help='Accession to work with. [Default: %default]')
	parser.add_option( '',	'--ploidy',			dest='ploidy',		default=None,			help='Accession ploidy')
	parser.add_option( '',	'--all',			dest='all',			default='n',			help='Allele should be present in all accessions to identify a group. [Default: %default]')
	(options, args) = parser.parse_args()
	
	if options.conf == None:
		sys.exit('Please provide a conf file to --conf argument')
	if options.origin == None:
		sys.exit('Please provide a origin file to --origin argument')
	if options.acc == None:
		sys.exit('Please provide a accession name to --acc argument')
	if options.ploidy == None:
		sys.exit('Please provide a ploidy level to --ploidy argument')
	# recording accession name
	ACCESS = options.acc
	
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
	
	outfile = open(ACCESS+'_AlleleOriginAndRatio.tab','w')
	outfile1 = open(ACCESS+'_stats.tab','w')
	
	# recording vcf file to work with
	file = open(options.conf)
	dico_vcf = set()
	for line in file:
		data = line.split()
		if data:
			dico_vcf.add(data[0])
	file.close()
	
	total_unmissing_sites = 0
	
	# For linear drawing
	dico_draw = {}
	dico_chr = {}
	total = []
	for vcf in dico_vcf:
		file = open(vcf)
		for line in file:
			data = line.split()
			if data:
				if data[0] == "#CHROM":
					header = data
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
							ACCESSION = data[header.index(ACC)].split(':')
							GENOTYPE = set(ACCESSION[FORMAT.index("GT")].replace('|','/').split('/'))
							for geno in GENOTYPE:
								dico_allele[gp].add(geno)
					
					# to manage missing data in groups
					Good = 1
					for gp in dico_allele:
						if options.all == 'n':
							if '.' in dico_allele[gp] and len(dico_allele[gp]) == 1:
								Good = 0
						elif options.all == 'y':
							if '.' in dico_allele[gp]:
								Good = 0
						else:
							sys.exit('Oups, their is a bug... either the vcf has a problem or I made a mistake in the programing')
					if Good:
						# recording accession 
						ACCESSION = data[header.index(ACCESS)].split(':')
						GENOTYPE = set(ACCESSION[FORMAT.index("GT")].replace('|','/').split('/'))
						COVERAGE = list(map(int, ACCESSION[FORMAT.index("AD")].split(',')))
						add_cov = 0
						if not ('.' in GENOTYPE):
							total_unmissing_sites += 1
						for geno in GENOTYPE:
							if geno != '.':
								group = []
								for gp in dico_allele:
									if geno in dico_allele[gp]:
										group.append(gp)
								# print('tutu', dico_allele, geno, group)
								RATIO = COVERAGE[int(geno)]/float(sum(COVERAGE))
								if len(group) == 1:
									add_cov = 1
									if not(CHR in dico_draw):
										dico_draw[CHR] = {}
									if not(POS in dico_draw[CHR]):
										dico_draw[CHR][POS] = {}
									dico_draw[CHR][POS][group[0]] = RATIO
									outfile.write('\t'.join([CHR, str(POS), VARIANT[int(geno)], group[0], str(RATIO)]))
									outfile.write('\n')
						if add_cov:
							dico_draw[CHR][POS]['cov'] = sum(COVERAGE)
							total.append(sum(COVERAGE))
				elif data[0][0:9] == '##contig=':
					chr_info = get_chr_size(data)
					dico_chr[chr_info[0]] = chr_info[1]
		file.close()
	outfile.close()
	
	# for global statistics
	dico_nb_group = {}
	dico_class_grouped = {}
	
	for gp in dico_group:
		dico_nb_group[gp] = 0
	
	for chr in dico_draw:
		for pos in dico_draw[chr]:
			liste_grouped = []
			for gp in dico_group:
				if gp in dico_draw[chr][pos]:
					liste_grouped.append(gp)
					dico_nb_group[gp] += 1
			liste_grouped.sort()
			couple = ':'.join(liste_grouped)
			if not(couple in dico_class_grouped):
				dico_class_grouped[couple] = 0
			dico_class_grouped[couple] += 1
	
	outfile1.write('\t'.join([ACCESS, 'Position_without_missing_data', str(total_unmissing_sites)]))
	outfile1.write('\n')
	for gp in dico_class_grouped:
		outfile1.write('\t'.join([ACCESS, 'Position_affected_to->'+gp, str(dico_class_grouped[gp])]))
		outfile1.write('\n')
	for gp in dico_nb_group:
		outfile1.write('\t'.join([ACCESS, 'Alleles_affected_to->'+gp, str(dico_nb_group[gp])]))
		outfile1.write('\n')
	outfile1.close()
	# It's time to draw
	draw_chr(dico_draw, dico_chr, dico_group, ACCESS, sum(total)/len(total), int(options.ploidy))
	
if __name__ == "__main__": __main__()