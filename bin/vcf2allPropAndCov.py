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
import gzip

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
from scipy.interpolate import spline
from textwrap import wrap


def draw_chr(DICO_INFO, CHR_INFO, DICO_GROUP, OUT, MEAN_COV, PLOIDY, DCURVE, halfWindow, PSIZE, LSIZE, VERT, VERTREG, DICOCOLOR):
	
	# getting chromosomes list
	chr_list_temp = sorted(list(CHR_INFO.keys()))
	chr_list = []
	for chr in chr_list_temp:
		if chr in DICO_INFO:
			chr_list.append(chr)
	
	# Calculating chromosome number
	NB = len(chr_list)
	
	# Calculating max chromosome size
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
	
	# Drawing figures
	POSSPAN = 0
	fig = plt.figure(figsize=(10.5, 14.85))
	fig.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)
	
	# Drawing per chromosome
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,1), colspan=14, rowspan=1)
		ax.set_ylim(0, 1.1)
		ax.set_xlim(0, MAXI)
		
		if chr in VERT:
			for pos in VERT[chr]:
				ax.axvline(x=pos, ymin=0, ymax = 2*MEAN_COV, linewidth=2, color='black')
		
		if chr in VERTREG:
			for pos in VERTREG[chr]:
				ax.fill_betweenx([0,2*MEAN_COV], int(pos[0]), int(pos[1]), color=(0.3,0.3,0.3), alpha=0.20)
		
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
			
			ax.plot(final_pos, value, 'o', ms=PSIZE, mew=0, mfc=DICOCOLOR[gp])
			
			if DCURVE == 'y':
				ValueNorm = []
				PosNorm = []
				Window = 1+2*halfWindow
				for n in range(len(value)-Window):
					IntermVal = value[n:n+Window]
					ValueNorm.append(sum(IntermVal)/len(IntermVal))
					PosNorm.append(final_pos[n+halfWindow])
				ax.plot(PosNorm, ValueNorm, color=DICOCOLOR[gp], lw=LSIZE)
			
			ax.axes.yaxis.set_ticklabels([])
			ax.axes.xaxis.set_ticklabels([])
		
		POSSPAN += 1
	ax.set_xticks(ticks_pos)
	ax.set_xticks(minor_ticks, minor=True)
	ax.set_xticklabels(ticks_labels)
	
	# Drawing chromosome name
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
	
	# Creating figure
	POSSPAN = 0
	fig = plt.figure(figsize=(10.5, 14.85))
	fig.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)
	
	# Drawing per chromosomes
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
		
		if chr in VERT:
			for pos in VERT[chr]:
				ax.axvline(x=pos, ymin=0, ymax = 2*MEAN_COV, linewidth=2, color='black')
		
		if chr in VERTREG:
			for pos in VERTREG[chr]:
				ax.fill_betweenx([0,2*MEAN_COV], int(pos[0]), int(pos[1]), color=(0.3,0.3,0.3), alpha=0.20)
		
		for pos in position:
			value.append(DICO_INFO[chr][pos]['cov'])
			final_pos.append(pos)
		ax.fill_between([0]+final_pos+[CHR_INFO[chr]], fill_plus1, fill_moins1, color=(1,0,0), alpha=0.20) ## For band
		ax.axhline(y=MEAN_COV, linewidth=0.5, color = (0,0,0))
		ax.axhline(y=(PLOIDY+1)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,1,0))
		ax.axhline(y=(PLOIDY-1)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,1,0))
		ax.axhline(y=(PLOIDY+2)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,0,1))
		ax.axhline(y=(PLOIDY-2)/PLOIDY*MEAN_COV, linewidth=0.5, color = (0,0,1))
		ax.plot(final_pos, value, 'o', ms=PSIZE, mew=0, mfc='black')
		ax.axes.yaxis.set_ticklabels([])
		ax.axes.xaxis.set_ticklabels([])
		
		POSSPAN += 1
	ax.set_xticks(ticks_pos)
	ax.set_xticks(minor_ticks, minor=True)
	ax.set_xticklabels(ticks_labels)
	
	# Adding chromosome name
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
	parser.add_option( '',	'--vcf',			dest='vcf',			default=None,			help='Path to single vcf file on vcf for all chromosomes is available. [Default: %default]')
	parser.add_option( '',	'--origin',			dest='origin',		default=None,			help='A 2 column file containing accession name (col1) origin/group (Col2). [Default: %default]')
	parser.add_option( '',	'--acc',			dest='acc',			default=None,			help='Accession to work with. [Default: %default]')
	parser.add_option( '',	'--ploidy',			dest='ploidy',		default=None,			help='Accession ploidy')
	parser.add_option( '',	'--NoMiss',			dest='NoMiss',		default='n',			help='No missing data are allowed in accessions used to group alleles. [Default: %default]')
	parser.add_option( '',	'--all',			dest='all',			default='n',			help='Allele should be present in all accessions of the group. [Default: %default]')
	parser.add_option( '',  '--dcurve',			dest='dcurve',		default='n',			help='Draw mean curve for ratio. Possible values: y or n [Default: %default]')
	parser.add_option( '',  '--psize',			dest='psize',		default='1.5',			help='Dot size in graph. [Default: %default]')
	parser.add_option( '',  '--lsize',			dest='lsize',		default='1',			help='Size of the line of the mean value curve. [Default: %default]')
	parser.add_option( '',  '--win',			dest='halfwin',		default='10',			help='Size of half sliding window that allow to draw mean value curve [Default: %default]')
	parser.add_option( '',  '--loc',			dest='loc',			default='',				help='Regions to locate by vertical line. This should be formated this way: Chromosome_name:position,chromosome_name:position, ... [Default: %default]')
	parser.add_option( '',  '--col',			dest='col',			default=None,			help='A color file with 4 columns: col1=group and the three last column corresponded to RGB code. [Default: %default]')
	parser.add_option( '',  '--excl',			dest='excl',		default=None,			help='A file containing region to exclude for ancestry attribution in some introgressed accessions. [Default: %default]')
	parser.add_option( '',	'--prefix',			dest='prefix',		default='',				help='Prefix for output files. Not required [Default: %default]')
	(options, args) = parser.parse_args()	
	
	if options.conf == None and options.vcf == None:
		sys.exit('Please provide either a conf file to --conf argument or a vcf file to --vcf argument')
	if options.origin == None:
		sys.exit('Please provide a origin file to --origin argument')
	if options.acc == None:
		sys.exit('Please provide a accession name to --acc argument')
	if options.ploidy == None:
		sys.exit('Please provide a ploidy level to --ploidy argument')
	
	VERT = {}
	VERTREG = {}
	for n in options.loc.split(':'):
		info = n.split(',')
		if len(info) == 1:
			pass
		else:
			chr = info[0]
			if not(chr in VERT):
				VERT[chr] = set()
				VERTREG[chr] = []
			if len(info) == 2:
				VERT[chr].add(int(info[1]))
			elif len(info) == 3:
				VERTREG[chr].append(info[1:])
			else:
				sys.stdout.write('Wrong format passed to --loc argument\n')
	
	# Recording regions to exclude if required
	sys.stdout.write("Recording excluded regions\n")
	sys.stdout.flush()
	DicoRegion = {}
	if options.excl != None:
		file = open(options.excl)
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
	
	# recording accession name
	ACCESS = options.acc
	PREFIX = options.prefix
	
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
	
	outfile = open(PREFIX+ACCESS+'_AlleleOriginAndRatio.tab','w')
	outfile1 = open(PREFIX+ACCESS+'_stats.tab','w')
	
	# Preparing color file
	DicoColor = {}
	if options.col == None:
		cmap = matplotlib.cm.get_cmap('gist_rainbow')
		groupList = sorted(list(dico_group.keys()))
		for i in range(len(groupList)):
			DicoColor[groupList[i]] = cmap(i/(len(groupList)-1))
	else:
		file = open(options.col)
		for line in file:
			data = line.split()
			if data:
				DicoColor[data[0]] = (int(data[1])/255.0,int(data[2])/255.0,int(data[3])/255.0)
	
	
	# recording vcf file to work with
	dico_vcf = set()
	if options.vcf == None:
		file = open(options.conf)
		for line in file:
			data = line.split()
			if data:
				dico_vcf.add(data[0])
		file.close()
	else:
		dico_vcf.add(options.vcf)
	
	total_unmissing_sites = 0
	
	
	sys.stdout.write("Working on the vcf\n")
	sys.stdout.flush()
	
	# For linear drawing
	dico_draw = {}
	dico_chr = {}
	total = []
	for vcf in dico_vcf:
		if vcf[-3:] == '.gz':
			file = gzip.open(vcf,'rt')
		else:
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
					
					# Recording non introgressed accessions at the position
					dicoIntro = set()
					for gp in dico_group:
						for ACC in dico_group[gp]:
							if ACC in DicoRegion:
								if CHR in DicoRegion[ACC]:
									for region in DicoRegion[ACC][CHR]:
										START = region[0]
										END = region[1]
										if START < POS and POS < END:
											dicoIntro.add(ACC)
					
					# recording group alleles
					dico_allele = {}
					for gp in dico_group:
						dico_allele[gp] = set()
						for ACC in dico_group[gp]:
							if not (ACC in dicoIntro):
								ACCESSION = data[header.index(ACC)].split(':')
								GENOTYPE = set(ACCESSION[FORMAT.index("GT")].replace('|','/').split('/'))
								for geno in GENOTYPE:
									dico_allele[gp].add(geno)
					
					# to manage missing data in groups
					Good = True
					for gp in dico_allele:
						if options.NoMiss == 'n':
							if '.' in dico_allele[gp] and len(dico_allele[gp]) == 1:
								Good = False
						elif options.NoMiss == 'y':
							if '.' in dico_allele[gp]:
								Good = False
						else:
							sys.exit('Oups, their is a bug... either the vcf has a problem or I made a mistake in the programing\n')
						if len(dico_allele[gp]) == 0:
							Good = False
							print('Unpredicted region because of introgression in all representative of a group :', CHR, POS, gp)
					if Good:
						# recording accession
						add_cov = 0
						if not ('.' in GENOTYPE):
							total_unmissing_sites += 1
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
									if not (ACC in dicoIntro):
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
						
						ACCESSION = data[header.index(ACCESS)].split(':')
						GENOTYPE = set(ACCESSION[FORMAT.index("GT")].replace('|','/').split('/'))
						COVERAGE = list(map(int, ACCESSION[FORMAT.index("AD")].split(',')))
						for geno in GENOTYPE:
							if geno != '.':
								group = []
								for gp in dico_allele:
									if geno in dico_allele[gp]:
										group.append(gp)
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
	draw_chr(dico_draw, dico_chr, dico_group, PREFIX+ACCESS, sum(total)/len(total), int(options.ploidy), options.dcurve, int(options.halfwin), float(options.psize), float(options.lsize), VERT, VERTREG, DicoColor)
	
if __name__ == "__main__": __main__()