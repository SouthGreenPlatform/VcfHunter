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
import random
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


def draw_chr(DICONUM, DICOCURVE, DICO_INFO, CHR_INFO, DICO_GROUP, OUT, PLOIDY, DCURVE, halfWindow, PSIZE, LSIZE, VERT, VERTREG, DICOCOLOR, FIGSIZE):
	
	# obtaining figure size
	FigureSize = list(map(float, FIGSIZE.split(',')))
	
	# getting chromosomes list
	chr_list_temp = sorted(list(CHR_INFO.keys()))
	chr_list = []
	for chr in chr_list_temp:
		if chr in DICO_INFO:
			chr_list.append(chr)
	
	# Calculating chromosome number
	if DICONUM:
		NB = len(chr_list)*2
		YNUMMAX = 0
		for chr in DICONUM:
			for gp in DICONUM[chr][0]:
				YNUMMAX = max(max(DICONUM[chr][0][gp]),YNUMMAX)
	else:
		NB = len(chr_list)
	
	# Calculating max chromosome size
	MAXI = 0
	for n in CHR_INFO:
		if MAXI < CHR_INFO[n]:
			MAXI = CHR_INFO[n]
	
	# Getting group order
	group = sorted(list(DICO_GROUP))
	
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
	fig = plt.figure(figsize=FigureSize)
	fig.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)
	
	# Drawing per chromosome
	for chr in chr_list:
		ax = plt.subplot2grid((NB,15),(POSSPAN,1), colspan=14, rowspan=1)
		ax.set_ylim(0, 1.1)
		ax.set_xlim(0, MAXI)
		
		if chr in VERT:
			for pos in VERT[chr]:
				ax.axvline(x=pos, ymin=0, ymax = 1.1, linewidth=2, color='black')
		
		if chr in VERTREG:
			for pos in VERTREG[chr]:
				ax.fill_betweenx([0,1.1], int(pos[0]), int(pos[1]), color=(0.3,0.3,0.3), alpha=0.20)
		
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
			
			if chr in DICOCURVE:
				if gp in DICOCURVE[chr][0]:
					ax.plot(DICOCURVE[chr][1], DICOCURVE[chr][0][gp], color=DICOCOLOR[gp], lw=LSIZE)
					if len(DICOCURVE[chr][1]) != len(DICOCURVE[chr][0][gp]):
						sys.exit('There is a bug...\n')
					
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
		
		if DICONUM:
			POSSPAN += 1
			ax = plt.subplot2grid((NB,15),(POSSPAN,1), colspan=14, rowspan=1)
			ax.set_ylim(0, YNUMMAX)
			ax.set_xlim(0, MAXI)
			# ax.axes.yaxis.set_ticklabels([])
			ax.axes.xaxis.set_ticklabels([])
			if chr in DICONUM:
				for gp in DICONUM[chr][0]:
					if chr in VERT:
						for pos in VERT[chr]:
							ax.axvline(x=pos, ymin=0, ymax = YNUMMAX, linewidth=2, color='black')
					
					if chr in VERTREG:
						for pos in VERTREG[chr]:
							ax.fill_betweenx([0,YNUMMAX], int(pos[0]), int(pos[1]), color=(0.3,0.3,0.3), alpha=0.20)
					
					if gp in DICOCOLOR:
						ax.plot(DICONUM[chr][1], DICONUM[chr][0][gp], color=DICOCOLOR[gp], lw=LSIZE)
					else:
						ax.plot(DICONUM[chr][1], DICONUM[chr][0][gp], lw=LSIZE)
		
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
		if DICONUM:
			POSSPAN += 1
			ax = plt.subplot2grid((NB,15),(POSSPAN,0), colspan=1, rowspan=1)
			ax.axis('off')	
			ax.axis([0, 1, 0, 1])
			ax.text(0, 0.5, chr, size=12, va='center', fontweight='bold')
		
		POSSPAN += 1
	
	fig.savefig(OUT+'Ratio.png')
	plt.close(fig)

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n This program draw allele ratios obtained by allele_ratio_per_acc.py program")
	# Wrapper options. 
	parser.add_option( '',	'--chr',			dest='chr',			default=None,			help='Path to single vcf file on vcf for all chromosomes is available. [Default: %default]')
	parser.add_option( '',	'--origin',			dest='origin',		default=None,			help='A 5 column file containing chromosome (col1), position (Col2), allele (Col3), origin (Col4) and allele ratio (Col 5). [Default: %default]')
	parser.add_option( '',	'--NormOri',		dest='NormOri',		default=None,			help='A 4+n column file containing chromosome name chromosome (CHR), median position (POS), start region (Start), end region (End) and after one column for each origin.')
	parser.add_option( '',	'--Unknown',		dest='Unknown',		default=None,			help='A 4 column file containing chromosome (col1), position (Col2), allele (Col3), origin (Col4). This file will be used to draw density along chromosomes.')
	parser.add_option( '',	'--acc',			dest='acc',			default=None,			help='Accession to work with. [Default: %default]')
	parser.add_option( '',	'--ploidy',			dest='ploidy',		default=None,			help='Accession ploidy')
	parser.add_option( '',  '--dcurve',			dest='dcurve',		default='n',			help='Draw mean curve for ratio. Possible values: y or n [Default: %default]')
	parser.add_option( '',  '--psize',			dest='psize',		default='1.5',			help='Dot size in graph. [Default: %default]')
	parser.add_option( '',  '--lsize',			dest='lsize',		default='1',			help='Size of the line of the mean value curve. [Default: %default]')
	parser.add_option( '',  '--win',			dest='halfwin',		default='10',			help='Size of half sliding window that allow to draw mean value curve [Default: %default]')
	parser.add_option( '',  '--loc',			dest='loc',			default='',				help='Regions to locate by vertical line. This should be formated this way: Chromosome_name,position:chromosome_name,position: ... [Default: %default]')
	parser.add_option( '',  '--col',			dest='col',			default=None,			help='A color file with 4 columns: col1=group and the three last column corresponded to RGB code. [Default: %default]')
	parser.add_option( '',  '--figsize',		dest='figsize',		default='10.5,14.85',	help='Figure size (in inches). [Default: %default]')
	parser.add_option( '',  '--MinMax',			dest='MinMax',		default='1,1',			help='Additional variations added to values equal to 1 to have a better idea of dot density. Do not exceed 1.1 in maximal value. Example: 0.9,1.1 [Default: %default]')
	parser.add_option( '',	'--prefix',			dest='prefix',		default='',				help='Prefix for output files. Not required [Default: %default]')
	(options, args) = parser.parse_args()	
	
	if options.chr == None:
		sys.exit('Please provide a chr file to --chr argument')
	if options.origin == None:
		sys.exit('Please provide a origin file to --origin argument')
	if options.acc == None:
		sys.exit('Please provide a accession name to --acc argument')
	if options.ploidy == None:
		sys.exit('Please provide a ploidy level to --ploidy argument')
	
	VARIATION = list(map(float, options.MinMax.split(',')))
	
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
	
	# recording accession name
	ACCESS = options.acc
	PREFIX = options.prefix
	
	# recording allele group
	dico_origin = {}
	dico_group = set()
	if options.origin[-3:] == '.gz':
		file = gzip.open(options.origin, 'rt')
	else:
		file = open(options.origin, 'r')
	header = False
	for line in file:
		data = line.split()
		if data:
			if data == ['chr', 'pos', 'allele', 'obs_ratio', 'exp_ratio', 'grp']:
				header = True
			else:
				if header:
					chr = data[0]
					pos = int(data[1])
					allele = data[2]
					Origin = data[5]
				else:
					chr = data[0]
					pos = int(data[1])
					allele = data[2]
					Origin = data[3]
				if not(chr in dico_origin):
					dico_origin[chr] = {}
				if not(pos in dico_origin[chr]):
					dico_origin[chr][pos] = {}
				if allele in dico_origin[chr][pos]:
					sys.exit('Oups, there is a problem, the allele '+allele+' at position '+str(pos)+' on chromosome '+chr+' has more than one origin... this is not allowed\n')
				dico_origin[chr][pos][allele] = Origin
				dico_group.add(Origin)
	file.close()
	
	# Preparing color file
	DicoColor = {}
	if options.col == None:
		cmap = matplotlib.cm.get_cmap('gist_rainbow')
		groupList = sorted(list(dico_group))
		for i in range(len(groupList)):
			DicoColor[groupList[i]] = cmap(i/(len(groupList)-1))
	else:
		file = open(options.col)
		for line in file:
			data = line.split()
			if data:
				DicoColor[data[0]] = (int(data[1])/255.0,int(data[2])/255.0,int(data[3])/255.0)
		file.close()
	
	# recording chromosome information
	dico_chr = {}
	file = open(options.chr)
	for line in file:
		data = line.split()
		if data:
			dico_chr[data[0]] = int(data[1])
	file.close()
	
	# recording allele information
	dico_draw = {}
	if options.origin[-3:] == '.gz':
		file = gzip.open(options.origin, 'rt')
	else:
		file = open(options.origin, 'r')
	header = False
	for line in file:
		data = line.split()
		if data:
			if data == ['chr', 'pos', 'allele', 'obs_ratio', 'exp_ratio', 'grp']:
				header = True
			else:
				if header:
					chr = data[0]
					pos = int(data[1])
					ori = data[5]
					value = float(data[3])
				else:
					chr = data[0]
					pos = int(data[1])
					ori = data[3]
					value = float(data[4])
				if value == 1:
					value = random.uniform(VARIATION[0], VARIATION[1])
				if not (chr in dico_draw):
					dico_draw[chr] = {}
				if not(int(pos) in dico_draw[chr]):
					dico_draw[chr][pos] = {}
				dico_draw[chr][pos][ori] = value
	file.close()
	
	# Recording curve to draw if any
	DicoCurve = {}
	DicoCurve[ACCESS] = {}
	DicoWin = {}
	DicoWin[ACCESS] = {}
	if options.NormOri != None:
		if options.NormOri[-3:] == '.gz':
			file = gzip.open(options.NormOri,'rt')
		else:
			file = open(options.NormOri, 'r')
		for chr in dico_chr:
			DicoCurve[ACCESS][chr] = [{},[]]
			DicoWin[ACCESS][chr] = []
			for n in DicoColor:
				DicoCurve[ACCESS][chr][0][n] = []
		
		header = file.readline().split()
		chr = ""
		CHRPos = header.index("CHR")
		POSPos = header.index("POS")
		STARTPos = header.index("Start")
		ENDPos = header.index("End")
		DicoColPos = {}
		for col in DicoColor:
			DicoColPos[col] = header.index(col)
		for line in file:
			data = line.split()
			if data:
				DicoWin[ACCESS][data[CHRPos]].append((int(data[STARTPos]), int(data[ENDPos])))
				if chr != data[CHRPos]: # We are on a new chromosome
					if chr: # This is not the first chromosome so we should add the last end window of former chromosomes
						DicoCurve[ACCESS][chr][1].append(END)
						for col in DicoColPos:
							DicoCurve[ACCESS][chr][0][col].append(DicoCurve[ACCESS][chr][0][col][-1])
				
					DicoCurve[ACCESS][data[CHRPos]][1].append(int(data[STARTPos])) # We add the first start window of new chromosome
					for col in DicoColPos:
						value = data[DicoColPos[col]]
						if value == "NA":
							DicoCurve[ACCESS][data[CHRPos]][0][col].append(np.nan)
						else:
							DicoCurve[ACCESS][data[CHRPos]][0][col].append(float(value))
					chr = data[CHRPos]
					
				else:
					chr = data[CHRPos]
					POS = int(data[POSPos])
					START = int(data[STARTPos])
					END = int(data[ENDPos])
					DicoCurve[ACCESS][chr][1].append(POS)
					for col in DicoColPos:
						value = data[DicoColPos[col]]
						if value == 'NA':
							DicoCurve[ACCESS][chr][0][col].append(np.nan)
						else:
							DicoCurve[ACCESS][chr][0][col].append(float(value))
		file.close()
		
		DicoCurve[ACCESS][chr][1].append(END)
		for col in DicoColPos:
			DicoCurve[ACCESS][chr][0][col].append(DicoCurve[ACCESS][chr][0][col][-1])
	
	# Recording allele number per windows if file provided
	DicoNum = {}
	GPset = set()
	ChrSet = set()
	if options.Unknown != None:
		if options.Unknown[-3:] == '.gz':
			file = gzip.open(options.Unknown,'rt')
		else:
			file = open(options.Unknown, 'r')
		for line in file:
			data = line.split()
			ChrSet.add(data[0])
			GPset.add(data[3])
		file.close()
		
		for chr in ChrSet:
			for gp in GPset:
				ListValues = [0]*dico_chr[chr]
				if options.Unknown[-3:] == '.gz':
					file = gzip.open(options.Unknown,'rt')
				else:
					file = open(options.Unknown, 'r')
				for line in file:
					data = line.split()
					if data:
						if data[0] == chr:
							if data[3] == gp:
								ListValues[int(data[1])] += 1
				file.close()
				# Filling the Values
				if not(chr in DicoNum):
					DicoNum[chr] = [{}, []]
					for i in range(len(DicoWin[ACCESS][chr])):
						if i == 0:
							DicoNum[chr][1].append(DicoWin[ACCESS][chr][i][0])
						DicoNum[chr][1].append(np.mean(DicoWin[ACCESS][chr][i]))
					DicoNum[chr][1].append(DicoWin[ACCESS][chr][i][1])
				DicoNum[chr][0][gp] = []
				for i in range(len(DicoWin[ACCESS][chr])):
					start, end = DicoWin[ACCESS][chr][i]
					if i == 0:
						DicoNum[chr][0][gp].append(sum(ListValues[start-1:end]))
					DicoNum[chr][0][gp].append(sum(ListValues[start-1:end]))
				DicoNum[chr][0][gp].append(sum(ListValues[start-1:end]))
	
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
	
	# It's time to draw
	draw_chr(DicoNum, DicoCurve[ACCESS], dico_draw, dico_chr, dico_group, PREFIX+ACCESS, int(options.ploidy), options.dcurve, int(options.halfwin), float(options.psize), float(options.lsize), VERT, VERTREG, DicoColor, options.figsize)
	
if __name__ == "__main__": __main__()