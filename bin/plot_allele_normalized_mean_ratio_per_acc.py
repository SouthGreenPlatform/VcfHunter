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
import os

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


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(
		description="This program draw normalized origin ratio along chromosomes obtained after PaintArp.",
		epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr)")
	# Wrapper options. 
	parser.add_option( '-c',	'--color',			dest='color',		default=None,			help='A color file, five columns: code, complete name, red, green, blue. [Default: %default]')
	parser.add_option( '-C',	'--chr',			dest='chr',			default=None,			help='A chromosome file, two columns: chromosome, size. [Default: %default]')
	parser.add_option( '-a',	'--acc',			dest='acc',			default=None,			help='An accession file, one columns: accession name. [Default: %default]')
	parser.add_option( '-r',	'--ratio',			dest='ratio',		default=None,			help='The ratio files containing on the four first column: chromosome (CHR), median position (POS), start region (Start), end region (End) and after one column for each origin. [Default: %default]')
	parser.add_option( '-l',	'--loc',			dest='loc',			default='',				help='Regions to locate by vertical line. This should be formated this way: Chromosome_name:position,chromosome_name:position, ... [Default: %default]')
	parser.add_option( '-p',	'--prefix',			dest='prefix',		default='Out',			help='Prefix for output files. [Default: %default]')
	parser.add_option( '-g',	'--graph',			dest='graph',		default='png',			help='graphic output. possible options: pdf, png, svg [Default: %default]')
	(options, args) = parser.parse_args()	
	
	if options.color == None:
		sys.exit('Please provide a color file to --color argument')
	if options.chr == None:
		sys.exit('Please provide a chr file to --chr argument')
	if options.acc == None:
		sys.exit('Please provide a acc file to --acc argument')
	
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
	
	file = open(options.color,'r')
	DicoColor = {}
	header = file.readline().split()
	for line in file:
		data = line.split()
		if data:
			DicoColor[data[0]] = (int(data[2])/255.0, int(data[3])/255.0, int(data[4])/255.0)
	file.close()
	
	MAXI = 0
	dicoChrSize = {}
	file = open(options.chr, 'r')
	for line in file:
		data = line.split()
		if data:
			dicoChrSize[data[0]] = int(data[1])
			MAXI = max(MAXI, int(data[1]))
	
	AllCHR = sorted(list(dicoChrSize.keys()))
	
	dicoAcc = {}
	acc = options.acc
	dicoAcc[acc] = {}
	for chr in AllCHR:
		dicoAcc[acc][chr] = [{},[]]
		for n in DicoColor:
			dicoAcc[acc][chr][0][n] = []
	
	
	if options.ratio[-3:] == '.gz':
		file = gzip.open(options.ratio,'rt')
	else:
		file = open(options.ratio, 'r')
	
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
			if chr != data[CHRPos]: # We are on a new chromosome
				if chr: # This is not the first chromosome so we should add the last end window of former chromosomes
					dicoAcc[acc][chr][1].append(END)
					for col in DicoColPos:
						dicoAcc[acc][chr][0][col].append(dicoAcc[acc][chr][0][col][-1])
			
				dicoAcc[acc][data[CHRPos]][1].append(int(data[STARTPos])) # We add the first start window of new chromosome
				for col in DicoColPos:
					value = data[DicoColPos[col]]
					if value == 'NA':
						dicoAcc[acc][data[CHRPos]][0][col].append(np.nan)
					else:
						dicoAcc[acc][data[CHRPos]][0][col].append(float(value))
				chr = data[CHRPos]
				
			else:
				chr = data[CHRPos]
				POS = int(data[POSPos])
				START = int(data[STARTPos])
				END = int(data[ENDPos])
				dicoAcc[acc][chr][1].append(POS)
				for col in DicoColPos:
					value = data[DicoColPos[col]]
					if value == 'NA':
						dicoAcc[acc][chr][0][col].append(np.nan)
					else:
						dicoAcc[acc][chr][0][col].append(float(value))
	dicoAcc[acc][chr][1].append(END)
	for col in DicoColPos:
		dicoAcc[acc][chr][0][col].append(dicoAcc[acc][chr][0][col][-1])
	file.close()
	
	
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
	
	ACCName = sorted(list(dicoAcc.keys()))
	
	NB = len(AllCHR)
	for acc in ACCName:
		
		fig = plt.figure(figsize=(21, 29.7))
		fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
		
		POSSPAN = 0
		
		for chr in AllCHR:
			ax = plt.subplot2grid((NB,30),(POSSPAN,0), colspan=29, rowspan=1)
			if AllCHR.index(chr) == 0:
				ax.set_title(acc, fontweight='bold', position=(0.5, 1.2), size=30)
			ax.set_ylim(-0.05, 1.05)
			ax.set_xlim(0, MAXI)
			
			if chr in VERT:
				for pos in VERT[chr]:
					ax.axvline(x=pos, ymin=-0.05, ymax = 1.05, linewidth=2, color='black')
			
			if chr in VERTREG:
				for pos in VERTREG[chr]:
					ax.fill_betweenx([-0.05,1.05], int(pos[0]), int(pos[1]), color=(0.3,0.3,0.3), alpha=0.20)
			
			for gp  in DicoColor:
				for i in range(9):
					ax.axhline(y=i/8, linewidth=0.1, color = (0,0,0.5))
				
				for i in range(int(MAXI/5000000)+1):
					ax.axvline(x=i*5000000, ymin=-0.05, ymax = 1.05, linewidth=0.1, color = (0,0,0.5))
				
				# Getting position
				ax.plot(dicoAcc[acc][chr][1], dicoAcc[acc][chr][0][gp], color=DicoColor[gp], linewidth=3, label=gp)
				
				if chr != AllCHR[-1]:
					ax.tick_params(axis='x', which='both',  top=False, labeltop=False, labelbottom=False)
				else:
					ax.tick_params(axis='x', which='both',  top=False, labeltop=False)
				ax.tick_params(axis='y', which='both', right=False, labelright=False)
				
				ax.set_xticks(ticks_pos)
				ax.set_xticks(minor_ticks, minor=True)
				
				ax.yaxis.set_ticks(np.arange(0, 1.25, 0.25))
				ax.xaxis.set_tick_params(labelsize=10)
			
			POSSPAN += 1
		ax.set_xticklabels(ticks_labels, fontsize=10)
		
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=min(8, len(DicoColor)), fontsize=25)
		
		# Drawing chromosome name
		POSSPAN = 0
		for chr in AllCHR:
			ax = plt.subplot2grid((NB,30),(POSSPAN,29), colspan=1, rowspan=1, facecolor=(0.75,0.75,0.75))
			ax.axis([0, 1, 0, 1])
			ax.text(0.5, 0.5, chr, size=20, va='center', ha='center', rotation=-90)
			ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False, top=False, labeltop=False)
			ax.tick_params(axis='y', which='both', right=False, labelright=False, left=False, labelleft=False)
			
			POSSPAN += 1
		
		fig.savefig(options.prefix+'.'+options.graph)
		plt.close(fig)

if __name__ == "__main__": __main__()
