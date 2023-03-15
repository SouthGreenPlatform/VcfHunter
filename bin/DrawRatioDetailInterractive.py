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
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
from scipy.interpolate import spline
from textwrap import wrap

def draw_chr(DICOCURVE, DICO_INFO, CHR_INFO, DICO_GROUP, PLOIDY, DCURVE, halfWindow, PSIZE, LSIZE, VERT, VERTREG, DICOCOLOR, FIGSIZE, CHR, START, END, VARIATION):
	
	# obtaining figure size
	FigureSize = list(map(float, FIGSIZE.split(',')))
	
	# Calculating group number plus 1
	NB = len(DICO_GROUP) + 1
	
	# Getting group order
	group = sorted(list(DICO_GROUP))
	
	# Calculating Ymax
	INTERVAL = VARIATION[1]+0.2
	YMAX = (len(group)+1)*INTERVAL
	
	
	######Drawing allele ratio######
	# Drawing a first graph with all groups
	fig,ax = plt.subplots(figsize=FigureSize)
	fig.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)
	ax.set_ylim(0, YMAX)
	ax.set_xlim(START, END)
	ax.axes.yaxis.set_ticklabels([])
	
	# Drawing horizontal lines separating each origin
	for i in range(len(group)+1):
		ax.axhline(y=i*INTERVAL, linewidth=1,  color='black')
	
	# Adding chromosome name
	for i in range(len(group)):
		ax.text(START-(END-START)*0.01, (i*INTERVAL)+0.5, group[i], size=15, va='center', ha='right', fontweight='bold')
	ax.text(START-(END-START)*0.01, (len(group)*INTERVAL)+0.5, 'All', size=15, va='center', ha='right', fontweight='bold')
	
	# Adding vertical regions if any
	if CHR in VERT:
		for pos in VERT[CHR]:
			ax.axvline(x=pos, ymin=0, ymax = YMAX, linewidth=2, color='black')
	if CHR in VERTREG:
		for pos in VERTREG[CHR]:
			ax.fill_betweenx([0,YMAX], int(pos[0]), int(pos[1]), color=(0.3,0.3,0.3), alpha=0.20)
	
	# Curve drawing if any
	if CHR in DICOCURVE:
		for i in range(len(group)):
			gp = group[i]
			if gp in DICOCURVE[CHR]:
				ax.plot(DICOCURVE[CHR]['Position'], np.array(DICOCURVE[CHR][gp])+(i*INTERVAL), color=DICOCOLOR[gp], lw=LSIZE)
				if len(DICOCURVE[CHR]['Position']) != len(DICOCURVE[CHR][gp]):
					sys.exit('There is a bug...\n')
	
	# Preparing for dot drawing
	XVal = []
	YVal = []
	Color = []
	names = []
	for i in range(len(group)):
		gp = group[i]
		position = sorted(list(DICO_INFO[CHR].keys()))
		for pos in position:
			if gp in DICO_INFO[CHR][pos]:
				YVal.append(DICO_INFO[CHR][pos][gp]+(i*INTERVAL))
				XVal.append(pos)
				names.append(str(pos))
				Color.append(DICOCOLOR[gp])
	sc = plt.scatter(XVal, YVal, c=Color, s=PSIZE, edgecolors='face')
	
	# Initiating annotation
	annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points", bbox=dict(boxstyle="round", fc="w"), arrowprops=dict(arrowstyle="->"))
	annot.set_visible(False)
	
	
	def update_annot(ind):
		pos = sc.get_offsets()[ind["ind"][0]]
		annot.xy = pos
		text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
		print(text)
		annot.set_text(text)

	def hover(event):
		vis = annot.get_visible()
		if event.inaxes == ax:
			cont, ind = sc.contains(event)
			if cont:
				update_annot(ind)
				annot.set_visible(True)
				fig.canvas.draw_idle()
			else:
				if vis:
					annot.set_visible(False)
					fig.canvas.draw_idle()
	
	fig.canvas.mpl_connect("button_press_event", hover)


	
	plt.show()
	sys.exit()
		
		
		
	
	
	
	# fig, ax = plt.subplot2grid((1,15),(0,1), colspan=14, rowspan=1)
	# ax.set_ylim(0, YMAX)
	# ax.set_xlim(START, END)
	# if CHR in VERT:
		# for pos in VERT[CHR]:
			# ax.axvline(x=pos, ymin=0, ymax = 1.1, linewidth=2, color='black')
	# if CHR in VERTREG:
		# for pos in VERTREG[CHR]:
			# ax.fill_betweenx([0,1.1], int(pos[0]), int(pos[1]), color=(0.3,0.3,0.3), alpha=0.20)
	
	# Drawing vertical line separating each origin
	# for i in range(len(group)+1):
		# ax.axvline(x=pos, ymin=0, ymax = 1.1, linewidth=2, color='black')
	
	
	for gp in group:
		# Getting position
		position = sorted(list(DICO_INFO[CHR].keys()))
		value =  []
		final_pos = []
		ax.axhline(y=0.25, linewidth=0.1, color = (0,0,0.5))
		ax.axhline(y=0.75, linewidth=0.1, color = (0,0,0.5))
		ax.axhline(y=0.5, linewidth=0.1, color = (0,0,0.5))
		ax.axhline(y=0.33, linewidth=0.1, color = (0.5,0.5,0.5))
		ax.axhline(y=0.66, linewidth=0.1, color = (0.5,0.5,0.5))
		for pos in position:
			if gp in DICO_INFO[CHR][pos]:
				value.append(DICO_INFO[CHR][pos][gp])
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
	
	POSSPAN += 1
	
	# Drawing per group
	for gp in group:
		ax = plt.subplot2grid((NB,15),(POSSPAN,1), colspan=14, rowspan=1)
		ax.set_ylim(0, 1.1)
		ax.set_xlim(START, END)
		
		if CHR in VERT:
			for pos in VERT[CHR]:
				ax.axvline(x=pos, ymin=0, ymax = 1.1, linewidth=2, color='black')
		
		if CHR in VERTREG:
			for pos in VERTREG[CHR]:
				ax.fill_betweenx([0,1.1], int(pos[0]), int(pos[1]), color=(0.3,0.3,0.3), alpha=0.20)
		
		# Getting position
		position = sorted(list(DICO_INFO[CHR].keys()))
		value =  []
		final_pos = []
		ax.axhline(y=0.25, linewidth=0.1, color = (0,0,0.5))
		ax.axhline(y=0.75, linewidth=0.1, color = (0,0,0.5))
		ax.axhline(y=0.5, linewidth=0.1, color = (0,0,0.5))
		ax.axhline(y=0.33, linewidth=0.1, color = (0.5,0.5,0.5))
		ax.axhline(y=0.66, linewidth=0.1, color = (0.5,0.5,0.5))
		for pos in position:
			if gp in DICO_INFO[CHR][pos]:
				value.append(DICO_INFO[CHR][pos][gp])
				final_pos.append(pos)
		
		if CHR in DICOCURVE:
			if gp in DICOCURVE[CHR]:
				ax.plot(DICOCURVE[CHR]['Position'], DICOCURVE[CHR][gp], color=DICOCOLOR[gp], lw=LSIZE)
				if len(DICOCURVE[CHR]['Position']) != len(DICOCURVE[CHR][gp]):
					sys.exit('There is a bug...\n')
		
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
	
	# Drawing chromosome name
	POSSPAN = 1
	for gp in group:
		ax = plt.subplot2grid((NB,15),(POSSPAN,0), colspan=1, rowspan=1)
		ax.axis('off')	
		ax.axis([0, 1, 0, 1])
		ax.text(0, 0.5, gp, size=20, va='center', fontweight='bold')
		POSSPAN += 1
	
	fig.canvas.mpl_connect('button_press_event', onclick)
	plt.show()
	# fig.savefig(OUT+'Ratio.png')
	# plt.close(fig)

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr\n This program draw allele ratios obtained by allele_ratio_per_acc.py program")
	# Wrapper options. 
	parser.add_option( '',	'--chr',			dest='chr',			default=None,			help='Path to a file containing chromosome information. Two columns are required: col1 -> chromosome, col2 -> size.')
	parser.add_option( '',	'--origin',			dest='origin',		default=None,			help='Either a 5 column file containing chromosome (col1), position (Col2), allele (Col3), origin (Col4) and allele ratio (Col 5) or the output file of allele_ratio_per_acc.py.')
	parser.add_option( '',	'--NormOri',		dest='NormOri',		default=None,			help='A 4+n column file containing chromosome name chromosome (CHR), median position (POS), start region (Start), end region (End) and after one column for each origin.')
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
	parser.add_option( '',	'--reg',			dest='reg',			default=None,			help='Region to draw. This could be a chromosome or a chromosome region. To specify a chromosome, just put the name. If it is a region, it should be formated as follows: name,start,end. [Default: %default]')
	# parser.add_option( '',	'--prefix',			dest='prefix',		default='',				help='Prefix for output files. Not required [Default: %default]')
	(options, args) = parser.parse_args()	
	
	if options.chr == None:
		sys.exit('Please provide a chr file to --chr argument')
	if options.origin == None:
		sys.exit('Please provide a origin file to --origin argument')
	if options.acc == None:
		sys.exit('Please provide a accession name to --acc argument')
	if options.ploidy == None:
		sys.exit('Please provide a ploidy level to --ploidy argument')
	if options.reg == None:
		sys.exit('Please provide a region to draw to --reg argument')
	
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
	# PREFIX = options.prefix
	
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
				dico_group.add(data[0])
		file.close()
	
	# recording chromosome information
	dico_chr = {}
	file = open(options.chr)
	for line in file:
		data = line.split()
		if data:
			dico_chr[data[0]] = int(data[1])
	file.close()
	
	# Identification of region to draw
	Region = options.reg.split(',')
	if len(Region) == 1:
		CHR = Region[0]
		START = 0
		END = dico_chr[CHR]
	elif len(Region) == 3:
		CHR = Region[0]
		START = int(Region[1])
		END = int(Region[2])
	else:
		sys.exit("Oups, the information passed to --reg option is not understood... The program exited without finishing\n")
	
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
	if options.NormOri != None:
		if options.NormOri[-3:] == '.gz':
			file = gzip.open(options.NormOri,'rt')
		else:
			file = open(options.NormOri, 'r')
		
		header = file.readline().split()
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
				if data[CHRPos] == CHR: # We are on a new chromosome
					if not (CHR in DicoCurve): # Management of first position
						DicoCurve[CHR] = {}
						DicoCurve[CHR]['Position'] = [int(data[STARTPos])]
						for col in DicoColor:
							value = data[DicoColPos[col]]
							if value == 'NA':
								DicoCurve[CHR][col] = [np.nan]
							else:
								DicoCurve[CHR][col] = [float(value)]
					# Management of "normal" line
					DicoCurve[CHR]['Position'].append(int(data[POSPos]))
					for col in DicoColor:
						value = data[DicoColPos[col]]
						if value == 'NA':
							DicoCurve[CHR][col].append(np.nan)
						else:
							DicoCurve[CHR][col].append(float(value))
					ENDCurve = int(data[ENDPos]) # For management of last line
		# Management of last line
		DicoCurve[CHR]['Position'].append(ENDCurve)
		for col in DicoColor:
			DicoCurve[CHR][col].append(DicoCurve[CHR][col][-1])
		file.close()
	
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
	draw_chr(DicoCurve, dico_draw, dico_chr, dico_group, int(options.ploidy), options.dcurve, int(options.halfwin), float(options.psize), float(options.lsize), VERT, VERTREG, DicoColor, options.figsize, CHR, START, END, VARIATION)
	
if __name__ == "__main__": __main__()