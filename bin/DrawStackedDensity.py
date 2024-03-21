#!/usr/bin/env python
#
#  Copyright 2019 CIRAD
# Guillaume Martin Franc-Christophe Baurens
#
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
#V0 = all genome by chr or by genome
#V1 = graphical output options added figure files are  named as 'genome'_'fasta'.'graphical format'

# -*- coding: utf-8 -*-

import optparse, os, shutil, subprocess, sys, tempfile, fileinput, operator, time, math, gzip
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#fonction de dessin ref1 en rouge
def draw_plot_per_chr(ListID, Y, X, OUT, CHR, FILLUNDER, YLIM, COLOR, NEGATIVE, VERT, VERTREG, FIGURESIZE): 
	
	# obtaining figure size
	FigureSize = list(map(float, FIGURESIZE.split(',')))
	
	fig = plt.figure(figsize=FigureSize)
	ax = plt.subplot2grid((1,1),(0,0))
	
	if NEGATIVE == 'y':
		for ID in reversed(ListID):
			ax.plot(X[CHR], Y[ID][CHR], linewidth=0.5, color=COLOR[ID], label=ID)
			ax.plot(X[CHR], np.array(Y[ID][CHR])*-1, linewidth=0.5, color=COLOR[ID])
			if FILLUNDER == 'y':
				ax.fill_between(X[CHR],np.array(Y[ID][CHR])*-1, np.array(Y[ID][CHR]), color=COLOR[ID])
	else:
		for ID in reversed(ListID):
			ax.plot(X[CHR], Y[ID][CHR], linewidth=0.5, color=COLOR[ID], label=ID)
			if FILLUNDER == 'y':
				ax.fill_between(X[CHR],0, np.array(Y[ID][CHR]), color=COLOR[ID])
	
	if NEGATIVE == 'y':
		if YLIM != None:
			ax.set_ylim(float(YLIM*-1), float(YLIM))
		else:
			YLIM = max(np.array(Y[ListID[-1]][CHR]))
			ax.set_ylim(YLIM*-1, YLIM)
	else:
		if YLIM != None:
			ax.set_ylim(0, float(YLIM))
		else:
			YLIM = max(np.array(Y[ListID[-1]][CHR]))
			ax.set_ylim(0, YLIM)
	
	if CHR in VERT:
		for pos in VERT[CHR]:
			ax.axvline(x=pos, ymin=0, ymax = float(YLIM), linewidth=2, color='black')
	
	if CHR in VERTREG:
		for pos in VERTREG[CHR]:
			ax.fill_betweenx([0,float(YLIM)], int(pos[0]), int(pos[1]), color=(0.3,0.3,0.3), alpha=0.20)
	
	ax.set_xlim(0, X[CHR][-1])
	ax.set_title(CHR, fontweight='bold', position=(0.5, 1))
	plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.03), fancybox=True, shadow=True, ncol=len(COLOR), fontsize=20)
	
	fig.savefig(OUT)
	plt.close(fig)

def draw_plot(DICOCHR, LISTCHR, ListID, Y, X, OUT, FILLUNDER, YLIM, COLOR, NEGATIVE, SCALE, VERT, VERTREG, FIGURESIZE):
	
	# obtaining figure size
	FigureSize = list(map(float, FIGURESIZE.split(',')))
	
	# Calculating chromosome number
	NB = len(LISTCHR)
	
	# Calculating max chromosome length
	XLIM = 0
	for chr in LISTCHR:
		XLIM = max(XLIM, max(X[chr]))
	
	# Calculating max Y value if not fixed
	if YLIM == None:
		YLIM = 0
		for chr in LISTCHR:
			for ID in ListID:
				YLIM = max(YLIM, max(Y[ID][chr]))
	else:
		YLIM = float(YLIM)
	
	# Drawing figures
	POSSPAN = 0
	fig = plt.figure(figsize=FigureSize)
	fig.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)
	
	# Drawing graph
	for CHR in LISTCHR:
		ax = plt.subplot2grid((NB+1,15),(POSSPAN,0), colspan=14, rowspan=1)
		
		if NEGATIVE == 'y':
			for ID in reversed(ListID):
				ax.plot(X[CHR], Y[ID][CHR], linewidth=0.5, color=COLOR[ID], label=ID)
				ax.plot(X[CHR], np.array(Y[ID][CHR])*-1, linewidth=0.5, color=COLOR[ID])
				if FILLUNDER == 'y':
					ax.fill_between(X[CHR],np.array(Y[ID][CHR])*-1, np.array(Y[ID][CHR]), color=COLOR[ID])
		else:
			for ID in reversed(ListID):
				ax.plot(X[CHR], Y[ID][CHR], linewidth=0.5, color=COLOR[ID], label=ID)
				if FILLUNDER == 'y':
					ax.fill_between(X[CHR],0, np.array(Y[ID][CHR]), color=COLOR[ID])
		
		if NEGATIVE == 'y':
			if YLIM != None:
				ax.set_ylim(float(YLIM*-1), float(YLIM))
			else:
				ax.set_ylim(max(np.array(Y[ListID[-1]][CHR]))*-1, max(np.array(Y[ListID[-1]][CHR])))
		else:
			if YLIM != None:
				ax.set_ylim(0, float(YLIM))
			else:
				ax.set_ylim(0, max(np.array(Y[ListID[-1]][CHR])))
		
		if 's' in SCALE:
			ax.set_xlim(0, XLIM)
		else:
			ax.set_xlim(0, X[CHR][-1])
		ax.axvline(x=DICOCHR[CHR],linewidth=2, color='black')
		
		if CHR in VERT:
			for pos in VERT[CHR]:
				ax.axvline(x=pos, ymin=0, ymax = float(YLIM), linewidth=2, color='black')
		
		if CHR in VERTREG:
			for pos in VERTREG[CHR]:
				ax.fill_betweenx([0,float(YLIM)], int(pos[0]), int(pos[1]), color=(0.3,0.3,0.3), alpha=0.20)
		
		POSSPAN += 1
	plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=len(COLOR), fontsize=20)
	
	# Drawing chromosome name
	POSSPAN = 0
	for chr in LISTCHR:
		ax = plt.subplot2grid((NB+1,15),(POSSPAN,14), colspan=1, rowspan=1)
		ax.axis('off')	
		ax.axis([0, 1, 0, 1])
		ax.text(0, 0.5, chr, size=12, va='center', ha='left', fontweight='bold')
		
		POSSPAN += 1
	
	fig.savefig(OUT)
	plt.close(fig)

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser()
	# Wrapper options. 
	parser.add_option( '-f', '--files',		dest='files',		default=None,			help='List of density files separated by ",". These files should be formated as follows: ID1:File1,ID2:File2,...')
	parser.add_option( '-F', '--FillUnder',	dest='FillUnder',	default='y',			help='Fill the curve under, (y or n)')
	parser.add_option( '-c', '--chrToRm',	dest='chrToRm',		default='',				help='List of chromosome to exclude from drawing separated with ","')
	parser.add_option( '-C', '--color',		dest='color',		default=None,			help='A color file in tabulated format (in RGB). Col1: RefxName, Col2: Red, Col3: Green, Col4: Blue')
	parser.add_option( '-y', '--Ylim',		dest='Ylim',		default=None,			help='Set Y limit for all plots. Recommended to compare genomes with the same graphic scale --> if omitted each graph may have different Yscale (adjusted from each chromosome)')
	parser.add_option( '-d', '--draw',		dest='draw',		default='i',			help='Drawing output. Possible options: One figure per genome (g), one figure per chromosomes (i), scaled image (s) - only with (g), in this case chromosomes are drawn scaled. Options can be combined')
	parser.add_option( '-n', '--negative',	dest='negative',	default='n',			help='Draw also negative curve: y or n')
	parser.add_option( '-l',  '--loc',		dest='loc',			default='',				help='Regions to locate by vertical line. This should be formated this way: Chromosome_name,position:chromosome_name,position: ... [Default: %default]')
	parser.add_option( '-g', '--graph',		dest='graph',		default='svg',			help='graphic output. possible options: pdf, png, svg')
	parser.add_option( '-s',  '--figsize',	dest='figsize',		default='21,29.7',		help='Figure size (in inches). [Default: %default]')
	parser.add_option( '-p', '--prefix',	dest='prefix',		default='Stacked',		help='Prefix for the output file(s)')


	(options, args) = parser.parse_args()
	
	# Recording regions to locate
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
	
	# Recording chromosomes to exclude
	ChrToExclude=options.chrToRm.split(',')
	
	# Recording chromosome information
	FilesList = options.files.split(',')
	dicoInfo={}
	ListID = []
	DicoPos = {}
	listchr = set()
	chrSize = {}
	for File in FilesList:
		ID, filename = File.split(":")
		ListID.append(ID)
		dicoInfo[ID] = {}
		if filename[-3:] == '.gz':
			file = gzip.open(filename, 'rt')
		else:
			file = open(filename, 'r')
		for line in file:
			data = line.split()
			if data:
				chr, debut, fin, value = data[0:4]
				debut = int(debut)
				fin = int(fin)
				pos = (debut+fin)/2
				if not (chr in ChrToExclude):
					if not(chr in chrSize):
						chrSize[chr] = 0
					chrSize[chr] = max(chrSize[chr], fin)
					listchr.add(chr)
					if len(ListID) == 1:
						if not (chr in dicoInfo[ID]):
							dicoInfo[ID][chr] = []
							DicoPos[chr] = []
						DicoPos[chr].append(pos)
						dicoInfo[ID][chr].append(float(value))
					else:
						if not (chr in dicoInfo[ID]):
							dicoInfo[ID][chr] = []
						precID = ListID[-2]
						if not(pos in DicoPos[chr]):
							sys.exit("There is a problem, windows are not the same in files\n")
						else:
							VALUE = float(value) + dicoInfo[precID][chr][DicoPos[chr].index(pos)]
							dicoInfo[ID][chr].append(float(VALUE))
		file.close()
	
	# Adding values at beginning and end of chromosomes
	for ID in ListID:
		for chr in dicoInfo[ID]:
			dicoInfo[ID][chr].insert(0, dicoInfo[ID][chr][0])
			dicoInfo[ID][chr].append(dicoInfo[ID][chr][-1])
	for chr in DicoPos:
		DicoPos[chr].insert(0, 1)
		DicoPos[chr].append(chrSize[chr])
	
	# Obtaining colors if any
	DicoCol = {}
	if options.color == None:
		cmap = matplotlib.cm.get_cmap('gist_rainbow')
		for i in range(len(ListID)):
			if i%2 == 0:
				DicoCol[ListID[i]] = cmap((i/2)/(len(ListID)))
			else:
				DicoCol[ListID[i]] = cmap((((i-1)/2)+math.ceil(len(ListID)/2))/(len(ListID)))
	else:
		file = open(options.color)
		for line in file:
			data = line.split()
			if data:
				DicoCol[data[0]] = (int(data[1])/255, int(data[2])/255, int(data[3])/255)
		file.close()
		for ref in ListID:
			if not (ref in DicoCol):
				print(ref+' was not found in color file, an arbitrary color was attributed')
				DicoCol[ref] = (238/255, 16/255, 16/255)
	
	# Drawing figures
	listchr = sorted(listchr)
	if 'g' in options.draw:
		draw_plot(chrSize, listchr, ListID, dicoInfo, DicoPos, options.prefix+"."+options.graph, options.FillUnder, options.Ylim, DicoCol, options.negative, options.draw, VERT, VERTREG, options.figsize)
	if 'i' in options.draw:
		for chr in listchr:
			draw_plot_per_chr(ListID, dicoInfo, DicoPos, options.prefix+'_'+chr+"."+options.graph, chr, options.FillUnder, options.Ylim, DicoCol, options.negative, VERT, VERTREG, options.figsize)


if __name__ == "__main__": __main__()


