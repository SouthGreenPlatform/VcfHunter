#!/usr/bin/env python
#
#  Copyright 2017 CIRAD
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

sys.stdout.write('modules loaded\n')

def get_color(GCOL, GROUP_TO_DRAW, dico_color):
	if GCOL == None:
		sys.stdout.write('No group color file was provided in --gCol argument. Color will be chosen randomly. If there is more than 7 groups, severals groups will have the same color.\n')
		color = [(0,0,0,0.7),(1,0,0,0.7),(0,1,0,0.7),(0,0,1,0.7),(1,1,0,0.7),(0,1,1,0.7),(1,0,1,0.7)]
		if len(GROUP_TO_DRAW) > len(color):
			sys.stdout.write('In the graphe some groups will have the same color because there is not enough color defined in the script.\nPlease contact Guillaume MARTIN (guillaume.martin@cirad.fr) for adding colors.\n')
		for i in range(len(GROUP_TO_DRAW)):
			dico_color[GROUP_TO_DRAW[i]] = color[i%len(color)]
			if int(i/len(color)) > 0:
				sys.stdout.write('In the graphe groups '+GROUP_TO_DRAW[i]+' will have the same color as the '+str(GROUP_TO_DRAW[i%len(color)])+' group.\n')
	else:
		file = open(GCOL)
		sur_group = False
		sur_color = False
		for line in file:
			data = line.split()
			if data:
				if data[0] == '[group]':
					sur_group = True
					sur_color = False
				elif data[0] == '[color]':
					sur_group = False
					sur_color = True
				elif sur_color:
					if data[0] in GROUP_TO_DRAW:
						color_list = data[1].split(':')
						dic = {}
						for k in color_list:
							dic[k.split('=')[0]] = float(k.split('=')[1])
						if not('red' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: red color proportion is missing')
						if not('green' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: green color proportion is missing')
						if not('blue' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: blue color proportion is missing')
						if not('alpha' in dic):
							sys.exit('There is a problem in the file passed to --gCol argument: alpha color proportion is missing')
						dico_color[data[0]] = (dic['red'], dic['green'], dic['blue'], dic['alpha'])
		# checking all groups have a color
		for n in GROUP_TO_DRAW:
			if not(n in dico_color):
				sys.exit('There is a problem. The group '+n+' have no color. the program exited without finishing')

def get_centro(CENTRO, dico_centro):
	
	file = open(CENTRO)
	for line in file:
		data = line.split()
		if data:
			dico_centro[data[0]] = [int(data[1]), int(data[2])]

def draw_plot(DICO, GCOL, GROUPS, ACC, LIST_CHR, CENTRO, PLOIDY, GRAPH):
	
	# recording centromere position
	dico_centro = {}
	get_centro(CENTRO, dico_centro)
	
	# recording color to draw
	dico_color = {}
	get_color(GCOL, GROUPS, dico_color)
	
	# Getting total chromosome number
	total_chr = 0
	total_gp = 0
	max_length = 0
	for chr in DICO:
		if chr in LIST_CHR:
			total_gp += 1
			for hap in DICO[chr]:
				if hap != 'length':
					total_chr += 1
			max_length = max(max_length, DICO[chr]['length'])
	
	# Calculating plot metrics
	short_space = 0.4
	large_space = 0.6
	ymax = (total_chr*1) + (total_chr*large_space/PLOIDY) + (total_chr*short_space/PLOIDY*(PLOIDY-1)) + large_space
	
	
	# Calculating ticks
	i = 5000000
	ticks_pos = [0]
	ticks_labels = [0]
	
	while i < max_length:
		ticks_pos.append(i)
		ticks_labels.append(str(int(i/1000000)))
		i += 5000000
	i = 1000000
	
	minor_ticks = []
	while i < max_length:
		minor_ticks.append(i)
		i += 1000000
	
	# Drawing plot
	fig = plt.figure(figsize=(14.85,10.5))
	fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.1)
	ax = plt.subplot2grid((1,15),(0,1), colspan=14, rowspan=1)
	
	ax.set_frame_on(False) # remove axis
	ax.get_xaxis().tick_bottom() # only draw bottom ticks
	ax.axes.get_yaxis().set_visible(False) # remove y ticks
	ax.add_artist(matplotlib.lines.Line2D((0, max_length), (0, 0), color='black', linewidth=2)) # on ajoute l'axe horizontal
	ax.set_xticks(ticks_pos)
	ax.set_xticks(minor_ticks, minor=True)
	ax.set_xticklabels(ticks_labels)
	
	ax.set_title(ACC, fontweight='bold', position=(0.5, 1.005), size=20)
	ax.set_ylim(0, ymax)
	ax.set_xlim(-max_length*0.01, max_length*1.01)
	
	# Drawing legend
	for gp in GROUPS:
		ax.plot([],[], color=dico_color[gp][0:3], label=gp, linewidth=10)
	ax.plot([],[], color=(0.6,0.6,0.6), label='unknown', linewidth=10)
	ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=len(GROUPS)+1, fontsize=15)
	
	# Drawing chromosomes shape
	y0 = ymax-short_space
	for chr in LIST_CHR:
		for hap in DICO[chr]:
			if hap != 'length':
				y0 -= 1
				x0 = 1
				x_size = DICO[chr]['length']
				ax.add_patch(mpatches.Rectangle((x0,y0) , x_size, 1, ec='black', fc='black', linewidth=3))
				y0 -= short_space
		y0 -= (large_space - short_space)
	
	# Locating pericentromeric regions
	y0 = ymax-short_space
	for chr in LIST_CHR:
		for hap in DICO[chr]:
			if hap != 'length':
				y0 -= 1
				x0 = dico_centro[chr][0]
				x_size = dico_centro[chr][1] - x0
				ax.add_patch(mpatches.Rectangle((x0,y0) , x_size, 1, ec='grey', fc='black', linewidth=3))
				y0 -= short_space
		y0 -= (large_space - short_space)
	
	# to fill blancs by grey
	y0 = ymax-short_space
	for chr in LIST_CHR:
		for hap in DICO[chr]:
			if hap != 'length':
				y0 -= 1
				x0 = 1
				x_size = DICO[chr]['length']
				ax.add_patch(mpatches.Rectangle((x0,y0) , x_size, 1, ec=(0.6,0.6,0.6), fc=(0.6,0.6,0.6)))
				y0 -= short_space
		y0 -= (large_space - short_space)
	
	y0 = ymax-short_space
	# Drawing chromosomes
	for chr in LIST_CHR:
		for hap in DICO[chr]:
			if hap != 'length':
				y0 -= 1
				for n in range(len(DICO[chr][hap]['pos'])):
					x0 = DICO[chr][hap]['pos'][n][0]
					x_size = (DICO[chr][hap]['pos'][n][1] - DICO[chr][hap]['pos'][n][0]) + 1
					gp = DICO[chr][hap]['gp'][n]
					if gp in dico_color:
						col = dico_color[gp][0:3]
					else:
						col = (0.6,0.6,0.6)
					ax.add_patch(mpatches.Rectangle((x0,y0), x_size, 1, color=col))
				y0 -= short_space
		y0 -= (large_space - short_space)
	
	
	# Adding chromosome names
	ax = plt.subplot2grid((1,15),(0,0), colspan=1, rowspan=1)
	ax.axis('off')	
	ax.axis([0, 1, 0, ymax])
	y0 = ymax-short_space
	for chr in LIST_CHR:
		for hap in DICO[chr]:
			if hap != 'length':
				y0 -= 1
				ax.text(1, y0+0.5, "\n".join(wrap(chr, 20)), size=10, va='center', ha='right', fontweight='bold')
				y0 -= short_space
		y0 -= (large_space - short_space)
	
	fig.savefig(ACC+'.'+GRAPH)
	plt.close(fig)
	sys.exit()

def get_haplotype(ACC, LIST_CHR, dico_haplotype, ploidy):
	
	i = 0
	while i < ploidy:
		i += 1
		for chr in LIST_CHR:
			file = open(ACC+'_'+chr+'_haplo'+str(i)+'.tab')
			if not(chr in dico_haplotype):
				dico_haplotype[chr] = {}
				dico_haplotype[chr]['length'] = 0
			dico_haplotype[chr][i] = {}
			dico_haplotype[chr][i]['pos'] = []
			dico_haplotype[chr][i]['gp'] = []
			for line in file:
				data = line.split()
				if data:
					dico_haplotype[chr][i]['pos'].append(list(map(int, [data[2], data[3]])))
					dico_haplotype[chr][i]['gp'].append(data[4])
					dico_haplotype[chr]['length'] = int(data[3])
			file.close()

def draw_chromosome(ACC, CHR, GCOL, GP, CENTRO, PLOIDY, GRAPH):
	
	# recording chromosome to draw
	list_chr = CHR.split(':')
	
	# recording groups to work  with
	groups = GP.split(':')
	
	# recording haplotypes
	dico_haplotype = {}
	get_haplotype(ACC, list_chr, dico_haplotype, PLOIDY)
	
	draw_plot(dico_haplotype, GCOL, groups, ACC, list_chr, CENTRO, PLOIDY, GRAPH)

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '',	'--acc',			dest='acc',			default=None,			help='The Accession name. [Default: %default]')
	parser.add_option( '',	'--chr',			dest='chr',			default=None, 			help='Chromosomes to draw, separated by ":". [Default: %default]')
	parser.add_option( '',	'--gcol',			dest='gcol',		default=None, 			help='Groups color. [Default: %default]')
	parser.add_option( '',	'--dg',				dest='dg',			default=None, 			help='Groups to draw, separated by ":". [Default: %default]')
	parser.add_option( '',	'--centro',			dest='centro',		default=None, 			help='(Peri)centromere positions. [Default: %default]')
	parser.add_option( '',	'--ploidy',			dest='ploidy',		default='2', 			help='Ploidy level. Possible values: 2, 3 [Default: %default]')
	parser.add_option( '',  '--graph',			dest='graph',		default='pdf',			help='graphic output. possible options: pdf, png, svg [Default: %default]')

	(options, args) = parser.parse_args()
	
	
	# Identify genome blocs in accessions and draw a circos by accessions
	if options.acc == None:
		sys.exit('Please provide a acc name to --acc argument')
	if options.chr == None:
		sys.exit('Please provide a chr list to --chr argument')
	if options.gcol == None:
		sys.exit('Please provide a gcol file to --gcol argument')
	if options.dg == None:
		sys.exit('Please provide a dg string to --dg argument')
	if options.centro == None:
		sys.exit('Please provide a centro file to --centro argument')
	draw_chromosome(options.acc, options.chr, options.gcol, options.dg, options.centro, int(options.ploidy), options.graph)
		
if __name__ == "__main__": __main__()