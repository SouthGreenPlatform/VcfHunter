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
import os
import sys
import numpy
import optparse
import datetime
import tempfile
import operator
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
from scipy.interpolate import spline
from textwrap import wrap

def couleur(VALUE):
	
	from scipy.stats import norm
	
	v03 = norm(0.05, 0.03)
	v06 = norm(0.09, 0.03)
	v1 = norm(0.16, 0.03)
	
	rouge = min(1, v03.pdf(VALUE)/v03.pdf(0.05))
	verte = min(1, v06.pdf(VALUE)/v06.pdf(0.09))
	bleu = min(1, v1.pdf(VALUE)/v1.pdf(0.16))
	
	if VALUE > 0.16:
		rouge = 1
		verte = 1
		bleu = 1
	if VALUE == 0:
		rouge = 0
		verte = 0
		bleu = 0
		
	
	
	return (rouge, verte, bleu)

def record_loci(loc, chr, DICO_LOCI):
	
	file = open(loc)
	if chr == None:
		for line in file:
			data = line.split()
			if data:
				if not (data[1] in DICO_LOCI):
					DICO_LOCI[data[1]] = []
				DICO_LOCI[data[1]].append([data[0],int(data[2])])
	else:
		list_chr = chr.split(':')
		for line in file:
			data = line.split()
			if data:
				if data[1] in list_chr:
					if not (data[1] in DICO_LOCI):
						DICO_LOCI[data[1]] = []
					DICO_LOCI[data[1]].append([data[0],int(data[2])])
	file.close()

def record_ordered_loci(DICO_LOCI, CHR, AGP, DICO_ORDERED_LOCI, PHYSICAL):
	
	# recording chromosome order
	if CHR == None:
		list_chr = sorted(DICO_LOCI.keys())
	else:
		list_chr = CHR.split(':')
	
	# Calculating mean marker distance
	total = []
	for chr in list_chr:
		list_pos = sorted(DICO_LOCI[chr], key=operator.itemgetter(1))
		if PHYSICAL == 'n':
			marker_order = 1
			for i in range(len(list_pos)):
				list_pos[i][1] = marker_order
				marker_order += 1
		elif PHYSICAL != 'y':
			sys.exit("Wrong value passed to --pysical argument. Only possible values are 'y' or 'n'")
		for i in range(len(list_pos)):
			if i > 0:
				total.append(list_pos[i][1]-list_pos[i-1][1]+1)
	Mean_dist = int(sum(total)/float(len(total)))
	
	# Creating marker order and position
	increment = 0
	if AGP == None:
		for chr in list_chr:
			DICO_ORDERED_LOCI['chr'][chr] = []
			DICO_ORDERED_LOCI['chr'][chr].append(increment)
			DICO_ORDERED_LOCI['marker'].append(chr)
			DICO_ORDERED_LOCI['pos'].append(increment)
			
			list_pos = sorted(DICO_LOCI[chr], key=operator.itemgetter(1))
			if PHYSICAL == 'n':
				marker_order = 1
				for i in range(len(list_pos)):
					list_pos[i][1] = marker_order
					marker_order += 1
			elif PHYSICAL != 'y':
				sys.exit("Wrong value passed to --pysical argument. Only possible values are 'y' or 'n'")
			
			for i in range(len(list_pos)):
				DICO_ORDERED_LOCI['marker'].append(list_pos[i][0])
				DICO_ORDERED_LOCI['pos'].append(list_pos[i][1] + increment)
				last_pos = list_pos[i][1] + increment
			DICO_ORDERED_LOCI['chr'][chr].append(last_pos)
			increment = last_pos + Mean_dist
	else:
		# recording chromosomes size in agp file
		file = open(AGP)
		dico_agp = {}
		for line in file:
			data = line.split()
			if data:
				if data[0] in list_chr:
					if not(data[0] in dico_agp):
						dico_agp[data[0]] = 0
					if dico_agp[data[0]] < int(data[2]):
						dico_agp[data[0]] = int(data[2])
		file.close()
		for chr in list_chr:
			DICO_ORDERED_LOCI['chr'][chr] = []
			DICO_ORDERED_LOCI['chr'][chr].append(increment)
			DICO_ORDERED_LOCI['marker'].append(chr+'-start')
			DICO_ORDERED_LOCI['pos'].append(increment)
			
			list_pos = sorted(DICO_LOCI[chr], key=operator.itemgetter(1))
			if PHYSICAL == 'n':
				marker_order = 1
				for i in range(len(list_pos)):
					list_pos[i][1] = marker_order
					marker_order += 1
			elif PHYSICAL != 'y':
				sys.exit("Wrong value passed to --pysical argument. Only possible values are 'y' or 'n'")
			
			for i in range(len(list_pos)):
				DICO_ORDERED_LOCI['marker'].append(list_pos[i][0])
				DICO_ORDERED_LOCI['pos'].append(list_pos[i][1] + increment)
			last_pos = dico_agp[chr] + increment
			DICO_ORDERED_LOCI['chr'][chr].append(last_pos)
			DICO_ORDERED_LOCI['marker'].append(chr+'-fin')
			DICO_ORDERED_LOCI['pos'].append(last_pos)
			increment = last_pos + 1

def fill_the_matrix(MATRIX, DICO_ORDERED_LOCI, MAT):
	
	file = open(MATRIX)
	header = file.readline().split()
	file.close()
	
	set_mark = set(DICO_ORDERED_LOCI['marker'])
	
	file = open(MATRIX, 'r')
	csvfile = csv.reader(file, delimiter='\t')
	for data in csvfile:
		if data:
			if data[0] in set_mark:
				dico_ordered_loci_index = DICO_ORDERED_LOCI['marker'].index(data[0])
				print(dico_ordered_loci_index)
				for i in range(len(DICO_ORDERED_LOCI['marker'])):
					if DICO_ORDERED_LOCI['marker'][i] in header:
						header_index = header.index(DICO_ORDERED_LOCI['marker'][i])
						value = float(data[header_index])
						MAT[i][dico_ordered_loci_index] = value
						MAT[dico_ordered_loci_index][i] = value
	file.close()

def draw_plot(MAT, DICO_ORDERED_LOCI, AGP, CHR, OUT, STAT):
	
	## Obtaining xmax and ymax
	xmax = DICO_ORDERED_LOCI['pos'][-1]
	ymax = xmax
	
	## done
	
	fig = plt.figure(figsize=(15,15))
	fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
	if STAT == None:
		ax = plt.subplot2grid((30,30),(1,1), colspan=28, rowspan=28)
	else:
		# recording statistics information
		list_pos = []
		list_stat = []
		file = open(STAT)
		for line in file:
			data = line.split()
			if data:
				if data[0] in DICO_ORDERED_LOCI['marker']:
					list_stat.append(float(data[1]))
					list_pos.append(DICO_ORDERED_LOCI['pos'][DICO_ORDERED_LOCI['marker'].index(data[0])])
		file.close()
		stat_max = max(list_stat)*1.1
		
		ax = plt.subplot2grid((30,30),(1,28), colspan=2, rowspan=27)
		if AGP == None:
			ax.set_ylim(0, ymax)
		else:
			falseYmin = -0.028*ymax
			ax.set_ylim(falseYmin, ymax)
		ax.set_xlim(0, stat_max)
		
		ax.axes.get_yaxis().set_visible(False) # remove y ticks
		for label in ax.get_xticklabels():
			label.set_rotation(-90)
			label.set_fontsize(10)
		# ax.axes.get_xaxis().set_visible(False) # remove y ticks
		
		ax.plot(list_stat, list_pos, 'o', mew=0, mfc='black', ms=4)
		
		ax = plt.subplot2grid((30,30),(1,1), colspan=27, rowspan=27)
	
	
	if AGP == None:
		ax.set_xlim(0, xmax)
		ax.set_ylim(0, ymax)
	else:
		falseXmax = xmax+(0.028*xmax)
		falseYmin = -0.028*ymax
		ax.set_xlim(0, falseXmax)
		ax.set_ylim(falseYmin, ymax)
	
	ax.axes.get_yaxis().set_visible(False) # remove y ticks
	ax.axes.get_xaxis().set_visible(False) # remove y ticks
	
	## drawing boxes
	for i in range(len(DICO_ORDERED_LOCI['pos'])):
		for j in range(len(DICO_ORDERED_LOCI['pos'])):
			if i < j:
				posi = DICO_ORDERED_LOCI['pos'][i]
				posj = DICO_ORDERED_LOCI['pos'][j]
				posiplus1 = DICO_ORDERED_LOCI['pos'][i+1]
				posjmoins1 = DICO_ORDERED_LOCI['pos'][j-1]
				val = MAT[i][j]
				if val <= 0.16:
					ax.add_patch(mpatches.Rectangle((posjmoins1,posi) , posj-posjmoins1, posiplus1-posi, fc=couleur(val), ec='black', linewidth=0))
	## done
	
	## Drawing a white triangle for the upper diagonal
	points = [[0, 0], [0, ymax], [xmax, ymax]]
	ax.add_patch(plt.Polygon(points, fc=(1,1,1), ec='red', linewidth=0))
	## done
	
	## drawing chromosome separation
	#	-recording chromosome order
	if CHR == None:
		list_chr = sorted(DICO_ORDERED_LOCI['chr'].keys())
	else:
		list_chr = CHR.split(':')
	#	-drawing separations
	for chr in list_chr:
		ax.add_artist(matplotlib.lines.Line2D((0, xmax), (DICO_ORDERED_LOCI['chr'][chr][1], DICO_ORDERED_LOCI['chr'][chr][1]), color='black', linewidth=1))
		ax.add_artist(matplotlib.lines.Line2D((DICO_ORDERED_LOCI['chr'][chr][1], DICO_ORDERED_LOCI['chr'][chr][1]), (0, ymax), color='black', linewidth=1))
	if AGP != None:
		ax.add_artist(matplotlib.lines.Line2D((0, xmax), (0, 0), color='black', linewidth=1))
	## done
	
	## It's time to draw the scaffolds
	# -recording scaffold position in the AGP
	dico_agp = {}
	if AGP != None:
		file = open(AGP)
		dico_agp = {}
		for line in file:
			data = line.split()
			if data:
				if data[0] in list_chr:
					if data[4] != 'N' and data[4] != 'U':
						if not(data[0] in dico_agp):
							dico_agp[data[0]] = []
						dico_agp[data[0]].append([int(data[1]), int(data[2]), data[8]])
		# -Now we need to draw the shape
		seuil_value = 0.01*ymax
		Ymin = -0.024*ymax
		Ymax = -0.004*ymax
		Xmin = xmax+(0.004*xmax)
		Xmax = xmax+(0.024*xmax)
		for chr in dico_agp:
			for poly in dico_agp[chr]:
				if poly[2] == '+':
					if poly[1] - poly[0] < seuil_value:
						p1 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[0], Ymin]
						p2 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[0], Ymax]
						p3 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[1], (falseYmin/2.0)]
						Xpoints = [p1, p2, p3]
						p1 = [Xmin, DICO_ORDERED_LOCI['chr'][chr][0]+poly[0]]
						p2 = [Xmax, DICO_ORDERED_LOCI['chr'][chr][0]+poly[0]]
						p3 = [((Xmax+Xmin)/2.0), DICO_ORDERED_LOCI['chr'][chr][0]+poly[1]]
						Ypoints = [p1, p2, p3]
					else:
						p1 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[0], Ymin]
						p2 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[0], Ymax]
						p3 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[1]-seuil_value, Ymax]
						p4 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[1], (falseYmin/2.0)]
						p5 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[1]-seuil_value, Ymin]
						Xpoints = [p1, p2, p3, p4, p5]
						p1 = [Xmin, DICO_ORDERED_LOCI['chr'][chr][0]+poly[0]]
						p2 = [Xmax, DICO_ORDERED_LOCI['chr'][chr][0]+poly[0]]
						p3 = [Xmax, DICO_ORDERED_LOCI['chr'][chr][0]+poly[1]-seuil_value]
						p4 = [((Xmax+Xmin)/2.0), DICO_ORDERED_LOCI['chr'][chr][0]+poly[1]]
						p5 = [Xmin, DICO_ORDERED_LOCI['chr'][chr][0]+poly[1]-seuil_value]
						Ypoints = [p1, p2, p3, p4, p5]
				else:
					if poly[1] - poly[0] < seuil_value:
						p1 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[1], Ymin]
						p2 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[1], Ymax]
						p3 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[0], (falseYmin/2.0)]
						Xpoints = [p1, p2, p3]
						p1 = [Xmin, DICO_ORDERED_LOCI['chr'][chr][0]+poly[1]]
						p2 = [Xmax, DICO_ORDERED_LOCI['chr'][chr][0]+poly[1]]
						p3 = [((Xmax+Xmin)/2.0), DICO_ORDERED_LOCI['chr'][chr][0]+poly[0]]
						Ypoints = [p1, p2, p3]
					else:
						p1 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[1], Ymin]
						p2 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[1], Ymax]
						p3 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[0]+seuil_value, Ymax]
						p4 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[0], (falseYmin/2.0)]
						p5 = [DICO_ORDERED_LOCI['chr'][chr][0]+poly[0]+seuil_value, Ymin]
						Xpoints = [p1, p2, p3, p4, p5]
						p1 = [Xmin, DICO_ORDERED_LOCI['chr'][chr][0]+poly[1]]
						p2 = [Xmax, DICO_ORDERED_LOCI['chr'][chr][0]+poly[1]]
						p3 = [Xmax, DICO_ORDERED_LOCI['chr'][chr][0]+poly[0]+seuil_value]
						p4 = [((Xmax+Xmin)/2.0), DICO_ORDERED_LOCI['chr'][chr][0]+poly[0]]
						p5 = [Xmin, DICO_ORDERED_LOCI['chr'][chr][0]+poly[0]+seuil_value]
						Ypoints = [p1, p2, p3, p4, p5]
				ax.add_patch(plt.Polygon(Xpoints, fc=(0.5,0.5,0.5), ec='black', linewidth=1))
				ax.add_patch(plt.Polygon(Ypoints, fc=(0.5,0.5,0.5), ec='black', linewidth=1))
		ax.add_patch(mpatches.Rectangle((0,-0.0032*ymax) , xmax, 0.0032*ymax, fc=(1,1,1), ec='black', linewidth=0))
		ax.add_patch(mpatches.Rectangle((0,falseYmin) , xmax, 0.0032*ymax, fc=(1,1,1), ec='black', linewidth=0))
		ax.add_patch(mpatches.Rectangle((xmax,0) , 0.0032*xmax, ymax, fc=(1,1,1), ec='black', linewidth=0))
		ax.add_patch(mpatches.Rectangle((falseXmax,0) , -0.0032*xmax, ymax, fc=(1,1,1), ec='black', linewidth=0))
		## done
	
	## Drawing chromosome names
	for chr in list_chr:
		start = DICO_ORDERED_LOCI['chr'][chr][0]
		end = DICO_ORDERED_LOCI['chr'][chr][1]
		# print(start, end)
		ax.text(end-(xmax*0.02), end-(xmax*0.004), chr, size=10, va='top', ha='right', fontweight='bold')
	## done
	
	fig.savefig(OUT)
	# fig.savefig(OUT+'.pdf')
	# fig.savefig(OUT+'.svg')
	plt.close(fig)

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr"
	"\n\n This program draw dot plot from a data matrix.")
	# Wrapper options. 
	parser.add_option( '-m', '--matrix',	dest='matrix',	default=None, 	help='The matrix marker file')
	parser.add_option( '-l', '--loc', 		dest='loc',		default=None, 	help='Loci to plot with their locations')
	parser.add_option( '-c', '--chr', 		dest='chr',		default=None, 	help='List of chromosomes to draw in this order (separated by ":")')
	parser.add_option( '-a', '--agp',		dest='agp',		default=None, 	help='Agp file locating scaffolds in the reference sequence')
	parser.add_option( '-s', '--stat',		dest='stat',	default=None, 	help='A two column file with column 1: marker name, column 2: statistics')
	parser.add_option( '-p', '--phys',		dest='phys',	default='y', 	help='A value specifying if the marker position should be defined based on physical position or not. Possible values: "y" or "n", [Default: %default]')
	parser.add_option( '-o', '--output',	dest='output',	default='test', help='Output file name')
	
	(options, args) = parser.parse_args()
	V = os.getpid()
	
	if options.matrix == None:
		sys.exit('Please provide a pairwise matrix file to --matrix argument')
	if options.loc == None:
		sys.exit('Please provide a loci location file to --loc argument')
	
	if not ('.' in options.output):
		sys.exit('Please provide an extension to --output argument. ex: toto.png or toto.svg or toto.pdf')
	
	# 0- Plotting heatmap color legend
	fig = plt.figure(figsize=(15,5))
	fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
	ax = plt.subplot2grid((30,30),(1,1), colspan=28, rowspan=28)
	ax.set_xlim(0, 0.5)
	ax.set_ylim(0, 1)
	ax.axes.get_yaxis().set_visible(False) # remove y ticks
	for i in range(500):
		ax.add_patch(mpatches.Rectangle((i/1000,0) , 1/1000, 1, fc=couleur(i/1000), ec='black', linewidth=0))
	
	fig.savefig('.'.join(options.output.split('.')[0:-1])+'_heatmap.'+options.output.split('.')[-1])
	plt.close(fig)
	
	# 1- Recording chromosomes and loci to draw
	dico_loci = {}
	record_loci(options.loc, options.chr, dico_loci)
	toto = os.system("pmap %s | grep total" % V)
	
	# 2- Recording marker order and their final position
	dico_ordered_loci = {}
	dico_ordered_loci['marker'] = []
	dico_ordered_loci['pos'] = []
	dico_ordered_loci['chr'] = {}
	record_ordered_loci(dico_loci, options.chr, options.agp, dico_ordered_loci, options.phys)
	toto = os.system("pmap %s | grep total" % V)
	
	# 3- It's time to create matrix
	list_4_matrix = list([999999] * len(dico_ordered_loci['marker']))
	mat = []
	for n in range(len(dico_ordered_loci['marker'])):
		mat.append(list(list_4_matrix))
	
	# 3- Now we fill the matrix
	fill_the_matrix(options.matrix, dico_ordered_loci, mat)
	
	# 4- It's time to plot the matrix
	draw_plot(mat, dico_ordered_loci, options.agp, options.chr, options.output, options.stat)
	
	toto = os.system("pmap %s | grep total" % V)

	
if __name__ == "__main__": __main__()