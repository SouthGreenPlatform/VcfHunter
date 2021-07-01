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
import configparser
import multiprocessing as mp
from inspect import currentframe, getframeinfo
from operator import itemgetter

sys.stdout.write('modules loaded\n')

def get_color(GCOL, GROUP_TO_DRAW, dico_color):
	if GCOL == None:
		sys.stdout.write('No group color file was provided in --gCol argument. Color will be chosen randomly. If there is more than 7 groups, severals groups will have the same color.\n')
		color = [(0*255,0*255,0*255,0.7),(1*255,0*255,0*255,0.7),(0*255,1*255,0*255,0.7),(0*255,0*255,1*255,0.7),(1*255,1*255,0*255,0.7),(0*255,1*255,1*255,0.7),(1*255,0*255,1*255,0.7)]
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
						dico_color[data[0]] = (dic['red']*255, dic['green']*255, dic['blue']*255, dic['alpha'])
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

def get_chromosome_info(LISTE_ACC, CHR, DICO_HAPLOTYPE):
	
	for acc in LISTE_ACC:
		i = 0
		DICO_HAPLOTYPE[CHR] = 0
		while i < acc[2]:
			i += 1
			file = open(acc[1]+'_'+CHR+'_haplo'+str(i)+'.tab')
			for line in file:
				data = line.split()
				if data:
					DICO_HAPLOTYPE[CHR] = max(int(data[3]), DICO_HAPLOTYPE[CHR])
			file.close()

def GetAcc(ACC, LISTE_ACC):
	file = open(ACC)
	for line in file:
		data = line.split()
		if data:
			LISTE_ACC.append([data[0].split('/')[-1], data[0], int(data[1])])

def convert2tile(FILE, OUT):
	
	file = open(FILE)
	outfile = open(OUT,'a')
	PrecPos = 1
	list2print = []
	for line in file:
		data = line.split()
		if data:
			NewPos = int(data[2])
			if NewPos - PrecPos > 1:
				if list2print[-1][3] == 'color=un':
					list2print[-1][2] = str(NewPos-1)
					list2print.append(data[1:4]+['color='+data[4]])
				elif data[4] == 'un':
					list2print.append([data[1]]+[str(PrecPos+1)]+[data[3]]+['color='+data[4]])
				else:
					list2print.append([data[1]]+[str(PrecPos+1)]+[str(NewPos-1)]+['color=un'])
					list2print.append(data[1:4]+['color='+data[4]])
			else:
				list2print.append(data[1:4]+['color='+data[4]])
			PrecPos = int(data[3])
	file.close()
	for n in list2print:
		outfile.write('\t'.join(n))
		outfile.write('\n')
				
	outfile.close()

def createKar(PREFIX, liste_chr, dico_centro, dico_haplotype):
	outfile = open(PREFIX+'.kar','w')
	for chr in liste_chr:
		outfile.write('\t'.join(['chr','-',chr,chr,'0',str(dico_haplotype[chr]),'white']))
		outfile.write('\n')
	
	for chr in dico_centro:
		if chr in liste_chr:
			outfile.write('\t'.join(['band', chr, 'p36.33', 'p36.33', str(dico_centro[chr][0]), str(dico_centro[chr][1]), 'black']))
			outfile.write('\n')
	outfile.close()

def create_housekeeping(fileName):

	"""
		Create the housekepping.conf file (needed if many point are drawn)

		:param PREFIX: Prefix of the output file.
		:type PREFIX: str
		:return: the housekeeping.conf file
		:rtype: void
	"""

	outfile = open(fileName,'w')
	outfile.write('anglestep       = 0.5\n')
	outfile.write('minslicestep    = 10\n')
	outfile.write('beziersamples   = 40\n')
	outfile.write('debug           = no\n')
	outfile.write('warnings        = no\n')
	outfile.write('imagemap        = no\n')
	outfile.write('paranoid        = yes\n')
	outfile.write('units_ok        = bupr\n')
	outfile.write('units_nounit    = n\n')
	outfile.write('file_delim = \s\n')
	outfile.write('file_delim_collapse = yes\n')
	outfile.write('list_record_delim = \s*[;,]\s*\n')
	outfile.write('list_field_delim  = \s*[:=]\s*\n')
	outfile.write('options_record_delim = [,;]\n')
	outfile.write('options_field_delim  = =\n')
	outfile.write('skip_missing_expression_vars = no\n')
	outfile.write('legacy_underline_expression_syntax = no\n')
	outfile.write('svg_font_scale = 1.3\n')
	outfile.write('sup_baseline_shift = 40\n')
	outfile.write('sub_baseline_shift = -40\n')
	outfile.write('sup_fontsize = 90\n')
	outfile.write('sub_fontsize = 90\n')
	outfile.write('default_font   = default\n')
	outfile.write('default_font_name  = Arial\n')
	outfile.write('default_font_color = black\n')
	outfile.write('default_color  = black\n')
	outfile.write('<guides>\n')
	outfile.write('thickness      = 1\n')
	outfile.write('size           = 5\n')
	outfile.write('type           = outline\n')
	outfile.write('<object>\n')
	outfile.write('all            = no\n')
	outfile.write('ideogram       = no\n')
	outfile.write('ideogram_label = no\n')
	outfile.write('</object>\n')
	outfile.write('<color>\n')
	outfile.write('default = lblue\n')
	outfile.write('text    = red\n')
	outfile.write('</color>\n')
	outfile.write('</guides>\n')
	outfile.write('debug_group = summary,output\n')
	outfile.write('debug_auto_timer_report = 30\n')
	outfile.write('debug_word_separator = " "\n')
	outfile.write('debug_undef_text     = _undef_\n')
	outfile.write('debug_empty_text     = _emptylist_\n')
	outfile.write('debug_validate       = yes\n')
	outfile.write('debug_output_tidy    = no\n')
	outfile.write('text_pixel_subsampling = 1\n')
	outfile.write('text_snuggle_method    = array\n')
	outfile.write('restrict_parameter_names = no\n')
	outfile.write('case_sensitive_parameter_names = no\n')
	outfile.write('calculate_track_statistics = yes\n')
	outfile.write('color_cache_static = yes\n')
	outfile.write('color_cache_file   = circos.colorlist\n')
	outfile.write('color_lists_use    = yes\n')
	outfile.write('memoize = yes\n')
	outfile.write('quit_on_dump = yes\n')
	outfile.write('offsets = 0,0\n')
	outfile.write('max_ticks            = 5000\n')
	outfile.write('max_ideograms        = 200\n')
	outfile.write('max_links            = 1000000\n')
	outfile.write('max_points_per_track = 1000000\n')
	outfile.write('undefined_ideogram = skip\n')
	outfile.write('relative_scale_iterations = 10\n')
	outfile.write('relative_scale_spacing    = mode\n')
	outfile.write('data_out_of_range = trim,warn #\n')
	outfile.write('track_defaults = etc/tracks\n')
	outfile.write('round_brush_use           = yes\n')
	outfile.write('round_brush_min_thickness = 5\n')
	outfile.write('anti_aliasing = yes\n')
	outfile.write('housekeeping = yes\n')
	outfile.write('auto_eval = no\n')
	outfile.close()

def draw_chromosome(ACC, CHR, GCOL, GP, CENTRO, PREFIX):
	
	# recording acc to draw and ploidy level
	liste_acc = []
	GetAcc(ACC, liste_acc)
	
	# recording chromosomes to draw
	liste_chr = CHR.split(':')
	
	# recording groups to work  with
	groups = GP.split(':')
	
	# converting the color file
	dico_color = {}
	get_color(GCOL, groups, dico_color)
	
	# recording centromere position
	dico_centro = {}
	get_centro(CENTRO, dico_centro)
	
	# converting the haplo files to circos tile file
	for acc in liste_acc:
		i = 0
		while i < acc[2]:
			i += 1
			outfile = open(acc[1]+'_haplo'+str(i)+'.tab','w')
			outfile.close()
			for chr in liste_chr:
				convert2tile(acc[1]+'_'+chr+'_haplo'+str(i)+'.tab', acc[1]+'_haplo'+str(i)+'.tab')
	
	# recording haplotypes
	dico_haplotype = {}
	for chr in liste_chr:
		get_chromosome_info(liste_acc, chr, dico_haplotype)
	
	# Creating the karyotype file
	createKar(PREFIX, liste_chr, dico_centro, dico_haplotype)
	
	# Creating the housekeeping.conf file
	create_housekeeping(PREFIX+'_housekeeping.conf')
	
	# Creating the configuration file
	
	RADIUS = 3000
	IDEORAD = 0.90
	
	outfile = open(PREFIX+'.conf','w')
	outfile.write('<<include '+os.getcwd()+'/'+PREFIX+'_housekeeping.conf'+'>>\n')
	outfile.write('<<include etc/colors_fonts_patterns.conf>>\n')
	outfile.write('<colors>\n')
	outfile.write('<<include etc/colors.conf>>\n')
	outfile.write('un = 180,180,180,1\n')
	for col in dico_color:
		outfile.write(' = '.join([col]+[','.join(list(map(str,map(int,dico_color[col]))))]))
		outfile.write('\n')
	outfile.write('</colors>\n')
	outfile.write('karyotype = %s\n' % os.path.realpath(PREFIX+'.kar'))
	outfile.write('<image>\n')
	outfile.write('<<include etc/image.conf>>\n')
	outfile.write('file* = '+PREFIX+'.png\n')
	outfile.write('radius* = '+str(RADIUS)+'p\n')
	outfile.write('angle_offset* = -90\n')
	outfile.write('svg* = no\n')
	outfile.write('</image>\n')
	outfile.write('chromosomes_units = %s\n' % 10000000)
	outfile.write('chromosomes = %s\n' % ';'.join(liste_chr))
	################################################
	outfile.write('<plots>\n')
	nb_sep = len(liste_acc)
	nbhaplo = 0
	MinRadius = 0.75
	for acc in liste_acc:
		nbhaplo += acc[2]
	rP = MinRadius/((4*nbhaplo)+nb_sep)
	rN = 4*rP
	ThicknessValue = (RADIUS-25)*IDEORAD
	pP = ThicknessValue * rP
	pN = ThicknessValue * rN
	print(rP, pP, int(pP))
	print(rN, pN, int(pN))
	
	start = 1
	
	for acc in liste_acc:
		i = 0
		start = start-rP
		while i < acc[2]:
			i += 1
			end = start - rN
			outfile.write('<plot>\n')
			outfile.write('file    = '+acc[1]+'_haplo'+str(i)+'.tab\n')
			outfile.write('show    = yes\n')
			outfile.write('type    = tile\n')
			outfile.write('layers = 15\n')
			outfile.write('margin = 0u\n')
			outfile.write('thickness = '+str(int(pN))+'\n')
			outfile.write('padding = 8\n')
			outfile.write('orientation = out\n')
			outfile.write('stroke_thickness = 1\n')
			outfile.write('stroke_color     = black\n')
			outfile.write('r0    = '+str(end)+'r\n')
			outfile.write('r1    = '+str(start)+'r\n')
			outfile.write('background       = no\n')
			outfile.write('</plot>\n')
			start = end
	outfile.write('</plots>\n')
	
	
	
	################################################
	outfile.write('<ideogram>\n')
	outfile.write('<spacing>\n')
	outfile.write('default = 0.01r\n')
	outfile.write('</spacing>\n')
	outfile.write('thickness        = 50p\n')
	outfile.write('stroke_thickness = 2\n')
	outfile.write('stroke_color     = black\n')
	outfile.write('fill             = yes\n')
	outfile.write('fill_color       = black\n')
	outfile.write('radius         = '+str(IDEORAD)+'r\n')
	outfile.write('show_label     = yes\n')
	outfile.write('label_font     = condensedbold\n')
	outfile.write('label_radius   = dims(ideogram,radius) + 0.05r\n')
	outfile.write('label_size     = 60\n')
	outfile.write('label_parallel = yes\n')
	outfile.write('band_stroke_thickness = 0\n')
	outfile.write('show_bands            = yes\n')
	outfile.write('fill_bands            = yes\n')
	outfile.write('band_transparency     = 1\n')
	outfile.write('</ideogram>\n')


	outfile.write('show_ticks          = yes\n')
	outfile.write('show_tick_labels    = yes\n')
	outfile.write('chrticklabels       = yes\n')
	outfile.write('chrticklabelfont    = default\n')
	outfile.write('grid_start         = dims(ideogram,radius_inner)-0.5r\n')
	outfile.write('grid_end           = dims(ideogram,radius_outer)+100\n')
	outfile.write('<ticks>\n')
	outfile.write('skip_first_label     = no\n')
	outfile.write('skip_last_label      = no\n')
	outfile.write('radius               = dims(ideogram,radius_outer)\n')
	outfile.write('tick_separation      = 2p\n')
	outfile.write('min_label_distance_to_edge = 0p\n')
	outfile.write('label_separation = 5p\n')
	outfile.write('label_offset     = 2p\n')
	outfile.write('label_size = 20p\n')
	outfile.write('multiplier = 1e-7\n')
	outfile.write('color = black\n')
	outfile.write('<tick>\n')
	outfile.write('spacing        = 0.5u\n')
	outfile.write('size           = 5p\n')
	outfile.write('thickness      = 2p\n')
	outfile.write('color          = black\n')
	outfile.write('show_label     = no\n')
	outfile.write('label_size     = 8p\n')
	outfile.write('label_offset   = 0p\n')
	outfile.write('format         = %.0f\n')
	outfile.write('grid           = yes\n')
	outfile.write('grid_color     = grey\n')
	outfile.write('grid_thickness = 1p\n')
	outfile.write('</tick>\n')
	outfile.write('<tick>\n')
	outfile.write('spacing        = 1u\n')
	outfile.write('size           = 8p\n')
	outfile.write('thickness      = 2p\n')
	outfile.write('color          = black\n')
	outfile.write('show_label     = yes\n')
	outfile.write('label_size     = 12p\n')
	outfile.write('label_offset   = 0p\n')
	outfile.write('format         = %.0f\n')
	outfile.write('grid           = yes\n')
	outfile.write('grid_color     = dgrey\n')
	outfile.write('grid_thickness = 1p\n')
	outfile.write('</tick>\n')
	outfile.write('</ticks>\n')
	outfile.close()

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN : guillaume.martin@cirad.fr")
	# Wrapper options. 
	parser.add_option( '',	'--acc',			dest='acc',			default=None,			help='The Accessions names configuration file. Column1: accession name, column2: ploidy level. [Default: %default]')
	parser.add_option( '',	'--chr',			dest='chr',			default=None, 			help='Chromosome(s) to draw, separated by ":". [Default: %default]')
	parser.add_option( '',	'--gcol',			dest='gcol',		default=None, 			help='Groups color. [Default: %default]')
	parser.add_option( '',	'--dg',				dest='dg',			default=None, 			help='Groups to draw, separated by ":". [Default: %default]')
	parser.add_option( '',	'--centro',			dest='centro',		default=None, 			help='(Peri)centromere positions. [Default: %default]')
	parser.add_option( '',	'--prefix',			dest='prefix',		default='All_acc', 		help='Prefix for output file. [Default: %default]')

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
	draw_chromosome(options.acc, options.chr, options.gcol, options.dg, options.centro, options.prefix)
	
	pathname = os.path.dirname(sys.argv[0])
	loca_programs = configparser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')
	CIRCOS = loca_programs.get('Programs','circos')
	
	qs=os.popen(CIRCOS+' -conf '+options.prefix+'.conf -noparanoid')
	
	for n in qs:
		sys.stdout.write(n.replace("\n","\r\n"))
	
if __name__ == "__main__": __main__()