#!/usr/bin/env python
#
#  Copyright 2023 CIRAD
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

import subprocess
import argparse
import tempfile
import fnmatch
import gzip
import sys
import os

from inspect import currentframe, getframeinfo
from Bio import SeqIO

def run_job (frameinfo, cmd_line, ERROR):
	try:
		tmp = tempfile.NamedTemporaryFile().name
		error = open(tmp, 'w')
		proc = subprocess.Popen( args=cmd_line, shell=True, stderr=error)
		returncode = proc.wait()
		error.close()
		error = open( tmp, 'r' )
		stderr = ''
		buffsize = 1048576
		try:
			while True:
				stderr += error.read( buffsize )
				if not stderr or len( stderr ) % buffsize != 0:
					break
		except OverflowError:
			pass
		error.close()
		os.remove(tmp)
		if returncode != 0:
			raise Exception
		return 0
	except Exception:
		return 'Line : '+str(frameinfo.lineno)+' - '+ERROR + str( stderr )

def record_color(FILE):
	DiCol = {}
	file  = open(FILE)
	header = file.readline().split()
	for line in file:
		data = line.split()
		if data:
			DiCol[data[header.index('group')]] = [data[header.index('name')], data[header.index('r')], data[header.index('g')], data[header.index('b')]]
	file.close()
	return DiCol

def format_haplo(CHILD, FOLDER, COLOR, CROSS, IND, PLO, CHR):
	for i in range(PLO):
		Indice = str(i+1)
		outfile = open(CHILD+'/'+CROSS+'.circos.'+IND+'.haplo'+str(Indice)+'.tab','w')
		for chro in CHR:
			fileName = FOLDER+'/'+IND+'_'+chro+'_haplo'+str(Indice)+'.tab'
			if os.path.exists(fileName):
				file = open(fileName, 'r')
				for line in file:
					data = line.split()
					if data[4] != 'un':
						data[4] = 'color='+data[4]
						outfile.write('\t'.join((data[1:])))
						outfile.write('\n')
			else:
				print('Warning: the file '+fileName+' was not found, this is not necessarily a problem, but if you expect to have painting on chromosome '+chro+' this means that there is a problem.')
		outfile.close()

def create_housekeeping(CN, CROSS):
	outfile = open(CN+'/'+CROSS+'_housekeeping.conf','w')
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
	outfile.write('max_links            = 100000000\n')
	outfile.write('max_points_per_track = 100000000\n')
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

def add_ideogram_section(FILE):
	FILE.write('<ideogram>\n')
	FILE.write('<spacing>\n')
	FILE.write('default = 0.01r\n')
	FILE.write('</spacing>\n')
	FILE.write('thickness        = 50p\n')
	FILE.write('stroke_thickness = 2\n')
	FILE.write('stroke_color     = black\n')
	FILE.write('fill             = yes\n')
	FILE.write('fill_color       = black\n')
	FILE.write('radius         = 0.9r\n')
	FILE.write('show_label     = yes\n')
	FILE.write('label_font     = condensedbold\n')
	FILE.write('label_radius   = dims(ideogram,radius) + 0.05r\n')
	FILE.write('label_size     = 90\n')
	FILE.write('label_parallel = yes\n')
	FILE.write('band_stroke_thickness = 0\n')
	FILE.write('show_bands            = yes\n')
	FILE.write('fill_bands            = yes\n')
	FILE.write('band_transparency     = 1\n')
	FILE.write('</ideogram>\n')
	FILE.write('show_ticks          = yes\n')
	FILE.write('show_tick_labels    = yes\n')
	FILE.write('chrticklabels       = yes\n')
	FILE.write('chrticklabelfont    = default\n')
	FILE.write('grid_start         = dims(ideogram,radius_inner)-0.5r\n')
	FILE.write('grid_end           = dims(ideogram,radius_outer)+100\n')
	FILE.write('<ticks>\n')
	FILE.write('skip_first_label     = no\n')
	FILE.write('skip_last_label      = no\n')
	FILE.write('radius               = dims(ideogram,radius_outer)\n')
	FILE.write('tick_separation      = 2p\n')
	FILE.write('min_label_distance_to_edge = 0p\n')
	FILE.write('label_separation = 5p\n')
	FILE.write('label_offset     = 2p\n')
	FILE.write('label_size = 20p\n')
	FILE.write('multiplier = 1e-6\n')
	FILE.write('color = black\n')
	FILE.write('<tick>\n')
	FILE.write('spacing        = 1u\n')
	FILE.write('size           = 10p\n')
	FILE.write('thickness      = 3p\n')
	FILE.write('color          = black\n')
	FILE.write('show_label     = no\n')
	FILE.write('label_size     = 8p\n')
	FILE.write('label_offset   = 0p\n')
	FILE.write('format         = %.0f\n')
	FILE.write('grid           = yes\n')
	FILE.write('grid_color     = grey\n')
	FILE.write('grid_thickness = 1p\n')
	FILE.write('</tick>\n')
	FILE.write('<tick>\n')
	FILE.write('spacing        = 10u\n')
	FILE.write('size           = 20p\n')
	FILE.write('thickness      = 5p\n')
	FILE.write('color          = black\n')
	FILE.write('show_label     = yes\n')
	FILE.write('label_size     = 35p\n')
	FILE.write('label_offset   = 5p\n')
	FILE.write('format         = %.0f\n')
	FILE.write('grid           = yes\n')
	FILE.write('grid_color     = dgrey\n')
	FILE.write('grid_thickness = 1p\n')
	FILE.write('</tick>\n')
	FILE.write('</ticks>\n')

def add_start_file(OUT, KAR, HOUSE, CHR_LIST, DICO_COLOR, FILE):
	wkdir = os.getcwd()
	FILE.write('<<include '+wkdir+'/'+HOUSE+'>>\n')
	FILE.write('<<include etc/colors_fonts_patterns.conf>>\n')
	FILE.write('<colors>\n')
	FILE.write('<<include etc/colors.conf>>\n')
	for gp in DICO_COLOR:
		FILE.write(gp+'='+','.join(DICO_COLOR[gp][1:]+['1'])+'\n')
	FILE.write('un=125,125,125,1\n')
	FILE.write('vvgrey=230,230,230,1\n')
	FILE.write('</colors>\n')
	FILE.write('karyotype = '+wkdir+'/'+KAR+'\n')
	FILE.write('<image>\n')
	FILE.write('<<include etc/image.conf>>\n')
	FILE.write('file* = '+OUT+'\n')
	FILE.write('radius* = 3000p\n')
	FILE.write('angle_offset* = -90\n')
	FILE.write('svg* = yes\n')
	FILE.write('background* = transparent\n')
	FILE.write('</image>\n')
	FILE.write('chromosomes_units = 1000000\n')
	FILE.write('chromosomes = '+';'.join(['ll','lll']+CHR_LIST+['l'])+'\n')

def add_tile(INFILE, R1, FILE):
	R0 = R1-0.03
	FILE.write('<plot>\n')
	FILE.write('file    = '+INFILE+'\n')
	FILE.write('show    = yes\n')
	FILE.write('type    = tile\n')
	FILE.write('layers = 15\n')
	FILE.write('margin = 0u\n')
	FILE.write('thickness = 80\n')
	FILE.write('padding = 8\n')
	FILE.write('orientation = out\n')
	FILE.write('stroke_thickness = 1\n')
	FILE.write('stroke_color     = black\n')
	FILE.write('r0    = '+str(R0)+'r\n')
	FILE.write('r1    = '+str(R1)+'r\n')
	FILE.write('<backgrounds>\n')
	FILE.write(' <background>\n')
	FILE.write(' color = un\n')
	FILE.write(' </background>\n')
	FILE.write('</backgrounds>\n')
	FILE.write('</plot>\n')
	FILE.write('\n')
	R1 = R1 - 0.04
	return R1

def add_single_curve_fill_under(INFILE, R1, FILE, FILL):
	R0 = R1-0.07
	FILE.write('<plot>\n')
	FILE.write('z = 10\n')
	FILE.write('type      = line\n')
	FILE.write('thickness = 10\n')
	FILE.write('max_gap = 10u\n')
	FILE.write('file    = '+INFILE+'\n')
	FILE.write('color   = dgrey\n')
	FILE.write('fill_color = dgrey\n')
	FILE.write('min     = 0\n')
	FILE.write('max     = 1\n')
	FILE.write('r0    = '+str(R0)+'r\n')
	FILE.write('r1    = '+str(R1)+'r\n')
	FILE.write('<backgrounds>\n')
	FILE.write(' <background>\n')
	FILE.write(' color = vvgrey\n')
	FILE.write(' </background>\n')
	FILE.write('</backgrounds>\n')
	FILE.write('<axes>\n')
	FILE.write('\n')
	FILE.write('<axis>\n')
	FILE.write('color     = lgrey_a2\n')
	FILE.write('thickness = 5\n')
	FILE.write('spacing   = 1.0r\n')
	FILE.write('</axis>\n')
	FILE.write('</axes>\n')
	FILE.write('</plot>\n')
	FILE.write('\n')
	R1 = R1 - 0.08
	return R1

def __main__():
	#Parse Command Line
	parser = argparse.ArgumentParser(
		description="Takes as input tile and curve files as well as karyotype file "
		"and generate configuration files required to draw circos representation",
		epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr)")
	# Wrapper options.
	parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')
	required.add_argument( '-f',	'--Files',		dest='Files',		required=True,	default=False,		help='Files used in the circos, the type of layer their order from outer to inner they should appear in the circos. Argument should be formated as follows: File1,type,order:File2,type,order:... Where possible types are "tile" or "line" and order is an interger 1, 2, 3, ... Two files can have the same order, they will thus be drawn over each other.')
	required.add_argument( '-o',	'--outfolder',	dest='outfolder',	required=True,	default=False,		help='The output folder that will contain circos files and pictures.')
	required.add_argument( '-p',	'--prefix',		dest='prefix',		required=True,	default=False,		help='Prefix for output files generated: "_housekeeping.conf", ".conf", ".png" and ".svg" files will be generated.')
	required.add_argument( '-k',	'--karyotype',	dest='karyotype',	required=True,	default=False,		help='Path to karyotype file')
	required.add_argument( '-P',	'--color',		dest='color',		required=False,	default=False,		help='A file containing color code for painting.')
	options = parser.parse_args()
	
	# Creating the outfolder if it does not exist
	print('\nCreating the outfolder if it does not exist.')
	os.makedirs(options.outfolder, exist_ok=True)
	
	print('\n Obtaining information')
	
	# Obtaining chromosomes and order from karyotype file
	print('-Obtaining chromosomes and order from karyotype file.')
	chr_list = []
	file = open(options.karyotype, 'r')
	for line in file:
		data = line.split()
		if data:
			chr_list.append(data[2])
	file.close()
	
	# Obtaining layers to draw
	print('-Obtaining layers to draw')
	dicoLayers = {}
	for n in options.Files.split(':'):
		f, t, l = n.split(',')
		l = int(l)
		if not(os.path.isfile(f)):
			sys.exit('ERROR: This is ambarrasing, it seems that the file '+f+' does not exists. The program exited without finishing.')
		if not (l in dicoLayers):
			dicoLayers[l] = []
		dicoLayers[l].append((f, t))
	MinL = min(dicoLayers.keys())
	MaxL = max(dicoLayers.keys())
	for n in range(MinL, MaxL+1):
		if not(n in dicoLayers):
			dicoLayers[n] = []
	
	# Recording color
	print('-Recording color.')
	dico_color = record_color(options.color)
	
	print('\nCreating files required for circos.')
	# Creating the housekeeping.conf file
	print('-Creating the housekeeping.conf file.')
	create_housekeeping(options.outfolder, options.prefix)
	
	# Creating circos file
	print('-Creating circos file.')
	outfile = open(options.outfolder+'/'+options.prefix+'.conf', 'w')
	add_start_file(options.outfolder+'/'+options.prefix+'.png', options.karyotype, options.outfolder+'/'+options.prefix+'_housekeeping.conf', chr_list, dico_color, outfile)
	outfile.write('<plots>\n')
	
	r1 = 0.99
	for n in dicoLayers:
		if len(dicoLayers[n]) == 0:
			r1 = r1 - 0.01
		else:
			if len(dicoLayers[n]) > 1:
				print('--Warning several files have the same order ('+str(n)+'), they will thus be drawn over each other')
				R1 = r1
				for f in dicoLayers[n]:
					if f[1] == 'line':
						r2 = add_single_curve_fill_under(f[0], r1, outfile, 'fill')
						R1 = min(r2, R1)
					elif f[1] == 'tile':
						r2 = add_tile(f[0], r1, outfile)
						R1 = min(r2, R1)
					else:
						sys.exit('--Error type of layer not recognised by this program :'+f[1]+'. The program exited without finishing...\n')
				r1 = R1
			else:
				for f in dicoLayers[n]:
					if f[1] == 'line':
						r1 = add_single_curve_fill_under(f[0], r1, outfile, 'fill')
					elif f[1] == 'tile':
						r1 = add_tile(f[0], r1, outfile)
					else:
						sys.exit('--Error type of layer not recognised by this program :'+f[1]+'. The program exited without finishing...\n')
	
	outfile.write('</plots>\n')
	add_ideogram_section(outfile)
	outfile.close()
	
if __name__ == '__main__': __main__()
