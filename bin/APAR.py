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

def create_karyotype(CN, CROSS, CHR_LIST, DICO_CHR_LENGTH):
	outfile = open(CN+'/'+CROSS+'.karyotype.tab', 'w')
	outfile.write('\t'.join(['chr', '-', 'll', 'll', '0', '1', 'white']))
	outfile.write('\n')
	outfile.write('\t'.join(['chr', '-', 'lll', 'lll', '0', '1', 'white']))
	outfile.write('\n')
	for chro in CHR_LIST:
		outfile.write('\t'.join(['chr', '-', chro, chro, '0', str(DICO_CHR_LENGTH[chro]), 'white']))
		outfile.write('\n')
	outfile.write('\t'.join(['chr', '-', 'l', 'l', '0', '1', 'white']))
	outfile.write('\n')
	outfile.close()

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

def add_tile(CN, IND, CP, CROSS, R1, FILE):
	R0 = R1-0.03
	for i in range(CP):
		indice = str(i+1)
		FILE.write('<plot>\n')
		FILE.write('file    = '+CN+'/'+CROSS+'.circos.'+IND+'.haplo'+indice+'.tab\n')
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
		R0 = R1 - 0.03
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

def circos_both_parents_comb(CN, CP, P1N, P1P, P2N, P2P, CROSS, DICO_COLOR, CHR_LIST):
	outfile = open(CN+'/'+CROSS+'.circos.ParentComb.conf', 'w')
	add_start_file(CN+'/'+CROSS+'-circos-ParentComb.png', CN+'/'+CROSS+'.karyotype.tab', CN+'/'+CROSS+'_housekeeping.conf', CHR_LIST, DICO_COLOR, outfile)
	outfile.write('<plots>\n')
	
	r1 = 0.99
	r1 = add_tile(CN, CN, int(CP), CROSS, r1, outfile)
	r1 = add_single_curve_fill_under(CN+'/'+CROSS+'_2Parents_NoOK_prop_SNPwin.tab', r1, outfile, 'fill')
	r1 = add_tile(CN, P1N, int(P1P), CROSS, r1, outfile)
	r1 = r1 - 0.02
	r1 = add_tile(CN, P2N, int(P2P), CROSS, r1, outfile)
	
	outfile.write('</plots>\n')
	add_ideogram_section(outfile)
	outfile.close()

def circos_one_parent_complete(CN, CP, PN, PP, CROSS, DICO_COLOR, CHR_LIST):
	outfile = open(CN+'/'+CROSS+'.circos.'+PN+'.conf', 'w')
	add_start_file(CN+'/'+CROSS+'-circos-'+PN+'.png', CN+'/'+CROSS+'.karyotype.tab', CN+'/'+CROSS+'_housekeeping.conf', CHR_LIST, DICO_COLOR, outfile)
	outfile.write('<plots>\n')
	
	r1 = 0.99
	r1 = add_tile(CN, CN, int(CP), CROSS, r1, outfile)
	r1 = add_single_curve_fill_under(CN+'/'+CROSS+'_1Parent-'+PN+'_NoOK_prop_SNPwin.tab', r1, outfile, 'fill')
	r1 = add_single_curve_fill_under(CN+'/'+CROSS+'_allele-'+PN+'_NoOK_prop.tab', r1, outfile, 'fill')
	r1 = add_tile(CN, PN, int(PP), CROSS, r1, outfile)
	
	outfile.write('</plots>\n')
	add_ideogram_section(outfile)
	outfile.close()

def circos_one_parent(CN, CP, PN, PP, CROSS, DICO_COLOR, CHR_LIST):
	outfile = open(CN+'/'+CROSS+'.circos.'+PN+'.conf', 'w')
	add_start_file(CN+'/'+CROSS+'-circos-'+PN+'.png', CN+'/'+CROSS+'.karyotype.tab', CN+'/'+CROSS+'_housekeeping.conf', CHR_LIST, DICO_COLOR, outfile)
	outfile.write('<plots>\n')
	
	r1 = 0.99
	r1 = add_tile(CN, CN, int(CP), CROSS, r1, outfile)
	r1 = add_single_curve_fill_under(CN+'/'+CROSS+'_1Parent-'+PN+'_NoOK_prop_SNPwin.tab', r1, outfile, 'fill')
	r1 = add_tile(CN, PN, int(PP), CROSS, r1, outfile)
	
	outfile.write('</plots>\n')
	add_ideogram_section(outfile)
	outfile.close()
	
def circos_one_parent_against_haplo(CN, CP, PN, PP, CROSS, DICO_COLOR, CHR_LIST):
	outfile = open(CN+'/'+CROSS+'.circos.haplo.'+PN+'.conf', 'w')
	add_start_file(CN+'/'+CROSS+'-circos-'+CN+'-haplo-'+PN+'.png', CN+'/'+CROSS+'.karyotype.tab', CN+'/'+CROSS+'_housekeeping.conf', CHR_LIST, DICO_COLOR, outfile)
	outfile.write('<plots>\n')
	
	r1 = 0.99
	r1 = add_tile(CN, CN, int(CP), CROSS, r1, outfile)
	r1 = add_single_curve_fill_under(CN+'/'+CROSS+'.'+PN+'.haplo_NoOK_prop_SNPwin.tab', r1, outfile, 'fill')
	r1 = add_tile(CN, PN, int(PP), CROSS, r1, outfile)
	
	outfile.write('</plots>\n')
	add_ideogram_section(outfile)
	outfile.close()

def __main__():
	#Parse Command Line
	parser = argparse.ArgumentParser(
		description="This programme goes through all the steps of parent-child trio analysis. "
		"It includes a chromosomal analysis of SNPs validating the trio and the duos (child, "
		"parent), a search for the complete genome of the parent in the child if appropriate "
		"(ploidy of the parent gamete is equal to the ploidy of the parent), an attempt to "
		"remove the complete genome of the parent from the child if appropriate (ploidy of "
		"the parent gamete is equal to the ploidy of the parent) and a file preparation for "
		"circos visualisation.",
		epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr)")
	# Wrapper options.
	parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')
	required.add_argument( '-1',	'--parent1',	dest='parent1',		required=True,	default=False,		help='Potential parent 1 its ploidy and ploidy of the generated gamete. Format should be as follows: ParentName,ParentPloidy,GametePloidy')
	required.add_argument( '-2',	'--parent2',	dest='parent2',		required=True,	default=False,		help='Potential parent 2 its ploidy and ploidy of the generated gamete. Format should be as follows: ParentName,ParentPloidy,GametePloidy')
	required.add_argument( '-C',	'--child',		dest='child',		required=True,	default=False,		help='Potential child its ploidy. Format should be as follows: ParentName,ParentPloidy. A folder with the acession name will be generated.')
	required.add_argument( '-n',	'--crossname',	dest='crossname',	required=True,	default=False,		help='Name of the cross. All generated files will be prefixed with this name.')
	required.add_argument( '-v',	'--vcf',		dest='vcf',			required=True,	default=False,		help='Path to the vcf file containing parents and child studied.')
	required.add_argument( '-f',	'--fasta',		dest='fasta',		required=True,	default=False,		help='The multifasta reference file.')
	required.add_argument( '-p',	'--painting',	dest='painting',	required=False,	default=False,		help='A folder containing chromosome painting for all accessions studied from vcfhunter tools.')
	required.add_argument( '-P',	'--color',		dest='color',		required=False,	default=False,		help='A file containing color code for painting.')
	optional.add_argument( '-c',	'--chr',		dest='chr',			required=False,	default=False,		help='Chromosome to work with. They should be separated by ":". If omitted, all chromosomes will be used.')
	optional.add_argument( '-w',	'--window',		dest='window',		required=False,	default='100',		help='Half window to calculate site proportion in accordance with tested parentage (unit = variant site). [Default: 100]')
	optional.add_argument( '-W',	'--WINDOW',		dest='WINDOW',		required=False,	default='100000',	help='The window to calculate site proportion in accordance with tested parentage (unit base pair). [Default: 100000]')
	options = parser.parse_args()
	
	# Parsing options
	P1n, P1p, P1g = options.parent1.split(',')
	P2n, P2p, P2g = options.parent2.split(',')
	Cn, Cp = options.child.split(',')
	
	# Creating folder if it does not exist
	os.makedirs(Cn, exist_ok=True)
	
	# Obtaining chromosome names and additional information
	sys.stdout.write('\nObtaining chromosome name and length\n')
	sys.stdout.flush()
	dico_chr_length = {}
	if options.vcf[-3:] == '.gz':
		file = gzip.open(options.vcf,'rt')
	else:
		file = open(options.vcf)
	for line in file:
		data = line.split()
		if data:
			if data[0][0] == "#":
				if data[0][0:13] == "##contig=<ID=":
					chrlen = int(data[0].split('=')[3].split('>')[0])
					chrname = data[0].split('=')[2].split(',')[0]
					dico_chr_length[chrname] = chrlen
			else:
				break
	file.close()
	
	if options.chr:
		chr_list = options.chr.split(':')
	else:
		chr_list = sorted(dico_chr_length.keys())
	
	
	# Creating accession name files required for filtering
	outfile=open(Cn+'/'+options.crossname+'.name.tab','w')
	outfile.write('\n'.join((P1n, P2n, Cn)))
	outfile.write('\t')
	outfile.close()
	outfile=open(Cn+'/'+options.crossname+'.'+P1n+'.name.tab','w')
	outfile.write('\n'.join((P1n, Cn)))
	outfile.write('\t')
	outfile.close()
	outfile=open(Cn+'/'+options.crossname+'.'+P2n+'.name.tab','w')
	outfile.write('\n'.join((P2n, Cn)))
	outfile.write('\t')
	outfile.close()
	
	# Running vcf filtering for selected parents-child trios
	sys.stdout.write('\nFiltering initial vcf file\n')
	sys.stdout.flush()
	run_job_error = run_job(getframeinfo(currentframe()), 'vcfFilter.1.0.py --vcf '+options.vcf+' --prefix '+Cn+'/'+options.crossname+' -g y --MinCov 10 --MinAl 3 --nMiss 0 --MaxCov 10000 --MinFreq 0.1 --RmAlAlt 1:3:4:5:6:7:8:9:10 --names  '+Cn+'/'+options.crossname+'.name.tab', 'Error occured in vcf filtering step\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Running only parent1 and child vcf filtering
	sys.stdout.write('\nFiltering vcf file with '+Cn+' and '+P1n+'\n')
	sys.stdout.flush()
	run_job_error = run_job(getframeinfo(currentframe()), 'vcfFilter.1.0.py --vcf '+Cn+'/'+options.crossname+'_filt.vcf.gz --prefix '+Cn+'/'+options.crossname+'.'+P1n+' -g y --MinCov 10 --MinAl 3 --nMiss 0 --MaxCov 10000 --MinFreq 0.1 --RmAlAlt 1:3:4:5:6 --names  '+Cn+'/'+options.crossname+'.'+P1n+'.name.tab', 'Error occured in vcf filtering step of parent1\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Running only parent2 and child vcf filtering
	sys.stdout.write('\nFiltering vcf file with '+Cn+' and '+P2n+'\n')
	sys.stdout.flush()
	run_job_error = run_job(getframeinfo(currentframe()), 'vcfFilter.1.0.py --vcf '+Cn+'/'+options.crossname+'_filt.vcf.gz --prefix '+Cn+'/'+options.crossname+'.'+P2n+' -g y --MinCov 10 --MinAl 3 --nMiss 0 --MaxCov 10000 --MinFreq 0.1 --RmAlAlt 1:3:4:5:6 --names  '+Cn+'/'+options.crossname+'.'+P2n+'.name.tab', 'Error occured in vcf filtering step of parent2\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Testing parentage with both parents
	sys.stdout.write('\nTesting parentage with both parents\n')
	sys.stdout.flush()
	if options.chr:
		run_job_error = run_job(getframeinfo(currentframe()), 'ACRO.py -w 100 -1 '+P1n+','+P1g+' -2 '+P2n+','+P2g+' -a '+Cn+' -v '+Cn+'/'+options.crossname+'_filt.vcf.gz -p '+Cn+'/'+options.crossname+'_2Parents --fasta '+options.fasta+' --chr '+options.chr, 'Error occured during testing of both parents\n')
	else:
		run_job_error = run_job(getframeinfo(currentframe()), 'ACRO.py -w 100 -1 '+P1n+','+P1g+' -2 '+P2n+','+P2g+' -a '+Cn+' -v '+Cn+'/'+options.crossname+'_filt.vcf.gz -p '+Cn+'/'+options.crossname+'_2Parents --fasta '+options.fasta, 'Error occured during testing of both parents\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Testing parentage with parent 1
	sys.stdout.write('\nTesting parentage with '+P1n+'\n')
	sys.stdout.flush()
	if options.chr:
		run_job_error = run_job(getframeinfo(currentframe()), 'ACRO.py -w 100 -1 '+P1n+','+P1g+' -a '+Cn+' -v '+Cn+'/'+options.crossname+'.'+P1n+'_filt.vcf.gz -p '+Cn+'/'+options.crossname+'_1Parent-'+P1n+' --fasta '+options.fasta+' --chr '+options.chr, 'Error occured during testing gamete of '+P1n+'\n')
	else:
		run_job_error = run_job(getframeinfo(currentframe()), 'ACRO.py -w 100 -1 '+P1n+','+P1g+' -a '+Cn+' -v '+Cn+'/'+options.crossname+'.'+P1n+'_filt.vcf.gz -p '+Cn+'/'+options.crossname+'_1Parent-'+P1n+' --fasta '+options.fasta, 'Error occured during testing gamete of '+P1n+'\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Testing parentage with parent 2
	sys.stdout.write('\nTesting parentage with '+P2n+'\n')
	sys.stdout.flush()
	if options.chr:
		run_job_error = run_job(getframeinfo(currentframe()), 'ACRO.py -w 100 -1 '+P2n+','+P2g+' -a '+Cn+' -v '+Cn+'/'+options.crossname+'.'+P2n+'_filt.vcf.gz -p '+Cn+'/'+options.crossname+'_1Parent-'+P2n+' --fasta '+options.fasta+' --chr '+options.chr, 'Error occured during testing gamete of '+P2n+'\n')
	else:
		run_job_error = run_job(getframeinfo(currentframe()), 'ACRO.py -w 100 -1 '+P2n+','+P2g+' -a '+Cn+' -v '+Cn+'/'+options.crossname+'.'+P2n+'_filt.vcf.gz -p '+Cn+'/'+options.crossname+'_1Parent-'+P2n+' --fasta '+options.fasta, 'Error occured during testing gamete of '+P2n+'\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Testing number of shared alleles between parent 1 and child
	sys.stdout.write('\nTesting number of shared alleles between '+P1n+' and '+Cn+'. This test alowed to look for complete genome restitution.\n')
	sys.stdout.flush()
	run_job_error = run_job(getframeinfo(currentframe()), 'vcfIdent.1.0.py -v '+Cn+'/'+options.crossname+'.'+P1n+'_filt.vcf.gz -a '+Cn+' -c '+P1n+' -o '+Cn+'/'+options.crossname+'_allele -w 100 -d n', 'Error occured during testing of '+P1n+' alleles\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Testing number of shared alleles between parent 2 and child
	sys.stdout.write('\nTesting number of shared alleles between '+P2n+' and '+Cn+'. This test alowed to look for complete genome restitution.\n')
	sys.stdout.flush()
	run_job_error = run_job(getframeinfo(currentframe()), 'vcfIdent.1.0.py -v '+Cn+'/'+options.crossname+'.'+P2n+'_filt.vcf.gz -a '+Cn+' -c '+P2n+' -o '+Cn+'/'+options.crossname+'_allele -w 100 -d n', 'Error occured during testing of '+P2n+' alleles\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Removing alleles from parents that have the same ploidy than the gamete they gives.
	if P1p == P1g:
		sys.stdout.write('\nRemoving alleles from parent '+P1n+' that have the same ploidy that the gamete it gives. In case of complete gamete restitution this gives acces to the gamete composition of the second parent.\n')
		sys.stdout.flush()
		run_job_error = run_job(getframeinfo(currentframe()), 'vcfRemove.1.0.py -v '+Cn+'/'+options.crossname+'_filt.vcf.gz -a '+Cn+' -r '+P1n+' -o '+Cn+'/'+options.crossname+'.ss.'+P1n, 'Error occured during testing of '+P1n+' alleles\n')
		if run_job_error:
			print(run_job_error)
			sys.exit()
	if P2p == P2g:
		sys.stdout.write('\nRemoving alleles from parent '+P2n+' that have the same ploidy that the gamete it gives. In case of complete gamete restitution this gives acces to the gamete composition of the second parent.\n')
		sys.stdout.flush()
		run_job_error = run_job(getframeinfo(currentframe()), 'vcfRemove.1.0.py -v '+Cn+'/'+options.crossname+'_filt.vcf.gz -a '+Cn+' -r '+P2n+' -o '+Cn+'/'+options.crossname+'.ss.'+P2n, 'Error occured during testing of '+P2n+' alleles\n')
		if run_job_error:
			print(run_job_error)
			sys.exit()
	
	# Analysing remaining haplotype(s) after removal of one parent
	if P1p == P1g:
		sys.stdout.write('\nAnalysing remaining haplotype(s) of '+Cn+' with '+P2n+' after removal of alleles from '+P1n+'. This makes sens only if '+P1n+' as performed a full genome restitution.\n')
		sys.stdout.flush()
		if options.chr:
			run_job_error = run_job(getframeinfo(currentframe()), 'ACRO.py -1 '+P2n+','+P2g+' -a '+Cn+' -v '+Cn+'/'+options.crossname+'.ss.'+P1n+'_haplo.vcf.gz -p '+Cn+'/'+options.crossname+'.'+P2n+'.haplo --fasta '+options.fasta+' --chr '+options.chr, 'Error occured during testing gamete of '+P2n+' with remaining haplotypes of '+Cn+'\n')
		else:
			run_job_error = run_job(getframeinfo(currentframe()), 'ACRO.py -1 '+P2n+','+P2g+' -a '+Cn+' -v '+Cn+'/'+options.crossname+'.ss.'+P1n+'_haplo.vcf.gz -p '+Cn+'/'+options.crossname+'.'+P2n+'.haplo --fasta '+options.fasta, 'Error occured during testing gamete of '+P2n+' with remaining haplotypes of '+Cn+'\n')
		if run_job_error:
			print(run_job_error)
			sys.exit()
	if P2p == P2g:
		sys.stdout.write('\nAnalysing remaining haplotype(s) of '+Cn+' with '+P1n+' after removal of alleles from '+P2n+'. This makes sens only if '+P2n+' as performed a full genome restitution.\n')
		sys.stdout.flush()
		if options.chr:
			run_job_error = run_job(getframeinfo(currentframe()), 'ACRO.py -1 '+P1n+','+P1g+' -a '+Cn+' -v '+Cn+'/'+options.crossname+'.ss.'+P2n+'_haplo.vcf.gz -p '+Cn+'/'+options.crossname+'.'+P1n+'.haplo --fasta '+options.fasta+' --chr '+options.chr, 'Error occured during testing gamete of '+P1n+' with remaining haplotypes of '+Cn+'\n')
		else:
			run_job_error = run_job(getframeinfo(currentframe()), 'ACRO.py -1 '+P1n+','+P1g+' -a '+Cn+' -v '+Cn+'/'+options.crossname+'.ss.'+P2n+'_haplo.vcf.gz -p '+Cn+'/'+options.crossname+'.'+P1n+'.haplo --fasta '+options.fasta, 'Error occured during testing gamete of '+P1n+' with remaining haplotypes of '+Cn+'\n')
		if run_job_error:
			print(run_job_error)
			sys.exit()
	
	# Preparing circos files if a folder containing accession chromosome painting is provided as well as a file with color code.
	if options.painting:
		if options.color:
			sys.stdout.write('\nPreparing files for circos drawing.\n')
			sys.stdout.flush()
			
			sys.stdout.write('-Recording color.\n')
			sys.stdout.flush()
			dico_color = record_color(options.color)
			
			sys.stdout.write('-Formatting chromosome painting files for circos representation.\n')
			sys.stdout.flush()
			sys.stdout.write('--'+P1n+'\n')
			sys.stdout.flush()
			format_haplo(Cn, options.painting, dico_color, options.crossname, P1n, int(P1p), chr_list)
			
			sys.stdout.write('--'+P2n+'\n')
			sys.stdout.flush()
			format_haplo(Cn, options.painting, dico_color, options.crossname, P2n, int(P2p), chr_list)
			
			sys.stdout.write('--'+Cn+'\n')
			sys.stdout.flush()
			format_haplo(Cn, options.painting, dico_color, options.crossname, Cn, int(Cp), chr_list)
			
			sys.stdout.write('-Creating the circos files.\n')
			sys.stdout.flush()
			
			sys.stdout.write('--Creating the karyotype file.\n')
			sys.stdout.flush()
			create_karyotype(Cn, options.crossname, chr_list, dico_chr_length)
			
			sys.stdout.write('--Creating the housekeeping.conf file.\n')
			sys.stdout.flush()
			create_housekeeping(Cn, options.crossname)
			
			sys.stdout.write('--Creating the circos file of data testing both parent combination.\n')
			sys.stdout.flush()
			circos_both_parents_comb(Cn, Cp, P1n, P1p, P2n, P2p, options.crossname, dico_color, chr_list)
			
			sys.stdout.write('--Creating the circos file of data testing '+P1n+'.\n')
			sys.stdout.flush()
			if P1p == P1g:
				circos_one_parent_complete(Cn, Cp, P1n, P1p, options.crossname, dico_color, chr_list)
			else:
				circos_one_parent(Cn, Cp, P1n, P1p, options.crossname, dico_color, chr_list)
			
			sys.stdout.write('--Creating the circos file of data testing '+P2n+'.\n')
			sys.stdout.flush()
			if P2p == P2g:
				circos_one_parent_complete(Cn, Cp, P2n, P2p, options.crossname, dico_color, chr_list)
			else:
				circos_one_parent(Cn, Cp, P2n, P2p, options.crossname, dico_color, chr_list)
			
			if P1p == P1g:
				sys.stdout.write('--Creating the circos file of data testing '+P2n+' with remaining haplotype(s) of the child '+Cn+'. This makes sense only if '+P1n+' makes a complete genome restitution.\n')
				sys.stdout.flush()
				circos_one_parent_against_haplo(Cn, Cp, P2n, P2p, options.crossname, dico_color, chr_list)
			else:
				pass
			
			if P2p == P2g:
				sys.stdout.write('--Creating the circos file of data testing '+P1n+' with remaining haplotype(s) of the child '+Cn+'. This makes sense only if '+P2n+' makes a complete genome restitution.\n')
				sys.stdout.flush()
				circos_one_parent_against_haplo(Cn, Cp, P1n, P1p, options.crossname, dico_color, chr_list)
			else:
				pass
			
		else:
			sys.stdout.write('\nA painting file is passed but no color file provided. The program will thus stop here and not go until circos drawing.\n')
			sys.stdout.flush()
	
if __name__ == '__main__': __main__()
