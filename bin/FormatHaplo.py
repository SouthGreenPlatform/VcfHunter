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

def format_haplo(OUTFOLDER, FOLDER, COLOR, IND, PLO, CHR):
	for i in range(PLO):
		Indice = str(i+1)
		outfile = open(OUTFOLDER+'/circos.'+IND+'.haplo'+str(Indice)+'.tab','w')
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

def __main__():
	#Parse Command Line
	parser = argparse.ArgumentParser(
		description="This program takes the mosaic painting from PaintArp and convert them to tile files for circos representation.",
		epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr)")
	# Wrapper options.
	parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')
	required.add_argument( '-i',	'--indiv',		dest='indiv',		required=True,	default=False,		help='Individual for which the haplotypes should be reformated for Circos drawing ans its ploidy. Format should be as follows: Name,Ploidy.')
	required.add_argument( '-v',	'--vcf',		dest='vcf',			required=True,	default=False,		help='Path to the vcf file containing individual studied.')
	required.add_argument( '-p',	'--painting',	dest='painting',	required=False,	default=False,		help='A folder containing chromosome painting for all accessions studied from vcfhunter tools.')
	required.add_argument( '-P',	'--color',		dest='color',		required=False,	default=False,		help='A file containing color code for painting.')
	optional.add_argument( '-c',	'--chr',		dest='chr',			required=False,	default=False,		help='Chromosome to work with. They should be separated by ":". If omitted, all chromosomes will be used.')
	optional.add_argument( '-o',	'--outfolder',	dest='outfolder',	required=False,	default=False,		help='Folder in which output should be written.')
	options = parser.parse_args()
	
	# Parsing options
	P1n, P1p = options.indiv.split(',')
	
	# Creating folder if it does not exist
	os.makedirs(options.outfolder, exist_ok=True)
	
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
	
	dico_color = record_color(options.color)
			
	format_haplo(options.outfolder, options.painting, dico_color, P1n, int(P1p), chr_list)
			
if __name__ == '__main__': __main__()
