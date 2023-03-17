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

def recode_color(FILE, OUT):
	GROUPS = []
	file  = open(FILE,'r')
	header = file.readline().split()
	outfile = open(OUT, 'w')
	outfile.write('[color]\n')
	for line in file:
		data = line.split()
		if data:
			GROUPS.append(data[header.index('group')])
			outfile.write(data[header.index('group')]+'\tred='+str(int(data[header.index('r')])/255)+':green='+str(int(data[header.index('g')])/255)+':blue='+str(int(data[header.index('b')])/255)+':alpha=1\n')
	outfile.close()
	file.close()
	
	return ':'.join(GROUPS)

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

def __main__():
	#Parse Command Line
	parser = argparse.ArgumentParser(
		description="This program remove the complete genome of a parent from a child, paint the "
		"remaining haplotype/genotype and prepare files for circos visualization.",
		epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr)")
	# Wrapper options.
	parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')
	required.add_argument( '-1',	'--parent1',	dest='parent1',		required=True,	default=False,		help='Potential parent 1 its ploidy and ploidy of the generated gamete. Format should be as follows: ParentName,ParentPloidy,GametePloidy')
	required.add_argument( '-C',	'--child',		dest='child',		required=True,	default=False,		help='Potential child its ploidy. Format should be as follows: ParentName,ParentPloidy. A folder with the acession name will be generated.')
	required.add_argument( '-O',	'--outfolder',	dest='outfolder',	required=True,	default=False,		help='Outfolder in which vcf and painting and circos will be put.')
	required.add_argument( '-v',	'--vcf',		dest='vcf',			required=True,	default=False,		help='Path to the vcf file containing parents and child studied.')
	required.add_argument( '-g',	'--groupFile',	dest='groupFile',	required=True,	default=False,		help='Group file, with a pair accession - group per line.')
	required.add_argument( '-i',	'--inputDir',	dest='inputDir',	required=True,	default=False,		help='Input directory containing grouped alleles and their expected ratio.')
	required.add_argument( '-c',	'--color',		dest='color',		required=True,	default=False,		help='Color file name. Tabulated file with 5 columns with header. Col1: group, col2: name, col3: r, col4: g, col5: b')
	required.add_argument( '-s',	'--size',		dest='size',		required=True,	default=False,		help='A file containing chromosome size. 2 columns are required: col1 : chromosome name, col2 : chromosome size.')
	required.add_argument( '-w',	'--win',		dest='win',			required=True,	default=False,		help='--win option of PaintArp.py: Half window size of markers to use for origin region attribution')
	required.add_argument( '-N',	'--noise',		dest='noise',		required=True,	default=False,		help='--noise option of PaintArp.py: Maximal mean read proportion threshold in which absence of haplotype probability is put to 1.')
	required.add_argument( '-n',	'--threshold',	dest='threshold',	required=True,	default=False,		help='--threshold option of PaintArp.py: Minimal mean read proportion threshold relative to the expected proportion in which haplotype probability is put to 1.')
	required.add_argument( '-r',	'--chro',		dest='chro',		required=False,	default=False,		help='Chromosome to work with. They should be separated by ":".')
	required.add_argument( '-E',	'--centro',		dest='centro',		required=True,	default=False,		help='File locating (peri)centromeric positions. It contained 3 columns: col1 -> chromosome, col2 -> start, col3 -> end')
	options = parser.parse_args()
	
	# Parsing options
	P1n, P1p = options.parent1.split(',')
	Cn, Cp = options.child.split(',')
	
	# Creating folder if it does not exist
	os.makedirs(options.outfolder, exist_ok=True)
	
	# Removing allele from the complete parent restitutor
	sys.stdout.write('\nRemoving allele from the complete parent restitutor\n')
	sys.stdout.flush()
	run_job_error = run_job(getframeinfo(currentframe()), 'vcfRemove.1.0.py -v '+options.vcf+' -A y -r '+P1n+' -a '+Cn+' -o '+options.outfolder+'/'+Cn, 'Error occured during complete restitutor removal step\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Removing supernumerary files that takes many place
	sys.stdout.write('\nRemoving unnecessary files\n')
	sys.stdout.flush()
	os.remove(options.outfolder+'/'+Cn+'_MoreThanHaplo.vcf.gz')
	os.remove(options.outfolder+'/'+Cn+'.vcf.gz')
	
	#Running allele ratio calculation
	sys.stdout.write('\nCreating configuration file for painting\n')
	sys.stdout.flush()
	outfile = open(options.outfolder+'/'+Cn+'.conf', 'w')
	outfile.write(options.outfolder+'/'+Cn+'_haplo.vcf.gz\n')
	outfile.close()
	
	sys.stdout.write('\nRunning allele ratio calculation\n')
	sys.stdout.flush()
	run_job_error = run_job(getframeinfo(currentframe()), 'allele_ratio_per_acc.py -c '+options.outfolder+'/'+Cn+'.conf -g '+options.groupFile+' -a '+Cn+' -i '+options.inputDir+' -o '+options.outfolder, 'Error occured during calculation of ancestral allele ratio in individual step\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Calculating blocks
	sys.stdout.write('\nCalculating ancestral blocks from ratio\n')
	sys.stdout.flush()
	run_job_error = run_job(getframeinfo(currentframe()), 'PaintArp.py  -a '+Cn+' -r '+options.outfolder+'/'+Cn+'_ratio.tab.gz -c '+options.color+' -o '+options.outfolder+'/'+Cn+' -w '+options.win+' -N '+options.noise+' -n '+options.threshold+' -p '+str(int(Cp)-int(P1p))+' -s '+options.size, 'Error occured during calculation of ancestral blocks step\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Drawing curves
	sys.stdout.write('\nDrawing normalized ratios\n')
	sys.stdout.flush()
	run_job_error = run_job(getframeinfo(currentframe()), 'plot_allele_normalized_mean_ratio_per_acc.py -a '+Cn+' -p '+options.outfolder+'/'+Cn+'_curves -r '+options.outfolder+'/'+Cn+'_win_ratio.tab.gz -g png -c '+options.color+' -C '+options.size, 'Error occured during ancestral ratio drawing step\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Drawing blocks
	sys.stdout.write('\nDrawing deduced ancestral blocks\n')
	sys.stdout.flush()
	groups = recode_color(options.color, options.outfolder+'/'+Cn+'.color.conf')
	run_job_error = run_job(getframeinfo(currentframe()), 'haplo2kar.1.0.py --acc '+options.outfolder+'/'+Cn+' --chr '+options.chro+' --gcol '+options.outfolder+'/'+Cn+'.color.conf --dg '+groups+' --ploidy '+str(int(Cp)-int(P1p))+' --centro '+options.centro+' --graph svg', 'Error occured during ancestral block drawing step\n')
	if run_job_error:
		print(run_job_error)
		sys.exit()
	
	# Formatting haplotype file for circos
	sys.stdout.write('\nFormatting haplotype file for circos\n')
	sys.stdout.write('-Recording color.\n')
	sys.stdout.flush()
	dico_color = record_color(options.color)
	sys.stdout.write('-Obtaining chromosome list.\n')
	sys.stdout.flush()
	chr_list = options.chro.split(':')
	sys.stdout.write('-Creating haplotype file(s).\n')
	sys.stdout.flush()
	format_haplo(options.outfolder, options.outfolder, dico_color, 'Removed', Cn, int(Cp)-int(P1p), chr_list)
	
if __name__ == '__main__': __main__()
