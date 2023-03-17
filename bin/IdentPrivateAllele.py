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
import multiprocessing as mp
import subprocess
import optparse
import gzip
import sys
import os

def concat(OUT, FILES):
	print('Concatenated files in ', OUT, ':', '\t'.join(FILES))
	outfile = gzip.open(OUT, 'wt')
	for file in FILES:
		f = gzip.open(file, 'rt')
		for line in f:
			outfile.write(line)
		f.close()
	outfile.close()
	
	for file in FILES:
		os.remove(file) 

def Parallelize(job):
	print('cmd: ', job)
	pipe = subprocess.Popen(job, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	res = pipe.communicate()
	return (pipe.returncode, res[1])

def Parallel(job):
	
	try:
		rslt = concat(job[0],job[1])
	except Exception as e:
		print (e)
		rslt = 1
	finally:
		return rslt


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN (guillaume.martin@cirad.fr)"
	"\n\nThis program run 'vcf2allPropAndCov.py' in chain in order to identify"
	"\nspecific alleles of genetic groups.")
	# Wrapper options.
	
	parser.add_option( '-c', '--conf',					dest='conf',					default=None,		help='Path to a file containing path to one or multiple vcf files (one per line)')
	parser.add_option( '-g', '--group',					dest='group',					default=None,		help='A two column file with accession in the first column and group tag (i.e. origin) in the second column')
	parser.add_option( '-o', '--outdir',				dest='outdir',					default='step1',	help='Path to the output directory, where the program will put the subdirectories per accession')
	parser.add_option( '-t', '--thread',				dest='thread',					default='1',		help='Number of processors available. [Default: %default]')
	parser.add_option( '-a', '--param_v2apac_all',		dest='param_v2apac_all',		default='n',		help='vcf2allPropAndCov parameter --all: allele should be present in all accessions of the group. Possible values "y", or "n". [Default: %default]')
	parser.add_option( '-i', '--param_v2apac_introg',	dest='param_v2apac_introg',		default=None,		help='vcf2allPropAndCov parameter --excl: a tabulated file locating introgression in ancestral accessions')
	parser.add_option( '-p', '--param_v2apac_prop',		dest='param_v2apac_prop',		default='n',		help='vcf2allPropAndCov parameter --prop: allele proportion in ancestral accessions. Value comprised between 0 and 1 or "n" if not using this parameter. [Default: %default]')
	parser.add_option( '-m', '--param_v2apac_NoMiss',	dest='param_v2apac_NoMiss',		default='y',		help='vcf2allPropAndCov parameter --NoMiss: No missing data are allowed in accessions used to group alleles. Value "y" for not allowing missing data, "n" for allowing missing data. [Default: %default]')
	(options, args) = parser.parse_args()
	
	# Verification of required argument are passed
	if options.conf == None:
		sys.exit('Please provide a configuration file to --conf argument')
	if options.group == None:
		sys.exit('Please provide a configuration file to --group argument')
	
	# Recording number of available processors
	nbProcs = int(options.thread)
	
	# Looking for vcf2allPropAndCov.py
	pathname = os.path.dirname(os.path.abspath(sys.argv[0]))
	
	# Creating the directory if it does not exists
	os.makedirs(options.outdir, exist_ok=True)
	
	# Recording group information
	dicoAcc = set()
	file = open(options.group)
	for line in file:
		data = line.split()
		if data:
			dicoAcc.add(data[0])
	file.close()
	
	# Recording vcf information
	ListVcf = []
	file = open(options.conf)
	for line in file:
		data = line.split()
		if data:
			ListVcf.append(data[0])
	file.close()
	
	# Creating subdirectories if they does not exist and preparing data for multiprocessing
	i = 0
	listJobs = []
	for vcf in ListVcf:
		i+=1
		for acc in dicoAcc:
			os.makedirs(options.outdir+'/'+acc, exist_ok=True)
			if options.param_v2apac_introg == None:
				cmd = pathname+'/vcf2allPropAndCov.py --vcf '+vcf+' --origin '+options.group+' --acc '+acc+' --ploidy 2 --NoMiss '+options.param_v2apac_NoMiss+' --prop '+options.param_v2apac_prop+' --all '+options.param_v2apac_all+' --prefix '+options.outdir+'/'+acc+'/tmp_'+str(i)+'_'
			else:
				cmd = pathname+'/vcf2allPropAndCov.py --vcf '+vcf+' --origin '+options.group+' --acc '+acc+' --ploidy 2 --NoMiss '+options.param_v2apac_NoMiss+' --prop '+options.param_v2apac_prop+' --all '+options.param_v2apac_all+' --prefix '+options.outdir+'/'+acc+'/tmp_'+str(i)+'_ --excl '+options.param_v2apac_introg
			print(cmd)
			listJobs.append(cmd)
	
	pool = mp.Pool(processes=nbProcs)
	results = pool.map(Parallelize, listJobs)	
	
	for n in results:
		if n[0] != 0:
			for line in n[1].decode(encoding='utf-8').split('\n'):
				print(line)
	
	# Merging files
	listJobs = []
	for acc in dicoAcc:
		i = 0
		listFiles = []
		for vcf in ListVcf:
			i+=1
			listFiles.append(options.outdir+'/'+acc+'/tmp_'+str(i)+'_'+acc+'_AlleleOriginAndRatio.tab.gz')
		listJobs.append((options.outdir+'/'+acc+'/'+acc+'_ratio.tab.gz', listFiles))
	
	pool = mp.Pool(processes=nbProcs)
	results = pool.map(Parallel, listJobs)

if __name__ == "__main__": __main__()
