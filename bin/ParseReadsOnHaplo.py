#!/usr/bin/env python
#
#  Copyright 2024 CIRAD
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

import optparse
import sys
import pysam
import gzip
import os


def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(usage="python %prog [options]\n\nProgram designed by Guillaume MARTIN (guillaume.martin@cirad.fr)"
	"\n\nThis program use a bi-allelic haplotype file to sort aligned reads according"
	"\nto haplotype information.")
	# Wrapper options.
	parser.add_option( '-H', '--haplo',		dest='haplo',		default=None,		help='The haplotype file')
	parser.add_option( '-c', '--conf',		dest='conf',		default=None,		help='Configuration file locating bam files to work with')
	parser.add_option( '-p', '--prefix',	dest='prefix',		default=None,		help='Prefix for output files')
	(options, args) = parser.parse_args()
	
	HAPLO = options.haplo
	CONF = options.conf
	PREFIX = options.prefix
	
	# Recording bam to work with
	SetBam = set()
	missing_files = []
	file = open(options.conf)
	for line in file:
		data = line.split()
		if data:
			SetBam.add(data[0])
			if not (os.path.exists(data[0])):
				missing_files.append(data[0])
	file.close()
	
	if missing_files:
		sys.exit('The following file(s) were not found:\n'+'\n'.join(missing_files))
	
	# Recording the haplotypes
	DicoRegions = {}
	file = open(HAPLO)
	header = file.readline().split()
	MarkerPos = header.index('Marker')
	PhasePos = header.index('FinalPhase')
	for line in file:
		data = line.split()
		if data:
			if len(set(data[PhasePos].split('|'))) == 2:
				DicoRegions[data[MarkerPos]] = [data[PhasePos].split('|'), 0, 0]
	file.close()
	
	# Working Bam files per Bam files
	DHaplo1 = {}
	DHaplo2 = {}
	AmbiguousReads = set()
	for bam in SetBam:
		MyBam = pysam.AlignmentFile(bam)
		HEAD = MyBam.header
		ChrAvail = set()
		for CHROM in HEAD['SQ']:
			ChrAvail.add(CHROM['SN'])
		for Marker in DicoRegions:
			PrecPos = {}
			CHR = Marker.split('M')[0]
			START = int(Marker.split('M')[1])-1
			END = START+1
			if CHR in ChrAvail:
				for pileupcolumn in MyBam.pileup(reference=CHR, start=START, end=END):
					# Recording last position before deletion in query
					for pileupread in pileupcolumn.pileups:
						if pileupread.alignment.is_paired:
							if pileupread.alignment.is_read1:
								Rname = pileupread.alignment.query_name + '/1'
							elif pileupread.alignment.is_read2:
								Rname = pileupread.alignment.query_name + '/2'
							else:
								Rname = pileupread.alignment.query_name
								sys.exit('There is a bug... Read '+Rname+' is paired but neither read1 or read2...\n')
						else:
							Rname = pileupread.alignment.query_name
						if pileupread.is_refskip:
							sys.exit('A refskip information has been found in read named: '+pileupread.alignment.query_name+'. I do not know what it is... So the program exited without finishing...\n')
						elif not pileupread.is_del:
							PrecPos[Rname] = pileupread.query_position
					# Working on the studied position
					if pileupcolumn.pos == START:
						for pileupread in pileupcolumn.pileups:
							if pileupread.alignment.is_paired:
								if pileupread.alignment.is_read1:
									Rname = pileupread.alignment.query_name + '/1'
								elif pileupread.alignment.is_read2:
									Rname = pileupread.alignment.query_name + '/2'
								else:
									Rname = pileupread.alignment.query_name
									sys.exit('There is a bug... Read '+Rname+' is paired but neither read1 or read2...\n')
							else:
								Rname = pileupread.alignment.query_name
							Rseq = pileupread.alignment.query_sequence
							RQual = ''.join(list(map(chr, map(lambda x:x+33, list(pileupread.alignment.query_qualities)))))
							RMQual = pileupread.alignment.mapping_quality
							if RMQual != 0:
								OK = True
								if pileupread.is_del:
									base = '*'
									if Rname in PrecPos:
										pos = PrecPos[Rname]
									else:
										OK = False
										sys.stdout.write('Ambigous read: this position is on a deletion and their is no former aligned position of the read. This is why it is removed: '+Rname+'\n')
								elif not pileupread.is_del and not pileupread.is_refskip:
									# query position is None if is_del or is_refskip is set.
									base = pileupread.alignment.query_sequence[pileupread.query_position]
									pos = pileupread.query_position
								else:
									sys.stdout.write('Ambigous read and I do not know why, this is why it is removed: '+Rname+'\n')
									OK = False
								if OK:
									if base in DicoRegions[Marker][0]:
										if DicoRegions[Marker][0].index(base) == 0:
											# print ('\tbase in read %s = %s %s %s %s' % (Rname, base, pos, Rseq, 'haplo1'))
											# outfile1.write('>'+Rname+'|'+Marker+'|'+PREFIX+'-haplo1|'+str(pos)+'\n')
											# outfile1.write(Rseq+'\n+\n')
											# outfile1.write(RQual+'\n')
											
											DicoRegions[Marker][1] += 1
											
											if not (Rname in DHaplo1):
												DHaplo1[Rname] = [Rseq, RQual, [], []]
											DHaplo1[Rname][2].append(Marker)
											DHaplo1[Rname][3].append(str(pos))
												
										elif DicoRegions[Marker][0].index(base) == 1:
											# print ('\tbase in read %s = %s %s %s %s' % (Rname, base, pos, Rseq, 'haplo2'))
											# outfile2.write('>'+Rname+'|'+Marker+'|'+PREFIX+'-haplo2|'+str(pos)+'\n')
											# outfile2.write(Rseq+'\n+\n')
											# outfile2.write(RQual+'\n')
											
											DicoRegions[Marker][2] += 1
											
											if not (Rname in DHaplo2):
												DHaplo2[Rname] = [Rseq, RQual, [], []]
											DHaplo2[Rname][2].append(Marker)
											DHaplo2[Rname][3].append(str(pos))
										else:
											sys.stdout.write('Ambigous read and I do not know why, this is why it is removed: '+Rname+'\n')
									else:
										sys.stdout.write('Ambigous read: '+Rname+' It has a SNP ('+base+') not attributed ('+','.join(DicoRegions[Marker][0])+') at the position '+Marker+'\n')
										AmbiguousReads.add(Rname)
		MyBam.close()
		# print('Ambigous reads:', len(AmbiguousReads))
		# print('Reads found in both haplo', len(DHaplo1&DHaplo2))
		# print('Reads in haplo1:', len(DHaplo1))
		# print('Reads in haplo2:', len(DHaplo2))
	# outfile1.close()
	# outfile2.close()
	
	# Identification of reads in the two haplotypes
	HAPLO1 = set(DHaplo1.keys())
	HAPLO2 = set(DHaplo2.keys())
	NoGoodReads = HAPLO1&HAPLO2
	print('Ambigous reads:', len(AmbiguousReads))
	print('Reads found in both haplo', len(NoGoodReads))
	print('Reads in haplo1:', len(HAPLO1))
	print('Reads in haplo2:', len(HAPLO2))
	
	# printing statistics
	# outfile = open(PREFIX+'.sum', 'w')
	# outfile.write('Marker\tAlleleHaplo1\tAlleleHaplo2\tNumReadHaplo1\tNumReadHaplo2\n')
	# for n in DicoRegions:
		# outfile.write('\t'.join([n]+DicoRegions[n][0]+[str(DicoRegions[n][1]),str(DicoRegions[n][2])]))
		# outfile.write('\n')
	
	# Printing reads grouped in haplotypes
	outfile = gzip.open(PREFIX+'haplo1.fastq.gz', 'wt')
	outfile1 = gzip.open(PREFIX+'haplo1.gz', 'wt')
	for Rname in HAPLO1:
		if not(Rname in NoGoodReads):
			outfile.write('@'+Rname+'\n')
			outfile.write(DHaplo1[Rname][0]+'\n')
			outfile.write('+\n')
			outfile.write(DHaplo1[Rname][1]+'\n')
			outfile1.write('>'+Rname+'|'+':'.join(DHaplo1[Rname][2])+'|'+':'.join(DHaplo1[Rname][3])+'\n')
			outfile1.write(DHaplo1[Rname][0]+'\n')
	outfile.close()
	outfile1.close()
	
	outfile = gzip.open(PREFIX+'haplo2.fastq.gz', 'wt')
	outfile1 = gzip.open(PREFIX+'haplo2.gz', 'wt')
	for Rname in HAPLO2:
		if not(Rname in NoGoodReads):
			outfile.write('@'+Rname+'\n')
			outfile.write(DHaplo2[Rname][0]+'\n')
			outfile.write('+\n')
			outfile.write(DHaplo2[Rname][1]+'\n')
			outfile1.write('>'+Rname+'|'+':'.join(DHaplo2[Rname][2])+'|'+':'.join(DHaplo2[Rname][3])+'\n')
			outfile1.write(DHaplo2[Rname][0]+'\n')
	outfile.close()
	outfile1.close()
	
if __name__ == "__main__": __main__()
