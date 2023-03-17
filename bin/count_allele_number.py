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

import configparser
import optparse
import os
import sys
import utilsSR.utilsSR as utils

def __main__():
	#Parse Command Line
	parser = optparse.OptionParser(description="This program count the number or read supporting alleles of each covered bases.",
									epilog="Program designed by Guillaume MARTIN (guillaume.martin@cirad.fr)")
	# Wrapper options.
	parser.add_option( '-r', '--ref', dest='ref', default=None, help='Reference multifasta file, [default: %default]')
	parser.add_option( '-b', '--bam', dest='bam', default=None, help='The bam file, [default: %default]')
	parser.add_option( '-o', '--out', dest='out', default=None, help='The output file name, [default: %default]')
	(options, args) = parser.parse_args()
	
	# locating loca_programs.conf file
	pathname = os.path.dirname(sys.argv[0])
	loca_programs = configparser.RawConfigParser()
	loca_programs.read(pathname+'/loca_programs.conf')
	bamreadcount = loca_programs.get('Programs','bamreadcount')
	
	# running bamreadcount
	utils.perform_count(bamreadcount, options.ref, options.bam, options.out+'.temp')
	
	# parsing bamreadcount file
	utils.generate_pseudo_vcf(options.out+'.temp', options.ref, options.out)
	

	
	
if __name__ == "__main__": __main__()
