#!/usr/bin/env bash
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

#-----------
# Parameters

TEMP=`getopt -o r:d:f:w:i:e:l:p:s:h: --long reference:,directory:,read-folder:,window:,include-pattern:,exclude-pattern:,color-file:,prefix:,steps:,help -- "$@"`
eval set -- "$TEMP"

help_line="Help: Paint assembly based on origin specific reads\n
-h|--help \n
\tprint this help\n

-r|--reference \n
\t Path to the reference fasta file\n

-d|--directory \n
\tA directory that will contains outputs\n

-f|--read-folder \n
\tA folder containing origin reads. A file per origin named ORIGIN.fastq.gz\n

-w|--window \n
\tWindows in which the number of mean number of reads hits will be counted\n

-i|--include-pattern \n
\tPattern(s) to select chromosomes into the figure. Each patterns should be separated by '|'\n

-e|--exclude-pattern \n
\tPattern(s) to exclude some chromosomes. Each patterns should be separated by '|'\n

-l|--color-file \n
\tA color file with 4 columns: col 1 = group name and col 2 to 4 = RGB color code\n

-p|--prefix \n
\tPrefix for the output figure and intermediate output files\n

-s|--steps \n
\tSteps of the analysis to perform:\n
\t1: Working folder creation, copy of the reference and indexation\n
\t2: Origin read libraries mapping\n
\t3: Mapping read filtration\n
\t4: Calculating the number of hits on non overlapping sliding windows\n
\t5: Selecting chromosomes on pattern and drawing curve figures\n
\t6: Attributing origin based on majority rule\n
\t7: Data formating for GEMO analysis\n"


# Default variable 
Ref="NoNaMe"
Dir="NoNaMe"
ReadFolder="NoNaMe"
WIN="NoNaMe"
CHRPAT="NoNaMe"
EXCCHRPAT="NoNaMe"
COLFILE="NoNaMe"
OUTFIGURE="NoNaMe"
STEPS="12345"

# Extract options and their arguments into variables.

while true ; do
    case "$1" in
        -r|--reference)
            Ref=$2 ; shift 2;;
        -d|--directory)
            Dir=$2 ; shift 2;;
        -f|--read-folder)
            ReadFolder=$2 ; shift 2;;
        -w|--window)
            WIN=$2 ; shift 2;;
        -i|--include-pattern)
            CHRPAT=$2 ; shift 2;;
        -e|--exclude-pattern)
            EXCCHRPAT=$2 ; shift 2;;
        -l|--color-file)
            COLFILE=$2 ; shift 2;;
        -p|--prefix)
            OUTFIGURE=$2 ; shift 2;;
        -s|--steps)
            STEPS=$2 ; shift 2;;
        -h|--help)
            echo -e $help_line ; exit 0;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1;;
    esac
done

if [ "$Ref" = "NoNaMe" ]; then
 echo -e "\nMissing required argument --reference\n"; echo -e $help_line; exit 0;
fi

if [ "$Dir" = "NoNaMe" ]; then
 echo -e "\nMissing required argument --directory\n"; echo -e $help_line; exit 0;
fi

if [ "$ReadFolder" = "NoNaMe" ]; then
 echo -e "\nMissing required argument --read-folder\n"; echo -e $help_line; exit 0;
fi

if [ "$WIN" = "NoNaMe" ]; then
 echo -e "\nMissing required argument --window\n"; echo -e $help_line;  exit 0;
fi

if [ "$CHRPAT" = "NoNaMe" ]; then
 echo -e "\nMissing required argument --include-pattern\n"; echo -e $help_line; exit 0;
fi

if [ "$EXCCHRPAT" = "NoNaMe" ]; then
 echo -e "\nMissing required argument --exclude-pattern\n"; echo -e $help_line;  exit 0;
fi

if [ "$COLFILE" = "NoNaMe" ]; then
 echo -e "\nMissing required argument --color-file\n"; echo -e $help_line;  exit 0;
fi

if [ "$OUTFIGURE" = "NoNaMe" ]; then
 echo -e "\nMissing required argument --prefix\n"; echo -e $help_line;  exit 0;
fi

# Obtaining multifasta file name
FileN=$(basename ${Ref})

# Obtaining script working directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Step1: Copying assembly and indexing
if [ $(grep 1 <<<$STEPS) ]; then
 echo "Step1: Copying assembly and indexing"
 mkdir -p ${Dir}
 cp $Ref ${Dir}/.
 bwa index ${Dir}/${FileN}
fi


#Step2: Mapping reads against assembly and removing unmapped reads
if [ $(grep 2 <<<$STEPS) ]; then
 echo "Step2: Mapping reads against assembly and removing unmapped reads"
 for j in ${ReadFolder}/*
 do
  bank=$(basename ${j} | cut -f 1 -d ".")
  bwa mem -t 1 -M ${Dir}/${FileN} ${j} | samtools view -bS -uF 4 - > ${Dir}/${OUTFIGURE}_${bank}.bam
 done
fi

#Step3: Counting read number. Multiple hit reads were removed as well as secondary alignments. Exact matches were scored 1 while other matches were scored 0.
if [ $(grep 3 <<<$STEPS) ]; then
 echo "Step3: Counting read number. Multiple hit reads were removed as well as secondary alignments. Exact matches were scored 1 while other matches were scored 0."
 for j in ${ReadFolder}/*
 do
  bank=$(basename ${j} | cut -f 1 -d ".")
  samtools view -h -q 2 ${Dir}/${OUTFIGURE}_${bank}.bam | samtools view -S -F 256 - | grep -w 'NM:i:0' | awk '($6 !~ /D/ && $6 !~ /I/ && $6 !~ /S/ && $6 !~ /H/)' | cut -f 3,4 | sed 's/$/\t1/' > ${Dir}/${OUTFIGURE}_${bank}.list ; samtools view -h -q 2 ${Dir}/${OUTFIGURE}_${bank}.bam | samtools view -S -F 256 - | grep -w 'NM:i:0' | awk '($6 ~ /D/ || $6 ~ /I/ || $6 ~ /S/ || $6 ~ /H/)' | cut -f 3,4 | sed 's/$/\t0/' >> ${Dir}/${OUTFIGURE}_${bank}.list ; samtools view -h -q 2 ${Dir}/${OUTFIGURE}_${bank}.bam | samtools view -S -F 256 - | grep -v 'NM:i:0' | cut -f 3,4 | sed 's/$/\t0/' >> ${Dir}/${OUTFIGURE}_${bank}.list
 done
fi

# Step4: Calculating the number of hits on sliding windows
if [ $(grep 4 <<<$STEPS) ]; then
 echo "Step4: Calculating the number of hits on sliding windows"
 for j in ${ReadFolder}/*
 do
  bank=$(basename ${j} | cut -f 1 -d ".")
  python $SCRIPT_DIR/calcul_pileup_count.py -p ${Dir}/${OUTFIGURE}_${bank}.list -w ${WIN} -f ${Dir}/${FileN} -o ${Dir}/${OUTFIGURE}_${bank}_${WIN}_count -d n
  python $SCRIPT_DIR/calcul_pileup_mean.py -p ${Dir}/${OUTFIGURE}_${bank}.list -w ${WIN} -f ${Dir}/${FileN} -o ${Dir}/${OUTFIGURE}_${bank}_${WIN}_mean -d n
 done
fi

# Step 5: Excluding chromosomes and drawing figure
if [ $(grep 5 <<<$STEPS) ]; then
 echo "Step 5: Excluding chromosomes and drawing figure"
 MOT=START
 for j in ${ReadFolder}/*
 do
  bank=$(basename ${j} | cut -f 1 -d ".")
  inp=${Dir}/${OUTFIGURE}_${bank}_${WIN}_count.tab.gz
  out=${Dir}/${OUTFIGURE}_${bank}_${WIN}_count.chrSelect.tab.gz
  zgrep -E ${CHRPAT} ${inp} | grep -E -v ${EXCCHRPAT} | gzip > ${out}
  if [ "${MOT}" = "START" ]; then
   MOT=${bank}":"${out}
  else
   MOT=$MOT","${bank}":"${out}
  fi
 done
 python $SCRIPT_DIR/DrawStackedDensity.py -f ${MOT} -g png -p ${Dir}/${OUTFIGURE}_count_${WIN} -d gs -C ${COLFILE}
 
 MOT=START
 for j in ${ReadFolder}/*
 do
  bank=$(basename ${j} | cut -f 1 -d ".")
  inp=${Dir}/${OUTFIGURE}_${bank}_${WIN}_mean.tab.gz
  out=${Dir}/${OUTFIGURE}_${bank}_${WIN}_mean.chrSelect.tab.gz
  zgrep -E ${CHRPAT} ${inp} | grep -E -v ${EXCCHRPAT} | gzip > ${out}
  if [ "${MOT}" = "START" ]; then
   MOT=${bank}":"${out}
  else
   MOT=$MOT","${bank}":"${out}
  fi
 done
 python $SCRIPT_DIR/DrawStackedDensity.py -f ${MOT} -g png -p ${Dir}/${OUTFIGURE}_mean_${WIN} -d gs -C ${COLFILE}
fi

# Step 6: Attributing origin based on majority rule on count and drawing

if [ $(grep 6 <<<$STEPS) ]; then
 echo "Step 6: Attributing origin based on majority rule on count and drawing"
 MOT=START
 out=${Dir}/${OUTFIGURE}_${WIN}_count.MajR.tab
 for j in ${ReadFolder}/*
 do
  bank=$(basename ${j} | cut -f 1 -d ".")
  inp=${Dir}/${OUTFIGURE}_${bank}_${WIN}_count.chrSelect.tab.gz
  if [ "${MOT}" = "START" ]; then
   MOT=${bank}":"${inp}
  else
   MOT=$MOT","${bank}":"${inp}
  fi
 done
 python $SCRIPT_DIR/MajorityRuleAttribution.py -f ${MOT} -c ${COLFILE} -o ${out} -m 10 -t 0.6
 awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' ${Dir}/${FileN} | paste - - | sed 's/>//' | grep -E ${CHRPAT} | grep -E -v ${EXCCHRPAT} > ${Dir}/${OUTFIGURE}_${WIN}.chrom
 perl $SCRIPT_DIR/Draw_ideograms_v3_Chr-Sorted.pl ${out} ${Dir}/${OUTFIGURE}_${WIN}.chrom ${Dir}/${OUTFIGURE}_${WIN}_count.MajR 15000 png
fi

# Step 7: Data formating for GEMO analysis
if [ $(grep 7 <<<$STEPS) ]; then
 echo "Step 7: Data formating for GEMO analysis"
 echo "chr len centromereInf centromereSup label" | sed 's: :\t:g'> ${Dir}/${OUTFIGURE}_${WIN}_GEMO_chrom.tab
 awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' ${Dir}/${FileN} | paste - - | sed 's/>//' | grep chr | grep -v random | while read ligne; do
  CHR=$(echo ${ligne} | cut -f 1 -d " ")
  LEN=$(echo ${ligne} | cut -f 2 -d " ")
  echo $CHR $LEN $(($LEN / 2)) $(($LEN / 2 + 2)) A | sed 's: :\t:g' >> ${Dir}/${OUTFIGURE}_${WIN}_GEMO_chrom.tab
 done
 echo "chr start end" | sed 's: :\t:g' > ${Dir}/${OUTFIGURE}_${WIN}_GEMO_header.tab
 > ${Dir}/${OUTFIGURE}_${WIN}_GEMO_body.tab
 i=0
 for j in ${ReadFolder}/*
 do
  bank=$(basename ${j} | cut -f 1 -d ".")
  if [ "$i" = 0 ]; then 
   echo ${ass}
   sed -i 's:$:\t'${bank}':' ${Dir}/${OUTFIGURE}_${WIN}_GEMO_header.tab 
   zcat ${Dir}/${OUTFIGURE}_${bank}_${WIN}_mean.chrSelect.tab.gz | cut -f 1,2,3,4 >> ${Dir}/${OUTFIGURE}_${WIN}_GEMO_body.tab
  else
   sed -i 's:$:\t'${bank}':' ${Dir}/${OUTFIGURE}_${WIN}_GEMO_header.tab 
   zcat ${Dir}/${OUTFIGURE}_${bank}_${WIN}_mean.chrSelect.tab.gz | cut -f 4 > ${Dir}/${OUTFIGURE}_${WIN}_GEMO_body-toto.tab
   paste ${Dir}/${OUTFIGURE}_${WIN}_GEMO_body.tab ${Dir}/${OUTFIGURE}_${WIN}_GEMO_body-toto.tab > ${Dir}/${OUTFIGURE}_${WIN}_GEMO_body-tutu.tab
   mv ${Dir}/${OUTFIGURE}_${WIN}_GEMO_body-tutu.tab ${Dir}/${OUTFIGURE}_${WIN}_GEMO_body.tab
   rm ${Dir}/${OUTFIGURE}_${WIN}_GEMO_body-toto.tab
  fi
  ((i++))
 done
 cat ${Dir}/${OUTFIGURE}_${WIN}_GEMO_header.tab ${Dir}/${OUTFIGURE}_${WIN}_GEMO_body.tab > ${Dir}/${OUTFIGURE}_${WIN}_GEMO_value.tab
 rm ${Dir}/${OUTFIGURE}_${WIN}_GEMO_header.tab
 rm ${Dir}/${OUTFIGURE}_${WIN}_GEMO_body.tab
 echo "group name hex" | sed 's: :\t:g' > ${Dir}/${OUTFIGURE}_${WIN}_GEMO_color.tab
 cat ${COLFILE} | while read line
 do
 NAME=$(echo ${line} | awk '{print $1}')
 RED=$(echo ${line} | awk '{print $2}')
 GREEN=$(echo ${line} | awk '{print $3}')
 BLUE=$(echo ${line} | awk '{print $4}')
 printf "${NAME}\t${NAME}\t#%02x%02x%02x\n" ${RED} ${GREEN} ${BLUE} >> ${Dir}/${OUTFIGURE}_${WIN}_GEMO_color.tab
 done
fi

