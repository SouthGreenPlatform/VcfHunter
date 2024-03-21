#!/usr/bin/perl
#
#  Copyright 2024 CIRAD
#
#  Olivier Garsmeur
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

use strict;
use warnings;
use GD::Simple;
use GD::SVG;


my @names = GD::Simple->color_names;
	# foreach my $item (@names){
		# print $item,"\n";
	# }

print "\n-----------------------------------------------------------------------------------------------------------------------------\n";
print "\n#####   usage:     #####"."\n\n";
print "perl Draw_ideogram_v3_Chr-Sorted.pl [ideogram.txt] [chr.txt] [prefix_for_output] [scale=echelle] [output format png or svg] \n\n";
print "[ideogram.txt] = tabular file: marker/chr/start/end/label/color"."\n"."\n"; 
print "[chr.txt] = tabular file : chr/length\n\n";
print "[scale] = image scale (ex: use 50000 for chromosomes or genome. small values increase image size \n\n";
print "output format = png or svg\n\n";
print "with svg format, provide colors as RGB code ex: for red color: 255,0,0","\n","\n";
print "-----------------------------------------------------------------------------------------------------------------------------\n\n";


my $file = $ARGV[0]; #fichier tabule : marker/chr/start/end/label/color trie par chr et start

# ATTTTTCTAAAAATAACCTATCATAGG_323 chr-10_1        19422   45075   S.officinarum   green
# CAGTTGGTTGATCGGACGCGCAGCGCG_442 chr-10_1        46208   46239   undetermined    gray
# ACTCTAACCGGACGCGGGCACTGTCGC_925 chr-10_1        46284   46314   undetermined    gray
# CAGCGCGTCCGGTCCGTTACGGGACCG_303 chr-10_1        46371   62649   S.spontaneum    red
# CCTAGGGTGAGGCCTATAGTTCTACAA_323 chr-10_1        66021   66047   undetermined    gray
# CACCGGACGCGGTCACTGTGCGACCGG_1526        chr-10_1        76687   135845  S.spontaneum    red
# ACAAAACATCACCTCAAGGCTCCCCGT_683 chr-10_1        136450  136480  undetermined    gray
# GCCAACGTCGCTTTCCCGCCTTTTTCC_425 chr-10_1        136870  151514  S.spontaneum    red



my $file2 = $ARGV[1]; #fichier taille des chromosomes

# chr-10-10	67335906
# chr10-2	66740959
#etc ...

my $prefix = $ARGV[2];
my $scale=$ARGV[3]; # passe à l'echelle pour passer des bp à une mesure gérable pour le graphique par exemple on divise par 50000 les coordonnées 


#prefix for the 2output files
my $format_image=$ARGV[4];

my $fileimg;
my $filelegend;

if ($format_image eq "png"){
	$fileimg=$prefix.".png";
	$filelegend=$prefix."-legend.png";
	print "2 files created : $prefix.png and $prefix-legend.png"."\n"."\n";
}

if ($format_image eq "svg"){
	$fileimg=$prefix.".svg";
	$filelegend=$prefix."-legend.svg";
	print "2 files created : $prefix.svg and $prefix-legend.svg"."\n"."\n";
}



#recupere la coordonnee la plus elevee pour dimensionner la hauteur de l'image
#=============================================================================

my $Hscale=0; # taille en cM du groupe le plus grand
my $nb_chromosomes_to_draw=0; #nombre de chromosomes à dessiner
my %chromosomeend; #taille de chaque chromosome;
my $longuest_name_of_chromosome=0;

# Contiendra l'ordre des chromosomes du file2 pour l'affichage
my @order = ();

open(SIZEFILE,"$file2");
while (my $line=<SIZEFILE>){
	
	if ($line=~ /^$/){
	}
	
	else{
	
		my ($chromosome_name,$chromosome_size)=split("\t",$line);

		#Ajout du chromosome pour l'ordre d'affichage
		push(@order,$chromosome_name); 

		my @short_chrom_name=split(" ",$chromosome_name);
		$chromosome_name=$short_chrom_name[0];
			$chromosome_size=$chromosome_size/$scale;
		
		#taille du plus grand chromosome
		if ($Hscale<$chromosome_size){
			$Hscale=$chromosome_size; 
		}
		
		#nombre de chromosomes
		if ($chromosome_name ne ""){
			$nb_chromosomes_to_draw++;
			$chromosomeend{$chromosome_name}=$chromosome_size;
		}
		
		#chromosome dont le nom est le plus long;
		if (length($chromosome_name) >$longuest_name_of_chromosome){
			$longuest_name_of_chromosome=length($chromosome_name);
		}
	}
}	

print "Number of chromosomes to draw = $nb_chromosomes_to_draw","\n";
#print,"Length of longuest chromosome name = $longuest_name_of_chromosome characteres","\n";

# create a new image (width, height) dimension en fonction du nombre de chr et taille des chr
#===========================================================================================
my $largeur; #(1 chromosome = 200 pixel pour 1 chromosome (espace interchromosome inclu)
my $hauteur; #(1 cM = 5 pixel)

my $marge_hauteur=($longuest_name_of_chromosome*20)+150;
my $marge_largeur=100;

$largeur = ($nb_chromosomes_to_draw*($marge_largeur*2))+$marge_largeur;
$hauteur = $Hscale+$marge_hauteur+50;

#create_image;

my $img;

if ($format_image eq "png"){
	$img = GD::Simple->new($largeur, $hauteur);
	#print "Image_size($largeur, $hauteur)","\n\n";
}

if ($format_image eq "svg"){
	GD::Simple->class('GD::SVG');
	$img = GD::Simple->new($largeur, $hauteur);
}


#draw body_of_chromosomes:
#============================

my %chromosome_coordinates;
my $topx_body_chromosome=$marge_largeur;
my $topy_body_chromosome=$marge_hauteur;
my $botx_body_chromosome=$topx_body_chromosome+$marge_largeur;
my $boty_body_chromosome;
my $taille_de_police=40;
my $chromosome_width;
my $xtxt;
my $ytxt;
my $color_text;
my %chromosome_to_draw;


# Parcours la liste des chromosomes en fonction de l'ordre indiqué dans le file2
foreach my $chromosomes_to_include(@order) {
	my $length_chromosome=$chromosomeend{$chromosomes_to_include};
	$boty_body_chromosome=$length_chromosome+$marge_hauteur;
	$chromosome_width=($botx_body_chromosome-$topx_body_chromosome);
	$chromosome_to_draw{$chromosomes_to_include}=$chromosomes_to_include;
	$img->fgcolor(225,225,225);
	$img->bgcolor(225,225,225);
	
	# body of chromosome
	$img->rectangle($topx_body_chromosome,$topy_body_chromosome,$botx_body_chromosome,$boty_body_chromosome);
	$chromosome_coordinates{$chromosomes_to_include}=$topx_body_chromosome."\t".$topy_body_chromosome."\t".$botx_body_chromosome."\t".$boty_body_chromosome;
	
	#dessine des bords arrondis en haut et bas du chromosome
	# $img ->moveTo($topx_body_chromosome+($chromosome_width/2),$topy_body_chromosome);
	# $img->ellipse($chromosome_width,$chromosome_width);
	# $img ->moveTo($botx_body_chromosome-($chromosome_width/2),$boty_body_chromosome);
	# $img->ellipse($chromosome_width,$chromosome_width);	
	
	# inscript nom du marqueur
	#si bord arrondis#   $img->moveTo($topx_body_chromosome+(($botx_body_chromosome-$topx_body_chromosome)/2)+($taille_de_police/2),$topy_body_chromosome-($chromosome_width/2)-10);
	
	if($format_image eq "png"){
		$img->moveTo($topx_body_chromosome+(($botx_body_chromosome-$topx_body_chromosome)/2)+($taille_de_police/2),$topy_body_chromosome-10);
		$img->angle(-90);
		$img->fgcolor('black');
		$img->font('Times');
		$img->fontsize($taille_de_police);
		$img->string($chromosomes_to_include);
	}
	
	if($format_image eq "svg"){
		$img->angle(-90);
		$xtxt=$topx_body_chromosome+(($botx_body_chromosome-$topx_body_chromosome)/2)+($taille_de_police/2);
		$ytxt=$topy_body_chromosome-10;

		$color_text=(0,0,0);
		my ($red,$green,$blue);
		if ($color_text =~/,/) {
			($red,$green,$blue) = (split(/\,/,$color_text));
		}
		
		$img->stringUp(gdLargeFont,$xtxt-25,$ytxt,$chromosomes_to_include,$color_text);
	}

	$topx_body_chromosome=$topx_body_chromosome+($marge_largeur*2);
	$botx_body_chromosome=$topx_body_chromosome+$marge_largeur;
}
	

#place markers on ideogram
#==========================

open (FILE, "$file");

my $topx=0;
my $topy;
my $botx;
my $boty;
my %legend;
my $marker_nb=0;
my $new_start;
my $new_end;
my($R,$G,$B);
my $colorname;


while (my $line = <FILE>){
	
	if ($line=~ /^$/){
	}
	
	else{
	
		if ($line!~ /^#/){
			chomp $line;
			my($marker,$chr,$start,$end,$label,$color)=split("\t",$line);
			
			if (exists($chromosome_to_draw{$chr})){
				
				if ($start > $end){
					$new_start=$end/$scale;
					$new_end=$start/$scale;
				}
				else{
					$new_start=$start/$scale;
					$new_end=$end/$scale;
				}
				
				#au cas ou la couleur ne soit pas defini,
				if (not (defined $color)){
						$color="black";
				}

				#draw markers
				if (exists($chromosome_coordinates{$chr})){

					#coordonnees du chromosome:
					my @coordinates_chrom=split("\t",$chromosome_coordinates{$chr}); 
					
					#place le marker sur le chromosome;
					$topx = $coordinates_chrom[0];
					$topy = $new_start+$coordinates_chrom[1];
					$botx = $coordinates_chrom[2];
					$boty = $new_end+$coordinates_chrom[1];
					
					
					if ($format_image eq "svg"){
						if($color=~/,/){			
							($R,$G,$B)=(split(/\,/,$color));
							$color=$img->colorAllocate($R,$G,$B);
							$img->bgcolor($color);
							$img->fgcolor($color);						
							$img->rectangle($topx,$topy,$botx,$boty);
							
							# for legend-color
							$legend{$label}=$color;
							
						}
						else {
							print "\nError! format svg requires colors as RGB code ex: 255,0,0\n";
							die;
						}
					}
					
					else{
						if($color=~/,/){ #convertit le code RGB en colorname
							($R,$G,$B)=split('\,',$color);
							$img->bgcolor($R,$G,$B);
							$img->fgcolor($R,$G,$B);
							$legend{$label}=$color;
						}
						 else{ # utilise des codes coleurs du style 'black' 'red' 'blue' etc ...
							$img->bgcolor($color);
							$img->fgcolor($color);
							$legend{$label}=$color;
						 }
						$img->rectangle($topx,$topy,$botx,$boty);
						
						# for legend-color
						
					}
					
					$marker_nb++;
				}
				else{
					#print "chromosome $chr does not exists in $file2","\n";
				}
			}
				
			else {
				#ne dessine pas ce chromosome car il n'est pas dans le fichier chr/taille
			}
		}
	}
}

print "Total number of hits: $marker_nb","\n\n";

#creer une image legende
#==========================
my $nbcolor;
$nbcolor =scalar(keys %legend);
my $sizelegend = $nbcolor*22;
my $ctx=10;
my $cty=10;
my $cbx=20;
my $cby=20;
my $key;
my $legend_color_box;
my $imglegend;


if ($format_image eq "svg"){
	GD::Simple->class('GD::SVG');
	$imglegend= GD::Simple->new(200,$sizelegend);
	
	foreach $key (sort(keys %legend)){
		$legend_color_box=$legend{$key};
		($R,$G,$B)=split("\,",$legend_color_box);
		$imglegend->bgcolor($R,$G,$B);
		$imglegend->fgcolor($R,$G,$B);		
		$imglegend->rectangle($ctx,$cty,$cbx,$cby);
		$legend_color_box=(0,0,0);
		#$legend_color_box=split("\,",$legend_color_box);
		$imglegend->string(gdLargeFont,$ctx+40,$cty-4,$key,$legend_color_box);
		$cty=$cty+20;
		$cby=$cby+20;
	}
}


else{
	$imglegend=GD::Simple->new(200,$sizelegend);

	#while (my ($keys, $value)=each(%legend)){
	foreach $key (sort(keys %legend)){
		my $color=$legend{$key};
		
		if($color=~/,/){
		($R,$G,$B)=split('\,',$color);
			$imglegend->fgcolor($R,$G,$B);
			$imglegend->bgcolor($R,$G,$B);
			$imglegend->rectangle($ctx,$cty,$cbx,$cby);
			$imglegend->moveTo($ctx+18,$cty+10);
			$imglegend->fgcolor('black');
			$imglegend->font('Arial');
			$imglegend->fontsize(12);
			$imglegend->string($key.' ('.$color.')');
			$cty=$cty+20;
			$cby=$cby+20;
		}
		else{
			$imglegend->fgcolor($color);
			$imglegend->bgcolor($color);
			$imglegend->rectangle($ctx,$cty,$cbx,$cby);
			$imglegend->moveTo($ctx+18,$cty+10);
			$imglegend->fgcolor('black');
			$imglegend->font('Arial');
			$imglegend->fontsize(12);
			$imglegend->string($key.' ('.$color.')');
			$cty=$cty+20;
			$cby=$cby+20;
		}
	}
}


#creer le fichier ideogram en png ou svg 


if ($format_image eq "png"){
	open my $out, '>', $fileimg or die;
	binmode $out;
	print $out $img->png;
}	


if ($format_image eq "svg"){
	open my $out, '>', $fileimg or die;
	print $out $img->svg;
}


#creer le fichier legend en png ou svg

if ($format_image eq "png"){
	open my $out2, '>', $filelegend or die;
	binmode $out2;
	print $out2 $imglegend->png;
}

if ($format_image eq "svg"){
	open my $out2, '>', $filelegend or die;
	binmode $out2;
	print $out2 $imglegend->svg;
}


system "rm -fr $file.temp*";
