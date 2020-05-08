#!/bin/bash

echo "Welcome to MSIsensor pipeline"
echo "WARNING: You are using female_cnv_calling_regions.mmus.38.bed file to delimit the analysis"
#echo "Reference genome?"
#read -e ref_genome
#echo "Reference genome directory?"
#read -e genome_direc
echo "Normal bam?"
read -e normal_bam_path
echo "Somatic bam?"
read -e somatic_bam_path
echo "Output name?"
read -e output_name
echo "Initializing..."

msi=msisensor
interval_list="/home/mav/genomes/mm10/female_cnv_calling_regions.mmus.38.bed"
file2="/home/mav/genomes/mm10/mmus.38.microsatellites.default.settings.list"

#Scan microsatellites from reference genome_genera el output pero vacio, da un error

#file1=$genome_direc"/microsatellites.default.settings.list"

#if [[ -e "$file1" ]]
#then
#echo -e "Microsatellites list exists"
#else
#$msi scan -d $ref_genome -o $genome_direc"/microsatellites.default.settings.list"
#fi

#MSI scoring

$msi msi -d $file2 -n $normal_bam_path -t $somatic_bam_path -e $interval_list -o $output_name".msi"  
