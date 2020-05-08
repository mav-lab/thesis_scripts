#!/bin/bash

echo "Welcome to somatic calling pipeline. We recommend using a gatk docker container to run this script"
echo "WARNING: You are using mouseBM_3refWT.pon.vcf.gz file as a panel of normals"
echo "Reference genome?"
read -e ref_genome
echo "Normal bam?"
read -e normal_bam_path
echo "Normal sample name?(RG)"
read -e normal_ref_name
echo "Somatic bam?"
read -e somatic_bam_path
echo "Output name?"
read -e output_name
echo "Minimum allele frequency?"
read AF
echo "Initializing..."

gatk="/home/mav/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar"


#Somatic calling
java -Xmx6G -jar $gatk Mutect2 -R $ref_genome -I $normal_bam_path -I $somatic_bam_path -normal $normal_ref_name \
--panel-of-normals /home/mav/alig/mouse_single_cell/pon_db/mouseBM_3refWT.pon.vcf.gz -O $output_name.somatic.vcf.gz

#Filtering
java -jar $gatk FilterMutectCalls -R $ref_genome -V $output_name.somatic.vcf.gz --min-allele-fraction $AF -O $AF"VAF.filtered.somatic.vcf.gz"

zcat $AF"VAF.filtered.somatic.vcf.gz" | grep '#' > $output_name"_"$AF"VAF.filtered.somatic.vcf" ; zcat $AF"VAF.filtered.somatic.vcf.gz" \
|  grep -v '#' | awk '$7=="PASS"' >>$output_name"_"$AF"VAF.filtered.somatic.vcf"

#Extract number of SNVs and Indels
SNVs=$(cat $output_name"_"$AF"VAF.filtered.somatic.vcf" | grep -v '#' | awk 'length($4) == 1 && length($5) == 1' | wc -l)
Del_1_3bp=$(cat $output_name"_"$AF"VAF.filtered.somatic.vcf" | grep -v '#' | awk 'length($4) > length($5) && (length($4)-length($5)) <= 3' | wc -l)
Del_3bp=$(cat $output_name"_"$AF"VAF.filtered.somatic.vcf" | grep -v '#' | awk 'length($4) > length($5) && (length($4)-length($5)) > 3' | wc -l)
Ins_1_3bp=$(cat $output_name"_"$AF"VAF.filtered.somatic.vcf" | grep -v '#' | awk 'length($5) > length($4) && (length($5)-length($4)) <= 3' | wc -l)
Ins_3bp=$(cat $output_name"_"$AF"VAF.filtered.somatic.vcf" | grep -v '#' | awk 'length($5) > length($4) && (length($5)-length($4)) > 3' | wc -l)
Indels=$(($Del_1_3bp+$Del_3bp+$Ins_1_3bp+$Ins_3bp))

echo -e "SNVs\t$SNVs\nIndels\t$Indels\nDel1-3bp\t$Del_1_3bp\nDel>3bp\t$Del_3bp\nIns1-3bp\t$Ins_1_3bp\nIns>3bp\t$Ins_3bp" > $output_name"_"$AF"VAF.SNV.Indels.summary.txt"

#Callable stats

java -jar /home/mav/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T CallableLoci -R $ref_genome \
-I $somatic_bam_path -summary $output_name"_callable_table.txt" -o $output_name"_callable_status.bed"

#Get AFs and plot them

perl /home/mav/software/getAF.pl $output_name"_"$AF"VAF.filtered.somatic.vcf" > passAFs_"$AF"

##Rscript /home/mav/software/hist.r passAFs_"$AF" passAFs_"$AF" ##Necesito poder responder a las preguntas primero
