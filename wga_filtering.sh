#!/bin/bash
echo "Reference genome?"
read -e ref_genome
echo "Input name?"
read -e output_name
echo "Minimum allele frequency?"
read AF
echo "Initializing..."

gatk="/home/mav/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar"

#Filtering
java -jar $gatk FilterMutectCalls -R $ref_genome -V $output_name.somatic.vcf.gz --min-allele-fraction $AF -O $AF"VAF.filtered.somatic.vcf.gz"

zcat $AF"VAF.filtered.somatic.vcf.gz" | perl /home/vqf/software/wga.pl > $output_name"_2d3_"$AF"VAF.filtered.somatic.vcf"

#Extract number of SNVs and Indels
SNVs=$(cat $output_name"_2d3_"$AF"VAF.filtered.somatic.vcf" | grep -v '#' | awk 'length($4) == 1 && length($5) == 1' | wc -l)
Del_1_3bp=$(cat $output_name"_2d3_"$AF"VAF.filtered.somatic.vcf" | grep -v '#' | awk 'length($4) > length($5) && (length($4)-length($5)) <= 3' | wc -l)
Del_3bp=$(cat $output_name"_2d3_"$AF"VAF.filtered.somatic.vcf" | grep -v '#' | awk 'length($4) > length($5) && (length($4)-length($5)) > 3' | wc -l)
Ins_1_3bp=$(cat $output_name"_2d3_"$AF"VAF.filtered.somatic.vcf" | grep -v '#' | awk 'length($5) > length($4) && (length($5)-length($4)) <= 3' | wc -l)
Ins_3bp=$(cat $output_name"_2d3_"$AF"VAF.filtered.somatic.vcf" | grep -v '#' | awk 'length($5) > length($4) && (length($5)-length($4)) > 3' | wc -l)
Indels=$(($Del_1_3bp+$Del_3bp+$Ins_1_3bp+$Ins_3bp))

echo -e "SNVs\t$SNVs\nIndels\t$Indels\nDel1-3bp\t$Del_1_3bp\nDel>3bp\t$Del_3bp\nIns1-3bp\t$Ins_1_3bp\nIns>3bp\t$Ins_3bp" > $output_name"_2d3_"$AF"VAF.SNV.Indels.summary.txt"

#Extract coordinates (nVenn)

cat $output_name"_2d3_"$AF"VAF.filtered.somatic.vcf" | grep -v '#' | awk '{print $1":"$2":"$4">"$5}' > $output_name"_2d3_"$AF.coords

#Get AFs

perl /home/mav/software/getAF.pl $output_name"_2d3_"$AF"VAF.filtered.somatic.vcf" > passAFs_2d3_"$AF"

##Rscript /home/mav/software/hist.r passAFs_"$AF" passAFs_"$AF" ##Necesito poder responder a las preguntas primero
