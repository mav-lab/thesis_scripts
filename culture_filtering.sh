#!/bin/bash
echo "Reference genome?"
read -e ref_genome
echo "File pattern?"
read -e output_name

echo "Initializing..."

gatk="/home/mav/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar"
AF="0.0"

#Filtering
java -jar $gatk FilterMutectCalls -R $ref_genome -V $output_name.somatic.vcf.gz --min-allele-fraction $AF -O $AF"VAF.filtered.somatic.vcf.gz"

zcat $AF"VAF.filtered.somatic.vcf.gz" | grep '#' > $output_name"_"$AF"VAF.filtered.somatic.vcf" ; zcat $AF"VAF.filtered.somatic.vcf.gz" \
|  grep -v '#' | awk '$7=="PASS"' >>$output_name"_"$AF"VAF.filtered.somatic.vcf"


#Culture filtering

cat $output_name"_"$AF"VAF.filtered.somatic.vcf" | perl /home/mav/software/cutAF.pl 0.3 > $output_name"_"$AF"_0.3VAF.filtered.somatic.vcf" 
cultureMut=$output_name"_"$AF"_0.3VAF.filtered.somatic.vcf"
#Extract number of SNVs and Indels
SNVs=$(cat $cultureMut | grep -v '#' | awk 'length($4) == 1 && length($5) == 1' | wc -l)
Del_1_3bp=$(cat $cultureMut | grep -v '#' | awk 'length($4) > length($5) && (length($4)-length($5)) <= 3' | wc -l)
Del_3bp=$(cat $cultureMut | grep -v '#' | awk 'length($4) > length($5) && (length($4)-length($5)) > 3' | wc -l)
Ins_1_3bp=$(cat $cultureMut | grep -v '#' | awk 'length($5) > length($4) && (length($5)-length($4)) <= 3' | wc -l)
Ins_3bp=$(cat $cultureMut | grep -v '#' | awk 'length($5) > length($4) && (length($5)-length($4)) > 3' | wc -l)
Indels=$(($Del_1_3bp+$Del_3bp+$Ins_1_3bp+$Ins_3bp))

echo -e "SNVs\t$SNVs\nIndels\t$Indels\nDel1-3bp\t$Del_1_3bp\nDel>3bp\t$Del_3bp\nIns1-3bp\t$Ins_1_3bp\nIns>3bp\t$Ins_3bp" > $output_name"_"$AF"_0.3VAF.SNV.Indels.summary.txt"
