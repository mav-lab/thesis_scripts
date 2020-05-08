#!/bin/bash

echo "Welcome to cnv calling pipeline. We recommend using a gatk docker container to run this script"
#Poner algo de que es necesario procesar antes un intervalo de regiones y un panel de normales
echo "WARNING: You are using ? file as a panel of normals"
echo "Reference genome?"
read -e ref_genome
echo "Normal bam?"
read -e normal_bam_path
echo "Normal path?"
read -e normal_path
echo "Normal panel count?"
read -e count_panel_of_normals
echo "Somatic bam?"
read -e somatic_bam_path
echo "Output name?"
read -e output_name

echo "Initializing..."

gatk="/home/mav/software/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar"
preprocessed_interval_list=/home/mav/genomes/mm10/female_cnv_calling_regions.mmus.38.preprocessed.interval_list
snp=/home/mav/genomes/mm10/snp_mm10_new.vcf
lower="0.9"
upper="1.1"
dict=/home/mav/genomes/mm10/mmus.38.female.chrs.noWG.dict
mContigL="61431566"

#CollectReadCounts

java -jar $gatk CollectReadCounts -I $somatic_bam_path -L $preprocessed_interval_list --interval-merging-rule OVERLAPPING_ONLY \
-O $output_name".cnv_calling_regions.counts.hdf5"

#DenoiseReadCounts

java -Xmx15G -jar $gatk DenoiseReadCounts -I $output_name".cnv_calling_regions.counts.hdf5" --count-panel-of-normals $count_panel_of_normals \
--standardized-copy-ratios $output_name".cnv_calling_regions.clean.standardizedCR.tsv" --denoised-copy-ratios $output_name".cnv_calling_regions.clean.denoisedCR.tsv"

#CollectAllelicCounts (somatic)

java -Xmx15G -jar $gatk CollectAllelicCounts -L $snp -I $somatic_bam_path -R $ref_genome -O $output_name".cnv_calling_regions.allelicCounts.tsv"

#CollectAllelicCounts (normal)

file1=$normal_path"/normal.cnv_calling_regions.allelicCounts.tsv"

if [[ -e "$file1" ]]
then
echo -e "Normal allelic count exists"
else
java -Xmx15G -jar $gatk CollectAllelicCounts -L $snp -I $normal_bam_path -R $ref_genome -O $file1
fi

#ModelSegments

java -Xmx30G -jar $gatk ModelSegments --denoised-copy-ratios $output_name".cnv_calling_regions.clean.denoisedCR.tsv" --allelic-counts $output_name".cnv_calling_regions.allelicCounts.tsv" \
--normal-allelic-counts $file1 --output . --output-prefix $output_name".cnv_calling_regions"

#Filtering

segfile=$output_name"_"$lower"_"$upper.called.seg

java -jar $gatk CallCopyRatioSegments --input $output_name".cnv_calling_regions.cr.seg" --output $segfile \
--neutral-segment-copy-ratio-lower-bound $lower --neutral-segment-copy-ratio-upper-bound $upper

Amp1=$(awk '$6 == "+" && ($3-$2) >= 0 && ($3-$2) < 2500000' $segfile | wc -l)
Amp2=$(awk '$6 == "+" && ($3-$2) >= 2500000 && ($3-$2) < 5000000' $segfile | wc -l)
Amp3=$(awk '$6 == "+" && ($3-$2) >= 5000000 && ($3-$2) < 7500000' $segfile | wc -l)
Amp4=$(awk '$6 == "+" && ($3-$2) >= 7500000 && ($3-$2) < 10000000' $segfile | wc -l)
Amp5=$(awk '$6 == "+" && ($3-$2) >= 10000000' $segfile | wc -l)

Del1=$(awk '$6 == "-" && ($3-$2) >= 0 && ($3-$2) < 25000' $segfile | wc -l)
Del2=$(awk '$6 == "-" && ($3-$2) >= 25000 && ($3-$2) < 50000' $segfile | wc -l)
Del3=$(awk '$6 == "-" && ($3-$2) >= 50000 && ($3-$2) < 75000' $segfile | wc -l)
Del4=$(awk '$6 == "-" && ($3-$2) >= 75000 && ($3-$2) < 100000' $segfile | wc -l)
Del5=$(awk '$6 == "-" && ($3-$2) >= 100000' $segfile | wc -l)

echo -e "Amp0_2500kb\t$Amp1\nAmp2500_5000kb\t$Amp2\nAmp5000_7500kb\t$Amp3\nAmp7500_10000kb\t$Amp4\nAmp>10000kb\t$Amp5\nDel0_25kb\t$Del1\nDel25_50kb\t$Del2\nDel50_75kb\t$Del3\nDel75_100kb\t$Del4\nDel>100kb\t$Del5" > $output_name"_"$lower"_"$upper".called.seg.summary.txt"

#Extract coordinates (biomart)

awk '$6 == "+"' $segfile | awk '{print $1":"$2":"$3":+1";print $1":"$2":"$3":-1"}' > $output_name"_"$lower"_"$upper"_gene_gain_regions.biomart.coords"

awk '$6 == "-"' $segfile | awk '{print $1":"$2":"$3":+1";print $1":"$2":"$3":-1"}' > $output_name"_"$lower"_"$upper"_gene_lost_regions.biomart.coords"

#Extract lenghts (violin plot)

echo -e "#SAMPLE\t#CALL\t#LENGTH" > $output_name"_"$lower"_"$upper.violin.plot.length

awk '$6 == "+"' $segfile | awk -v output_name="$output_name" '{print output_name"\t+\t"($3-$2)}' >> $output_name"_"$lower"_"$upper.violin.plot.length

awk '$6 == "-"' $segfile | awk -v output_name="$output_name" '{print output_name"\t-\t"($3-$2)}' >> $output_name"_"$lower"_"$upper.violin.plot.length

#Plot Modeled Segments

java -jar $gatk PlotModeledSegments --denoised-copy-ratios $output_name".cnv_calling_regions.clean.denoisedCR.tsv" --allelic-counts $output_name".cnv_calling_regions.hets.tsv" \
--segments $output_name".cnv_calling_regions.modelFinal.seg" --sequence-dictionary $dict --minimum-contig-length $mContigL --output . --output-prefix $output_name"_"$lower"_"$upper

#Callable stats

file2=$output_name"_callable_status.bed"

if [[ -e "$file2" ]]
then
echo -e "Callable file  exists"
else
java -jar /home/mav/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T CallableLoci -R $ref_genome \
-I $somatic_bam_path -summary $output_name"_callable_table.txt" -o $file2
fi


#Scan microsatellites - MSI scoring

msi=/home/mav/software/msisensor/msisensor
interval_list="/home/mav/genomes/mm10/female_cnv_calling_regions.mmus.38.bed"
file3="/home/mav/genomes/mm10/mmus.38.microsatellites.default.settings.list"

$msi msi -d $file3 -n $normal_bam_path -t $somatic_bam_path -e $interval_list -o $output_name".msi"
