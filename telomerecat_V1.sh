#!/bin/bash

echo "Welcome to telomerecat pipeline"

echo "Somatic bam?"
read -e somatic_bam
echo "Number of cores?"
read  cores
echo "Output name?"
read -e output_name
echo "Initializing..."

telomerecat=telomerecat


#Telomere length measurement
echo "#!/bin/bash" > $output_name.tel.sh
echo "$telomerecat bam2length -p $cores -v 2 $somatic_bam" >>  $output_name.tel.sh

chmod 755 $output_name.tel.sh
