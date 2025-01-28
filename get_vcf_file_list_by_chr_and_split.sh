#!/bin/bash
# Andrew Wood
# No warranty.

# Set some global variables here:
DIR="/Bulk/GATK and GraphTyper WGS/GraphTyper population level WGS variants, pVCF format [500k release]"
SPLITS=500

# Define function to split file where suffix not fixed size
split() {

  # Get chr file count
  FILEC=`wc -l < "$1" | awk '{print $1}'`

  # GET SIZE OF LISTS WE WILL REQUIRE BASED ON N_LINES / N_SPLITS
  LISTCOUNT1=$(($FILEC / $SPLITS))
  LISTCOUNT2=$(($LISTCOUNT1+1))
  # GET NUMBER OF FILES FOR LISTCOUNT1 AND LISTCOUNT2
  COUNT2=`expr "$FILEC" % "$SPLITS" | bc`
  COUNT1=`expr $SPLITS - $COUNT2`

  # CYCLE THROUGH ALL FILES NEEDED TO GENERATE
  OFFSET=1
  for (( i=1; i<=$SPLITS; i++ ))
  do
    if [ $i -le $COUNT1 ]
    then
      tail -n +$OFFSET $1 | head -$LISTCOUNT1 > $1_$i
      OFFSET=$(($OFFSET + $LISTCOUNT1))
    else
      tail -n +$OFFSET $1 | head -$LISTCOUNT2 > $1_$i
      OFFSET=$(($OFFSET + $LISTCOUNT2))
    fi
  done

}

# Get autosomal chromosome vcf file list ordered by block
for i in {1..22}
do
  # Generate file containing full list of VCFs per chromosome, ordered by block:
  dx ls "${DIR}/chr${i}/ukb23374_c${i}_b*_v1.vcf.gz" | sort -t"b" -k3.1 -n | awk -v d="$DIR/chr${i}/" '{print d $0}' > chr${i}_vcf_list
  # Split lists out
  split chr${i}_vcf_list
done

# Get X vcf file list ordered by block
dx ls "${DIR}/chrX/ukb23374_cX_b*_v1.vcf.gz" | sort -t"b" -k3.1 -n | awk -v d="$DIR/chrX/" '{print d $0}' > chrX_vcf_list
# Split lists out
split chrX_vcf_list

