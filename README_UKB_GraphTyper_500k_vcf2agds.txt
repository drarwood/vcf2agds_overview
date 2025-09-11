##########################################################
# UK Biobank GraphTyper 500k Whole-Genome Sequencing Data
# vcf2agds Pipeline on UKB RAP (DNAnexus Platform)
# Xihao Li, Andrew R. Wood, Yuxin Yuan,
# Manrui Zhang, Yushu Huang, Gareth Hawkes,
# Robin N. Beaumont, Michael N. Weedon,
# Wenyuan Li, Xiaoyu Li, Xihong Lin, Zilin Li
##########################################################

##### Install DNAnexus Platform SDK
# Open anaconda prompt:
pip3 install dxpy
#pip3 install --upgrade dxpy

##### Log in
cd /path_to_dx_toolkit/
source dx-toolkit-0.381.0/environment
cd dx-toolkit-0.381.0/bin
dx login
#dx ls (here we assume to have a project named UKB_500k_WGS)


##### Step 0 (Processing each chromosome over 500 jobs)
[link: https://github.com/drarwood/vcf_to_gds_overview/blob/master/get_vcf_file_list_by_chr_and_split.sh]

#!/bin/bash

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


### Move all the generated files from local (dx-toolkit/bin) to an online subfolder (UKB_500k_WGS:/UKB_500k_WGS_vcf_list) on DNAnexus


##### Step 1: Trimming down data in the VCFs (vcf_trimmer)
# Clone this github repo to a local directory:
git clone https://github.com/drarwood/vcf_trimmer

# Navigate to a relevant directory within the project directory on the DNAnexus platform
dx cd UKB_500k_WGS:/

# Now you are ready to build and upload the applet to the DNAnexus platform directory
#dx build -f vcf_trimmer

# Create a new folder under the project directory (Step1_vcf_trimmed_500k)

for CHR in {1..22} 
do
  for i in {1..500}
  do
    dx run UKB_500k_WGS:/vcf_trimmer \
      -ivcf_file_list=UKB_500k_WGS:/UKB_500k_WGS_vcf_list/chr${CHR}_vcf_list_${i} \
      -ifile_label=trimmed \
      -ioutput_dir=Step1_vcf_trimmed_500k \
      -iqc_thresholds="INFO/AAScore>=0.5" \
      -ifields_to_remove="FORMAT/FT,FORMAT/AD,FORMAT/MD,FORMAT/DP,FORMAT/RA,FORMAT/PP,FORMAT/GQ,FORMAT/PL" \
      --instance-type="mem2_ssd1_v2_x32"\
      -y
  done
done


##### Step 2: Merging VCFs (vcf_merger)
# Clone this github repo to a local directory:
git clone https://github.com/drarwood/vcf_merger

# Navigate to a relevant directory within the project directory on the DNAnexus platform
dx cd UKB_500k_WGS:/

# Now you are ready to build and upload the applet to the DNAnexus platform directory
dx build -f vcf_merger

for CHR in {1..22} 
do
  dx ls "UKB_500k_WGS:/Step1_vcf_trimmed_500k/ukb23374_c${CHR}_b*_v1_trimmed.vcf.gz" | sort -t"b" -k3.1 -n | awk -v d="Step1_vcf_trimmed_500k/" '{print d $0}' > chr${CHR}_vcfs_to_merge
done

# Put the chr${CHR}_vcfs_to_merge files to a local folder (UKB_500k_WGS_trimmed_vcf_list)

# Before step 2, need to check files in UKB_500k_WGS_trimmed_vcf_list/chr${CHR}_vcfs_to_merge. If there are duplicates, we need to remove the duplicated files and regenerate the chr${CHR}_vcfs_to_merge files

-------------------------------------------
--- R scripts to check duplicated files ---
-------------------------------------------

files <- read.table("/path_to_dx_toolkit/bin/UKB_500k_WGS_trimmed_vcf_list/chr1_vcfs_to_merge")
head(files)
dim(files)

files[duplicated(files$V1),]
-------------------------------------------

### Move all the generated files from local (dx-toolkit/bin) to an online subfolder (UKB_500k_WGS:/UKB_500k_WGS_trimmed_vcf_list) on DNAnexus

# Create a new folder under the project directory (Step2_vcf_merged_500k)

dx cd UKB_500k_WGS:/Step2_vcf_merged_500k

for CHR in {1..8}
do
  dx run UKB_500k_WGS:/vcf_merger \
    -ivcf_file_list=UKB_500k_WGS:/UKB_500k_WGS_trimmed_vcf_list/chr${CHR}_vcfs_to_merge \
    -imerged_vcf_filename=chr${CHR}_merged.vcf.gz \
    --priority high \
    --instance-type="mem2_ssd2_v2_x16"\
    -y
done

for CHR in {9..22}
do
  dx run UKB_500k_WGS:/vcf_merger \
    -ivcf_file_list=UKB_500k_WGS:/UKB_500k_WGS_trimmed_vcf_list/chr${CHR}_vcfs_to_merge \
    -imerged_vcf_filename=chr${CHR}_merged.vcf.gz \
    --priority high \
    --instance-type="mem2_ssd2_v2_x8"\
    -y
done


##### Step 3: Converting VCF to SeqArray GDS for subsequent annotation in the FAVOR database (vcf2gds)
# Clone this github repo to a local directory:
git clone https://github.com/drarwood/vcf2gds

# Navigate to a relevant directory within the project directory on the DNAnexus platform
dx cd UKB_500k_WGS:/

# Now you are ready to build and upload the applet to the DNAnexus platform directory
dx build -f vcf2gds

# Create a new folder under the project directory (Step3_gds_500k)

dx cd UKB_500k_WGS:/Step3_gds_500k

for CHR in {1..2} {4..4} {6..6} {8..8}
do
  dx run UKB_500k_WGS:/vcf2gds \
    -ivcf_file=UKB_500k_WGS:/Step2_vcf_merged_500k/chr${CHR}_merged.vcf.gz \
    -igds_filename=chr${CHR}.gds \
    -iparallel=30 \
    --priority high \
    --instance-type="mem1_ssd1_v2_x36"\
    -y
done

for CHR in {3..3} {5..5} {7..7} {9..22}
do
  dx run UKB_500k_WGS:/vcf2gds \
    -ivcf_file=UKB_500k_WGS:/Step2_vcf_merged_500k/chr${CHR}_merged.vcf.gz \
    -igds_filename=chr${CHR}.gds \
    -iparallel=16 \
    --priority high \
    --instance-type="mem2_ssd1_v2_x16"\
    -y
done


##### Step 4: Annotating GDS to aGDS for subsequent association analysis in the STAARpipeline (favorannotator)
# Clone this github repo to some directory:
git clone https://github.com/li-lab-genetics/favorannotator-rap.git

# Navigate to a relevant directory within the project directory on the DNAnexus platform
dx cd UKB_500k_WGS:/

# Compile the source code:
dx build -f favorannotator-rap

# Create a new folder under the project directory (UKB_500k_WGS_aGDS)
for CHR in {1..22}
do
  dx run UKB_500k_WGS:/favorannotator \
  -igds_file=UKB_500k_WGS:/Step3_gds_500k/chr${CHR}.gds \
  -ichromosome=${CHR} \
  -iuse_compression=NO \
  -ioutfile=ukb.500k.wgs.chr${CHR}.pass.annotated \
  --destination=UKB_500k_WGS:/UKB_500k_WGS_aGDS --yes
done

