###################################################################
# UK Biobank ML-Corrected DRAGEN 500k Whole-Genome Sequencing Data
# vcf2agds Pipeline on UKB RAP (DNAnexus Platform)
# Xihao Li, Andrew R. Wood, Yuxin Yuan,
# Manrui Zhang, Yushu Huang, Gareth Hawkes,
# Robin N. Beaumont, Michael N. Weedon,
# Wenyuan Li, Xiaoyu Li, Xihong Lin, Zilin Li
###################################################################

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

##### Step 2: Merging VCFs (vcf_merger)
# Clone this github repo to a local directory:
git clone https://github.com/drarwood/vcf_merger

# Note: change the timeoutPolicy parameter from 120 (5 days) to 600 (25 days) in the dxapp.json file

# Navigate to a relevant directory within the project directory on the DNAnexus platform
dx cd UKB_500k_WGS:/

# Now you are ready to build and upload the applet to the DNAnexus platform directory
# Note: change the timeoutPolicy parameter from 48 (2 days) to 600 (25 days) in the dxapp.json file
dx build -f vcf_merger

# Get autosomal chromosome vcf file list
for i in {1..22}
do
  DIR_chr="${DIR}chr${i}/"
  # Processing each chromosome
  # Generate file containing full list of VCFs per chromosome:
  # list of vcf.gz files
  dx ls "${DIR_chr}/ukb24311_c${i}_b*_v1.vcf.gz" | sort -t"b" -k3.1 -n | awk -v d="$DIR_chr" '{print d $0}' > chr${i}_vcfs_to_merge
done

# Put the chr${CHR}_vcfs_to_merge files to a local folder (UKB_DRAGEN_500K_WGS_vcf_list)

# Before step 2, need to check files in UKB_DRAGEN_500K_WGS_vcf_list/chr${CHR}_vcfs_to_merge. If there are duplicates, we need to remove the duplicated files and regenerate the chr${CHR}_vcfs_to_merge files

-------------------------------------------
--- R scripts to check duplicated files ---
-------------------------------------------

files <- read.table("~/Dropbox/STAAR/Packages/Analysis_Commons/MacOSM/dx-toolkit-0.381.0/bin/UKB_DRAGEN_500K_WGS_vcf_list/chr1_vcfs_to_merge")
head(files)
dim(files)

files[duplicated(files$V1),]
-------------------------------------------


### Move all the generated files from local (dx-toolkit/bin) to an online subfolder (UKB_500k_WGS:/UKB_DRAGEN_500K_WGS_vcf_list) on DNAnexus

# Create a new folder under the project directory (Step2_vcf_merged_500K)

dx cd UKB_500k_WGS:/Step2_vcf_merged_500K

############
for CHR in {1..4} {7..7} {9..9} {11..12} {18..18}
do
  dx run UKB_500k_WGS:/vcf_merger \
    -ivcf_file_list=UKB_500k_WGS:/UKB_DRAGEN_500K_WGS_vcf_list/chr${CHR}_vcfs_to_merge \
    -imerged_vcf_filename=chr${CHR}_merged.vcf.gz \
    --priority high \
    --instance-type="mem3_ssd2_v2_x16"\
    -y
done

for CHR in {5..6} {8..8} {10..10} {13..17} {19..22}
do
  dx run UKB_500k_WGS:/vcf_merger \
    -ivcf_file_list=UKB_500k_WGS:/UKB_DRAGEN_500K_WGS_vcf_list/chr${CHR}_vcfs_to_merge \
    -imerged_vcf_filename=chr${CHR}_merged.vcf.gz \
    --priority high \
    --instance-type="mem3_ssd2_v2_x8"\
    -y
done

##### Step 3: Converting VCF to GDS for subsequent annotation/analysis in the STAAR pipeline (vcf2gds)
# Clone this github repo to a local directory:
git clone https://github.com/drarwood/vcf2gds

# Navigate to a relevant directory within the project directory on the DNAnexus platform
dx cd UKB_500k_WGS:/

# Now you are ready to build and upload the applet to the DNAnexus platform directory
# Note: change the timeoutPolicy parameter from 72 (3 days) to 600 (25 days) in the dxapp.json file
dx build -f vcf2gds

# Create a new folder under the project directory (Step3_gds_500K)

dx cd UKB_500k_WGS:/Step3_gds_500K

for CHR in {1..2} {7..7} {11..11}
do
  dx run UKB_500k_WGS:/vcf2gds \
    -ivcf_file=UKB_500k_WGS:/Step2_vcf_merged_500K/chr${CHR}_merged.vcf.gz \
    -igds_filename=chr${CHR}.gds \
    -iparallel=32 \
    --priority high \
    --instance-type="mem2_ssd2_v2_x32"\
    -y
done

for CHR in {3..6} {8..10} {12..22}
do
  dx run UKB_500k_WGS:/vcf2gds \
    -ivcf_file=UKB_500k_WGS:/Step2_vcf_merged_500K/chr${CHR}_merged.vcf.gz \
    -igds_filename=chr${CHR}.gds \
    -iparallel=16 \
    --priority high \
    --instance-type="mem2_ssd2_v2_x16"\
    -y
done


##### Step 4: Annotate GDS to AGDS by incorporating variant functional annotations (favorannotator)

# Clone this github repo to some directory:
git clone https://github.com/li-lab-genetics/favorannotator-rap.git

# Navigate to a relevant directory within the project directory on the DNAnexus platform
dx cd UKB_500k_WGS:/

# Compile the source code:
dx build -f favorannotator-rap


### Uncompressed AGDS
# Create a new folder under the project directory (UKB_500K_WGS_AGDS_uncompressed)
for CHR in {1..22}
do
  dx run UKB_500k_WGS:/favorannotator \
  -igds_file=UKB_500k_WGS:/Step3_gds_500K/chr${CHR}.gds \
  -ichromosome=${CHR} \
  -iuse_compression=NO \
  -ioutfile=ukb.500k.wgs.chr${CHR}.pass.annotated \
  --destination=UKB_500k_WGS:/UKB_500K_WGS_AGDS_uncompressed --yes
done
