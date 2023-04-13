# VCF to GDS Conversion
### Description of VCF to GDS conversion on the UK Biobank RAP for STAAR processing

Due to the vcf file sizes on the DNAnexus platform (especially WGS), it may be necessary to reduce the amount of data contained within VCFs if wanting to merge all VCFs associated with a chromosome on a workstation for subsequent STAAR annotation and processing. Some applets are provided here that may help facilitate this process.

### Step 1: Trimming down data in the VCFs
See [here](https://github.com/drarwood/vcf_trimmer) for applet that removes fields and performs required filtering through bcftools.
This applet takes in a list of files as an input to process. This applet could be used across multiple jobs submitted on the DNAnexus platform which would require unique lists of VCFs to be split. 
#### Example: Processing chromosome 17 over 100 jobs
##### Generating the input file lists:
If you wanted to process the 200K WGS data release for chromosome 17 then `vcf_trimmer` would require a list of VCFs associated with that chromosome. 
Furthermore, you may also want to submit 100 jobs whereby the list of VCFs associated with chromosome 17 is split into 100 unique and equally sized VCF lists.
Running the `get_vcf_file_list_by_chr_and_split.sh` bash script in this repo will produce all the relevant input files for chromosome 17 (`chr17_vcf_list_[1:100]`)
as well as for all other chromosomes for the 60,648 VCFs currently available.
This script has been provided as a starting point that can be amended to suit your requirements.
Note, the file lists generated by this bash script will need to be subsequently uploaded to a project folder on the DNAnexus platform.
##### Submitting the jobs
Once you have a set of 100 files listing the VCFs, you can run the `vcf_trimmer` with each of the files generated above. 
For example, if we wanted to run `vcf_trimmer` on chromosome 17 to remove `FORMAT` fields except `GT` and keep only variants with `AAscore >= 0.5`, 
then the command would look something like this (note priority set to high):

```
for i in {1..100}
do
  dx run /path/to/vcf_trimmer \
    -ivcf_file_list=/path/to/chr17_vcf_list_${i} \
    -ifile_label=trimmed \
    -ioutput_dir=/path/to/output/dir \
    -iqc_thresholds="INFO/AAScore>=0.5" \
    -ifields_to_remove="FORMAT/FT,FORMAT/AD,FORMAT/MD,FORMAT/DP,FORMAT/RA,FORMAT/PP,FORMAT/GQ,FORMAT/PL" \
    --priority high \
    -y
done
```

### Step 2: Merging VCFs
See [here](https://github.com/drarwood/vcf_merger). Another applet wrapper for bcftools.

### Step 3: Converting VCF to GDS for subsequent annotation/analysis in the STAAR pipeline.
See [here](https://github.com/drarwood/vcf2gds). This applet comes with an R library that will be unpacked during runtime and used for data conversion.
