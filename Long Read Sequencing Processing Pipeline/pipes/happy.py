#---------------------
## Run Hap.py (Optional Benchmarking against WGS of the same sample)
## -> The *.summary.csv will tell how effective this tool is for calling
##    SNP and INDEL when comapred to WGS (Illumina short read w/ GATK)
#---------------------

## This 'Truth set' is specific to SVI_0034
TRU_DIR="/data/svi/WGS_SVI/sample_vcfs"
TRU="/data/svi/WGS_SVI/sample_vcfs/${NWGC_id}_sorted.vcf.gz"
TRU_i="${TRU}.tbi"

## Set BED directory
BED_DIR="/data/ref/ConfidentRegions"

## Reference Bed file for chr20
#BED="HG003_GRCh38_chr20_v4.2.1_benchmark_noinconsistent.bed"
BED="HG001_GRCh38_1_22_v4.2.1_benchmark.bed"


## Output VCF file from Pepper-Margin-DeepVariant
OUTPUT_VCF="${PREFIX}.phased.vcf.gz"

docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-v "${TRU_DIR}":"${TRU_DIR}" \
-v "${REF_DIR}":"${REF_DIR}" \
-v "${BED_DIR}":"${BED_DIR}" \
--user $(id -u):$(id -g) \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${TRU} \
${OUTPUT_DIR}/${OUTPUT_VCF} \
-f "${BED_DIR}/${BED}" \
-r "${REF_DIR}/${HGREF}" \
-o "${OUTPUT_DIR}/${PREFIX}" \
--engine=vcfeval \
--threads="${THREADS}" \
--pass-only

