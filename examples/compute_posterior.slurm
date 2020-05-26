#!/bin/bash

# 2018 Jeff Vierstra
# jvierstra@altius.org

set -o pipefail

metadata_file=$1 # Metadata file
interval_file=$2 # File (BED-format) containing regions where footprint detection will occur (i.e., hotspots)
output_dir=$3 # Path to directory where output will be stored (created by script if non-existent)

chunksize=2500

prefix="fps"

####

mkdir -p ${output_dir}/logs

# load modules
module load bedops/2.4.3

cat <<__SCRIPT__ > ${output_dir}/slurm.compute_posterior_chunk
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

set -e -o pipefail

INPUT_FILES=(\`cat ${output_dir}/slurm.compute_posterior_chunk.params | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1\`)

ftd-compute-posterior \
	${base_dir}/metadata.txt \${INPUT_FILES[0]} \
| sort --buffer-size=16G -k1,1 -k2,2n > \${INPUT_FILES[0]}.out
__SCRIPT__

cat ${interval_file} | split -l ${chunksize} -a 4 -d - ${output_dir}/interval.
ls ${output_dir}/interval.* > ${output_dir}/slurm.compute_posterior_chunk.params

JOB0=$(sbatch --export=ALL \
	--job-name=${prefix}.compute_posterior_chunk \
	--array=1-$(wc -l < ${output_dir}/slurm.compute_posterior_chunk.params) \
	${output_dir}/slurm.compute_posterior_chunk)

echo $JOB0

cat <<__SCRIPT__ > ${output_dir}/slurm.compute_posterior_merge_sort
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -e -o pipefail

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

INPUT_FILES=(\`cat ${output_dir}/slurm.compute_posterior_chunk.params\`)
OUTPUT_FILES=("\${INPUT_FILES[@]/%/.out}")

sort -k1,1 -k2,2n -S 32G -m \${OUTPUT_FILES[@]} > \${TMPDIR}/interval.all.bed

bgzip -c \${TMPDIR}/interval.all.bed > ${output_dir}/interval.all.bed.gz
tabix -0 -p bed ${output_dir}/interval.all.bed.gz

rm -rf \${TMPDIR}
__SCRIPT__

JOB1=$(sbatch --export=ALL \
	--job-name=${prefix}.compute_posterior_merge_sort \
	--depend=afterok:${JOB0##* }\
	${output_dir}/slurm.compute_posterior_merge_sort)
