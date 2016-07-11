#!/bin/bash

# Copyright 2015 Jeff Vierstra

set -o pipefail

interval_file=/home/jvierstra/proj/ftd/results.slurm/regions.bed
posterior_file=/home/jvierstra/proj/ftd/results.slurm/posteriors/interval.all.bed.gz
output_dir=/home/jvierstra/proj/ftd/results.slurm/posteriors-mi

chunksize=2500

prefix="fps"

####

mkdir -p ${output_dir}/logs

# load modules
module load bedops/2.4.3
module load python/2.7.11
module load gcc/5.3.0

cat <<__SCRIPT__ > ${output_dir}/slurm.compute_mutual_information_chunk
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=32G
#SBATCH --cpus-per-task=12

set -e -o pipefail

INPUT_FILES=(\`cat ${output_dir}/slurm.compute_mutual_information_chunk.params | head -n \${SLURM_ARRAY_TASK_ID} | tail -n 1\`)

ftd-compute-mutual-information \
	--processors \${SLURM_CPUS_PER_TASK} \
	\${INPUT_FILES[0]} ${posterior_file} \
> \${INPUT_FILES[0]}.out
__SCRIPT__

cat ${interval_file} | split -l ${chunksize} -a 4 -d - ${output_dir}/interval.
ls ${output_dir}/interval.* > ${output_dir}/slurm.compute_mutual_information_chunk.params

JOB0=$(sbatch --export=ALL \
	--job-name=${prefix}.compute_mutual_information_chunk \
	--array=1-$(wc -l < ${output_dir}/slurm.compute_mutual_information_chunk.params) \
	${output_dir}/slurm.compute_mutual_information_chunk)

echo $JOB0

cat <<__SCRIPT__ > ${output_dir}/slurm.compute_mutual_information_merge_sort
#!/bin/bash
#
#SBATCH --output=${output_dir}/logs/%J.out
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -e -o pipefail

TMPDIR=/tmp/\$SLURM_JOB_ID
mkdir -p \${TMPDIR}

INPUT_FILES=(\`cat ${output_dir}/slurm.compute_mutual_information_chunk.params\`)
OUTPUT_FILES=("\${INPUT_FILES[@]/%/.out}")

sort -k1,1 -k2,2n -S 32G \${OUTPUT_FILES[@]} > \${TMPDIR}/interval.all.bed

bgzip -c \${TMPDIR}/interval.all.bed > ${output_dir}/interval.all.bed.gz
tabix -0 -p bed ${output_dir}/interval.all.bed.gz

rm -rf \${TMPDIR}
__SCRIPT__

JOB1=$(sbatch --export=ALL \
	--job-name=${prefix}.compute_mutual_information_merge_sort \
	--depend=afterok:${JOB0##* }\
	${output_dir}/slurm.compute_mutual_information_merge_sort)
