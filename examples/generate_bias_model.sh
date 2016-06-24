#!/bin/bash

# Copyright 2016 Jeff Vierstra
#

#bam_filepath="/home/jvierstra/proj/dnase-perspective/cleavage_model/reads.filtered.bam"
#genome_mappability_filepath=
#fasta_filepath="/home/jvierstra/data/genomes/hg19/hg.ribo.all.fa"

bam_filters="-q 1 -F 512"
filtered_contigs="chrX,chrY,chrM"
max_mem="16G"

tmpdir=$(mktemp -u)

usage() {
	echo "Generate a cleavage bias model using a reference dataset and a genome mappability file"
	echo -e "\nUsage: $0 [options] bam_file mappability_file fasta_file bias_model_file"
}


TEMP=`getopt -o f:m:b:t:h --long filtered-contigs:,max-mem:,bam-filters:,temporary-dir:,help -n 'test.sh' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -f|--filtered-contigs)
            case "$2" in
                "") shift 2;;
                *) filtered_contigs=$2 ; shift 2 ;;
            esac ;;
        -m|--max-mem) 
			case "$2" in
				"") shift 2;;
				*) max_mem=$2 ; shift 2;;
			esac ;;
        -b|--bam-arguments)
            case "$2" in
                "") shift 2;;
                *) bam_filters=$2 ; shift 2 ;;
            esac ;;
        -t|--temporary-dir)
			case "$2" in
				"") shift 2;;
				*) tmpdir=$2; shift 2;;
			esac ;;
		-h|--help)
			usage; exit 0 ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

if [ "$#" -ne 4 ]; then
	echo -e "Error: Not enough arguments!\n"
	usage
	exit 1;
fi

bam_filepath=$1
genome_mappability_filepath=$2
fasta_filepath=$3
model_filepath=$4

rm -rf ${tmpdir} && mkdir -p ${tmpdir}

#
echo -e -n "+ Creating a cleavage file from BAM file..."
#
start=$(date +%s.%N)
#
samtools view ${bam_args} ${bam_filepath} \
| awk -v OFS="\t" '{ \
	if (and($2, 16)) { \
		start = $4 + length($10) - 1; \
		strand = "-"; \
	} else { \
		start = $4; \
		strand = "+"; \
	} \
	print $3, start, start + 1, ".", ".", strand; \
}' \
| grep -v -E $(echo ${filtered_contigs} | tr "," "|") \
| sort-bed --max-mem ${max_mem} - \
| tee ${tmpdir}/cuts.bed \
| bedops -m - | bedops --chop 1 - \
> ${tmpdir}/positions.bed
# 
end=$(date +%s.%N)
dur=$(echo "$end-$start" | bc)
#
echo "Done! ($dur secs)"

echo -e -n "+ Computing sequence context for each cleavage..."
#
start=$(date +%s.%N)
#
mkfifo ${tmpdir}/pos.bed && awk '$6 == "+" { print; }' ${tmpdir}/cuts.bed > ${tmpdir}/pos.bed &
mkfifo ${tmpdir}/neg.bed && awk '$6 == "-" { print; }' ${tmpdir}/cuts.bed > ${tmpdir}/neg.bed &
#
cat > ${tmpdir}/observed_context.py <<__EOF__
from __future__ import print_function
import sys
import pyfaidx

fasta = pyfaidx.Fasta("${fasta_filepath}", sequence_always_upper=True)

for line in sys.stdin:
	fields = line.strip().split('\t')
	chrom = fields[0]
	start = int(fields[1])
	end = int(fields[2])
	
	try:
		print("%s\t%d\t%d\t%s\t%s" % (chrom, start, end, fasta[chrom][start-3:end+2], -fasta[chrom][start-2:end+3]))
	except:
		pass
__EOF__
#
cat ${tmpdir}/positions.bed \
| python ${tmpdir}/observed_context.py \
| bedmap --faster --ec --delim "\t" --echo --count - ${tmpdir}/pos.bed \
| bedmap --faster --ec --delim "\t" --echo --count - ${tmpdir}/neg.bed \
| awk -v OFS="\t" ' \
	{ cnts[$4] += $6; cnts[$5] += $7; } \
	END { \
		for(k in cnts) { \
			print k, cnts[k]; \
		} \
	}' \
| sort -k1,1 - \
> ${tmpdir}/observed.hexamers.txt
#
end=$(date +%s.%N)
dur=$(echo "$end-$start" | bc)
#
echo "Done! ($dur secs)"

#
echo -e -n "+ Computing background sequence context..."
#
start=$(date +%s.%N)
#
cat > ${tmpdir}/expected_context.py <<__EOF__
from __future__ import print_function
import sys
import pyfaidx

fasta = pyfaidx.Fasta("${fasta_filepath}", sequence_always_upper=True)

cnts = {}

for line in sys.stdin:
	fields = line.strip().split('\t')
	chrom = fields[0]
	start = int(fields[1])
	end = int(fields[2])
	strand = fields[5]

	try:
		region = fasta[chrom][start-3:end+2] if strand == '+' else -fasta[chrom][start-2:end+3]
	
		for i in range(3, len(region)-2):
			kmer = region.seq[i-3:i+3]
			cnts[kmer] = cnts.get(kmer, 0) + 1
	except:
		pass	

for key, value in cnts.items():
	print("%s\t%d" % (key, value))
__EOF__

cat ${genome_mappability_filepath} \
| grep -v -E $(echo ${filtered_contigs} | tr "," "|") \
| python ${tmpdir}/expected_context.py \
| sort -k1,1 - \
> ${tmpdir}/expected.hexamers.txt
#
end=$(date +%s.%N)
dur=$(echo "$end-$start" | bc)
#
echo "Done! ($dur secs)"

echo -e -n "+ Completing..."

join -j 1 ${tmpdir}/observed.hexamers.txt ${tmpdir}/expected.hexamers.txt \
| tr " " "\t" | grep -v "N" \
| awk -v OFS="\t" '{ print $0, $2/$3; }' \
> ${model_filepath}

echo "Done!"
