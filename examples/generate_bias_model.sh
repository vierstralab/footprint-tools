#!/bin/bash

# Copyright 2016-21 Jeff Vierstra
# Altius Institute for Biomedical Sciences

mapq=1
filtered_contigs="chrX,chrY,chrM"
max_mem="16G"

left_offset=3
right_offset=2

tmpdir=$(mktemp -u)

usage() {
	echo "Generate a cleavage bias model using a reference dataset and a genome mappability file"
	echo -e "\nUsage: $0 [options] bam_file mappability_file fasta_file bias_model_file"
	echo -e "\n\tOptions:"
}


TEMP=`getopt -o l:r:f:m:b:t:h --long left:,right:,filtered-contigs:,max-mem:,bam-filters:,temporary-dir:,help -n 'test.sh' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
	-l|--left)
	    case "$2" in
		"") shift 2;;
		*) left_offset=$2 ; shift 2 ;;
	    esac ;;
	-r|--right)
	    case "$2" in
		"") shift 2;;
		*) right_offset=$2 ; shift 2 ;;
	    esac ;; 
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
        -q|--mapq)
            case "$2" in
                "") shift 2;;
                *) mapq=$2 ; shift 2 ;;
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

echo -e "+ Temporary directory: ${tmpdir}"
#
echo -e -n "+ Creating a cleavage file from BAM file..."
#
start=$(date +%s.%N)
#
cat > ${tmpdir}/cuts.py <<__EOF__
from __future__ import print_function
import sys
import pysam
import pyfaidx

reads = pysam.AlignmentFile(sys.argv[1], "rb")
read = None

while(1):
	try:
		read = reads.next()
	except:
		break

	if read.mapq < ${mapq}:
		continue

	if read.is_qcfail:
		continue

	if read.is_reverse:
		i = int(read.aend) - 1
	else:
		i = int(read.pos)

	chrom = reads.getrname(read.reference_id)

	print("%s\t%d\t%d\t.\t.\t%s" % (chrom, i, i + 1, '-' if read.is_reverse else '+'))
__EOF__

python ${tmpdir}/cuts.py ${bam_filepath} \
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
awk '$6 == "+" { print; }' ${tmpdir}/cuts.bed > ${tmpdir}/pos.bed &
awk '$6 == "-" { print; }' ${tmpdir}/cuts.bed > ${tmpdir}/neg.bed &

wait; wait;

cat >> ${tmpdir}/observed_context.py <<__EOF__
from __future__ import print_function
import sys
import pyfaidx

fasta = pyfaidx.Fasta("${fasta_filepath}", one_based_attributes=False, sequence_always_upper=True)

loffset = ${left_offset}
roffset = ${right_offset}

for line in sys.stdin:
	fields = line.strip().split('\t')
	chrom = fields[0]
	start = int(fields[1])
	end = int(fields[2])
	
	try:
		pos = fasta[chrom][start-loffset:end+roffset]
		neg = -fasta[chrom][start-roffset:end+loffset]
		print("%s\t%d\t%d\t%s\t%s" % (chrom, start, end, pos, neg))
	except:
		pass
__EOF__

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
> ${tmpdir}/observed.kmers.txt

end=$(date +%s.%N)
dur=$(echo "$end-$start" | bc)

echo "Done! ($dur secs)"

echo -e -n "+ Computing background sequence context..."

start=$(date +%s.%N)

cat > ${tmpdir}/expected_context.py <<__EOF__
from __future__ import print_function
import sys
import pyfaidx

fasta = pyfaidx.Fasta("${fasta_filepath}", one_based_attributes=False, sequence_always_upper=True)

cnts = {}

loffset = ${left_offset}
roffset = ${right_offset}

for line in sys.stdin:
	fields = line.strip().split('\t')
	chrom = fields[0]
	start = int(fields[1])
	end = int(fields[2])
	strand = fields[5]

	try:
		region = fasta[chrom][start-loffset:end+roffset] if strand == '+' else -fasta[chrom][start-roffset:end+loffset]
	
		for i in range(loffset, len(region)-roffset):
			kmer = region.seq[i-loffset:i+roffset+1]
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
> ${tmpdir}/expected.kmers.txt

end=$(date +%s.%N)
dur=$(echo "$end-$start" | bc)

echo "Done! ($dur secs)"

echo -e -n "+ Completing..."

cat > ${tmpdir}/finalize.py <<__EOF__
from __future__ import print_function
import sys
from itertools import product

loffset = ${left_offset}
roffset = ${right_offset}

exp = {}
obs = {}

with open(sys.argv[1]) as f:
	for line in f:
		(kmer, val) = line.strip().split("\t")
		obs[kmer] = int(val)

with open(sys.argv[2]) as f:
	for line in f:
		(kmer, val) = line.strip().split("\t")
		exp[kmer] = int(val)

kmers = set(map(lambda x: ''.join(x), list(product('ACGT', repeat = loffset + roffset + 1))))

for kmer in kmers:
	o = obs.get(kmer, 1)
	e = exp.get(kmer, 1)
	print("%s\t%d\t%d\t%f" % (kmer, o, e, float(o)/float(e)))

__EOF__
 
python ${tmpdir}/finalize.py ${tmpdir}/observed.kmers.txt ${tmpdir}/expected.kmers.txt > ${model_filepath}

echo "Done!"
