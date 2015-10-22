import sys
 
sys.path.append('/home/jvierstra/proj/code/footprint-tools')

#import footprint_tools
from footprint_tools import modeling, stats, genomic_interval, cutcounts

import numpy as np
import pysam

intervals_file = sys.argv[1]

tabix = pysam.TabixFile("/home/jvierstra/proj/dnase-perspective/test/interval.chisq.bedgraph.gz")

intervals = genomic_interval.genomic_interval_set(intervals_file)

def parse(interval, tabix_handle):

	x = np.zeros(len(interval))

	for entry in tabix_handle.fetch(interval.chrom, interval.start, interval.end, parser = pysam.asTuple()):

		i = int(entry[1]) - interval.start

		x[i] = np.float(entry[3])

	return x

data = np.vstack( [parse(interval, tabix) for interval in intervals] )

prof = np.mean(data, axis = 0)

vals = np.max(data, axis = 1)

order = np.argsort(vals)[::-1]

print len(prof)

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

from pylab import rcParams

rcParams['pdf.fonttype'] = 42
rcParams['figure.figsize'] = 3, 10

#plt.plot(prof)
plt.imshow(data[order,], aspect = "auto", cmap = "bwr", vmin = 0, vmax = 15, interpolation = "bicubic")

plt.savefig("out.pdf")
