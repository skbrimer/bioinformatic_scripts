# simple script for histograms of any file modifiied from Biopython tutorial

from Bio import SeqIO
import pylab
import sys

# read in file
seqs = sys.argv[1]
_fileType = seqs.split(".")

# binning the lengths
sizes = [len(rec) for rec in SeqIO.parse(seqs,_fileType[1])]
# number of bins, set at 20% of the number of seqs
numberOfbins = len(sizes)*0.2

# making plot
pylab.hist(sizes, bins=int(numberOfbins))
pylab.title("%i in file\nLengths from %i to %i"\
		%(len(sizes),min(sizes),max(sizes)))
pylab.xlabel("Sequence length (bp)")
pylab.ylabel("Count")
pylab.show()

