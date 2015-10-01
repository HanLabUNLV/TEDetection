#!/usr/bin/env python
import os.path
import sys
import pysam
import itertools

def main(args):
  basename = os.path.splitext(os.path.basename(args[1]))[0]
  discfile = pysam.Samfile(args[1], "rb")
  readdepth = int(args[2])
  minclustsize = int(args[3])
  isfile = open(args[4], 'r')
  ismean = float(isfile.readline().strip('\n'))
  isfile.close()

  ## Magic Numbers
  minqual = 5 # minimum quality to filter out multi-mapped reads
  readlen = 100 # read length
  bufrange = int(ismean) # buffer range between discordant reads
  windowsize = int(ismean) # window size for finding highest density of reads
  # This is to prevent large repeat regions from clogging file writing pipeline
  #maxclustsize = (readdepth*ismean*15)/readlen

  rangesfile = open("Results/" + basename + ".clusters.ranges.txt", 'w')
  rangesfile.write("chromosome\tcluster\tstart\tend\tnum_of_reads\n")
  clusfile = open("Results/" + basename + ".clusters.txt", 'w')

  cluster = {}
  clusternum = 0
  totcluster = 0
  chrom = 0
  leftpos = 0
  rightpos = 0
  bam_iter = discfile.__iter__()

  def writeCluster():
    if (clusternum >= minclustsize): #and clusternum <= maxclustsize):
      clusfile.write(">%i\n" % (totcluster))

      bins = []

      for i in xrange(((rightpos - leftpos) / windowsize) + 1):
        bins.append([])

      for value in cluster.itervalues():
        pos = value.pos
        if (read.cigar[0][0] == 4):
          pos -= read.cigar[0][1]
        index = (pos - leftpos) / windowsize
        bins[index].append(value)

      max_bin = bins.index(max(bins, key=lambda x: len(x))) # get bin with largest number of reads

      count = 0
      for i in xrange(max(max_bin - 2, 0), min(max_bin + 3, len(bins))): # get 2 bins to the left and right of max
        for disc_read in bins[i]:
          clusfile.write("%s\\%i\n" % (disc_read.qname, disc_read.is_read1))
          count += 1

      clusleftpos = leftpos + (max_bin - 2)*windowsize # 2 bins before max
      clusrightpos = clusleftpos + 5*windowsize # 2 bins after max (2 before + 4 + 1 to reach end)

      rangesfile.write("%s\t%i\t%i\t%i\t%i\n" % (discfile.getrname(chrom), totcluster, clusleftpos+1, clusrightpos+1, count))
      return True
    return False


  for read in bam_iter:
    if (read.mapq > 5): # Do not cluster multi-mapping reads
      read_pos = read.pos
      if (read.cigar[0][0] == 4):
        read_pos -= read.cigar[0][1]
      if (read.tid == chrom):
        if (clusternum == 0):
          leftpos = read_pos - bufrange
          rightpos = read_pos + read.rlen + bufrange
          cluster[clusternum] = read
          clusternum += 1
        elif (read_pos > leftpos and read_pos <= rightpos):
          cluster[clusternum] = read
          rightpos = read_pos + read.rlen + bufrange
          clusternum += 1
        else:
          if (writeCluster()):
            totcluster += 1
          clusternum = 0
          cluster.clear()
          cluster[0] = read
          clusternum += 1
          leftpos = read_pos - bufrange
          rightpos = read_pos + read.rlen + bufrange

      else:
        if (writeCluster()):
          totcluster += 1
        clusternum = 0
        cluster.clear()
        cluster[0] = read
        clusternum += 1
        leftpos = read_pos - bufrange
        rightpos = read_pos + read.rlen + bufrange
        chrom = read.tid

  writeCluster()

  discfile.close()
  rangesfile.close()
  clusfile.close()

if (__name__ == "__main__"):
  main(sys.argv)
