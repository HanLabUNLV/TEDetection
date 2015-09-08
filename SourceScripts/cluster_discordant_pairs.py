#!/usr/bin/env python
import os.path
import sys
import pysam
import itertools

def main(args):
  basename = os.path.splitext(os.path.basename(args[1]))[0]
  discfile = pysam.Samfile(args[1], "rb")
  readdepth = int(args[2])
  isfile = open(args[3], 'r')
  ismean = float(isfile.readline().strip('\n'))
  isfile.close()

  ## Magic Numbers
  minqual = 5 # minimum quality to filter out multi-mapped reads
  readlen = 100 # read length
  bufrange = int(ismean) # buffer range between discordant reads
  # This is to prevent large repeat regions from clogging file writing pipeline
  #maxclustsize = (readdepth*ismean*15)/readlen

  rangesfile = open("Results/" + basename + ".clusters.ranges.txt", 'w')
  rangesfile.write("chromosome\tcluster\tstart\tend\tnum_of_reads\n")
  clusfile = open("Results/" + basename + ".clusters.txt", 'w')

  cluster = {}
  clusternum = 0
  totcluster = 0
  chrom = 0
  bam_iter = discfile.__iter__()

  def writeCluster():
    #if (clusternum >= minclustsize): #and clusternum <= maxclustsize):
    clusfile.write(">%i\n" % (totcluster))

    for value in cluster.itervalues():
      clusfile.write("%s\\%i\n" % (value.qname, value.is_read1))

    rangesfile.write("%s\t%i\t%i\t%i\t%i\n" % (discfile.getrname(chrom), totcluster, leftpos+1, rightpos+1, clusternum))
    return True
    #return False


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
