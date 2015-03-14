#!/usr/bin/env python
import os.path
import sys
import pysam
import itertools

def main(args):
  basename = os.path.splitext(os.path.basename(args[1]))[0]
  discfile = pysam.Samfile(args[1], "rb")
  minclustsize = int(args[2])
  readdepth = int(args[3])
  isfile = open(args[4], 'r')
  ismean = float(isfile.readline().strip('\n'))
  isfile.close()

  ## Magic Numbers
  minqual = 5 # minimum quality to filter out multi-mapped reads
  readlen = 100 # read length
  # This is to prevent large repeat regions from clogging file writing pipeline
  maxclustsize = (readdepth*ismean)/readlen

  rangesfile = open("Results/" + basename + ".clusters.ranges.txt", 'w')
  rangesfile.write("chromosome\tcluster\tstart\tend\tnum_of_reads\tstrand\n")
  clusfile = open("Results/" + basename + ".clusters.txt", 'w')

  poscluster = {}
  negcluster = {}
  posclusternum = 0
  negclusternum = 0
  negtot = 0
  postot = 0
  totcluster = 0
  chrom = 0
  bam_iter = discfile.__iter__()

  for read in bam_iter:
    if (read.mapq > 5): # Do not cluster multi-mapping reads
      if (read.tid == chrom):
        if (read.is_reverse):
          if (negclusternum == 0):
            negleftpos = read.pos
            negrightpos = read.pos + read.rlen
            negcluster[negclusternum] = read
            negclusternum += 1
          elif (read.pos > negleftpos and read.pos <= negrightpos):
            negcluster[negclusternum] = read
            if (read.cigar[-1][0] != 4):
              negrightpos = read.pos + read.rlen
            else:
              negrightpos = read.pos + read.rlen - read.cigar[-1][1]
            negclusternum += 1
          else:
            if (negclusternum >= minclustsize and negclusternum <= maxclustsize):
              if (negtot == postot):
                totcluster += 1
                negtot = totcluster
              clusfile.write(">%i\n" % (negtot))

              for value in negcluster.itervalues():
                clusfile.write("%s\\%i\n" % (value.qname, value.is_read1))

              rangesfile.write("%s\t%i\t%i\t%i\t%i\t%s\n" % (discfile.getrname(chrom), negtot, negleftpos+1, negrightpos+1, negclusternum, '-'))
              totcluster += 1
              negtot = totcluster

            negclusternum = 0
            negcluster.clear()
            negcluster[0] = read
            negclusternum += 1
            negleftpos = read.pos
            negrightpos = read.pos + read.rlen

        else:
          if (posclusternum == 0):
            posleftpos = read.pos
            posrightpos = read.pos + read.rlen
            poscluster[posclusternum] = read
            posclusternum += 1
          elif (read.pos > posleftpos and read.pos <= posrightpos):
            poscluster[posclusternum] = read
            if (read.cigar[-1][0] != 4):
              posrightpos = read.pos + read.rlen
            else:
              posrightpos = read.pos + read.rlen - read.cigar[-1][1]
            posclusternum += 1
          else:
            if (posclusternum >= minclustsize and posclusternum <= maxclustsize):
              if (negtot == postot):
                totcluster += 1
                postot = totcluster
              clusfile.write(">%i\n" % (postot))

              for value in poscluster.itervalues():
                clusfile.write("%s\\%i\n" % (value.qname, value.is_read1))

              rangesfile.write("%s\t%i\t%i\t%i\t%i\t%s\n" % (discfile.getrname(chrom), postot, posleftpos+1, posrightpos+1, posclusternum, '+'))
              totcluster += 1
              postot = totcluster

            posclusternum = 0
            poscluster.clear()
            poscluster[0] = read
            posclusternum += 1
            posleftpos = read.pos
            posrightpos = read.pos + read.rlen

      else:
        if (negclusternum >= minclustsize and negclusternum <= maxclustsize):
          clusfile.write(">%i\n" % (negtot))

          for value in negcluster.itervalues():
            clusfile.write("%s\\%i\n" % (value.qname, value.is_read1))

          rangesfile.write("%s\t%i\t%i\t%i\t%i\t%s\n" % (discfile.getrname(chrom), negtot, negleftpos+1, negrightpos+1, negclusternum+1, '-'))
          totcluster += 1
          negtot = totcluster

        negclusternum = 0
        negcluster.clear()

        if (read.is_reverse):
          negcluster[0] = read
          chrom = read.tid
          negleftpos = read.pos
          negrightpos = read.pos + read.rlen

        if (posclusternum >= minclustsize and posclusternum <= maxclustsize):
          clusfile.write(">%i\n" % (postot))

          for value in poscluster.itervalues():
            clusfile.write("%s\\%i\n" % (value.qname, value.is_read1))

          rangesfile.write("%s\t%i\t%i\t%i\t%i\t%s\n" % (discfile.getrname(chrom), postot, posleftpos+1, posrightpos+1, posclusternum+1, '+'))
          totcluster += 1
          postot = totcluster

        posclusternum = 0
        poscluster.clear()

        if (not read.is_reverse):
          poscluster[0] = read
          chrom = read.tid
          posleftpos = read.pos
          posrightpos = read.pos + read.rlen

  discfile.close()
  rangesfile.close()
  clusfile.close()

if (__name__ == "__main__"):
  main(sys.argv)
