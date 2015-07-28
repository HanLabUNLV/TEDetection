#!/usr/bin/env python

def main(args):
  # Insert size file must contain median and std
  # Assumes bam file is sorted by query name!!
  import pysam
  import sys
  import os.path
  import itertools

  ## Magic Numbers
  minqual = 5 # minimum quality to filter out multi-mapping reads
  Ychromnum = 24 # chromosome number for Y, not looking at any other chromosomes after Y

  bamfile = pysam.Samfile(args[1], "rb")
  discfilename = "Scratch/" + os.path.splitext(os.path.basename(args[1]))[0]+".disc.bam"
  discfile = pysam.Samfile(discfilename, "wb", template=bamfile)
  isfile = open(args[2], "r")
  ismean = float(isfile.readline().strip('\n'))
  isstd = float(isfile.readline().strip('\n'))
  isfile.close()

  discindexfile = open("Scratch/" + os.path.splitext(os.path.basename(args[1]))[0] + ".disc.index", 'w')

  bam_iter = bamfile.__iter__()
  first_read = bam_iter.next()
  for second_read in bam_iter:
    if (first_read.qname == second_read.qname):
      if (first_read.is_paired and second_read.is_paired):
        if (first_read.mapq > minqual or second_read.mapq > minqual):
          if (first_read.tid < Ychromnum and second_read.tid < Ychromnum): # Only look at chr 1-22,X,Y
            if ((first_read.tid != second_read.tid) \
                or (first_read.is_unmapped and not second_read.is_unmapped) \
                or (not first_read.is_unmapped and second_read.is_unmapped) \
                or abs(first_read.tlen) > (ismean + 3*isstd)):
              #read1pos = basepos
              offset1 = discfile.write(first_read)

              #read2pos = basepos + offset1
              offset2 = discfile.write(second_read)

              #discindexfile.write("%s\t%i\t%i\n" % (first_read.qname, first_read.is_read1, read1pos))
              #discindexfile.write("%s\t%i\t%i\n" % (second_read.qname, second_read.is_read1, read2pos))

      try:
        first_read = bam_iter.next()
      except StopIteration:
        continue
        # End of file
    else:
      # Error in file sort
      first_read = second_read

  bamfile.close()
  discfile.close()

  discfile = pysam.Samfile(discfilename, "rb")
  disc_iter = discfile.__iter__()
  readpos = discfile.tell()
  for read in disc_iter:
    discindexfile.write("%s\t%i\t%i\n" % (read.qname, read.is_read1, readpos))
    readpos = discfile.tell()

  discindexfile.close()

if (__name__ == "__main__"):
  main(sys.argv)
