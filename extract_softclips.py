#!/usr/bin/env python
import pysam
import os.path
import sys
import itertools
import subprocess

def main(args):
  bamfile = pysam.Samfile(args[1], 'rb')
  basename = os.path.splitext(os.path.basename(args[1]))[0]
  rangefilename = args[2]
  searchrange = int(args[3])
  avereaddepth = int(args[4])
  TEReffilename = args[5]
  numthreads = args[6]
  bpfile = open("Results/" + basename + ".breakpoints.txt", 'w')
  bpfile.write("Chromosome\tCluster\tPosition\tSupportingReads\tNonsupportingReads\tSideOfSoftclip\n")
  skippedclusters = open("Results/" + basename + ".skipped.clusters.txt", 'w')
  scfilename = "Scratch/" + basename + ".softclips.fasta"
  scindexfilename = "Scratch/" + basename + ".softclips.sai"
  scsamfilename = "Scratch/" + basename + ".softclips.sam"
  scfile = open(scfilename, 'w')

  ## Magic Numbers
  minqual = 5 # minimum quality
  softclipID = 4 # pysam ID for softclipped region in cigar string
  minphredqual = 10 # minimum phred quality for softclipped sequences
  phredscaleoffset = 33 # offset to calculate phred score of a nucleotide
  maxreaddepthmult = 15 # read depth multiplier to filter out large repeat regions
  minsoftcliplen = 7 # minimum length of softclipped sequence to be considered
  minnumsoftclips = 2 # minimum number of softclipped reads need to support a breakpoint
                      # This is different than the user input minimum number of softclipped reads
                      # User input is used in filtering a later step

  leftscseqs = {}
  rightscseqs = {}

  with open(rangefilename, 'r') as rangefile:
    header = rangefile.readline() # first line is a header
    for line in rangefile:
      leftscseqs.clear()
      rightscseqs.clear()

      ranges = line.split('\t')
      chrom = ranges[0]
      clusnum = ranges[1]
      readstrand = ranges[5].rstrip('\n')
      if (readstrand == '+'):
        rstart = int(ranges[3]) - searchrange
        rend = int(ranges[3]) + searchrange
      elif (readstrand == '-'):
        rstart = int(ranges[2]) - searchrange
        rend = int(ranges[2]) + searchrange

      range_iter = bamfile.fetch(chrom, rstart - 1, rend - 1).__iter__()

      count = 0
      largedepth = False
      for read in range_iter:
        if (not read.is_duplicate and read.mapq >= minqual):
          if (read.cigar[0][0] == softclipID):
            sclen = read.cigar[0][1]
            scpos = read.pos
            if (ord(read.qual[sclen-1]) - phredscaleoffset > minphredqual): # check phred quality of softclipped region
              if (scpos in leftscseqs):
                leftscseqs[scpos].append(read.seq[:sclen])
              else:
                leftscseqs[scpos] = []
                leftscseqs[scpos].append(read.seq[:sclen])
          if (read.cigar[-1][0] == softclipID):
            sclen = read.cigar[-1][1]
            scpos = read.pos + read.rlen - sclen + 1
            if (ord(read.qual[-sclen-1]) - phredscaleoffset > minphredqual): # check phred quality of softclipped region
              if (scpos in rightscseqs):
                rightscseqs[scpos].append(read.seq[-sclen:])
              else:
                rightscseqs[scpos] = []
                rightscseqs[scpos].append(read.seq[-sclen:])
        count += 1
        if (count >= avereaddepth*maxreaddepthmult):
          largedepth = True
          break

      if (largedepth):
        skippedclusters.write("%s:%i-%i\t%s\n" % (chrom, rstart, rend, clusnum))
        continue

      longscs = 0
      leftmaxseqs = 0
      leftscpos = 0
      for pos, seqs in leftscseqs.iteritems():
        longscs = 0
        for i in xrange(len(seqs)):
          if (len(seqs[i]) > minsoftcliplen):
            longscs += 1
        if (longscs > leftmaxseqs):
          leftmaxseqs = longscs
          leftscpos = pos

      longscs = 0
      rightmaxseqs = 0
      rightscpos = 0
      for pos, seqs in rightscseqs.iteritems():
        longscs = 0
        for i in xrange(len(seqs)):
          if (len(seqs[i]) > minsoftcliplen):
            longscs += 1
        if (longscs > rightmaxseqs):
          rightmaxseqs = longscs
          rightscpos = pos

      if (leftmaxseqs == 0 and rightmaxseqs == 0):
        continue

      maxseqs = leftmaxseqs + rightmaxseqs
      if (maxseqs >= minnumsoftclips):
        coverage = 0
        if (leftmaxseqs > 0):
          coverage_iter = bamfile.fetch(chrom, leftscpos-2, leftscpos).__iter__() # positions are 0-based, -2 to +0 is -1 to +1
          for read in coverage_iter:
            if (read.cigarstring == "100M" and read.qual > minqual): # Count unique reads mapped through the breakpoint
              coverage += 1
        if (rightmaxseqs > 0):
          coverage_iter = bamfile.fetch(chrom, rightscpos-2, rightscpos).__iter__()
          for read in coverage_iter:
            if (read.cigarstring == "100M" and read.qual > minqual):
              coverage += 1

        if (leftmaxseqs > 0):
          seqs = leftscseqs[leftscpos]
          for i in xrange(len(seqs)):
            if (len(seqs[i]) > minsoftcliplen): # length of the softclip region is greater than 7
              scfile.write(">leftsc%i:%s\n%s\n" % (i, clusnum, seqs[i]))
          bpfile.write(chrom + '\t' + str(ranges[1]) + '\t' + str(leftscpos) + '\t' + str(maxseqs) + '\t' + str(coverage) + '\tleft\n')

        if (rightmaxseqs > 0):
          seqs = rightscseqs[rightscpos]
          for i in xrange(len(seqs)):
            if (len(seqs[i]) > minsoftcliplen): # length of the softclip region is greater than 7
              scfile.write(">rightsc%i:%s\n%s\n" % (i, clusnum, seqs[i]))
          bpfile.write(chrom + '\t' + str(ranges[1]) + '\t' + str(rightscpos) + '\t' + str(maxseqs) + '\t' + str(coverage) + '\tright\n')

  bpfile.close()
  bamfile.close()
  skippedclusters.close()
  scfile.close()

  scalncommand = "bwa aln -t %s %s %s > %s 2>> bwa.err" % (numthreads, TEReffilename, scfilename, scindexfilename)
  scalnproc = subprocess.Popen(scalncommand, stdout=subprocess.PIPE, shell=True)
  scalnproc.wait()

  scmapcommand = "bwa samse %s %s %s > %s 2>> bwa.err" % (TEReffilename, scindexfilename, scfilename, scsamfilename)
  scmapproc = subprocess.Popen(scmapcommand, stdout=subprocess.PIPE, shell=True)
  scmapproc.wait()

if (__name__ == "__main__"):
  main(sys.argv)
