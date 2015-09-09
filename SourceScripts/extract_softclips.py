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
  avereaddepth = int(args[3])
  TEReffilename = args[4]
  numthreads = args[5]
  minclustsize = int(args[6])
  bpfile = open("Results/" + basename + ".breakpoints.txt", 'w')
  bpfile.write("Chromosome\tCluster\tPosition\tSupportingReads\tNonsupportingReads\tSideOfSoftclip\n")
  skippedclusters = open("Results/" + basename + ".skipped.clusters.txt", 'w')
  scfilename = "Scratch/" + basename + ".softclips.fasta"
  scindexfilename = "Scratch/" + basename + ".softclips.sai"
  scsamfilename = "Scratch/" + basename + ".softclips.sam"
  scfile = open(scfilename, 'w')

  ## Magic Numbers
  readlen = 100 # estimated read length for search range
  minqual = 5 # minimum quality
  minsearchrange = 500 # minimum range to search from clusters
  softclipID = 4 # pysam ID for softclipped region in cigar string
  minphredqual = 30.0 # minimum phred quality for softclipped sequences
  phredscaleoffset = 33 # offset to calculate phred score of a nucleotide
  minsoftcliplen = 10 # minimum length of softclipped sequence to be considered
  minnumsoftclips = 5 # minimum number of softclipped reads need to support a breakpoint
                      # This is different than the user input minimum support
                      # User input is used in filtering a later step

  leftscseqs = {}
  rightscseqs = {}

  with open(rangefilename, 'r') as rangefile:
    header = rangefile.readline() # first line is a header
    for line in rangefile:
      leftscseqs.clear()
      rightscseqs.clear()

      ranges = line.split('\t')
      clustsize = int(ranges[4])
      if (clustsize < minclustsize):
        continue
      chrom = ranges[0]
      clusnum = ranges[1]
      rstart = int(ranges[2]) - readlen
      rend = int(ranges[3]) + readlen
      searchrange = rend - rstart
      if (searchrange < minsearchrange):
        rstart -= (minsearchrange - searchrange)/2
        rend += (minsearchrange - searchrange)/2

      rstart = max(rstart, 1) # make sure start range is not negative

      maxreaddepthmult = ((searchrange*2)/readlen) * 5 # read depth multiplier to filter out large repeat regions
                                                       # will filter anything with >5x average read depth
      maxdepth = avereaddepth*maxreaddepthmult

      range_iter = bamfile.fetch(chrom, rstart - 1, rend - 1).__iter__()

      count = 0
      largedepth = False
      for read in range_iter:
        if (not read.is_duplicate and read.mapq >= minqual):
          mismatches = [tag for tag in read.tags if tag[0] == 'NM']
          if (len(mismatches) == 0 or mismatches[0][1] <= 2): # ignore reads with more than 2 mismatches in unique region or continue of no NM tag is found
            if (read.cigar[0][0] == softclipID):
              sclen = read.cigar[0][1]
              scpos = read.pos
              scqual = [ord(read.qual[i]) - phredscaleoffset for i in range(0, sclen)]
              if (float(sum(scqual)/len(scqual)) > minphredqual): # check phred quality of softclipped region
                if (len(read.seq[:sclen]) >= minsoftcliplen):
                  if (scpos in leftscseqs):
                    leftscseqs[scpos].append(read.seq[:sclen])
                  else:
                    leftscseqs[scpos] = []
                    leftscseqs[scpos].append(read.seq[:sclen])
            if (read.cigar[-1][0] == softclipID):
              sclen = read.cigar[-1][1]
              scpos = read.pos + read.rlen - sclen + 1
              if (read.cigar[0][0] == softclipID): # samfiles don't count leftside softclips in the position
                scpos -= read.cigar[0][1]
              scqual = [ord(read.qual[i]) - phredscaleoffset for i in range(-sclen, 0)]
              if (float(sum(scqual)/len(scqual)) > minphredqual): # check phred quality of softclipped region
                if (len(read.seq[-sclen:]) >= minsoftcliplen):
                  if (scpos in rightscseqs):
                    rightscseqs[scpos].append(read.seq[-sclen:])
                  else:
                    rightscseqs[scpos] = []
                    rightscseqs[scpos].append(read.seq[-sclen:])
        count += 1
        if (count >= maxdepth):
          largedepth = True
          break

      if (largedepth):
        skippedclusters.write("%s:%i-%i\t%s\n" % (chrom, rstart, rend, clusnum))
        continue

      leftmaxseqs = 0
      leftscpos = 0
      for pos, seqs in leftscseqs.iteritems():
        if (len(seqs) > leftmaxseqs):
          leftmaxseqs = len(seqs)
          leftscpos = pos

      # Add reads in +-5 buffer range to the max
      for i in xrange(leftscpos - 5, leftscpos + 5):
        if (i != leftscpos and i in leftscseqs):
          for j in xrange(len(leftscseqs[i])):
            leftscseqs[leftscpos].append(leftscseqs[i][j])
            leftmaxseqs += 1

      rightmaxseqs = 0
      rightscpos = 0
      for pos, seqs in rightscseqs.iteritems():
        if (len(seqs) > rightmaxseqs):
          rightmaxseqs = len(seqs)
          rightscpos = pos

      for i in xrange(rightscpos - 5, rightscpos + 5):
        if (i != rightscpos and i in rightscseqs):
          for j in xrange(len(rightscseqs[i])):
            rightscseqs[rightscpos].append(rightscseqs[i][j])
            rightmaxseqs += 1

      maxseqs = leftmaxseqs + rightmaxseqs
      if (maxseqs >= minnumsoftclips):
        coverage = 0
        if (leftmaxseqs > 0):
          coverage_iter = bamfile.fetch(chrom, max(leftscpos-2, 0), leftscpos).__iter__() # positions are 0-based, -2 to +0 is -1 to +1
          for read in coverage_iter:
            if (len(read.cigarstring.split('M')) == 2 and read.mapq > minqual): # Count unique reads mapped through the breakpoint
              coverage += 1
        if (rightmaxseqs > 0):
          coverage_iter = bamfile.fetch(chrom, max(rightscpos-2, 0), rightscpos).__iter__()
          for read in coverage_iter:
            if (len(read.cigarstring.split('M')) == 2 and read.mapq > minqual): # Count unique reads mapped through the breakpoint
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
