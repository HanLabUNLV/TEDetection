#!/usr/bin/env python
import pysam
import os.path
import sys
import itertools
import subprocess
from multiprocessing import Process, Queue

def main(args):
  basename = os.path.splitext(os.path.basename(args[1]))[0]
  rangefilename = args[2]
  avereaddepth = int(args[3])
  TEReffilename = args[4]
  numthreads = args[5]
  minclustsize = int(args[6])
  runaligner = int(args[7])
  minnumsoftclips = int(args[8])
  minphredqual = float(args[9]) # minimum phred quality for softclipped sequences
  minsoftcliplen = int(args[10]) # minimum length of softclipped sequence to be considered
  maxnummismatches = int(args[11]) # max mismatches in read to be considered
  reportallbps = int(args[12])
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
  softclipID = 4 # pysam ID for softclipped region in cigar string
  phredscaleoffset = 33 # offset to calculate phred score of a nucleotide
  #minsoftcliplen = 10 # minimum length of softclipped sequence to be considered
  #minnumsoftclips = 2 # minimum number of softclipped reads need to support a breakpoint
                      # This is different than the user input minimum support
                      # User input is used in filtering a later step

  minsearchrange = avereaddepth # minimum range to search from clusters
  maxsearchrange = avereaddepth*6 # maximum search range from clusters (will use middle range if too large)

  def threadedExtract(range_lines, results_queue, thread_num, bamfile, reportallbps):
    results = ["", "", ""]
    range_lines = range_lines.split('\n')[:-1]
    for line in range_lines:
      leftscseqs = {}
      rightscseqs = {}
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
  #      elif (searchrange > maxsearchrange):
  #        mid = (rend + rstart)/2
  #        rstart = mid - maxsearchrange/2
  #        rend = mid + maxsearchrange/2

      rstart = max(rstart, 1) # make sure start range is not negative
      searchrange = rend - rstart
      maxreaddepthmult = (searchrange/readlen) * 5 # read depth multiplier to filter out large repeat regions
                                                       # will filter anything with >5x average read depth
      maxdepth = avereaddepth*maxreaddepthmult

      range_iter = bamfile.fetch(chrom, rstart - 1, rend - 1).__iter__()

      count = 0
      largedepth = False
      for read in range_iter:
        if (not read.is_duplicate and read.mapq >= minqual):
          mismatches = filter(lambda x: x[0] == "NM", read.tags)
          if (mismatches == [] or mismatches[0][1] <= maxnummismatches): # ignore reads with more than 2 mismatches in unique region or continue of no NM tag is found
            if (read.cigar[0][0] == softclipID):
              sclen = read.cigar[0][1]
              scpos = read.pos
              if (sclen >= minsoftcliplen):
                ave = (sum(map(ord, read.qual[sclen - minsoftcliplen:sclen])) - phredscaleoffset*minsoftcliplen)/minsoftcliplen
                if (ave > minphredqual): # check phred quality of softclipped region
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
              if (sclen >= minsoftcliplen):
                ave = (sum(map(ord, read.qual[-sclen:-sclen+minsoftcliplen])) - phredscaleoffset*minsoftcliplen)/minsoftcliplen
                if (ave > minphredqual): # check phred quality of softclipped region
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
        results[2] += "%s:%i-%i\t%s\n" % (chrom, rstart, rend, clusnum)
        continue

      if (not reportallbps):
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
  #        if (leftmaxseqs > 0):
  #          coverage_iter = bamfile.fetch(chrom, max(leftscpos-2, 0), leftscpos).__iter__() # positions are 0-based, -2 to +0 is -1 to +1
  #          for read in coverage_iter:
  #            if (read.mapq > minqual and len(read.cigarstring.split('M')) == 2): # Count unique reads mapped through the breakpoint
  #              coverage += 1
  #        if (rightmaxseqs > 0):
  #          coverage_iter = bamfile.fetch(chrom, max(rightscpos-2, 0), rightscpos).__iter__()
  #          for read in coverage_iter:
  #            if (read.mapq > minqual and len(read.cigarstring.split('M')) == 2): # Count unique reads mapped through the breakpoint
  #              coverage += 1

          if (leftmaxseqs > 0):
            seqs = leftscseqs[leftscpos]
            for i in xrange(len(seqs)):
              results[0] += ">leftsc%i:%s\n%s\n" % (i, clusnum, seqs[i])
            results[1] += chrom + '\t' + str(clusnum) + '\t' + str(leftscpos) + '\t' + str(leftmaxseqs) + '\t' + str(coverage) + '\tleft\n'

          if (rightmaxseqs > 0):
            seqs = rightscseqs[rightscpos]
            for i in xrange(len(seqs)):
              results[0] += ">rightsc%i:%s\n%s\n" % (i, clusnum, seqs[i])
            results[1] += chrom + '\t' + str(clusnum) + '\t' + str(rightscpos) + '\t' + str(rightmaxseqs) + '\t' + str(coverage) + '\tright\n'
      # if all the breakpoints are to be reported..
      else:
        coverage = 0
        j = 0
        for pos, seqs in leftscseqs.iteritems():
          for i in xrange(len(seqs)):
            results[0] += ">leftsc%i:%s.%i\n%s\n" % (i, clusnum, j, seqs[i])
          results[1] += chrom + '\t' + str(clusnum) + '.' + str(j) + '\t' + str(pos) + '\t' + str(len(seqs)) + '\t' + str(coverage) + '\tleft\n'
          j += 1

        for pos, seqs in rightscseqs.iteritems():
          for i in xrange(len(seqs)):
            results[0] += ">rightsc%i:%s.%i\n%s\n" % (i, clusnum, j, seqs[i])
          results[1] += chrom + '\t' + str(clusnum) + '.' + str(j) + '\t' + str(pos) + '\t' + str(len(seqs)) + '\t' + str(coverage) + '\tright\n'
          j += 1

    results_queue.put(results)

  bamfile_handles = []
  for i in xrange(int(numthreads)):
    bamfile_handles.append(pysam.Samfile(args[1], 'rb'))

  with open(rangefilename, 'r') as rangefile:
    for i, l in enumerate(rangefile):
      pass
  num_lines = i + 1

  with open(rangefilename, 'r') as rangefile:
    header = rangefile.readline() # first line is a header
    eof = False
    line_num = 0
    lines = [""] * int(numthreads)
    for i in xrange(int(numthreads)):
      try:
        for j in xrange(num_lines / (int(numthreads)) + (num_lines % int(numthreads)) / int(numthreads) + 1):
          lines[i] += rangefile.next()
          line_num += 1
      except StopIteration:
        eof = True
        break

    threads = []
    results_queue = Queue()

    #print "Starting threads %i-%i" % (line_num - int(numthreads), line_num)
    for i in xrange(int(numthreads)):
      threads.append(Process(target=threadedExtract, args=(lines[i], results_queue, i, bamfile_handles[i], reportallbps)))
      threads[i].start()

    scfile_to_write = ""
    bpfile_to_write = ""
    skipped_to_write = ""
    for i in xrange(int(numthreads)):
      results = results_queue.get()
      scfile_to_write += results[0]
      bpfile_to_write += results[1]
      skipped_to_write += results[2]

    for i in xrange(int(numthreads)):
      threads[i].join()

    #print "Joined threads %i-%i" % (line_num - int(numthreads), line_num)

    scfile.write(scfile_to_write)
    bpfile.write(bpfile_to_write)
    skippedclusters.write(skipped_to_write)

  bpfile.close()
  for i in xrange(int(numthreads)):
    bamfile_handles[i].close()
  skippedclusters.close()
  scfile.close()

  if (runaligner):
    scalncommand = "bwa aln -t %s %s %s > %s 2>> bwa.err" % (numthreads, TEReffilename, scfilename, scindexfilename)
    scalnproc = subprocess.Popen(scalncommand, stdout=subprocess.PIPE, shell=True)
    scalnproc.wait()

    scmapcommand = "bwa samse %s %s %s > %s 2>> bwa.err" % (TEReffilename, scindexfilename, scfilename, scsamfilename)
    scmapproc = subprocess.Popen(scmapcommand, stdout=subprocess.PIPE, shell=True)
    scmapproc.wait()

if (__name__ == "__main__"):
  main(sys.argv)
