#!/usr/bin/env python
import sys
import pysam
import glob
import os.path
import itertools
import subprocess

def main(args):
  bamfilename = args[1]
  discfilename = args[2] # disc pair index file
  clusfilename = args[3] # cluster qname file
  TEReffilename = args[4]
  numthreads = args[5]

  clusfastqfilename = "Scratch/" + os.path.splitext(os.path.basename(discfilename))[0] + ".clusters.fastq"
  clusindexfilename = "Scratch/" + os.path.splitext(os.path.basename(discfilename))[0] + ".clusters.sai"
  clussamfilename = "Scratch/" + os.path.splitext(os.path.basename(discfilename))[0] + ".clusters.sam"
  clusfastqfile = open(clusfastqfilename, 'w')

  discseqs = {}

  with open(discfilename, 'r') as discfile:
    for line in discfile:
      line_sp = line.rstrip('\n').split('\t')
      qname = line_sp[0]
      readnum = line_sp[1]
      filepos = line_sp[2]
      discseqs[(qname, readnum)] = int(filepos)

  clusters = {}
  curclus = 0

  with open(clusfilename, 'r') as clusfile:
    for line in clusfile:
      stripline = line.rstrip('\n')
      if (line[0] == ">"):
        curclus = stripline[1:]
        clusters[curclus] = []
      else:
        qname = stripline[0:len(stripline)-2] # -2 to get rid of read number (/1 /2)
        readnum = stripline[-1]
        clusters[curclus].append((qname, readnum))

  bamfile = pysam.Samfile(bamfilename, "rb")

  for key, value in clusters.iteritems():
    for i in xrange(len(value)):
      if (value[i][1] == "0"):
        qname = value[i][0]
        if ((qname, "1") in discseqs):
          try:
            bamfile.seek(discseqs[(qname, "1")])
            mate = bamfile.next()
            # cluster number will be last number after colon in qname
            clusfastqfile.write("@%s:%s\n%s\n+\n%s\n" % (qname, key, mate.seq, mate.qual))
          except StopIteration:
            print qname, discseqs[(qname, "1")]
        else:
          print "Could not find %s mate" % (qname)
      else:
        qname = value[i][0]
        if ((qname, "0") in discseqs):
          try:
            bamfile.seek(discseqs[(qname, "0")])
            mate = bamfile.next()
            # cluster number will be last number after colon in qname
            clusfastqfile.write("@%s:%s\n%s\n+\n%s\n" % (qname, key, mate.seq, mate.qual))
          except StopIteration:
            print qname, discseqs[(qname, "0")]
        else:
          print "Could not find %s mate" % (qname)

  clusfastqfile.close()
  bamfile.close()

  clusalncommand = "bwa aln -t %s %s %s > %s 2>> bwa.err" % (numthreads, TEReffilename, clusfastqfilename, clusindexfilename)
  clusalnproc = subprocess.Popen(clusalncommand, stdout=subprocess.PIPE, shell=True)
  clusalnproc.wait()

  clusmapcommand = "bwa samse %s %s %s > %s 2>> bwa.err" % (TEReffilename, clusindexfilename, clusfastqfilename, clussamfilename)
  clusmapproc = subprocess.Popen(clusmapcommand, stdout=subprocess.PIPE, shell=True)
  clusmapproc.wait()

if (__name__ == "__main__"):
  main(sys.argv)
