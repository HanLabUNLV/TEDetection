#!/usr/bin/env python
import sys
import os.path

def main(args):
  patID = args[1]
  cancerext = args[2]
  normalext = args[3]
  srange = int(args[4])
  minscsup = int(args[5])
  polymorphfilename = args[6]

  cancerRBP = []
  cancerBP = {}
  normalRBP = []
  normalBP = {}
  polymorphBP = {}

  cancerBPfilename = patID + cancerext + ".breakpoints.txt"
  normalBPfilename = patID + normalext + ".breakpoints.txt"
  cancerRBPfilename = patID + cancerext + ".refined.breakpoints.txt"
  normalRBPfilename = patID + normalext + ".refined.breakpoints.txt"

  with open(cancerRBPfilename, 'r') as cancerRBPfile:
    header = cancerRBPfile.readline()
    for line in cancerRBPfile:
      line_sp = line.rstrip('\n').split('\t')
      #if (line_sp[4] != "NA" and line_sp[5] != "NA" and line_sp[-1] == "Yes" and int(line_sp[2]) >= minscsup \
      if (line_sp[-1] == "Yes" and int(line_sp[2]) >= minscsup \
          and (line_sp[3] == "SINE1/7SL" or line_sp[3] == "L1")):
        cancerRBP.append(line)

  with open(normalRBPfilename, 'r') as normalRBPfile:
    header = normalRBPfile.readline()
    for line in normalRBPfile:
      line_sp = line.rstrip('\n').split('\t')
      #if (line_sp[4] != "NA" and line_sp[5] != "NA" and line_sp[-1] == "Yes" and int(line_sp[2]) >= minscsup \
      if (line_sp[-1] == "Yes" and int(line_sp[2]) >= minscsup \
          and (line_sp[3] == "SINE1/7SL" or line_sp[3] == "L1")):
        normalRBP.append(line)

  with open(cancerBPfilename, 'r') as cancerBPfile:
    header = cancerBPfile.readline()
    for line in cancerBPfile:
      line_sp = line.rstrip('\n').split('\t')
      chrom = line_sp[0]
      for i in xrange(int(line_sp[2]) - srange, int(line_sp[2]) + srange):
        cancerBP[(chrom, i)] = 1 # value doesn't matter, just adding key to hash table

  with open(normalBPfilename, 'r') as normalBPfile:
    header = normalBPfile.readline()
    for line in normalBPfile:
      line_sp = line.rstrip('\n').split('\t')
      chrom = line_sp[0]
      for i in xrange(int(line_sp[2]) - srange, int(line_sp[2]) + srange):
        normalBP[(chrom, i)] = 1

  if (os.path.isfile(polymorphfilename)):
    with open(polymorphfilename, 'r') as polymorphfile:
      for line in polymorphfile:
        line_sp = line.rstrip('\n').split('\t')
        chrom = line_sp[0]
        polymorphBP[(chrom, line_sp[4])] = 1
        polymorphBP[(chrom, line_sp[5])] = 1

  overlapfile = open("%s.overlaps.txt" % (patID), 'w+')
  overlapfile.write("Chromosome\tCluster\tSupportingReads\tTEFamily\tLeftBP\tRightBP\tTEMatch\n")
  normalonlyfile = open("%s.normalonly.txt" % (patID), 'w+')
  normalonlyfile.write("Chromosome\tCluster\tSupportingReads\tTEFamily\tLeftBP\tRightBP\tTEMatch\n")
  canceronlyfile = open("%s.canceronly.txt" % (patID), 'w+')
  canceronlyfile.write("Chromosome\tCluster\tSupportingReads\tTEFamily\tLeftBP\tRightBP\tTEMatch\n")

  overlapfiletowrite = {}
  appendpolymorphfile = open(polymorphfilename, 'a+')

  for i in xrange(len(cancerRBP)):
    line_sp = cancerRBP[i].rstrip('\n').split('\t')
    chrom = line_sp[0]
    if (line_sp[4] != "NA" and ((chrom, int(line_sp[4])) in normalBP or (chrom, line_sp[4]) in polymorphBP)):
      overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = cancerRBP[i]
      if (not (chrom, line_sp[4]) in polymorphBP or not (chrom, line_sp[5]) in polymorphBP):
        appendpolymorphfile.write(cancerRBP[i])
        polymorphBP[(chrom, line_sp[4])] = 1
        polymorphBP[(chrom, line_sp[5])] = 1
    elif (line_sp[5] != "NA" and ((chrom, int(line_sp[5])) in normalBP or (chrom, line_sp[5]) in polymorphBP)):
      overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = cancerRBP[i]
      if (not (chrom, line_sp[4]) in polymorphBP or not (chrom, line_sp[5]) in polymorphBP):
        appendpolymorphfile.write(cancerRBP[i])
        polymorphBP[(chrom, line_sp[4])] = 1
        polymorphBP[(chrom, line_sp[5])] = 1
    else:
      canceronlyfile.write(cancerRBP[i])

  for i in xrange(len(normalRBP)):
    line_sp = normalRBP[i].rstrip('\n').split('\t')
    chrom = line_sp[0]
    if (line_sp[4] != "NA" and ((chrom, int(line_sp[4])) in cancerBP or (chrom, line_sp[4]) in polymorphBP)):
      overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = cancerRBP[i]
      if (not (chrom, line_sp[4]) in polymorphBP or not (chrom, line_sp[5]) in polymorphBP):
        appendpolymorphfile.write(cancerRBP[i])
        polymorphBP[(chrom, line_sp[4])] = 1
        polymorphBP[(chrom, line_sp[5])] = 1
    elif (line_sp[5] != "NA" and ((chrom, int(line_sp[5])) in cancerBP or (chrom, line_sp[5]) in polymorphBP)):
      overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = cancerRBP[i]
      if (not (chrom, line_sp[4]) in polymorphBP or not (chrom, line_sp[5]) in polymorphBP):
        appendpolymorphfile.write(cancerRBP[i])
        polymorphBP[(chrom, line_sp[4])] = 1
        polymorphBP[(chrom, line_sp[5])] = 1
    else:
      normalonlyfile.write(cancerRBP[i])

  for key, value in overlapfiletowrite.iteritems():
    overlapfile.write(value)

  overlapfile.close()
  normalonlyfile.close()
  canceronlyfile.close()
  appendpolymorphfile.close()

if (__name__ == "__main__"):
  main(sys.argv)
