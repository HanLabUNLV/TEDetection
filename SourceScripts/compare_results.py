#!/usr/bin/env python
import sys
import os.path
import os
import datetime

def main(args):
  patID = args[1]
  cancerext = args[2]
  normalext = args[3]
  scsrange = int(args[4])
  discsrange = int(args[5])
  minscsup = int(args[6])
  polymorphfilename = args[7]

  cancerRBP = []
  cancerBP = {}
  normalRBP = []
  normalBP = {}
  polymorphBP = {}

  cancerBPfilename = patID + cancerext + ".breakpoints.txt"
  normalBPfilename = patID + normalext + ".breakpoints.txt"
  cancerRBPfilename = patID + cancerext + ".refined.breakpoints.txt"
  normalRBPfilename = patID + normalext + ".refined.breakpoints.txt"
  cancerCRfilename = patID + cancerext + ".disc.clusters.ranges.txt"
  normalCRfilename = patID + normalext + ".disc.clusters.ranges.txt"

  print("\nReading %s..." % (cancerRBPfilename))
  print(str(datetime.datetime.today()))
  with open(cancerRBPfilename, 'r') as cancerRBPfile:
    header = cancerRBPfile.readline()
    file_size = os.stat(cancerRBPfilename).st_size
    count = 0
    for line in cancerRBPfile:
      count += len(line)
      sys.stdout.write("\r%i%%" % (int((count*100)/file_size)))
      sys.stdout.flush()
      line_sp = line.rstrip('\n').split('\t')
      TE_list = line_sp[3].split(',')
      chrom = line_sp[0]
      if (chrom not in cancerBP):
        cancerBP[chrom] = set()
      #if (line_sp[4] != "NA" and line_sp[5] != "NA" and line_sp[-1] == "Yes" and int(line_sp[2]) >= minscsup \
      if (int(line_sp[2]) >= minscsup and line_sp[-2] == "Yes"):
        if ("SINE1/7SL" in TE_list or "L1" in TE_list):
          cancerRBP.append(line)
          #if (line_sp[-1] == "No" and line_sp[-2] == "No"):
          if (line_sp[4] != "NA" and line_sp[5] != "NA"):
            cancerBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[5]) + scsrange))
          elif (line_sp[4] != "NA"):
            cancerBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[4]) + scsrange))
          elif (line_sp[5] != "NA"):
            cancerBP[chrom].add((int(line_sp[5]) - scsrange, int(line_sp[5]) + scsrange))

  print("\nReading %s..." % (normalRBPfilename))
  print(str(datetime.datetime.today()))
  with open(normalRBPfilename, 'r') as normalRBPfile:
    header = normalRBPfile.readline()
    file_size = os.stat(normalRBPfilename).st_size
    count = 0
    for line in normalRBPfile:
      count += len(line)
      sys.stdout.write("\r%i%%" % (int((count*100)/file_size)))
      sys.stdout.flush()
      line_sp = line.rstrip('\n').split('\t')
      TE_list = line_sp[3].split(',')
      chrom = line_sp[0]
      if (chrom not in normalBP):
        normalBP[chrom] = set()
      #if (line_sp[4] != "NA" and line_sp[5] != "NA" and line_sp[-1] == "Yes" and int(line_sp[2]) >= minscsup \
      if (int(line_sp[2]) >= minscsup and line_sp[-2] == "Yes"):
        if ("SINE1/7SL" in TE_list or "L1" in TE_list):
          normalRBP.append(line)
          #if (line_sp[-1] == "No" and line_sp[-2] == "No"):
          if (line_sp[4] != "NA" and line_sp[5] != "NA"):
            normalBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[5]) + scsrange))
          elif (line_sp[4] != "NA"):
            normalBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[4]) + scsrange))
          elif (line_sp[5] != "NA"):
            normalBP[chrom].add((int(line_sp[5]) - scsrange, int(line_sp[5]) + scsrange))

  print("\nReading %s..." % (cancerBPfilename))
  print(str(datetime.datetime.today()))
  with open(cancerBPfilename, 'r') as cancerBPfile:
    header = cancerBPfile.readline()
    file_size = os.stat(cancerBPfilename).st_size
    count = 0
    for line in cancerBPfile:
      count += len(line)
      sys.stdout.write("\r%i%%" % (int((count*100)/file_size)))
      sys.stdout.flush()
      line_sp = line.rstrip('\n').split('\t')
      chrom = line_sp[0]
      if (chrom not in cancerBP):
        cancerBP[chrom] = set()
      cancerBP[chrom].add((int(line_sp[2]) - scsrange, int(line_sp[2]) + scsrange))

  print("\nReading %s..." % (normalBPfilename))
  print(str(datetime.datetime.today()))
  with open(normalBPfilename, 'r') as normalBPfile:
    header = normalBPfile.readline()
    file_size = os.stat(normalBPfilename).st_size
    count = 0
    for line in normalBPfile:
      count += len(line)
      sys.stdout.write("\r%i%%" % (int((count*100)/file_size)))
      sys.stdout.flush()
      line_sp = line.rstrip('\n').split('\t')
      chrom = line_sp[0]
      if (chrom not in normalBP):
        normalBP[chrom] = set()
      normalBP[chrom].add((int(line_sp[2]) - scsrange, int(line_sp[2]) + scsrange))

  print("\nReading %s..." % (cancerCRfilename))
  print(str(datetime.datetime.today()))
  with open(cancerCRfilename, 'r') as cancerCRfile:
    header = cancerCRfile.readline()
    file_size = os.stat(cancerCRfilename).st_size
    count = 0
    for line in cancerCRfile:
      count += len(line)
      sys.stdout.write("\r%i%%" % (int((count*100)/file_size)))
      sys.stdout.flush()
      line_sp = line.rstrip('\n').split('\t')
      chrom = line_sp[0]
      if (chrom not in cancerBP):
        cancerBP[chrom] = set()
      left_side = int(line_sp[2])
      right_side = int(line_sp[3])
      #center = (right_side + left_side)/2
      cancerBP[chrom].add((left_side - discsrange, right_side + discsrange))

  print("\nReading %s..." % (normalCRfilename))
  print(str(datetime.datetime.today()))
  with open(normalCRfilename, 'r') as normalCRfile:
    header = normalCRfile.readline()
    file_size = os.stat(normalCRfilename).st_size
    count = 0
    for line in normalCRfile:
      count += len(line)
      sys.stdout.write("\r%i%%" % (int((count*100)/file_size)))
      sys.stdout.flush()
      line_sp = line.rstrip('\n').split('\t')
      chrom = line_sp[0]
      if (chrom not in normalBP):
        normalBP[chrom] = set()
      left_side = int(line_sp[2])
      right_side = int(line_sp[3])
      #center = (right_side + left_side)/2
      normalBP[chrom].add((left_side - discsrange, right_side + discsrange))

  if (os.path.isfile(polymorphfilename)):
    print("\nReading %s..." % (polymorphfilename))
    print(str(datetime.datetime.today()))
    with open(polymorphfilename, 'r') as polymorphfile:
      for line in polymorphfile:
        line_sp = line.rstrip('\n').split('\t')
        chrom = line_sp[0]
        if (chrom not in normalBP):
          normalBP[chrom] = set()
        if (chrom not in cancerBP):
          cancerBP[chrom] = set()
        if (line_sp[4] != "NA" and line_sp[5] != "NA"):
          normalBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[5]) + scsrange))
          cancerBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[5]) + scsrange))
          for i in xrange(int(line_sp[4]) - scsrange, int(line_sp[5]) + scsrange):
            polymorphBP[(chrom, i)] = 1
        elif (line_sp[4] != "NA"):
          normalBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[4]) + scsrange))
          cancerBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[4]) + scsrange))
          for i in xrange(int(line_sp[4]) - scsrange, int(line_sp[4]) + scsrange):
            polymorphBP[(chrom, i)] = 1
        elif (line_sp[5] != "NA"):
          normalBP[chrom].add((int(line_sp[5]) - scsrange, int(line_sp[5]) + scsrange))
          cancerBP[chrom].add((int(line_sp[5]) - scsrange, int(line_sp[5]) + scsrange))
          for i in xrange(int(line_sp[5]) - scsrange, int(line_sp[5]) + scsrange):
            polymorphBP[(chrom, i)] = 1


  overlapfile = open("%s.overlaps.txt" % (patID), 'w+')
  overlapfile.write("Chromosome\tCluster\tSupportingReads\tTEFamily\tLeftBP\tRightBP\tBoth_Align\tTEMatch\n")
  normalonlyfile = open("%s.normalonly.txt" % (patID), 'w+')
  normalonlyfile.write("Chromosome\tCluster\tSupportingReads\tTEFamily\tLeftBP\tRightBP\tBoth_Align\tTEMatch\n")
  canceronlyfile = open("%s.canceronly.txt" % (patID), 'w+')
  canceronlyfile.write("Chromosome\tCluster\tSupportingReads\tTEFamily\tLeftBP\tRightBP\tBoth_Align\tTEMatch\n")

  overlapfiletowrite = {}
  appendpolymorphfile = open(polymorphfilename, 'a+')

  print("\nComparing cancer breakpoints...")
  print(str(datetime.datetime.today()))
  for i in xrange(len(cancerRBP)):
    sys.stdout.write("\r%i%%" % (int((i*100)/len(cancerRBP))))
    sys.stdout.flush()
    line_sp = cancerRBP[i].rstrip('\n').split('\t')
    chrom = line_sp[0]
    found = False
    bp_list = list(normalBP[chrom])
    for j in xrange(len(bp_list)):
      normal_left_bp = bp_list[j][0]
      normal_right_bp = bp_list[j][1]
      if (line_sp[4] != "NA" and int(line_sp[4]) >= normal_left_bp and int(line_sp[4]) <= normal_right_bp):
        overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = cancerRBP[i]
        if (not (chrom, int(line_sp[4])) in polymorphBP):
          appendpolymorphfile.write(cancerRBP[i])
          polymorphBP[(chrom, int(line_sp[4]))] = 1
        found = True
        break
      elif (line_sp[5] != "NA" and int(line_sp[5]) >= normal_left_bp and int(line_sp[5]) <= normal_right_bp):
        overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = cancerRBP[i]
        if (not (chrom, int(line_sp[5])) in polymorphBP):
          appendpolymorphfile.write(cancerRBP[i])
          polymorphBP[(chrom, int(line_sp[5]))] = 1
        found = True
        break
    if (not found):
      canceronlyfile.write(cancerRBP[i])

  print("\nComparing normal breakpoints...")
  print(str(datetime.datetime.today()))
  for i in xrange(len(normalRBP)):
    sys.stdout.write("\r%i%%" % (int((i*100)/len(normalRBP))))
    sys.stdout.flush()
    line_sp = normalRBP[i].rstrip('\n').split('\t')
    chrom = line_sp[0]
    found = False
    bp_list = list(cancerBP[chrom])
    for j in xrange(len(cancerBP[chrom])):
      cancer_left_bp = bp_list[j][0]
      cancer_right_bp = bp_list[j][1]
      if (line_sp[4] != "NA" and int(line_sp[4]) >= cancer_left_bp and int(line_sp[4]) <= cancer_right_bp):
        overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = normalRBP[i]
        if (not (chrom, int(line_sp[4])) in polymorphBP):
          appendpolymorphfile.write(normalRBP[i])
          polymorphBP[(chrom, int(line_sp[4]))] = 1
        found = True
        break
      elif (line_sp[5] != "NA" and int(line_sp[5]) >= cancer_left_bp and int(line_sp[5]) <= cancer_right_bp):
        overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = normalRBP[i]
        if (not (chrom, int(line_sp[5])) in polymorphBP):
          appendpolymorphfile.write(normalRBP[i])
          polymorphBP[(chrom, int(line_sp[5]))] = 1
        found = True
        break
    if (not found):
      normalonlyfile.write(normalRBP[i])

  print("\nWriting overlaps to file...")
  print(str(datetime.datetime.today()))
  for key, value in overlapfiletowrite.iteritems():
    overlapfile.write(value)

  overlapfile.close()
  normalonlyfile.close()
  canceronlyfile.close()
  appendpolymorphfile.close()

  print("Finished comparing results.")
  print(str(datetime.datetime.today()))

if (__name__ == "__main__"):
  main(sys.argv)
