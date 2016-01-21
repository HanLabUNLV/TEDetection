#!/usr/bin/env python
import sys
import os.path
import os
import datetime
import subprocess

def main(args):
  #sys.path.append(args[8])
  import extract_softclips
  patID = args[1]
  cancerext = args[2]
  normalext = args[3]
  scsrange = int(args[4])
  discsrange = int(args[5])
  minscsup = int(args[6])
  polymorphfilename = args[7]
  cancer_bam = args[8]
  norm_bam = args[9]
  TEReffilename = args[10]
  avereaddepth = args[11]
  allow_disconly = int(args[12])

  cancerRBP = []
  cancerBP = {}
  cancerBPclusters = set()
  normalRBP = []
  normalBP = {}
  normalBPclusters = set()
  polymorphBP = {}

  cancerBPfilename = "Results/" + patID + cancerext + ".breakpoints.txt"
  normalBPfilename = "Results/" + patID + normalext + ".breakpoints.txt"
  cancerRBPfilename ="Results/" +  patID + cancerext + ".refined.breakpoints.txt"
  normalRBPfilename ="Results/" +  patID + normalext + ".refined.breakpoints.txt"
  cancerCRfilename = "Results/" + patID + cancerext + ".disc.clusters.ranges.txt"
  normalCRfilename = "Results/" + patID + normalext + ".disc.clusters.ranges.txt"

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
      cluster = line_sp[1]
      cancerBPclusters.add(cluster)
      if (chrom not in cancerBP):
        cancerBP[chrom] = set()
      cancerBP[chrom].add((int(line_sp[2]) - scsrange, int(line_sp[2]) + scsrange))
  sys.stdout.write("\r100%")
  sys.stdout.flush()

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
      cluster = line_sp[1]
      normalBPclusters.add(cluster)
      if (chrom not in normalBP):
        normalBP[chrom] = set()
      normalBP[chrom].add((int(line_sp[2]) - scsrange, int(line_sp[2]) + scsrange))
  sys.stdout.write("\r100%")
  sys.stdout.flush()

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
      cluster = line_sp[1]
      chrom = line_sp[0]
      has_bp = line_sp[6]
      if (chrom not in cancerBP):
        cancerBP[chrom] = set()
      #if (line_sp[4] != "NA" and line_sp[5] != "NA" and line_sp[-1] == "Yes" and int(line_sp[2]) >= minscsup \
      if (has_bp == "Yes" or allow_disconly):
        if ("SINE1/7SL" in TE_list or "L1" in TE_list):
          if (int(line_sp[2]) >= minscsup): #and cluster in cancerBPclusters):
            cancerRBP.append(line)
          #if (line_sp[-1] == "No" and line_sp[-2] == "No"):
          if (line_sp[4] != "NA" and line_sp[5] != "NA"):
            cancerBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[5]) + scsrange))
          elif (line_sp[4] != "NA"):
            cancerBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[4]) + scsrange))
          elif (line_sp[5] != "NA"):
            cancerBP[chrom].add((int(line_sp[5]) - scsrange, int(line_sp[5]) + scsrange))
  sys.stdout.write("\r100%")
  sys.stdout.flush()

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
      cluster = line_sp[1]
      chrom = line_sp[0]
      has_bp = line_sp[6]
      if (chrom not in normalBP):
        normalBP[chrom] = set()
      #if (line_sp[4] != "NA" and line_sp[5] != "NA" and line_sp[-1] == "Yes" and int(line_sp[2]) >= minscsup \
      if (has_bp == "Yes" or allow_disconly):
        if ("SINE1/7SL" in TE_list or "L1" in TE_list):
          if (int(line_sp[2]) >= minscsup): #and cluster in normalBPclusters):
            normalRBP.append(line)
          #if (line_sp[-1] == "No" and line_sp[-2] == "No"):
          if (line_sp[4] != "NA" and line_sp[5] != "NA"):
            normalBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[5]) + scsrange))
          elif (line_sp[4] != "NA"):
            normalBP[chrom].add((int(line_sp[4]) - scsrange, int(line_sp[4]) + scsrange))
          elif (line_sp[5] != "NA"):
            normalBP[chrom].add((int(line_sp[5]) - scsrange, int(line_sp[5]) + scsrange))
  sys.stdout.write("\r100%")
  sys.stdout.flush()

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


  overlapfile = open("Results/%s.overlaps.txt" % (patID), 'w+')
  overlapfile.write("Chromosome\tCluster\tSupportingReads\tTEFamily\tLeftBP\tRightBP\tHas_BP\tSoftclip_Align\tDiscordant_Align\tTEMatch\n")
  normalonlyfile = open("Results/%s.normalonly.txt" % (patID), 'w+')
  normalonlyfile.write("Chromosome\tCluster\tSupportingReads\tTEFamily\tLeftBP\tRightBP\tHas_BP\tSoftclip_Align\tDiscordant_Align\tTEMatch\n")
  canceronlyfile = open("Results/%s.canceronly.txt" % (patID), 'w+')
  canceronlyfile.write("Chromosome\tCluster\tSupportingReads\tTEFamily\tLeftBP\tRightBP\tHas_BP\tSoftclip_Align\tDiscordant_Align\tTEMatch\n")

  overlapfiletowrite = {}
  appendpolymorphfile = open(polymorphfilename, 'a+')

  print("\nComparing cancer breakpoints to normal breakpoints...")
  print(str(datetime.datetime.today()))
  for i in xrange(len(cancerRBP)):
    sys.stdout.write("\r%i%%" % (int((i*100)/len(cancerRBP))))
    sys.stdout.flush()
    line_sp = cancerRBP[i].rstrip('\n').split('\t')
    chrom = line_sp[0]
    found = False
    if (chrom in normalBP):
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
        if (line_sp[4] != "NA" and line_sp[5] != "NA"):
          if ((normal_left_bp >= int(line_sp[4]) and normal_left_bp <= int(line_sp[5]))\
              or (normal_right_bp >= int(line_sp[4]) and normal_right_bp <= int(line_sp[5]))):
            overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = cancerRBP[i]
            if (not (chrom, int(line_sp[4])) in polymorphBP):
              appendpolymorphfile.write(cancerRBP[i])
              polymorphBP[(chrom, int(line_sp[4]))] = 1
            found = True
            break
    if (not found):
      # Check for breakpoints in normal file
      if (line_sp[4] != "NA"):
        rstart = int(line_sp[4])
        if (line_sp[5] != "NA"):
          rend = int(line_sp[5])
        else:
          rend = int(line_sp[4])
      else:
        rstart = int(line_sp[5])
        rend = int(line_sp[5])

      if (rstart > rend):
        temp = rstart
        rstart = rend
        rend = temp

      rstart -= scsrange
      rend += scsrange

      temp_bam_proc = subprocess.Popen("samtools view -hb %s %s:%i-%i > temp.bam && samtools index temp.bam" % (norm_bam, chrom, rstart - int(avereaddepth), rend + int(avereaddepth)), shell=True, stdout=subprocess.PIPE)
      temp_bam_proc.wait()

      temp_ranges_file = open("temp_ranges.txt", 'w+')
      temp_ranges_file.write("Header\n%s\t0\t%i\t%i\t1\n" % (chrom, rstart, rend))
      temp_ranges_file.close()
      extract_softclips.main(["", "temp.bam", "temp_ranges.txt", avereaddepth, TEReffilename, "1", 0, 0, 2, 15, 5, 20, 1])
      temp_bps = []

      with open("Results/temp.breakpoints.txt", 'r') as temp_bp_file:
        temp_bp_file.readline()
        for line in temp_bp_file:
          split_line = line.rstrip('\n').split('\t')
          chr = split_line[0]
          pos = int(split_line[2])
          sup = int(split_line[3])
          temp_bps.append((chr, pos, sup))

      poly = False
      for k in xrange(len(temp_bps)):
        #if (temp_bps[k][0] == chrom and temp_bps[k][2] > 1):
        if (line_sp[4] != "NA" and int(line_sp[4]) >= temp_bps[k][1] - 5 \
            and int(line_sp[4]) <= temp_bps[k][1] + 5):
          poly = True
        elif (line_sp[5] != "NA" and int(line_sp[5]) >= temp_bps[k][1] - 5 \
            and int(line_sp[5]) <= temp_bps[k][1] + 5):
          poly = True
        if (line_sp[4] != "NA" and line_sp[5] != "NA"):
          if (temp_bps[k][1] >= line_sp[4] and temp_bps[k][1] <= line_sp[5]):
            poly = True

      if (not poly):
        canceronlyfile.write(cancerRBP[i])
      else:
        overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = cancerRBP[i]
  sys.stdout.write("\r100%")
  sys.stdout.flush()

  print("\nComparing normal breakpoints to cancer breakpoints...")
  print(str(datetime.datetime.today()))
  for i in xrange(len(normalRBP)):
    sys.stdout.write("\r%i%%" % (int((i*100)/len(normalRBP))))
    sys.stdout.flush()
    line_sp = normalRBP[i].rstrip('\n').split('\t')
    chrom = line_sp[0]
    found = False
    if (chrom in cancerBP):
      bp_list = list(cancerBP[chrom])
      for j in xrange(len(bp_list)):
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
        # Check if breakpoints are between
        if (line_sp[4] != "NA" and line_sp[5] != "NA"):
          if ((cancer_left_bp >= int(line_sp[4]) and cancer_left_bp <= int(line_sp[5]))\
              or (cancer_right_bp >= int(line_sp[4]) and cancer_right_bp <= int(line_sp[5]))):
            overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = normalRBP[i]
            if (not (chrom, int(line_sp[4])) in polymorphBP):
              appendpolymorphfile.write(normalRBP[i])
              polymorphBP[(chrom, int(line_sp[4]))] = 1
            found = True
            break
    if (not found):
      # Check for breakpoints in cancer file
      if (line_sp[4] != "NA"):
        rstart = int(line_sp[4])
        if (line_sp[5] != "NA"):
          rend = int(line_sp[5])
        else:
          rend = int(line_sp[4])
      else:
        rstart = int(line_sp[5])
        rend = int(line_sp[5])

      if (rstart > rend):
        temp = rstart
        rstart = rend
        rend = temp

      rstart -= scsrange
      rend += scsrange

      temp_bam_proc = subprocess.Popen("samtools view -hb %s %s:%i-%i > temp.bam 2> sam.err && samtools index temp.bam" % (cancer_bam, chrom, rstart - int(avereaddepth), rend + int(avereaddepth)), shell=True, stdout=subprocess.PIPE)
      temp_bam_proc.wait()

      temp_ranges_file = open("temp_ranges.txt", 'w+')
      temp_ranges_file.write("Header\n%s\t0\t%i\t%i\t1\n" % (chrom, rstart, rend))
      temp_ranges_file.close()
      extract_softclips.main(["", "temp.bam", "temp_ranges.txt", avereaddepth, TEReffilename, "1", 0, 0, 2, 15, 5, 20, 1])
      temp_bps = []

      with open("Results/temp.breakpoints.txt", 'r') as temp_bp_file:
        temp_bp_file.readline()
        for line in temp_bp_file:
          split_line = line.rstrip('\n').split('\t')
          chr = split_line[0]
          pos = int(split_line[2])
          sup = int(split_line[3])
          temp_bps.append((chr, pos, sup))

      poly = False
      for k in xrange(len(temp_bps)):
        #if (temp_bps[k][0] == chrom and temp_bps[k][2] > 1):
        if (line_sp[4] != "NA" and int(line_sp[4]) >= temp_bps[k][1] - 5 \
            and int(line_sp[4]) <= temp_bps[k][1] + 5):
          poly = True
        elif (line_sp[5] != "NA" and int(line_sp[5]) >= temp_bps[k][1] - 5 \
            and int(line_sp[5]) <= temp_bps[k][1] + 5):
          poly = True
        if (line_sp[4] != "NA" and line_sp[5] != "NA"):
          if (temp_bps[k][1] >= line_sp[4] and temp_bps[k][1] <= line_sp[5]):
            poly = True

      if (not poly):
        normalonlyfile.write(normalRBP[i])
      else:
        overlapfiletowrite[(chrom, line_sp[4], line_sp[5])] = normalRBP[i]
  sys.stdout.write("\r100%")
  sys.stdout.flush()

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
