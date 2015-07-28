#!/usr/bin/env python
import sys
import scipy.stats

def main(args):
  discfilename = args[1]
  softclipfilename = args[2]
  breakpointfilename = args[3]
  outputfilename = args[4]
  groupedTEsfilename = args[5]

  dfmapped = {}
  scmapped = {}
  breakpoints = {}
  TEgroups = {}

  # Checks for a given bit in a flag
  # Used to check sam flag for the bit
  def testBit(int_type, offset):
    mask = 1 << offset
    return(int_type & mask)

  if (groupedTEsfilename != "none"):
    with open(groupedTEsfilename, 'r') as groupedTEsfile:
      for line in groupedTEsfile:
        line_sp = line.rstrip('\n').split('\t')
        for i in xrange(1, len(line_sp)):
          TEgroups[line_sp[i]] = line_sp[0]

  with open(discfilename, 'r') as discfile:
    for line in discfile:
      line_sp = line.split('\t')
      if (line[0:3] != "@SQ" and line[0:3] != "@PG"):
        flag = int(line_sp[1])
        if (testBit(flag, 2) == 0):
          curclus = line_sp[0].split(':')[-1]
          TE = line_sp[2]
          if (curclus in dfmapped):
            dfmapped[curclus].add(TE)
          else:
            dfmapped[curclus] = set()
            dfmapped[curclus].add(TE)
          if (line_sp[-1][:5] == "XA:Z:"):
            XAMatches = line_sp[-1].rstrip('\n').split(';')
            XAMatches[0] = XAMatches[0][5:] # To get rid of the XA:Z:
            for i in xrange(len(XAMatches)-1): # There is an extra ; at the end
              comsplit = XAMatches[i].split(',')
              dfmapped[curclus].add(comsplit[0])

  with open(softclipfilename, 'r') as softclipfile:
    for line in softclipfile:
      line_sp = line.split('\t')
      if (line[0:3] != "@SQ" and line[0:3] != "@PG"):
        flag = int(line_sp[1])
        if (testBit(flag, 2) == 0):
          curclus = line_sp[0].split(':')[1]
          TE = line_sp[2]
          if (curclus in scmapped):
            scmapped[curclus].add(TE)
          else:
            scmapped[curclus] = set()
            scmapped[curclus].add(TE)
          if (line_sp[-1][:5] == "XA:Z:"):
            XAMatches = line_sp[-1].rstrip('\n').split(';')
            XAMatches[0] = XAMatches[0][5:] # To get rid of the XA:Z:
            for i in xrange(len(XAMatches)-1):
              comsplit = XAMatches[i].split(',')
              scmapped[curclus].add(comsplit[0])

  with open(breakpointfilename, 'r') as breakpointfile:
    header = breakpointfile.readline()
    for line in breakpointfile:
      line_sp = line.rstrip('\n').split('\t')
      breakpoints[(line_sp[1], line_sp[-1])] = line.rstrip('\n')

  outputfile = open(outputfilename, 'w')
  header = "Chromosome\tCluster\tSupportingReads\tTEFamily\tLeftBP\tRightBP\tTEMatch\n"
  outputfile.write(header)
  outputfiletowrite = []
  bps = {}

  for key,value in scmapped.iteritems():
    if (key in dfmapped):
      TEmatch = "Yes"
      if (groupedTEsfilename == "none"):
        possibleTEs = ','.join(value.intersection(dfmapped[key]))
        if (possibleTEs == ''):
          possibleTEs = ','.join(value.union(dfmapped[key]))
          TEmatch = "No"
      else:
        sclist = list(value)
        dflist = list(dfmapped[key])
        scgroups = set()
        dfgroups = set()
        for i in xrange(len(sclist)):
          scgroups.add(TEgroups[sclist[i]])
        for i in xrange(len(dflist)):
          dfgroups.add(TEgroups[dflist[i]])

        possibleTEs = ','.join(scgroups.intersection(dfgroups))
        if (possibleTEs == ''):
          possibleTEs = ','.join(value.union(dfmapped[key]))
          TEmatch = "No"

      if ((key, "left") in breakpoints):
        left_bp_split = breakpoints[(key, "left")].split('\t')
        left_bp = left_bp_split[2]
        chrom = left_bp_split[0]
        clus = left_bp_split[1]
        support = left_bp_split[3]
      else:
        left_bp = "NA"

      if ((key, "right") in breakpoints):
        right_bp_split = breakpoints[(key, "right")].split('\t')
        right_bp = right_bp_split[2]
        chrom = right_bp_split[0]
        clus = right_bp_split[1]
        support = right_bp_split[3]
      else:
        right_bp = "NA"

      if (not (left_bp, right_bp, chrom) in bps):
        outputfiletowrite.append(chrom + '\t' + clus + '\t' + str(support) + '\t' + possibleTEs + '\t' + left_bp + '\t' + right_bp + '\t' + TEmatch + '\n')
        bps[(left_bp, right_bp, chrom)] = 1
    else: # If no discordant reads mapped to anything
      if ((key, "left") in breakpoints):
        left_bp_split = breakpoints[(key, "left")].split('\t')
        left_bp = left_bp_split[2]
        chrom = left_bp_split[0]
        clus = left_bp_split[1]
        support = left_bp_split[3]
      else:
        left_bp = "NA"

      if ((key, "right") in breakpoints):
        right_bp_split = breakpoints[(key, "right")].split('\t')
        right_bp = right_bp_split[2]
        chrom = right_bp_split[0]
        clus = right_bp_split[1]
        support = right_bp_split[3]
      else:
        right_bp = "NA"

      if (groupedTEsfilename == "none"):
        possibleTEs = ','.join(value)
      else:
        sclist = list(value)
        scgroups = set()
        for i in xrange(len(sclist)):
          scgroups.add(TEgroups[sclist[i]])
        possibleTEs = ','.join(scgroups)
      TEmatch = "Yes"
      if (not (left_bp, right_bp, chrom) in bps):
          outputfiletowrite.append(chrom + '\t' + clus + '\t' + str(support) + '\t' + possibleTEs + '\t' + left_bp + '\t' + right_bp + '\t' + TEmatch + '\n')
          bps[(left_bp, right_bp, chrom)] = 1

  for key, value in dfmapped.iteritems(): # call reads with only discordant pairs aligning to TE
    if (key not in scmapped):
      if ((key, "left") in breakpoints):
        left_bp_split = breakpoints[(key, "left")].split('\t')
        left_bp = left_bp_split[2]
        chrom = left_bp_split[0]
        clus = left_bp_split[1]
        support = left_bp_split[3]
      else:
        left_bp = "NA"

      if ((key, "right") in breakpoints):
        right_bp_split = breakpoints[(key, "right")].split('\t')
        right_bp = right_bp_split[2]
        chrom = right_bp_split[0]
        clus = right_bp_split[1]
        support = right_bp_split[3]
      else:
        right_bp = "NA"

      if (groupedTEsfilename == "none"):
        possibleTEs = ','.join(value)
      else:
        dflist = list(value)
        dfgroups = set()
        for i in xrange(len(dflist)):
          dfgroups.add(TEgroups[dflist[i]])
        possibleTEs = ','.join(dfgroups)
      TEmatch = "Yes"
      if (not (left_bp, right_bp, chrom) in bps):
          outputfiletowrite.append(chrom + '\t' + clus + '\t' + str(support) + '\t' + possibleTEs + '\t' + left_bp + '\t' + right_bp + '\t' + TEmatch + '\n')
          bps[(left_bp, right_bp, chrom)] = 1

  outputfiletowrite = sorted(list(set(outputfiletowrite)), key=lambda x: int(x.split('\t')[1]))

  i=0
  while (i + 1 < len(outputfiletowrite)):
    first_sp = outputfiletowrite[i].rstrip('\n').split('\t')
    second_sp = outputfiletowrite[i+1].rstrip('\n').split('\t')
    if (first_sp[0] == second_sp[0]):
      if (first_sp[4] == second_sp[4] and (first_sp[4] != "NA")) \
          or (first_sp[5] == second_sp[5] and (first_sp[5] != "NA")):
        if (first_sp[4] == "NA" or first_sp[5] == "NA"):
          outputfiletowrite.remove(outputfiletowrite[i])
        elif (second_sp[4] == "NA" or second_sp[5] == "NA"):
          outputfiletowrite.remove(outputfiletowrite[i+1])
    i += 1

  for line in outputfiletowrite:
    outputfile.write(line)

  outputfile.close()

if (__name__ == "__main__"):
  main(sys.argv)
