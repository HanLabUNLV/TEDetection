#!/usr/bin/env python
import sys
import scipy.stats
import itertools
import operator

def main(args):
  discfilename = args[1]
  softclipfilename = args[2]
  breakpointfilename = args[3]
  outputfilename = args[4]
  groupedTEsfilename = args[5]
  clusterrangesfilename = args[6]
  mappeddiscfilename = args[7]
  mappedscfilename = args[8]
  minsupport = int(args[9])
  calldisconly = int(args[10])

  dfmapped = {}
  scmapped = {}
  breakpoints = {}
  TEgroups = {}
  cluster_ranges = {}

  sc_align = "Yes"
  dc_align = "Yes"

  # Checks for a given bit in a flag
  # Used to check sam flag for the bit
  def testBit(int_type, offset):
    mask = 1 << offset
    return(int_type & mask)

  """
  # Count the most common TE
  def most_common(L):
    # get sorted list with iterable pairs
    SL = sorted((x, i) for i, x in enumerate(L))
    # separate into groups
    groups = itertools.groupby(SL, key=operator.itemgetter(0))
    counts = {}
    # auxillary function to get count of an item
    def _auxfun(g):
      item, iterable = g
      count = 0
      min_index = len(L)
      for _, where in iterable:
        count += 1
        min_index = min(min_index, where)
      counts[item] = count
      return count, -min_index
    max_item = max(groups, key=_auxfun)[0]
    count = counts[max_item]
    return max_item, count
  """

  def count(L):
    return sorted([(x, L.count(x)) for x in set(L)], key=lambda x: x[1], reverse=True)

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
            dfmapped[curclus].append(TE)
          else:
            dfmapped[curclus] = []
            dfmapped[curclus].append(TE)
          if (line_sp[-1][:5] == "XA:Z:"):
            XAMatches = line_sp[-1].rstrip('\n').split(';')
            XAMatches[0] = XAMatches[0][5:] # To get rid of the XA:Z:
            for i in xrange(len(XAMatches)-1): # There is an extra ; at the end
              comsplit = XAMatches[i].split(',')
              dfmapped[curclus].append(comsplit[0])

  with open(softclipfilename, 'r') as softclipfile:
    for line in softclipfile:
      line_sp = line.split('\t')
      if (line[0:3] != "@SQ" and line[0:3] != "@PG"):
        flag = int(line_sp[1])
        if (testBit(flag, 2) == 0):
          curclus = line_sp[0].split(':')[1]
          TE = line_sp[2]
          if (curclus in scmapped):
            scmapped[curclus].append(TE)
          else:
            scmapped[curclus] = []
            scmapped[curclus].append(TE)
          if (line_sp[-1][:5] == "XA:Z:"):
            XAMatches = line_sp[-1].rstrip('\n').split(';')
            XAMatches[0] = XAMatches[0][5:] # To get rid of the XA:Z:
            for i in xrange(len(XAMatches)-1):
              comsplit = XAMatches[i].split(',')
              scmapped[curclus].append(comsplit[0])

  with open(breakpointfilename, 'r') as breakpointfile:
    header = breakpointfile.readline()
    for line in breakpointfile:
      line_sp = line.rstrip('\n').split('\t')
      breakpoints[(line_sp[1], line_sp[-1])] = line.rstrip('\n')

  with open(clusterrangesfilename, 'r') as clusterrangesfile:
    header = clusterrangesfile.readline()
    for line in clusterrangesfile:
      line_sp = line.rstrip('\n').split('\t')
      cluster_ranges[line_sp[1]] = line.rstrip('\n')

  outputfile = open(outputfilename, 'w')
  header = "Chromosome\tCluster\tSupportingReads\tTEFamily\tLeftBP\tRightBP\tHas_BP\tSoftclip_Align\tDiscordant_Align\tTEMatch\tSubFamily\n"
  outputfile.write(header)
  outputfiletowrite = []
  bps = {}
  sc_TE = {}
  sc_count = {}
  dc_TE = {}
  dc_count = {}

  for key,value in scmapped.iteritems():
    if (key in dfmapped):
      has_bp = "No"
      grouped_sc_TEs = value
      grouped_dc_TEs = dfmapped[key]
      sc_align = "Yes"
      dc_align = "Yes"
      TEmatch = "Yes"
      #subFam, subCount = most_common(dfmapped[key])
      subFam = ','.join(['{0} {1}'.format(x,y) for x,y in count(dfmapped[key])])
      if (groupedTEsfilename != "none"):
        #for i in xrange(len(value)):
          #grouped_sc_TEs[i] = TEgroups[value[i]]
        grouped_sc_TEs = [TEgroups[i] for i in value]
        #for i in xrange(len(dfmapped[key])):
          #grouped_dc_TEs[i] = TEgroups[dfmapped[key][i]]
        grouped_dc_TEs = [TEgroups[i] for i in dfmapped[key]]
      sc_common = count(grouped_sc_TEs)
      dc_common = count(grouped_dc_TEs)
      sc_TE[key], sc_count[key] = sc_common[0]
      dc_TE[key], dc_count[key] = dc_common[0]
      if (sc_TE[key] != dc_TE[key]):
        TEmatch = "No"
      support = dc_count[key]
      ins_TE = dc_TE[key]

      if ((key, "left") in breakpoints):
        left_bp_split = breakpoints[(key, "left")].split('\t')
        left_bp = left_bp_split[2]
        chrom = left_bp_split[0]
        support += int(left_bp_split[3])
        has_bp = "Yes"
      else:
        left_bp = "NA"

      if ((key, "right") in breakpoints):
        right_bp_split = breakpoints[(key, "right")].split('\t')
        right_bp = right_bp_split[2]
        chrom = right_bp_split[0]
        support += int(right_bp_split[3])
        has_bp = "Yes"
      else:
        right_bp = "NA"

      clus = key

      if (not (left_bp, right_bp, chrom) in bps and support >= minsupport):
        outputfiletowrite.append(chrom + '\t' + clus + '\t' + str(support) + '\t' + ins_TE + '\t' + left_bp + '\t' + right_bp + '\t' + has_bp + '\t' + sc_align + '\t' + dc_align + '\t' + TEmatch + '\t' + subFam + '\n')
        bps[(left_bp, right_bp, chrom)] = 1

    else: # If no discordant reads mapped to anything
      has_bp = "No"
      dc_align = "No"
      sc_align = "Yes"
      grouped_sc_TEs = value
      #subFam, subCount = most_common(value)
      subFam = ','.join(['{0} {1}'.format(x,y) for x,y in count(value)])
      if (groupedTEsfilename != "none"):
        #for i in xrange(len(value)):
          #grouped_sc_TEs[i] = TEgroups[value[i]]
        grouped_sc_TEs = [TEgroups[i] for i in value]
      sc_common = count(grouped_sc_TEs)
      sc_TE[key], sc_count[key] = sc_common[0]
      ins_TE = sc_TE[key]
      support = 0

      if ((key, "left") in breakpoints):
        left_bp_split = breakpoints[(key, "left")].split('\t')
        left_bp = left_bp_split[2]
        chrom = left_bp_split[0]
        support += int(left_bp_split[3])
        has_bp = "Yes"
      else:
        left_bp = "NA"

      if ((key, "right") in breakpoints):
        right_bp_split = breakpoints[(key, "right")].split('\t')
        right_bp = right_bp_split[2]
        chrom = right_bp_split[0]
        support += int(right_bp_split[3])
        has_bp = "Yes"
      else:
        right_bp = "NA"

      clus = key
      TEmatch = "No"
      if (not (left_bp, right_bp, chrom) in bps and support >= minsupport):
          outputfiletowrite.append(chrom + '\t' + clus + '\t' + str(support) + '\t' + ins_TE + '\t' + left_bp + '\t' + right_bp + '\t' + has_bp + '\t' + sc_align + '\t' + dc_align + '\t' + TEmatch + '\t' + subFam + '\n')
          bps[(left_bp, right_bp, chrom)] = 1

  if (calldisconly):
    for key, value in dfmapped.iteritems(): # call reads with only discordant pairs aligning to TE
      if (key not in scmapped):
        has_bp = "No"
        sc_align = "No"
        dc_align = "Yes"
        #subFam, subCount = most_common(value)
        subFam = ','.join(['{0} {1}'.format(x,y) for x,y in count(value)])
        grouped_dc_TEs = value
        if (groupedTEsfilename != "none"):
          #for i in xrange(len(value)):
            #grouped_dc_TEs[i] = TEgroups[value[i]]
          grouped_dc_TEs = [TEgroups[i] for i in value]
        dc_common = count(grouped_dc_TEs)
        int_TE, support = dc_common[0]

        cluster_range = cluster_ranges[key].split('\t')
        chrom = cluster_range[0]
        clus = cluster_range[1]
        left_bp = cluster_range[2]
        right_bp = cluster_range[3]

        if ((key, "left") in breakpoints):
          left_bp_split = breakpoints[(key, "left")].split('\t')
          #left_bp = left_bp_split[2]
          support += int(left_bp_split[3])
          has_bp = "Yes"

        if ((key, "right") in breakpoints):
          right_bp_split = breakpoints[(key, "right")].split('\t')
          #right_bp = right_bp_split[2]
          support += int(right_bp_split[3])
          has_bp = "Yes"

        TEmatch = "No"
        if (not (left_bp, right_bp, chrom) in bps and support >= minsupport):
            outputfiletowrite.append(chrom + '\t' + clus + '\t' + str(support) + '\t' + ins_TE + '\t' + left_bp + '\t' + right_bp + '\t' + has_bp + '\t' + sc_align + '\t' + dc_align + '\t' + TEmatch + '\t' + subFam + '\n')
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

  mapped_clust_file = open(mappeddiscfilename, 'w+')
  mapped_clust_file.write("Cluster\tTE\tCount\n")
  for clusnum, cluster_TE in dc_TE.iteritems():
    cluster_count = dc_count[clusnum]
    mapped_clust_file.write("%s\t%s\t%i\n" %(clusnum, cluster_TE, cluster_count))

  mapped_clust_file.close()

  mapped_sc_file = open(mappedscfilename, 'w+')
  mapped_sc_file.write("Cluster\tTE\tCount\n")
  for clusnum, cluster_TE in sc_TE.iteritems():
    cluster_count = sc_count[clusnum]
    mapped_sc_file.write("%s\t%s\t%i\n" %(clusnum, cluster_TE, cluster_count))

  mapped_sc_file.close()

if (__name__ == "__main__"):
  main(sys.argv)
