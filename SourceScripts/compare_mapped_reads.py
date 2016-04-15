#!/usr/bin/env python
import sys
import scipy.stats
import itertools
import operator
import re

def TestBit(int_type, offset):
    """
    Checks for a given bit in a int flag
    Returns the value of the bit in the flag
    """
    mask = 1 << offset
    return(int_type & mask)

def Count(L, BPs, cluster):
    """
    Counts the number of occurrences of each element in a list L
    Returns a sorted list of tuples with each element and their count
    from highest to lowest
    """
    return sorted([(x, L.count(x), BPs[(cluster, x)]) for x in set(L)], key=lambda x: x[2], reverse=True)

def GetMatchingBPCount(cigar):
    """
    Gets all the matching BPs (M in the cigar string) and sums their BP count
    """
    return sum([int(x[:-1]) for x in re.findall("\d*M", cigar)])

def GetSubfamilyCounts(filename, clusters): #, subfamilyCounts, subfamilyBPs):
    """
    Reads the generated *.sam to get all of the subfamilies the aligned reads mapped to
    Returns a dictionary with the cluster number as the key and the tuple list of counts
    for each subfamily
    """
    subfamilyCounts = {}
    currentCluster = 0
    clusterFamilies = []
    subfamilyBPs = {}
    with open(filename, 'r') as file:
        for line in file:
            if (line[0:3] != "@SQ" and line[0:3] != "@PG"):
                linesp = line.rstrip('\n').split('\t')
                flag = int(linesp[1])
                if (TestBit(flag, 2) == 0):
                    clusterNum = linesp[0].split(':')[-1]
                    if (currentCluster == 0):
                        currentCluster = clusterNum
                    subfamilies = [linesp[2]]
                    subfamilyBPs[(clusterNum, linesp[2])] = GetMatchingBPCount(linesp[5])
                    if (linesp[-1][:5] == "XA:Z:"):
                        XAMatches = linesp[-1].rstrip('\n').split(';')
                        XAMatches[0] = XAMatches[0][5:]
                        for match in XAMatches[:-1]:
                            matchsp = match.split(',')
                            subfamilies += [matchsp[0]]
                            if ((clusterNum, matchsp[0]) not in subfamilyBPs):
                                subfamilyBPs[(clusterNum, matchsp[0])] = 0
                            # Extra matches add less support
                            subfamilyBPs[(clusterNum, matchsp[0])] += float(GetMatchingBPCount(matchsp[2]) / (len(XAMatches[:-1]) + 1))
                        subfamilyBPs[(clusterNum, linesp[2])] /= float((len(XAMatches[:-1]) + 1))
                    if (currentCluster == clusterNum):
                        clusterFamilies += subfamilies
                    else:
                        subfamilyCounts[currentCluster] = Count(clusterFamilies, subfamilyBPs, currentCluster)
                        clusterFamilies = subfamilies
                        clusters.add(clusterNum)
                        currentCluster = clusterNum
        # For the last cluster
        subfamilyCounts[clusterNum] = Count(clusterFamilies, subfamilyBPs, clusterNum)
        clusterFamilies = subfamilies
        clusters.add(clusterNum)
        currentCluster = clusterNum

    return subfamilyCounts, clusters

def GetTEGroups(filename):
    """
    Gets the list of TE groups for each subfamily and organizes it into a dictionary
    """
    TEGroups = {}
    with open(filename, 'r') as file:
        for line in file:
            linesp = line.rstrip('\n').split('\t')
            group = linesp[0]
            for subfamily in linesp[1:]:
                TEGroups[subfamily] = group
    return TEGroups

def GetBreakpoints(filename, clusters):
    """
    Reads the *.breakpoints file to get the cluster and left/right positions of the breakpoints
    """
    breakpoints = {}
    with open(filename, 'r') as breakpointFile:
        header = breakpointFile.readline()
        for line in breakpointFile:
            line_sp = line.rstrip('\n').split('\t')
            if (line_sp[1] in clusters):
                breakpoints[(line_sp[1], line_sp[-1])] = line.rstrip('\n')

    return breakpoints

def GetClusterRanges(filename, clusters):
    """
    Reads the *.ranges.txt file to get the cluster ranges
    """
    clusterRanges = {}
    with open(filename, 'r') as clusterRangesFile:
        header = clusterRangesFile.readline()
        for line in clusterRangesFile:
            line_sp = line.rstrip('\n').split('\t')
            if (line_sp[1] in clusters):
                clusterRanges[line_sp[1]] = line.rstrip('\n')

    return clusterRanges

def AddSupport(discSubfamilies, softclipSubfamilies, readLength):
    """
    Adds the support values for the same subfamilies between the disc and softclip subfamilies
    """
    counts = {}
    families = set()
    for family, count, bpCount in discSubfamilies:
        counts[family] = bpCount
        families.add(family)

    for family, count, bpCount in softclipSubfamilies:
        if (family not in counts):
            counts[family] = 0
            families.add(family)
        counts[family] += bpCount

    return sorted([(family, counts[family] / readLength) for family in families], key=lambda x: x[1], reverse=True)

def GroupTEs(TEGroups, support):
    """
    Combines all the TE subfamilies into their groups and adds support values
    """
    counts = {}
    families = set()
    for family, count in support:
        if (TEGroups[family] not in counts):
            counts[TEGroups[family]] = 0
            families.add(TEGroups[family])
        counts[TEGroups[family]] += count
    return sorted([(family, counts[family]) for family in families], key=lambda x: x[1], reverse=True)

def WriteCallFile(clusterRanges, breakpoints, discSubfamilies, softclipSubfamilies, clusters, minSupport, readLength, callDiscOnly, TEGroups, filename):
    """
    Iterate over the clusters and write final output to filename. This also calculates the
    added support for softclips and discordant reads.
    """
    outputFile = open(filename, 'w+')
    outputFile.write("Chromosome\tCluster\tSupport\tTEFamily\tLeftBP\t\
            RightBP\tHasBP\tSoftclipAlign\tDiscordantAlign\tTEMatch\tSubFamilyCounts\n")
    for cluster in sorted([int(cluster) for cluster in clusters]):
        cluster = str(cluster)
        if (cluster in softclipSubfamilies):
            softclipAligned = "Yes"
        else:
            softclipAligned = "No"
            softclipSubfamilies[cluster] = []

        if (cluster in discSubfamilies):
            discordantAligned = "Yes"
        else:
            discordantAligned = "No"
            discSubfamilies[cluster] = []

        clusterRange = clusterRanges[cluster].split('\t')
        if ((cluster, "left") in breakpoints):
            leftBP = breakpoints[(cluster, "left")].split('\t')[2]
        else:
            leftBP = "NA"
        if ((cluster, "right") in breakpoints):
            rightBP = breakpoints[(cluster, "right")].split('\t')[2]
        else:
            rightBP = "NA"

        hasBP = "Yes"
        if (leftBP == "NA" and rightBP == "NA"):
            hasBP = "No"
            if (callDiscOnly):
                leftBP = clusterRange[2]
                rightBP = clusterRange[3]
            else:
                continue

        if (discordantAligned == "Yes" and (hasBP == "Yes" or callDiscOnly)):
            subfamilySupports = AddSupport(discSubfamilies[cluster], softclipSubfamilies[cluster], readLength)
            if (TEGroups != None):
                if (discSubfamilies[cluster] != [] and softclipSubfamilies[cluster] != []):
                    TEMatch = "Yes" if TEGroups[discSubfamilies[cluster][0][0]] == TEGroups[softclipSubfamilies[cluster][0][0]] else "No"
                else:
                    TEMatch = "No"
                groupedSupports = GroupTEs(TEGroups, subfamilySupports)
            else:
                if (discSubfamilies[cluster] != [] and softclipSubfamilies[cluster] != []):
                    TEMatch = "Yes" if discSubfamilies[cluster][0][0] == softclipSubfamilies[cluster][0][0] else "No"
                else:
                    TEMatch = "No"
                groupedSupports = subfamilySupports

            if (groupedSupports[0][1] >= minSupport):
                chromosome = clusterRange[0]
                support = "%.2f" % groupedSupports[0][1]
                family = groupedSupports[0][0]
                subfamilyCounts = ','.join(['%s %.2f' % (x,y) for x,y in subfamilySupports])

                outputFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format( \
                        chromosome, cluster, support, family, leftBP, rightBP, hasBP, \
                        softclipAligned, discordantAligned, TEMatch, subfamilyCounts))

    outputFile.close()

def main(args):
    discSamFilename = args[1]
    softclipSamFilename = args[2]
    breakpointFilename = args[3]
    outputFilename = args[4]
    groupedTEsFilename = args[5]
    clusterRangesFilename = args[6]
    minSupport = float(args[7])
    readLength = int(args[8])
    callDiscOnly = int(args[9])

    # Calculate the subfamily support for each call
    discSubfamilies, clusters = GetSubfamilyCounts(discSamFilename, set())
    softclipSubfamilies, clusters = GetSubfamilyCounts(softclipSamFilename, clusters)

    # Grab all the relevant data
    breakpoints = GetBreakpoints(breakpointFilename, clusters)
    clusterRanges = GetClusterRanges(clusterRangesFilename, clusters)
    if (groupedTEsFilename != "none"):
        TEGroups = GetTEGroups(groupedTEsFilename)
    else:
        TEGroups = None

    WriteCallFile(clusterRanges, breakpoints, discSubfamilies, softclipSubfamilies, clusters, \
            minSupport, readLength, callDiscOnly, TEGroups, outputFilename)

if (__name__ == "__main__"):
  main(sys.argv)
