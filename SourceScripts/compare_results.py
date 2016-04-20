#!/usr/bin/env python
import sys
import argparse
import pysam
import os.path

def ParseArgs():
    """
    Sets up the argument parser and returns the parsed command line arguments
    """
    parser = argparse.ArgumentParser(description="Compare the breakpoints between the cancer and normal patient files and identify polymorphisms.", \
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-sr', dest='searchRange', \
            type=int, default=100, help='Search range to extend breakpoints searching for polymorphisms')
    parser.add_argument('-id', dest='patientID', \
            help='ID of the patient in the .bam file without extensions', required=True)
    parser.add_argument('-ce', dest='cancerExt', \
            help='File name extension of the cancer .bam file', required=True)
    parser.add_argument('-ne', dest='normalExt', \
            help='File name extension of the normal .bam file', required=True)
    parser.add_argument('-pf', dest='polymorphBasename', required=True, \
            help='File path and basename (without the chromosome #) for the polymorphisms files')
    parser.add_argument('-cb', dest='cancerBam', \
            help='File path to cancer .bam file', required=True)
    parser.add_argument('-nb', dest='normalBam', \
            help='File path to normal .bam file', required=True)
    parser.add_argument('-rf', dest="resultsFolder", \
            help="File path to Results/ folder from TEDetection", default="Results/")

    return parser.parse_args()

def GetPolymorphisms(polymorphFilenames):
    """
    Reads all the polymorphisms*.txt files to get all the polymorph breakpoints.
    Returns a dictionary of polymorphisms for each chromosome with the list of breakpoints
    """
    polymorphisms = {}
    for chromosome, filename in polymorphFilenames:
        if (os.path.isfile(filename)):
            polymorphisms[chromosome] = []
            with open(filename, 'r') as polyFile:
                for line in polyFile:
                    linesp = line.rstrip('\n').split('\t')
                    leftBP = linesp[4]
                    rightBP = linesp[5]

                    if (leftBP == "NA"):
                        leftBP = int(rightBP)
                    if (rightBP == "NA"):
                        rightBP = int(leftBP)

                    leftBP = int(leftBP)
                    rightBP = int(rightBP)

                    if (leftBP > rightBP):
                        # Swap the variables
                        leftBP, rightBP = rightBP, leftBP

                    polymorphisms[chromosome].append((leftBP, rightBP))

    return polymorphisms

def GetBreakpoints(filename):
    """
    Reads a *.refined.breakpoints.txt file and returns a list of tuples with each breakpoint
    """
    breakpoints = []
    fileLines = {}
    with open(filename, 'r') as file:
        header = file.readline()
        for line in file:
            linesp = line.rstrip('\n').split('\t')
            chromosome = linesp[0]
            cluster = linesp[1]
            leftBP = linesp[4]
            rightBP = linesp[5]
            fileLines[cluster] = line

            if (leftBP == "NA"):
                leftBP = int(rightBP)
            if (rightBP == "NA"):
                rightBP = int(leftBP)

            leftBP = int(leftBP)
            rightBP = int(rightBP)

            if (leftBP > rightBP):
                # Swap the variables
                leftBP, rightBP = rightBP, leftBP

            breakpoints.append((chromosome, leftBP, rightBP, cluster))

    return fileLines, breakpoints, header

def GetDiscRanges(filename, searchRange, mappedClusters):
    """
    Reads the *.disc.clusters.ranges.txt files and returns a list of tuples with the discordant
    cluster ranges. This also adds the search range to each range.
    """
    discRanges = {}
    with open(filename, 'r') as file:
        header = file.readline()
        for line in file:
            linesp = line.rstrip('\n').split('\t')
            chromosome = linesp[0]
            if (chromosome not in discRanges):
                discRanges[chromosome] = []
            cluster = linesp[1]
            if cluster in mappedClusters:
                leftBP = int(linesp[2]) - searchRange
                rightBP = int(linesp[3]) + searchRange
                discRanges[chromosome].append((leftBP, rightBP, cluster))

    return discRanges

def CalculateQuality(qualitySeq, softclipLength, minSoftclipLength, phredOffset, side):
    """
    Calculates the average quality of the softclip sequence up to the minimum softclip length.
    The side indicates which side the softclip is on
    """
    if (side == 0):
        return (sum(map(ord, qualitySeq[softclipLength - \
                minSoftclipLength:softclipLength])) - \
                phredOffset*minSoftclipLength)/minSoftclipLength
    else:
        return (sum(map(ord, qualitySeq[-softclipLength: \
                -softclipLength + minSoftclipLength])) - \
                phredOffset*minSoftclipLength)/minSoftclipLength

def CheckSoftclips(bamFile, chromosome, position, side):
    """
    Extracts the softclips around a breakpoint and returns true if at least 2
    are found around the position. Returns false otherwise
    """
    minQuality = 5
    maxMismatches = 2
    softclipID = 4
    minSoftclipLength = 5
    minPhredQuality = 20
    minNumberOfSoftclips = 1
    phredOffset = 33

    numberOfSoftclips = 0
    if (side == 0):
        # Left side has softclip position closer
        bamReadIter = bamFile.fetch(chromosome, position - 5, position + 5).__iter__()
    else:
        # Right side can have read position far from softclip
        bamReadIter = bamFile.fetch(chromosome, position - 100, position + 100).__iter__()
    for read in bamReadIter:
        if (read.cigar and read.cigar[side][0] == softclipID and \
                not read.is_duplicate and read.mapq >= minQuality):
            mismatches = filter(lambda x: x[0] == "NM", read.tags)
            if (mismatches or mismatches[0][1] <= maxMismatches):
                softclipLength = read.cigar[side][1]
                if (side == 0):
                    softclipPosition = read.pos
                else:
                    softclipPosition = read.pos + read.rlen - softclipLength + 1
                    if (read.cigar[0][0] == softclipID):
                        softclipPosition -= read.cigar[0][1]

                if (softclipLength >= minSoftclipLength):
                    averageQuality = CalculateQuality(read.qual, softclipLength, \
                            minSoftclipLength, phredOffset, side)
                    if (averageQuality >= minPhredQuality and \
                            position - 5 <= softclipPosition <= position + 5):
                        numberOfSoftclips += 1
                        if (numberOfSoftclips >= minNumberOfSoftclips):
                            return True

    return False

def CompareBreakpoints(cancerLines, cancerBPs, cancerDiscRanges, normalLines, normalBPs, normalDiscRanges, polymorphisms, searchRange, polymorphFilenames, cancerBamFilename, normalBamFilename, resultsFolder, patientID, header):
    """
    Compares the cancer breakpoints to the normal breakpoints and polymorphisms.
    Once the cancerBPs and normalBPs have been checked here, they don't need to be
    checked again for the normal.
    """
    cancerBamFile = pysam.Samfile(cancerBamFilename, 'rb')
    normalBamFile = pysam.Samfile(normalBamFilename, 'rb')
    polymorphFiles = {}
    for chrom, filename in polymorphFilenames:
        polymorphFiles[chrom] = open(filename, 'a+')

    cancerOnlyFile = open(resultsFolder + patientID + ".canceronly.txt", 'w+')
    normalOnlyFile = open(resultsFolder + patientID + ".normalonly.txt", 'w+')
    cancerOnlyFile.write(header)
    normalOnlyFile.write(header)

    for cancerChrom, cancerLeftBP, cancerRightBP, cancerCluster in cancerBPs:
        found = False
        # Check against the normal breakpoints
        for normalChrom, normalLeftBP, normalRightBP, normalCluster in normalBPs:
            adjustedLeftBP = normalLeftBP - searchRange
            adjustedRightBP = normalRightBP + searchRange
            if (cancerChrom == normalChrom) and (adjustedLeftBP <= cancerLeftBP <= adjustedRightBP or \
                    adjustedLeftBP <= cancerRightBP <= adjustedRightBP):
                found = True
                # Remove the subfamilies from polymorph file line
                polymorphLine = '\t'.join(cancerLines[cancerCluster].split('\t')[:-1]) + '\n'
                polymorphFiles[cancerChrom].write(polymorphLine)
                cancerLines.pop(cancerCluster, 0)
                normalLines.pop(normalCluster, 0)
                break

        if found:
            continue

        # Check against the normal discordant ranges
        for discLeftBP, discRightBP, discCluster in normalDiscRanges[cancerChrom]:
            if (discLeftBP <= cancerLeftBP <= discRightBP or \
                    discLeftBP <= cancerRightBP <= discRightBP):
                found = True
                # Remove the subfamilies from polymorph file line
                polymorphLine = '\t'.join(cancerLines[cancerCluster].split('\t')[:-1]) + '\n'
                polymorphFiles[cancerChrom].write(polymorphLine)
                cancerLines.pop(cancerCluster, 0)
                normalLines.pop(discCluster, 0)
                break

        if found:
            continue

        # Check against the polymorphisms
        if (cancerChrom in polymorphisms):
            for polyLeftBP, polyRightBP in polymorphisms[cancerChrom]:
                if (polyLeftBP <= cancerLeftBP <= polyRightBP or \
                        polyLeftBP <= cancerRightBP <= polyRightBP):
                    found = True
                    cancerLines.pop(cancerCluster, 0)
                    break

        if found:
            continue

        # Final check for softclips in normal before calling insertion
        splitLine = cancerLines[cancerCluster].split('\t')
        if ((splitLine[4] != "NA" and CheckSoftclips(normalBamFile, cancerChrom, int(splitLine[4]), 0)) \
                or (splitLine[5] != "NA" and CheckSoftclips(normalBamFile, cancerChrom, int(splitLine[5]), -1))):
            found = True
            # Remove the subfamilies from polymorph file line
            polymorphLine = '\t'.join(cancerLines[cancerCluster].split('\t')[:-1]) + '\n'
            polymorphFiles[cancerChrom].write(polymorphLine)
            cancerLines.pop(cancerCluster, 0)
            continue

        # Passed all checks, add to calls
        cancerOnlyFile.write(cancerLines[cancerCluster])

    for normalChrom, normalLeftBP, normalRightBP, normalCluster in normalBPs:
        found = False
        # Cluster could have been removed while checking cancer breakpoints
        if (normalCluster in normalLines):
            # Check against the cancer discordant ranges
            for discLeftBP, discRightBP, discCluster in cancerDiscRanges[normalChrom]:
                if (discLeftBP <= normalLeftBP <= discRightBP or \
                        discLeftBP <= normalRightBP <= discRightBP):
                    found = True
                    # Remove the subfamilies from polymorph file line
                    polymorphLine = '\t'.join(normalLines[normalCluster].split('\t')[:-1]) + '\n'
                    polymorphFiles[normalChrom].write(polymorphLine)
                    cancerLines.pop(normalCluster, 0)
                    normalLines.pop(discCluster, 0)
                    break

            if found:
                continue

            # Check against the polymorphisms
            if (normalChrom in polymorphisms):
                for polyLeftBP, polyRightBP in polymorphisms[normalChrom]:
                    if (polyLeftBP <= normalLeftBP <= polyRightBP or \
                            polyLeftBP <= normalRightBP <= polyRightBP):
                        found = True
                        normalLines.pop(normalCluster, 0)
                        break

            if found:
                continue

            # Final check for softclips in cancer before calling insertion
            splitLine = normalLines[normalCluster].split('\t')
            if ((splitLine[4] != "NA" and CheckSoftclips(cancerBamFile, normalChrom, int(splitLine[4]), 0)) \
                    or (splitLine[5] != "NA" and CheckSoftclips(cancerBamFile, normalChrom, int(splitLine[5]), -1))):
                found = True
                # Remove the subfamilies from polymorph file line
                polymorphLine = '\t'.join(normalLines[normalCluster].split('\t')[:-1]) + '\n'
                polymorphFiles[normalChrom].write(polymorphLine)
                normalLines.pop(discCluster, 0)
                continue

            # Passed all checks, add to calls
            normalOnlyFile.write(normalLines[normalCluster])

    cancerOnlyFile.close()
    normalOnlyFile.close()
    cancerBamFile.close()
    normalBamFile.close()

def main():
    args = ParseArgs()

    cancerRefinedBPFilename = args.resultsFolder + args.patientID + args.cancerExt + ".refined.breakpoints.txt"
    normalRefinedBPFilename = args.resultsFolder + args.patientID + args.normalExt + ".refined.breakpoints.txt"
    cancerMappedClustersFilename = args.resultsFolder + args.patientID + args.cancerExt + ".refined.breakpoints.mappedreads.txt"
    normalMappedClustersFilename = args.resultsFolder + args.patientID + args.normalExt + ".refined.breakpoints.mappedreads.txt"
    cancerClustersFilename = args.resultsFolder + args.patientID + args.cancerExt + ".disc.clusters.ranges.txt"
    normalClustersFilename = args.resultsFolder + args.patientID + args.normalExt + ".disc.clusters.ranges.txt"

    chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']
    polymorphFilenames = [(chrom, args.polymorphBasename + chrom + ".txt") for chrom in chromosomes]

    polymorphisms = GetPolymorphisms(polymorphFilenames)

    cancerLines, cancerBPs, header = GetBreakpoints(cancerRefinedBPFilename)
    normalLines, normalBPs, header = GetBreakpoints(normalRefinedBPFilename)

    cancerMappedClusters = set()
    normalMappedClusters = set()
    with open(cancerMappedClustersFilename, 'r') as cancerMappedFile:
        for line in cancerMappedFile:
            cancerMappedClusters.add(line.rstrip('\n'))
    with open(normalMappedClustersFilename, 'r') as normalMappedFile:
        for line in normalMappedFile:
            normalMappedClusters.add(line.rstrip('\n'))

    cancerDiscRanges = GetDiscRanges(cancerClustersFilename, args.searchRange, cancerMappedClusters)
    normalDiscRanges = GetDiscRanges(normalClustersFilename, args.searchRange, normalMappedClusters)

    CompareBreakpoints(cancerLines, cancerBPs, cancerDiscRanges, normalLines, normalBPs, normalDiscRanges, polymorphisms, args.searchRange, polymorphFilenames, args.cancerBam, args.normalBam, args.resultsFolder, args.patientID, header)

if __name__ == '__main__':
    main()
