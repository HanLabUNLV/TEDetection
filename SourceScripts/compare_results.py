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
    parser.add_argument('-ce', dest='cancerExt', default = ".cancer", \
            help='File name extension of the cancer .bam file')
    parser.add_argument('-ne', dest='normalExt', default = ".normal", \
            help='File name extension of the normal .bam file')
    parser.add_argument('-pf', dest='polymorphBasename', required=True, \
            help='File path and basename (without the chromosome #) for the polymorphisms files')
    parser.add_argument('-rf', dest="resultsFolder", \
            help="File path to Results/ folder from TEDetection", default="Results/")
    parser.add_argument('-ex', dest="callFileExt", \
            help="Extension for the call file", default = ".refined.breakpoints.txt")
    parser.add_argument('-mr', dest="maxRange", \
            help="Maximum range between breakpoints for an insertion", type=int, default = 4000)

    return parser.parse_args()

def GetPolymorphisms(polymorphFilenames, MAXRANGE):
    """
    Reads all the polymorphisms*.txt files to get all the polymorph breakpoints.
    Returns a dictionary of polymorphisms for each chromosome with the list of breakpoints
    """
    polymorphisms = {}
    for chromosome, filename in polymorphFilenames:
        if (os.path.isfile(filename)):
            polymorphisms[chromosome] = []
            with open(filename, 'r') as polyFile:
                lineNumber = 0
                for line in polyFile:
                    lineNumber += 1
                    linesp = line.rstrip('\n').split('\t')
                    TEFamily = linesp[3]
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

                    if (rightBP - leftBP > MAXRANGE):
                        sys.stderr.write("Insertion BPs greater than max range at line: " + str(lineNumber) + " in polymorphisms file at chromosome: " + str(chromosome) + "\n")
                        continue

                    polymorphisms[chromosome].append((leftBP, rightBP, TEFamily))

    return polymorphisms

def GetBreakpoints(filename, MAXRANGE):
    """
    Reads a *.refined.breakpoints.txt file and returns a list of tuples with each breakpoint
    """
    breakpoints = []
    fileLines = {}
    with open(filename, 'r') as file:
        header = file.readline()
        lineNumber = 1
        for line in file:
            lineNumber += 1
            linesp = line.rstrip('\n').split('\t')
            chromosome = linesp[0]
            cluster = linesp[1]
            TEFamily = linesp[3]
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

            if (rightBP - leftBP > MAXRANGE):
                sys.stderr.write("Insertion BPs greater than max range at line: " + str(lineNumber) + " in " + filename + "\n")
                continue

            breakpoints.append((chromosome, leftBP, rightBP, cluster, TEFamily))

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
                TEFamily = mappedClusters[cluster]
                leftBP = int(linesp[2]) - searchRange
                rightBP = int(linesp[3]) + searchRange
                discRanges[chromosome].append((leftBP, rightBP, cluster, TEFamily))

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

def CompareBreakpoints(cancerLines, cancerBPs, cancerDiscRanges, normalLines, normalBPs, normalDiscRanges, polymorphisms, searchRange, polymorphFilenames, resultsFolder, patientID, header):
    """
    Compares the cancer breakpoints to the normal breakpoints and polymorphisms.
    Once the cancerBPs and normalBPs have been checked here, they don't need to be
    checked again for the normal.
    """
    polymorphFiles = {}
    for chrom, filename in polymorphFilenames:
        polymorphFiles[chrom] = open(filename, 'a+')

    cancerOnlyFile = open(resultsFolder + patientID + ".canceronly.txt", 'w+')
    normalOnlyFile = open(resultsFolder + patientID + ".normalonly.txt", 'w+')
    cancerOnlyFile.write(header)
    normalOnlyFile.write(header)

    for cancerChrom, cancerLeftBP, cancerRightBP, cancerCluster, cancerTEFamily in cancerBPs:
        found = False
        # Check against the normal breakpoints
        for normalChrom, normalLeftBP, normalRightBP, normalCluster, normalTEFamily in normalBPs:
            adjustedLeftBP = normalLeftBP - searchRange
            adjustedRightBP = normalRightBP + searchRange
            if (cancerChrom == normalChrom) and \
                    (cancerTEFamily == normalTEFamily) and \
                    (adjustedLeftBP <= cancerLeftBP <= adjustedRightBP or \
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
        for discLeftBP, discRightBP, discCluster, discTEFamily in normalDiscRanges[cancerChrom]:
            if (cancerTEFamily == discTEFamily) and \
                    (discLeftBP <= cancerLeftBP <= discRightBP or \
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
            for polyLeftBP, polyRightBP, polyTEFamily in polymorphisms[cancerChrom]:
                if (cancerTEFamily == polyTEFamily) and \
                        (polyLeftBP <= cancerLeftBP <= polyRightBP or \
                        polyLeftBP <= cancerRightBP <= polyRightBP):
                    found = True
                    cancerLines.pop(cancerCluster, 0)
                    break

        if found:
            continue

        # Passed all checks, add to calls
        cancerOnlyFile.write(cancerLines[cancerCluster])

    for normalChrom, normalLeftBP, normalRightBP, normalCluster, normalTEFamily in normalBPs:
        found = False
        # Cluster could have been removed while checking cancer breakpoints
        if (normalCluster in normalLines):
            # Check against the cancer discordant ranges
            for discLeftBP, discRightBP, discCluster, discTEFamily in cancerDiscRanges[normalChrom]:
                if (normalTEFamily == discTEFamily) and \
                        (discLeftBP <= normalLeftBP <= discRightBP or \
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
                for polyLeftBP, polyRightBP, polyTEFamily in polymorphisms[normalChrom]:
                    if (normalTEFamily == polyTEFamily) and \
                            (polyLeftBP <= normalLeftBP <= polyRightBP or \
                            polyLeftBP <= normalRightBP <= polyRightBP):
                        found = True
                        normalLines.pop(normalCluster, 0)
                        break

            if found:
                continue

            # Passed all checks, add to calls
            normalOnlyFile.write(normalLines[normalCluster])

    cancerOnlyFile.close()
    normalOnlyFile.close()

def main():
    args = ParseArgs()

    cancerRefinedBPFilename = args.resultsFolder + args.patientID + args.cancerExt + args.callFileExt
    normalRefinedBPFilename = args.resultsFolder + args.patientID + args.normalExt + args.callFileExt
    cancerMappedClustersFilename = args.resultsFolder + args.patientID + args.cancerExt + ".refined.breakpoints.mappedreads.txt"
    normalMappedClustersFilename = args.resultsFolder + args.patientID + args.normalExt + ".refined.breakpoints.mappedreads.txt"
    cancerClustersFilename = args.resultsFolder + args.patientID + args.cancerExt + ".disc.clusters.ranges.txt"
    normalClustersFilename = args.resultsFolder + args.patientID + args.normalExt + ".disc.clusters.ranges.txt"

    chromosomes = [str(x) for x in range(1, 23)] + ['X', 'Y']
    polymorphFilenames = [(chrom, args.polymorphBasename + chrom + ".txt") for chrom in chromosomes]

    polymorphisms = GetPolymorphisms(polymorphFilenames, args.maxRange)

    cancerLines, cancerBPs, header = GetBreakpoints(cancerRefinedBPFilename, args.maxRange)
    normalLines, normalBPs, header = GetBreakpoints(normalRefinedBPFilename, args.maxRange)

    cancerMappedClusters = {}
    normalMappedClusters = {}
    with open(cancerMappedClustersFilename, 'r') as cancerMappedFile:
        for line in cancerMappedFile:
            linesp = line.rstrip('\n').split('\t')
            cancerMappedClusters[linesp[0]] = linesp[1]
    with open(normalMappedClustersFilename, 'r') as normalMappedFile:
        for line in normalMappedFile:
            linesp = line.rstrip('\n').split('\t')
            normalMappedClusters[linesp[0]] = linesp[1]

    cancerDiscRanges = GetDiscRanges(cancerClustersFilename, args.searchRange, cancerMappedClusters)
    normalDiscRanges = GetDiscRanges(normalClustersFilename, args.searchRange, normalMappedClusters)

    CompareBreakpoints(cancerLines, cancerBPs, cancerDiscRanges, normalLines, normalBPs, normalDiscRanges, polymorphisms, args.searchRange, polymorphFilenames, args.resultsFolder, args.patientID, header)

if __name__ == '__main__':
    main()
