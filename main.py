#!/usr/bin/env python
import os
import os.path
import sys
import datetime
import subprocess

class switch(object):
  def __init__(self, value):
    self.value = value
    self.fall = False

  def __iter__(self):
    yield self.match
    raise StopIteration

  def match(self, *args):
    if self.fall or not args:
      return True
    elif self.value in args:
      self.fall = True
      return True
    else:
      return False

def parseArgs(sysargs):
  usage = "Usage: ./TEDetection\n \
            -sd [Path/to/SourceScripts]\n \
            -bf [Path/to/BamFile.bam]\n \
            -tr [Path/to/TErefgen.fasta]\n \
            -bb [Path/to/SortedByNameBamfile.bam] (will be created if left empty)\n \
            -is [Path/to/isinfo.txt]  mean + std dev of insert sizes (will be created if left empty)\n \
            -db [Path/to/DiscordantPairs.bam] coordinate sorted file of discordant pairs (will be created if left empty)\n \
            -dq [Path/to/DiscordantBynamePairs.bam] queryname sorted file of discordant pairs (will be created if left empty)\n \
            -di [Path/to/DiscordantPairs.index] index file for discpairs.bam (will be created if left empty)\n \
            -cf [Path/to/ClustersFile.txt] (will be created if left empty)\n \
            -cr [Path/to/ClusterRanges.txt] (will be created if left empty)\n \
            -cs [Path/to/ClusterMapped.sam] (will be created if left empty)\n \
            -ss [Path/to/SoftclipsMapped.sam] (will be created if left empty)\n \
            -bp [Path/to/Breakpoints.txt] (will be created if left empty)\n \
            -tg [Path/to/TEGroups.txt]  File with tab seperate groups of similar TEs (optional)\n \
            -rd <int> default:50   Average read depth (used for filtering sections of abnormally large read depth)\n \
            -mc <int> default:5    Minimum cluster size for discordant reads\n \
            -ms <int> default:3    Minimum softclips supporting a breakpoint\n \
            -sr <int> default:75  Search range for breakpoints\n \
            -st <int> default:1    Number of threads for sorting\n \
            -sm <int>(bytes) default: 1073741824 (1 GB) Amount of memory per thread used in sorting\n \
            -d  Discard intermediate files when finished\n \
            -h/--help This usage output"

  args = {}
  # Default values
  args["-rd"] = 50
  args["-mc"] = 5
  args["-ms"] = 3
  args["-sr"] = 75
  args["-st"] = 1
  args["-sm"] = 1073741824
  args["-d"] = False
  i = 1
  default = True
  if (len(sysargs) == 1):
    print usage
    quit()
  while (i < len(sysargs)):
    for case in switch(sysargs[i]):
      if case("-d"):
        opt = sysargs[i]
        args[opt] = True
        i += 1
        break
      if case("-h") or case("--help") or case (""):
        print usage
        quit()

      if (default):
        opt = sysargs[i]
        i += 1
        args[opt] = sysargs[i]
        i += 1
        break

  return args

def main(args):
  argdict = parseArgs(args)
  sys.path.append(argdict["-sd"])
  import cluster_discordant_pairs
  import compare_mapped_reads
  import extract_discordant_pairs

  try:
    os.makedirs("Scratch")
  except OSError:
    print "Scratch folder already exists"

  try:
    os.makedirs("Results")
  except OSError:
    print "Results folder already exists"

  if ("-bf" in argdict):
    basename = os.path.splitext(os.path.basename(argdict["-bf"]))[0]
  else:
    print "Bamfile option required for program (-bf bamfile.bam)"
    print "Exiting now.."
    quit()

  buffersize = 0
  outputfile = open("%s.output" % (basename), 'w+', buffersize)

  outputfile.write("Starting TEDetection program.." + '\n')
  outputfile.write(str(datetime.datetime.today()) + '\n')

  collateran = False

  if (not "-bb" in argdict):
    outputfile.write("Sorting bamfile by query name.." + '\n')
    outputfile.write(str(datetime.datetime.today()) + '\n')
    collateran = True
    collatecommand = "bamcollate2 filename=%s O=Scratch/%s.byname.bam outputthreads=%s inputbuffersize=%s 2>> bamcollate.err" % (argdict["-bf"], basename, argdict["-st"], argdict["-sm"])
    collateproc = subprocess.Popen(collatecommand, stdout=subprocess.PIPE, shell=True)
    argdict["-bb"] = "Scratch/%s.byname.bam" % (basename)

  isinforan = False

  if (not "-is" in argdict):
    outputfile.write("Calculating mean and stddev of insert sizes.." + '\n')
    outputfile.write(str(datetime.datetime.today()) + '\n')
    isinforan = True
    getiscommand = "python %s/extract_insert_info.py %s" % (argdict["-sd"], argdict["-bf"])
    getisproc = subprocess.Popen(getiscommand, stdout=subprocess.PIPE, shell=True)
    argdict["-is"] = "%s.isinfo.txt" % (basename)

  if (collateran):
    collateproc.wait()
  if (isinforan):
    getisproc.wait()

  if (not "-db" in argdict or not "-di" in argdict or not "-dq" in argdict):
    outputfile.write("Extracting discordant read pairs.." + '\n')
    outputfile.write(str(datetime.datetime.today()) + '\n')

    extract_discordant_pairs.main(["none", argdict["-bb"], argdict["-is"]])

    discsortcommand = "samtools sort Scratch/%s.disc.bam Scratch/%s.disc" % (os.path.splitext(os.path.basename(argdict["-bb"]))[0], basename)
    discsortproc = subprocess.Popen(discsortcommand, stdout=subprocess.PIPE, shell=True)
    discsortproc.wait()

    discindexcommand = "samtools index Scratch/%s.disc.bam" % (basename)
    discindexproc = subprocess.Popen(discindexcommand, stdout=subprocess.PIPE, shell=True)
    discindexproc.wait()
    argdict["-db"] = "Scratch/%s.disc.bam" % (basename)
    argdict["-di"] = "Scratch/%s.disc.index" % (os.path.splitext(os.path.basename(argdict["-bb"]))[0])
    argdict["-dq"] = "Scratch/%s.disc.bam" % (os.path.splitext(os.path.basename(argdict["-bb"]))[0])

  if (not "-cf" in argdict or not "-cr" in argdict):
    outputfile.write("Clustering discordant reads.." + '\n')
    outputfile.write(str(datetime.datetime.today()) + '\n')

    cluster_discordant_pairs.main(["none", argdict["-db"], argdict["-mc"], argdict["-rd"], argdict["-is"]])
    argdict["-cf"] = "Results/%s.clusters.txt" % (os.path.splitext(os.path.basename(argdict["-db"]))[0])
    argdict["-cr"] = "Results/%s.clusters.ranges.txt" % (os.path.splitext(os.path.basename(argdict["-db"]))[0])

  mapclusran = False
  if (not "-cs" in argdict):
    outputfile.write("Mapping cluster mates to TE reference genome.." + '\n')
    outputfile.write(str(datetime.datetime.today()) + '\n')
    genmatefastqcommand = "python %s/generate_mate_fastq.py %s %s %s %s %s" % (argdict["-sd"], argdict["-dq"], argdict["-di"], argdict["-cf"], argdict["-tr"], argdict["-st"])
    genmatefastqproc = subprocess.Popen(genmatefastqcommand, stdout=subprocess.PIPE, shell=True)
    argdict["-cs"] = "Scratch/%s.clusters.sam" % (os.path.splitext(os.path.basename(argdict["-di"]))[0])
    mapclusran = True

  mapscran = False
  if (not "-ss" in argdict or not "-bp" in argdict):
    outputfile.write("Mapping softclip reads to TE reference genome.." + '\n')
    outputfile.write(str(datetime.datetime.today()) + '\n')
    extractsccommand = "python  %s/extract_softclips.py %s %s %s %s %s %s" % (argdict["-sd"], argdict["-bf"], argdict["-cr"], argdict["-sr"], argdict["-rd"], argdict["-tr"], argdict["-st"])
    extractscproc = subprocess.Popen(extractsccommand, stdout=subprocess.PIPE, shell=True)
    argdict["-ss"] = "Scratch/%s.softclips.sam" % (basename)
    argdict["-bp"] = "Results/%s.breakpoints.txt" % (basename)
    mapscran = True

  if (mapclusran):
    genmatefastqproc.wait()
  if (mapscran):
    extractscproc.wait()

  refinedoutputfilename = "Results/%s.refined.breakpoints.txt" % (basename)
  outputfile.write("Comparing mapped cluster mates to mapped softclip sequences.." + '\n')
  outputfile.write(str(datetime.datetime.today()) + '\n')
  compare_mapped_reads.main(["none", argdict["-cs"], argdict["-ss"], argdict["-bp"], refinedoutputfilename, argdict["-tg"]])

  if (argdict["-d"] == True):
    outputfile.write("Removing temp data from scratch directory.." + '\n')
    rmcommand = "rm Scratch/%s*" % (basename)
    rmproc = subprocess.Popen(rmcommand, stdout=subprocess.PIPE, shell=True)
    rmproc.wait()

  outputfile.write("TEDetection program finished" + '\n')
  outputfile.write(str(datetime.datetime.today()) + '\n')

if (__name__ == "__main__"):
  main(sys.argv)
