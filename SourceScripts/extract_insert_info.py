#!/usr/bin/env python

# pipe insert sizes into this program from stdin
import os
import sys
import pysam
import math

def main(args):
  bamfile = pysam.Samfile(args[1], "rb")
  basename = os.path.splitext(os.path.basename(args[1]))[0]
  maxis = 10000 # To prevent artifacts from skewing values
  count = 0
  mean = 0.0
  M2 = 0.0
  for read in bamfile:
    insert_size = int(read.tlen)
    if (insert_size > 0 and insert_size < maxis):
      count = count + 1
      delta = insert_size - mean
      mean = mean + (delta/count)
      M2 = M2 + delta*(insert_size - mean)

  isfile = open("%s.isinfo.txt" % (basename), "w+")
  stddev = math.sqrt(M2/(count - 1))
  isfile.write("%f\n%f" % (mean, stddev))
  isfile.close()

if (__name__ == "__main__"):
  main(sys.argv)
