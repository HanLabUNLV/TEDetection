#!/usr/bin/env python
import sys

reffilename = sys.argv[1]

TEgroups = {}

with open(reffilename, 'r') as reffile:
  for line in reffile:
    if (line[0] == '>'):
      line_sp = line.split('\t')
      line_sp[0] = line_sp[0][1:]
      if (line_sp[1] in TEgroups):
        TEgroups[line_sp[1]].append(line_sp[0])
      else:
        TEgroups[line_sp[1]] = []
        TEgroups[line_sp[1]].append(line_sp[0])

outfile = open("grouped.TEs.txt", 'w')
for key,value in TEgroups.iteritems():
  group = '\t'.join(value)
  outfile.write(key + '\t' + group + '\n')

outfile.close()
