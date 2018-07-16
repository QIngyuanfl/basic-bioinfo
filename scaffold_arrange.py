#! /usr/bin/env python

import sys
import re

def main(fasta):
    newfasta = {}
    for line in open(fasta, 'r'):
        line = line.rstrip()
        if re.match(r'^>', line):
            seqname = line
            newfasta[seqname] = []
        else:
            newfasta[seqname].append(line)

    for line in newfasta:
        print '%s\n%s' % (line, ''.join(newfasta[line]))

if len(sys.argv) > 1:
    main(sys.argv[1])
