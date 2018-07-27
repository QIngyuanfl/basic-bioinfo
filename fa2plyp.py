#! /usr/bin/env python
import os
import sys
import re
(xx,a,b)=sys.argv
with open(a ,'r') as fin:
    sequences = [(m.group(1), ''.join(m.group(2).split()))
    for m in re.finditer(r'(?m)^>([^ \n]+)[^\n]*([^>]*)', fin.read())]
with open(b, 'w') as fout:
    fout.write('%d %d\n' % (len(sequences), len(sequences[0][1])))
    for item in sequences:
        fout.write('%-14s %s\n' % item)
