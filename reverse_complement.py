#! /usr/bin/env python
import sys
from Bio import SeqIO
with open(sys.argv[1]) as f:
    seq_record = SeqIO.read(f,'fasta')
    title=seq_record.id
    seq=seq_record.seq
with open(sys.argv[2],'w') as f:
    f.write('>%s\n' % title)
    f.write('%s' % seq.reverse_complement())

    
