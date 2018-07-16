#! /usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
fa=SeqIO.read(sys.argv[1],'fasta')
start=int(sys.argv[3])
end=int(sys.argv[4])
if start < end :
    seq=str(fa.seq[start-1:end])
if start > end  :
    seq=str(fa.seq[end-1:start].reverse_complement())
seq1=Seq(str(seq),IUPAC.unambiguous_dna)
trans=seq1.translate(table=int(sys.argv[2]))
if len(sys.argv) == 5 :
    print trans
else:
    print 'usage:python %s <fa> <trans_table(int)> <start> <end>' % (sys.argv[0])
