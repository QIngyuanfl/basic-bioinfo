#!/usr/bin/env python
import sys
from Bio import SeqIO
def read_fa():
    seq={}
    with open(sys.argv[1]) as f:
        seq_record=SeqIO.read(f,'fasta')
        seq[seq_record.id] = seq_record.seq
    return seq

def parse_fa(seq):
    start=int(sys.argv[2])-1
    ID=seq.keys()[0]
    sequence=seq.values()[0]
    length=len(sequence)
    with open ('new_%s'% sys.argv[1],'w') as f:
        f.write('>%s' %ID)
        f.write('\n')
        f.write('%s' % (sequence[start:length]))
        f.write('%s' % sequence[0:start])
        
    
def main():
    seq=read_fa()
    parse_fa(seq)
if __name__ == '__main__':
    main()
