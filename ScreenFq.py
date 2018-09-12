#!/usr/bin/env python

import sys
import gzip
from Bio import SeqIO
from functools import partial
from multiprocessing import Pool

def ScreenFq(fq):
    ID = []
    n = 0
    with gzip.open(fq) as f:
        for seq_record in SeqIO.parse(f, 'fastq'):
            n += 1
            if 'N' in seq_record.seq:
                continue
            count = 0
            for i in seq_record.letter_annotations['phred_quality']:
                if i < 20:
                    count += 1
                if count >= 3:
                    continue
            ID.append(seq_record.id)
    return ID

def PassFq(fq, reads):
    name = fq.split('.')[0]+'.fq'
    with open(name, 'w') as o:
        with gzip.open(fq) as f:
            for seq_record in SeqIO.parse(f, 'fastq'):
                if seq_record.id in reads:
                    SeqIO.write(seq_record, o ,'fastq')

def main(fq1, fq2):
    FqFile = []
    for i in sys.argv[1:]:
        FqFile.append(i)
    p = Pool(2)
    ID1, ID2 = p.map(ScreenFq,FqFile)
    CommonID = list(set(ID1) & set(ID2))
    del ID1; del ID2
    p = Pool(16)
    split_length = len(CommonID)/128
    reads = [CommonID[x*split_length: (x+1)*(split_lenth-1)] for x in xrange(128)]
    p.map(partial(PassFq, reads = CommonID), FqFile)
    p.terminate()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'usage python ScreenFq.py <fq1> <fq2>'
    else:
        main(sys.argv[1], sys.argv[2])
