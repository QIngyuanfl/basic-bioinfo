#!/usr/bin/env python
from Bio import SeqIO
from Bio import SeqUtils
import sys
from argparse import ArgumentParser
__author__ = 'Qingyuan Zhang(zhangqingyuan@scgene.com)'
__version__ = '0.1.0'
__date__ = '18 July 2018'

def read_params(argvs):
    parser = ArgumentParser(description = 'chlorplast genome V2 parser' +__version__+ '(' +__date__ + ')'
                                            'AUTHOR:'+__author__)
    arg = parser.add_argument
    arg ('fasta', type=str, default=None, metavar='input nucl fasta', help='Fasta file with nucleotide and tRNA sequence')
    arg ('-seqtype', type=str, choices=['tRNA','CDS'], default='CDS', metavar='type of sequence', help='The default sequence type is CDS, tRNAs and rRNAs are supported')
    return vars(parser.parse_args())
    

def seq_parse(fasta, seqtype):
    seg={}
    with open(fasta) as nucl:
        for i in SeqIO.parse(nucl,'fasta'):
            seg[i.id]=i.seq
    if seqtype == 'CDS':
        for key in list(seg):
            if 'tRNA' in key or 'rrn' in key:
                seg.pop(key)
    if seqtype == 'tRNA':
        for key in list(seg):
            if not 'tRNA' in key:
                seg.pop(key)
    return seg

def gcskew(seg):
    print 'PE' + '\t' + 'GC_skew' + '\t' + 'AC_skew' + '\t'
    for i in seg:
        gc_skew=SeqUtils.GC_skew(seg[i],window=100)
        print i , '\t' , gc_skew[0] , '\t'
    
    
                
def main():
    args = read_params(sys.argv)
    fasta=seq_parse(args['fasta'], args['seqtype'])
    gcskew(fasta)
if __name__ == "__main__":
    main()
