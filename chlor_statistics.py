#!/usr/bin/env python

import os
import sys
from Bio import SeqIO
from argparse import ArgumentParser

__author__ = 'Qingyuan Zhang(zhangqingyuan@scgene.com)'
__version__ = '0.1.0'
__date__ = '18 July 2018'

def read_params(argvs):
    parser = ArgumentParser(description = 'chlorplast genome V2 parser' +__version__+ '(' +__date__ + ')'
                                            'AUTHOR:'+__author__)
    arg = parser.add_argument
    arg ('fasta', type=str, default=None, metavar='input nucl fasta', help='Fasta file with nucleotide and tRNA sequence')
    arg ('-seqtype', type=str, choices=['tRNA','CDS'], default='None', metavar='type of sequence', help='CDSs, tRNAs and rRNAs are supported')
    arg('-task', type=str, choices=['cs'], default=None, metavar='type of analyse', help="cs\t calculate and draw CDSs' gc skew" )
    return vars(parser.parse_args())
    

def seq_parse(fasta, seqtype):
    seg={}
    with open(fasta) as nucl:
        for i in SeqIO.parse(nucl,'fasta'):
            try:
                while i.id in seg:
                    i.id += '+' 
            except KeyError:
                pass
            seg[i.id] = i.seq
    if seqtype == 'CDS':
        for key in list(seg):
            if 'tRNA' in key or 'rrn' in key:
                seg.pop(key)
    if seqtype == 'tRNA':
        for key in list(seg):
            if not 'tRNA' in key:
                seg.pop(key)

    return seg

def gcskew(sequence,args):
    if args['seqtype'] == 'CDS':
        f=open('gc_skew.csv','w')
        f.write('Pc' + '\t' + 'GC_skew' + '\t' + 'AT_skew' + '\n')
    G = 0; C = 0; A = 0; T = 0
    for i in sorted(sequence):
        g=sequence[i].count('G')
        c=sequence[i].count('C')
        a=sequence[i].count('A')
        t=sequence[i].count('T')
        GC_skew = float(g-c)/(g+c)
        AT_skew = float(a-t)/(a+t)
        G += g; C += c; A += a; T += t
        try:
            f.write('%s\t%s\t%s\n' % (i, GC_skew, AT_skew))
        except UnboundLocalError:
            pass
    try:
        f.close()
    except UnboundLocalError:
        pass
    GC_skew = float(G-C)/(G+C)
    AT_skew = float(A-T)/(A+T)
    if args['seqtype'] == 'tRNA':
        print 'Total', '\t', GC_skew, '\t', AT_skew

def skew_draw():
    csv = {}
    n = 0
    breaks = ''

    with open('gc_skew.csv') as f:
        f.next()
        for line in f:
            line=line.strip().split()
            line[0]='"'+line[0]+'"'
            csv[line[0]]=[line[1],line[2]]
            n += 1
    Pc = ', '.join(sorted(csv))
    for i in range (1,n):
        breaks += ('%s,' % i)
    breaks += str(n)
        
    with open('gc_skew.R','w') as f:
        f.write('library(ggplot2)\n')
        f.write('library(reshape2)\n')
        f.write('gc <- read.csv(file="gc_skew.csv", header=TRUE, sep="\t")\n')
        f.write('pdf("gc_skew.pdf",width =35,height=15)\n')
        f.write('gc_long <- melt(gc, id = "Pc", value= "skew")\n')
        f.write('ggplot(gc_long,aes(x=as.numeric(Pc),y=value,colour=variable)) + geom_line() +geom_point()+labs (y= "value", x= "Gene", colour="skew") +scale_x_continuous(limits=c(1,%s),breaks=c(%s),labels=c(%s))' % (n, breaks, Pc))

    os.system('Rscript gc_skew.R')

#def gene_classfication(sequence):
    
    

def main():
    args = read_params(sys.argv)
    fasta = seq_parse(args['fasta'], args['seqtype'])
    if args['task'] == 'cs':
        gcskew(fasta,args)
        skew_draw()
if __name__ == "__main__":
    main()
