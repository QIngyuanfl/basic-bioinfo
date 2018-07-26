#!/usr/bin/env python

import os
import re
import sys
from Bio import SeqIO
from argparse import ArgumentParser

__author__ = 'Qingyuan Zhang(zhangqingyuan@scgene.com)'
__version__ = 'V1.0'
__date__ = '18 July 2018'

def read_params(argvs):
    parser = ArgumentParser(description = 'chlorplast genome V2 parser' +__version__+ '(' +__date__ + ')'
                                            'AUTHOR:'+__author__)
    arg = parser.add_argument
    arg ('genbank', type=str, default=None, metavar='input query genbank', help='Genbank file with nucleotide and tRNA sequence')
    arg ('-seqtype', type=str, choices=['tRNA','CDS', 'all'], default='all', metavar = ['tRNA', 'CDS', 'all'], help='CDSs, tRNAs and rRNAs are supported')
    arg('-task', type=str, choices=['cs','cl','bs', 'f', 'l'], default=None, metavar=['cs','cl','bs', 'f', 'l'], help="cs -> calculate gc skew and at skew; cl -> gene classification forms; bs -> base statistics(seqtype CDS is required, f -> fracture gene statistics, l->chlorplast genome, protein coding gene, tRNAs base statistics)" )
    return vars(parser.parse_args())
    
def base_ratio(seq):
    A = seq.count('A')
    T = seq.count('T')
    G = seq.count('G')
    C = seq.count('C')
    length = A + T + G +C
    a = float(A)/length*100
    t = float(T)/length*100
    g = float(G)/length*100
    c = float(C)/length*100
    return a, t, g, c, length

def skew(seq):
    g=seq.count('G')
    c=seq.count('C')
    a=seq.count('A')
    t=seq.count('T')
    GC_skew = float(g-c)/(g+c)
    AT_skew = float(a-t)/(a+t)
    return GC_skew, AT_skew

def genbank_parse(genbank, seqtype):
    seg = {}
    misc = {}
    frac = {}
    single = {}
    repeat = []
    single_CDS = ''
    single_tRNA = ''
    for i in SeqIO.parse(genbank,'genbank'):
        single[i.id] = i.seq
        for seq_feature in i.features:
            if seq_feature.type == 'tRNA':
                single_tRNA += seq_feature.extract(i.seq)
            if seq_feature.type == 'CDS':
                single_CDS += seq_feature.extract(i.seq)
            if seq_feature.type == 'misc_feature':
                misc_feature = seq_feature.qualifiers.values()[0][0]
                geneSeq = seq_feature.extract(i.seq)
                misc[misc_feature] = geneSeq
            if 'join' in str(seq_feature.location):
                frac[seq_feature.qualifiers['gene'][0]] = seq_feature.location
            if 'gene' in seq_feature.qualifiers and len(seq_feature.qualifiers) > 1:
                if 'intron' not in seq_feature.type and 'exon' not in seq_feature.type:
                    gene = seq_feature.qualifiers['gene'][0]
                    gene_location = seq_feature.location
                    geneSeq = seq_feature.extract(i.seq)
                    try:
                        if gene in seg:
                            repeat.append(gene) 
                    except KeyError:
                        pass    
                    seg[gene] = geneSeq
    if seqtype == 'CDS':
        for key in list (seg):
            if 'tRNA' in key or 'rrn' in key:
                seg.pop(key)
    if seqtype == 'tRNA':
        for key in list(seg):
            if not 'tRNA' in key:
                seg.pop(key)
    return seg, misc, repeat, frac, single, single_CDS, single_tRNA     

def gcskew(sequence, args):
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

def gene_classfication(sequence, repeat):
    n_tRNA = 0; n_rRNA =0
    rrn = ''; rps = ''; rpl = ''; rpo = ''; ndh = ''; psa = ''; psb = ''; pet = ''; atp = ''; ycf = ''
    for i in sequence:
        if 'tRNA' in i:
            n_tRNA += 1
        elif 'rrn' in i and i not in repeat:
            rrn += (i + ', ')
        elif 'rps' in i and i not in repeat:
            rps += (i + ', ')
        elif 'rpl' in i and i not in repeat:
            rpl += (i + ', ')
        elif 'rpo' in i and i not in repeat:
            rpo += (i + ', ')
        elif 'ndh' in i and i not in repeat:
            ndh += (i + ', ')
        elif 'psa' in i and i not in repeat:
            psa += (i + ', ')
        elif 'psb' in i and i not in repeat:
            psb += (i + ', ')
        elif 'pet' in i and i not in repeat:
            pet += (i + ', ')
        elif 'atp' in i and i not in repeat:
            atp += (i + ', ')
        elif 'ycf' in i and i not in repeat:
            ycf += (i + ', ')
        else :
            if i not in repeat:
                print i
    for j in repeat:
        if 'tRNA' in j:
            n_tRNA += 1
        if 'rrn' in j:
            rrn += ('%s(x2), ' % (j))
        if 'rps' in j:
            rps += ('%s(x2), ' % (j))
        if 'rpl' in j:
            rpl += ('%s(x2), ' % (j))
        if 'rpo' in j:
            rpo += ('%s(x2), ' % (j))
        if 'ndh' in j:
            ndh += ('%s(x2), ' % (j))
        if 'psa' in j:
            psa += ('%s(x2), ' % (j))
        if 'psb' in j:
            psb += ('%s(x2), ' % (j))
        if 'pet' in j:
            pet += ('%s(x2), ' % (j))
        if 'atp' in j:
            atp += ('%s(x2), ' % (j))
        if 'ycf' in j:
            ycf += ('%s(x2), ' % (j))
    print '%s tRNA genes' % (n_tRNA)
    print rrn
    print rps
    print rpl
    print rpo
    print ndh
    print psa
    print psb
    print pet
    print atp
    print ycf

def section_statistics(misc, seg, cds):
    ta = 0; tt = 0; tc = 0; tg = 0; sum_all = 0
    print 'region', 'A(U)%', 'T(%)', 'G(%)', 'C(%)' , 'length(bp)'
    for i in misc:
        a, t, g, c, sum_base = base_ratio(misc[i])
        print i, '\t', a,'\t', t,'\t', g,'\t', c,'\t', sum_base
        ta += a; tt += t; tg += g; tc += c
        sum_all += sum_base
    a = float(ta)/sum_all*100; t = float(tt)/sum_all*100
    g = float(tg)/sum_all*100; c = float(tc)/sum_all*100
    print 'Total', '\t', a,'\t',t,'\t', g,'\t', c,'\t', sum_all
    p0 = ''; p1 = ''; p2 = ''
    sequence = cds
    A, T, G, C, length = base_ratio(sequence)
    print 'CDS', '\t',A,'\t', T,'\t', G,'\t', C,'\t', length
    for b in range(len(sequence)):
        if b % 3 == 0:
            p0 += sequence[b]
        if b % 3 == 1:
            p1 += sequence[b]
        if b % 3 == 2:
            p2 += sequence[b]
    a0, t0, g0, c0, length = base_ratio(p0)
    print '1st position', '\t', a0,'\t', t0,'\t', g0,'\t', c0,'\t', length
    a1, t1, g1, c1, length = base_ratio(p1)
    print '2ed position', '\t', a1,'\t',t1,'\t', g1,'\t', c1,'\t', length
    a2, t2, g2, c2, length = base_ratio(p2)
    print '3rd position', '\t',a2,'\t', t2,'\t', g2,'\t', c2,'\t', length

def frac_info(genbank):
    gb = SeqIO.read(genbank,'genbank')
    cds = {}
    start = 0; end = 0
    repeat = False
    for i in gb.features:
        if i.type == 'CDS' and 'join' in str(i.location):
            if i.location.start > start and i.location.end < end:
                length = []
                loci = sorted(re.findall(r'\d+\.?\d*' , str(i.location)))
                for j in range(1, len(loci)):
                    length.append(str(int(loci[j]) - int(loci[j-1])))
                print i.qualifiers['gene'][0], '\t', misc, '\t' ,'\t'.join((length))
            else:
                if repeat == False:
                    loci = sorted(re.findall(r'\d+\.?\d*' , str(i.location)))
                    length = []
                    for j in range(1, len(loci), 2):
                        length.append(str(int(loci[j]) - int(loci[j-1])))
                    print i.qualifiers['gene'][0], '\t', 'LSC/IR', '\t','\t'.join((length))
                    repeat = True
        if i.type == 'misc_feature':
            start = i.location.start
            end = i.location.end
            misc = i.qualifiers['note'][0]

def base_statistics(single, CDS, tRNA):
    sequence = str(single.values()[0])
    AT_percent = base_ratio(sequence)[0] + base_ratio(sequence)[1]
    GC_skew, AT_skew = skew(sequence)
    print 'Entire genome, length:%s, AT:%s, AT_skew:%s, GC_skew:%s' % (len(sequence), AT_percent, AT_skew, GC_skew)
    sequence = CDS
    AT_percent = base_ratio(sequence)[0] + base_ratio(sequence)[1]
    p = ''
    for i in range(len(sequence)):
        if i % 3 == 2:
            p += sequence[i]
    AT_percent = base_ratio(sequence)[0] + base_ratio(sequence)[1]
    th_percent = base_ratio(p)[0] + base_ratio(p)[1]
    print 'Protein coding gene, length(aa):%s, AT(all):%s, AT(3rd):%s' %(len(sequence)/3, AT_percent, th_percent)
    sequence = tRNA
    AT_percent = base_ratio(sequence)[0] + base_ratio(sequence)[1]
    GC_skew, AT_skew = skew(sequence)
    print 'tRNA, length:%s, AT:%s, AT-skew:%s, GC-skew:%s' % (len(sequence), AT_percent, AT_skew, GC_skew)
    
def main():
    args = read_params(sys.argv)
    fasta, misc, repeat, frac, single, CDS, tRNA = genbank_parse(args['genbank'], args['seqtype'])
    if args['task'] == 'cs':
        gcskew(fasta, args)
        skew_draw()
    if args['task'] == 'cl':
        gene_classfication(fasta, repeat)
    if args['task'] == 'bs':
        section_statistics(misc, fasta, CDS)
    if args['task'] == 'f':
        frac_info(args['genbank'])
    if args['task'] == 'l':
        base_statistics(single, CDS, tRNA)
if __name__ == "__main__":
    main()
