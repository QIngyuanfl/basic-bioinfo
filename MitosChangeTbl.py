#! /usr/bin/env python

import sys
def MitosChangeTbl(trna):
    change={'A':'T','T':'A','C':'G','G':'C'}
    name=trna[3].split('-')[-1].capitalize()[0:3]
    Code=trna[1].split('(')[-1].split(')')[0].upper()
    revCode=''
    for i in Code:
        revCode+=change[i]
    revCode=revCode[::-1]
    GeneName='tRNA-%s(%s)' % (name,Code)
    revGeneName='tRNA-%s(%s)' % (name,revCode)
    trna[1]='gene\t%s' % (revGeneName)
    trna[3]='product\t%s' % (GeneName)
    return trna
def main(tbl,out,table):
    info={}
    gene=[]
    tblf=open(tbl,'r')
    outf=open(out,'w')
    for Line in tblf:
        if Line[0] !='>':
            Split=Line.strip().split()
            gene.append(Line.strip())
            if Split[0] == 'product':
                start=int(gene[0].split()[0])
                if gene[2].split()[-1] =='tRNA':
                    gene=MitosChangeTbl(gene)
                if gene[-2].split()[-1] == 'rRNA':
                    gene=[gene[0],'gene\tRNA',gene[1],gene[2]]
                if gene[2].split()[-1] =='CDS':
                    gene.append('transl_table\t%s' % table)
                info[start]=gene
                gene=[]
        if Line[0] == '>':
            outf.write('%s\n' % Line.strip())
    size=info.keys()
    size.sort()
    for i in size:
        if len(info[i]) ==4:
            outf.write('%s\n\t\t\t%s\n%s\n\t\t\t%s\n' % (info[i][0],info[i][1],info[i][2],info[i][3]))
        if len(info[i]) ==5:
            outf.write('%s\n\t\t\t%s\n%s\n\t\t\t%s\n\t\t\t%s\n' % (info[i][0],info[i][1],info[i][2],info[i][3],info[i][4]))
if len(sys.argv) ==4 :
    main(sys.argv[1],sys.argv[2],sys.argv[3])
if len(sys.argv) !=4 :
    print 'Usage:python %s MitosTbl output transl_table' % sys.argv[0]
