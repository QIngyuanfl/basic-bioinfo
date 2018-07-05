#! /usr/bin/env python

import sys
from Bio import SeqIO

Seq=SeqIO.read(sys.argv[1],'fasta')
fa=Seq.seq
Gene=0
f=open('ref.fasta','w')
for i in open(sys.argv[2],'r'):
    line=i.strip().split('\t')
    if len(line)!=0 and line[0][0]!='#' and line[2] != 'gene' and line[2] != 'exon':
        start=int(line[3])
        end=int(line[4])
        strange=line[6]
        seq=fa[start-1:end]
        if strange=='+':
            result_seq=seq
        if strange=='-':
            result_seq=seq.reverse_complement()
        info=line[8]
        ID=info.split(';')[0].split('=')[1]
        if line[2] =='CDS':
            gene=info.split('gene=')[-1].split(';')[0]
            product=info.split('product=')[-1].split(';')[0]
            if product != 'hypothetical protein':
                if "pseudo=true" in i.strip():
                    note=info.split('Note=')[-1].split(';')[0]
                    if 'ID=' not in gene:
                        f.write('>%s\t[gene=%s]\t[note=%s;product=%s]\ttype=CDS\t%s\n%s' % (ID,gene,note,product,len(result_seq),result_seq))
                    else:
                        Gene+=1
                    
                        f.write('>%s\t[gene=Gene%s]\t[note=%s;product=%s]\ttype=CDS\t%s\n%s' % (ID,Gene,note,product,len(result_seq),result_seq))
                if "pseudo=true" not in i.strip():
                    if 'ID=' not in gene:
                        f.write('>%s\t[gene=%s]\t[product=%s]\ttype=CDS\t%s\n%s' % (ID,gene,product,len(result_seq),result_seq))
                    else:
                        Gene+=1
                        f.write('>%s\t[gene=Gene%s]\t[product=%s]\ttype=CDS\t%s\n%s' % (ID,Gene,product,len(result_seq),result_seq))
        if line[2] =='tRNA':
            product=info.split('product=')[-1].split(';')[0]
            annit=info.split('Note=')[-1].split(';')[0]
            f.write('>%s-tRNA\t[gene=%s]\t[product=%s]\ttype=tRNA\n%s' % (ID,product,product,result_seq))
        if line[2] =='rRNA':
            gene=info.split('gene=')[-1].split(';')[0]
            product=info.split('product=')[-1].split(';')[0]
            if 'ID=' not in gene:
                f.write('>%s-rRNA:%s\t[gene=%s]\t[product=%s]\ttype=rRNA\n%s' % (ID,len(result_seq),gene,product,result_seq))
            else:
                Gene+=1
                f.write('>%s-rRNA:%s\t[gene=Gene]\t[product=%s]\ttype=rRNA\n%s' % (ID,len(result_seq),product,result_seq))
        '''
        if line[2] == 'repeat_region':
            gene='repeat_region'
            product=info.split('Note=')[-1].split(';')[0]
            f.write( '>%s\t[gene=%s]\t[note=%s]\ttype=repeat_region\n%s' % (ID,gene,product,result_seq)
        if line[2] =='ncRNA' or line[2] == 'antisense_RNA': 
            gene=info.split('gene=')[-1].split(';')[0]
            product=info.split('product=')[-1].split(';')[0]
            f.write( '>%s\t[gene=%s]\t[product=%s]\ttype=ncRNA\n%s' % (ID,gene,product,result_seq)
        if line[2] =='mobile_genetic_element':
            gene='mobile_element'
            product=info.split('=')[-1]
            f.write( '>%s\t[gene=%s]\t[mobile_element_type=%s]\ttype=mobile_element\n%s' % (ID,gene,product,result_seq)
        '''
