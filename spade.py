#! /usr/bin/env python

import os
import sys
from Bio import SeqIO

def spade(bam,fa,fq1,fq2,start,end):
    fa_dic={}
    os.system("samtools view %s > bam.txt" % (bam))
    fq_list={}
    for i in SeqIO.parse(fa,'fasta'):
        fa_dic[i.id]=len(i.seq)
    if start == 1 :
        for i in open("bam.txt",'r'):
            line=i.strip().split()
            if int(line[3]) < 500:
                fq_list[line[0]]=0
            for k in fa_dic.keys():
                if line[2] == k :
                    if int(line[3]) > (fa_dic[k] - 700):
                        fq_list[line[0]]=0
    if start != 1:
        if len(fa_dic) >1 :
            print "If the number of scaffold bigger than 1 ,this program only extend start and end"
            print "While the number of scaffold is 1 ,this program will repair the size you choice"
            sys.exit()
        else:
            for i in open("bam.txt",'r'):
                line=i.strip().split()
                if  int(start) < int(line[3]) < int(end):
                    fq_list[line[0]]=0
    Out=open("list.txt",'w')
    for i in fq_list.keys():
        Out.write('%s\n' % (i))
    Out.close()
    os.system("seqtk subseq %s list.txt > R1.fq" % (fq1))
    os.system("seqtk subseq %s list.txt > R2.fq" % (fq2))
    #os.system("rm -f bam.txt list.txt ")
    Out1=open("spade.sh",'w')
    Out1.write("#$ -S /bin/sh\n/hellogene/scgene01/opt/bin/python /hellogene/scgene01/bio/pip/ass/SPAdes-3.5.0-Linux/bin/spades.py -k 79,97 --only-assembler -1 R1.fq -2 R2.fq -t 15 -m 30 -o spadeAss")
    Out1.close()
    os.system("sh spade.sh")
    Out2=open("scaffold_long.fa",'w')
    for i in SeqIO.parse("spadeAss/scaffolds.fasta",'fasta'):
        Seq=i.seq
        revSeq=i.seq.reverse_complement()
        Out2.write(">%s\n%s\n>%s_rev\n%s" % (i.id,Seq,i.id,revSeq))
    Out.close()
if len(sys.argv) == 1 :
    print "usage:python %s <bam> <fa> <fq1> <fq2> (start) (end)"
if len(sys.argv) ==5:
    start=1
    end=1
    spade(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],start,end)
if len(sys.argv) ==7:
    spade(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
