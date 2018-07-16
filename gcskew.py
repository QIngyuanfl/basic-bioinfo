#! /usr/bin/env python 

import sys
from Bio import SeqIO
from Bio.SeqUtils import GC

def main(fasta):
    a_num = 0
    t_num = 0
    g_num = 0
    c_num = 0
    p1=''
    p2=''
    p3=''
    for i in SeqIO.parse(open(fasta), 'fasta'):
        length = len(i.seq)
        for j in range(0,len(i.seq)):
            if i.seq[j] == 'G':
                g_num += 1.00
            elif i.seq[j] == 'C':
                c_num += 1.00
            elif i.seq[j] == 'T':
                t_num += 1.00
            elif i.seq[j] == 'A':
                a_num += 1.00
            if j %3 ==1:
                p2+=i.seq[j]        
            if j %3 ==2:
                p3+=i.seq[j]
            if j %3 ==0:
                p1+=i.seq[j]
    #A=p3.count('A')
    #T=p3.count('G')
    #C=p3.count('C')
    #G=p3.count('G')
    #sum=A+T+C+G
    #print (float(G)/float(sum))
        GCskew = round(float(g_num - c_num) / float(g_num + c_num), 4)
        ATskew = round(float(a_num - t_num) / float(a_num + t_num), 4)
        seq = i.id + '\t' + str(GCskew) + '\t' + str(ATskew)
        print seq

if len(sys.argv) > 1:
    main(sys.argv[1])  
