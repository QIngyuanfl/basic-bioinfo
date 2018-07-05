#!/usr/bin/env python
import sys
import os
from Bio import SeqIO

(xx,fname,a,b) = sys.argv
with open(fname) as f:
    read=f.read()
if read.count('>') ==1:
    fp=SeqIO.read(fname,'fasta')
    data=fp.seq
if read.count('>') >1:
    title=raw_input('input query title:') 
    for seq_record in SeqIO.parse(fname,'fasta'):
        if title==seq_record.id:
            data=seq_record.seq
if int(b) > int(a):
		start=int(a)-1
		stop=int(b)
		print (str(data[start:stop]),start+1,stop,'seqlen%%3=%d' %((stop-start)%3))
if int(a) > int(b):
		start=int(b)-1
		stop=int(a)
		print (str(data[start:stop].reverse_complement()),start+1,stop,'seqlen%%3=%d' %((stop-start)%3))
# print data[a-100:a],'***100nt Before start***'  
# print data[b:b+100],'***100nt After end***'

