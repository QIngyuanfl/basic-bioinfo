#! /usr/bin/env python
# -*- coding: utf-8 -*-
'''''''''''''''''''''''''''''''''''''''''''''
此程序是为了提取出断裂区域的序列，可把出来的文件直接覆盖V1的Feature的
CDS和nucl两个文件。
'''''''''''''''''''''''''''''''''''''''''''''
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
scgene,faf,tblf=sys.argv
def main(faf,tblf):
	fa=SeqIO.read(faf,'fasta')
	seq=str(fa.seq)
	tbl=open(tblf,'r')
	out=open("CDS.tbl",'w')
	out1=open("nucl.tbl",'w')
	tlist=[]
	for line in tbl:
		tlist.append(line.strip().split())
	tlist=tlist[1:]
	list=[]
	for i in tlist:
		if i != []:
			list.append(i)
	a=[i for i,x in enumerate(list) if list[i][0]=='gene']
	b=[i for i ,x in enumerate(list) if list[i][0]=='product']
	for i in range(len(b)):
		c=int(a[i])
		d=int(b[i])
		type=list[c+1][2]
		gene=list[c][1]
		product="%s=%s" %(list[d][0],' '.join(list[d][1:]))
		if int(d)-int(c)<=2:
			if int(list[c-1][0]) <int(list[c-1][1]):
				start=int(list[c-1][0])-1
				end=int(list[c-1][1])
				seq2=str(seq[start:end])
				start1=start+1
				end1=end
			if int(list[c-1][1]) <int(list[c-1][0]):
				start=int(list[c-1][1])-1
				end=int(list[c-1][0])
				seq1=Seq(seq[start:end],IUPAC.unambiguous_dna)
				seq2=seq1.reverse_complement()
				start1="complment[("+str(start+1)+".."
				end1=str(end)+")]"
			product="%s=%s" %(list[d][0],' '.join(list[d][1:]))
			out1.write(">%s [%s] [(%s..%s)]\n%s\n" %(gene,product,start1,end1,seq2))
			if type=="CDS":
				out.write(">%s [%s] [(%s..%s)]\n%s\n" %(gene,product,start1,end1,seq2))
		if int(d)-int(c)>2:
			e=[]
			f=[]
			z=[]
			g="no"
			type=list[c+1][2]
			gene=list[c][1]
			product="%s=%s" %(list[d][0],' '.join(list[d][1:]))
			for x in range(c+1,d,1):
				if int(list[x][0]) < int(list[x][1]):
					start=int(list[x][0])-1
					end=int(list[x][1])
					f.append(seq[start:end])
					e.append("(%s..%s)"  % (start+1,end))
				if int(list[x][0]) > int(list[x][1]):
					start=int(list[x][1])-1
					end=int(list[x][0])
					seq1=Seq(seq[start:end],IUPAC.unambiguous_dna)
					seq2=seq1.reverse_complement()
					f.append(str(seq2))
					z.append([start+1,end])
					z.sort()
					g='yes'
			for y in z :
				e.append('(%s..%s)' % (y[0],y[1]))
			seq3=str(''.join(f))
			sub=' '.join(e)
			if g =="no":
				out1.write(">%s [%s] [%s]\n%s\n" %(gene,product,sub,seq3))
				if type=="CDS":
					out.write(">%s [%s] [%s]\n%s\n" %(gene,product,sub,seq3))
			if g =="yes":
				out1.write(">%s [%s] complment[%s]\n%s\n" %(gene,product,sub,seq3))
				if type=="CDS":
					out.write(">%s [%s] complment[%s]\n%s\n" %(gene,product,sub,seq3))
	out.close()
	out1.close()
main(faf,tblf)
