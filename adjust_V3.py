# -*- coding: utf-8 -*-
#!/usr/bin/env python

#使用原始fasta文件找到调整区cox1位置后（将cox1作为fasta文件起始），
#此程序可以将原始fasta文件切成需要的fasta文件, 将tbl文件作相应调整
#使用方法：1.找到cox1的起始位点，未来第1点的位置（如果cox1为正向则是最小值，如果cox1是反向则为最大值）
#	  2.在工作目录下建立一个input文件，将fasta文件和tbl文件移到input文件中
#	   3.在input文件夹外（工作目录下），直接运行python Adjust_v2.py会有提示
#	    4.按照提示输入后会在input文件夹中产生*.out文件，即为需要文件
#V3在V2的基础上修改了start的算法，稍微优化代码，结果不变
#last edit: 2015.12.16
#author: lurui@scgene.com

import re
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
def AdjustSeq(fsafile, start, reverse):
        ''' output fsafile which adjusted. at the same time,
                record the whole length of seq of fsafile'''
        fp = SeqIO.read(fsafile, 'fasta')
        seqlen = len(fp.seq)        
        if re.match(r'^[nN]', reverse):
                newseq = fp.seq[(start-1):] + fp.seq[:(start-1)]
                seq = newseq
        elif re.match(r'^[yY]', reverse):
                newseq = fp.seq[start:] + fp.seq[:start]
                seq = newseq.reverse_complement()
        else:
                print '-*- 请输入yes或no !! -*-'
                sys.exit(1)
        newfsa = open(fsafile+'.out', 'w')
        newfsa.write('>'+fp.description+'\n')        
        newfsa.write(str(seq))
        newfsa.close()
        return seqlen
        
def AdjustTbl(tblfile, start, seqlen, reverse):
        '''可改写为读一行写一行'''
        newtbl = open(tblfile+'.out', 'w')
        infile = open(tblfile)        
        origintbl = infile.readlines()
        for eachline in origintbl:
                newline = ''
                llist = []
                llist = eachline.strip().split('\t')
                if len(llist)!=0 and re.match(r'^[0-9]', llist[0]):
                        if re.match(r'^[nN]', reverse):
                                t = start - 1 
                                if int(llist[0]) > t and int(llist[1]) > t:
                                        llist[0] = str(int(llist[0]) - t)
                                        llist[1] = str(int(llist[1]) - t)
                                elif int(llist[0]) < t and int(llist[1]) < t:
                                        llist[0] = str(seqlen - t + int(llist[0]))
                                        llist[1] = str(seqlen - t + int(llist[1]))
                                else:
                                        print '-*-起始区cox1有重叠，请先删除重叠CDS，再运行-*-'
                                        sys.exit(1)
                        elif re.match(r'^[yY]', reverse):
                                s = start + 1
                                if int(llist[0]) > s and int(llist[1]) > s:
                                        llist[0] = str(seqlen + s - int(llist[0]))
                                        llist[1] = str(seqlen + s - int(llist[1]))
                                elif int(llist[0]) < s and int(llist[1]) < s:
                                        llist[0] = str(s - int(llist[0]))
                                        llist[1] = str(s - int(llist[1]))
                                else:
                                        print '-*-起始区cox1有重叠，请先删除重叠CDS，再运行-*-'
                                        sys.exit(1)                                
                        else:
                                print '-*- 请输入yes或no !! -*-'
                                sys.exit(1)
                        newline = '\t'.join(llist) + '\n'
                        newtbl.write(newline)
                else:
                        newtbl.write(eachline)
        infile.close()
        newtbl.close()

def main():
        fname = os.listdir(os.getcwd())[:]
        fsafile = ''
        tblfile = ''
        for i in fname:
                if re.search(r'(.fsa)|(.fasta)|(.fa)',i):
                        fsafile = i
                elif re.search(r'(.tbl)',i):
                        tblfile = i
                else:
                        continue
        if fsafile == '' or tblfile == '':
                print '请将原始fasta和tbl文件放在工作目录下'
                sys.exit(1)                
        start = int(raw_input('请输入COX1的起始位点(M端)：'))
#        start = start -1
        reverse = raw_input('轻重链是否翻转(yes/no)?')
#        fsafile = 'input/' + fsafile
#        tblfile = 'input/' + tblfile
        seqlen = AdjustSeq(fsafile, start, reverse)
        AdjustTbl(tblfile, start, seqlen, reverse)

if __name__ == '__main__':
        main()
                
