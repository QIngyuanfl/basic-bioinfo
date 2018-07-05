#! /usr/bin/env python
import sys
from Bio import SeqIO
def read_tbl():
    loci=[]
    kind=[]
    info=[]
    gene=[]
    with open(sys.argv[1]) as f:
        count=0
        for line in f:
            if '0'<= line[0]<='9':
                line=line.strip().split()
                count+=1
                if count %2 ==0:
                    loci.append([int(line[0]),int(line[1])])
                if line[2] !='gene':
                    kind.append('type=%s'%line[2])
            else:
                line=line.strip()
                if 'a'<= line[0]<= 'z':
                    line=line.split('\t')
                    if line[0]=='gene':
                        gene.append('[gene=%s]'%line[1])
                    else:
                        if line[0]=='product':
                            info.append('[product=%s]'%line[1])
    return loci,kind,info,gene

def parse_fa(loci,kind,info,gene):
    start=[]
    end=[]
    with open(sys.argv[2]) as f:
        seq_record=SeqIO.read(f,'fasta')
        seq=seq_record.seq
        title=seq_record.id
    start=[x[0] for x in loci]
    end=[x[1] for x in loci]
    with open(sys.argv[3],'w+') as f:
        for i in range (len(kind)):
            if kind[i]== 'type=CDS':
                if start[i] < end[i]:
                    segment=seq[start[i]-1:end[i]]
                    f.write('>%s\t%s\t%s\n%s\n'%(gene[i],info[i],kind[i],segment))
                if start[i] > end[i]:
                    segment=seq[end[i]-1:start[i]].reverse_complement()
                    f.write('>%s\t%s\t%s\n%s\n'%(gene[i],info[i],kind[i],segment))

                    
            
def main():
    loci,kind,info,gene=read_tbl()
    parse_fa(loci,kind,info,gene)
if __name__ == '__main__':
    main()
