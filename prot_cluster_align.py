#! /usr/bin/env python
import sys
import os
def get_ref_ID():
    with open(sys.argv[1]) as f:
        ID=[]
        for line in f:
            if line[0]=='>':
                ID.append(line.strip().lstrip('>').split('_i')[0])
        ID=list(set(ID))
    return ID
def find_blastp(ID):
    bp=os.listdir('blastp/')
    with open('prot_alignments.xls','w') as f:
        f.write('\t%s\n'% ('\t'.join(ID)))
        for i in bp:
            with open('blastp/%s'%i) as blastp:
                blastp.seek(0)
                f.write('%s\t' % i.split('.')[0])
                align={}
                for x in ID:
                    align.setdefault(x,'#N/A')
                for lines in blastp:
                    lines=lines.strip().split()
                    if float(lines[2])>=90:
                        query=lines[0].split('_i')[0]
                        if align[query] != '#N/A' and lines[1] not in align[query]:
                            align[query]+=(';%s'%lines[1])
                        if align[query]=='#N/A':
                            align[query]=lines[1]
                for x in ID:
                    f.write('%s\t'%align[x])
                f.write('\n')
                            
                            
                            
def main():
    ID=get_ref_ID()
    find_blastp(ID)

if __name__ =='__main__':
    main()
