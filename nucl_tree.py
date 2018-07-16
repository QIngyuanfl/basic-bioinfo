#!/usr/bin/env python
import os
import sys
import numpy as np

def parse():
    f=open(sys.argv[1],'r')
    data=[]
    dt=[]
    for lines in f.readlines():
        line=lines.rstrip()
        data.append(line)
    for i in data:
        a=i.split()
        dt.append(a)
    column=np.array(dt)
    asn=[a[0] for a in column]
    f.close()
    return asn

def downloadsASN():
    ID=parse()
    if not os.path.exists('./cdsbank'):
		os.system('mkdir cdsbank')
    fs=os.listdir('./cdsbank')
    for i in range(0,len(ID)):
        if (ID[i]+'.txt') not in fs:
            os.system('asn2all -r -A %s -f d > ./cdsbank/%s.txt' % (ID[i],ID[i]))
        else:
            print ('%s.txt existed' % ID[i])
    for j in range (0,len(ID)):
        f=open('./cdsbank/%s.txt'%ID[j],'r')
        words=f.readline()
        if words =='':
            os.system('asn2all -r -A %s -f d > ./cdsbank/%s.txt' % (ID[j],ID[j]))
            print ('%s.txt has been downloaded' % ID[j])
        f.close()
    print ('Mission Complete')

def common_gene():
    fs=os.listdir('./cdsbank')
    text=[]
    msg=[]
    genes=[]
    common_gene=[]
    for i in fs:
        handle=open('./cdsbank/%s' %i)
        for lines in handle:
            if lines[0] == '>':
                line=lines.rstrip()
                text.append(line)
        handle.close()
    for info in text:
        a=info.split()
        msg.append(a)
    table=np.array(msg)
    genes=[b[1] for b in table]
    gene=set(genes)
    for x in gene:
        ngenes=genes.count(x)
        if ngenes >= len(fs):
            common_gene.append(x)
    return common_gene

def verify_gene():
    cmg=common_gene()
    fs=os.listdir('./cdsbank')
    vrg=[]
    vg=[]
    for i in fs:
        handle=open('./cdsbank/%s' %i)
        words=handle.read()
        vrg.append(words)
        handle.close()
    for cg in cmg:
        n=0
        for x in range (len(vrg)):
            if cg in vrg[x]:
                n+=1
        if n == len(fs):
            vg.append(cg)
    
    return vg
        


def combine_gene():
    fs=os.listdir('./cdsbank')
    if not os.path.exists('samegene'):
		os.system('mkdir samegene')
    cmg=verify_gene()
    for cds in cmg:
        a=cds.split('[gene=')[1]
        b=a.strip(']')
        f=open('./samegene/%s.txt' % b ,'w+')
        for i in fs:
            handle=open('./cdsbank/%s' % i)
            words=handle.read()
            record=words.split('>')
            handle.close()
            cg=[]
            for j in range (0,len(record)):
                if cds in record[j]:
                    c=record[j].split('|')
                    d=c[1].split('_cds_')
                    accession=(d[0])
                    d=record[j].split('\n',1)
                    sequence=d[1]
                    cg+=record[j]
                    if cds not in cg:
                        f.write('>')
                        f.write(accession+'\n')
                        f.writelines(sequence)
                        cg=''.join(cg)

def muscle():
    if not os.path.exists('muscle'):
		os.system('mkdir muscle')
    f=os.listdir('./muscle')
    fs=os.listdir('./samegene')
    for i in fs:
        if (i+'.nucl') not in f:
            os.system('muscle -in ./samegene/%s -out ./muscle/%s.nucl' % (i,i)) 
        else:
            print (i+'.nucl existed')
            continue

def combine_asn():
    fs=os.listdir('./muscle')
    asn=parse()
    f=open('combine.fa','w')
    for num in asn:
        f.write('>')
        f.write(num + '\n')
        for i in fs:
            handle=open('./muscle/%s' % i)
            words=handle.read()
            records=words.split('>')
            for j in range(len(records)):
                if num in records[j]:
                    seq=records[j].split('\n',1)
                    sequence=''.join(seq[1].split('\n'))
                    f.write('%s'%sequence)
            handle.close()
        f.write('\n')
    f.close()
                    
def main():
    downloadsASN()
    combine_gene()
    muscle()
    combine_asn()

if __name__=="__main__":
    main()
