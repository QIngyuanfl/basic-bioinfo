#!/usr/bin/env  python
import sys
from Bio import SeqIO
def read_blastn():
    with open (sys.argv[1]) as f:
        loci={}
        info=f.readlines()
        length= len(info)
        d={}
        init=info[0].strip().split()
        if 'RNA' in init[0]:
            if 'tRNA' in init[0]:
                loci[int(init[8]),int(init[9])]=[init[0]]
            if 'rRNA' in init[0]:
                length=abs(int(init[8])-int(init[9]))
                distance=int(init[0].split(':')[-1])
                if abs(length-distance) <= 100:
                    loci[int(init[8]),int(init[9])]=[line[0]]
        if 'cds' in init[0]:
            if init[6] == '1' :
                d[init[0]] = init[7]
                loci[int(init[8]),int(init[9])]=[init[0]]
        for i in range(1,length):
            line=info[i].strip().split()
            if 'RNA' in line[0]:
                if 'rRNA' in line[0]:
                    length=abs(int(line[8])-int(line[9]))
                    distance=int(line[0].split(':')[-1])
                    if abs(length-distance) <= 100:
                        loci[int(line[8]),int(line[9])]=[line[0]]
                if 'tRNA' in line[0]:
                    loci[int(line[8]),int(line[9])]=[line[0]]
            if 'cds' in line[0]:
                if line[0] == info[i-1].strip().split()[0]:
                    info[i]==0
                if line[6] == '1' and info[i] !=0:
                    d[line[0]] = line[7]
                    loci[int(line[8]),int(line[9])]=[line[0]]
        return loci,d

def write_tbl(loci,d):
    with open(sys.argv[2]) as f:
        for lines in f:
            if lines[0] == '>':
                line=lines.lstrip('>').strip().split('\t')
                for i,j in loci.items():
                    if line[0] ==j[0]: 
                        if len(line)==5:
                            j.append([line[1],line[2],line[3],line[4]])
                        else:
                            j.append([line[1],line[2],line[3]])
    site=list(sorted(loci))
    with open(sys.argv[3],'w') as f:
        f.write('>Feature\n')
        Gene=0
        for i in site:
            for x,y in loci.items():
                gene=y[1][0].split('[gene=')[-1].split(']')[0]
                if 'product' in y[1][1]:
                    product=y[1][1].split('[product=')[-1].split(']')[0]
                    if '%3B' in product:
                        product=product.replace('%3B',';')
                    if '%2C' in product:
                        product=product.replace('%2C',',')
                if 'note' in y[1][1]:
                    note=y[1][1].split('[note=')[-1].split(']')[0]
                    if '%3B' in note:
                        note=note.replace('%3B',';')
                    if '%2C' in note:
                        note=note.replace('%2C',',')
                if i==x:
                    if y[1][2]=='type=CDS':
                        length=int(d[y[0]])
                        cmp_len= int(y[1][3])-length
                        if -3<=cmp_len<=3:
                            start=i[0]
                            stop=i[1]
                            if cmp_len >0:
                                if i[0]<i[1]:
                                    stop=i[1] +cmp_len
                                if i[0] >i[1]:
                                    stop=i[1]-cmp_len
                            if cmp_len <0:
                                print y[0]
                            if 'note' in y[1][1]:
                                note=y[1][1].split('[note=')[-1].split(';')[0]
                                product=y[1][1].split('product=')[-1].split(']')[0]
                                if 'Gene' not in gene and product !='hypothetical protein':
                                    f.write('%s\t%s\tgene\n\t\t\tgene\t%s\n%s\t%s\tCDS\n\t\t\tnote\tpseudo\n\t\t\tnote\t%s\n\t\t\tproduct\t%s\n\t\t\ttransl_table\t11\n'%(start,stop,gene,start,stop,note,product))
                                if 'Gene' in gene and product !='hypothetical protein':
                                    Gene+=1
                                    f.write('%s\t%s\tgene\n\t\t\tgene\tUniGene%s\n%s\t%s\tCDS\n\t\t\tnote\tpseudo\n\t\t\tnote\t%s\n\t\t\tproduct\t%s\n\t\t\ttransl_table\t11\n'%(start,stop,Gene,start,stop,note,product))
                            else:
                                if 'Gene' not in gene and product !='hypothetical protein':
                                    f.write('%s\t%s\tgene\n\t\t\tgene\t%s\n%s\t%s\tCDS\n\t\t\tproduct\t%s\n\t\t\ttransl_table\t11\n'%(start,stop,gene,start,stop,product))
                                if 'Gene' in gene and product !='hypothetical protein':
                                    Gene+=1
                                    f.write('%s\t%s\tgene\n\t\t\tgene\tUniGene%s\n%s\t%s\tCDS\n\t\t\tproduct\t%s\n\t\t\ttransl_table\t11\n'%(start,stop,Gene,start,stop,product)) 
                    if y[1][2]=='type=tRNA':
                        f.write('%s\t%s\tgene\n\t\t\tgene\t%s\n%s\t%s\ttRNA\n\t\t\tproduct\t%s\n'%(i[0],i[1],gene,i[0],i[1],product))
                    if y[1][2]=='type=rRNA':
                        if 'Gene' not in gene:
                            f.write('%s\t%s\tgene\n\t\t\tgene\t%s\n%s\t%s\trRNA\n\t\t\tproduct\t%s\n'%(i[0],i[1],gene,i[0],i[1],product))
                        if 'Gene' in gene:
                            Gene+=1
                            f.write('%s\t%s\tgene\n\t\t\tgene\tUniGene%s\n%s\t%s\trRNA\n\t\t\tproduct\t%s\n'%(i[0],i[1],Gene,i[0],i[1],product))

    
        
def main():
    if len(sys.argv) ==1:
        print 'Usage: split.py <blastn> <ref_fasta> <out_tbl>'
    else:
        loci,d=read_blastn()
        write_tbl(loci,d)

if __name__ =='__main__':
    main()
