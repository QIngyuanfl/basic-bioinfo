#! /usr/bin/env python
import sys
def complement_reverse(codon):
        paired={'A':'T','T':'A','C':'G','G':'C'}
        anticodon=''
        for i in codon:
            anticodon+=paired[i]
        anticodon=anticodon[::-1]
        return anticodon
def main():
    gene=[]
    product=[]
    upper=[]
    lower=[]
    with open (sys.argv[1]) as f:
        for line in f:
            if 'a'<=line.strip()[0]<='z':
                line=line.strip().split()
                if line[1].startswith('trn'):
                    gene.append(line[1])
                if line[1].startswith('tRNA'):
                    product.append(line[1])
    tRNA=zip(gene,product)
    for i in tRNA:
        codon=i[0].split('-')[1]
        product=i[1].split('-')[1]
        anticodon=complement_reverse(codon)
        upper.append(i[0].split('-')[0]+'-'+product+'('+codon+')')
        lower.append(i[1].split('-')[0]+'-'+product+'('+anticodon+')')
    with open(sys.argv[2],'w') as f:
        infile=open(sys.argv[1],'r')
        n=0
        m=0
        for line in infile:
            if 'a'<=line.strip()[0]<='z':
                lines=line.strip().split()
                if 'trn' in lines[1]:
                    f.write('\t\t\tgene\t%s\n'%upper[n])
                    n+=1
                if 'tRNA' in lines[1]:
                    f.write('\t\t\tproduct\t%s\n'%lower[m])
                    m+=1
                if 'trn' not in lines[1] and 'tRNA' not in lines[1]:
                    f.write(line)
            else:f.write(line)
        infile.close()
if __name__ =='__main__':
    main()
