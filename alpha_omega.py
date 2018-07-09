#! /usr/bin/env python
import sys
from Bio import SeqIO
def main():
    loci=[]
    gene=[]
    if sys.argv[3]==str(2):
        M=['ATT','ATA','ATC','ATG','GTG']
        star=['TAA','TAG','AGA','AGG']
    if sys.argv[3] ==str(4):
        M=['TTA','TTG','CTG','ATT','ATC','ATA','ATG','GTG']
        star=['TAA','TAG']
    if sys.argv[3] ==str(11):
        M=['TTG','GTG','ATT','ATC','ATA','ATG','CTG']
        star=['TAA','TAG','TGA']
    with open(sys.argv[1]) as f:
        n=3
        fa=open(sys.argv[2],'r')
        sec_record=SeqIO.read(fa,'fasta')
        title=sec_record.id
        sequence=sec_record.seq
        for line in f:
            if line.strip().startswith('gene'):
                line=line.strip().split()
                n_line=f.next()
                if 'CDS' in n_line:
                    gene.append(line[1])
                if '0'<=n_line[0]<='9':
                    n_line=n_line.strip().split()
                    if n_line[2] == 'CDS':
                        start=int(n_line[0])
                        end=int(n_line[1])
                        loci.append([start,end])
        for i in range (len(gene)):
            if loci[i][0] < loci[i][1]:
                seq=sequence[loci[i][0]-1:loci[i][1]]
            if loci[i][0] > loci[i][1]:
                seq=sequence[loci[i][1]-1:loci[i][0]].reverse_complement()
            length =len(seq)
            aa=[seq[n:n+3] for n in range(0,len(seq[:-4]),3)]
            if length % 3 !=0:
                print ('Warning: invalid length for gene %s' % gene[i])
            if seq[:3] not in M:
                print ('Warning: incorrect initiation codon %s for gene %s at %s' %(seq[:3],gene[i],loci[i][0]))
            if seq[-3:] not in star:
                print('Warning: incorrect stop codon %s for gene %s at %s' % (seq[-3:],gene[i],loci[i][1]))
            for y in aa:
                if y in [x for x in star]:
                    print('Warning: internal stop codon in gene %s'% (gene[i]))
                    break
                

        print ('Examine complete')
if __name__ =='__main__':
    if len(sys.argv) !=4:
        print ('-------------------------------------------------')
        print ('This tool aims to check the abnormal initiation codon and stop condon in a feature table')
        print ('Usage: python alpha_omega.py <tbl> <tbl_fasta> <transl_table>')
        print ('-------------------------------------------------')
    if len(sys.argv) ==4:
        main() 

