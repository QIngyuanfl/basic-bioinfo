#! /usr/bin/env/python
import sys
from Bio import SeqIO

def main(blastx, fa):
    bx = {}
    with open(blastx) as f:
        for line in f:
            line = line.rstrip().split()
            nucl = line[0]
            prot = line[1]
            nucl_start, nucl_end = int(line[6]), int(line[7])
            prot_start, prot_end = int(line[8]), int(line[9])
            if nucl not in bx:
                bx[nucl] = {}
            if prot not in bx[nucl]:
                bx[nucl][prot] = []
            bx[nucl][prot].append([nucl_start, nucl_end, prot_start, prot_end])

    seq = {}
    for seq_record in SeqIO.parse(fa, 'fasta'):
        seq[seq_record.id] = seq_record

    extract = {}
    for nucl in bx:
        for prot in bx[nucl]:
            pos = sorted(bx[nucl][prot], key = lambda x:x[2])
            prot_loci = []
            nucl_loci = []
            for i in pos:
                if i[0] < i[1]:
                    reverse = False
                if i[0] > i[1]:
                    reverse = True
                if len(i) == 1:
                    if reverse:
                        extract[prot] = seq[nucl][i[1] - 1:i[0]].reverse_completement()
                    if not reverse:
                        extract[prot] = seq[nucl][i[0] - 1:i[1]]
                if len(i) > 1:
                    nucl_loci.append(i[0])
                    nucl_loci.append(i[1])
                    prot_loci.append(i[2])
                    prot_loci.append(i[3])
            if len(pos) > 1:
                for j in range(len(prot_loci)):
                    print prot_loci[j]
                        

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2])
