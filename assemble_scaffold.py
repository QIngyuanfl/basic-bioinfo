import sys
from Bio import SeqIO
def read_list():
    info=[]
    with open(sys.argv[1]) as f:
        for line in f:
            line=line.strip().split()
            if len(line)==3:
                info.append([line[0],line[1],line[2]])
            if len(line)==2:
                info.append([line[0],line[1],line[2]])
    return info

def read_scaffolds():
    seq={}
    with open(sys.argv[2]) as f:
        for seq_record in SeqIO.parse(f,'fasta'):
            seq[seq_record.id]=seq_record.seq
    return seq

def write_fasta(info,seq):
    name=sys.argv[3].split('.')[0]
    with open(sys.argv[3],'w') as f:
        f.write('>%s\n'% name)
        for i in info:
            for key,value in seq.items():
                if i[0] in key:
                    if i[1] =='+':
                        f.write('%s' % value)
                    if i[1] =='-':
                        f.write('%s' % value.reverse_complement())
                    if len(i)==3:
                        f.write('N'*int(i[2]))
                    if len(i)==2:
                        f.write('N'*20)
def main():
    info=read_list()
    seq=read_scaffolds()
    write_fasta(info,seq)
if __name__ == '__main__':
    main()
