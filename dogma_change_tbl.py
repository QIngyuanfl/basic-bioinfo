#!/usr/bin/env python
import sys
ontology={'trnA-TGC':'tRNA-Ala','ycf68':'Ycf68','trnI-GAT':'tRNA-Gly','trnP-GGG':'tRNA-Gly','trnV-TAC':'tRNA-Val','psbG':'photosystem II reaction center protein G','trnL-TAA':'tRNA-Leu','trnfM-CAT':'tRNA-Met','trnG-TCC':'tRNA-Gly','trnK-TTT':'tRNA-Lys','trnQ-TTG': 'tRNA-Gln', 'clpP': 'ATP-dependent Clp protease proteolytic subunit', 'trnM-CAT': 'tRNA-Met-M', 'trnY-GTA': 'tRNA-Tyr', 'rrn4.5': '4.5S ribosomal RNA', 'trnH-GTG': 'tRNA-His', 'atpI': 'ATP synthase CF0 A subunit', 'atpH': 'ATP synthase CF0 C subunit', 'atpB': 'ATP synthase CF1 beta subunit', 'atpA': 'ATP synthase CF1 alpha subunit', 'atpF': 'ATP synthase CF0 B subunit', 'atpE': 'ATP synthase CF1 epsilon subunit', 'rbcL': 'ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit', 'petD': 'cytochrome b6/f complex subunit IV', 'petG': 'cytochrome b6/f complex subunit V', 'petA': 'cytochrome f', 'trnT-GGT': 'tRNA-Thr', 'petL': 'cytochrome b6/f complex subunit VI', 'petN': 'cytochrome b6/f complex subunit VIII', 'rpl20': 'ribosomal protein L20', 'trnL-CAA': 'tRNA-Leu', 'rpl22': 'ribosomal protein L22', 'psaI': 'photosystem I subunit VIII', 'psaJ': 'photosystem I subunit IX', 'rrn23': '23S ribosomal RNA', 'psaA': 'photosystem I P700 chlorophyll a apoprotein A1', 'psaB': 'photosystem I P700 chlorophyll a apoprotein A2', 'psaC': 'photosystem I subunit VII', 'trnP-TGG': 'tRNA-Pro', 'rpl36': 'ribosomal protein L36', 'rpl32': 'ribosomal protein L32', 'trnR-TCT': 'tRNA-Arg', 'rpl33': 'ribosomal protein L33', 'rpoB': 'RNA polymerase beta subunit', 'rpoA': 'RNA polymerase alpha subunit', 'rps12': 'ribosomal protein S12', 'rps11': 'ribosomal protein S11', 'rps16': 'ribosomal protein S16', 'rps15': 'ribosomal protein S15', 'rps14': 'ribosomal protein S14', 'rps19': 'ribosomal protein S19', 'rps18': 'ribosomal protein S18', 'rpl2': 'ribosomal protein L2', 'trnS-GCT': 'tRNA-Ser', 'trnL-TAG': 'tRNA-Leu', 'trnV-GAC': 'tRNA-Val', 'rrn5': '5S ribosomal RNA', 'trnG-GCC': 'tRNA-Gly', 'trnW-CCA': 'tRNA-Trp', 'trnT-TGT': 'tRNA-Thr', 'accD': 'acetyl-CoA carboxylase beta subunit', 'petB': 'cytochrome b6', 'trnN-GTT': 'tRNA-Asn', 'ccsA': 'cytochrome c biogenesis protein', 'trnS-TGA': 'tRNA-Ser', 'rpoC2': "RNA polymerase beta'' subunit", 'infA': 'translation initiation factor 1', 'ycf15': 'Ycf15', 'rpl23': 'ribosomal protein L23', 'cemA': 'envelope membrane protein', 'ndhB': 'NADH dehydrogenase subunit 2', 'ndhC': 'NADH dehydrogenase subunit 3', 'ndhA': 'NADH dehydrogenase subunit 1', 'ndhF': 'NADH dehydrogenase subunit 5', 'ndhG': 'NADH dehydrogenase subunit 6', 'ndhD': 'NADH dehydrogenase subunit 4', 'ndhE': 'NADH dehydrogenase subunit 4L', 'ndhJ': 'NADH dehydrogenase subunit J', 'ndhK': 'NADH dehydrogenase subunit K', 'ndhH': 'NADH dehydrogenase subunit 7', 'ndhI': 'NADH dehydrogenase subunit I', 'rps7': 'ribosomal protein S7', 'rps4': 'ribosomal protein S4', 'rps3': 'ribosomal protein S3', 'rps2': 'ribosomal protein S2', 'trnF-GAA': 'tRNA-Phe', 'rps8': 'ribosomal protein S8', 'psbI': 'photosystem II protein I', 'ycf4': 'photosystem I assembly protein Ycf4', 'ycf2': 'Ycf2', 'ycf3': 'photosystem I assembly protein Ycf3', 'ycf1': 'Ycf1', 'trnS-GGA': 'tRNA-Ser', 'trnC-GCA': 'tRNA-Cys', 'matK': 'maturase K', 'rpl14': 'ribosomal protein L14', 'rpoC1': "RNA polymerase beta' subunit", 'rpl16': 'ribosomal protein L16', 'trnD-GTC': 'tRNA-Asp', 'trnE-TTC': 'tRNA-Glu', 'psbE': 'photosystem II protein V', 'psbD': 'photosystem II protein D2', 'trnR-ACG': 'tRNA-Arg', 'psbF': 'photosystem II protein VI', 'psbA': 'photosystem II protein D1', 'psbC': 'photosystem II 44 kDa protein', 'psbB': 'photosystem II 47 kDa protein', 'psbM': 'photosystem II protein M', 'psbL': 'photosystem II protein L', 'psbN': 'photosystem II protein N', 'rrn16': '16S ribosomal RNA', 'psbH': 'photosystem II protein H', 'psbK': 'photosystem II protein K', 'psbJ': 'photosystem II protein J', 'psbT': 'photosystem II protein T', 'trnI-CAT': 'tRNA-Met-I', 'psbZ': 'photosystem II protein Z'}
def read_dogma():
    with open (sys.argv[1]) as f:
        dogma=[]
        for line in f:
            line=line.strip().split()
            if 'U' in line[2]:
                line[2]=line[2].replace('U','T')
            if line[3]=='-':
                line[0]=str(int(line[0])-3)
                dogma.append([line[1],line[0],line[2]])
            if line[3]=='+':
                line[1]=str(int(line[1])+3)
                dogma.append([line[0],line[1],line[2]])
    return dogma
def write_tbl(dogma):
    with open(sys.argv[2],'w') as f:
        f.write('>Features\n')
        for i in dogma:
            if 'trn' in i[2]:
                f.write('%s\t%s\tgene\n\t\t\tgene\t%s\n%s\t%s\t%s\n\t\t\tproduct\t%s\n' %(i[0],i[1],i[2],i[0],i[1],'tRNA',ontology[i[2]]))
            elif 'rrn' in i[2]:
                f.write('%s\t%s\tgene\n\t\t\tgene\t%s\n%s\t%s\t%s\n\t\t\tproduct\t%s\n' %(i[0],i[1],i[2],i[0],i[1],'rRNA',ontology[i[2]]))
            else:
                f.write('%s\t%s\tgene\n\t\t\tgene\t%s\n%s\t%s\t%s\n\t\t\tproduct\t%s\n\t\t\ttransl_table\t11\n' %(i[0],i[1],i[2],i[0],i[1],'CDS',ontology[i[2]]))
def main():
    dogma=read_dogma()
    write_tbl(dogma)

if __name__ =='__main__':
    main()
