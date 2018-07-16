# -*- coding: utf-8 -*-
'''*****************************************************************************		     
				   last edit: 2016.5.12
*****************************************************************************'''
from Bio import SeqIO
from Bio.SeqUtils import GC
import sys
import re
from xlwt import Workbook,XFStyle,Borders
import xlwt
(xx,gbf)=sys.argv

TRNA={"Ala":"A","Arg":"R","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","His":"H","Ile":"I","Gly":"G","Asn":"N","Leu":"L",
"Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Sec":"U"," ":""}
base={'A':'T','a':'T','T':'A','t':'A','G':'C','g':'C','C':'G','c':'G'}

class ReadGenbank:
        gbfile=''
        refile=None
        def __init__(self,gbfile):
                self.gbfile=gbfile
                self.refile=SeqIO.read(self.gbfile,'gb')
        
        def getStats(self):
                element=[]
                onelettercode=[]
                codetype=[]
                size=[]
                GCpercent=[]
                beginend=[]
                fromto=[]
                No_ofAA=[]
                initiation=[]
                termination=[]
                anticode=[]
                geneperiod=[]
                intervening=[]
                intergen=[]
                for feature in self.refile.features:
                        seq=''
                        gc=0
                        beginend=[]
                        figure=re.findall(r'\d+',str(feature.location))
                        symbol=re.findall(r'[\+\-]',str(feature.location))[0]
                        chaintype='H' if symbol=='+' else 'L'
                        
                        # 2015.12.11 add three lines below, for facing the problem if reverse translate, numbers in figure list disorganized
                        figure0=[int(i) for i in figure]
                        figure0.sort()
                        figure=figure0[:]                  
                        
                        for j in range(len(figure)/2):
                                beginend.append((int(figure[j*2])+1,int(figure[j*2+1])))
                                seq+=self.refile.seq[int(figure[j*2]):int(figure[j*2+1])]
                                        
                        gc=round(GC(seq),2)/100
                        
                        if feature.type=='gene':
                                geneperiod.append([int(figure[0]),int(figure[1])])
                                
                        if feature.type=='CDS':		
                                element.append(feature.qualifiers['gene'][0])
                                onelettercode.append(' ')
                                codetype.append(chaintype)
                                size.append(len(seq))
                                AAnum=len(seq)/3-1 if len(seq)%3==0 else len(seq)/3
                                No_ofAA.append(AAnum)
                                GCpercent.append(gc)
                                fromto.append(beginend)
                                initiation.append(self.extractCodon(seq,symbol)[0])
                                termination.append(self.extractCodon(seq,symbol)[1])
                                anticode.append(' ')

                        if feature.type=='tRNA':
                                element.append(feature.qualifiers['gene'][0][:8])
                                aa=re.findall(r'\w+',str(feature.qualifiers['gene'][0]))[1]
                                onelettercode.append(TRNA[aa])
                                codetype.append(chaintype)
                                size.append(len(seq))
                                No_ofAA.append(' ')
                                GCpercent.append(gc)
                                fromto.append(beginend)
                                initiation.append(' ')
                                termination.append(' ')
				# 2015.11.02 change feature.qualifiers['note'][0])>>>feature.qualifiers['gene'][0])
				genecode=re.findall(r'\w+',str(feature.qualifiers['gene'][0]))
				# 2015.10.23 change anticode.append(anti[-1])>>>anticode.append(genecode[-1])
                                # anti=re.findall(r'\w+',str(feature.qualifiers['note'][0]))
                                anticode.append(genecode[-1])
                                
                        if feature.type=='rRNA':
                                element.append(feature.qualifiers['gene'][0][:8])
                                onelettercode.append(' ')
                                codetype.append(chaintype)
                                size.append(len(seq))
                                No_ofAA.append(' ')
                                GCpercent.append(gc)
                                fromto.append(beginend)
                                initiation.append(' ')
                                termination.append(' ')
                                anticode.append(' ')

		# if intervening over the range of int, then throw a warn
                for i in range(1,len(geneperiod)):
			intervening.append(geneperiod[i][0]-geneperiod[i-1][1])
			if geneperiod[i][0]-geneperiod[i-1][1]<32767:
                                intergen.append(str(self.refile.seq[geneperiod[i-1][1]:geneperiod[i][0]]))                                
			else: intergen.append('too long, over the range of int')		
                intervening.append(len(self.refile.seq)-geneperiod[-1][1])
		if len(self.refile.seq)-geneperiod[-1][1]<32767:
	                intergen.append(str(self.refile.seq[geneperiod[-1][1]:]))
		else: intergen.append('too long, over the range of int')

                workbook=xlwt.Workbook()
                sheet1=workbook.add_sheet('Genome Base Content',cell_overwrite_ok=True)
                # 2015.10.23 exchange 'InterveningSequence','IntergenicNucleotides'>>>'IntergenicNucleotides','InterveningSequence'
                title=['Gene/element','OneLetterCode','From','To','Size','GC_Percent','Codetype','No. of AA',
                       'InferredInitiationCoden','InferredTermination','Anti-Coden','IntergenicNucleotides',
                       'InterveningSequence']
                title=['NO.','Gene/element', 'Strand', 'From', 'to', 'Size (bp)', 'GC_Percent',
                       'Amino Acids (aa)', 'Inferred Initiation Codon', 'Inferred Termination Codon',
                       'One Letter code', 'Anti-codon', 'Intergenic nucleotide*(bp)', 'Intervening sequence']
                fromcol=[]
                tocol=[]
                t=0
                for i in range(len(element)):
                        for j in range(len(fromto[i])):
                                fromcol.append(fromto[i][j][0])
                                tocol.append(fromto[i][j][1])

                                # 2015.10.22 (i+j)>>>t
                                if j!=0:
                                        element.insert(t,' ')
                                        onelettercode.insert(t,' ')
                                        size.insert(t,' ')
                                        GCpercent.insert(t,' ')
                                        codetype.insert(t,' ')
                                        No_ofAA.insert(t,' ')
                                        initiation.insert(t,' ')
                                        termination.insert(t,' ')
                                        anticode.insert(t,' ')
                                        intervening.insert(t,' ')
                                        intergen.insert(t,' ')                                        
                                t+=1
                # 2016.5.12 change order of col
                for i in range(t):
                        sheet1.write(i+1,1,element[i])
                        sheet1.write(i+1,10,onelettercode[i])
                        sheet1.write(i+1,3,fromcol[i])
                        sheet1.write(i+1,4,tocol[i])
                        sheet1.write(i+1,5,size[i])
                        sheet1.write(i+1,6,GCpercent[i])
                        sheet1.write(i+1,2,codetype[i])
                        sheet1.write(i+1,7,No_ofAA[i])
                        sheet1.write(i+1,8,initiation[i])
                        sheet1.write(i+1,9,termination[i])
                        sheet1.write(i+1,11,anticode[i])
                        sheet1.write(i+1,12,intervening[i])
                        sheet1.write(i+1,13,intergen[i])
                            
                for i in range(len(title)):
                        sheet1.write(0,i,title[i])                        
                workbook.save('Gene Information.xls')

        def extractCodon(self,seq='',symbol='-'):
                if symbol=='-':
                        sequence=''
                        inverseq=seq[::-1]
                        for i in inverseq:
                                if i in base:
                                        sequence+=base[i]
                                else:
                                        sequence+=i
                        start=str(sequence[:3])
                        end=str(sequence[-3:]) if len(sequence)%3==0 else str(sequence[-(len(sequence)%3):])
                else:
                        start=str(seq[:3])
                        end=str(seq[-3:]) if len(seq)%3==0 else str(seq[-(len(seq)%3):])
                return [start,end]

ReadGenbank(gbf).getStats()




