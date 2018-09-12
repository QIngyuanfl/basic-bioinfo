#! /usr/bin/env python
import re
import os
import sys
import time
import threading
import commands
import pickle
from Bio import SeqIO
from argparse import ArgumentParser

sys.path.append("/hellogene/scgene01/user/chenjiehu/bin/recovery/module/")
import pipeline

__author__ = 'Qingyuan Zhang(zhangqingyuan@scgene.com)'
__version__ = 'V1.1'
__date__ = 'August 22nd, 2018'

FASTQC = '/hellogene/scgene01/bio/software/FastQC/fastqc'
ITOOLS = '/hellogene/scgene01/bio/bin/iTools'

def read_params(argvs):
    parser = ArgumentParser(description = 'RAD pipline' +__version__+ '(' +__date__ + ')'
                                                'AUTHOR:'+__author__)
    arg = parser.add_argument
    arg('-t', choices =['x','dir'], type = str, default = None, metavar = 'task')
    arg('-pop', type = str, default = 'Unknown', metavar = 'The population name')
    arg('-ploidy', type = str, default = '2', metavar = 'The ploidy of samples')
    return vars(parser.parse_args())

def creat_workspace():
    if not os.path.exists('RAD'):
        os.mkdir('RAD')
        os.system('mkdir RAD/01assemble RAD/02join2ref RAD/03fqIndex RAD/04align RAD/05sampe RAD/06division RAD/07sort RAD/08ECSelect RAD/09Haplotype RAD/10Misa RAD/11SNP RAD/12SSR ')
        os.system('cp /hellogene/scgene01/bio/software/misa/misa* RAD/10Misa/')
        
def extract_fq_from_dir():
    file_name = os.listdir('RAD/03fqIndex')
    fq={}
    for f in file_name:
        if f.endswith('1.fq'):
            fq[f] = f.split('1.fq')[0] + ('2.fq')
    return fq

def getQ20(fq):
    if not os.path.exists('RAD/qc'):
        os.mkdir('RAD/qc')
    cmd = ITOOLS + ' Fqtools stat -MinBaseQ ! -CPU 4 ' 
    for f in sorted(fq):
        cmd += ' -InFq' + ' ' + 'RAD/03fqIndex/' + f + ' -InFq ' + ' ' + 'RAD/03fqIndex/' + fq[f]
    cmd += ' -OutStat RAD/qc/baseinfo.txt'
    if os.path.exists(ITOOLS):
        os.system(cmd)
    else:
        print '-- Cannot find iTools in path "/hellogene/scgene01/bio/bin/bwa"'
        exit(1)

def fastqc(fq):
    cmd = FASTQC + ' -o RAD/qc/fastqc --extract -t 16 -f fastq'
    if not os.path.exists('RAD/qc/fastqc'):
        os.mkdir('RAD/qc/fastqc')
    for f in sorted(fq):
        cmd += (' RAD/03fqIndex/' + f + ' ' + 'RAD/03fqIndex/' + fq[f])
    with open('fastqc.sh', 'w') as sh:
        sh.write(cmd)
    pipeline.man_jobs(['fastqc.sh'], '-cwd -l vf=8g,p=16')
    exit(1)

def statRepli(fq):
    START = 30
    END = 60
    LINE = 4*600000
    outfile = open('RAD/qc/Repli.txt', 'w')
    for f in sorted(fq):
        seledic = {}
        fq1 = open('RAD/03fqIndex/%s' % ( f))
        fq2 = open('RAD/03fqIndex/%s' % (fq[f]))
        for i in xrange(1,LINE):
            if (i+2) % 4 == 0:
                read1 = fq1.readline()
                read2 = fq2.readline()
                seq1 = read1[START: END]
                seq2 = read2[START: END]
                seq = "%s%s" % (seq1, seq2)
                if seq in seledic:
                    seledic[seq] += 1
                else:
                    seledic[seq] = 1
            else:
                fq1.readline()
                fq2.readline()
        rep = 0
        for i in seledic.values():
            if i!= 1:
                rep += (i-1)
        reprate = round(float(rep)*4*100/LINE, 2)
        outfile.write('%s\tRepetition Rate: %s%%\n' %(f.split('1.fq')[0],reprate))
    outfile.close()

def info2xls1(baseinfo, replic):
    outfile = open('RAD/qc/Basicinfo.xls', 'w')
    outfile.write('#%s \t%s \t%s \t%s \t%s \t%s \t%s' %('FastqFile',
        'ReadsNum', 'BaseNum', 'GC', 'Q20', 'Q30', 'RepRate') + '\n')
    infolist = {}
    repline = open(replic, 'r')
    for line in repline:
        sample_name = line.strip().split()[0]
        reprate = line.strip().split()[-1]
        infolist[sample_name] = [reprate]
    baseinfo = open(baseinfo, 'r')
    for i in baseinfo:
        if i.startswith('##'):
            samplename = i.strip().strip('#')
            infolist[re.split(r"[12].\w*fq", samplename)[0]].append(samplename)
        if i.startswith('#ReadNum'):
            numlist = re.findall(r'\d+', i)
            infolist[re.split(r"[12].\w*fq", samplename)[0]].append(numlist[0])
            infolist[re.split(r"[12].\w*fq", samplename)[0]].append(numlist[1])
        if i.startswith('#GC'):
            gc = i.split()[1]
            infolist[re.split(r"[12].\w*fq", samplename)[0]].append(gc)
        if i.startswith('#BaseQ') and re.search(r'Q20', i):
            Q20 = i.split()[-1]
            infolist[re.split(r"[12].\w*fq", samplename)[0]].append(Q20)
        if i.startswith('#BaseQ') and re.search(r'Q30', i):
            Q30 = i.split()[-1]
            infolist[re.split(r"[12].\w*fq", samplename)[0]].append(Q30)
    for x in sorted(infolist):
        outfile.write('\t'.join(infolist[x][1:7])+'\n')
        outfile.write('\t'.join(infolist[x][7:]) + '\t' + infolist[x][0]+'\n')
    outfile.close()

def join_survey_scaffolds():
    os.system('python /hellogene/scgene02/RD/RAD/RADPlus/bin/joinRef.py RAD/01assemble/* > RAD/02join2ref/genome.fa')
    os.system('bwa index RAD/02join2ref/genome.fa')
    os.system('java -jar /hellogene/scgene02/RD/resequencing/bin/picard-tools-1.119/CreateSequenceDictionary.jar R=RAD/02join2ref/genome.fa O=RAD/02join2ref/genome.dict')

def ala_align(fq):
    jobs = []
    for f in fq:
        fq_name = f.split('1.fq')[0]
        sample_name = f.split('_')[0]
        with open('RAD/04align/aln_%s.sh' % fq_name, 'w') as sh:
            sh.write('bwa aln -t 16 -f RAD/04align/%s.sai RAD/02join2ref/genome.fa RAD/03fqIndex/%s\n'%(f,f))
            sh.write('bwa aln -t 16 -f RAD/04align/%s.sai RAD/02join2ref/genome.fa RAD/03fqIndex/%s\n'%(fq[f], fq[f]))
            sh.write("bwa sampe -r '@RG\\tID:%s\\tSM:%s' -f RAD/05sampe/%s.sam RAD/02join2ref/genome.fa RAD/04align/%s.sai RAD/04align/%s.sai RAD/03fqIndex/%s RAD/03fqIndex/%s\n" % (sample_name, sample_name, fq_name, f, fq[f], f, fq[f]))
            sh.write('sleep 1s')
        jobs.append('RAD/04align/aln_%s.sh' % fq_name)
    pipeline.man_jobs(jobs, '-cwd -l vf=8g,p=16')

def cut_filter_sort():
    os.system('python /hellogene/scgene02/RD/RAD/RADPlus/RADSeq/splitSamPE.py RAD/05sampe/ RAD/06division/')
    os.system('python /hellogene/scgene01/user/chenjiehu/bin/recovery/resequencing/resequencing.py sort RAD/06division/ RAD/07sort/')
    exist = False
    if os.path.exists('RAD/08ECSelect/EClib.group') and os.path.exists('RAD/08ECSelect/EClib.list'):
        exist = True
    else:
        print 'Suspend! Please check 08ECSelect directory and replenish the configuration file.'
        while not exist:
            time.sleep(60)
            if os.path.exists('RAD/08ECSelect/EClib.group') and os.path.exists('RAD/08ECSelect/EClib.list'):
                exist = True
    os.system('python /hellogene/scgene02/RD/RAD/RADPlus/RADSeq/ECSelect2.py RAD/08ECSelect/EClib.list RAD/08ECSelect/EClib.group RAD/02join2ref/genome.fa RAD/08ECSelect/')
    with open('RAD/08ECSelect/EClib.group') as f0:
        values = f0.read().strip().split()
        enzyme = values[0]
        n1 = values[1]
        n2 = values[2]
    with open('RAD/08ECSelect/True-MseI.fa', 'w') as f1:
        for seq_record in SeqIO.parse('RAD/08ECSelect/%s-%s-%s.fa' % (enzyme, n1, n2), 'fasta'):
            if not 'N' in seq_record.seq:
                f1.write('>%s\n%s\n' % (seq_record.id, seq_record.seq))
        
    os.system('python /hellogene/scgene02/RD/RAD/RADPlus/RADSeq/targetIntervals.py RAD/08ECSelect/True-%s.fa > RAD/09Haplotype/target.list' % enzyme)

def HaplotypeCaller_GATK():
    bam = os.listdir('RAD/07sort')
    commands.getstatusoutput('samtools faidx RAD/02join2ref/genome.fa')
    index = []
    GATK = []
    cmd = 'java -Xmx10G -jar /hellogene/scgene02/RD/resequencing/bin/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T HaplotypeCaller -R RAD/02join2ref/genome.fa'
    gvcfs= 'java -Xmx10G -jar /hellogene/scgene02/RD/resequencing/bin/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R RAD/02join2ref/genome.fa'
    for fa in bam:
        if fa.endswith('_map.sort.bam'):
            sample_name = fa.split('_map.sort.bam')[0]
            with open('RAD/07sort/index_%s.sh' % fa,'w') as sh:
                sh.write('samtools index RAD/07sort/%s' % fa)
            with open('RAD/09Haplotype/GATK_%s.sh' % fa, 'w') as out:
                out.write(cmd + ' -I RAD/07sort/%s -o RAD/09Haplotype/%s.g.vcf --emitRefConfidence GVCF -L RAD/09Haplotype/target.list -nct 4' % (fa,sample_name))
            index.append('RAD/07sort/index_%s.sh' % fa)
            GATK.append('RAD/09Haplotype/GATK_%s.sh' % fa)
            gvcfs += ' -V RAD/09Haplotype/%s.g.vcf' %(sample_name)
    #pipeline.man_jobs(index, '-cwd -l vf=4g,p=1')
    #pipeline.man_jobs(GATK, '-cwd -l vf=4g,p=4')
    gvcfs += ' -o RAD/09Haplotype/All.vcf'
    with open('RAD/09Haplotype/combine.sh', 'w') as f:
        f.write('%s\n' % gvcfs)
    with open('RAD/09Haplotype/selectVar.sh', 'w') as f:
        #f.write('java -Xmx10G -jar /hellogene/scgene02/RD/resequencing/bin/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R RAD/02join2ref/genome.fa -V RAD/09Haplotype/All.vcf -selectType SNP -o RAD/09Haplotype/raw_snps.vcf\n')
        f.write('java -Xmx10G -jar /hellogene/scgene02/RD/resequencing/bin/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T VariantFiltration -R RAD/02join2ref/genome.fa -V RAD/09Haplotype/raw_snps.vcf --filterExpression "QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o RAD/09Haplotype/filtered_snps.vcf')
    #pipeline.man_jobs(['RAD/09Haplotype/combine.sh'], '-cwd -l vf=10g,p=1')
    pipeline.man_jobs(['RAD/09Haplotype/selectVar.sh'], '-cwd -l vf=10g,p=1')

def misa():
    os.system('ln -s ../02join2ref/genome.fa RAD/10Misa/genome.fa')
    os.chdir('RAD/10Misa/')
    os.system('perl misa.pl genome.fa')
    os.chdir('../../')

def SNPpattern():
    repeat_dic = {}
    with open('RAD/10Misa/genome.fa.misa') as f:
        for line in f:
            if re.search("^ID", line):
                continue
            lines = line.rstrip().split()
            chrom = lines[0]
            repeat_type = lines[2]
            if not re.search("^p", repeat_type):
                continue
            repeat_ele = lines[3].split(')')[-1]
            start = int(lines[-2])
            end = int(lines[-1])
            if chrom not in repeat_dic:
                repeat_dic[chrom] = {}
            for i in xrange(start, end + 1):
                repeat_dic[chrom][i] = repeat_ele
    
    scale = {}
    with open('RAD/09Haplotype/target.list') as f:
        for line in f:
            line = line.rstrip().split(":")
            chrom = line[0]
            start, end = line[1].split("-")
            if chrom not in scale:
                scale[chrom] = {}
            if (start, end) not in scale[chrom]:
                scale[chrom][start, end] = 0 
    
    SNP = {}
    with open('RAD/09Haplotype/filtered_snps.vcf') as f:
        for line in f:
            if re.search(r'^##', line):
                continue
            if line[0] == '#':
                sample_name = line.strip().split()[9:]
                continue
            line = line.strip().split()
            ref = line[3]
            alt = line[4]
            if len(ref) > 1 or len(alt) >1:
                continue
            Filter = line[6]
            if Filter != 'PASS':
                continue
            sample = line[9:]
            F = 0
            for i in sample:
                if i[0] == '.':
                    continue
                DP = i.split(':')[2]
                if int(DP) < 4:
                    F += 1
            if float(F)/len(sample) > 0.2:
                continue
            chrom = line[0]
            pos = int(line[1])
            if pos in repeat_dic[chrom]:
                continue
            
            target = 0
            for m in sorted(scale[chrom]):
                if pos > m[1]:
                    break
                if pos >= m[0] and pos <= m[1]:
                    scale[chrom][m] += 1
                if scale[chrom][m] > 1:
                    target = 1

            if target == 1:
               continue

            if chrom not in SNP:
                SNP[chrom] = {}
            if pos not in SNP[chrom]:
                SNP[chrom][pos] = []
            SNP[chrom][pos].append([ref, alt, sample])
    
    with open('RAD/11SNP/SNP.xls', 'w') as f:
        f.write('#Genome\tPos\tSNP_Type\tSNP\tSample_Info\n')
        for i in sorted(SNP, key = lambda x:int(x.split('genome')[-1])):
            for j in sorted(SNP[i]):
                print SNP[i][j]
                ref = SNP[i][j][0][0]
                alt = SNP[i][j][0][1]
                SNP_Type = SNPtype(ref, alt)
                snp = ref + '/' + alt
                f.write('%s\t%s\t%s\t%s\t' % (i, j, SNP_Type, snp))
                n = 0
                for y in SNP[i][j][0][2]:
                    n += 1
                    if SNP[i][j][0][2][0] == '.':
                        f.write('%s:NA\t' % sample_name[n])
                    if SNP[i][j][0][2][:2] == '0/0':
                        f.write('%s:%s\t' % (sample_name[n], ref*2))
                    if SNP[i][j][0][2][:2] == '1/1':
                        f.write('%s:%s\t' % (sample_name[n], alt*2))
                    if SNP[i][j][0][2][:2] == '0/1':
                        f.write('%s:%s\t' % (sample_name[n], ''.join(sorted(ref + alt))))
                f.write('%s\n')
                
                        

def SSR():
    sort_bam = os.listdir('RAD/07sort')
    with open('RAD/07sort/sort_bam.list', 'w') as f:
        for b in sort_bam:
            if b.endswith('_map.sort.bam'):
                f.write('RAD/07sort/%s\n' % b)
    with open('RAD/07sort/merge_bam.sh', 'w') as f:
        f.write('samtools merge -c -p -@ 32 -b RAD/07sort/sort_bam.list RAD/07sort/All_map.sort.bam\n')
        f.write('samtools index -@ 32 RAD/07sort/All_map.sort.bam')
    pipeline.man_jobs(['RAD/07sort/merge_bam.sh'], '-cwd -l vf=8g,p=32')
    os.system('python /hellogene/scgene02/RD/RAD/RADPlus/RADSeq/SSRFindPerSample.py RAD/02join2ref/genome.fa RAD/10Misa/genome.fa.misa RAD/07sort/All_map.sort.bam 2 > RAD/12SSR/All.SSR.txt')
    os.system('python /hellogene/scgene02/RD/RAD/RADPlus/RADSeq/PrimerDesign.py RAD/12SSR/All.SSR.txt RAD/12SSR/Primer')

def analyse_SNP(ploidy):
    table = {}
    SNP_info = {}
    SNP_file = os.listdir('RAD/11SNP')
    sample_name = []
    location = []
    loc_dic = {}
    for doc in SNP_file:
        if doc.endswith('.xls'):
            sample_name.append(doc.split('.xls')[0])
            location.append(re.split(r'[0-9]', doc)[0])
    location = list(set(sorted(location)))
    count = 0
    for point in location:
        count += 1
        loc_dic[point] = str(count)
    with open('RAD/11SNP/Primer/SNP.xls') as f:
        f.next()
        for line in f:
            line = line.strip().split()
            SNP_ID = line[0]
            SNP = line[2]
            n_Sample = int(line[3])
            Sample_info = line[4]
            if n_Sample >= 0.8*(len(sample_name)):
                hetero = 0
                SNP_info[SNP_ID] = Sample_info, SNP
                b1 = SNP[0]; b2 = SNP[2]
                pb1 = 0; pb2 = 0
                for i in re.findall('[AGCT]/[AGCT]', Sample_info):
                    pair = i.split('/')
                    if pair[0] != pair[1]:
                        hetero += 1
                    pb1 += i.count(b1)
                    pb2 += i.count(b2)
                pb1 = float(pb1)/ (2*len(re.findall('[AGCT]/[AGCT]', Sample_info)))
                pb2 = float(pb2)/ (2*len(re.findall('[AGCT]/[AGCT]', Sample_info)))
                observed_heterozygosity = float(hetero)/n_Sample
                expected_heterozygosity = 1-(observed_heterozygosity*observed_heterozygosity+(1-observed_heterozygosity)*(1-observed_heterozygosity))
                PIC = 1 - (pb1*pb1+pb2*pb2) - 2*(pb1*pb1*pb2*pb2)
                table[SNP_ID] = [str(observed_heterozygosity), str(expected_heterozygosity), str(PIC)]
    base_dic = {"A":'11', "C":'12', "G":'13', "T":'14'}
    with open('RAD/11SNP/SNP_Structure.txt', 'w') as s:
        s.write('\t\t%s\n' % '\t'.join(sorted(SNP_info)))
        for ind in sample_name:
            loc = re.split(r'[0-9]', ind)[0]
            mark1 = []; mark2 = []
            for x in sorted(SNP_info):
                if ind in SNP_info[x][0]:
                    allele = SNP_info[x][0].split(ind)[-1].split(';')[0].split(':')[-1]
                    mark1.append(base_dic[allele[0]])
                    mark2.append(base_dic[allele[2]])
                if ind not in SNP_info[x][0]:
                    mark1.append('-9')
                    mark2.append('-9')
            s.write('%s\t%s\t%s\n' % (ind, loc_dic[loc], '\t'.join(mark1)))
            s.write('%s\t%s\t%s\n' % (ind, loc_dic[loc], '\t'.join(mark2)))
                    
                     

            
    with open('RAD/11SNP/SNP.csv', 'w') as w:
        w.write('Sample_ID\tLocation\t%s\n' % '\t'.join([str(x) for x in xrange(1, len(SNP_info)+1)]))
        for sample in sorted(sample_name):
            location = re.split(r'[0-9]', sample)[0]
            w.write('%s\t%s\t' % (sample, location))
            for j in sorted(SNP_info):
                if sample in SNP_info[j]:
                    allele = sorted(SNP_info[j].split(sample)[-1].split(';')[0].split(':')[-1].split('/'))
                    allele = allele[0] + allele[1]
                if sample not in SNP_info[j]:
                    allele = 'NA'
                w.write('%s\t' % allele)
            w.write('\n')
    '''
    with open('RAD/11SNP/SNP.R', 'w') as r:
        r.write('library("poppr")\n')
        r.write('mydata <- read.table("SNP.csv", header = TRUE, check.names = FALSE)\n\
        dim(mydata)\n\
        ind <- mydata$Sample_ID\n\
        dim(mydata)\n\
        locus <- mydata[, -c(1, %s:ncol(mydata))]\n' % len(SNP_info))
        r.write('data <- df2genind(locus, ploidy = %s, ind.names = ind, sep = "\\t")\n' % ploidy)
        r.write('dt <- missingno(data, type = "loci", cutoff = 0.05, quiet = FALSE)\n')
        r.write('pop <- summary(dt)\n')
    #os.system('Rscript RAD/11SNP/SNP.R > RAD/11SNP/expected_heterozygosity.txt')
    with open('RAD/11SNP/expected_heterozygosity.txt') as e:
        for line in e:
            if line.startswith('// Expected heterozygosity:'):
                line = line.strip().split(':')[-1].split()
                print len(line )
    '''

    with open("RAD/11SNP/SNP_Analyse.txt", 'w') as o:
        o.write('SNP_ID\tobserved_heterozygosity\texpected_heterozygosity\tPIC\n')
        for tab in sorted(table):
            o.write('%s\t%s\n' % (tab, '\t'.join((table[tab]))))
    

def analyse_SSR(ploidy):
    pop = []
    sample = []
    SSR_info = {}
    pop = []
    pop_dic = {}
    Primer = os.listdir('RAD/12SSR/Primer')
    sample_file = os.listdir('RAD/11SNP/')
    for i in sample_file:
        if i.endswith('.xls'):
            sample.append(i.split('.')[0])
            pop.append(re.split(r'[0-9]', i)[0])
    pop = list(set(sorted(pop)))
    count = 0
    for i in pop:
        count += 1
        pop_dic[i] = count
    '''  
    for j in sample:
        pop.append(re.split(r'[0-9]', j)[0])
        for f in sorted(Primer):
            if re.match(r'Polymorphism_p[1-9].xls', f)  != None:
                ssr_file = re.match(r'Polymorphism_p[1-9].xls', f).group()
                with open('RAD/12SSR/Primer/%s' % ssr_file) as ssr:
                    for line in ssr:
                        if line[0] != '#':
                            line = line.strip().split()
                            l_sample = re.findall(r'[0-9]*/[0-9]*', line[6])
                            if len(l_sample) > 0.8*(len(sample)):
                                Sample_info = line[6].split(';')
                                for x in Sample_info:
                                    if j in x:
                                        if j not in SSR_info:
                                            SSR_info[j] = []
                                        SSR_info[j].append([line[0],re.split(r'[:/;?]',x)])
                                    if j not in ''.join(Sample_info):
                                        if j not in SSR_info:
                                            SSR_info[j] = []
                                        try:
                                            if SSR_info[j][-1] != [line[0], [j, '0', '0']]:
                                                SSR_info[j].append([line[0],[j,'0','0']])
                                        except IndexError:
                                            pass
    with open('RAD/12SSR/SSR.csv', 'w') as p:
        loci = []
        for x in sorted(SSR_info):
            for y in SSR_info[x]:
                loci.append(y[0].replace('_', '-'))
            break
        p.write('%s\t%s\t%s\n' % (len(loci), len(SSR_info), len(pop)))
        p.write('ind\tPop\t%s\n' % '\t'.join(loci*2))
        for m in sorted(SSR_info):
            location = re.split(r'[0-9]', m)[0]
            p.write('%s\t%s\t' % (m, location))
            for n in SSR_info[m]:
                p.write('%s\t%s\t' % (n[1][1], n[1][2]))
            p.write('\n')
    with open('RAD/12SSR/SSR.R', 'w') as r:
        r.write('library("poppr")\n')
        r.write('Mydata <- read.genalex("SSR.csv", ploidy = %s, geo = FALSE, region = FALSE, genclone = TRUE, sep = ",", recode = FALSE)\n' % ploidy)
        r.write('popdata <- summary(mydata)\n')
        r.write('pop <- missingno(popdata, type ="loci", cutoff = 0.05, quiet = FALSE)')
    '''

    table = {}
    for f in sorted(Primer):
        if re.match(r'Polymorphism_p[1-9].xls', f)  != None:
            ssr_file = re.match(r'Polymorphism_p[1-9].xls', f).group()
            with open('RAD/12SSR/Primer/%s' % ssr_file) as ssr:
                for line in ssr:
                    if line[0] != '#':
                        line = line.strip().split()
                        SSR_ID = line[0]
                        Sample_Info = line[6]
                        l_sample = re.findall(r'[0-9]*/[0-9]*', line[6])
                        if len(l_sample) > 0.8*(len(sample)):
                            SSR_info[SSR_ID] = Sample_Info
                            poly = []
                            frequency = []
                            hetero = 0
                            for repeat in l_sample:
                                rr = repeat.split('/')
                                poly.append(rr[0])
                                poly.append(rr[1])
                                if rr[0] != rr[1]:
                                    hetero += 1
                            observed_heterozygosity = float(hetero)/len(l_sample)
                            
                            for ele in set(sorted(poly)):
                                n = poly.count(ele)    
                                frequency.append(float(n)/len(poly))
                            pic = PIC(len(set(sorted(poly))), frequency)
                            He = expected_heterozygosity(len(set(sorted(poly))), frequency)
                            table[SSR_ID] = [str(observed_heterozygosity), str(He), str(pic)]
    with open('RAD/12SSR/SSR_Analyse.txt', 'w') as o:
        o.write('SSR_ID\tobserved_heterozygosity\texpected_heterozygosity\tPIC\n')
        for tab in sorted(table):
            o.write('%s\t%s\n' %(tab, '\t'.join(table[tab])))
    with open('RAD/12SSR/SSR_Structure.txt', 'w') as s:
        s.write('\t\t%s\n' % '\t'.join(sorted(table)))
        for ind in (sample):
            pop = re.split(r'[0-9]', ind)[0]
            mark1 = []; mark2 = []
            for i in sorted(table):
                if ind in SSR_info[i]:
                    score = SSR_info[i].split(ind)[-1].split(';')[0].split(':')[-1]
                    score = score.split('/')
                    mark1.append(score[0])
                    mark2.append(score[1])
                if ind not in SSR_info[i]:
                    mark1.append('-9')
                    mark2.append('-9')
            s.write('%s\t%s\t%s\n' % (ind, pop_dic[pop], '\t'.join(mark1)))
            s.write('%s\t%s\t%s\n' % (ind, pop_dic[pop], '\t'.join(mark2)))
       
                    
def SNPtype(ref, alt):
    if len(ref) == 1 and len(alt) == 1 :
        co = ''.join(sorted(ref + alt))
        if co == 'AG' or co == 'CT':
            mutation = "Transitions"
        else:
            mutation = "Transversions"
    return mutation

def PIC(n, frequency):
    f2 = 0; f3 = 0
    for i in range(1,n+1):
        f2 += frequency[i-1]*frequency[i-1]
        try:
            f3 += 2*frequency[i-1]*frequency[i-1]*frequency[i]*frequency[i]
        except IndexError:
            pass
    PIC = 1 - f2 - f3
    return PIC

def expected_heterozygosity(n,frequency):
    f1 = 0
    for i in xrange(n):
      f1 += frequency[i]*frequency[i]
    He = 1-f1
    return He

def qc(fq):
    getQ20(fq)
    fastqc(fq)
    statRepli(fq)
    info2xls1('RAD/qc/baseinfo.txt', 'RAD/qc/Repli.txt')

def survey():
    join_survey_scaffolds()
    return 0

def RAD(fq):
    #ala_align(fq)
    #cut_filter_sort()

    #HaplotypeCaller_GATK()
    #misa()
    SNPpattern()
    #SSR()
    
def main():
    args = read_params(sys.argv)
    if args['t'] == 'dir':
        creat_workspace()
    if args['t'] == 'x':
        fq = extract_fq_from_dir()
        threads = []
        #t1 = threading.Thread(target = qc, args = (fq,))
        #threads.append(t1)
        #t2 = threading.Thread(target = survey)
        #threads.append(t2)
        #for t in threads:
        #    t.start()
        #for t in threads:
        #    t.join()
        RAD(fq)
        #analyse_SNP(args['ploidy'])
        analyse_SSR(args['ploidy'])

if __name__ == '__main__':
    main()
