#!/usr/bin/env python

import os
import re
import sys
import commands

sys.path.append("/hellogene/scgene01/user/chenjiehu/bin/recovery/module/")
import pipeline

def parseSortBam(sam, split):
    name = sam.split('.')[0]
    output = commands.getoutput('samtools view -H %s >header' % sam)
    h = open('header')
    header = h.read()
    h.close()
    header_size = os.path.getsize('header')
    sam_size = os.path.getsize(sam)
    reads_size = 460
    reads_num = (sam_size - header_size)/reads_size
    base_num = reads_num*150
    reads_per_sam = reads_num/split
    if reads_per_sam % 2 != 0:
        reads_per_sam += 1
    count = 0
    n = 1
    if not os.path.exists('sam'):
        os.mkdir('sam')
    init = open('sam/%s1.sam' % name, 'w')
    init.write(header)
    with open(sam) as f:
        for line in f:
            if line[0] == '@':
                continue
            count += 1
            if n == 1:
                init.write(line)
                if count > reads_per_sam:
                    n += 1
                    init.close()
                    
            if n > 1:
                if count <= reads_per_sam:
                    out.write(line)
                if count > reads_per_sam:
                    count = 0
                    try:
                        out.close()
                    except:
                        pass
                    out = open('sam/%s%d.sam' % (name, n), 'w')
                    n += 1
                    out.write(header)
    return base_num

def sort(sam, split):
    name = sam.split('.')[0]
    shell = []
    if not os.path.exists('sort'):
        os.mkdir('sort')
    for i in xrange(1, int(split)+1):
        with open('sort/sort_%s%d.sh' % (name, i), 'w') as sh:
            sh.write('samtools sort %s%d.sam -o sort/%s%d.sort.bam -@ 4' % (name, i, name, i))
        shell.append('sort/sort_%s%d.sh' % (name, i))
    pipeline.man_jobs(shell, '-cwd -l vf=8g,p=4')

def depth():
    sort_dir = os.listdir('sort')
    bam = []
    with open ('sort/bam_list', 'w') as f:
        for i in sort_dir:
            if i.endswith('bam'):
                f.write('sort/%s\n' % i)
    os.system('samtools depth -f sort/bam_list -a > all.dp')

def covVSdepth(gc): # base_num as a parameter
    length500 = {}
    with open(gc) as f:
        for line in f:
            if line[0] == '#':
                continue
            line = line.rstrip().split()
            length = int(line[-2])
            if length <= 500:
                continue
            ID = line[0].split('>')[-1].split('_')[0]
            length500[ID] = 0
    names = {}
    first = True
    genome_len = 0
    with open('all.dp') as f:
        for line in f:
            
            lines = line.rstrip().split()[2:]
            ID = line.rstrip().split()[0]
            if ID not in length500:
                continue
            genome_len += 1
            counter = 0
            for i in xrange(len(lines)):
                if first:
                    names['dp%s' %(i+1)] = 0
                    if lines[i] == '0':
                        pass
                    if lines[i] != '0':
                        counter = 1
                    if counter == 1:
                        names['dp%s' %(i+1)] += 1
                if not first:
                    if lines[i] == '0':
                        pass
                    if lines[i] != '0':
                        counter = 1
                    if counter == 1:
                        names['dp%s' % (i+1)] += 1
            first = False

    with open('covVSdepth.txt', 'w') as f:
        f.write('Depth\tCoverage\n')
        for i in sorted(names, key = lambda x:int(x.split('dp')[-1])):
            f.write('%s\t%s\n' % (i.split('dp')[-1], float(names[i])/(genome_len)))

def main():
    parseSortBam(sys.argv[1], int(sys.argv[2]))
    sort(sys.argv[1], sys.argv[2])
    depth()
    covVSdepth(sys.argv[3])


if __name__ == "__main__":
    main()
