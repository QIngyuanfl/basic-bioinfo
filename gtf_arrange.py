#! /usr/bin/env python

import sys
import re
import os
import math
import numpy

def readgtf():
    files = os.listdir(sys.argv[1])
    sample = {}
    sample_log10 = {}
    for name in files:
        if name.endswith('.gtf'):
            gtf = open(name, 'r')
            name = name.split('.')[0]
            sample[name] = []
            sample_log10[name] = []
            for line in gtf:
                line = line.rstrip()
                if not re.match(r'#', line):
                    newline = line.split('\t')
                    attributes = newline[8].split(';')
                    if newline[2] == 'transcript':
                        fpkm = float(attributes[4].split()[1][1:-1])
                        if fpkm > 0:
                            sample[name].append(fpkm)
                            fpkm = math.log10(fpkm)
                            if fpkm >= -2:
                                sample_log10[name].append(fpkm)
            gtf.close()
    return sample, sample_log10

def Boxplot(sample_log10):
    boxplot = open('boxplot.txt', 'w')
    boxplot.write('Sample\tlog10FPKM\n')
    for sample_name in sample_log10:
        for fpkm in sample_log10[sample_name]:
                boxplot.write('%s\t%s\n' % (sample_name, fpkm))
    boxplot.close()

def Density(sample_log10):
    density = open('density.txt', 'w')
    density.write('Sample\tlog10FPKM\tdensity\n')
    for sample_name in sample_log10:
        fpkm = sorted(sample_log10[sample_name])
        total_expr_gene = len(fpkm)
        for i in numpy.arange(fpkm[0], fpkm[-1] + 0.0001, 0.0001):
            a = 0
            newfpkm = {}
            for line in fpkm:
                if i <= line < i + 0.0001:
                    fpkm.remove(line)
                    fpkm_mean = i+0.00005 
                    a += 1
                    newfpkm[fpkm_mean] = a
            if newfpkm:
                (key, value), = newfpkm.items()
                value = float(value) / float(total_expr_gene)
                density.write('%s\t%s\t%s\n' % (sample_name, key, value))
    density.close()

def GeneNumber(sample):
    genenumber = open('genenumber.txt', 'w')
    genenumber.write('Sample\tFPKM_Range\tGeneNumber\n')
    for sample_name in sample:
        verylow = 0
        low = 0
        midhigh = 0
        for fpkm in sample[sample_name]:
            if fpkm <= 1:
                verylow += 1
            elif 1 < fpkm < 10:
                low += 1
            else:
                midhigh += 1
        genenumber.write('%s\tFPKM <=1\t%s\n%s\tFPKM 1~10\t%s\n%s\tFPKM >=10\t%s\n' % (sample_name, verylow, sample_name, low, sample_name, midhigh))
    genenumber.close()

if len(sys.argv) > 1:
    sample, sample_log10 = readgtf()
    #Boxplot(sample_log10)
    Density(sample_log10)
    #GeneNumber(sample)
else:
    print >>sys.stderr, "python %s gtf_path" % (sys.argv[0])
