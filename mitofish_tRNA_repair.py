#!/usr/bin/env python
# coding: utf-8

__author__ = "Yonv1943"
__version__ = "2018-03-12 08:00"
"""
default:
python D:/Code/PyTest/mitofish_tRNA_repair.py -i mito.tbl mito.arwen -o mito_repair.tbl
"""
import os
import sys
import numpy as np
from functools import reduce


class DataIO:
    def __init__(self):
        self.infile, self.outfile = self.read_parameter()
        self.arwen_data = self.read_arwen()

        self.base_parting = dict([('A', 'T'), ('T', 'A'),
                                  ('C', 'G'), ('G', 'C')])

    def test(self):
        pass

    def read_parameter(self):
        """
        get parameter from command line
        :return: infile, outfile
        """

        # default parameter
        if len(sys.argv) == 1:
            print("\n||| Python Running", sys.argv,
                  "\n||| need to tell me the [inputData, outputData]"
                  "\n||| or i would used default parameter"
                  "\n||| -i *.tbl\t*.arwen\t[inputData]"
                  "\n||| -o repair.tbl\t\t[outputData]")

            cwd = os.path.dirname(__file__)
            listdir = os.listdir(cwd)

            file_tbl = [x for x in listdir if (x.find(".tbl") >= 0) and (x.find("_repair") < 0)][0]
            file_arwen = [x for x in listdir if x.find(".arwen") >= 0][0]

            sys.argv.extend([
                "-i %s %s" % (file_tbl, file_arwen),
                "-o %s_repair.tbl" % file_tbl[:-4],
            ])

        command = reduce(lambda a, b: a + b + " ", sys.argv[1:])
        command = command.split("-")[1:]

        infile = []
        outfile = []
        for item in command:
            if item[0] == 'i':
                infile = item[1:].split()
            elif item[0] == 'o':
                outfile = item[1:].split()

        return infile, outfile

    def read_arwen(self):
        arwen_pwd = [x for x in self.infile if x.find(".arwen") >= 0][0]
        # arwen_data = [ begin-end, protein(condon)+ ]
        arwen_data = []
        with open(arwen_pwd, 'r') as f:
            line = True
            while line:
                line = f.readline()
                find_end = line.find("tRNA-") + 5
                if find_end == 5:
                    arwen_data.append([line[5:-1], ])
                elif find_end > 5:
                    protein_condon = line[find_end:-1].split()[0]
                    if protein_condon[0] == '(':
                        seprarator_id = protein_condon.find("|")
                        arwen_data[-1].extend([protein_condon[1:seprarator_id] + protein_condon[-5:],
                                               protein_condon[seprarator_id+1:-6] + protein_condon[-5:], ])
                    else:
                        arwen_data[-1].append(protein_condon)

        # for i in arwen_data:
        #     print(i)
        # print(arwen_data[0])

        return arwen_data

    def repair_tbl(self):
        # read_tbl
        tbl_pwd = [x for x in self.infile if x.find(".tbl") >= 0][0]
        # tbl_repair_list = [ index_of_file, begin-end, protein ]
        tbl_repair_list = []
        with open(tbl_pwd, 'r') as f:
            indx = -1
            line = True
            while line:
                indx += 1
                line = f.readline()

                if line.find("\ttRNA\n") >= 0:
                    line_ary = line[:-1].split()
                    tbl_repair_list.append([indx, "%s-%s" % (line_ary[0], line_ary[1])])

                    indx += 1
                    line = f.readline()
                    tbl_repair_list[-1].append(line[line.find("-")+1:-1])

                    beg_end = 0

        tbl_repair_list = np.array(tbl_repair_list)
        # for i in tbl_repair_list:
        #     print(i)
        # print(tbl_repair_list[0])

        tbl_file = []
        with open(tbl_pwd, 'r') as f:
            tbl_file = list(f.readlines())

        # repair_tbl
        for tbl_file_id, beg_end, protein in tbl_repair_list:
            tbl_file_id = int(tbl_file_id)
            arwen_id = int(np.where(tbl_repair_list[:, 1] == beg_end)[0])

            protein_condon = self.arwen_data[arwen_id][1:]
            condon = [x for x in protein_condon if x.find(protein) == 0]

            if condon:
                condon = condon[0][-4:-1].upper()
            else:
                print("|||Outlier ", beg_end, protein, protein_condon)
                continue
            # print(protein, condon.upper(), tbl_file_id - 1)
            tbl_file[tbl_file_id - 1] = "%s(%s)\n" % (tbl_file[tbl_file_id - 1][:-1], condon)
            tbl_file[tbl_file_id + 1] = "%s(%s)\n" % (tbl_file[tbl_file_id + 1][:-1], self.reversed_pairing(condon))

        repair_tbl_pwd = [x for x in self.outfile if x.find(".tbl") >= 0][0]
        with open(repair_tbl_pwd, 'w') as f:
            f.writelines(tbl_file)

    def reversed_pairing(self, condon):
        return ''.join([self.base_parting[x] for x in condon[::-1]])


if __name__ == '__main__':
    dataIO = DataIO()
    dataIO.test()
    dataIO.repair_tbl()

pass
