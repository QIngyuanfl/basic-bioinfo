#! /usr/bin/env python
import sys
import numpy as np
import os
def txt():
    info={}
    tox={}
    with open(sys.argv[1]) as f:
        for line in f:
            line=line.strip().split('\t')
            if len(line)>=3:
                info[line[1]]=[line[2].replace(',','.').lstrip('.')]
    with open(sys.argv[2]) as f:
        for line in f:
            line=line.strip().split()
            if line[0] in info.keys():
                info[line[0]].append(float(line[1]))
    for k,v in info.items():
        v[0]=v[0].split('.')
        del(v[0][-1])
        tox[k]=v
    return tox

def annot(tox):
    name=sys.argv[1].split('.')[0]
    out=open('%s.txt'% name,'w')
    with open('%s.annot'%name,'w') as f:
        phylo=[x[0] for x in tox.values()]
        length=[]
        for i in phylo:
            length.append(len(i))
        l=min(length)
        arr=np.array(phylo)
        f.write('clade_separation\t0.7\nbranch_thickness\t1.5\nbranch_bracket_depth\t1\nclade_marker_edge_color\t#555555\nclade_marker_edge_width\t0.5\nbranch_bracket_width\t1\n')
        f.write('clade_marker_size\t10\nclade_marker_edge_color\t#000000\n')
        total_plotted_degrees=input('Please adjust the total_plotted_degrees:(0~360):')
        annotation_font_size=input('please type in the font size,5 as default(no/3~12):')
        if annotation_font_size != 'no':
            f.write('annotation_font_size\t%s\n'% annotation_font_size)
        if annotation_font_size == 'no':
             f.write('annotation_font_size\t5\n')
        f.write('total_plotted_degrees\t%s\nstart_rotation\t270\n'%total_plotted_degrees)
        rm=[]
        for i in range (0,l-1):
            a=[x[i] for x in arr]
            if len(set(a)) == 1:
                b=a[0]
                rm.append(b)
        for x in phylo:
            for y in rm:
                x.remove(y)
        max_abundance=max([x[1] for x in tox.values()])
        m_abundance=max_abundance
        count =0
        while max_abundance >=50:
            max_abundance=max_abundance/2
            count+=1
            outer_layer=list(set([x[0][0] for x in tox.values()]))
        for x in outer_layer:
            color=raw_input('please type in a favor color for %s (none/color code):'%x)
            if color != 'none':
                f.write('%s\tannotation\t%s\n' % (x,x))
            for k,v in tox.items():
                for i in v[0]:
                    out.write('%s.'% i)
                out.write('%s\n' % k)
                f.write('%s\tclade_marker_size\t%f\n' %(k,float(v[1])/(3*count)))
                if v[0][0] == x and color != 'none':
                    if v[1] == m_abundance:
                        f.write('%s\tclade_marker_color\t%s\n'%(k,color))
                    f.write('%s\tannotation_background_color\t%s\n' % (v[0][0],color))
                    f.write('%s\tannotation_background_color\t%s\n' % (k,color))
                    f.write('%s\tannotation_background_color\t%s\n' % (v[0][-1],color))
                    f.write('%s\tannotation\t%s\n' % (k,k))
                    f.write('%s\tannotation\t%s\n' %(v[0][-1],v[0][-1]))
    out.close()
def shell():
    name=name=sys.argv[1].split('.')[0]
    with open('%s.sh'%name,'w') as f:
        f.write('graphlan_annotate.py --annot %s.annot %s.txt %s.xml\n'%(name,name,name))
        f.write('graphlan.py %s.xml %s.png --dpi 300 --size 3.5\n'%(name,name))
    os.system('sh %s.sh'%name)
    os.system('display %s.png'%name)
def main():
    if len(sys.argv)!=3:
        print 'Usage:format.py <tox> <gene abundance>'
    else:
        tox=txt()
        annot(tox)
        shell()
if __name__ == '__main__':
    main()
