#!/bin/env python3
# -*- coding: utf-8 -*-

import sys
   
# process psf
try:
    infile = open(sys.argv[1],'r').readlines()
except:
    print('First command-line argument is not a proper psf file.')
    sys.exit()
    
numlp=0
numlph_list=[]
outlines = []

for curline in infile:
    curline_list = curline.split()
    if len(curline_list) <=6:
        outlines.append(curline)
        continue
    typename = curline_list[5]
    if typename == "OT":
        numlp+=1
        numlph_list.append(curline_list[0])
    outlines.append(curline)

outlines.append("    {}   {}   !{}  {}{}".format(numlp,numlp*3,"NUMLP","NUMLPH",'\n'))

index = 1
for lphindex in numlph_list:
    outlines.append('    3   '+str(index)+'    F    '+'-0.15000    '+'0.0000    '+'0.0000'+'\n')
    index = index+4;
    
index = 1
for lphindex in numlph_list:
    outlines.append("    "+str(int(lphindex)+3)+"    "+lphindex+"    "+str(int(lphindex)+1)+"    "+str(int(lphindex)+2)+'\n')

# write a new file
outfile = open(sys.argv[1][:-4]+'_lph_modified.psf','w')
for line in outlines:
    outfile.write(line)
outfile.close()

print('A file with modified lphost entries was written: {0}_lph_modified.psf'.format(sys.argv[1][:-4]))
