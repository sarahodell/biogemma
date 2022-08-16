#!/usr/bin/env python

from subprocess import Popen,PIPE
import re

axiomfile='MAGICSim_062018/Axiom600KSNPs.txt'
biogemmafile='/group/jrigrp/Share/genotypes/biogemma/Geno_BALANCE_Axiom_600K_548475_SNP_Imputes_v.haploid_26.01.16.txt'


biogemma={}
with open(biogemmafile,'r') as bfile:
    count=1
    for line in bfile:
        if re.search('^BGA_ID',line)==None:
            rs = line.split()[0]
            biogemma[rs]=count
        else:
            header=('\t').join(line.split('\t')[1:])
        count+=1

axiom={}
with open(axiomfile,'r') as afile:
    count=1
    for line in afile:
        if re.search('^#',line)==None and re.search('^ID',line)==None:
            rs = line.split()[3]
            ax = line.split()[0]
            axiom[rs]=count
            axiom[ax]=count
        count+=1

index_list = []
for b in biogemma:
    if b in axiom:
        index_list.append([biogemma[b],axiom[b]])

txt='rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t'
txt+=header
for i in index_list:
    bnum=i[0]
    anum=i[1]
    q = '{0}p'.format(bnum)
    p = '{0}p'.format(anum)
    b_process = Popen(['sed','-n',q,biogemmafile],stdout=PIPE,stderr=PIPE)
    stdout,stderr=b_process.communicate()
    bline = stdout.split('\t')
    a_process = Popen(['sed','-n',p,axiomfile],stdout=PIPE,stderr=PIPE)
    stdout,stderr=a_process.communicate()
    ainfo= stdout.split('\t')
    rs=bline[0]
    snps=('\t').join(bline[1:])
    strand=['+' if ainfo[8]=='f' else '-' if ainfo[8]=='r' else '.'][0]
    rest = '.\t.\t.\t.\t.\t.\t'
    txt+='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(rs,ainfo[9],ainfo[5],ainfo[6],strand,rest,snps)

with open('biogemma_600K_genotypes.hmp.txt','w') as outfile:
    outfile.write(txt)
