#!/usr/bin/env python

import sys
import argparse
from subprocess import Popen,PIPE
import string
import pandas as pd


def parse_args():
    """ -h for info on arguments
    """
    parser = argparse.ArgumentParser(description="""Program description:
                                     Takes a vcf file and converts it to a csv format with founder as rows and marker genotypes as columns
                                     with nucleotide information encoded as A for reference allele and B for alternate alleles. Founder lines will be given letter codes based on the order that they are listed in the vcf file ("A" for the first sample and so on).
                                     This csv file is formatted for use with R/qtl2. It requires that bcftools be installed.
                                     """)
    parser.add_argument("infile",type=str,help="""The input vcf file""")
    parser.add_argument("outfile",type=str,help="""The output csv filename""")
    args=parser.parse_args()
    return args

def get_founders(vcf):
    """Gets a list of the samples in the vcf file"""
    process=Popen(["bcftools", "query","-l",vcf],stdout=PIPE,stderr=PIPE)
    stdout,stderr=process.communicate()
    print(stderr)
    founders=stdout.split('\n')[:-1]    
    return founders


def discordance(samples):
    discord={}
    for s in range(1,len(samples)+1):
        discord[s]={}
        df=pd.read_table('biogemma/founders/{0}.diff.sites'.format(samples[s-1]),sep='\t')
        discordant=df[df['N_DISCORD']!= 0]
        for index,row in discordant.iterrows():
            chrom=int(row['CHROM'])
            pos=int(row['POS'])
            if chrom in discord[s]:
                discord[s][chrom][pos]=True
            else:
                discord[s][chrom]={pos:True}
    return discord


def get_tmp(vcf):
    process=Popen(['bcftools','query','-f','%ID\t%CHROM\t%POS[\tGT=%GT]\n',vcf],stdout=PIPE,stderr=PIPE)
    stdout,stderr=process.communicate()
    info = stdout.split('\n')[:-1]
    print(stderr)
    return info

def geno_info(info,discord):
    positions={}
    for i in info:
        split=i.split('\t')
        chrom=int(split[1])
        pos=int(split[2])
        marker=split[0]
        if chrom in positions:
            positions[chrom][pos]={'marker':marker}
        else:
            positions[chrom]={pos:{'marker':marker}}
        count=1
        for n in split[3:]:
            if './.' in n:
                n='-'
            elif '1' in n:
                n='A'
            elif '2' in n:
                n='B'
            else:
                n='-'
                positions[chrom][pos][count]=n
                count+=1
    return positions

                
def get_foundergenofile():
    args=parse_args()
    founders=get_founders(args.infile)
    discord=discordance(founders)
    numf=len(founders)
    print(founders)
    f=''
    lcode=string.ascii_uppercase
    for s in range(1,numf+1):
        f+='{0},{1}\n'.format(founders[s-1],lcode[s-1])
    print "Writing founder codes to FounderCodes.csv"
    with open('FounderCodes.csv','w') as ffile:
        ffile.write(f)
    info=get_tmp(args.infile)
    array=get_tmp(args.arrayfile)
    positions=geno_info(info,discord)
    positions=geno_info(array)
    for c in sorted(positions.keys()):
        geno={}
        txt='ind'
        for s in range(1,numf+1):
            geno[s]=[founders[s-1]]
        for p in sorted(positions[c].keys()):
            txt+=','+'S'+str(c)+'_'+str(p)
            for s in range(1,numf+1):
                geno[s].append(positions[c][p][s])
        for j in geno.keys():
            txt+='\n'+ ','.join(geno[j])
        print "Writing out chromosome {0}".format(str(c))
        with open(args.outfile+'_chr'+str(c)+'.csv','w') as outfile:
            outfile.write(txt)




if __name__ == "__main__":
    get_foundergenofile()


