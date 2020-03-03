#!/usr/bin/env python

"""
Takes a vcf file and converts it to a csv format with markers as columns and samples as rows
with nucleotide information in IUPAC format.
This csv file is formatted for use with R/qtl2
"""

import pandas as pd
from subprocess import PIPE,Popen,STDOUT
import sys
import argparse

def parse_args():
    """ -h for info on arguments
    """
    parser = argparse.ArgumentParser(description="""Program description""")
    parser.add_argument("infile",type=str,help="""The input vcf file""")
    parser.add_argument("outfile",type=str,help="""The output csv filename""")
    args=parser.parse_args()
    return args


def call_bcftools(infile):
    """Reads vcf file and extracts data on samples and genotypes for each marker
    Returns a list of samples and list of markers and genotypes (info)
    """
    process1 = Popen(['bcftools','query','-l',infile],stdout=PIPE,stderr=STDOUT)
    stdout,stderr=process1.communicate()
    print(stderr)
    samples = stdout.split('\n')[:-1]
    process2 = Popen(['bcftools','query','-f','%ID[\tGT=%GT]\n',infile],stdout=PIPE,stderr=STDOUT)
    stdout,stderr=process2.communicate()
    info = stdout.split('\n')[:-1]
    print(stderr)
    print "Read vcf file: Contains info on {0} samples and {1} markers".format(len(samples),len(info))
    return samples,info

def get_genofile():
    args=parse_args()
    samples,info = call_bcftools(args.infile)
    txt='ind'
    geno={}
    for s in range(1,len(samples)+1):
        geno[s] = [samples[s-1]]
    for i in info:
        split=i.split('\t')
        marker = split[0]
        txt+=','+marker
        count=1
        for n in split[1:]:
            if './.' in n:
                n='NA'
            elif '0' in n:
                n='A'
            else:
                n='B'
            geno[count].append(n)
            count+=1
    for j in geno.keys():
        txt+='\n'+ ','.join(geno[j])
    print "Writing out to {0}".format(args.outfile)
    with open(args.outfile,'w') as outfile:
        outfile.write(txt)



if __name__ == "__main__":
    get_genofile()


