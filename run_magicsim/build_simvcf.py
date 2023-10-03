#!/usr/bin/env python

from subprocess import Popen, PIPE
import numpy as np
import pandas as pd
import sys
import argparse

# 9/22/23 Updated format and Popen to work in Python 3.9
def get_args():
    parser=argparse.ArgumentParser(description="""Program description""")
    parser.add_argument("infile",type=str,help="""The input vcf file""")
    parser.add_argument("outfile",type=str,help="""The output vcf file""")
    parser.add_argument("headerfile",type=str,help="""vcf file with header used for output file""")
    parser.add_argument("donorpath",type=str,help="""The path to the location of the donor files""")
    parser.add_argument("--markerfile",type=str,help=""""File with list of marker postions""")
    parser.add_argument("--all",type=bool,help="""If True, use all of the marker positions available in the donor files""")
    args = parser.parse_args()
    return args


def co_loc(sample,bedfile):
    """Returns regions of a chromosome received from particular parents in an individual RIL chromosome
    Input:
    breaks: pandas df of format: "start" "donor"
    c: chromosome number (int)

    Output:
    2-D list of form at [[chromsome,start,end,donor],...]
    """
    s = bedfile[bedfile['sample']==sample]
    locs=[]
    parents = s['donor'].unique()
    for index,row in s.iterrows():
        locs.append([row['chr'],int(row['start']),int(row['end']),row['donor']])
    return locs,parents


def marker_regions(pbreaks,markerfile,rfile,c=10):
    markers=[]
    with open(markerfile,'r') as infile:
        for line in infile:
            markers.append(line[:-1])
    regions=""
    for m in markers:
        for i in pbreaks:
            chrom=c
            start=i[1]
            end=i[2]
            if int(m) >= start and int(m) <=end:
                regions+=f"{chrom}\t{m}\n"
    with open(rfile,'w') as outfile:
    	outfile.write(regions)


def all_regions(pbreaks,rfile):
    txt=""
    for i in pbreaks:
        txt+=f"{i[0]}\t{i[1]}\t{i[2]}\n"
    with open(rfile,'w') as outfile:
        outfile.write(txt)


def bcftools_view(donorfile,regionsfile=None,header=False):
    if header ==True:
        process = Popen(['bcftools','view','-Ov','-h',donorfile],stdout=PIPE,stderr=PIPE)
    else:
        process = Popen(['bcftools','view','-Ov','-H','-R',regionsfile,donorfile],stdout=PIPE,stderr=PIPE)
    stdout,stderr = process.communicate()
    return stdout,stderr


def main():
    args=get_args()
    bedfile = pd.read_csv(args.infile,sep='\t')
    samples = bedfile['sample'].unique()
    header,stderr=bcftools_view(donorfile=f"{args.headerfile}",header=True)
    #print(str(stderr,'utf-8'))
    for sample in samples:
        vcf = str(header,'utf-8')
        breaks,parents = co_loc(sample,bedfile)
        for i in parents:
            #print(i)
            pbreaks = [j for j in breaks if j[3]==i]
            regionsfile=f"{i}_{sample}_regions.txt"
            if args.all == True:
                all_regions(pbreaks=pbreaks,rfile=regionsfile)
            else:
                marker_regions(pbreaks=pbreaks,markerfile=args.markerfile,rfile=regionsfile,c=10)
            positions,stderr=bcftools_view(donorfile=f"{args.donorpath}/{i}_600K.vcf.gz",regionsfile=regionsfile)
            #print(str(stderr,'utf-8'))
            vcf+=str(positions,'utf-8')
        with open(f"{sample}_{args.outfile}",'w') as outfile:
            outfile.write(vcf)


if __name__=="__main__":
        main()
