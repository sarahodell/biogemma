#!/usr/bin/env python

from subprocess import Popen, PIPE
import numpy as np
import pandas as pd
import sys
import argparse


def get_args():
    parser=argparse.ArgumentParser(description="""Program description""")
    parser.add_argument("infile",type=str,help="""The input vcf file""")
    parser.add_argument("outfile",type=str,help="""The output vcf file""")
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
    parents = s['donor1'].unique()
    for index,row in s.iterrows():
        locs.append([row['chr'],int(row['start']),int(row['end']),row['donor1']])
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
                regions+='{0}\t{1}\n'.format(chrom,m)
    with open(rfile,'w') as outfile:
    	outfile.write(regions)

        
def all_regions(pbreaks,rfile):
    txt=""
    for i in pbreaks:
        txt+='{0}\t{1}\t{2}\n'.format(i[0],i[1],i[2])
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
    header,stderr=bcftools_view(donorfile='{0}/Biogemma_Founders_600K_Genotypes_AGPv4_no_tester_chr10.vcf.gz'.format(args.donorpath),header=True)
    print(stderr)
    for sample in samples:
        vcf = str(header)
        breaks,parents = co_loc(sample,bedfile)
        for i in parents:
            pbreaks = [j for j in breaks if j[3]==i]
            regionsfile='{0}_{1}_regions.txt'.format(i,sample)
            if args.all == True:
                all_regions(pbreaks=pbreaks,rfile=regionsfile)
            else:
                marker_regions(pbreaks=pbreaks,markerfile=args.markerfile,rfile=regionsfile,c=10)
            positions,stderr=bcftools_view(donorfile='{0}/{1}_600K.vcf.gz'.format(args.donorpath,i),regionsfile=regionsfile)
            print(stderr)
            vcf+=str(positions)
        with open('tmp/{0}_{1}'.format(sample,args.outfile),'w') as outfile:
            outfile.write(vcf)


if __name__=="__main__":
        main()
