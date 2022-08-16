#!/usr/bin/env python

import pandas as pd
import sys
import argparse
import numpy as np

def get_args():
    parser=argparse.ArgumentParser(description="""Program description""")
    parser.add_argument("infile",type=str,help="""The input ld file""")
    parser.add_argument("outfile",type=str,help="""The output ld file""")
    #parser.add_argument("markerfile",type=str,help=""""File with list of marker postions""")
    args = parser.parse_args()
    return args


def main():
    args=get_args()
    ldfile = pd.read_csv('{0}'.format(args.infile),sep='\t')
    samples = lddfile.columns[3:]
    header,stderr=bcftools_view(donorfile='hmp3_founders2/hmp3_founders_final.vcf.gz',header=True)
    print(stderr)
    for sample in samples:
        breaks = co_loc(bedfile[["start",sample]])
        vcf = header
        parents = bedfile[sample].unique()
        for i in parents:
            pbreaks=[j for j in breaks if j[3]==i]
            regionsfile='{0}_regions.txt'.format(i)
            marker_regions(pbreaks=pbreaks,markerfile=args.markerfile,rfile=regionsfile,c=10)
            positions,stderr=bcftools_view(donorfile='hmp3_founders2/{0}_c10_hmp321_final.vcf.gz'.format(i),regionsfile=regionsfile)
            print(stderr)
	    vcf+=positions
        with open('{0}_{1}'.format(sample,args.outfile),'w') as outfile:
            outfile.write(vcf)


if __name__=="__main__":
        main()
