#!/usr/bin/env python

import numpy as np
import argparse
import sys
from subprocess import Popen, PIPE

def get_args():
    parser=argparse.ArgumentParser(description="""Program description""")
    parser.add_argument("infile",type=str,help="""The input vcf file""")
    parser.add_argument("outfile",type=str,help="""The output regions file""")
    parser.add_argument("m",type=int,help="""The number of markers to select""")
    args=parser.parse_args()
    return args

def get_regions(infile):
    process=Popen(['bcftools','view', '-f','%POS\n',infile],stdout=PIPE,stderr=PIPE)
    stdout,stderr=process.communicate()
    print(stderr)
    full = stdout.split('\n')
    return full

def main():
    args=get_args()
    full = get_regions(args.infile)
    markers = np.random.choice(full,args.m,replace=False)
    markers = sorted(map(int,markers))
    with open(args.outfile,'w') as outfile:
        outfile.write(markers)




if __name__ == "__main__":
    main()
