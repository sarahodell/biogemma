#!usr/bin/env python

### Assess FILLIN Projection Alignment
### Created by: Sarah Odell
### Date Created: 07/06/18
###

from subprocess import Popen, PIPE
import pandas as pd
import argparse

def arg_parse():
    parser=argparse.ArgumentParser(description="""Program description""")
    parser.add_argument("proja",type=str,help="""Projection alignment text file""")
    parser.add_argument("outfile",type=str,help="""Outfile name""")
    args = parser.parse_args()
    return args


def read_proja(proja):
    process = Popen(['awk','/#Donor Haplotypes/{flag=1;next}/#Taxa Breakpoints/{flag=0}flag',proja],stdout=PIPE,stderr=PIPE)
    stdout,sterr=process.communicate()
    donors={}
    for i in stdout.split('\n')[:-1]:
        line = i.split('\t')
        donors[line[0]]=line[1]
    process = Popen(['awk','/#Block are defined chr:startPos:endPos:donor1:donor2/{flag=1;next}/,0/{flag=0}flag',proja],stdout=PIPE,stderr=PIPE)
    stdout,sterr=process.communicate()
    breakpoints=[]
    for i in stdout.split('\n')[:-1]:
        line = i.split('\t')
        brks={}
        sample=line[0]
        for j in line[1:-1]:
            info=j.split(':')
            chrom=info[0]
            start=info[1]
            end=info[2]
            donor1=info[3]
            donor2=info[4]
            brks = {'sample':sample, 'chr':chrom, 'start':start, 'end':end, 'donor1': donors[donor1], 'donor2':donors[donor2]}
            breakpoints.append(brks)
    df = pd.DataFrame(breakpoints)
    projadf = df[['sample','chr','start','end','donor1','donor2']]
    return projadf

            
def main():
    args=arg_parse()
    breakpoints=read_proja(args.proja)
    breakpoints.to_csv(args.outfile,sep='\t',index=False)


if __name__ == "__main__":
    main()
