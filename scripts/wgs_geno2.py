#!/usr/bin/env python

"""
Takes a vcf file and converts it to a csv format with markers as columns and samples as rows
with nucleotide information in IUPAC format. Creates a file filled with NA values for WGS data to be used for imputation. Ensure that the samples within the input vcf file are in alphabetical order.
This csv file is formatted for use with R/qtl2
"""

import pandas as pd
from subprocess import PIPE,Popen,STDOUT
import sys
import argparse

def parse_args():
    """ -h for info on arguments
    """
    parser = argparse.ArgumentParser(description="""Program description: Takes a vcf file and converts it to a csv format with markers as columns and samples as rows             
with nucleotide information in IUPAC format Creates a file filled with NA values for WGS data in between genotype vlaues to be used for imputation. Ensure that the samples within the input vcf file are in alphabetical order.                                                             
This csv file is formatted for use with R/qtl2""")
    parser.add_argument("vcffile",type=str,help="""The vcf file containing genotype data for the 600K array""")
    parser.add_argument("outfile",type=str,help="""The output csv filename""")
    args=parser.parse_args()
    return args


#def discordance(founders):
#    discord={}
#    samples=[]
#    with open(founders,'r') as founderfile:
#        for line in founderfile:
#            samples.append(line[:-1])
#    for s in samples:
#        df=pd.read_table('biogemma/founders/{0}.diff.sites'.format(s),sep='\t')
#        discordant=df[df['N_DISCORD']!= 0]
#        for index,row in discordant.iterrows():
#            chrom=int(row['CHROM'])
#            pos=int(row['POS'])
#            if chrom in discord:
#                discord[chrom][pos]=True
#            else:
#                discord[chrom]={pos:True}
#    return discord


def call_bcftools(infile):
    """Reads vcf file and extracts data on samples and genotypes for each marker
    Returns a list of samples and list of markers and genotypes (info)
    """
    process1 = Popen(['bcftools','query','-l',infile],stdout=PIPE,stderr=STDOUT)
    stdout,stderr=process1.communicate()
    print(stderr)
    samples = stdout.split('\n')[:-1]
    process2 = Popen(['bcftools','query','-f','%ID\t%CHROM\t%POS[\tGT=%GT]\n',infile],stdout=PIPE,stderr=STDOUT)
    stdout,stderr=process2.communicate()
    output = stdout.split('\n')[:-1]
    print "Read vcf file: Contains info on {0} samples and {1} markers".format(len(samples),len(output))

    return samples,output

def get_positions(positions,output):
    """Reads in the wgs and array position files a returns a dictionary of dictionaries for each chromosome and each position within the chromosome. If a position is True, then it is within the 600K array, else is False.
    """
#    with open(wgsfile,'r') as file1:
#       for line in file1:
#            split=line.split('\t')
#            chrom=int(split[0])
#            pos=int(split[1][:-1])
#            marker='S{0}_{1}'.format(str(chrom),str(pos))
#            if chrom in positions:
#                positions[chrom][pos]={'in_array':False,'in_vcf':False,'marker':marker}
#            else:
#                positions[chrom]={pos:{'in_array':False,'in_vcf':False,'marker':marker}}
#    with open(arrayfile,'r') as file2:
#        for line in file2:
#            split=line.split('\t')
#            chrom=int(split[0])
#            pos=int(split[1][:-1])
#            marker='S{0}_{1}'.format(str(chrom),str(pos))
#            if chrom in positions:
#                positions[chrom][pos]={'in_array':True,'in_vcf':False,'marker':marker}
#            else:
#                positions[chrom]={pos:{'in_array':True,'in_vcf':False,'marker':marker}}
    for i in output:
        split=i.split('\t')
        marker=split[0]
        chrom=int(split[1])
        pos=int(split[2])
        if chrom in positions:
            positions[chrom][pos]={'marker':marker}
        else:
            positions[chrom]={pos:{'marker':marker}}
        count=1
        if '2/2' in split[3:]:
            for n in split[3:]:
                if './.' in n:
                    n='NA'
                elif '1' in n:
                    n='A'
                elif '2' in n:
                    n='B'
                else:
                    n='NA'
                positions[chrom][pos][count]=n
                count+=1
        else:
            for n in split[3:]:
                if './.' in n:
                    n='NA'
                elif '0' in n:
                    n='A'
                elif '1' in n:
                    n='B'
                else:
                    n='NA'
                positions[chrom][pos][count]=n
                count+=1
    return positions


def get_genofile():
    args=parse_args()
    positions={}
    samples,output=call_bcftools(args.vcffile)
    positions=get_positions(positions,output)
    for c in sorted(positions.keys()):
        geno={}
        txt='ind'
        for s in range(1,len(samples)+1):
            geno[s]=[samples[s-1]]
        for p in sorted(positions[c].keys()):
            txt+=','+positions[c][p]['marker']
            for s in range(1,len(samples)+1):
                geno[s].append(positions[c][p][s])
        print "Writing out chromosome {0}".format(c)            
        for j in geno.keys():
            txt+='\n'+ ','.join(geno[j])
        with open(args.outfile+'_chr'+str(c)+'.csv','w') as outfile:
            outfile.write(txt)



if __name__ == "__main__":
    get_genofile()


