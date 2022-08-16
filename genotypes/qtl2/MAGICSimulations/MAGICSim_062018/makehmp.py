#!/usr/bin/env python


import sqlite3

db = sqlite3.connect('axiom600K.db')
cursor=db.cursor()
no_match=[]
txt='rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t'
with open('founder_genos.txt','r') as infile:
    for line in infile:
        l = line[:-1].split('\t')
        if l[0]=='BGA_ID':
            txt+=('\t').join(l[1:])+'\n'
        else:
            query = "SELECT chromosome,physical_position, tile_strand, alleles_f FROM snps WHERE cust_id='{0}'".format(l[0])
            cursor.execute(query)
            result = cursor.fetchall()
            if len(result) != 0:
                a = result[0]
                alleles=str(a[3])
                chrom=[a[0] if a[0] != -1 else '.'][0]
                pos=[a[1] if a[1] != -1 else '.'][0]
                strand=['+' if a[2]=='f' else '-' if a[2]=='r' else '.'][0]
                rest = '.\t.\t.\t.\t.\t.\t'
                geno=('\t').join(l[1:])
                txt+='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(l[0],alleles,chrom,pos,strand,rest,geno)
            else:
                no_match.append(b)

                
print("Writing to file")
db.close()

with open('biogemma_founder.hmp.txt','w') as outfile:
    outfile.write(txt)

if len(no_match) != 0:
    m=''
    print("{0} markers are not found in database".format(len(no_match)))
    for n in no_match:
        m+=n+'\n'
    with open('not_axiom600K.txt','w') as nofile:
        nofile.write(m)
