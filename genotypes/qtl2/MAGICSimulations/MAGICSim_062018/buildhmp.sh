#!/bin/bash


axiom=$1
geno=$2

head=$(echo -e "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t")
head=$head$(head -n1 $geno | cut -f2-)
echo $head > 'biogemma_600K_genotypes.hmp.txt'

sed 1d $geno | while read line; do
    rs=$(echo $line | cut -d ' ' -f1)
    snps=$(echo $line | cut -d ' ' -f2-)
    ax="$(grep $rs $axiom)"
    info=$(echo $ax | cut -d ' ' -f 12,7,8)
    newline=$(echo -e "\n$rs\t$info\t.\t.\t.\t.\t.\t.\t.\t$snps")
    echo $newline >> 'biogemma_600K_genotypes.hmp.txt'
done


