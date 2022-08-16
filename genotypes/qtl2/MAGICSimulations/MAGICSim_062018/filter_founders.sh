#!/bin/bash


file=$1
founders=$2
cols=$(head -n1 "$file" | tr \\t \\n | egrep -nf $founders | cut -d: -f1 | paste -sd,)
echo $cols
cut -f1,$cols "$file" > founder_genos.txt
