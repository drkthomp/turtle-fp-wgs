#!/bin/bash
for SAMPLE in "$@"; do
    sort ${SAMPLE}_cnvs.txt | uniq | grep -v 'TRNA' > tmp.txt
    while read line; do
        echo "${line}" | grep -f - -w -m 1 ProteinTable_13308_1483792_Cm_NEW_NCBI.txt | awk 'BEGIN {FS="\t";OFS="\t"} {print $7,$11}';
    done < tmp.txt > ${SAMPLE}_cnvs_w_gene_descriptions.txt
done
rm tmp.txt
