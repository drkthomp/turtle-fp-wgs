find -name 'bwa_*_CheMyd_sorted_MITO.bam' -exec sh consensus.sh "{}" \;
find -name 'bwa_*_CheMyd_sorted_MITO.fq' -exec echo "{}" \; -exec sh remove.sh "{}" \;
find -name 'bwa_*_CheMyd_sorted_MITO.fq' -exec echo "{}" \; -exec sh fastqtofasta.sh "{}" \;
find -name '*.fasta' -exec echo "{}" \; -exec sh fasta.sh "{}" \;
rm combined.fasta
cat *.fasta > combined.fasta

