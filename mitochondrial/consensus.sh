echo $1
base=$(basename "$1" .bam)
bcftools mpileup -f ../../ref/GCF_015237465.1_rCheMyd1.pri_mitochondrial.fasta -d 3459810 -q 60 -Q 40 --threads 2 -o ${base}_mpileup.vcf $1 
bcftools call -c -o ${base}_call.vcf --threads 2 ${base}_mpileup.vcf
vcfutils.pl vcf2fq ${base}_call.vcf > ${base}.fq
