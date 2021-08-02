while IFS= read -r line; do
    echo $line 
    bcftools view -r $line $1 > $(basename $1 .vcf.gz)_$line.vcf;
done < chromosomes.txt
