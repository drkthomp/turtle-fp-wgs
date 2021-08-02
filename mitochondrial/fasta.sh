base=$(basename $1 .fasta)
echo $base
sed -i "1 s/.*/>${base}/" $1
