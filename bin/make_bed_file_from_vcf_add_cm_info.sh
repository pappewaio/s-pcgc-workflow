input_vcf=${1}
input_cm=${2}
output=${3}

# Convert from vcf to bed
plink2 --vcf ${input_vcf} --make-bed --max-alleles 2 --out ${output}.tmp

# Add Centimorgan info
plink --bfile ${output}.tmp --make-bed --cm-map ${input_cm}--out ${output}


