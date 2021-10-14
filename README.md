# s-pcgc-workflow
An attempt to run s-pcgc from scratch using 1000 Genomes data


## Installing the environment

Use conda/mamba to install the prerequisites to run S-PCGC

```
# Install environment for s-pcgc
mamba create -n s-pcgc python=3.9 numpy pandas scikit-learn scipy tqdm bitarray pandas-plink plink==1.90b6.18 plink2==2.00a2.3 --channel bioconda

# Install environmane for ldsc (needed to generate some reference data)
# Reading env from file using mamba requires an empty env and then updating
mamba create -n ldsc
mamba env update -n ldsc --file repos/ldsc/environment.yml

```

## Using interactive node on HPC using slurm
When testing the code on a HPC, remeber to start an interactive node, like this:

```
projectname="some_name"
srun --mem=40g --ntasks1 --cpus-per-task 6 --time=5:00:00 --account ${projectname} --pty /bin/bash
```
## Run the basic example described in S-PCGC readme

```
## Run toy example 

spcgc="repos/S-PCGC"
mkdir temp_results

#create a sync file (a one-time offline operation)
python ${spcgc}/pcgc_sync.py \
--annot-chr ${spcgc}/example/model. \
--frqfile-chr ${spcgc}/example/model. \
--out temp_results/model

#compute cross-r^2 between functional annotations
python ${spcgc}/pcgc_r2.py \
--bfile ${spcgc}/example/ref_panel \
--annot-chr ${spcgc}/example/model. \
--sync temp_results/model. \
--out temp_results/prodr2

#Compute summary statistics for study 1
python pcgc_sumstats_creator.py \
--bfile example/s1 \
--pheno example/s1.phe \
--covar example/s1.cov \
--frqfile-chr example/model. \
--annot-chr example/model. \
--sync temp_results/model. \
--prev 0.01 \
--out temp_results/s1

#Compute summary statistics for study 2
python pcgc_sumstats_creator.py \
--bfile example/s2 \
--pheno example/s2.phe \
--covar example/s2.cov \
--frqfile-chr example/model. \
--annot-chr example/model. \
--sync temp_results/model. \
--prev 0.01 \
--out temp_results/s2

#Run S-PCGC to estimate h^2, rg and functional enrichment
python pcgc_main.py \
--annot-chr example/model. \
--sync temp_results/model. \
--sumstats temp_results/s1.,temp_results/s2. \
--prodr2 temp_results/prodr2. \
--out temp_results/results

#view heritability estimates for the two studies
cat temp_results/results.s1.output | column -t
cat temp_results/results.s2.output | column -t

#view a table of genetic correlation estimates between the studies
cat temp_results/results.rg

#view the functional enrichment estimates of the two studies
cat temp_results/results.s1.results | column -t
cat temp_results/results.s2.results | column -t


```

## Prepare 1000 Genomes data, which is going to be used in the real annotations example
This conversion will make hard calls, which can have side-effects for example, the allele frequency if you use imputed/dosage data. In this case the 1000G VCFs are not imputed so the conversion will be ok.

```
# Download and place in data/vcf (or symlink from a common repository)

```

```
# Create a run file, indicaing ,which chromosome each file has

```
Convert from vcf to bed format and add centimorgan info

```
# Run using example data
nextflow insert_map_in_vcf.nf --input data/runfile/runfile.txt

```

## Prepare the "good set of snps" using the recommended list built from hapmap3 snps
```
mkdir data/hapmap3
cd data/hapmap3
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/hapmap3_snps.tgz
tar -zxvf hapmap3_snps.tgz
```

## Run the real annotations example using the 1000G data

### prepare 1000G a bit further

```
# QC filtering
Might be worth doing some QC here(more than removing duplicates)

## drop duplicated IDs
# Use suppress-first, if you wnat to keep the first found element in each duplicated variant
# Here we do the conservative filtering removing indels overlapping the position, it can be wise to filter out indels first. 
# NOTE: it might be that the old plink bfile format removes indels by default (but that has to be checked)
# A full description of plink2 variant filtering can be found here: https://www.cog-genomics.org/plink/2.0/filter#variant

# remove all non 'snps'
PLINK_BFILE_folder=""
mkdir -p ${out_KGP}/snpsonly
for chr1-22
  do
    plink2 \
    --bfile <file> \
    --snps-only \
    --make-bed
    --out ${out_KGP}/snpsonly
done

# remove duplicates
mkdir -p  ${out_KGP}/duprm
for chr1-22
  do
    plink2 \
    --bfile ${out_KGP}/snpsonly/1000G.EUR.${chr} \
    --rm-dup excluding-all list \
    --out ${out_KGP}/duprm/1000G.EUR.${chr}
done

rm "${out_KGP}/duprm/1000G.EUR.all.rmdup.list"
for chr1-22
do
    if[ -f "${out_KGP}/duprm/1000G.EUR.${chr}.rmdup.list" ]
      cat "${out_KGP}/duprm/1000G.EUR.${chr}.rmdup.list" >> "${out_KGP}/duprm/1000G.EUR.all.rmdup.list"
    fi
    # add a dot for all missing rsids
    echo "." >> "${out_KGP}/duprm/1000G.EUR.all.rmdup.list"
done

# check for duplicates across all chromosomes
for 1-22
do
plink2 \
    --bfile ${out_KGP}/snpsonly/1000G.EUR.${chr} \
    --rm-dup excluding-all list \
    --out ${out_KGP}/duprm/1000G.EUR.${chr}


```

### Create annotation files to compute ld scores


```
# Annotation source premade by ldsc-team
mkdir -p data/annotation
cd data/annotation
wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/make_annot_sample_files/ENSG_coord.txt

```

I recommend to not use the annotation file creator, as it appears to create the wrong result. The reason is that it doesn't support chr, bp and cm info. For now I only want to use all snps, so only using the base case annotation setting 1 for all rows would create that. And if there are any partitioned estimates we are interested, then it is better anyways that we create them ourselves than relying on the provided script (make_annot.py)

So instead we will add annotation to the extracted fields of the bim file. As previously mentionen, the simplest annotation is the base case, where everything is just 1. So in the new annot file we will have this header "#CHR BP SNP CM ANNOT".

for 1-22;
do 
  awk -vOFS="\t" -vheader="CHR	BP	SNP	CM	ANNOT" 'BEGIN{print header}; {print $1, ${out_KGP}/duprm/1000G.EUR.${chr}}' | gzip -c > out/KGP_preparation/full_annotation_LD/baselineLD.${chr}.annot.gz"
done


### Compute LD scores using an annot file 
Takes time Chr22 ~ 22min. Previous success using a walltime of 15h.

```

# Generate example data
```
# Download and take out a random subset from of 1000G, one portion from each chromosome

# Download HapMap3 genetic maps



```



