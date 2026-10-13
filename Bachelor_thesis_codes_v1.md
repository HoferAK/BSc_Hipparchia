# POPULATION GENOMICS *Hipparchia* Butterflies Material and Methods
This script contains codes and some notes that were used during the course of my Bachelor thesis. 

The sources used to create this code (apart from assistance of [ChatGPT](https://chatgpt.com) and [Claude](https://claude.ai/)):
> [POPULATION GENOMICS Heliconius](https://hackmd.io/@EnCpOBucSKGlaCt5-riPXQ/r15hEyBx1g#2-Mapping)
> [Speciation & Population Genomics: a how-to-guide](https://speciationgenomics.github.io/)
> [KAY Lucek Pipeline](https://hackmd.io/UbLfs1AfSTqCDmIOjWhm6Q#GENERATE-VCF-FILES)

## Cleaning files with Fastp
The fist thing to do is to clean and trim the files with fastp, later on multiqc is used to summarize the overall qualtiy of the files. Note that unless specified [Metacentrum](https://metavo.metacentrum.cz/) is used to run the codes.
In this code a Poly-G text file is used to get rid of the adapter, since the polyG filter was not enough (based on [POPULATION GENOMICS Heliconius](https://hackmd.io/@EnCpOBucSKGlaCt5-riPXQ/r15hEyBx1g#2-Mapping)).
```
#!/bin/bash
#PBS -N fastp_all
#PBS -l select=1:ncpus=1:mem=10gb:scratch_local=50gb
#PBS -l walltime=10:00:00
#PBS -e PBS/
#PBS -o PBS/
#PBS -J 1-60

DATADIR=/storage/frontends/home/username/bachelor_thesis/files_raw/hipparchia_fastq

#clean scratch after the end
trap 'clean_scratch' TERM EXIT

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# copy input file to scratch directory
# if the copy operation fails, issue error message and exit

module load fastp

# go to  scratch directory
cd $SCRATCHDIR || exit 1
####

## Defining the file to process using the Array ID
file=$(ls $DATADIR/*_R1_001.fastq.gz | sed -n "${PBS_ARRAY_INDEX}p")

f=${file/_R1_001.fastq.gz/}
name=$(basename ${f})

## Forward and reverse
FORWARD=$file
REVERSE=${f}_R2_001.fastq.gz

out=/storage/frontends/home/username/metacentrum_output_fastp

fastp --thread 1 --in1 $FORWARD --in2 $REVERSE \
-l 50 --detect_adapter_for_pe \
--trim_poly_g \
--adapter_fasta /storage/frontends/home/username/my_code1/adapter.txt \
--out1 ${out}/${name}_1R.fastq.gz --out2 ${out}/${name}_2R.fastq.gz \
-h /storage/frontends/home/username/metacentrum_output_fastp/html/${name}.html \
&> /storage/frontends/home/username/metacentrum_output_fastp/log/${name}.log  

## Run FastQC
out1=${out}/${name}_1R.fastq.gz
out2=${out}/${name}_2R.fastq.gz

module load fastqc
fastqc -t 1 $out1 $out2 -o /storage/frontends/home/username/metacentrum_output_fastp/FastQC
```
Since the used memory was not sufficient and some files were not finished and vanished, another code was rerun for the files that were corrupted and could not be used for fastqc.
The missing files were found and extracted to a text file:
```
for f in "$folderA"/*.fastq.gz; do     base=$(basename "$f" .fastq.gz);      zipfile="$folderB/${base}_fastqc.zip";      if [[ -f "$zipfile" ]]; then         :;     else         echo "Missing: $zipfile" >> empty2.txt; fi; done 
```
Then I created a folder with the missing files and ran the same code as above with increased memory settings for the missing files.

## MultiQC
With this code from [POPULATION GENOMICS Heliconius](https://hackmd.io/@EnCpOBucSKGlaCt5-riPXQ/r15hEyBx1g#2-Mapping) the multiqc was run directly on the frontend in the folder with the final files for the *.zip files:
```
module load mambaforge
mamba activate multiqc_v1.12_py3.7
multiqc .
mamba deactivate
```

## Download Reference Genome
The reference Genomes of H. semele and H. statilinus were downloaded and unzipped:
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/933/228/805/GCA_933228805.2_ilHipSeme1.2/GCA_933228805.2_ilHipSeme1.2_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/971/852/105/GCA_971852105.1_ilHipStat1.hap1.1/GCA_971852105.1_ilHipStat1.hap1.1_genomic.fna.gz

# navigate in the folder with the genomes on the frontend

gunzip GCA_933228805.2_ilHipSeme1.2_genomic.fna.gz

gunzip GCA_971852105.1_ilHipStat1.hap1.1_genomic.fna.gz
```


### Index the Reference Genome
The indexing was fast using screen with the following code ([POPULATION GENOMICS Heliconius](https://hackmd.io/@EnCpOBucSKGlaCt5-riPXQ/r15hEyBx1g#2-Mapping) and [Speciation & Population Genomics: a how-to-guide](https://speciationgenomics.github.io/)):
```
module load bwa-mem2 
screen -S indexstatilinus
bwa-mem2 index GCA*.fna
ctrl+a+d
screen -rd indexstatilinus #to go back to the screen
```
the same was done for H. semele in a seperate screen.


## Mapping Reads
Now it is crucial to seperate the two species in seperate folders and work on the two species individually. Only one of the two codes is written here for better readability, since the two codes are essentially the same, only the file paths need to be adjusted. 
The mapping was done following [POPULATION GENOMICS Heliconius](https://hackmd.io/@EnCpOBucSKGlaCt5-riPXQ/r15hEyBx1g#2-Mapping).
```
#!/bin/bash
#PBS -N MappReadssemele
#PBS -l select=1:ncpus=2:mem=50gb:scratch_local=50gb
#PBS -l walltime=24:00:00
#PBS -e PBS/
#PBS -o PBS/
#PBS -J 0-44

# Input data directory
DATADIR="/storage/frontends/home/username/allign_semele/clean_files"

# Defining the file to process using the Array ID
file=$(ls ${DATADIR}/*_1R.fastq.gz | sed -n "${PBS_ARRAY_INDEX}p")

# Load required programs
module load bwa-mem2
module load samtools

# Output directory
output="/storage/frontends/home/username/allign_semele/fastq_mapping/bam_2.5"

# Reference genome (already indexed)
ref="/storage/frontends/home/username/reference_genomes/H_semele/GCA_933228805.2_ilHipSeme1.2_genomic.fna"

# Aligning individual sequences to the reference
f=${file/_1R.fastq.gz/}
name=$(basename ${f})
reverse=${file/_1R.fastq.gz/_2R.fastq.gz}

bwa-mem2 mem -t 2 -M -a ${ref} ${file} ${reverse} > ${output}/${name}.sam

# Convert sam to bam
samtools view -bS ${output}/${name}.sam -o ${output}/${name}.bam
rm ${output}/${name}.sam

# Sort the bam file
samtools sort ${output}/${name}.bam -o ${output}/${name}.sorted.bam
rm ${output}/${name}.bam

# Clean scratch after the end
trap 'clean_scratch' TERM EXIT
```


### Mapping Statistics
To run the Stats [POPULATION GENOMICS Heliconius](https://hackmd.io/@EnCpOBucSKGlaCt5-riPXQ/r15hEyBx1g#2-Mapping) was used.
```
#!/bin/bash
#PBS -N StatsMappingReadsSemele
#PBS -l select=1:ncpus=1:mem=50gb:scratch_local=100gb
#PBS -l walltime=24:00:00


DATADIR=/storage/frontends/home/username/allign_semele/

#clean scratch after the end
trap 'clean_scratch' TERM EXIT

module load samtools

# go to  scratch directory
cd $SCRATCHDIR || exit 1

output="$DATADIR/Stats"

touch ${output}/res.aln.reads.out
touch ${output}/reads_mapping_Qscores.txt
touch ${output}/mean_coverage.txt
touch ${output}/breadth_cov.txt

for file in $DATADIR/all_mapped/*.sorted.bam
 do
        f=${file/.sorted.bam/}
        name=`basename ${f}`

        # Number of reads with R1 and R2 mapping close to each other and other parameters
        stats=`samtools flagstat ${file} | awk '{print $1}' | column -x`
        echo ${name} ${stats} >> ${output}/res.aln.reads.out \
        && \
        # Mapping quality
        qual=`samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'`
        echo ${name} ${qual} >> ${output}/reads_mapping_Qscores.txt \
        && \
        # Nucleotide coverage
        cov=`samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'`
        echo ${name} ${cov} >> ${output}/mean_coverage.txt \
        && \
        # Breadth of coverage
        breadthcov=`samtools depth -a ${file} | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'`
        echo ${name} ${breadthcov} >> ${output}/breadth_cov.txt
done
```

### Removal of Duplicates
Instead of using a list and looping over each file individually like in [POPULATION GENOMICS Heliconius](https://hackmd.io/@EnCpOBucSKGlaCt5-riPXQ/r15hEyBx1g#2-Mapping), I used the PBS batch system in order to be faster. 
```
#!/bin/bash
#PBS -N RemovePCRduplicatessemele
#PBS -l select=1:ncpus=1:mem=65gb:scratch_local=200gb
#PBS -l walltime=10:00:00
#PBS -e /storage/frontends/home/username/PBS/
#PBS -o /storage/frontends/home/username/PBS/
#PBS -J 0-45

# Directories
DATADIR="/storage/frontends/home/username/allign_semele/all_mapped"
OUTPUT="/storage/frontends/home/username/allign_semele/finalbam"

# Load modules
module load samtools
module load picard/2.27.5-gcc-10.2.1-zpit36s

# Check scratch
test -n "$SCRATCHDIR" || { echo >&2 "SCRATCHDIR not set"; exit 1; }
cd "$SCRATCHDIR"

mkdir tmp.out
export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=$SCRATCHDIR"

# Select BAM file for this array job
file=$(ls ${DATADIR}/*.sorted.bam | sed -n "${PBS_ARRAY_INDEX}p")

# Safety check
[ -z "$file" ] && { echo "No file for index ${PBS_ARRAY_INDEX}"; exit 2; }

name=$(basename "${file/.sorted.bam/}")

# Copy input BAM to scratch
cp "$file" . || { echo "Failed to copy $file"; exit 3; }

# Remove PCR duplicates
picard MarkDuplicates \
    -REMOVE_DUPLICATES true \
    -ASSUME_SORTED true \
    -VALIDATION_STRINGENCY SILENT \
    -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
    -INPUT ${name}.sorted.bam \
    -OUTPUT tmp.out/${name}.rmd.bam \
    -METRICS_FILE tmp.out/${name}.rmd.metrics.txt

# Sort final BAM
samtools sort tmp.out/${name}.rmd.bam -o tmp.out/${name}.final.bam

# Index final BAM
samtools index tmp.out/${name}.final.bam

# Copy results back
cp tmp.out/${name}.final.bam* "$OUTPUT" || { echo "Copy failed"; exit 4; }

# Clean scratch
trap 'clean_scratch' TERM EXIT
```


### Stats after Removal of Duplicates
```
#!/bin/bash
#PBS -N StatsMappingReadsSemelermdup
#PBS -l select=1:ncpus=1:mem=50gb:scratch_local=100gb
#PBS -l walltime=24:00:00


DATADIR=/storage/frontends/home/username/allign_semele/

#clean scratch after the end
trap 'clean_scratch' TERM EXIT

module load samtools

# go to  scratch directory
cd $SCRATCHDIR || exit 1

output="$DATADIR/Stats2rmdup"

touch ${output}/res.aln.reads.out
touch ${output}/reads_mapping_Qscores.txt
touch ${output}/mean_coverage.txt
touch ${output}/breadth_cov.txt

for file in $DATADIR/finalbam/*.final.bam
 do
        f=${file/.final.bam/}
        name=`basename ${f}`

        # Number of reads with R1 and R2 mapping close to each other and other parameters
        stats=`samtools flagstat ${file} | awk '{print $1}' | column -x`
        echo ${name} ${stats} >> ${output}/res.aln.reads.out \
        && \
        # Mapping quality
        qual=`samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'`
        echo ${name} ${qual} >> ${output}/reads_mapping_Qscores.txt \
        && \
        # Nucleotide coverage
        cov=`samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'`
        echo ${name} ${cov} >> ${output}/mean_coverage.txt \
        && \
        # Breadth of coverage
        breadthcov=`samtools depth -a ${file} | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'`
        echo ${name} ${breadthcov} >> ${output}/breadth_cov.txt
done
```

## Variant Calling
Here I worked with bcf-tools mainly following [Speciation & Population Genomics: a how-to-guide](https://speciationgenomics.github.io/) while modifying some parts with AI.

### Indexing with Samtools
For bcftools another index is necessary, the previously used index does not work here:
```
module load samtools
samtools faidx ${ref}
```
To get the naming of the chromosomes:
```
grep ">" ${ref}
```

Create a textdocument with line one corresponding to chromosome one etc. This will be used for calling SNPs in a bit.
It is important to also create a bam_file.txt with all the paths and names of the files:
```
ls /storage/frontends/home/username/allign_statilinus/finalbam/*.bam > /storage/frontends/home/username/allign_statilinus/finalbam/bam_list.txt
```

### SNP calling
Now the SNP calling following [KAY Lucek Pipeline](https://hackmd.io/UbLfs1AfSTqCDmIOjWhm6Q#GENERATE-VCF-FILES)
```
#!/bin/bash
#PBS -N SNPcallingsemele
#PBS -l select=1:ncpus=1:mem=60gb:scratch_local=150gb
#PBS -l walltime=48:00:00
#PBS -e /storage/frontends/home/username/PBS/
#PBS -o /storage/frontends/home/username/PBS/
#PBS -J 1-30

DATADIR=/storage/frontends/home/username/reference_genomes

## if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
## clean scratch after the end
trap 'clean_scratch' TERM EXIT
cd $SCRATCHDIR || exit 1

## Set scaffold to run with PBS_ARRAY_INDEX
chr_file="$DATADIR/textbcf/semele.txt"
chr=$(cat $chr_file | sed -n "${PBS_ARRAY_INDEX}p")

## cp ref genome to scratchdir
cp /storage/frontends/home/username/reference_genomes/samtools/semele/GCA_933228805.2_ilHipSeme1.2_genomic.* .

## reference genome
ref="GCA_933228805.2_ilHipSeme1.2_genomic.fna" 

## modules
module load bcftools

## SNP calling with BCFtools. Allsites calling (no `bcftools call -v`)
bcftools mpileup -Ou -f ${ref} -b /storage/frontends/home/username/allign_semele/finalbam/bam_list.txt -r ${chr} --threads 1 -q 20 -a AD,DP,SP | bcftools call --threads 1 -m -Oz -a GQ,GP -o Semele_${chr}.vcf.gz

cp Semele_${chr}.vcf.gz /storage/frontends/home/username/SNP/bcf/semele

# Clean scratch
trap 'clean_scratch' TERM EXIT
```

### Filtering with VCFtools
Filtering was done following [KAY Lucek Pipeline](https://hackmd.io/UbLfs1AfSTqCDmIOjWhm6Q#GENERATE-VCF-FILES) and [speciationgenomics.github](https://speciationgenomics.github.io/variant_calling/). 
First a text file with the name of the Chromosomes was created using `grep ">" $file`. The the names of the Chromosomes were copied to a new textfile. Additionally a textfile with the paths to the final bam files was created with `ls /storage/frontends/home/username/allign_semele/finalbam/*.bam > /storage/frontends/home/username/allign_semele/finalbam/bam_list.txt` (AI assisted). These textfiles were used in the following code:
```
#!/bin/bash
#PBS -N Filtering_OnlyVariantsSemele
#PBS -l select=1:ncpus=1:mem=30gb:scratch_local=100gb
#PBS -l walltime=24:00:00
#PBS -e /storage/frontends/home/username/PBS/
#PBS -o /storage/frontends/home/username/PBS/
#PBS -J 1-30

DATADIR=/storage/frontends/home/username/reference_genomes


test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
## clean scratch after the end
trap 'clean_scratch' TERM EXIT
cd $SCRATCHDIR

## Set scaffold to run with PBS_ARRAY_INDEX
chr_file="$DATADIR/textbcf/semele.txt"
chr=$(cat $chr_file | sed -n "${PBS_ARRAY_INDEX}p")

in="/storage/frontends/home/username/SNP/bcf/semele"
out="/storage/frontends/home/username/SNP/bcf/semele_filtered" 
#repeats="/storage/brno12-cerit/home/repulles/03.Heliconius/02.1.fastqMapping/00.Hmel2.5_lepbase_repeats/repeats_all.bed"

## modules
module load vcftools

### First filtering keeping only the big scaffolds and with minGQ 20 (problems with it?)
vcftools --gzvcf ${in}/Semele_${chr}.vcf.gz \
	--minGQ 20 \
	--minDP 4 --maxDP 100 --max-missing 0.9 \
	--remove-indels --max-alleles 2 \
	--min-alleles 2 --maf 0.03 \
	--recode --recode-INFO-all \
	--out ${out}/Semele03_${chr}
gzip ${out}/Semele03_${chr}.recode.vcf 

# Clean scratch
trap 'clean_scratch' TERM EXIT
```

### VCF Stats
This step is based on [speciationgenomics.github](https://speciationgenomics.github.io/variant_calling/) and [KAY Lucek Pipeline](https://hackmd.io/UbLfs1AfSTqCDmIOjWhm6Q#GENERATE-VCF-FILES).
To have an insight into the quality of the SNPs the VCF Stats were run in MetaCentrum:
```
#run the stats on the frontend, not too demanding only a few seconds

DATADIR="/storage/frontends/home/username/SNP/bcf/statilinus_autosomes/stats"
OUT="/storage/frontends/home/username/SNP/bcf/statilinus_autosomes/vcftools"
VCF_FILE="${DATADIR}/SNPstatilinus_autosomes.vcf.gz"

module load vcftools

#calculate allele frequency:
vcftools --gzvcf $VCF_FILE --freq2 --out $OUT --max-alleles 2

#calculate mean depth per individual:
vcftools --gzvcf $VCF_FILE --depth --out $OUT

#calculate mean depth per site:
vcftools --gzvcf $VCF_FILE --site-mean-depth --out $OUT

#calculate site quality:
vcftools --gzvcf $VCF_FILE --site-quality --out $OUT

#calculate proportion of missing data per individual:
vcftools --gzvcf $VCF_FILE --missing-indv --out $OUT

#calculate proportion of missing data per site:
vcftools --gzvcf $VCF_FILE --missing-site --out $OUT


#Calculate heterozygosity and inbreeding coefficient per individual:
vcftools --gzvcf $VCF_FILE --het --out $OUT
```
and now the visualization of the Stats using R and AI:
```

# ===============================
# VCFtools QC plotting script
# ===============================

# 1) Load libraries
library(tidyverse)

# 2) Set working directory to where your VCFtools outputs are
setwd("C:/Users/username1/Documents/BT_R/vcf_stats2/stats_vcf_semele/")


# 3) Read VCFtools output files

# Site quality
var_qual <- read_delim("stats_vcf_semele.lqual", delim = "\t",
                       skip = 1,
                       col_names = c("chr","pos","qual"))

# Mean depth per site
site_depth <- read_delim("stats_vcf_semele.ldepth.mean", delim = "\t",
                         skip = 1,
                         col_names = c("chr","pos","mean_depth"))

# Depth per individual
ind_depth <- read_delim("stats_vcf_semele.idepth", delim = "\t",
                        skip = 1,
                        col_names = c("sample","N_called","mean_depth","sum_depth"))

# Missing per individual
ind_missing <- read_delim("stats_vcf_semele.imiss", delim = "\t",
                          skip = 1,
                          col_names = c("sample","N_called","N_GENOTYPES_FILTERED", "N_missing","frac_missing"))


var_miss <- read_delim("stats_vcf_semele.lmiss", delim = "\t",
                       
                       col_names = c("chr","pos","nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)


ind_het <- read_delim("stats_vcf_semele.het", delim = "\t",
                      skip = 1,
                      col_names = c("ind","ho","he","nsites","f"))

var_freq <- read_delim("stats_vcf_semele.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))



# 4) Create density plots

p1 <- ggplot(var_qual, aes(qual)) +
  geom_density(fill = "dodgerblue1", alpha = 0.5) +
  coord_cartesian(xlim = c(0, 3000)) +
  labs(title = "Site Quality Distribution") +
  theme_light()

p2 <- ggplot(site_depth, aes(mean_depth)) +
  geom_density(fill = "forestgreen", alpha = 0.5) +
  labs(title = "Mean Depth per Site") +
  theme_light()

p3 <- ggplot(ind_depth, aes(mean_depth)) +
  geom_density(fill = "orange", alpha = 0.5) +
  labs(title = "Mean Depth per Individual") +
  theme_light()

p4 <- ggplot(ind_missing, aes(frac_missing)) +
  geom_density(fill = "firebrick", alpha = 0.5) +
  labs(title = "Missing Data per Individual") +
  theme_light()

p55 <- ggplot(var_miss, aes(fmiss)) +
  geom_density(fill = "purple", alpha = 0.5) +
  labs(title = "Missing Data per Site") +
  theme_light()

p66 <- ggplot(ind_het, aes(f)) +
  geom_histogram(fill = "darkcyan", alpha = 0.5) +
  labs(title = "Inbreeding Coefficient (F)") +
  theme_light()

p7 <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(title = "MAF") 
p7 + theme_light()

# 5) Display all plots
print(p1)
print(p2)
print(p3)
print(p4)
print(p55)
print(p66)
print(p7)

# 6) Optional: save all plots to pngs
ggsave("01_site_quality.png", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("02_mean_depth_per_site.png", plot = p2, width = 8, height = 6, dpi = 300)
ggsave("03_mean_depth_per_individual.png", plot = p3, width = 8, height = 6, dpi = 300)
ggsave("04_missing_per_individual.png", plot = p4, width = 8, height = 6, dpi = 300)
ggsave("05_missing_per_site.png", plot = p55, width = 8, height = 6, dpi = 300)
ggsave("06_inbreeding_coefficient.png", plot = p66, width = 8, height = 6, dpi = 300)
ggsave("07_maf.png", plot = p7, width = 8, height = 6, dpi = 300)
```



## Plink/1.90 and PCA before LD Pruning as a Guideline
Plink is used here to convert the VCF data into formats that can be used for for instance PCA, following [speciationgenomics.github](https://speciationgenomics.github.io/variant_calling/). 
Note that we agreed on using plink/1.90 since the newer version would not have worked on this small dataset.
We tried to run the PCA on the individual chromosomes at first but since our goal is the whole genome (EXCLUDING the sex chromosomes!) it is necessary to concatenate the autosomes first.

Here a text file with the names of the filtered VCF files is needed. It was necessary to to gunzip *.vcf.gz after copying the original files to a different folder and then:
```
for f in *.vcf; do
    bgzip "$f"
done
```
BCFtools need a *.bzip file, hence a new index is needed. I renamed one file for H. semele since a "," breaks the concat command. this was done with `mv "Semele03_OW121722.1,CAKOFV020000007.1.recode.vcf.gz" "Semele03_OW121722.1_CAKOFV020000007.1.recode.vcf.gz"` with this the naming from the text file has to be changed accordingly to the rename. 
```
module load bcftools
cd ~/SNP/bcf/semele_filtered/bzip

for f in *.vcf.gz; do
    bcftools index --tbi "$f"
done
```
now the following code should work for H. semele:
```
module load bcftools
bcftools concat \
--allow-overlaps \
--file-list /storage/frontends/home/username/SNP/bcf/semele_filtered/autosomes_semele.txt \
-O z -o /storage/frontends/home/username/SNP/bcf/semele_autosomes/SNPsemele_autosomes.vcf.gz
```
and for H. statilinus:
```
bcftools concat \
--allow-overlaps \
--file-list /storage/frontends/home/username/SNP/bcf/statilinus_filtered/autosomes_statilinus.txt \
-O z -o /storage/frontends/home/username/SNP/bcf/statilinus_autosomes/SNPstatilinus_autosomes.vcf.gz
```
### Plink
```
#!/bin/bash
#PBS -N PlinkSemeleAuto
#PBS -l select=1:ncpus=1:mem=15gb:scratch_local=100gb
#PBS -l walltime=00:10:00
#PBS -e /storage/frontends/home/username/PBS/
#PBS -o /storage/frontends/home/username/PBS/

module load plink/1.90

in="/storage/frontends/home/username/SNP/bcf/semele_autosomes"
out="/storage/frontends/home/username/SNP/bcf/semele_autosomes"
VCF="${in}/SNPsemele_autosomes.vcf.gz"


# perform linkage pruning - i.e. identify prune sites
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out ${out}/SNPsemele_autosomes_plink


# prune and create pca
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ${out}/SNPsemele_autosomes_plink.prune.in \
--make-bed --pca --out ${out}/SNPsemele_autosomes_plink

# Clean scratch
trap 'clean_scratch' TERM EXIT
```


### PCA
This code for visualization is computed in R after downloading the *.eigenval and *.eigenvec files from MetaCentrum.
So now the id+loc were renamed according to country. A textfile was created, where old sample id's are replaced with new names (only number and year).
```
setwd("C:/Users/username1/Documents/BT_R/semele_autosomes")

library(tidyverse)
library(ggrepel)

# ===============================
# USER INPUTS – EDIT THESE
# ===============================
base_dir       <- "."  # current working directory
eigenvec_file  <- "SNPsemele_autosomes_plink.eigenvec"  # PLINK concatenated eigenvec
eigenval_file  <- "SNPsemele_autosomes_plink.eigenval"  # PLINK concatenated eigenval
meta_file      <- "semele_id+loc.txt"                    # sample ↔ locality
rename_file    <- "sample_rename.txt"                        # old_name → new_name for plotting
out_dir        <- "."                                        # output PDF folder
plot_title     <- "Semele PCA - concatenated autosomes"  # plot title

# ===============================
# read metadata
# ===============================
metadata <- read_tsv(file.path(base_dir, meta_file), show_col_types = FALSE) %>%
  rename(sample = id, locality = loc)

# read sample rename table
rename_table <- read_tsv(file.path(base_dir, rename_file), show_col_types = FALSE) %>%
  rename(old_name = old_name, new_name = new_name)

# ===============================
# read PCA data
# ===============================
pca <- read_table(file.path(base_dir, eigenvec_file), col_names = FALSE)
eigenval <- scan(file.path(base_dir, eigenval_file))

# remove FID column (column 1)
pca <- pca[, -1]

# rename columns
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# clean sample names
pca$sample <- basename(pca$ind)
pca$sample <- sub("\\.final\\.bam$", "", pca$sample)

# ===============================
# add metadata and rename info
# ===============================
pca <- left_join(pca, metadata, by = "sample")
pca <- left_join(pca, rename_table, by = c("sample" = "old_name"))

# create plot names: use new_name if available, otherwise keep original sample
pca$plot_name <- ifelse(is.na(pca$new_name), pca$sample, pca$new_name)

# warn if any metadata missing
if (any(is.na(pca$locality))) {
  warning("Some samples missing locality info in metadata")
}

# ===============================
# percentage variance explained
# ===============================
pve <- eigenval / sum(eigenval) * 100

# ===============================
# create PCA plot
# ===============================
p <- ggplot(pca, aes(PC1, PC2, color = locality)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = plot_name), size = 3) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve[2], 3), "%)")) +
  ggtitle(plot_title)

# ===============================
# save plot
# ===============================
ggsave(
  filename = file.path(out_dir, paste0("PCA_", plot_title, ".pdf")),
  plot = p,
  width = 8,
  height = 6
)

message("PCA plot saved to: ", file.path(out_dir, paste0("PCA_", plot_title, ".pdf")))
```


## Optional: Barcoding
Since the PCA revealed to outliers in H. semele there was the suspicion that they could belong to another species.
The pipeline to extract the mitochondrial COI barcode region for species identification was guided by AI.

I identify the COI coordinates in the mitochondrial genome since the mitogenome is unannotated, first i find where COI is located. The fastest way is to BLAST a known COI sequence from a related species against the mitochondrial genome.
```
module load blast
makeblastdb -in mitochondrial_DNA_semele.fasta -dbtype nucl

blastn \
 -query COI_semele.fasta \
 -db mitochondrial_DNA_semele.fasta \
 -outfmt 6
```
The output will give something like:
```
COI_ref   mito_genome   99.2   658   ...   5321   5978
DQ338868.1      OW121739.2      98.966  1450    15      0       1       1450    4777    3328    0.0     2595
#for statilinus
DQ338596.1      OZ310758.1      94.000  1450    87      0       1       1450    4791    3342    0.0     2217
```
that means: sstart = 4777 send = 3328

and for statilinus sstart = 4791 send = 3342

So COI coordinates are: 4777-3328 (typical barcode region ~658 bp). This only needs to be done once per reference genome.

Now I map WGS reads to the mitochondrial genome (used the mapping code from above, only modified to map to mitochondrial DNA)
```
#index genome
module load bwa-mem2 
bwa-mem2 index mit*.fasta


#script for mapping
#!/bin/bash
#PBS -N MappReadssemele
#PBS -l select=1:ncpus=2:mem=50gb:scratch_local=50gb
#PBS -l walltime=10:00:00
#PBS -e /storage/frontends/home/username/PBS/
#PBS -o /storage/frontends/home/username/PBS/
#PBS -J 0-45

# Input data directory
DATADIR="/storage/frontends/home/username/allign_semele/clean_files"

# Defining the file to process using the Array ID
file=$(ls ${DATADIR}/*_1R.fastq.gz | sed -n "${PBS_ARRAY_INDEX}p")

# Load required programs
module load bwa-mem2
module load samtools

# Output directory
output="/storage/frontends/home/username/semele_mitochondria"

# Reference genome (already indexed)
ref="/storage/frontends/home/username/reference_genomes/weird_semele_ref/mitochondrial_DNA_semele.fasta"

# Aligning individual sequences to the reference
f=${file/_1R.fastq.gz/}
name=$(basename ${f})
reverse=${file/_1R.fastq.gz/_2R.fastq.gz}

bwa-mem2 mem -t 2 -M -a ${ref} ${file} ${reverse} > ${output}/${name}.sam

# Convert sam to bam
samtools view -bS ${output}/${name}.sam -o ${output}/${name}.bam
rm ${output}/${name}.sam

# Sort the bam file
samtools sort ${output}/${name}.bam -o ${output}/${name}.sorted.bam
rm ${output}/${name}.bam

# Clean scratch after the end
trap 'clean_scratch' TERM EXIT
```
Index:
```
module load samtools
samtools index *.sorted.bam
```
Because mitochondrial DNA has very high coverage in WGS, this mapping step is usually fast.
Extract reads overlapping the COI region using the coordinates from BLAST:
```
module load samtools
module load bcftools
module load seqtk

for bam in *.sorted.bam; do
    name=$(basename "$bam" .sorted.bam)
    samtools view -b "$bam" OW121739.2:3328-4777 > "${name}_COI.bam"
	#Convert to FASTQ if needed:
    samtools fastq "${name}_COI.bam" > "${name}_COI.fastq"


	#5. (Optional but recommended) Generate a COI consensus
	#This gives you a barcode sequence for each sample.

	samtools mpileup -uf mitochondrial_DNA_semele.fasta "${name}_COI.bam" | \
    bcftools call -c | \
    vcfutils.pl vcf2fq > "${name}_COI_consensus.fastq"
	#Convert to FASTA:

    seqtk seq -a "${name}_COI_consensus.fastq" > "${name}_COI_consensus.fasta"
	done
```
Since the outlier samples did fit the barcode sequence, they are likely H. semele samples. We then tried to find out if H. statilinus would map to H. semele and also show a barcode.
```
#!/bin/bash
#PBS -N MappReadssemelestatilinus
#PBS -l select=1:ncpus=2:mem=50gb:scratch_local=50gb
#PBS -l walltime=03:00:00
#PBS -e /storage/frontends/home/username/PBS/
#PBS -o /storage/frontends/home/username/PBS/
#PBS -J 0-4

# Input data directory
DATADIR="/storage/frontends/home/username/semele_mitochondria/statilinus_try_neg"

# Defining the file to process using the Array ID
file=$(ls ${DATADIR}/*_1R.fastq.gz | sed -n "${PBS_ARRAY_INDEX}p")

# Load required programs
module load bwa-mem2
module load samtools

# Output directory
output="/storage/frontends/home/username/semele_mitochondria/indexedstat"

# Reference genome (already indexed)
ref="/storage/frontends/home/username/reference_genomes/weird_semele_ref/mitochondrial_DNA_semele.fasta"

# Aligning individual sequences to the reference
f=${file/_1R.fastq.gz/}
name=$(basename ${f})
reverse=${file/_1R.fastq.gz/_2R.fastq.gz}

bwa-mem2 mem -t 2 -M -a ${ref} ${file} ${reverse} > ${output}/${name}.sam

# Convert sam to bam
samtools view -bS ${output}/${name}.sam -o ${output}/${name}.bam
rm ${output}/${name}.sam

# Sort the bam file
samtools sort ${output}/${name}.bam -o ${output}/${name}.sorted.bam
rm ${output}/${name}.bam

# Clean scratch after the end
trap 'clean_scratch' TERM EXIT

after submitting themapping step redo the search for the barcode sequence

module load samtools
samtools index *.sorted.bam

module load bcftools
module load seqtk


for bam in *.sorted.bam; do
    name=$(basename "$bam" .sorted.bam)
    samtools view -b "$bam" OW121739.2:3328-4777 > "${name}_COI.bam"
    samtools fastq "${name}_COI.bam" > "${name}_COI.fastq"

    # Generate consensus
    samtools mpileup -uf mitochondrial_DNA_semele.fasta "${name}_COI.bam" | \
    bcftools call -c | \
    vcfutils.pl vcf2fq > "${name}_COI_consensus.fastq"
    seqtk seq -a "${name}_COI_consensus.fastq" > "${name}_COI_consensus.fasta"
	done
```
The mapping of H. statilinus on H. semele did not give any output, since we now want to try if the mapping of the COI works for H. statilinus on H. statilinus the same procedure was repeated.
```
#script for mapping, everything else is basically identical to the previous steps
#!/bin/bash
#PBS -N MappReadsMTstat
#PBS -l select=1:ncpus=2:mem=50gb:scratch_local=50gb
#PBS -l walltime=04:00:00
#PBS -e /storage/frontends/home/username/PBS/
#PBS -o /storage/frontends/home/username/PBS/
#PBS -J 0-4

# Input data directory
DATADIR="/storage/frontends/home/username/semele_mitochondria/statilinus_try_neg"

# Defining the file to process using the Array ID
file=$(ls ${DATADIR}/*_1R.fastq.gz | sed -n "${PBS_ARRAY_INDEX}p")

# Load required programs
module load bwa-mem2
module load samtools

# Output directory
output="/storage/frontends/home/username/semele_mitochondria/statilinus_ref_COI"

# Reference genome (already indexed)
ref="/storage/frontends/home/username/semele_mitochondria/statilinus_ref_COI/MT_statilinus.fasta"

# Aligning individual sequences to the reference
f=${file/_1R.fastq.gz/}
name=$(basename ${f})
reverse=${file/_1R.fastq.gz/_2R.fastq.gz}

bwa-mem2 mem -t 2 -M -a ${ref} ${file} ${reverse} > ${output}/${name}.sam

# Convert sam to bam
samtools view -bS ${output}/${name}.sam -o ${output}/${name}.bam
rm ${output}/${name}.sam

# Sort the bam file
samtools sort ${output}/${name}.bam -o ${output}/${name}.sorted.bam
rm ${output}/${name}.bam

samtools index ${output}/${name}.sorted.bam

# Clean scratch after the end
trap 'clean_scratch' TERM EXIT

and here an aditional debug with help of ai

for bam in *.sorted.bam; do
    name=$(basename "$bam" .sorted.bam)
    
    # Extract COI region
    samtools view -b "$bam" OZ310758.1:3342-4791 > "${name}_COI.bam"
    
    # Skip empty BAMs
    mapped=$(samtools view -c -F 0x4 "${name}_COI.bam")
    if [ "$mapped" -eq 0 ]; then
        echo "No reads for $name, skipping"
        continue
    fi

    # Optional FASTQ
    samtools fastq "${name}_COI.bam" > "${name}_COI.fastq"

    # Generate consensus using the SAME reference used for mapping
    samtools mpileup -uf MT_statilinus.fasta "${name}_COI.bam" | \
    bcftools call -c | \
    vcfutils.pl vcf2fq > "${name}_COI_consensus.fastq"

    seqtk seq -a "${name}_COI_consensus.fastq" > "${name}_COI_consensus.fasta"
done
```

The two oultliers likely are H. semele butterflies while showing good quality. Therefore they are not excluded from further analyses.


## Excluding Files and new Filters
Since the H. statilinus dataset is rather small we decided to modify the filters a bit in order to get less but higher quality SNPs. Also one individual of each species was excluded due to a too high % of missingness (~40% and ~30%).
/storage/frontends/home/username/allign_semele/finalbam/P33009_1059_S273_L006.final.bam   0.3288
/storage/frontends/home/username/allign_statilinus/finalbam/P33009_1044_S258_L006.final.bam   0.4319

And for H. statilinus max depth was set to 40, maf was set to 0.1 (2 individuals, 30 alleles) and the 40% missingness sample was removed.


### Preparation
For the removal of files I first created a text file where I entered e.g. 
`/storage/frontends/home/username/SNP/bcf/statilinus/P33009_1044_S258_L006.final.bam`
to exclude the corresponding samples.

Here it is important that this is not really a filepath anymore, but the name that is stored in the files (that could be confused with filepaths) such as:
/storage/frontends/home/username/allign_statilinus/finalbam/P33009_1044_S258_L006.final.bam

### VCF Filtering
[KAY Lucek Pipeline](https://hackmd.io/UbLfs1AfSTqCDmIOjWhm6Q#GENERATE-VCF-FILES) was used to conduct the VCF filtering. When having the textfile with the samples to exclude from the analyses, the code for VCF filtering was repeated. Here I depicted both codes since there are different filter adjustments. (H. statilinus had adjusted filters and one sample was excluded):
```
#!/bin/bash
#PBS -N Filtering_OnlyVariantsStatilinus2
#PBS -l select=1:ncpus=1:mem=30gb:scratch_local=100gb
#PBS -l walltime=24:00:00
#PBS -e /storage/frontends/home/username/PBS/
#PBS -o /storage/frontends/home/username/PBS/
#PBS -J 1-30

DATADIR=/storage/frontends/home/username/reference_genomes


test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
## clean scratch after the end
trap 'clean_scratch' TERM EXIT
cd $SCRATCHDIR

## Set scaffold to run with PBS_ARRAY_INDEX
chr_file="$DATADIR/textbcf/statilinus.txt"
chr=$(cat $chr_file | sed -n "${PBS_ARRAY_INDEX}p")

in="/storage/frontends/home/username/SNP/bcf/statilinus"
out="/storage/frontends/home/username/refiltered_vcf/statilinus"
##make outdir if not existent!!
#repeats="/storage/brno12-cerit/home/repulles/03.Heliconius/02.1.fastqMapping/00.Hmel2.5_lepbase_repeats/repeats_all.bed"

## modules
module load vcftools

### First filtering keeping only the big scaffolds and with minGQ 20 (problems with it?)
vcftools --gzvcf ${in}/Statilinus_${chr}.vcf.gz \
	--minGQ 20 \
	--minDP 4 --maxDP 40 --max-missing 0.9 \
	--remove /storage/frontends/home/username/refiltered_vcf/remtext/remove_vcf2_statilinus.txt \
	--remove-indels --max-alleles 2 \
	--min-alleles 2 --maf 0.1 \
	--recode --recode-INFO-all \
	--out ${out}/Statilinus03_${chr}
gzip ${out}/Statilinus03_${chr}.recode.vcf #this was added by chatgpt since streaming does not work for vcf


# Clean scratch
trap 'clean_scratch' TERM EXIT
```

H. semele with same filters applied, but with one sample excluded:
```
#!/bin/bash
#PBS -N Filtering_OnlyVariantsSemele2
#PBS -l select=1:ncpus=1:mem=30gb:scratch_local=100gb
#PBS -l walltime=24:00:00
#PBS -e /storage/frontends/home/username/PBS/
#PBS -o /storage/frontends/home/username/PBS/
#PBS -J 1-30

DATADIR=/storage/frontends/home/username/reference_genomes


test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
## clean scratch after the end
trap 'clean_scratch' TERM EXIT
cd $SCRATCHDIR

## Set scaffold to run with PBS_ARRAY_INDEX
chr_file="$DATADIR/textbcf/semele.txt"
chr=$(cat $chr_file | sed -n "${PBS_ARRAY_INDEX}p")

in="/storage/frontends/home/username/SNP/bcf/semele"
out="/storage/frontends/home/username/refiltered_vcf/semele"
##make outdir if not existent!!
#repeats="/storage/brno12-cerit/home/repulles/03.Heliconius/02.1.fastqMapping/00.Hmel2.5_lepbase_repeats/repeats_all.bed"

## modules
module load vcftools

### First filtering keeping only the big scaffolds and with minGQ 20 (problems with it?)
vcftools --gzvcf ${in}/Semele_${chr}.vcf.gz \
	--minGQ 20 \
	--minDP 4 --maxDP 100 --max-missing 0.9 \
	--remove /storage/frontends/home/username/refiltered_vcf/remtext/remove_vcf2_semele.txt \
	--remove-indels --max-alleles 2 \
	--min-alleles 2 --maf 0.03 \
	--recode --recode-INFO-all \
	--out ${out}/Semele03_${chr}
gzip ${out}/Semele03_${chr}.recode.vcf 

# Clean scratch
trap 'clean_scratch' TERM EXIT
```

With the line below one can check the samples that are stored here. Do this in the console beforehand to make sure the samples are really excluded.
```
zgrep "^#CHROM" Statilinus_${chr}.vcf.gz
```

### VCF Stats
This step is based on [speciationgenomics.github](https://speciationgenomics.github.io/variant_calling/) and [KAY Lucek Pipeline](https://hackmd.io/UbLfs1AfSTqCDmIOjWhm6Q#GENERATE-VCF-FILES). To have an insight into the quality of the SNPs the VCF stats are run in MetaCentrum:
```
#run the stats on the frontend, not too demanding only a few seconds

DATADIR="/storage/frontends/home/username/refiltered_vcf/autosomes_statilinus"
OUT="/storage/frontends/home/username/refiltered_vcf/stats_vcf_statilinus"
VCF_FILE="${DATADIR}/SNPstatilinus_autosomes.vcf.gz"

module load vcftools

#calculate allele frequency:
vcftools --gzvcf $VCF_FILE --freq2 --out $OUT --max-alleles 2

#calculate mean depth per individual:
vcftools --gzvcf $VCF_FILE --depth --out $OUT

#calculate mean depth per site:
vcftools --gzvcf $VCF_FILE --site-mean-depth --out $OUT

#calculate site quality:
vcftools --gzvcf $VCF_FILE --site-quality --out $OUT

#calculate proportion of missing data per individual:
vcftools --gzvcf $VCF_FILE --missing-indv --out $OUT

#calculate proportion of missing data per site:
vcftools --gzvcf $VCF_FILE --missing-site --out $OUT


#Calculate heterozygosity and inbreeding coefficient per individual:
vcftools --gzvcf $VCF_FILE --het --out $OUT
```
and here the visualization of the Stats using R and AI:
```
# ===============================
# VCFtools QC plotting script
# ===============================

# 1) Load libraries
library(tidyverse)

# 2) Set working directory to where your VCFtools outputs are
setwd("C:/Users/username1/Documents/BT_R/vcf_stats2/stats_vcf_semele/")


# 3) Read VCFtools output files

# Site quality
var_qual <- read_delim("stats_vcf_semele.lqual", delim = "\t",
                       skip = 1,
                       col_names = c("chr","pos","qual"))

# Mean depth per site
site_depth <- read_delim("stats_vcf_semele.ldepth.mean", delim = "\t",
                         skip = 1,
                         col_names = c("chr","pos","mean_depth"))

# Depth per individual
ind_depth <- read_delim("stats_vcf_semele.idepth", delim = "\t",
                        skip = 1,
                        col_names = c("sample","N_called","mean_depth","sum_depth"))

# Missing per individual
ind_missing <- read_delim("stats_vcf_semele.imiss", delim = "\t",
                          skip = 1,
                          col_names = c("sample","N_called","N_missing","frac_missing"))


var_miss <- read_delim("stats_vcf_semele.lmiss", delim = "\t",
                       
                       col_names = c("chr","pos","nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)


ind_het <- read_delim("stats_vcf_semele.het", delim = "\t",
                      skip = 1,
                      col_names = c("ind","ho","he","nsites","f"))

var_freq <- read_delim("stats_vcf_semele.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))



# 4) Create density plots

p1 <- ggplot(var_qual, aes(qual)) +
  geom_density(fill = "dodgerblue1", alpha = 0.5) +
  xlim(0, 100) +
labs(title = "Site Quality Distribution") +
  theme_light()

p2 <- ggplot(site_depth, aes(mean_depth)) +
  geom_density(fill = "forestgreen", alpha = 0.5) +
  labs(title = "Mean Depth per Site") +
  theme_light()

p3 <- ggplot(ind_depth, aes(mean_depth)) +
  geom_density(fill = "orange", alpha = 0.5) +
  labs(title = "Mean Depth per Individual") +
  theme_light()

p4 <- ggplot(ind_missing, aes(frac_missing)) +
  geom_density(fill = "firebrick", alpha = 0.5) +
  labs(title = "Missing Data per Individual") +
  theme_light()

p55 <- ggplot(var_miss, aes(fmiss)) +
  geom_density(fill = "purple", alpha = 0.5) +
  labs(title = "Missing Data per Site") +
  theme_light()

p66 <- ggplot(ind_het, aes(f)) +
  geom_histogram(fill = "darkcyan", alpha = 0.5) +
  labs(title = "Inbreeding Coefficient (F)") +
  theme_light()

p7 <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  labs(title = "MAF") 
p7 + theme_light()

# 5) Display all plots
print(p1)
print(p2)
print(p3)
print(p4)
print(p55)
print(p66)
print(p7)

# 6) Optional: save all plots to a PDF
pdf("vcftools_QC_summary_plots.pdf", width = 8, height = 6)
print(p1)
print(p2)
print(p3)
print(p4)
print(p55)
print(p66)
print(p7)
dev.off()
```
A summary table can be created manually by importing the csv files to R and taking min(), max(), and mean() of the corresponding columns.



### Concatenate Autosomes
Since the files need to be zipped in another way, the zipping had to be repeated (the following steps are all computed using the console, but navigate into the right folder first:
H. statilinus:
```
gunzip *.vcf.gz
module load samtools
for f in *.vcf; do
    bgzip "$f"
done
```
H. semele (since the file had to be renamed (see above)):
```
mv "Semele03_OW121722.1,CAKOFV020000007.1.recode.vcf.gz" "Semele03_OW121722.1_CAKOFV020000007.1.recode.vcf.gz"
gunzip *.vcf.gz
module load samtools
for f in *.vcf; do
    bgzip "$f"
done
```

For both species seperately:
```
module load bcftools

for f in *.vcf.gz; do
    bcftools index --tbi "$f"
done
```

The following code should work for H. semele now, navigation into the folder with the newly created recode.vcf files beforehand is necessary.
```
module load bcftools
bcftools concat \
--allow-overlaps \
--file-list /storage/frontends/home/username/SNP/bcf/semele_filtered/autosomes_semele.txt \
-O z -o /storage/frontends/home/username/refiltered_vcf/autosomes_semele/SNPsemele_autosomes.vcf.gz
```

and for H. statilinus:
```
bcftools concat \
--allow-overlaps \
--file-list /storage/frontends/home/username/SNP/bcf/statilinus_filtered/autosomes_statilinus.txt \
-O z -o /storage/frontends/home/username/refiltered_vcf/autosomes_statilinus/SNPstatilinus_autosomes.vcf.gz
```

Plink, PCA and the Stats were rerun in order to see what changed. But since in the next step LD pruning is applied, a new PCA had to be done. The exact codes are not depicted here, since its repetitive and can be conducted from the codes above/below. 

## LD Pruning and PCA
We now decided on taking maf = 0.1 for H. statilinus and the maf = 0.3 for H. semele with one sample removed in each species (refiltered_vcf in metacentrum). LD pruning was performed using [KAY Lucek Pipeline](https://hackmd.io/UbLfs1AfSTqCDmIOjWhm6Q#GENERATE-VCF-FILES), directly on the frontend.
```
module load plink/1.90
#semele
plink --vcf SNPsemele_autosomes.vcf.gz --allow-extra-chr --double-id --set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --out Variants_AllChr
#statilinus
plink --vcf SNPstatilinus_autosomes.vcf.gz --allow-extra-chr --double-id --set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --out Variants_AllChr
```
Another PCA:
```
module load plink/1.90
#semele
plink --vcf SNPsemele_autosomes.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --extract Variants_AllChr.prune.in --pca --out ../PCA_plink_semele/Variants_AllChr_pruned_semele
#statilinus
plink --vcf SNPstatilinus_autosomes.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --extract Variants_AllChr.prune.in --pca --out ../PCA_plink_statilinus/Variants_AllChr_pruned_statilinus
```
And here the new PCA plot in R:
```
setwd("C:/Users/username1/Documents/BT_R/01_ld_pca/PCA_plink_statilinus")

library(tidyverse)
library(ggrepel)

# ===============================
# USER INPUTS – EDIT THESE
# ===============================
base_dir       <- "."  # current working directory
eigenvec_file  <- "Variants_AllChr_pruned_statilinus.eigenvec"  # PLINK concatenated eigenvec
eigenval_file  <- "Variants_AllChr_pruned_statilinus.eigenval"  # PLINK concatenated eigenval
meta_file      <- "statilinus_id+loc.txt"                    # sample ↔ locality
rename_file    <- "sample_rename.txt"                        # old_name → new_name for plotting
out_dir        <- "."                                        # output PDF folder
plot_title     <- "Statilinus PCA - concatenated autosomes"  # plot title

# ===============================
# read metadata
# ===============================
metadata <- read_tsv(file.path(base_dir, meta_file), show_col_types = FALSE) %>%
  rename(sample = id, locality = loc)

# read sample rename table
rename_table <- read_tsv(file.path(base_dir, rename_file), show_col_types = FALSE) %>%
  rename(old_name = old_name, new_name = new_name)

# ===============================
# read PCA data
# ===============================
pca <- read_table(file.path(base_dir, eigenvec_file), col_names = FALSE)
eigenval <- scan(file.path(base_dir, eigenval_file))

# remove FID column (column 1)
pca <- pca[, -1]

# rename columns
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# clean sample names
pca$sample <- basename(pca$ind)
pca$sample <- sub("\\.final\\.bam$", "", pca$sample)

# ===============================
# add metadata and rename info
# ===============================
pca <- left_join(pca, metadata, by = "sample")
pca <- left_join(pca, rename_table, by = c("sample" = "old_name"))

# create plot names: use new_name if available, otherwise keep original sample
pca$plot_name <- ifelse(is.na(pca$new_name), pca$sample, pca$new_name)

# warn if any metadata missing
if (any(is.na(pca$locality))) {
  warning("Some samples missing locality info in metadata")
}

# ===============================
# percentage variance explained
# ===============================
pve <- eigenval / sum(eigenval) * 100

# ===============================
# create PCA plot
# ===============================
p <- ggplot(pca, aes(PC1, PC2, color = locality)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = plot_name), size = 3) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve[2], 3), "%)")) +
  ggtitle(plot_title)

# ===============================
# save plot
# ===============================
ggsave(
  filename = file.path(out_dir, paste0("PCA_", plot_title, ".png")),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

message("PCA plot saved to: ", file.path(out_dir, paste0("PCA_", plot_title, ".pdf")))
```

## Admixture
Here mostly [KAY Lucek Pipeline](https://hackmd.io/UbLfs1AfSTqCDmIOjWhm6Q#GENERATE-VCF-FILES), [Speciation & Population Genomics: a how-to-guide](https://speciationgenomics.github.io/) were used, and some code snippets were received from Mar personally.

First input files are created with:
```
plink --vcf SNPsemele_autosomes.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract /storage/frontends/home/username/01final_files_vcf/autosomes_semele/Variants_AllChr.prune.in --make-bed --out Variants_AllChr_LDpruned_admix_semele
```
```
plink --vcf SNPstatilinus_autosomes.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract /storage/frontends/home/username/01final_files_vcf/autosomes_statilinus/Variants_AllChr.prune.in --make-bed --out Variants_AllChr_LDpruned_admix_statilinus
```
Then the new files are taken to avoid chromosome errors.
```
FILE=Variants_AllChr_LDpruned_admix_semele

awk '{$1="0";print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim $FILE.bim.originalChrnames
mv $FILE.bim.tmp $FILE.bim
```
The different k cluster values and cv errors were extracted to achieve the best value of k clusters with low cv error by looping over several k.
```
module load admixture
for K in {1..6}; do
  admixture --cv $FILE.bed $K > log${K}.out
done
```
To make plotting easier, a file with individual names and species names in seperate columns was created.
```
grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > $FILE.cv.error
```
Extract the species name from the individual name:
```
awk '{split($1,name,"."); print $1,name[2]}' $FILE.nosex > $FILE.list
```

Using this files admixture can be plotted in R as followed, considering the id+loc.txt and sample_rename.txt tsv-files are prepared already.
```
############
# Plot admixt results following https://github.com/laurabenestan/Admixture
###########

#more like classification/summary than a classic ADMIXTURE visualization
#code from Mar; sent to Pavel



############
# ADMIXTURE plotting script (fixed)
############

setwd("C:/Users/username1/Documents/BT_R/01_admixture/semele/variants_admixture")

library(stringr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)

############################
# 1. LOAD METADATA
############################
meta <- read.delim("semele_id+loc.txt", header=TRUE, stringsAsFactors=FALSE)
colnames(meta) <- c("id","loc")
meta$spp <- "semele"

name_map <- read.delim("sample_rename.txt", header=TRUE, stringsAsFactors=FALSE)

############################
# 2. LOOP OVER K = 1–6
############################

for (k in 2:6) {
  
  # read Q file safely
  admixture <- read.table(
    paste0("Variants_AllChr_LDpruned_admix_semele.", k, ".Q"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  admixture <- read.table(
    paste0("Variants_AllChr_LDpruned_admix_semele.", k, ".Q"),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  id <- read.table(
    "Variants_AllChr_LDpruned_admix_semele.list",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  id <- id[,1]
  
  admixture$IND <- id
  
  colnames(admixture) <- c(paste0("K", 1:k), "IND")
  admixture <- admixture %>% relocate(IND)
  
  # name columns
  colnames(admixture) <- c("IND", paste0("K", 1:k))
  
  # clean IDs
  admixture$IND <- str_remove(basename(admixture$IND), "\\.final\\.bam$")
  
  # merge metadata
  admixture <- admixture %>%
    left_join(meta, by = c("IND" = "id")) %>%
    left_join(name_map, by = c("IND" = "old_name"))
  
  # apply rename
  admixture$IND <- ifelse(
    is.na(admixture$new_name),
    admixture$IND,
    admixture$new_name
  )
  
  ############################
  # ORDER INDIVIDUALS
  ############################
  
  admixture <- admixture %>%
    arrange(loc, IND)
  
  admixture$IND <- factor(admixture$IND, levels = unique(admixture$IND))
  
  ############################
  # reshape
  ############################
  
  # keep only real K columns
  admixture_clean <- admixture %>%
    select(IND, spp, loc, matches("^K[0-9]+$"))
  
  # melt
  admixture_long <- melt(
    admixture_clean,
    id.vars = c("IND","spp","loc"),
    variable.name = "ANCESTRY",
    value.name = "PERC"
  )
  
  # force numeric (CRITICAL FIX)
  admixture_long$PERC <- as.numeric(as.character(admixture_long$PERC))
  
  
  # filtering
  admixture_long_50 <- subset(admixture_long, PERC >= 0.50)
  
  ############################
  # plot
  ############################
  
  p <- ggplot(admixture_long, aes(x = IND, y = PERC, fill = ANCESTRY)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    labs(x = "Individual", y = "Ancestry") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(colour="grey", linetype="dashed"),
      axis.title.x = element_text(size=14, family="Helvetica", face="bold"),
      axis.text.x = element_text(size=8, angle=90, hjust=1),
      axis.title.y = element_text(size=14, family="Helvetica", face="bold"),
      axis.text.y = element_text(size=12)
    )
  
  print(p)
  
  # save
  ggsave(paste0("Admixture_K", k, ".png"), p, width=10, height=5, dpi=300)
}
```
### Plot the Admixture on a map
After the creation of a textfile where date and locations are specified (in tsv format, see scheme below), the following code can be used.
```
old_name	new_name	GPS_N	GPS_E	year
P33009_1015_S229_L006	15_2021_SK	48.5771222112687	17.183917938063807	2021
P33009_1016_S230_L006	16_2021_SK	48.5771222112687	17.183917938063807	2021
P33009_1017_S231_L006	17_2021_SK	48.5771222112687	17.183917938063807	2021
P33009_1043_S257_L006	43_1977_GER	52.34382965425453	13.877711291775757	1977
P33009_1044_S258_L006	44_1977_GER	52.34382965425453	13.877711291775757	1977
P33009_1045_S259_L006	45_2000_GER	52.047881795072705	13.435286724237844	2000
P33009_1046_S260_L006	46_2000_GER	52.047881795072705	13.435286724237844	2000
P33009_1047_S261_L006	47_2003_GER	52.1900000	13.4855556	2003
P33009_1048_S262_L006	48_2018_GER	52.03750	13.74861	2018
P33009_1049_S263_L006	49_2019_GER	52.03750	13.74861	2019
P33009_1050_S264_L006	50_1974_SK	48.5771222112687	17.183917938063807	1974
P33009_1051_S265_L006	51_1973_SK	48.5771222112687	17.183917938063807	1973
P33009_1052_S266_L006	52_1975_SK	48.5771222112687	17.183917938063807	1975
P33009_1053_S267_L006	53_1975_SK	48.5771222112687	17.183917938063807	1975
P33009_1060_S274_L006	60_1983_GER	52.395238608559445	13.08238583240185	1983
```
```
############
# ADMIXTURE PIE MAP – THREE MAPS: ALL / OLD (≤2000) / RECENT (>2000)
############

setwd("C:/Users/username1/Documents/BT_R/01_admixture/semele/variants_admixture")

library(stringr)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggforce)
library(tidyr)

############################
# PARAMETERS
############################

PIE_R         <- 0.15
SPREAD_FACTOR <- 1.6

############################
# LOAD DATA
############################

# Use the new TSV with the year column
meta     <- read.table("old_new_N_E_with_year.tsv", header = TRUE, stringsAsFactors = FALSE)
name_map <- read.table("sample_rename.txt",          header = TRUE, stringsAsFactors = FALSE)

k      <- 2
k_cols <- paste0("K", 1:k)

q   <- read.table(paste0("Variants_AllChr_LDpruned_admix_semele.", k, ".Q"))
q   <- as.data.frame(lapply(q, as.numeric))
colnames(q) <- k_cols

ids <- read.table("Variants_AllChr_LDpruned_admix_semele.list", header = FALSE)
ids <- data.frame(IND = str_remove(basename(ids$V1), "\\.final\\.bam$"))

df_all <- cbind(ids, q) %>%
  left_join(name_map, by = c("IND" = "old_name")) %>%
  left_join(meta,     by = "new_name")

df_all$longitude <- as.numeric(df_all$GPS_E)
df_all$latitude  <- as.numeric(df_all$GPS_N)
df_all$year      <- as.integer(df_all$year)

df_all <- df_all %>% filter(!is.na(longitude), !is.na(latitude), !is.na(year))

############################
# SHARED MAP EXTENT (same for all three maps)
############################

lon_range <- range(df_all$longitude)
lat_range <- range(df_all$latitude)
world     <- ne_countries(scale = "medium", returnclass = "sf")

############################
# HELPER FUNCTION
############################

make_pie_map <- function(df, title_label) {
  
  # circular spread
  df <- df %>%
    group_by(longitude, latitude) %>%
    mutate(
      n_at_site    = n(),
      rank_in_site = row_number(),
      spread_r     = ifelse(n_at_site == 1, 0,
                            PIE_R * SPREAD_FACTOR * (1 + log(n_at_site))),
      angle_j      = ifelse(n_at_site == 1, 0,
                            2 * pi * (rank_in_site - 1) / n_at_site),
      longitude_j  = longitude + spread_r * sin(angle_j),
      latitude_j   = latitude  + spread_r * cos(angle_j)
    ) %>%
    ungroup()
  
  df$id <- seq_len(nrow(df))
  
  coords_lookup <- df %>%
    select(id, longitude, latitude, longitude_j, latitude_j, n_at_site)
  
  df_pie <- df %>%
    select(id, all_of(k_cols)) %>%
    pivot_longer(cols = all_of(k_cols), names_to = "cluster", values_to = "value") %>%
    group_by(id) %>%
    arrange(cluster, .by_group = TRUE) %>%
    mutate(
      end   = cumsum(value),
      start = lag(end, default = 0)
    ) %>%
    ungroup() %>%
    left_join(coords_lookup, by = "id")
  
  leaders <- coords_lookup %>%
    filter(n_at_site > 1) %>%
    distinct()
  
  ggplot() +
    
    geom_sf(data = world, fill = "#c7d59f", color = "#6b6b4b", linewidth = 0.2) +
    
    geom_segment(
      data = leaders,
      aes(x = longitude, y = latitude, xend = longitude_j, yend = latitude_j),
      color = "grey25", linewidth = 0.35
    ) +
    
    geom_point(
      data  = leaders %>% distinct(longitude, latitude),
      aes(x = longitude, y = latitude),
      shape = 21, fill = "white", color = "grey20", size = 1.5, stroke = 0.6
    ) +
    
    geom_arc_bar(
      data = df_pie,
      aes(x0 = longitude_j, y0 = latitude_j,
          r0 = 0, r = PIE_R,
          start = 2 * pi * start, end = 2 * pi * end,
          fill = cluster),
      color = "black", linewidth = 0.25
    ) +
    
    scale_fill_manual(values = c(K1 = "#2c7bb6", K2 = "#d7191c")) +
    
    coord_sf(
      xlim   = c(lon_range[1] - 1.5, lon_range[2] + 1.5),
      ylim   = c(lat_range[1] - 1.5, lat_range[2] + 1.5),
      expand = FALSE
    ) +
    
    labs(fill = "cluster", x = NULL, y = NULL, title = title_label) +
    
    theme_minimal() +
    theme(
      panel.background  = element_rect(fill = "#dbe7c9", color = NA),
      panel.grid        = element_blank(),
      plot.title        = element_text(face = "bold", size = 13, hjust = 0.5,
                                       margin = margin(b = 6)),
      legend.position   = c(0.92, 0.20),
      legend.background = element_rect(fill = "white", color = "grey70",
                                       linewidth = 0.4),
      legend.margin     = margin(6, 8, 6, 8)
    )
}

############################
# GENERATE THREE MAPS
############################

p_all    <- make_pie_map(df_all,                          "All samples")
p_old    <- make_pie_map(df_all %>% filter(year <= 2000), "Historic samples (≤ 2000)")
p_recent <- make_pie_map(df_all %>% filter(year >  2000), "Recent samples (> 2000)")

############################
# SAVE
############################

ggsave("ADMIXTURE_PIE_all.png",    p_all,    width = 9, height = 11, dpi = 300)
ggsave("ADMIXTURE_PIE_old.png",    p_old,    width = 9, height = 11, dpi = 300)
ggsave("ADMIXTURE_PIE_recent.png", p_recent, width = 9, height = 11, dpi = 300)
```


## Pixy Nucleotide Diversity (pi), Dxy, Fst, ...
Since for this step another VCF filtered file is needed, Invariants and allsites had to be summarized prior to the analysis. The next steps are mostly based on [KAY Lucek Pipeline](https://hackmd.io/UbLfs1AfSTqCDmIOjWhm6Q#GENERATE-VCF-FILES).

### Filter Invariants
This code is similar to previous ones, only thing that is changed now is that *max-maf 0* is used
```
#!/bin/bash
#PBS -N Filtering_OnlyInvariants_Semele
#PBS -l select=1:ncpus=1:mem=30gb:scratch_local=100gb
#PBS -l walltime=24:00:00
#PBS -e /storage/frontends/home/username/PBS/
#PBS -o /storage/frontends/home/username/PBS/
#PBS -J 1-30

DATADIR=/storage/frontends/home/username/reference_genomes

test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
## clean scratch after the end
trap 'clean_scratch' TERM EXIT
cd $SCRATCHDIR

## Set scaffold to run with PBS_ARRAY_INDEX
chr_file="$DATADIR/textbcf/semele.txt"
chr=$(cat $chr_file | sed -n "${PBS_ARRAY_INDEX}p")

## files
in="/storage/frontends/home/username/SNP/bcf/semele"
out="/storage/frontends/home/username/invariants_for_pixy/semele"
#repeats="/storage/brno12-cerit/home/repulles/03.Heliconius/02.1.fastqMapping/00.Hmel2.5_lepbase_repeats/repeats_all.bed"

## modules
module load vcftools

### First filtering keeping only the big scaffolds and with minGQ 20 (problems with it?)
vcftools --gzvcf ${in}/Semele_${chr}.vcf.gz \
	--minGQ 20 \
	--minDP 4 --maxDP 100 --max-missing 0.9 \
	--remove /storage/frontends/home/username/refiltered_vcf/remtext/remove_vcf2_semele.txt \
	--remove-indels --max-alleles 2 \
	--max-maf 0 \
	--recode --recode-INFO-all \
	--out Semele_Invariants_${chr}

## Compress and index
module load bcftools
bcftools view -O z -o ${out}/Semele_Invariants_${chr}.vcf.gz Semele_Invariants_${chr}.recode.vcf
bcftools index -t ${out}/Semele_Invariants_${chr}.vcf.gz
```

### Concatenate Chromosomes
The files need to be zipped in another way the zipping is redone in the frontend, navigate in the right folder first: 
```
mv "Semele_Invariants_OW121722.1,CAKOFV020000007.1.vcf" "Semele_Invariants_OW121722.1_CAKOFV020000007.1.vcf" #thats just renaming since the "," breaks the command
gunzip *.vcf.gz
module load samtools
for f in *.vcf; do
    bgzip "$f"
done
```
Then another index:
```
module load bcftools

for f in *.vcf.gz; do
    bcftools index --tbi "$f"
done
```
Then concatenate again:
```
module load bcftools
bcftools concat \
--allow-overlaps \
--file-list semele_chr.txt \
-O z -o Invariants_semele_autosomes.vcf.gz
```
In order to extract the file names of the VCF files this command can be used. It is also helpful for the poplulation file that will be used soon.
```
bcftools query -l ${file}.vcf.gz > statilinus_pop.txt
```

So the text file that is used for the later R script and so on looks something like this after adding the countries:
```
/storage/frontends/home/username/allign_semele/finalbam/P33009_1001_S215_L006.final.bam	Slovakia
/storage/frontends/home/username/allign_semele/finalbam/P33009_1002_S216_L006.final.bam	Slovakia
/storage/frontends/home/username/allign_semele/finalbam/P33009_1003_S217_L006.final.bam	Slovakia
/storage/frontends/home/username/allign_semele/finalbam/P33009_1004_S218_L006.final.bam	Slovakia
/storage/frontends/home/username/allign_semele/finalbam/P33009_1005_S219_L006.final.bam	Czech Republic
/storage/frontends/home/username/allign_semele/finalbam/P33009_1006_S220_L006.final.bam	Austria
/storage/frontends/home/username/allign_semele/finalbam/P33009_1007_S221_L006.final.bam	Austria
/storage/frontends/home/username/allign_semele/finalbam/P33009_1008_S222_L006.final.bam	Austria
/storage/frontends/home/username/allign_semele/finalbam/P33009_1009_S223_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1010_S224_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1011_S225_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1012_S226_L006.final.bam	Germany
/storage/frontends/home/username/allign_semele/finalbam/P33009_1013_S227_L006.final.bam	Germany
/storage/frontends/home/username/allign_semele/finalbam/P33009_1014_S228_L006.final.bam	Germany
/storage/frontends/home/username/allign_semele/finalbam/P33009_1018_S232_L006.final.bam	Slovakia
/storage/frontends/home/username/allign_semele/finalbam/P33009_1019_S233_L006.final.bam	Slovakia
/storage/frontends/home/username/allign_semele/finalbam/P33009_1020_S234_L006.final.bam	Czech Republic
/storage/frontends/home/username/allign_semele/finalbam/P33009_1021_S235_L006.final.bam	Czech Republic
/storage/frontends/home/username/allign_semele/finalbam/P33009_1022_S236_L006.final.bam	Czech Republic
/storage/frontends/home/username/allign_semele/finalbam/P33009_1023_S237_L006.final.bam	Czech Republic
/storage/frontends/home/username/allign_semele/finalbam/P33009_1024_S238_L006.final.bam	Czech Republic
/storage/frontends/home/username/allign_semele/finalbam/P33009_1025_S239_L006.final.bam	Czech Republic
/storage/frontends/home/username/allign_semele/finalbam/P33009_1026_S240_L006.final.bam	Czech Republic
/storage/frontends/home/username/allign_semele/finalbam/P33009_1027_S241_L006.final.bam	Czech Republic
/storage/frontends/home/username/allign_semele/finalbam/P33009_1028_S242_L006.final.bam	Austria
/storage/frontends/home/username/allign_semele/finalbam/P33009_1029_S243_L006.final.bam	Austria
/storage/frontends/home/username/allign_semele/finalbam/P33009_1030_S244_L006.final.bam	Austria
/storage/frontends/home/username/allign_semele/finalbam/P33009_1031_S245_L006.final.bam	Austria
/storage/frontends/home/username/allign_semele/finalbam/P33009_1032_S246_L006.final.bam	Austria
/storage/frontends/home/username/allign_semele/finalbam/P33009_1033_S247_L006.final.bam	Austria
/storage/frontends/home/username/allign_semele/finalbam/P33009_1034_S248_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1035_S249_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1036_S250_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1037_S251_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1038_S252_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1039_S253_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1040_S254_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1041_S255_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1042_S256_L006.final.bam	Hungary
/storage/frontends/home/username/allign_semele/finalbam/P33009_1054_S268_L006.final.bam	Germany
/storage/frontends/home/username/allign_semele/finalbam/P33009_1055_S269_L006.final.bam	Germany
/storage/frontends/home/username/allign_semele/finalbam/P33009_1056_S270_L006.final.bam	Germany
/storage/frontends/home/username/allign_semele/finalbam/P33009_1057_S271_L006.final.bam	Germany
/storage/frontends/home/username/allign_semele/finalbam/P33009_1058_S272_L006.final.bam	Germany
```

### Combine Invariants and Variants
Firstly an index is needed:
```
bcftools index --tbi SNPsemele_autosomes.vcf.gz
bcftools index --tbi Invariants_semele_autosomes.vcf.gz
```
After this the files are concatenated:
```
bcftools concat --allow-overlaps -O z -o allsites_autosomes_semele.vcf.gz SNPsemele_autosomes.vcf.gz Invariants_semele_autosomes.vcf.gz
```

### Installation of Pixy
Since the download did not work on the tarkil node we switched to another frontend (...@zenith.cerit-sc.cz) and with these commands from Pedro, the installation worked. But first navigate inside the folder where you want the installation to be:

```
module add mambaforge
mamba create -n pixy -c conda-forge --override-channels pixy
#write "Y"
mamba activate pixy
mamba install -c bioconda -c conda-forge htslib samtools=1.21
#write "Y"
pixy --help
```
use `conda info --envs` to find the path where pixy is installed and use this in your script

/storage/brno12-cerit/home/username/.conda/envs/pixy

### Run Pixy
Now pixy is ready to run and for this a text file with the names of the samples as in the VCF files and tab separated locality (or anything to keep the populations apart) is needed and could look something like this:

```
/storage/frontends/home/username/allign_semele/finalbam/P33009_1001_S215_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1002_S216_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1003_S217_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1004_S218_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1005_S219_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1006_S220_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1007_S221_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1008_S222_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1009_S223_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1010_S224_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1011_S225_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1012_S226_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1013_S227_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1014_S228_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1018_S232_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1019_S233_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1020_S234_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1021_S235_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1022_S236_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1023_S237_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1024_S238_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1025_S239_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1026_S240_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1027_S241_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1028_S242_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1029_S243_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1030_S244_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1031_S245_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1032_S246_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1033_S247_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1034_S248_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1035_S249_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1036_S250_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1037_S251_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1038_S252_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1039_S253_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1040_S254_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1041_S255_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1042_S256_L006.final.bam	AUT_HUN_SK
/storage/frontends/home/username/allign_semele/finalbam/P33009_1054_S268_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1055_S269_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1056_S270_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1057_S271_L006.final.bam	CZ_GER
/storage/frontends/home/username/allign_semele/finalbam/P33009_1058_S272_L006.final.bam	CZ_GER
```

And this is the code for the analysis in metacentrum:
```
#!/bin/bash
#PBS -N pixy
#PBS -l select=1:ncpus=3:mem=15gb:scratch_local=100gb
#PBS -l walltime=02:00:00
#PBS -e /storage/frontends/home/username/PBS/
#PBS -o /storage/frontends/home/username/PBS/

DATADIR=/storage/brno12-cerit/home/username/files_pixy


test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
## clean scratch after the end
trap 'clean_scratch' TERM EXIT
cd $SCRATCHDIR

### Run pixy for pi and dxy

## load pixy
module add mambaforge
source activate /storage/brno12-cerit/home/username/.conda/envs/pixy

## files
cp $DATADIR/all_files_semele/allsites_autosomes_semele.vcf.gz* .
VCF=allsites_autosomes_semele.vcf.gz

## compress and index your vcf (if you didn't)
#bgzip [your.file.vcf]
#tabix [your.file.vcf.gz]

mkdir -p pixy_1_semele

## run pixy by chr
pixy --stats pi fst dxy \
--vcf $VCF \
--populations $DATADIR/pop_semele1.txt \
--window_size 10000 \
--n_cores 3 \
--output_folder pixy_1_semele

cp -r pixy_1_semele $DATADIR
```

And then another pixy was run with the same code as above. Only now "old" and "recent" were added to the second column of the tsv file to obtain 4 populations.


### Visualize pixy results and VCF-based metrics in R with se as error bars
For visualizing the metrics obtained above the AI generated code ([Claude](https://claude.ai/)) was used:
```
################################################################################
# Hipparchia semele Population Genomics: Pixy + VCF-based Diversity Analysis
# Populations:
#   2-pop: AUT_HUN_SK, CZ_GER
#   4-pop: AUT_HUN_SK_old, AUT_HUN_SK_recent, CZ_GER_old, CZ_GER_recent
#
# NOTE: Fst is intentionally excluded here (handled separately elsewhere).
# Error bars throughout are STANDARD ERROR (SE = SD / sqrt(n)), not SD.
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

# ── File paths ────────────────────────────────────────────────────────────────
dir_2pop <- "C:/Users/username1/Documents/BT_R/pixy/semele/pixy_2pop_semele/"
dir_4pop <- "C:/Users/username1/Documents/BT_R/pixy/semele/pixy_4pop_semele/"
vcf_file <- "C:/Users/username1/Documents/BT_R/pixy/semele/SNPsemele_autosomes.vcf.gz"
pop_2pop <- file.path(dir_2pop, "2pop_semele.txt")
pop_4pop <- file.path(dir_4pop, "4pop_semele.txt")
out_dir  <- "C:/Users/username1/Documents/BT_R/diversity_se/semele/"

# ── Colour palette (locality) ─────────────────────────────────────────────────
loc_colors <- c(
  "AUT_HUN_SK" = "#2196F3",
  "CZ_GER"     = "#E91E63"
)

# ── Helper: add locality + time metadata from pop name ───────────────────────
add_meta <- function(df, pop_col = "pop") {
  df %>%
    mutate(
      locality = case_when(
        grepl("AUT_HUN_SK", .data[[pop_col]]) ~ "AUT_HUN_SK",
        grepl("CZ_GER",     .data[[pop_col]]) ~ "CZ_GER"
      ),
      time = case_when(
        grepl("_old",    .data[[pop_col]]) ~ "Old",
        grepl("_recent", .data[[pop_col]]) ~ "Recent",
        TRUE ~ "Combined"
      )
    )
}

################################################################################
# 1.  PIXY pi — nucleotide diversity
################################################################################

pi2 <- read.table(file.path(dir_2pop, "pixy_pi.txt"), header = TRUE, sep = "\t")
pi4 <- read.table(file.path(dir_4pop, "pixy_pi.txt"), header = TRUE, sep = "\t")

pi2_summary <- pi2 %>%
  filter(no_sites > 50) %>%
  group_by(pop) %>%
  summarise(mean_pi = mean(avg_pi, na.rm = TRUE),
            sd_pi   = sd(avg_pi,   na.rm = TRUE),
            n       = n()) %>%
  mutate(se_pi = sd_pi / sqrt(n)) %>%
  add_meta("pop") %>%
  mutate(dataset = "2-pop")

pi4_summary <- pi4 %>%
  filter(no_sites > 50) %>%
  group_by(pop) %>%
  summarise(mean_pi = mean(avg_pi, na.rm = TRUE),
            sd_pi   = sd(avg_pi,   na.rm = TRUE),
            n       = n()) %>%
  mutate(se_pi = sd_pi / sqrt(n)) %>%
  add_meta("pop") %>%
  mutate(dataset = "4-pop")

# Plot: 4-pop Old->Recent per locality
p_pi <- ggplot(
  pi4_summary %>% mutate(time = factor(time, levels = c("Old", "Recent"))),
  aes(x = time, y = mean_pi, colour = locality, group = locality)
) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3.5, shape = 16) +
  geom_errorbar(aes(ymin = mean_pi - se_pi,
                    ymax = mean_pi + se_pi),
                width = 0.08, linewidth = 0.7) +
  geom_text(aes(label = n), vjust = -1.2, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = loc_colors) +
  labs(title  = expression(italic("H. semele") ~ "- Nucleotide Diversity (pi)"),
       x      = "Time period",
       y      = "Mean pi (+/-SE)",
       colour = "Locality") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right")

print(p_pi)
ggsave(file.path(out_dir, "semele_plot_pi.png"), p_pi, width = 6, height = 4.5, dpi = 300)

################################################################################
# 2.  PIXY Dxy — absolute divergence
################################################################################

dxy2 <- read.table(file.path(dir_2pop, "pixy_dxy.txt"), header = TRUE, sep = "\t")
dxy4 <- read.table(file.path(dir_4pop, "pixy_dxy.txt"), header = TRUE, sep = "\t")

dxy2_summary <- dxy2 %>%
  filter(no_sites > 50) %>%
  group_by(pop1, pop2) %>%
  summarise(mean_dxy = mean(avg_dxy, na.rm = TRUE),
            sd_dxy   = sd(avg_dxy,   na.rm = TRUE),
            n        = n(), .groups = "drop") %>%
  mutate(se_dxy = sd_dxy / sqrt(n),
         comparison = paste(pop1, "vs", pop2), dataset = "2-pop")

dxy4_summary <- dxy4 %>%
  filter(no_sites > 50) %>%
  group_by(pop1, pop2) %>%
  summarise(mean_dxy = mean(avg_dxy, na.rm = TRUE),
            sd_dxy   = sd(avg_dxy,   na.rm = TRUE),
            n        = n(), .groups = "drop") %>%
  mutate(se_dxy = sd_dxy / sqrt(n),
         comparison = paste(pop1, "vs", pop2), dataset = "4-pop")

dxy_all <- bind_rows(dxy2_summary, dxy4_summary)

p_dxy <- ggplot(dxy_all,
                aes(x = comparison, y = mean_dxy,
                    colour = dataset, shape = dataset)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = mean_dxy - se_dxy,
                    ymax = mean_dxy + se_dxy),
                width = 0.15, linewidth = 0.7,
                position = position_dodge(width = 0.3)) +
  geom_text(aes(label = n), vjust = -1.2, size = 2.8,
            position = position_dodge(width = 0.3),
            show.legend = FALSE) +
  scale_colour_manual(values = c("2-pop" = "#9C27B0", "4-pop" = "#4CAF50")) +
  scale_shape_manual(values  = c("2-pop" = 17,        "4-pop" = 16)) +
  labs(title  = expression(italic("H. semele") ~ "- Absolute Divergence (Dxy)"),
       x      = "Population comparison",
       y      = "Mean Dxy (+/-SE)",
       colour = "Dataset", shape = "Dataset") +
  theme_bw(base_size = 12) +
  theme(axis.text.x      = element_text(angle = 35, hjust = 1, vjust=1, size = 8),
        panel.grid.minor = element_blank(),
        legend.position  = "right",
        plot.margin      = margin(t = 5, r = 5, b = 5, l = 30, unit = "mm"))

print(p_dxy)
ggsave(file.path(out_dir, "semele_plot_dxy.png"), p_dxy, width = 8, height = 4.5, dpi = 300)

################################################################################
# 3.  VCF-based stats: Ho, Hs, Fis  (hierfstat)
################################################################################

library(vcfR)
library(adegenet)
library(hierfstat)

run_hierfstat <- function(vcf_path, pop_path, label) {
  message("Reading VCF for ", label, " ...")
  vcf    <- read.vcfR(vcf_path, verbose = FALSE)
  genind <- vcfR2genind(vcf)

  pop_info <- read.table(pop_path, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE)
  # Expected: V1 = SampleID (full BAM path), V2 = Population

  # Try an exact match first (VCF sample names == full BAM paths).
  matched <- match(indNames(genind), pop_info$V1)

  # If that fails entirely, the VCF likely only stores basenames
  # (e.g. "P33009_1001_S215_L006.final" or similar) -> fall back to
  # matching on basenames from the pop table instead.
  if (all(is.na(matched))) {
    message("  Exact path match failed for ", label,
            " - falling back to basename matching.")
    matched <- match(indNames(genind), basename(pop_info$V1))
  }

  pop_vec <- pop_info$V2[matched]

  if (any(is.na(pop_vec))) {
    warning("Sample IDs in VCF not found in pop table for ", label, ": ",
            paste(indNames(genind)[is.na(pop_vec)], collapse = ", "))
  }

  pop(genind) <- pop_vec
  gf    <- genind2hierfstat(genind)
  stats <- basic.stats(gf)

  # NOTE: n per statistic is tracked separately because Fis (and sometimes
  # Ho/Hs) can be NA at individual loci within a population (e.g. Fis is
  # undefined at monomorphic loci -> division by zero), so the number of
  # loci actually contributing to each mean can differ between statistics
  # and between populations.
  data.frame(
    Population = colnames(stats$Ho),
    Ho         = colMeans(stats$Ho,  na.rm = TRUE),
    Hs         = colMeans(stats$Hs,  na.rm = TRUE),
    Fis        = colMeans(stats$Fis, na.rm = TRUE),
    Ho_sd      = apply(stats$Ho,  2, sd, na.rm = TRUE),
    Hs_sd      = apply(stats$Hs,  2, sd, na.rm = TRUE),
    Fis_sd     = apply(stats$Fis, 2, sd, na.rm = TRUE),
    n_Ho       = colSums(!is.na(stats$Ho)),
    n_Hs       = colSums(!is.na(stats$Hs)),
    n_Fis      = colSums(!is.na(stats$Fis)),
    dataset    = label,
    row.names  = NULL
  )
}

vcf_stats_2pop <- run_hierfstat(vcf_file, pop_2pop, "2-pop")
vcf_stats_4pop <- run_hierfstat(vcf_file, pop_4pop, "4-pop")
vcf_stats_all  <- bind_rows(vcf_stats_2pop, vcf_stats_4pop)

# Plot Ho / Hs / Fis — 4-pop, Old->Recent per locality
vcf4_long <- vcf_stats_4pop %>%
  add_meta("Population") %>%
  mutate(time = factor(time, levels = c("Old", "Recent"))) %>%
  pivot_longer(cols = c(Ho, Hs, Fis),
               names_to  = "Statistic",
               values_to = "Value") %>%
  mutate(
    SD = case_when(Statistic == "Ho"  ~ Ho_sd,
                    Statistic == "Hs"  ~ Hs_sd,
                    Statistic == "Fis" ~ Fis_sd),
    N  = case_when(Statistic == "Ho"  ~ n_Ho,
                    Statistic == "Hs"  ~ n_Hs,
                    Statistic == "Fis" ~ n_Fis),
    SE = SD / sqrt(N)
  )

p_vcf <- ggplot(vcf4_long,
                aes(x = time, y = Value,
                    colour = locality, group = locality)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3.5, shape = 16) +
  geom_errorbar(aes(ymin = Value - SE, ymax = Value + SE),
                width = 0.08, linewidth = 0.7) +
  scale_colour_manual(values = loc_colors) +
  facet_wrap(~ Statistic, scales = "free_y") +
  labs(title  = expression(italic("H. semele") ~ "- VCF-based Statistics (Ho, Hs, Fis)"),
       x      = "Time period",
       y      = "Value (+/-SE across loci)",
       colour = "Locality") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "right")

print(p_vcf)
ggsave(file.path(out_dir, "semele_plot_vcf_stats.png"), p_vcf, width = 9, height = 4, dpi = 300)

################################################################################
# 4.  Export combined summary table
################################################################################

pi_combined <- bind_rows(pi2_summary, pi4_summary)

summary_table <- pi_combined %>%
  select(pop, locality, time, dataset, mean_pi, se_pi, n) %>%
  left_join(
    vcf_stats_all %>% select(Population, dataset, Ho, Hs, Fis),
    by = c("pop" = "Population", "dataset")
  )

write.table(summary_table,
            file.path(out_dir, "semele_diversity_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Done! Plots saved to: ", out_dir)
message("  semele_plot_pi.png | semele_plot_dxy.png | semele_plot_vcf_stats.png")
message("  (Fst excluded - handled separately)")
message("  semele_diversity_summary.txt")

```



### Assessing statistical significanve on the decline in pi
Using pi results and the emmeans R package also statistical significance on a linear model can be computed with the following script. In my thesis i only used the p-values and omitted the figures themselves.
```
## =============================================================
## Nucleotide diversity (pi) analysis: 2-pop and 4-pop comparisons
## Models: lmer with window as random effect
## Output: pi_analysis_report.txt
## =============================================================
setwd("C:/Users/username1/Documents/BT_R/pixy/statilinus/emmeans/")
library(tidyverse)
library(lme4)
library(emmeans)

# ---- EDIT THESE PATHS ----
file_2pop <- "pixy_pi_2pop.txt"   # Germany vs Slovakia
file_4pop <- "pixy_pi_4pop.txt"   # locality_old / locality_recent

out_file  <- "pi_analysis_report.txt"

## small helper: send a header + object's printed output into the report
log_section <- function(title, expr_result, con) {
  cat("\n", strrep("=", 70), "\n", sep = "", file = con, append = TRUE)
  cat(title, "\n", file = con, append = TRUE)
  cat(strrep("=", 70), "\n", file = con, append = TRUE)
  capture_text <- capture.output(print(expr_result))
  cat(capture_text, sep = "\n", file = con, append = TRUE)
}

## start fresh report
cat("Nucleotide diversity (pi) analysis report\n",
    "Generated:", as.character(Sys.time()), "\n",
    file = out_file)

## =============================================================
## PART 1: 2-population analysis (CZ_GER vs AUT_HUN_SK)
## =============================================================

pi2 <- read_tsv(file_2pop, show_col_types = FALSE)

pi2_f <- pi2 %>%
  filter(!is.na(avg_pi)) %>%
  mutate(window_id = paste(chromosome, window_pos_1, sep = "_"))

cat("\n", strrep("#", 70), "\n", "PART 1: 2-POPULATION ANALYSIS\n", strrep("#", 70), "\n",
    file = out_file, append = TRUE)

log_section("Sample sizes per population (after removing NA windows)",
            table(pi2_f$pop), out_file)

mod1 <- lmer(avg_pi ~ pop + (1 | window_id), data = pi2_f)

log_section("Model summary: avg_pi ~ pop + (1|window_id)",
            summary(mod1), out_file)

log_section("Overall F-test for 'pop' effect",
            joint_tests(mod1), out_file)

emm1 <- emmeans(mod1, ~ pop)
log_section("Estimated marginal means (pop)", emm1, out_file)

log_section("Pairwise comparison (Tukey-adjusted)",
            pairs(emm1, adjust = "tukey"), out_file)

## save diagnostic plots
png("diagnostics_2pop.png", width = 1000, height = 800)
par(mfrow = c(2,2))
plot(fitted(mod1), resid(mod1), main = "Residuals vs Fitted (2pop)")
abline(h = 0, col = "red")
qqnorm(resid(mod1), main = "QQ plot (2pop)"); qqline(resid(mod1))
hist(resid(mod1), main = "Residual histogram (2pop)", breaks = 40)
dev.off()


## =============================================================
## PART 2: 4-population analysis (locality x age)
## =============================================================

pi4 <- read_tsv(file_4pop, show_col_types = FALSE)

pi4_f <- pi4 %>%
  filter(!is.na(avg_pi)) %>%
  separate(pop, into = c("locality", "age"), sep = "_(?=old|recent)") %>%
  mutate(window_id = paste(chromosome, window_pos_1, sep = "_"))

cat("\n", strrep("#", 70), "\n", "PART 2: 4-POPULATION ANALYSIS (locality x age)\n", strrep("#", 70), "\n",
    file = out_file, append = TRUE)

log_section("Sample sizes per locality x age group (after removing NA windows)",
            table(pi4_f$locality, pi4_f$age), out_file)

mod2 <- lmer(avg_pi ~ locality * age + (1 | window_id), data = pi4_f)

log_section("Model summary: avg_pi ~ locality * age + (1|window_id)",
            summary(mod2), out_file)

log_section("Overall F-tests (locality, age, interaction)",
            joint_tests(mod2), out_file)

emm2 <- emmeans(mod2, ~ locality * age)
log_section("Estimated marginal means (locality x age)", emm2, out_file)

log_section("All pairwise comparisons (Tukey-adjusted)",
            pairs(emm2, adjust = "tukey"), out_file)

log_section("Simple effects: old vs recent, within each locality",
            emmeans(mod2, ~ age | locality) %>% pairs(adjust = "tukey"), out_file)

log_section("Simple effects: locality comparison, within each age class",
            emmeans(mod2, ~ locality | age) %>% pairs(adjust = "tukey"), out_file)

## save diagnostic plots
png("diagnostics_4pop.png", width = 1000, height = 800)
par(mfrow = c(2,2))
plot(fitted(mod2), resid(mod2), main = "Residuals vs Fitted (4pop)")
abline(h = 0, col = "red")
qqnorm(resid(mod2), main = "QQ plot (4pop)"); qqline(resid(mod2))
hist(resid(mod2), main = "Residual histogram (4pop)", breaks = 40)
dev.off()

emm2_df <- as.data.frame(emm2)
print(names(emm2_df))   # <-- add this temporarily

simple_age <- emmeans(mod2, ~ age | locality) %>% pairs(adjust = "tukey")
simple_age_df <- as.data.frame(simple_age)

ann_df <- simple_age_df %>%
  mutate(
    p_label = paste0("p = ", formatC(p.value, digits = 3, format = "f")),
    x = 1.5
  ) %>%
  left_join(
    emm2_df %>% group_by(locality) %>% summarise(y_max = max(asymp.UCL)),
    by = "locality"
  ) %>%
  mutate(y = y_max * 1.05)

p <- ggplot(emm2_df, aes(x = age, y = emmean, color = locality, group = locality)) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  geom_text(data = ann_df, aes(x = x, y = y, label = p_label, color = locality),
            show.legend = FALSE, size = 3.5, fontface = "italic") +
  labs(y = "Estimated mean pi (95% CI)", x = "Time period",
       title = expression("Model-estimated nucleotide diversity (pi) for " * italic("H. statilinus"))) +
  theme_minimal()

ggsave("pi_emmeans_plot.png", p, width = 7, height = 5)
```


## Fst with VCFtools
To obtain the Fst, the Invariants.vcf file is required and the population files are different from pixy. There should be population files with the names specified. Since in the .list files i had only the sample codes, i added the paths that are present in the vcf file in the script:
```
#!/bin/bash
#PBS -N Fst_sem_AUT_HUN_SK
#PBS -l select=1:ncpus=1:mem=15gb:scratch_local=100gb
#PBS -l walltime=24:00:00
#PBS -e /storage/brno12-cerit/home/username/PBS
#PBS -o /storage/brno12-cerit/home/username/PBS

DATADIR=/storage/brno12-cerit/home/username/files_pixy/all_files_semele
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
## clean scratch after the end
trap 'clean_scratch' TERM EXIT
cd $SCRATCHDIR

cp $DATADIR/SNPsemele_autosomes.vcf.gz* $SCRATCHDIR/
cp $DATADIR/*.list $SCRATCHDIR/


out="$DATADIR/Fst"
mkdir -p $out

## run vcftools
module load vcftools

## spp file lists (list of spp samples in vcf file)
PREFIX="/storage/frontends/home/username/allign_semele/finalbam/"

sed "s|^|${PREFIX}|; s|\$|.final.bam|" AUT_HUN_SK_old.list > AUT_HUN_SK_old.vcf.list

sed "s|^|${PREFIX}|; s|\$|.final.bam|" AUT_HUN_SK_recent.list > AUT_HUN_SK_recent.vcf.list

vcftools --gzvcf SNPsemele_autosomes.vcf.gz \
  --weir-fst-pop AUT_HUN_SK_old.vcf.list \
  --weir-fst-pop AUT_HUN_SK_recent.vcf.list \
  --out ${out}/Fst_AUT_old_vs_recent
  
########################################
#!/bin/bash
#PBS -N Fst_sem_CZ_GER
#PBS -l select=1:ncpus=1:mem=15gb:scratch_local=100gb
#PBS -l walltime=24:00:00
#PBS -e /storage/brno12-cerit/home/username/PBS
#PBS -o /storage/brno12-cerit/home/username/PBS

DATADIR=/storage/brno12-cerit/home/username/files_pixy/all_files_semele
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
## clean scratch after the end
trap 'clean_scratch' TERM EXIT
cd $SCRATCHDIR

cp $DATADIR/SNPsemele_autosomes.vcf.gz* $SCRATCHDIR/
cp $DATADIR/*.list $SCRATCHDIR/


out="$DATADIR/Fst"
mkdir -p $out

## run vcftools
module load vcftools

## spp file lists (list of spp samples in vcf file)
PREFIX="/storage/frontends/home/username/allign_semele/finalbam/"

sed "s|^|${PREFIX}|; s|\$|.final.bam|" CZ_GER_old.list > CZ_GER_old.vcf.list

sed "s|^|${PREFIX}|; s|\$|.final.bam|" CZ_GER_recent.list > CZ_GER_recent.vcf.list

vcftools --gzvcf SNPsemele_autosomes.vcf.gz \
  --weir-fst-pop CZ_GER_old.vcf.list \
  --weir-fst-pop CZ_GER_recent.vcf.list \
  --out ${out}/Fst_CZ_old_vs_recent
  
  
 ###########################
#!/bin/bash
#PBS -N Fst_sem_loc
#PBS -l select=1:ncpus=1:mem=15gb:scratch_local=100gb
#PBS -l walltime=24:00:00
#PBS -e /storage/brno12-cerit/home/username/PBS
#PBS -o /storage/brno12-cerit/home/username/PBS

DATADIR=/storage/brno12-cerit/home/username/files_pixy/all_files_semele
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
## clean scratch after the end
trap 'clean_scratch' TERM EXIT
cd $SCRATCHDIR

cp $DATADIR/SNPsemele_autosomes.vcf.gz* $SCRATCHDIR/
cp $DATADIR/*.list $SCRATCHDIR/


out="$DATADIR/Fst"
mkdir -p $out

## run vcftools
module load vcftools

## spp file lists (list of spp samples in vcf file)
PREFIX="/storage/frontends/home/username/allign_semele/finalbam/"

sed "s|^|${PREFIX}|; s|\$|.final.bam|" AUT_HUN_SK.list > AUT_HUN_SK.vcf.list

sed "s|^|${PREFIX}|; s|\$|.final.bam|" CZ_GER.list > CZ_GER.vcf.list

vcftools --gzvcf SNPsemele_autosomes.vcf.gz \
  --weir-fst-pop AUT_HUN_SK.vcf.list  \
  --weir-fst-pop CZ_GER.vcf.list \
  --out ${out}/Fst_AUT_HUN_SK_vs_CZ_GER
  
```


The output files are uploaded to R and the following code was run:
```
################################################################################
# Semele Population Genomics: vcftools Weir & Cockerham Fst Visualization
# Mirrors the pixy/plink Fst plot style (mean +/- SD per comparison)
#
# vcftools --weir-fst-pop output format (per SNP, not windowed):
#   CHROM   POS   WEIR_AND_COCKERHAM_FST
################################################################################

setwd("C:/Users/username1/Documents/BT_R/Fst")
library(dplyr)
library(ggplot2)

# ── File paths ────────────────────────────────────────────────────────────────
# One file per pairwise comparison, as produced by:
#   vcftools --weir-fst-pop popA.txt --weir-fst-pop popB.txt ...
fst_dir <- "C:/Users/username1/Documents/BT_R/Fst"

fst_files <- c(
  "Fst_CZ_old_vs_recent.weir.fst"      = "CZ_GER_old vs CZ_GER_recent",
  "Fst_AUT_old_vs_recent.weir.fst"     = "AUT_HUN_SK_old vs AUT_HUN_SK_recent",
  "Fst_AUT_HUN_SK_vs_CZ_GER.weir.fst"  = "CZ_GER vs AUT_HUN_SK"
)

out_dir <- "C:/Users/username1/Documents/BT_R/Fst"

# ── Helper: classify comparison type from its label ──────────────────────────
classify_type <- function(comparison) {
  parts <- strsplit(comparison, " vs ")[[1]]
  pop1 <- parts[1]; pop2 <- parts[2]
  
  loc_of <- function(x) {
    case_when(
      grepl("CZ_GER",     x) ~ "CZ_GER",
      grepl("AUT_HUN_SK", x) ~ "AUT_HUN_SK",
      TRUE ~ NA_character_
    )
  }
  loc1 <- loc_of(pop1); loc2 <- loc_of(pop2)
  
  if (!is.na(loc1) && !is.na(loc2) && loc1 == loc2) {
    "Within locality (temporal)"
  } else {
    "Between localities"
  }
}

# ── Read + summarise each per-site Fst file ───────────────────────────────────
read_vcftools_fst <- function(path, comparison) {
  df <- read.table(path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE,
                   col.names = c("CHROM", "POS", "fst"))
  
  # vcftools codes undefined Fst (zero variance) as -nan / NaN -> read as NA
  df <- df %>% filter(!is.na(fst))
  
  df %>%
    summarise(
      mean_fst = mean(fst, na.rm = TRUE),
      sd_fst   = sd(fst,   na.rm = TRUE),
      n        = n()
    ) %>%
    mutate(se_fst = sd_fst / sqrt(n)) %>%
    mutate(
      comparison = comparison,
      type       = classify_type(comparison),
      dataset    = "vcftools"
    )
}

fst_summary_list <- lapply(names(fst_files), function(fname) {
  read_vcftools_fst(file.path(fst_dir, fname), fst_files[[fname]])
})

fst_all <- bind_rows(fst_summary_list)

# ── Plot (same style as the pixy/plink Fst plot) ─────────────────────────────
p_fst_vcftools <- ggplot(fst_all,
                         aes(x = comparison, y = mean_fst, colour = type)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = mean_fst - se_fst,
                    ymax = mean_fst + se_fst),
                width = 0.15, linewidth = 0.7,
                position = position_dodge(width = 0.3)) +
  geom_text(aes(label = n), vjust = -1.5, size = 2.8,
            position = position_dodge(width = 0.3),
            show.legend = FALSE) +
  scale_colour_manual(values = c("Within locality (temporal)" = "#2196F3",
                                 "Between localities"         = "#E91E63")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.15))) +
  labs(title  = expression(italic("Semele") ~ "- Pairwise Fst (vcftools, Weir & Cockerham)"),
       x      = "Population comparison",
       y      = "Mean Fst (+/-SE across SNPs)",
       colour = "Comparison type") +
  theme_bw(base_size = 12) +
  theme(axis.text.x      = element_text(angle = 35, hjust = 1, size = 8),
        panel.grid.minor = element_blank(),
        legend.position  = "right",
        plot.margin      = margin(t = 15, r = 10, b = 40, l = 10))

print(p_fst_vcftools)
ggsave(file.path(out_dir, "plot_fst_vcftools.png"), p_fst_vcftools,
       width = 11, height = 5.5, dpi = 300)

# ── Export summary table ─────────────────────────────────────────────────────
write.table(fst_all,
            file.path(out_dir, "vcftools_fst_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Done! Saved: plot_fst_vcftools.png and vcftools_fst_summary.txt")
```


### Creating a table for the obtained diversity metrics
To create a table with the exact values use this R script:
```
################################################################################
# Hipparchia semele - Appendix table export (mean +/- SE only)
#
# Produces ONE clean, Excel-ready txt table (tab-separated) with just the
# columns you need for an appendix: which analysis/comparison, mean, SE, n.
# No SD, no extra internal columns - just what you'd want to paste into a
# Word/Excel appendix table.
#
# Open the output file in Excel via: Data > Get Data > From Text/CSV
# (or just double-click it - Excel will usually auto-detect the tab delimiter)
################################################################################

library(dplyr)
library(tidyr)

# ── File paths (keep in sync with your other scripts) ────────────────────────
dir_2pop <- "C:/Users/username1/Documents/BT_R/pixy/semele/pixy_2pop_semele/"
dir_4pop <- "C:/Users/username1/Documents/BT_R/pixy/semele/pixy_4pop_semele/"
vcf_file <- "C:/Users/username1/Documents/BT_R/pixy/semele/SNPsemele_autosomes.vcf.gz"
pop_2pop <- file.path(dir_2pop, "2pop_semele.txt")
pop_4pop <- file.path(dir_4pop, "4pop_semele.txt")

fst_dir <- "C:/Users/username1/Documents/BT_R/Fst"
fst_files <- c(
  "Fst_CZ_old_vs_recent.weir.fst"      = "CZ_GER_old vs CZ_GER_recent",
  "Fst_AUT_old_vs_recent.weir.fst"     = "AUT_HUN_SK_old vs AUT_HUN_SK_recent",
  "Fst_AUT_HUN_SK_vs_CZ_GER.weir.fst"  = "CZ_GER vs AUT_HUN_SK"
)

out_dir  <- "C:/Users/username1/Documents/BT_R/diversity_se/semele/"
out_file <- file.path(out_dir, "semele_appendix_table_mean_SE.txt")

################################################################################
# 1. PIXY pi
################################################################################
pi2 <- read.table(file.path(dir_2pop, "pixy_pi.txt"), header = TRUE, sep = "\t")
pi4 <- read.table(file.path(dir_4pop, "pixy_pi.txt"), header = TRUE, sep = "\t")

pi_table <- bind_rows(
  pi2 %>% mutate(Dataset = "2-pop"),
  pi4 %>% mutate(Dataset = "4-pop")
) %>%
  filter(no_sites > 50) %>%
  group_by(pop, Dataset) %>%
  summarise(mean = mean(avg_pi, na.rm = TRUE),
            sd   = sd(avg_pi,   na.rm = TRUE),
            n    = n(), .groups = "drop") %>%
  transmute(
    Statistic  = "pi",
    Comparison = pop,
    Dataset,
    Mean = round(mean, 5),
    SE   = round(sd / sqrt(n), 5),
    n
  )

################################################################################
# 2. PIXY Dxy
################################################################################
dxy2 <- read.table(file.path(dir_2pop, "pixy_dxy.txt"), header = TRUE, sep = "\t")
dxy4 <- read.table(file.path(dir_4pop, "pixy_dxy.txt"), header = TRUE, sep = "\t")

dxy_table <- bind_rows(
  dxy2 %>% mutate(Dataset = "2-pop"),
  dxy4 %>% mutate(Dataset = "4-pop")
) %>%
  filter(no_sites > 50) %>%
  group_by(pop1, pop2, Dataset) %>%
  summarise(mean = mean(avg_dxy, na.rm = TRUE),
            sd   = sd(avg_dxy,   na.rm = TRUE),
            n    = n(), .groups = "drop") %>%
  transmute(
    Statistic  = "Dxy",
    Comparison = paste(pop1, "vs", pop2),
    Dataset,
    Mean = round(mean, 5),
    SE   = round(sd / sqrt(n), 5),
    n
  )

################################################################################
# 3. VCF-based stats: Ho, Hs, Fis (hierfstat)
################################################################################
library(vcfR)
library(adegenet)
library(hierfstat)

run_hierfstat <- function(vcf_path, pop_path, label) {
  message("Reading VCF for ", label, " ...")
  vcf    <- read.vcfR(vcf_path, verbose = FALSE)
  genind <- vcfR2genind(vcf)

  pop_info <- read.table(pop_path, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE)

  matched <- match(indNames(genind), pop_info$V1)
  if (all(is.na(matched))) {
    message("  Exact path match failed for ", label,
            " - falling back to basename matching.")
    matched <- match(indNames(genind), basename(pop_info$V1))
  }

  pop_vec <- pop_info$V2[matched]
  if (any(is.na(pop_vec))) {
    warning("Sample IDs in VCF not found in pop table for ", label, ": ",
            paste(indNames(genind)[is.na(pop_vec)], collapse = ", "))
  }

  pop(genind) <- pop_vec
  gf    <- genind2hierfstat(genind)
  stats <- basic.stats(gf)

  data.frame(
    Population = colnames(stats$Ho),
    Ho         = colMeans(stats$Ho,  na.rm = TRUE),
    Hs         = colMeans(stats$Hs,  na.rm = TRUE),
    Fis        = colMeans(stats$Fis, na.rm = TRUE),
    Ho_sd      = apply(stats$Ho,  2, sd, na.rm = TRUE),
    Hs_sd      = apply(stats$Hs,  2, sd, na.rm = TRUE),
    Fis_sd     = apply(stats$Fis, 2, sd, na.rm = TRUE),
    n_Ho       = colSums(!is.na(stats$Ho)),
    n_Hs       = colSums(!is.na(stats$Hs)),
    n_Fis      = colSums(!is.na(stats$Fis)),
    Dataset    = label,
    row.names  = NULL
  )
}

vcf_stats_2pop <- run_hierfstat(vcf_file, pop_2pop, "2-pop")
vcf_stats_4pop <- run_hierfstat(vcf_file, pop_4pop, "4-pop")
vcf_stats_all  <- bind_rows(vcf_stats_2pop, vcf_stats_4pop)

vcf_table <- vcf_stats_all %>%
  pivot_longer(cols = c(Ho, Hs, Fis), names_to = "Statistic", values_to = "mean") %>%
  mutate(
    sd = case_when(Statistic == "Ho"  ~ Ho_sd,
                   Statistic == "Hs"  ~ Hs_sd,
                   Statistic == "Fis" ~ Fis_sd),
    n  = case_when(Statistic == "Ho"  ~ n_Ho,
                   Statistic == "Hs"  ~ n_Hs,
                   Statistic == "Fis" ~ n_Fis)
  ) %>%
  transmute(
    Statistic,
    Comparison = Population,
    Dataset,
    Mean = round(mean, 5),
    SE   = round(sd / sqrt(n), 5),
    n
  )

################################################################################
# 4. vcftools Weir & Cockerham Fst
################################################################################
read_vcftools_fst <- function(path, comparison) {
  df <- read.table(path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE,
                   col.names = c("CHROM", "POS", "fst"))
  df <- df %>% filter(!is.na(fst))
  df %>%
    summarise(mean = mean(fst, na.rm = TRUE),
              sd   = sd(fst,   na.rm = TRUE),
              n    = n()) %>%
    transmute(
      Statistic  = "Fst",
      Comparison = comparison,
      Dataset    = "vcftools",
      Mean = round(mean, 5),
      SE   = round(sd / sqrt(n), 5),
      n
    )
}

fst_table <- bind_rows(lapply(names(fst_files), function(fname) {
  read_vcftools_fst(file.path(fst_dir, fname), fst_files[[fname]])
}))

################################################################################
# 5. Combine and export ONE clean appendix table
################################################################################
appendix_table <- bind_rows(pi_table, dxy_table, vcf_table, fst_table)

write.table(appendix_table, out_file,
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Done! Appendix table (Statistic | Comparison | Dataset | Mean | SE | n) written to:")
message("  ", out_file)
message("Open in Excel: double-click the file, or Data > Get Data > From Text/CSV")
```





## Fitting a model over the Ho values

Run this command in metacentrum:
```
module load vcftools/
(BOOKWORM)username@zenith:~/files_pixy/all_files_semele$ vcftools --gzvcf allsites_autosomes_semele.vcf.gz --het --out semele_het_per_ind
```
Use the output in R and run the following script:
```
library(tidyverse)
library(emmeans)

log_section <- function(title, expr_result, con) {
  cat("\n", strrep("=", 70), "\n", sep = "", file = con, append = TRUE)
  cat(title, "\n", file = con, append = TRUE)
  cat(strrep("=", 70), "\n", file = con, append = TRUE)
  cat(capture.output(print(expr_result)), sep = "\n", file = con, append = TRUE)
}

analyze_pi_over_time <- function(het_file, meta_file, group_col, species_label, out_png, report_file) {
  
  cat("Pi analysis report for", species_label, "\nGenerated:", as.character(Sys.time()), "\n",
      file = report_file)
  
  het <- read_tsv(het_file, show_col_types = FALSE) %>%
    mutate(
      old_id = str_remove(basename(INDV), "\\.final\\.bam$"),
      heterozygosity = (N_SITES - `O(HOM)`) / N_SITES
    )
  
  meta <- read_tsv(meta_file, show_col_types = FALSE)
  
  pi_full <- het %>%
    left_join(meta, by = "old_id") %>%
    mutate(year = as.numeric(year))
  
  if (any(is.na(pi_full[[group_col]]))) {
    warning("Unmatched samples in ", species_label, ": ",
            paste(pi_full$old_id[is.na(pi_full[[group_col]])], collapse = ", "))
  }
  
  form <- as.formula(paste("heterozygosity ~", group_col, "* year"))
  mod  <- lm(form, data = pi_full)
  
  log_section(paste(species_label, "- model summary"), summary(mod), report_file)
  
  aov_tab <- anova(mod)
  log_section(paste(species_label, "- ANOVA"), aov_tab, report_file)
  
  ## ---- pull out p-values to annotate ----
  group_p <- aov_tab[group_col, "Pr(>F)"]
  year_p  <- aov_tab["year", "Pr(>F)"]
  inter_row <- rownames(aov_tab)[grepl(":", rownames(aov_tab))]
  inter_p <- aov_tab[inter_row, "Pr(>F)"]
  
  fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, digits = 3, format = "f")))
  
  trend_form <- as.formula(paste("~ year |", group_col))
  emt <- emtrends(mod, trend_form, var = "year")
  log_section(paste(species_label, "- estimated slopes"), emt, report_file)
  emt_pairs <- pairs(emt)
  log_section(paste(species_label, "- slope comparison between groups"), emt_pairs, report_file)
  
  ## per-group slope p-values, for labels placed at the end of each trend line
  emt <- emtrends(mod, trend_form, var = "year")
  log_section(paste(species_label, "- estimated slopes"), emt, report_file)
  
  emt_test <- test(emt)   # <-- this adds t.ratio and p.value (tests slope = 0)
  log_section(paste(species_label, "- per-group slope significance"), emt_test, report_file)
  
  emt_pairs <- pairs(emt)
  log_section(paste(species_label, "- slope comparison between groups"), emt_pairs, report_file)
  
  ## per-group slope p-values, for labels placed at the end of each trend line
  emt_df <- as.data.frame(emt_test) %>%
    rename(slope_p = p.value) %>%
    mutate(x = max(pi_full$year, na.rm = TRUE),
           label = paste0(.data[[group_col]], ": ", fmt_p(slope_p)))
  
  y_positions <- pi_full %>%
    group_by(.data[[group_col]]) %>%
    summarise(y = max(heterozygosity, na.rm = TRUE), .groups = "drop")
  
  emt_df <- emt_df %>% left_join(y_positions, by = group_col)
  
  subtitle_txt <- paste0(
    "Group ", fmt_p(group_p), "  |  Year ", fmt_p(year_p),
    "  |  Group\u00d7Year ", fmt_p(inter_p)
  )
  
  p <- ggplot(pi_full, aes(x = year, y = heterozygosity, colour = .data[[group_col]])) +
    geom_jitter(width = 0.6, height = 0, alpha = 0.75, size = 2.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
    geom_text(data = emt_df, aes(x = x, y = y, label = label, colour = .data[[group_col]]),
              hjust = 1, vjust = -0.5, size = 3.2, fontface = "italic", show.legend = FALSE) +
    labs(x = "Collection year", y = "Individual heterozygosity",
         colour = str_to_title(group_col),
         title = bquote(italic(.(species_label)) ~ "- individual diversity over time"),
         subtitle = subtitle_txt) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 10, colour = "grey30"))
  
  print(p)
  ggsave(out_png, p, width = 7.5, height = 5, dpi = 300)
  
  message("Report written to: ", report_file)
  list(data = pi_full, model = mod, plot = p)
}

semele_res <- analyze_pi_over_time(
  het_file      = "semele_het_per_ind.het",
  meta_file     = "individual_metadata_semele.txt",
  group_col     = "group",
  species_label = "H. semele",
  out_png       = "heterozygosityividual_over_time_semele.png",
  report_file   = "semele_pi_report.txt"
)

statilinus_res <- analyze_pi_over_time(
  het_file      = "statilinus_het_per_ind.het",
  meta_file     = "individual_metadata_statilinus.txt",
  group_col     = "locality",
  species_label = "H. statilinus",
  out_png       = "heterozygosityividual_over_time_statilinus.png",
  report_file   = "statilinus_pi_report.txt"
)
```
To eventually obtain the last text file for this thesis containing the numerical values this code can be used:
```
################################################################################
# Hipparchia semele + H. statilinus - Heterozygosity-over-time appendix table
#
# Same scheme as the other appendix tables (one clean, Excel-ready txt file),
# adapted for this analysis: instead of just group means, this model produces
# slopes and significance tests, so those are included too.
#
# Output columns: Species | Statistic | Group | Estimate | SE | p_value | n
#
#   Statistic = "Mean heterozygosity"      -> per-group raw mean +/- SE (n = individuals)
#   Statistic = "Slope (het ~ year)"       -> per-group regression slope +/- SE, with p-value
#   Statistic = "ANOVA: <term>"            -> F-test p-value for Group / Year / Group:Year
#   Statistic = "Slope difference"         -> pairwise comparison of slopes between groups
#
# This does NOT change your original script's plots/reports - it just also
# extracts the same numbers into one combined table for your appendix.
#
# Open the output file in Excel via: Data > Get Data > From Text/CSV
################################################################################

library(tidyverse)
library(emmeans)

out_file <- "C:/Users/username1/Documents/BT_R/diversity_se/het_over_time_appendix_table.txt"

fmt_p <- function(p) round(p, 5)

################################################################################
# ── Core analysis + table extraction (per species) ──────────────────────────
################################################################################

analyze_pi_over_time <- function(het_file, meta_file, group_col, species_label) {

  het <- read_tsv(het_file, show_col_types = FALSE) %>%
    mutate(
      old_id = str_remove(basename(INDV), "\\.final\\.bam$"),
      heterozygosity = (N_SITES - `O(HOM)`) / N_SITES
    )

  meta <- read_tsv(meta_file, show_col_types = FALSE)

  pi_full <- het %>%
    left_join(meta, by = "old_id") %>%
    mutate(year = as.numeric(year))

  if (any(is.na(pi_full[[group_col]]))) {
    warning("Unmatched samples in ", species_label, ": ",
            paste(pi_full$old_id[is.na(pi_full[[group_col]])], collapse = ", "))
  }

  form <- as.formula(paste("heterozygosity ~", group_col, "* year"))
  mod  <- lm(form, data = pi_full)
  aov_tab <- anova(mod)

  ## ---- 1. Raw group means +/- SE (n = individuals) ----
  mean_table <- pi_full %>%
    group_by(.data[[group_col]]) %>%
    summarise(mean = mean(heterozygosity, na.rm = TRUE),
              sd   = sd(heterozygosity,   na.rm = TRUE),
              n    = n(), .groups = "drop") %>%
    transmute(
      Species   = species_label,
      Statistic = "Mean heterozygosity",
      Group     = .data[[group_col]],
      Estimate  = round(mean, 5),
      SE        = round(sd / sqrt(n), 5),
      p_value   = NA_real_,
      n
    )

  ## ---- 2. Per-group slopes (heterozygosity ~ year) +/- SE, with p-value ----
  trend_form <- as.formula(paste("~ year |", group_col))
  emt      <- emtrends(mod, trend_form, var = "year")
  emt_test <- test(emt)

  slope_table <- as.data.frame(emt_test) %>%
    transmute(
      Species   = species_label,
      Statistic = "Slope (het ~ year)",
      Group     = .data[[group_col]],
      Estimate  = round(year.trend, 6),
      SE        = round(SE, 6),
      p_value   = fmt_p(p.value),
      n         = NA_integer_
    )

  ## ---- 3. ANOVA effect p-values (Group, Year, Group x Year) ----
  inter_row <- rownames(aov_tab)[grepl(":", rownames(aov_tab))]

  anova_table <- tibble(
    Species   = species_label,
    Statistic = paste0("ANOVA: ", c(group_col, "year", inter_row)),
    Group     = NA_character_,
    Estimate  = round(aov_tab[c(group_col, "year", inter_row), "F value"], 4),
    SE        = NA_real_,
    p_value   = fmt_p(aov_tab[c(group_col, "year", inter_row), "Pr(>F)"]),
    n         = NA_integer_
  )

  ## ---- 4. Pairwise slope comparisons between groups ----
  emt_pairs <- pairs(emt)
  pairs_table <- as.data.frame(emt_pairs) %>%
    transmute(
      Species   = species_label,
      Statistic = "Slope difference",
      Group     = contrast,
      Estimate  = round(estimate, 6),
      SE        = round(SE, 6),
      p_value   = fmt_p(p.value),
      n         = NA_integer_
    )

  bind_rows(mean_table, slope_table, anova_table, pairs_table)
}

################################################################################
# ── Run for both species and combine ─────────────────────────────────────────
################################################################################

semele_table <- analyze_pi_over_time(
  het_file      = "semele_het_per_ind.het",
  meta_file     = "individual_metadata_semele.txt",
  group_col     = "group",
  species_label = "H. semele"
)

statilinus_table <- analyze_pi_over_time(
  het_file      = "statilinus_het_per_ind.het",
  meta_file     = "individual_metadata_statilinus.txt",
  group_col     = "locality",
  species_label = "H. statilinus"
)

appendix_table <- bind_rows(semele_table, statilinus_table)

write.table(appendix_table, out_file,
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

message("Done! Heterozygosity-over-time appendix table written to:")
message("  ", out_file)
message("Columns: Species | Statistic | Group | Estimate | SE | p_value | n")
message("Open in Excel: double-click the file, or Data > Get Data > From Text/CSV")
```
