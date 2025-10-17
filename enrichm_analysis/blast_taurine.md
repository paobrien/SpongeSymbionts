# Blast taurine transporter against genomes to check if present

### Make a genome database for each symbiont group and blast tauABC genes


```bash
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --job-name=blast
#SBATCH --time=6:00:00
#SBATCH --partition=general
#SBATCH --account=a_ace
#SBATCH -o /home/uqpobri2/scripts/out/blast.output
#SBATCH -e /home/uqpobri2/scripts/error/blast.error

##########################

# use blast to search taurine transporter against symbiont genomes (reviewer comment)

# load program
module load blast

# assign variables
GenomeDir="/scratch/user/uqpobri2/data/sponge_symbionts"
DataDir="/scratch/user/uqpobri2/data/sponge_symbionts/taurine"
OutDir="/scratch/user/uqpobri2/analysis/sponge_symbionts/blast"

mkdir -p $OutDir/db
mkdir -p $OutDir/taurine

# check variables
echo ""
echo "#====================#"
echo "#===== JOB INFO =====#"
echo "#====================#"
echo ""
echo "    SLURM Job ID      :  ${SLURM_JOB_ID}"
echo "    SLURM Job Name    :  ${SLURM_JOB_NAME}"
echo ""
echo "    Input genome path :  ${GenomeDir}"
echo "    Input sequences   :  ${DataDir}" 
echo "    Output directory  :  ${OutDir}"
echo ""
echo "#====================#"
echo ""

## blast to genomes

# concatenate genomes from each family to one file
for Dir in $(ls $GenomeDir | grep -v taurine); do
    # get dir name
    echo $Dir
   # concatenate files in dir
    cat $GenomeDir/${Dir}/*.fna > $GenomeDir/${Dir}_all.fna
done

# make symbiont genome database
for File in $(ls $GenomeDir/*all.fna); do
    # check file name
    echo $File
    # get symbiont name
    Name=$(basename $File | cut -f 1 -d '_')
    # check name
    echo $Name
    #make blast db
   srun makeblastdb \
        -in $File \
        -dbtype nucl \
        -out $OutDir/db/$Name/$Name
done

# tblastx to genome database (translated protein search)
for db in $(ls $OutDir/db); do
    #check db
    echo $db
    # blast sequences to db 
    srun tblastx \ 
        -query $DataDir/tauABC.fasta \
        -db $OutDir/db/$db/$db \
        -out $OutDir/taurine/taurine_tblastx_results_${db}.txt \
        -outfmt 6
done

```
The columns in the tabular output file will include:

1 Query sequence ID
2 Subject sequence ID
3 Percent identity
4 Alignment length
5 Number of mismatches
6 Number of gap openings
7 Start position in query sequence
8 End position in query sequence
9 Start position in subject sequence
10 End position in subject sequence
11 E-value
12 Bit score

```bash
# to view: sort by bit score - 12,12nr - reflects overall quality - higher the better
# note: nr = numeric, reverse order
sort -k 12,12nr $OutDir/taurine/taurine_blast_results_micro.txt | tsv_view | less -S

```

### Find which genomes have contigs that have a signficant match to tauABC


```bash
# output significant results, in this case percent identity > 30, aligment >100, e-value < 1e-03, bit score >50
# based on Pearson 2013

# from output dir..
for File in *tblastx*; do
    # get sample id
    echo $File
    SampleID=$(echo $File | cut -d '_' -f 4 )
    echo $SampleID
   # get significant results
    awk '$3 > 30 && $4 > 100 && $11 < 1e-03 && $12 > 50 {print}' $File | sort -k 12,12nr > taurine_tblastx_significant_results_e03-100_${SampleID}
done

# get unique gene/contig matches
for File in *tblastx_significant_results*; do
    # get sample id
    echo $File
    SampleID=$(echo $File | cut -d '_' -f 6 )
    echo $SampleID
    # output unique matches
    awk '{print $1, $2}' $File | sort | uniq > taurine_tblastx_unique_matches_e03-100_${SampleID}
done

# get genomes for matched contigs
# running each separately to avoid looping issues
# endo
while read -r line; do
    col1=$(echo "$line" | awk '{print $1}')
    col2=$(echo "$line" | awk '{print $2}')
    match_found=false

    for file in /scratch/user/uqpobri2/data/sponge_symbionts/endo/*; do
        GenomeID=$(basename "$file")
        if grep -q "$col2" "$file"; then
            echo "$col1 $col2 $GenomeID" >> genomes_with_taurine_e03-100_endo.txt
            match_found=true
            break
        fi
    done

    if [ "$match_found" = false ]; then
        echo "$col1 $col2" >> genomes_with_taurine_e03-100_endo.txt
    fi
done < taurine_tblastx_unique_matches_e03-100_endo.txt 

# micro
while read -r line; do
    col1=$(echo "$line" | awk '{print $1}')
    col2=$(echo "$line" | awk '{print $2}')
    match_found=false

    for file in /scratch/user/uqpobri2/data/sponge_symbionts/micro/*; do
        GenomeID=$(basename "$file")
        if grep -q "$col2" "$file"; then
            echo "$col1 $col2 $GenomeID" >> genomes_with_taurine_e03-100_micro.txt
            match_found=true
            break
        fi
    done

    if [ "$match_found" = false ]; then
        echo "$col1 $col2" >> genomes_with_taurine_e03-100_micro.txt
    fi
done < taurine_tblastx_unique_matches_e03-100_micro.txt

# nitro
while read -r line; do
    col1=$(echo "$line" | awk '{print $1}')
    col2=$(echo "$line" | awk '{print $2}')
    match_found=false

    for file in /scratch/user/uqpobri2/data/sponge_symbionts/nitro/*; do
        GenomeID=$(basename "$file")
        if grep -q "$col2" "$file"; then
            echo "$col1 $col2 $GenomeID" >> genomes_with_taurine_e03-100_nitro.txt
            match_found=true
            break
        fi
    done

    if [ "$match_found" = false ]; then
        echo "$col1 $col2" >> genomes_with_taurine_e03-100_nitro.txt
    fi
done < taurine_tblastx_unique_matches_e03-100_nitro.txt

# spiro
while read -r line; do
    col1=$(echo "$line" | awk '{print $1}')
    col2=$(echo "$line" | awk '{print $2}')
    match_found=false

    for file in /scratch/user/uqpobri2/data/sponge_symbionts/spiro/*; do
        GenomeID=$(basename "$file")
        if grep -q "$col2" "$file"; then
            echo "$col1 $col2 $GenomeID" >> genomes_with_taurine_e03-100_spiro.txt
            match_found=true
            break
        fi
    done

    if [ "$match_found" = false ]; then
        echo "$col1 $col2" >> genomes_with_taurine_e03-100_spiro.txt
    fi
done < taurine_tblastx_unique_matches_e03-100_spiro.txt

# thermo
while read -r line; do
    col1=$(echo "$line" | awk '{print $1}')
    col2=$(echo "$line" | awk '{print $2}')
    match_found=false

    for file in /scratch/user/uqpobri2/data/sponge_symbionts/thermo/*; do
        GenomeID=$(basename "$file")
        if grep -q "$col2" "$file"; then
            echo "$col1 $col2 $GenomeID" >> genomes_with_taurine_e03-100_thermo.txt
            match_found=true
            break
        fi
    done

    if [ "$match_found" = false ]; then
        echo "$col1 $col2" >> genomes_with_taurine_e03-100_thermo.txt
    fi
done < taurine_tblastx_unique_matches_e03-100_thermo.txt

# get unique gene/genome matches
for File in *with_taurine*; do
    # get sample id
    echo $File
    SampleID=$(echo $File | cut -d '_' -f 5 )
    echo $SampleID
    # output unique matches
    awk '{print $1, $3}' $File | sort | uniq > genomes_with_taurine_unique_matches_e03-100_${SampleID}
done

```


```bash
# note: genes lengths for tauABC
# tauA = 980
# tauB = 781
# tauC = 842
```
