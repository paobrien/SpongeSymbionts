# Endozoicomonadaceae enrichment analysis

#### Enrichment analysis of Endozoicomonadaceae genomes isolated from sponges against other environments

## Infer phylogeny to assess which genomes to include


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate gtdbtk-1.4.0

#run GTDB-tk
gtdbtk de_novo_wf \
    --genome_dir ~/Data/MAGs/Endo \
    --bacteria \
    --outgroup_taxon g__Chromatocurvus \
    --out_dir ~/Metagenomics/GTDB/Endo/No_derep \
    -x fna \
    --taxa_filter f__Endozoicomonadaceae \
    --prefix Endo \
    --cpus 4

```

## Run checkM to get genome quality for dereplication


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate checkm-genome-1.1.3

# get checkm quality
checkm lineage_wf \
    ~/Data/MAGs/Endo/ \
    ~/Metagenomics/Checkm/Endo \
    -x fna -t 16

# create checkm file
checkm qa \
    ~/Metagenomics/Checkm/Endo/lineage.ms \
    ~/Metagenomics/Checkm/Endo \
    -o 1 -f ~/Metagenomics/Checkm/Endo/checkm_endo_short.tsv --tab_table -t 1

```

## Dereplicate and remove low quality genomes


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate coverm-0.6.0

# dereplicate with coverm
coverm cluster \
    --genome-fasta-directory ~/Data/MAGs/Endo \
    -x fna \
    --ani 95 \
    --checkm-tab-table ~/Metagenomics/Checkm/Endo/checkm_endo_short.tsv \
    --output-representative-fasta-directory ~/Data/MAGs/Endo/Endo_95_75comp \
    --precluster-method finch \
    --min-completeness 50 \
    --max-contamination 10 \
    -t 16

```

## Annotate genomes


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate enrichm_0.5.0rc1

# annotate genomes using enrichm with KO, Pfam and CAZy databases
enrichm annotate \
    --output ~/Metagenomics/Enrichm/Annotate/Endo_95/ \
    --genome_directory ~/Data/MAGs/Endo/Endo_95 \
    --force \
    --ko \
    --pfam \
    --cazy \
    --threads 16 \
    --suffix .fna 
```

## Redo phylogeny with dereplicated genomes


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate gtdbtk-1.4.0

# run gtdb-tk
gtdbtk de_novo_wf \
    --genome_dir ~/Data/MAGs/Endo/Endo_95 \
    --bacteria \
    --outgroup_taxon g__Chromatocurvus \
    --out_dir ~/Metagenomics/GTDB/Endo/95_ani \
    -x fna \
    --taxa_filter f__Endozoicomonadaceae \
    --prefix Endo \
    --cpus 4


```

## Run enrichment analysis


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate enrichm_0.5.0rc1

# CAZy analysis
enrichm enrichment \
    --output ~/Metagenomics/Enrichm/Enrichment/Endo_95/CAZY \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Endo_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Endo_95/metadata_endo95_83_sponge_v_nonsponge.txt \
    --cazy \
    --force

# Pfam analysis
enrichm enrichment \
    --output ~/Metagenomics/Enrichm/Enrichment/Endo_95/PFAM \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Endo_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Endo_95/metadata_endo95_83_sponge_v_nonsponge.txt \
    --pfam \
    --force

# KO analysis
enrichm enrichment \
    --output ~/Metagenomics/Enrichm/Enrichment/Endo_95/KO \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Endo_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Endo_95/metadata_endo95_83_sponge_v_nonsponge.txt \
    --ko \
    --force

```
