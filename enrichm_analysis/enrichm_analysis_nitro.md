# Nitrosopumiliaceae enrichment analysis

#### Enrichment analysis of nitrosopumiliaceae genomes isolated from sponges against other environments

## Infer phylogeny to assess which genomes to include


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate gtdbtk-1.4.0

#run GTDB-tk
gtdbtk de_novo_wf \
    --genome_dir ~/Data/MAGs/Nitro/95_ani \
    --archaea \
    --outgroup_taxon f__Nitrosocaldaceae \
    --out_dir ~/Metagenomics/GTDB/Nitro/No_derep \
    -x fna \
    --taxa_filter f__Nitrosopumilaceae \
    --prefix Nitro \
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
    ~/Data/MAGs/Nitro/95_ani/ \
    ~/Checkm/Nitro_95 \
    -x fna -t 16

# create checkm file
checkm qa \
    ~/Checkm/Nitro_95/lineage.ms \
    ~/Checkm/Nitro_95 \
    -o 1 -f ~/Checkm/Nitro_95/checkm_nitro95_short.tsv --tab_table -t 16

```

## Dereplicate and remove low quality genomes


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate coverm-0.6.0

# dereplicate with coverm
coverm cluster \
    --genome-fasta-directory ~/Data/MAGs/Nitro/95_ani/ \
    -x fna \
    --ani 95 \
    --checkm-tab-table ~/Metagenomics/Checkm/Nitro_95/checkm_nitro95_short.tsv \
    --output-representative-fasta-directory ~/Data/MAGs/Nitro/95_ani/Dereplicated \
    --precluster-method finch \
    --min-completeness 50 \
    --max-contamination 10 \
    -t 8

```

## Annotate genomes


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate enrichm_0.5.0rc1

# annotate genomes using enrichm with KO, Pfam and CAZy databases
enrichm annotate \
    --output ~/Metagenomics/Enrichm/Annotate/Nitro_95/ \
    --genome_directory ~/Data/MAGs/Nitro/95_ani \
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
    --genome_dir ~/Data/MAGs/Nitro/95_ani/Dereplicated/Enrichment_analysis \
    --archaea \
    --outgroup_taxon f__Nitrosocaldaceae \
    --out_dir ~/Metagenomics/GTDB/Nitro/Enrichment_analysis \
    -x fna \
    --taxa_filter f__Nitrosopumilaceae \
    --prefix nitro_95 \
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
    --output ~/Metagenomics/Enrichm/Enrichment/Nitro_95/CAZY \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Nitro_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Nitro_95/metadata_nitro95_85comp_sponge_vs_nonsponge.txt \
    --cazy \
    --force

# Pfam analysis
enrichm enrichment \
    --output ~/Metagenomics/Enrichm/Enrichment/Nitro_95/PFAM \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Nitro_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Nitro_95/metadata_nitro95_85comp_sponge_vs_nonsponge.txt \
    --pfam \
    --force

# KO analysis
enrichm enrichment \
    --output ~/Metagenomics/Enrichm/Enrichment/Nitro_95/KO \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Nitro_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Nitro_95/metadata_nitro95_85comp_sponge_vs_nonsponge.txt \
    --ko \
    --force

```
