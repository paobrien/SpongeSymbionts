# Microtrichaceae enrichment analysis

#### Enrichment analysis of microtrichaceae genomes isolated from sponges against other environments

## Infer phylogeny to assess which genomes to include


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate gtdbtk-1.4.0

#run GTDB-tk
gtdbtk de_novo_wf \
    --genome_dir ~/Data/MAGs/Micro \
    --bacteria \
    --outgroup_taxon g__Microthrix \
    --out_dir ~/Metagenomics/GTDB/Micro/No_derep \
    -x fa \
    --taxa_filter f__UBA11606,f__Bin134,f__TK06 \
    --prefix Micro \
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
    ~/Data/MAGs/Micro/ \
    ~/Metagenomics/Checkm/Micro \
    -x fna -t 16

# Create checkm file
checkm qa \
    ~/Metagenomics/Checkm/Micro/lineage.ms \
    ~/Metagenomics/Checkm/Micro \
    -o 1 -f ~/Metagenomics/Checkm/Micro/checkm_micro_short.tsv --tab_table -t 1
```

## Dereplicate and remove low quality genomes


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate coverm-0.6.0

# dereplicate with coverm
coverm cluster \
    --genome-fasta-directory ~/Data/MAGs/Micro \
    -x fna \
    --ani 95 \
    --checkm-tab-table ~/Metagenomics/Checkm/Micro/checkm_micro_short.tsv \
    --output-representative-fasta-directory ~/Data/MAGs/Micro/Micro_95 \
    --precluster-method finch \
    --min-completeness 85 \
    --max-contamination 10 \
    -t 1
```

## Annotate genomes


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate enrichm_0.5.0rc1

# annotate genomes using enrichm with KO, Pfam and CAZy databasesn
enrichm annotate \
    --output ~/Metagenomics/Enrichm/Annotate/Micro_95/ \
    --genome_directory ~/Data/MAGs/Micro/Micro_95 \
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
    --genome_dir ~/Data/MAGs/Micro/Micro_95/no_gtdb \
    --bacteria \
    --outgroup_taxon g__Microthrix \
    --out_dir ~/Metagenomics/GTDB/Micro/95_ani \
    -x fna \
    --taxa_filter f__UBA11606,f__Bin134,f__TK06 \
    --prefix Micro \
    --cpus 4
```

## Run enrichment analysis


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate enrichm_0.5.0rc1

# KO analysis
enrichm enrichment \
    --output ~/Metagenomics/Enrichm/Enrichment/Micro_95/KO \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Micro_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Micro_95/metadata_micro95_85_sponge_v_nonsponge.txt \
    --ko \
    --force

# Pfam analysis
enrichm enrichment \
    --output ~/Metagenomics/Enrichm/Enrichment/Micro_95/PFAM \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Micro_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Micro_95/metadata_micro95_85_sponge_v_nonsponge.txt \
    --pfam \
    --force

# CAZy analysis
enrichm enrichment \
    --output ~/Metagenomics/Enrichm/Enrichment/Micro_95/CAZY \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Micro_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Micro_95/metadata_micro95_85_sponge_v_nonsponge.txt \
    --cazy \
    --force
```
