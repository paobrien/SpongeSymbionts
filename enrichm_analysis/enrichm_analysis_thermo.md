# Thermoanaerobaculaceae enrichment analysis

#### Enrichment analysis of thermoanaerobaculaceae genomes isolated from sponges against other environments

## Infer phylogeny to assess which genomes to include


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate gtdbtk-1.4.0

#run GTDB-tk
gtdbtk de_novo_wf \
    --genome_dir ~/Data/MAGs/Thermo \
    --bacteria \
    --outgroup_taxon o__Fen-336 \
    --out_dir ~/Metagenomics/GTDB/Thermo/No_derep \
    -x fna \
    --taxa_filter c__Thermoanaerobaculia \
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
    ~/Data/MAGs/Thermo/ \
    ~/Metagenomics/Checkm/Thermo \
    -x fna -t 16

# create checkm file
checkm qa \
    ~/Metagenomics/Checkm/Thermo/lineage.ms \
    ~/Metagenomics/Checkm/Thermo \
    -o 1 -f ~/Metagenomics/Checkm/Thermo/checkm_thermo_short.tsv --tab_table -t 1
```

## Dereplicate and remove low quality genomes


```bash
#!/bin/bash

# load programs
module load miniconda3
conda activate coverm-0.6.0

# dereplicate with coverm
coverm cluster \
    --genome-fasta-directory ~/Data/MAGs/Thermo \
    -x fna \
    --ani 95 \
    --checkm-tab-table ~/Metagenomics/Checkm/Thermo/checkm_thermo_short.tsv \
    --output-representative-fasta-directory ~/Data/MAGs/Thermo/Thermo_95 \
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

# annotate genomes using enrichm with KO, Pfam and CAZy databases
enrichm annotate \
    --output ~/Metagenomics/Enrichm/Annotate/Thermo_95/ \
    --genome_directory ~/Data/MAGs/Thermo/Thermo_95 \
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
    --genome_dir ~/Data/MAGs/Thermo/Thermo_95/no_gtdb \
    --bacteria \
    --outgroup_taxon o__Fen-336 \
    --out_dir ~/Metagenomics/GTDB/Thermo/95_ani \
    -x fna \
    --taxa_filter c__Thermoanaerobaculia \
    --prefix Thermo \
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
    --output ~/Metagenomics/Enrichm/Enrichment/Thermo_95/CAZY \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Thermo_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Thermo_95/metadata_thermo95_85_sponge_v_nonsponge.txt \
    --cazy \
    --force

# Pfam analysis
enrichm enrichment \
    --output ~/Metagenomics/Enrichm/Enrichment/Thermo_95/PFAM \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Thermo_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Thermo_95/metadata_thermo95_85_sponge_v_nonsponge.txt \
    --pfam \
    --force

# KO analysis
enrichm enrichment \
    --output ~/Metagenomics/Enrichm/Enrichment/Thermo_95/KO \
    --annotate_output ~/Metagenomics/Enrichm/Annotate/Thermo_95 \
    --metadata ~/Metagenomics/Enrichm/Enrichment/Thermo_95/metadata_thermo95_85_sponge_v_nonsponge.txt \
    --ko \
    --force
```
