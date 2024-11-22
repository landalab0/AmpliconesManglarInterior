# 16S Mangrove Interior Sediments

## Project description

## Workspace

Directory structure and content information

**`data/`** Directory with raw sequences and symbolic links of sequences.

**`src/`** Scripts directory to perform each analysis

**`results/`** Contains results of each analysis.

**`outs/`** Contains the output reports of each analysis.

**`LICENSE`** File with license specifications.

**`README`** File with workflow and repository details.

## Workflow

### 01. Get sediment files

``` bash
# Create dir
mkdir -p data/{run_one,run_two}

# Create simbolic line of fastq files inside each run directory

cd data/run_one
bash get_fastq_files_run1.sh

cd data/run_two
bash get_fastq_files_run2.sh
```

### 02. Combine sequences

``` bash
bash src/01.combine_files.sh > outs/01.combine_files.out
```

### 03. Filter short reads and remove adapters

``` bash
nohup bash src/02.cutadapt.sh > outs/02.cutadapt.out &
[1] 1654676
```

### 04. Get ASVs

```bash
mkdir -p results/03.Dada2
mkdir -p results/plots
```

``` r
nohup Rscript src/03.dada.R > outs/03.dada2.out 2> outs/03.dada2.error &
[1] 1546439

```

```bash
#save PID
echo $! > outs/03.dada2_PID.txt
#check progress
ps -ef | grep Rscript
ps -f -p $(cat outs/03.dada2_PID.txt)
```

### 05. Taxonomy assignment

```bash
mkdir -p results/04.qiime

#Import dada R sequences and ASV table to QIIME2
bash src/04.import_data_to_qiime2.sh 1> outs/04.import.out

#Taxonomy classification
nohup bash src/05.tax_assign.sh 1> outs/05.tax_assign.out 2> outs/05.tax_assign.error &
#check progress
echo $! > outs/05.tax_PID.txt
ps -f -p $(cat outs/05.tax_PID.txt) 
```

### 06. Filters

```bash
bash src/06.filters.sh 1> outs/06.filters.out
```

### 07. Phylogeny

```bash
nohup bash src/07.phylogeny.sh 1> outs/07.phylogeny.out 2> outs/07.phylogeny.error &
echo $! > outs/07.phylogeny_PID.txt
ps -f -p $(cat outs/07.phylogeny_PID.txt)
```

### 08. Exports


```bash

bash src/08.exports.sh 1> outs/08.exports.out

```

### 09. Explore Diversity

check `**src/Interior_Mangroove_Sediments**`

