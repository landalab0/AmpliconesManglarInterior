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

