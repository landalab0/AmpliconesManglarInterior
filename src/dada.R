getwd()
library(dada2)
#Load trim fastq files and list fastq_path content
fastq_path <- "/axolote/diana/manglares/results/02.cutadapt"
list.files(fastq_path) 

#Sort file names
Fs <- sort(list.files(fastq_path, pattern="_1.fastq"))
Rs <- sort(list.files(fastq_path, pattern="_2.fastq"))

# Extract sample names
sampleNames <- sapply(strsplit(Fs, "_1"), `[`, 1)
sampleNames
# Add complete path to remove ambiguities errors
Fs <- file.path(fastq_path, Fs)
Rs <- file.path(fastq_path, Rs)

# Quality check plot with only the first fastq file
#pdf("results/03.Dada2/CHIH_A1_QualityProfile_plot.pdf")
plotQualityProfile(c(rbind(Fs[1],Rs[1])))
#dev.off()



# Create directory for clean reads
filt_path <- file.path("results/03.Dada2" , "filter_reads") 
if(!file_test("-d", filt_path)) 
  dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

# Quality control
# V1 based on quality plots truncLen=250,200 and permisive maxEE=5,5
out1 <- filterAndTrim(Fs, filtFs, Rs, filtRs,
                      truncLen=c(250,200),
                      maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) 

#V2 luigui tutorial aprox 0 0 5,5
out2 <- filterAndTrim(Fs, filtFs, Rs, filtRs,
                      truncLen=c(0,0),
                      maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) 


#V3 based on figaro test = 280,200 maxEE=2,2
out3 <- filterAndTrim(Fs, filtFs, Rs, filtRs,
                      truncLen=c(280,210),
                      maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) 

##v4
out4 <- filterAndTrim(Fs, filtFs, Rs, filtRs,
                      truncLen=c(0,200),
                      maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) 


v1 <- as.data.frame(out1)
v2 <- as.data.frame(out2)
v3 <- as.data.frame(out3)
v4 <- as.data.frame(out4)

# Calcula el porcentaje para cada uno y añade una columna de grupo
calculate_percentage <- function(df, group_name) {
  df$percentage <- (df$reads.out / df$reads.in) * 100
  df$version <- group_name
  return(df)
}

# Aplica la función a cada dataframe y almacena los resultados
out1_with_percentage <- calculate_percentage(v1, 'out1')
out2_with_percentage <- calculate_percentage(v2, 'out2')
out3_with_percentage <- calculate_percentage(v3, 'out3')
out4_with_percentage <- calculate_percentage(v4, 'out4')

# Combina los dataframes en uno solo
combined_data <- rbind(out1_with_percentage, out2_with_percentage, out3_with_percentage, out4_with_percentage)

# Ahora puedes crear un boxplot con ggplot2 o la base de gráficos en R
library(ggplot2)
boxplot_versions <- ggplot(combined_data, aes(x = version, y = percentage, fill = version)) +
  geom_boxplot() +
  labs(x = "Filter version", y = "Percentage of reads after filter") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2")

png("results/plots/summary/04.boxplot_trunc_versions.png")
boxplot_versions
dev.off()

#Save info of final version
write.table(out1, file="results/03.Dada2/Dada_clean.tsv", quote=F, sep="\t",col.names=NA) # Table with the totals before and after cleaning


## Error Model
#De-replicate to reduce redundance 

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Add names to de-rep object
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

#Generate Error model IMPORTANT
errF <- learnErrors(derepFs, multithread=TRUE, verbose = TRUE)
errR <- learnErrors(derepRs, multithread=TRUE, verbose=TRUE)

save.image(file = "Dada2.RData") # Save point to stop for now

## ASVs inference

dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = "pseudo", verbose=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = "pseudo", verbose = TRUE)

save.image(file = "Dada2.RData") # Save point to stop for now


# Merge pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 8, verbose=TRUE)

# Create ASVs table 
seqtabAll <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtabAll)))

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtabAll, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab_nochim)
sum(seqtab_nochim)/sum(seqtabAll)

# create a new table with each ASV number and its representative sequence
PE.table_tsv_output <- seqtab_nochim
PE.table_tsv_output[PE.table_tsv_output==1]=0 # Don't consider those values that have a single observation per sample, make them 0 (sample singletons)
PE.table_tsv_output <- PE.table_tsv_output[,colSums(PE.table_tsv_output)>1] # filter singleton ASVs across the table

# Export sequences as in fasta format
uniquesToFasta(PE.table_tsv_output, fout="results/03.Dada2/ASVs.fasta", ids=paste("ASV_",1:ncol(PE.table_tsv_output), sep=""))
nochim=PE.table_tsv_output
write.table(cbind("ASVs"=1:nrow(t(PE.table_tsv_output)),"rep_seq"=rownames(t(PE.table_tsv_output))), file="results/03.Dada2/ASV_to_seqs-nochim.tsv", quote=F, sep="\t",row.names=FALSE)

# replace the rep_seq with an incremental ASV number
PE.table_tsv_output <- t(PE.table_tsv_output)
rownames(PE.table_tsv_output) <- paste0("ASV_",1:nrow(PE.table_tsv_output))

# and print the output ASV table
write.table(PE.table_tsv_output, file="results/03.Dada2/ASV_to_seqs-nochim.tsv", quote=F, sep="\t",col.names=NA)

# evaluate the total table dimensions
dim(nochim)
# [1]   64 102068
table(nchar(getSequences(nochim))) 

###Track reads lost per step ###

# By using this, we can create a function to automate this for all samples in a set:
getN <- function(x) sum(getUniques(x)) # Where getUniques gets non-repeated sequences from a dada2 object or merger object (joined reads)
track <- cbind(out1, sapply(derepFs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN), rowSums(seqtabAll), rowSums(nochim))
colnames(track) <- c("Raw", "Qual_filter", "Derep", "ASVs R1", "ASVs R2", "Merged", "nonchim")
rownames(track) <- sampleNames
write.table(track, "results/03.Dada2/Seqs_lost_in_ASVs_processing.tsv", col.names=NA, sep="\t")


# Create a quick assesment of sequences lost throughout the process
pdf("results/03.Dada2/preview_reads_passing_ASV_processing.pdf")
#png("results/plots/summary/05.ASV_secuences_process.png") #to summary
# And same thing for the percentage remaining
matplot(t(track[,-5]/track[,1]*100),type='l',xaxt='n', main="Sequences remaining after each step  - R1 (%)", xlab="Step", ylab=" Percentage of Sequences remaining")
axis(1,at=1:ncol(track[,-5]),labels=colnames(track[,-5]))
# R2
matplot(t(track[,-4]/track[,1]*100),type='l',xaxt='n', main="Sequences remaining after each step  - R2 (%)", xlab="Step", ylab=" Percentage of Sequences remaining")
axis(1,at=1:ncol(track[,-4]),labels=colnames(track[,-4]))

dev.off()

##Add final table
track2 <- data.frame(track)
track2$percentage_used <-(track2$nonchim / track2$Raw) * 100
track2
write.table(track2, "results/03.Dada2/Seqs_lost_in_ASVs_processing_percentage.tsv", col.names=NA, sep="\t")

# Save work so far
save.image(file = "Dada2.RData") 
