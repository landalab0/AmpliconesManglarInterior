# Get path info
getwd()

# Load packages
library(Rcpp)
library(dada2)
library(ggplot2)
library(gridExtra)

##### 01. Prepare files -------------------------------------------------------

#Load trim fastq files and list fastq_path content
fastq_path <- "/axolote/diana/manglares/v0.2/results/02.cutadapt"
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

##### 02. Check Quality --------------------------------------------------------

# Quality check plot with only the first fastq file
QC_F1_15 <- plotQualityProfile(Fs[1:15], aggregate = TRUE)
QC_R1_15 <- plotQualityProfile(Rs[1:15], aggregate = TRUE)
QCsFR1_15 <- grid.arrange(QC_F1_15, QC_R1_15, nrow = 1)

#save in png format
ggsave("results/plots/01.QualityProfile_1-15.png", QCsFR1_15, width = 7, height = 3)

#save in pdf format
QCsFR1_15 <- grid.arrange(QC_F1_15, QC_R1_15, nrow = 1)
ggsave("results/plots/01.QualityProfile_1-15.pdf", QCsFR1_15, width = 7, height = 3)

##### 03.Quality control -------------------------------------------------------

# Create directory for clean reads
filt_path <- file.path("results/03.Dada2" , "01.filter_reads") 
if(!file_test("-d", filt_path)) 
  dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

# Filter versions

# V1 
out1 <- filterAndTrim(Fs, filtFs, Rs, filtRs,
                      truncLen=c(250,200),
                      maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) 
head(out1)

#V2 extra permissive
out2 <- filterAndTrim(Fs, filtFs, Rs, filtRs,
                      maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) 
head(out2)

#V3
out3 <- filterAndTrim(Fs, filtFs, Rs, filtRs,
                      truncLen=c(280,200),
                      maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) 
head(out3)

##v4
out4 <- filterAndTrim(Fs, filtFs, Rs, filtRs,
                      truncLen=c(0,200),
                      maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) 
head(out4)

##v5
out5 <- filterAndTrim(Fs, filtFs, Rs, filtRs,
                      truncLen=c(260,200),
                      maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE)
head(out5)

## compare trunc versions
v1 <- as.data.frame(out1)
v2 <- as.data.frame(out2)
v3 <- as.data.frame(out3)
v4 <- as.data.frame(out4)
v5 <- as.data.frame(out5)


# Percentage function
calculate_percentage <- function(df, group_name) {
  df$percentage <- (df$reads.out / df$reads.in) * 100
  df$version <- group_name
  return(df)
}

# Get percentage
out1_with_percentage <- calculate_percentage(v1, 'v1:250-200')
out2_with_percentage <- calculate_percentage(v2, 'v2:0-0')
out3_with_percentage <- calculate_percentage(v3, 'v3:280-200')
out4_with_percentage <- calculate_percentage(v4, 'v4:0-200')
out5_with_percentage <- calculate_percentage(v5, 'v5:260-200')

# Combine percentage versions
combined_data <- rbind(out1_with_percentage, out2_with_percentage, 
                       out3_with_percentage, out4_with_percentage,
                       out5_with_percentage)

# Compare plot
boxplot_versions <- ggplot(combined_data, aes(x = version, y = percentage, 
                    fill = version)) + geom_boxplot() + theme_bw() +
  labs(x = "Filter version", y = "Percentage of reads after filter") +
  scale_fill_brewer(palette = "Set2")

boxplot_versions

#save plot as png
ggsave("results/plots/02.boxplot_trunc_versions.png", boxplot_versions, width = 6)

#save plot as pdf
ggsave("results/plots/02.boxplot_trunc_versions.pdf", boxplot_versions, width = 6)

#Save info of final version
#We chose v5 
write.table(out5_with_percentage, file="results/03.Dada2/Dada_clean_reads.tsv", quote=F, sep="\t",col.names=NA) # Table with the totals before and after cleaning

##### 04.Error Model -----------------------------------------------------------

#De-replicate to reduce redundance 

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Add names to de-rep object
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

#Generate Error model IMPORTANT
errF <- learnErrors(derepFs, multithread=TRUE, verbose = TRUE)
errR <- learnErrors(derepRs, multithread=TRUE, verbose=TRUE)

save.image(file = "src/Dada2.RData") # Save point to stop for now

##### 0.5 Get ASVs -------------------------------------------------------------
# ASVs inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = "pseudo", verbose=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = "pseudo", verbose = TRUE)

save.image(file = "src/Dada2.RData") # Save point to stop for now

# Merge pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 8, verbose=TRUE)

# Create ASVs table 
seqtabAll <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtabAll)))

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtabAll, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab_nochim)
sum(seqtab_nochim)/sum(seqtabAll)

##### 0.6 info -----------------------------------------------------------------

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

table(nchar(getSequences(nochim))) 

###Track reads lost per step ###

# By using this, we can create a function to automate this for all samples in a set:
getN <- function(x) sum(getUniques(x)) # Where getUniques gets non-repeated sequences from a dada2 object or merger object (joined reads)
track <- cbind(out1, sapply(derepFs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN), rowSums(seqtabAll), rowSums(nochim))
colnames(track) <- c("Raw", "Qual_filter", "Derep", "ASVs R1", "ASVs R2", "Merged", "nonchim")
rownames(track) <- sampleNames
write.table(track, "results/03.Dada2/Seqs_lost_in_ASVs_processing.tsv", col.names=NA, sep="\t")


# Create a quick assesment of sequences lost throughout the process
png("results/plots/03.preview_reads_passing_ASV_processing.png")
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
save.image(file = "src/Dada2.RData") 
