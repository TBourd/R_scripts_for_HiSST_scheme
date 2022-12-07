#'dada2 for S. marcescens HiSST scheme : Illumina-sequences Management and Quality Control


#'This script originate from dada2 processing (https://benjjneb.github.io/dada2/index.html)
#'and was adapted for HiSST analysis .
#'
#'The following script is an example of raw sequencing reads processing for dhaM locus 
#'included primer sequences removal with the software Cutadapt v. 2.10, 
#'followed by quality control, paired ends merging and chimera check 
#'using the default parameters specified in the package dada2 v1.8.0 
#'that include packages ShortRead v1.48.0 and Biostrings v2.58.0
#'



library(dada2)
library(ShortRead)
library(Biostrings)

path = setwd("~/R/Sequencing_Results/Serratia_marcescens/dhaM/") #Set the working directory path of your project

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_RXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

#### Identify Primers ####
FWD = "GGCGTCCAGCATYGCCTT"
REV = "GACGTGCGCGACATGCTG"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#"pre-filter" the sequences just to remove those with Ns
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#count the number of times the primers appear
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#### Remove Primers ####

system2("cutadapt", args = "--version") # Run on Linux system

#cutadapt <- "C:/Users/bourd/AppData/Roaming/Python/Python39/Scripts/cutadapt.exe" # Run on Windows
#system2(cutadapt, args = "--version") # Run on Windows

#create output filenames for the cutadapt-ed files, and define the parameters we are going to give the cutadapt command.
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt on linux system - (On Windows: remove quotation marks on "cutadapt" in the function 'system2')
for(i in seq_along(fnFs)) {
  system2("cutadapt", args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               "--discard-untrimmed", #Discard reads in which no adapter was found
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


#_________________________________________________________#

# Remove reads with more than 'maxLen = x' nucleotides #
#fnFs.cutFilt <- file.path(path, "cutFilt", basename(fnFs.cut)) # Put N-filterd files in filtN/ subdirectory
#fnRs.cutFilt <- file.path(path, "cutFilt", basename(fnRs.cut))
#filterAndTrim(fnFs.cut, fnFs.cutFilt, fnRs.cut, fnRs.cutFilt, maxLen = 200, multithread = FALSE)

#path.cutFilt <- file.path(path, "cutFilt")

#_______________________________________________________#

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
cutFs.name=sub("S502.","x.",cutFs)
cutFs.name=sub("S503.","x.",cutFs.name)
cutFs.name=sub("S505.","x.",cutFs.name)
cutFs.name=sub("S506.","x.",cutFs.name)
cutFs.name=sub("S507.","x.",cutFs.name)
cutFs.name=sub("S510.","x.",cutFs.name)
cutFs.name=sub("S511.","x.",cutFs.name)
cutFs.name=sub("S513.","x.",cutFs.name)
cutFs.name=sub("S515.","x.",cutFs.name)
cutFs.name=sub("S516.","x.",cutFs.name)
cutFs.name=sub("S517.","x.",cutFs.name)
cutFs.name=sub("S518.","x.",cutFs.name)
cutFs.name=sub("S520.","x.",cutFs.name)
cutFs.name=sub("S521.","x.",cutFs.name)
cutFs.name=sub("S522.","x.",cutFs.name)

get.sample.name <- function(fname) strsplit(basename(fname), "x.")[[1]][2]
sample.names <- unname(sapply(cutFs.name, get.sample.name))
head(sample.names)


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(sample.names), "_R"), `[`, 1)



#### dada2 ####


plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))



out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(220,170),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)


# Learn the Error Rates
# Tool to visualize the frequency of error rate as a function of quality score.
# Necessary for the algortith - see the paper
errF <- learnErrors(filtFs, multithread=TRUE, MAX_CONSIST = 20)
errR <- learnErrors(filtRs, multithread=TRUE, MAX_CONSIST = 20)

plotErrors(errF, nominalQ=TRUE)


# Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding 
# “abundance” equal to the number of reads with that unique sequence. Dereplication substantially reduces 
# computation time by eliminating redundant comparisons.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#We are now ready to apply the core sample inference algorithm to the dereplicated data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#dadaFs[[1]]
# dadaRs[[1]]


# Merging paired ends
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
# head(mergers[[1]])

# We can now construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Optional : remove unexpected sequences regarding amplicon size
seqtab.filt <- seqtab[,nchar(as.character(colnames(seqtab)))==243] #Change value by the expected amplicon size whithout primers
dim(seqtab.filt)
table(nchar(getSequences(seqtab.filt)))

# Remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab.filt, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file = "Summary.csv")

taxa <- assignSpecies(seqtab.nochim, "ST_dhaM.fa", allowMultiple = TRUE)#Taxonomic assignment to the ST level by exact matching
colnames(taxa) <- c("Species", "ST")#Change column names for ST assignment

# taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa")

## making and writing out standard output files:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV_dhaM", i, sep="_")
}

# fasta:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "dhaM_ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "dhaM_ASVs_counts.txt", sep="\t", quote=F)

# taxa table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "dhaM_ASVs_taxonomy.txt", sep="\t", quote=F)



### _____ Export results in a matrix of samples (rows) and AVS counts (columns) ___ ####

#'The output is used in the script "Rscript_Step1_Create_dendrogram_and_binary-matrix.R"
#' available at: https://github.com/TBourd/R_scripts_HiSST_SM-outbreaks/tree/main/Step1_Create_dendrogram_and_binary-matrix


ASV.sample <- seqtab.nochim
colnames(ASV.sample) <- sub(">", "", asv_headers)

ASV.sample.sort <- ASV.sample
for (i in ncol(ASV.sample.sort):1) {
  ASV.sample.sort <- ASV.sample.sort[order(-ASV.sample.sort[,i]),]
}

ASV.sample.sort <- data.frame(Samples = row.names(ASV.sample.sort),ASV.sample.sort,row.names = NULL)

write.table(ASV.sample.sort, "dhaM_ASV_samples.txt", sep="\t", quote=F)



#### ______ Export a formatted taxon table _______ ####

#'The output is used in the script "ST-Assignation_Clinic-eDNA.R"
#' available at: https://github.com/TBourd/R_scripts_HiSST_SM-outbreaks/tree/main/Step2_Assign-ST-infos

ASV.sample.mat <- seqtab.nochim
colnames(ASV.sample.mat) <- sub(">", "", asv_headers)
name.ST.ASV <- data.frame(ST = asv_tax[,2], ASV = row.names(asv_tax), row.names = NULL)

# Change column names without any ST-ID by ASV names 
name.STASV <- name.ST.ASV
for (i in 1:nrow(name.STASV)) {
  ifelse(is.na(name.STASV[i,1]),name.STASV[i,1]<-name.STASV[i,2],name.STASV[i,1])
}
colnames(ASV.sample.mat) <- name.STASV$ST

# Compute a data frame with Sample names, and ST-ID or ASV-ID
df <- data.frame(Samples=rownames(ASV.sample.mat), ST=nrow(ASV.sample.mat))
for (j in 1:nrow(ASV.sample.mat)) {
  max(ASV.sample.mat[j,]) -> x
  for (i in 1:ncol(ASV.sample.mat)) {
    ifelse(ASV.sample.mat[j,i] == x,
           z <- colnames(ASV.sample.mat)[i],z <- "NA")
    ifelse(z == "NA",NA, df$ST[j] <- z)
  }
}

write.table(df, "dhaM_Samples_and_ASV-ST.txt", sep="\t", quote=F, row.names = F)