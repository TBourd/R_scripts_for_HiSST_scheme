#'SNP binary matrix R script for HiSST scheme


#'This R script provides binary matrix based on SNP at each nucleotide position of all sequence type 
#' compared to the first ST1 (i.e. ST of the reference strain S. marcescens DB11).
#'It helps to see the nucleotide polymorphisms of each HiSST locus and their relationship 
#' visualized by minimum spanning trees.
#' 
#'This R script uses dplyr package
#'
#'The data frame used for this script provided from BioEdit software (http://www.mbio.ncsu.edu/BioEdit/bioedit.html),
#' using the « True sequence positions of selected sequences » export option. 
#'Example files 'Seq_bssA.txt','Seq_gabR.txt', or 'Seq_dhaM.txt' adapted for this script are available 
#' at https://github.com/TBourd/Rscript_for_HiSST_scheme_of_S_marcescens. 



library(dplyr)

setwd("~/SNP_bssA_All_ST") #Set the working directory path of your project

Seq <- read.delim("~/Seq_bssA.txt",skip = 1)

Seq <- Seq %>% select(-contains('Pos')) #delete columns containing "Pos" characters
Seq <- Seq[, colSums(is.na(Seq)) != nrow(Seq)] #delete columns with NA value


# Function to create SNP binary matrix for a sequence type, based on the first column (ST 1)
rbind.FUN <- function(ST){
  SNPbin <- NULL
  l=0 #Row
  m=ST  #Column
  repeat{
    for (i in unique(1:nrow(Seq))){
      l = l+1
      tmp <- ifelse(Seq[l,m]==Seq[l,1], 1, 0)
      SNPbin <- rbind(SNPbin,tmp)
    }
    if (l == length(Seq[,1])) {break}
  }
  return(SNPbin)
}

# Compile SNP binary matrix of all STs
SNPbin.tot <- NULL
ST = 0
repeat {
  for (i in unique(1:ncol(Seq))){
    ST = ST+1
    SNPbin.tot <- cbind (SNPbin.tot,rbind.FUN(ST))
  }
  if (ST == length(Seq[1,])) {break}
}

# Reverse rows and columns with Transpose
inv_SNPbin <- t(SNPbin.tot)

# Customize columns and rows names
colnames(inv_SNPbin) <- paste("Pos",1:ncol(inv_SNPbin),sep = "_")
inv_SNPbin=data.frame(ST=paste("ST",1:nrow(inv_SNPbin)),inv_SNPbin)

# Save ANIb binary matrix
write.table(inv_SNPbin,file="SNPbin.txt", sep = "\t",row.names = FALSE)



