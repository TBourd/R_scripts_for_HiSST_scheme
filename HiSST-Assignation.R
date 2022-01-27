#'HiSST-Assignation


#'This R script allows to identify the HiSST-profile ID of the corresponding S. marcecens isolate(s),
#'by using the sequence type (ST) ID of each HiSST-locus (i.e. bssA, gabR and dhaM) of the isolate(s).
#'
#'The last section of this script allows to add a new HiSST profile to the database.
#'
#'The file 'HiSST_Database.txt' needed for this script is available 
#'at https://github.com/TBourd/R_scripts_for_HiSST_scheme/tree/main/Data.


  
setwd("~/HiSST_profiles") #Set the working directory path of your project
  
HiSST <- read.table("~/v01-2022_HiSST_database.txt",header = TRUE) #Import Database 


#'Combine loci-ST profiles of data base in 1 column  
combination_ST <- mapply(paste, sep = "-", HiSST$bssA, HiSST$gabR,HiSST$dhaM)
  
HiSST.cb <- data.frame(HiSST,combination_ST)

## -- Import data --

#### Without input file ####
#ID <- c("a","b","c") #'Enter the ID of the strain(s) : e.g. strain #a = "a"
#query <- data.frame(ID=ID,bssA=c(2,4,3),gabR=c(1,7,1),dhaM=c(1,4,1)) #'Enter ST of each HiSST-locus

#### Import .txt file ####
# These tables are obtained from the R script "dada2_for_HiSST_scheme.R" (lines 242-266).

bssA <- read.table("~/bssA_Samples_and_ASV-ST.txt",header = TRUE)
dhaM <- read.table("~/dhaM_Samples_and_ASV-ST.txt",header = TRUE)
gabR <- read.table("~/gabR_Samples_and_ASV-ST.txt",header = TRUE)

bssA <- data.frame(ID=bssA$Samples, bssA=bssA$ST)
bssA$ID <- gsub("-bssA", "", bssA$ID)
bssA$bssA <- gsub("ST", "", bssA$bssA)
dhaM <- data.frame(ID=dhaM$Samples, dhaM=dhaM$ST)
dhaM$ID <- gsub("-dhaM", "", dhaM$ID)
dhaM$dhaM <- gsub("ST", "", dhaM$dhaM)
gabR <- data.frame(ID=gabR$Samples, gabR=gabR$ST)
gabR$ID <- gsub("-gabR", "", gabR$ID)
gabR$gabR <- gsub("ST", "", gabR$gabR)

library(dplyr)
query <- full_join(bssA, dhaM, by = "ID")
query <- full_join(query, gabR, by = "ID")

#### Import .txt file ####
#'with the following format : Col = 'ID, bssA, gabR, dhaM' and Rows = loci-ST profiles of strains 
#strain <- read.table("~/HiSST_profiles/strain.txt",header = TRUE) #Replace the directory path by your file such as described above
#ID <- strain$ID
#query <- data.frame(ID=ID,bssA=strain$bssA,gabR=strain$gabR,dhaM=strain$dhaM)

## -- End Import data ----


#'Combine loci-ST profiles of strains in 1 column
combin_query <-setNames(as.data.frame(mapply(paste, sep = "-",query$bssA,query$gabR,query$dhaM)),
                             nm="combination_ST")
row.names(combin_query) <-1:nrow(combin_query)
df_query <- data.frame(ID=ID,combination_ST=combin_query)

#' Assignation of HiSST ID for each query Strains with the data base
assign_HiSST_profile=merge(df_query,HiSST.cb, by=c("combination_ST"),all.x = TRUE)

Strain.ST <- data.frame(ID=assign_HiSST_profile$ID,HiSST=assign_HiSST_profile$HiSST,
                        HiSST.profile=assign_HiSST_profile$combination_ST)
#Strain.ST <- with(Strain.ST,  Strain.ST[order(ID) , ]) #Sort data frame by ID 

write.table(Strain.ST,file="HiSST ID of strains.txt", sep="\t",row.names = FALSE, quote = F)


#### Add new HiSST profile to the database ####

#run this function
add.newST <- function(StrainID){
  a=query[grepl(StrainID,query$ID),]
  a$HiSST <- nrow(HiSST)+1
  HiSST.new <- merge(HiSST,a,all=TRUE)
  HiSST.new <- subset(HiSST.new,select = -ID)
  return(HiSST.new)
}

#'Write the Strain ID into add.newST function 
#'to add the new HiSST profile in the HiSST database

new.HiSST =add.newST("1-BB24") # Change "a" with the strain ID
write.table(new.HiSST,file="New HiSST database.txt", sep="\t",row.names = FALSE)



#_______________________________

####_______Add more new ST ________ ####

more.newST <- function(StrainID){
  a=query[grepl(StrainID,query$ID),]
  a$HiSST <- nrow(new.HiSST)+1
  HiSST.new <- merge(new.HiSST,a,all=TRUE)
  HiSST.new <- subset(HiSST.new,select = -ID)
  return(HiSST.new)
}

### Repeat for each new ST :

# STEP 1:____

#'Combine new loci-ST profiles of data base in 1 column  
combination_ST <- mapply(paste, sep = "-", new.HiSST$bssA, new.HiSST$gabR,new.HiSST$dhaM)

new.HiSST.cb <- data.frame(new.HiSST,combination_ST)
#' Assignation of HiSST ID for each query Strains with the NEW data base
new.HiSST.profile=merge(df_query,new.HiSST.cb, by=c("combination_ST"),all.x = TRUE)

new.Strain.ST <- data.frame(ID=new.HiSST.profile$ID,HiSST=new.HiSST.profile$HiSST)
#Strain.ST <- with(Strain.ST,  Strain.ST[order(ID) , ]) #Sort data frame by ID 


# STEP 2:____

##/!\ Check the "new.Strain.ST" file to find which strain doesn't have an assigned ST


# STEP 3:____

new.HiSST =more.newST("strainID") # Change "strainID" with the strain ID which doesn't have an assigned ST

### /!\ Repeat STEP 1 and 2 to check if all of strain have an assigned ST: 
# if not, repeat STEP 3, then STEP 1 and 2


write.table(new.HiSST,file="New HiSST database.txt", sep="\t",row.names = FALSE)
