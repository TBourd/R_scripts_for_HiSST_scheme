#'HiSST-Assignation


#'This R script allows to identify the HiSST-profile ID of the corresponding S. marcecens isolate(s),
#'by using the sequence type (ST) ID of each HiSST-locus (i.e. bssA, gabR and dhaM) of the isolate(s).
#'
#'The last section of this script allows to add a new HiSST profile to the database.
#'
#'The file 'HiSST_Database.txt' needed for this script is available 
#'at https://github.com/TBourd/R_scripts_for_HiSST_scheme/tree/main/Data.


  
setwd("~/HiSST_profiles") #Set the working directory path of your project
  
HiSST <- read.table("~/HiSST_Database.txt",header = TRUE) #Import Database 


#'Combine loci-ST profiles of data base in 1 column  
combination_ST <- mapply(paste, sep = "-", HiSST$bssA, HiSST$gabR,HiSST$dhaM)
  
HiSST.cb <- data.frame(HiSST,combination_ST)

## -- Import data --

#### Without input file ####
ID <- c("a","b","c") #'Enter the ID of the strain(s) : e.g. strain #a = "a"
query <- data.frame(ID=ID,bssA=c(2,4,3),gabR=c(1,7,1),dhaM=c(1,4,1)) #'Enter ST of each HiSST-locus

#### Import .txt file ####
#'with the following format : Col = 'ID, bssA, gabR, dhaM' and Rows = loci-ST profiles of strains 
strain <- read.table("~/HiSST_profiles/strain.txt",header = TRUE) #Replace the directory path by your file such as described above
ID <- strain$ID
query <- data.frame(ID=ID,bssA=strain$bssA,gabR=strain$gabR,dhaM=strain$dhaM)

## -- End Import data ----


#'Combine loci-ST profiles of strains in 1 column
combin_query <-setNames(as.data.frame(mapply(paste, sep = "-",query$bssA,query$gabR,query$dhaM)),
                             nm="combination_ST")
row.names(combin_query) <-1:nrow(combin_query)
df_query <- data.frame(ID=ID,combination_ST=combin_query)

#' Assignation of HiSST ID for each query Strains with the data base
assign_HiSST_profile=merge(df_query,HiSST.cb, by=c("combination_ST"),all.x = TRUE)

Strain.ST <- data.frame(ID=assign_HiSST_profile$ID,HiSST=assign_HiSST_profile$HiSST)
#Strain.ST <- with(Strain.ST,  Strain.ST[order(ID) , ]) #Sort data frame by ID 

write.table(Strain.ST,file="HiSST ID of strains.txt", sep="\t",row.names = FALSE)


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

new.HiSST =add.newST("a") # Change "a" with the strain ID
write.table(new.HiSST,file="New HiSST database.txt", sep="\t",row.names = FALSE)
