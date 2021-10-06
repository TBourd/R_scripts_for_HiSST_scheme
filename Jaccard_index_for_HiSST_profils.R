#'Jaccard distance for HiSST profils.
#'
#'See an example in the "Figure 6: Survey of Serratia marcescens in sink drains of a NICU.", in:
#' Bourdin T, Monnier A, Benoit MÈ, Bédard E, Prévost M, Quach C, Déziel E, Constant P. 
#' A High-Throughput Short Sequence Typing Scheme for Serratia marcescens Pure Culture 
#' and Environmental DNA. Appl Environ Microbiol. 2021 Sep 29:AEM0139921. doi: 10.1128/AEM.01399-21.


#'This script creates an UPGMA dendrogram based on Jaccard distance 
#'computed with the HiSST profile of gabR, bssA and dhaM loci amongst several environmental samples
#'
#'This script use the following packages : 
#', SnowballC, dplyr, proxy, dendextend. 


setwd("~/HiSST")

Compil_ST <- read.delim("bssA_ASV_samples.txt") #dada2 output : a matrix of samples (rows) and AVS counts (columns), see "dada2_for_HiSST_scheme.R" script.


## Delete query samples
#seqtab.nochim.env <- seqtab.nochim[ grep("001-", row.names(seqtab.nochim), invert = TRUE) , ] # Delete sample that contain '001-'

library(SnowballC)
library(dplyr)

# Create a function to merge replicate samples
Comp.merge <- function(x,df) {
  c <- (filter(Compil_ST,grepl(x,Compil_ST$Samples)))
  c <- as.data.frame(c(Samples=x,c[,2:ncol(c)]))
  Comp.comb <- rbind(df,c)
}

df <- data.frame()
df <- Comp.merge("sample-name_1",df)
df <- Comp.merge("sample-name_2",df)
df <- Comp.merge("sample-name_3",df)

# For sequential sample names. Ex.: "sample_01", "sample_02", ..., "sample_30" :
ls.BB <- c(sprintf("sample-name_%02d", seq(1,30)))
for (i in ls.BB) {
  a <- i
  df <- Comp.merge (a,df)
}


# Agregate with min, max or mean, and then compute a binary matrix
Compil.ST <- aggregate(x=df[,2:ncol(df)], by=list(Samples=df$Samples), 
                       min, na.rm = TRUE) #min ou max ou mean ?
Compil.ST.bin <- as.data.frame(Compil.ST[,2:ncol(Compil.ST)], row.names = Compil.ST[,1])
Compil.ST.bin[Compil.ST.bin > 0] <- 1

# OR
#'Aggregate with a cutoff and binary output :
Compil.ST <- aggregate(x=df[,2:ncol(df)], by=list(Samples=df$Samples),
                       function(x) ifelse(mean(x)<50,0,1)) # Ex. : Cutoff 50



## /!\ Change bssA or dhaM or gabR /!\
Compil.bssA <- Compil.ST.bin
#Compil.gabR <- Compil.ST.bin
#Compil.dhaM <- Compil.ST.bin

Compil.ALL <- cbind(Compil.bssA,Compil.gabR,Compil.dhaM)

### Jaccard Index #### 

#install.packages("proxy")


require(proxy)

Jac_dist = proxy::dist(Compil.ALL, by_rows = TRUE, method = "Jaccard")


#### Dendro ####
library(dendextend)


d2 <- hclust(Jac_dist,method = "average")



d1 = as.dendrogram(d2)


par(mar=c(2,2,2,3),mgp=c(1,1,0),cex.axis=0.5)
d1 = d1 %>% 
  set("labels_col", h=0.8) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color", h=0.8) %>%
  set("leaves_pch", 19) %>%
  set("leaves_col", "dark cyan") %>%
  set("labels_cex",0.75) 

plot(d1,horiz = TRUE, yaxt="n",)
axis(1,cex.axis=0.5, ylim = c(0,1))

title(main = "Sequence type profile of sinks")
#title(sub="Sink ID")

#click_rotate (d1)


# Put a line to the desired height #
abline(h=0.8,col="red")


d2$height= ifelse(d2$height< 0.01,NA,d2$height) # Condition of which bootstrap height to show

# Write bootstrap values on the dendrogram #
with(pvclust:::hc2axes(d2), 
     text(x.axis, y.axis, round(y.axis, 2), adj = c(.5, 1),cex = 0.3)
)

