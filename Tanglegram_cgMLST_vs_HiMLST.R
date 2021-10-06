#'Customize the Tanglegram of HiSST scheme and whole genome sequences.


#'This script was inspired from URL https://www.r-graph-gallery.com/340-custom-your-dendrogram-with-dendextend
#'and has been adapted for the development of the HiSST scheme, which helps to select the loci of the scheme.
#'
#'The following script draws dendrograms based on the ANIb score of concatenated loci 
#'selected for HiSST scheme and core genome MLST of WG to discriminate strains of Serratia marcescens;
#'using pvclust package for Hierchical Clustering.
#'Dendextend and tidyverse packages were used for the customized Tanglegram.
#'
#'It uses the files called "ANIb_S_marcescens_gabR-bssA-dhaM.tab" and "ANIb_WG_S_marcescens.tab",
#'available at https://github.com/TBourd/R_scripts_for_HiSST_scheme/tree/main/Data.

library(pvclust)

setwd("~/R/cgMLST_vs_HiMLST") #Set the working directory path of your project



#### PHYLOGENOMICS HiSST SNPs vs cgMLST ####

#### HiSST ##

HiSST.dist <- readDist("HiSST_diff-matrix.tab",format = "phylip") # Import a distance matrix of concatenated loci

compH <- hclust(HiSST.dist,method = "complete")

#library('TreeTools')
#setdiff(TipLabels(tree1), TipLabels(tree2)) # In tree1 but not tree2
labels(compH)=gsub("Serratia_marcescens_","",labels((compH)))
labels(compH)=gsub("_"," ",labels(compH)) #Delete '_' character in the sample names


#### WG ##

cgMLST.dist <- readDist("cgMLST_WGS.tab",format = "phylip") # Import a distance matrix based on the core-genome MLST, compute from http://wgmlstdb.imst.nsysu.edu.tw/index.php

compG <- hclust(cgMLST.dist,method = "complete")

labels(compG)=gsub("Serratia_marcescens_","",labels((compG)))
labels(compG)=gsub("_"," ",labels(compG)) #Delete '_' character in the sample names


#___________ Statistical test ______________####

# Robinson-Foulds distances 

#see Smith 2020, https://doi.org/10.1093/bioinformatics/btaa614

Tree1 = as.phylo(compH)

Tree2 = as.phylo(compG)

TreeDistance(Tree1,Tree2) #Calculate tree dissimilarity based generalized Robinson–Foulds distances 

RF.dist(Tree1,Tree2,normalize = TRUE)


#'  >>> Contingency table comparison using Pearson's Chi-squared Test for Count Data <<<<
#'  
#'  The chi-square test of independence is used to analyze the frequency table (i.e. contengency table) 
#'  formed by two categorical variables. The chi-square test evaluates whether there is a significant association 
#'  between the categories of the two variables. 
HiSST.dist <- as.dist(HiSST.dist)
cgMLST.dist <- as.dist(cgMLST.dist)
HiSST.dist.norm <-decostand(HiSST.dist, method="norm")
cgMLST.dist.norm <-decostand(cgMLST.dist, method="norm")

HiSST.km <- kmeans(HiSST.dist.norm,centers = 6,nstart = 10000)

cgMLST.km <- kmeans(cgMLST.dist.norm,centers = 6,nstart = 10000)

table(as.factor(cgMLST.km$cluster),as.factor(cgMLST.km$cluster))


chisq.test(as.factor(HiSST.km$cluster),as.factor(cgMLST.km$cluster),
           simulate.p.value = TRUE, B=10000) # H0 : the two typology are independent, no relationship between variable A and B


# Optimal number of clusters according to silhouette widths
library(cluster)

asw <- numeric(64) # In numeric(), put the number of organisms inputted

for (k in 2:(64-1)) { # Replace '64' by the corresponding number of organisms inputted
  sil <- silhouette(cutree(compG,k=k),cgMLST.dist)
  asw[k] <- summary(sil)$avg.width
}

k.best <- which.max(asw)
plot(asw, type = "h", xlab = "k", ylab = "Average silhouette width")
axis(1, k.best,paste("optimum",k.best,sep = "\n"),col = "red", col.axis="red")
points(k.best,max(asw),col="red")




#### Compilation ####

library(dendextend)

library(tidyverse)


d1 = as.dendrogram(compH)

d2 = as.dendrogram(compG)


dend12_corrected_1 <- untangle_random_search(d1,d2) #random shuffling the two untangled dendrogram and each time checking if their entanglement was improved
dend12_corrected_2 <- untangle_step_rotate_2side(dend12_corrected_1[[1]],dend12_corrected_1[[2]]) #rotate the tree and looking for a better match. Use it if you think the two trees could better match.


#Customize the first dendrogram corresponding to the HiSST scheme UPGMA
dend12_corrected_2[[1]] = dend12_corrected_2[[1]] %>% 
  set("labels_col", k=6, value = c("orange","bisque4","sky blue","navy","dark green","dark red")) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color",k=6, value = c("orange","bisque4","sky blue","navy","dark green","dark red")) %>%
  set("leaves_pch", 19) %>%
  set("leaves_col", "dark cyan")

#Customize the second dendrogram corresponding to the cgMLST
dend12_corrected_2[[2]] = dend12_corrected_2[[2]] %>% 
  set("labels_col",k=6, value = c("orange","bisque4","sky blue","navy","dark green","dark red")) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color",k=6, value = c("orange","bisque4","sky blue","navy","dark green","dark red")) %>%
  set("leaves_pch", 19) %>%
  set("leaves_col", "dark cyan")


#### Plot trees together ####

graph = tanglegram(dend12_corrected_2[[1]],dend12_corrected_2[[2]], 
           margin_inner=16,
           common_subtrees_color_lines = FALSE, common_subtrees_color_branches = FALSE,
           highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE,
           lwd=1,
           columns_width = c(5, 2, 5),
           margin_outer = 0.5,
           main_left = "                                     HiSST SNPs",  main_right = "cgMLST                               ")



# ________Correlation tests________ ####

# Cophenetic correlation matrix 
dl <- dendlist(d1,d2)
#cor.dendlist(dl, method = "cophenetic")
cor_cophenetic(dl)
#Baker’s Gamma Index
cor_bakers_gamma(dl) #The value can range between -1 to 1. With near 0 values meaning that the two trees are not statistically similar

