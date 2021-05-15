#'Customize the Tanglegram of HiSST scheme and whole genome sequences.


#'This script was inspired from URL https://www.r-graph-gallery.com/340-custom-your-dendrogram-with-dendextend
#'and has been adapted for the development of the HiSST scheme, which helps to select the loci of the scheme.
#'
#'The following script draws an UPGMA dendrogram based on the ANIb score of concatenated loci 
#'selected for HiSST scheme and genome similarity to discriminate strains of Serratia marcescens;
#'using pvclust package for Hierchical Clustering.
#'Dendextend and tidyverse packages were used for the customized Tanglegram.
#'
#'It uses the files called "ANIb_S_marcescens_gabR-bssA-dhaM.tab" and "ANIb_WG_S_marcescens.tab",
#'available at https://github.com/TBourd/Rscript_for_HiSST_scheme_of_S_marcescens

library(pvclust)

setwd("~/Tanglegram") #Set the working directory path of your project



#### HiMLST ####
Sm_UPGMA_HiSST=read.table("ANIb_S_marcescens_gabR-bssA-dhaM.tab") # Import ANIb scores of concatenated loci
Sm_UPGMA_HiSST=pvclust(Sm_UPGMA_HiSST, method.dist="cor", method.hclust="average", nboot=1000,store=TRUE) # Change nboot value for a faster run, ex.: 'nboot=100'


#### WG ####

Sm_UPGMA_WG=read.table("ANIb_WG_S_marcescens.tab") # Import ANIb scores based on whole genomes sequences
Sm_UPGMA_WG=pvclust(Sm_UPGMA_WG, method.dist="cor", method.hclust="average", nboot=1000,store=TRUE) # Change nboot value for a faster run, ex.: 'nboot=100'


#### Compilation ####

library(dendextend)

library(tidyverse)


d1 = as.dendrogram(Sm_UPGMA_HiSST)

labels(d1)=gsub("_"," ",labels(d1)) #Delete '_' character in the sample names


d2 = as.dendrogram(Sm_UPGMA_WG)

labels(d2)=gsub("_"," ",labels(d2))


dend12_corrected_1 <- untangle_random_search(d1,d2) #random shuffling the two untangled dendrogram and each time checking if their entanglement was improved
dend12_corrected_2 <- untangle_step_rotate_2side(dend12_corrected_1[[1]],dend12_corrected_1[[2]]) #rotate the tree and looking for a better match. Use it if you think the two trees could better match.


#Customize the first dendrogram corresponding to the HiSST scheme UPGMA
dend12_corrected_2[[1]] = dend12_corrected_2[[1]] %>% 
  set("labels_col",value = c("grey","orange","dark green","navy","sky blue"), k=5) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color",value = c("grey","orange","dark green","navy","sky blue"), k=5) %>%
  set("leaves_pch", 19) %>%
  set("leaves_col", "dark cyan")

#Customize the second dendrogram corresponding to the WG UPGMA
dend12_corrected_2[[2]] = dend12_corrected_2[[2]] %>% 
  set("labels_col",value = c("grey","orange","dark green","navy","sky blue"), k=5) %>%
  set("branches_lty", 1) %>%
  set("branches_k_color",value = c("grey","orange","dark green","navy","sky blue"), k=5) %>%
  set("leaves_pch", 19) %>%
  set("leaves_col", "dark cyan")


# Plot them together

graph = tanglegram(dend12_corrected_2[[1]],dend12_corrected_2[[2]], 
           margin_inner=16,
           common_subtrees_color_lines = FALSE, common_subtrees_color_branches = FALSE,
           highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE,
           lwd=1,
           columns_width = c(5, 2, 5),
           margin_outer = 0.5,
           main_left = "                                     HiSST scheme",  main_right = "WG                               ")


