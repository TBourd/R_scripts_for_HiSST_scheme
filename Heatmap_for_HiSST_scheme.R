#'Create a heatmap to visualize ANIb analysis of HiSST scheme or Whole genome


#'This script follows the instructions from URL https://jokergoo.github.io/ComplexHeatmap-reference/book/,
#'using ComplexHeatmap package.
#'
#'The following script customizes a heatmap adapted for ANIb analysis visualization of HiSST scheme or Whole genome
#'using S. marcescens strains (n= 60) and Serratia spp. strains (n=9) 
#'without the strain S_proteamaculans_336X, whose ANIb score is too dissimilar from S. marcescens strains.
#'
#'ggplot2, circlize and dendextend were used for the customized heatmap.
#'
#'It uses the files called "HiSST_ANIb_Serratia_sp.tab" and "WG_ANIb_Serratia_sp.tab",
#'available at https://github.com/TBourd/R_scripts_for_HiSST_scheme/tree/main/Data



library(ggplot2)

#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

library(circlize)
library(dendextend)

setwd("~/HiSST_Heatmap") #Set the working directory path of your project

#import data

ANIb = read.table("HiSST_ANIb_Serratia_sp.tab") # Import ANIb scores of concatenated loci

#Remove S_proteamaculans_336X : not any similarity with S. marcescens strains
ANIb = ANIb[!(row.names(ANIb) %in% c("S_proteamaculans_336X")),!(names(ANIb) %in% c("S_proteamaculans_336X"))]
str(ANIb2)


#Replace row and column names
names(ANIb) <- gsub(x = names(ANIb), pattern = "_", replacement = " ")
names(ANIb) <- gsub(x = names(ANIb), pattern = "S marcescens ", replacement = "")


row.names(ANIb) = gsub(x = row.names(ANIb), pattern = "_", replacement = " ")
row.names(ANIb) = gsub(x = row.names(ANIb), pattern = "S marcescens ", replacement = "")

#Order dataframe by names for standardize and compare heatmaps
ANIb.order = ANIb[order(rownames(ANIb)),order(names(ANIb))]

#For Heatmap, data must be numerical matrix
ANIb.mat = as.matrix(ANIb.order)

#Change clustering threshold by color
col_fun = colorRamp2(c(0.86, 0.94, 1), c("blue", "white", "red"))

#Cluster 
dend = hclust(dist(ANIb.mat))
dend = color_branches(dend, k = 2, col = c("red", "blue"))

#Heatmap
Heatmap(ANIb.mat, name = "mat",
        col = col_fun, #heatmap color
        heatmap_legend_param = list( #Legend parameters
          at = c(0, 0.25, 0.5, 0.75,0.90, 1),
          title = "ANIb percentage identity",
          legend_height = unit(4, "cm"),
          title_position = "lefttop-rot"),
        row_names_gp = gpar(fontsize = 5), #row names size
        column_names_gp = gpar(fontsize = 5), #column name size
        heatmap_width = unit(15, "cm"), heatmap_height = unit(15, "cm"), #heatmap size
        row_split = 2,
        column_split = 2,
        cluster_rows = dend,
        cluster_columns = dend,
        row_title_gp = gpar(col = c("red", "blue")),
        row_title = c(expression(italic(S.~marcescens)),expression(italic(Serratia)*' sp.')),
        column_title_gp = gpar(col = c("red", "blue")),
        column_title = c(expression(italic(S.~marcescens)),expression(italic(Serratia)*' sp.'))
        #cluster_rows = FALSE # turn off row clustering
)

#Save as '.svg' file : width = 700 , Height = 600  

