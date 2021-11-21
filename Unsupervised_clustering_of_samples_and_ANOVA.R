##############################
#
# Unsupervised clustering of the samples and ANOVA test
#
# Written by Beáta Szeitz
#
##############################


###########################
# Set wd
###########################

setwd("C:/Users/User/PhD/MM_pilot/Scripts")


###########################
# Load libraries and custom scripts
###########################
.packages = c("ggplot2","ggbiplot","dplyr","tibble", "proBatch", "gridExtra", 
              "cowplot", "ComplexHeatmap", "circlize", "WGCNA", "RColorBrewer",
              "ggpubr", "survival", "reshape2")
lapply(.packages, library, character.only=TRUE)
source("Utility_functions.R")
source("Color_list.R")


###########################
# Load tables
###########################

Expr.Table <- read.delim("Protein_intensity_table_batchCorr.txt", check.names = F)
Annotation <- read.delim("Annotation_table.txt", check.names = F)
Protein.Info <- read.delim("Protein_info_table.txt", check.names = F)
row.names(Expr.Table) <- Expr.Table$Accession
Expr.Table$Accession <- NULL
row.names(Annotation) <- Annotation$UniqueID
Annotation <- Annotation[colnames(Expr.Table),]



###########################
# Filter for missing values (<50%)
###########################

Expr.Table.F50 <- filter_missingvalues(Expr.Table, 50)


###########################
# Compare dynamic tree cut and constant height cut
###########################

distM <- dist(scale(t(Expr.Table.F50)),method = "euclidean")
hc <- hclust(distM, "complete")

set.seed(12345)
Sample.Clusters.C<- data.frame(Name = colnames(Expr.Table.F50),
                               Cluster = cutree(hc, k=12),
                               row.names = colnames(Expr.Table.F50))
Annotation$Constant.height.cut <- Sample.Clusters.C[Annotation$UniqueID,"Cluster"]

Annotation$Constant.height.cut <- factor(Annotation$Constant.height.cut, levels =
                                           c("1","2","3","4","5","6","7","8","9","10","11",
                                             "12"))

set.seed(12345)
Sample.Clusters.D<- data.frame(Name = colnames(Expr.Table.F50),
                             Cluster = cutreeDynamic(hc, minClusterSize=0, method="hybrid",
                                                     distM=as.matrix(distM, deepSplit=4, 
                                                                     maxCoreScatter=NULL, minGap=NULL, 
                                                                     maxAbsCoreScatter=NULL, minAbsGap=NULL)),
                             row.names = colnames(Expr.Table.F50))

Annotation$Dynamic.cut <- Sample.Clusters.D[Annotation$UniqueID,"Cluster"]

column.ha.cluster <- create_column_ha(Annotation, 
                                      c("Constant.height.cut","Dynamic.cut"), 
                                      colorlist)

draw(Heatmap(t(scale(t(Expr.Table.F50))), 
             name="Z-score",
             col=colorRamp2(c(-4, 0, 4), c("blue","white", "red")),
             cluster_rows = F, top_annotation = column.ha.cluster, 
             cluster_columns = T,
             show_row_names = F, 
             show_heatmap_legend =F, 
             show_column_names = T,
             row_title = NULL,
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 10),
             height  = unit(0.01, "cm")),
     annotation_legend_side="bottom")


pdf(file="Constant_vs_Dynamic_Cut.pdf",width = 15, height = 5)
draw(Heatmap(t(scale(t(Expr.Table.F50))), 
             name="Z-score",
             col=colorRamp2(c(-4, 0, 4), c("blue","white", "red")),
             cluster_rows = F, top_annotation = column.ha.cluster, 
             cluster_columns = T,
             show_row_names = F, 
             show_heatmap_legend =F, 
             show_column_names = T,
             row_title = NULL,
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 10),
             height  = unit(0.01, "cm")),
     annotation_legend_side="bottom")
dev.off()

###########################
# Unsupervised clustering of the samples with dynamic tree cut
###########################

set.seed(12345)
Sample.Clusters<- data.frame(Name = colnames(Expr.Table.F50),
                             Cluster = cutreeDynamic(hc, minClusterSize=0, method="hybrid",
                                                     distM=as.matrix(distM, deepSplit=4, 
                                                                     maxCoreScatter=NULL, minGap=NULL, 
                                                                     maxAbsCoreScatter=NULL, minAbsGap=NULL)),
                             row.names = colnames(Expr.Table.F50))

Annotation$SampleCluster <- Sample.Clusters[Annotation$UniqueID,"Cluster"]

summary(as.factor(Annotation$SampleCluster))



###########################
# Overrepresentation analysis of clinical and histopathological characteristics
###########################

clusters <- unique(Annotation$SampleCluster)
clusters <- clusters[order(clusters)]
clusters <- clusters[clusters!="0"]

cluster.enrichment.p0.2 <- list()

annot.columns <- c("Type of Samples", "Organ of the Samples", "Tumor content (%).Cat", "Therapy groups","BRAFstate",
                   "Age at sample collectin.Cat", "Sex", "DFS (m).Cat", "PFS (m).Cat", "OS (m).Cat",
                   "Live", "Loc_p", "AJCC8", "Type.Simplified", "T", "U", "Breslow (mm).Cat", "Regres", "BRAFstate_simple",
                   "MSOrder.Cat", "Year.of.collection.Cat", "No.MVs.Cat", "Stage","Progression.time.months.Cat")


Enrichment.Results.All <- matrix(nrow=0, ncol=6)
colnames(Enrichment.Results.All) <- c("PatientRatio","ConfInt","OddsRatio","P.value","Cluster","Trait")

for (i in 1:length(clusters)){
  
  Enrichment.Results.Sub <- matrix(nrow=0, ncol=4)
  colnames(Enrichment.Results.Sub) <- c("PatientRatio","ConfInt","OddsRatio","P.value")
  
  for (C in 1:length(annot.columns)){
    Enrichment.Results.Sub <- rbind(Enrichment.Results.Sub, 
                                    enrichment_test(Annotation, "SampleCluster",annot.columns[C], clusters[i], "fisher"))
  }
  
  colnames(Enrichment.Results.Sub)[4] <- "P.value"
  
  Enrichment.Results.Sub$Cluster <- clusters[i]
  Enrichment.Results.Sub$Trait <- row.names(Enrichment.Results.Sub)
  
  Enrichment.Results.All <- rbind(Enrichment.Results.All, Enrichment.Results.Sub)
  
  cluster.enrichment.p0.2[[i]] <- Enrichment.Results.Sub[as.numeric(Enrichment.Results.Sub$P.value) < 0.1,]
  names(cluster.enrichment.p0.2)[i] <- paste0("Cluster",clusters[i])
  
  
}

# Print results

for (i in 1:length(cluster.enrichment.p0.2)){
  if (nrow(cluster.enrichment.p0.2[[i]])>0){
    cluster.enrichment.p0.2[[i]]$P.value <- as.numeric(cluster.enrichment.p0.2[[i]]$P.value)
    printtable <- cluster.enrichment.p0.2[[i]][order(cluster.enrichment.p0.2[[i]]$P.value),]
    print(knitr::kable(printtable, 
                       caption = paste0(names(cluster.enrichment.p0.2)[i],
                                        ", Fisher test p < 0.1")))
  }
}


# Export results

Enrichment.Results.All.Export <- Enrichment.Results.All[,c("Cluster","Trait", "P.value",
                                                           "PatientRatio","OddsRatio")]

#write.table(Enrichment.Results.All.Export, "SampleClusters_Enrichment_Results.txt", 
#            sep="\t", quote = F, row.names = F)



###########################
# ANOVA with sample clusters
###########################

#ANOVA.Cluster <- ANOVAandTK_padj_noMVfilter(annotation = Annotation, colname_for_factor = "SampleCluster" ,
#                                            levelOrder = c("1", "2", "3", "4","5","6"), 
#                                            m = as.data.frame(Expr.Table), 
#                                            analysisName= "", filter=0.8, stars =F)

#save(ANOVA.Cluster, file="ANOVA.Cluster.RData")
load("ANOVA.Cluster.RData")



###########################
# Clustering of the top 1000 differentially expressed proteins
###########################


ANOVA.Cluster.Top1000 <- ANOVA.Cluster[order(ANOVA.Cluster$`p.adj SampleCluster`) ,]
ANOVA.Cluster.Top1000 <- ANOVA.Cluster.Top1000[1:1000,]

distM.ANOVA.Cluster <- dist(t(scale(t(Expr.Table[ANOVA.Cluster.Top1000$Accession,]))),method = "euclidean")
hc.ANOVA.Cluster <- hclust(distM.ANOVA.Cluster, "complete")

set.seed(12345)
ProteinClusters<- data.frame(Accession = ANOVA.Cluster.Top1000$Accession,
                             Top1000Protein.Cluster = cutreeDynamic(hc.ANOVA.Cluster, minClusterSize=50, method="hybrid",
                                                                    distM=as.matrix(distM.ANOVA.Cluster, deepSplit=0, 
                                                                                    maxCoreScatter=NULL, minGap=NULL, 
                                                                                    maxAbsCoreScatter=NULL, minAbsGap=NULL)))

summary(as.factor(ProteinClusters$Top1000Protein.Cluster))



###########################
# Visualize sample clusters with top1000 proteins
###########################


Expr.Table.sign <- Expr.Table[ProteinClusters$Accession,]

column_ha <- create_column_ha(Annotation, 
                               c("Type of Samples",
                                 "Organ of the Samples",
                                 "Therapy groups",
                                 "AJCC8",
                                 "T", "U", 
                                 "Breslow (mm).Cat",
                                 "Regres",
                                 "Loc_p", # Primary localization
                                 "BRAFstate",
                                 "Type.Simplified", # Melanoma type
                                 "Live",
                                 "DFS (m).Cat",
                                 "PFS (m).Cat",
                                 "OS (m).Cat",
                                 "Sex",
                                 "Age at sample collectin.Cat",
                                 "No.MVs.Cat", # No. of missing protein intensities
                                 "Year.of.collection.Cat",
                                 "Tumor content (%).Cat"
                               ), 
                               colorlist)


draw(Heatmap(t(scale(t(Expr.Table.sign))), 
             name="Z-score",
             col=colorRamp2(c(-2, 0, 2), c("blue","white", "red")),
             cluster_rows = T, top_annotation = column_ha, 
             cluster_columns = T,
             row_split = ProteinClusters$Top1000Protein.Cluster,
             #clustering_distance_rows = "pearson",
             column_split = Annotation$SampleCluster,
             #column_labels = columnnames,
             show_row_names = F, 
             #height  = unit(0.01, "cm"),
             width = unit(18, "cm"),
             show_heatmap_legend = T, 
             show_row_dend = F,
             show_column_names = F,
             #row_title = NULL,
             #row_title_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 10)),
     #row_split = rowsplits, gap = unit(5, "mm"),
     #height  = unit(0.01, "cm")),
     heatmap_legend_side ="left",
     annotation_legend_side="right")




pdf(file="SampleClusters_Top1000proteins.pdf", 
    width=20, 
    height=15, 
    pointsize=12)
draw(Heatmap(t(scale(t(Expr.Table.sign))), 
             name="Z-score",
             col=colorRamp2(c(-2, 0, 2), c("blue","white", "red")),
             cluster_rows = T, top_annotation = column_ha, 
             cluster_columns = T,
             row_split = ProteinClusters$Top1000Protein.Cluster,
             #clustering_distance_rows = "pearson",
             column_split = Annotation$SampleCluster,
             #column_labels = columnnames,
             show_row_names = F, 
             #height  = unit(0.01, "cm"),
             width = unit(18, "cm"),
             show_heatmap_legend = T, 
             show_row_dend = F,
             show_column_names = F,
             #row_title = NULL,
             #row_title_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 10)),
     #row_split = rowsplits, gap = unit(5, "mm"),
     #height  = unit(0.01, "cm")),
     heatmap_legend_side ="left",
     annotation_legend_side="right")
dev.off()



###########################
# Export results
###########################

#write.table(Annotation, "Annotation_table_SampleClusters.txt", 
#            sep="\t", quote = F, row.names = F)


Protein.Info.ANOVA <-  merge(Protein.Info, ANOVA.Cluster, by="Accession")
#write.table(Protein.Info.ANOVA, "Protein_info_table_SampleClusterANOVA.txt", 
#            sep="\t", quote = F, row.names = F)



