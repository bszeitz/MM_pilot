##############################
#
# Therapy response prediction with Cox regression analysis
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
Annotation <- read.delim("Annotation_table_SampleClusters.txt", check.names = F)
Protein.Info <- read.delim("Protein_info_table_SampleClusterANOVA.txt", check.names = F)
row.names(Expr.Table) <- Expr.Table$Accession
Expr.Table$Accession <- NULL
row.names(Annotation) <- Annotation$UniqueID
Annotation <- Annotation[colnames(Expr.Table),]


###########################
# Create time, status, and Met annotations, 
# and separate annotation tables for immuno- and targeted therapy
###########################


Annotation$time <- Annotation$Progression.time.months
Annotation$status <- Annotation$Progression.event
Annotation$Met <- as.factor(ifelse(Annotation$`Type of Samples`=="Met", 1, 0))

Annotations.Therapy <- list(Immun= Annotation[Annotation$Progression.therapyType =="Immun" & !is.na(Annotation$Progression.therapyType),],
                            Target= Annotation[Annotation$Progression.therapyType =="Target" & !is.na(Annotation$Progression.therapyType),])


###########################
# Perform Cox regression analysis
###########################

VVfilter <- 80
#Cox.Results.MetStrata <- list(Immun= Cox_withStrata(Expr.Table[,row.names(Annotations.Therapy$Immun)], 
#                                                    VVfilter, Annotations.Therapy$Immun, "ImmunoTh"),
#                              Target= Cox_withStrata(Expr.Table[,row.names(Annotations.Therapy$Target)], 
#                                                     VVfilter, Annotations.Therapy$Target, "TargetedTh"))
#save(Cox.Results.MetStrata, file="Cox.results.RData")
load("Cox.results.RData")


###########################
# Annotate proteins based on their suggested predictivity
###########################

Cox.Results.MetStrata$Immun$Prediction_ImmunoTh <- ""
Cox.Results.MetStrata$Immun$Prediction_ImmunoTh <- ifelse(Cox.Results.MetStrata$Immun$Cox_p.v_ImmunoTh < 0.05 & 
                                                  Cox.Results.MetStrata$Immun$Cox_HR_ImmunoTh < 1 & 
                                                  Cox.Results.MetStrata$Immun$Proport.Haz_p.v_ImmunoTh > 0.05, "Upregulated in patients with better response",
                                                Cox.Results.MetStrata$Immun$Prediction_ImmunoTh)
Cox.Results.MetStrata$Immun$Prediction_ImmunoTh <- ifelse(Cox.Results.MetStrata$Immun$Cox_p.v_ImmunoTh < 0.05 & 
                                                  Cox.Results.MetStrata$Immun$Cox_HR_ImmunoTh > 1 & 
                                                  Cox.Results.MetStrata$Immun$Proport.Haz_p.v_ImmunoTh > 0.05, "Downregulated in patients with better response",
                                                Cox.Results.MetStrata$Immun$Prediction_ImmunoTh)

Cox.Results.MetStrata$Target$Prediction_TargetedTh <- ""
Cox.Results.MetStrata$Target$Prediction_TargetedTh <- ifelse(Cox.Results.MetStrata$Target$Cox_p.v_TargetedTh < 0.05 & 
                                                   Cox.Results.MetStrata$Target$Cox_HR_TargetedTh < 1 & 
                                                   Cox.Results.MetStrata$Target$Proport.Haz_p.v_TargetedTh > 0.05, "Upregulated in patients with better response",
                                                 Cox.Results.MetStrata$Target$Prediction_TargetedTh)
Cox.Results.MetStrata$Target$Prediction_TargetedTh <- ifelse(Cox.Results.MetStrata$Target$Cox_p.v_TargetedTh < 0.05 & 
                                                   Cox.Results.MetStrata$Target$Cox_HR_TargetedTh > 1 & 
                                                   Cox.Results.MetStrata$Target$Proport.Haz_p.v_TargetedTh > 0.05, "Downregulated in patients with better response",
                                                 Cox.Results.MetStrata$Target$Prediction_TargetedTh)


###########################
# Merge results from immuno and targeted
###########################

Cox.Results.MetStrata.Both <- merge(Cox.Results.MetStrata$Immun, Cox.Results.MetStrata$Target, by="Accession",
                                    all=T)


immun.sign <- Cox.Results.MetStrata.Both[Cox.Results.MetStrata.Both$Prediction_ImmunoTh!="" & !is.na(Cox.Results.MetStrata.Both$Prediction_ImmunoTh),]
target.sign <- Cox.Results.MetStrata.Both[Cox.Results.MetStrata.Both$Prediction_TargetedTh!="" & !is.na(Cox.Results.MetStrata.Both$Prediction_TargetedTh),]


###########################
# Visualize results from the immunotherapy subgroup
###########################

Annotations.Therapy$Immun <- Annotations.Therapy$Immun[order(Annotations.Therapy$Immun$Progression.time.months),]

Annotations.Therapy$Immun$Progression.time.months.Cat <- factor(Annotations.Therapy$Immun$Progression.time.months.Cat,
                                                                levels=c("[0,3]","(3,7]","(7,12]","(12,60]"))

column_ha <- create_column_ha(Annotations.Therapy$Immun, c("Progression.time.months.Cat","Progression.event"), colorlist)

draw(Heatmap(t(scale(t(Expr.Table[immun.sign$Accession,row.names(Annotations.Therapy$Immun)]))), 
             name="Z-score",
             col=colorRamp2(c(-2, 0, 2), c("blue","white", "red")),
             cluster_rows = T, top_annotation = column_ha, 
             cluster_columns = F,
             row_split = immun.sign$Prediction_ImmunoTh,
             cluster_column_slices = F,
             clustering_distance_rows = "euclidean",
             #column_split = progression.time,
             #column_labels = columnnames,
             show_row_names = F, 
             show_heatmap_legend = T,
             width = unit(15, "cm"),
             #row_title = NULL,
             #row_title_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 10)),
     #row_split = rowsplits, gap = unit(5, "mm"),
     #height  = unit(20, "cm")),
     heatmap_legend_side="right",
     annotation_legend_side="bottom")

pdf(file="Immunotherapy_response_predicting_proteins.pdf",width = 15, height = 10)
draw(Heatmap(t(scale(t(Expr.Table[immun.sign$Accession,row.names(Annotations.Therapy$Immun)]))), 
             name="Z-score",
             col=colorRamp2(c(-2, 0, 2), c("blue","white", "red")),
             cluster_rows = T, top_annotation = column_ha, 
             cluster_columns = F,
             row_split = immun.sign$Prediction_ImmunoTh,
             cluster_column_slices = F,
             clustering_distance_rows = "euclidean",
             #column_split = progression.time,
             #column_labels = columnnames,
             show_row_names = F, 
             show_heatmap_legend = T,
             width = unit(15, "cm"),
             #row_title = NULL,
             #row_title_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 10)),
     #row_split = rowsplits, gap = unit(5, "mm"),
     #height  = unit(20, "cm")),
     heatmap_legend_side="right",
     annotation_legend_side="bottom")
dev.off()


###########################
# Visualize results from the targeted therapy subgroup
###########################

Annotations.Therapy$Target$Progression.time.months.Cat <- factor(Annotations.Therapy$Target$Progression.time.months.Cat,
                                                                levels=c("[0,3]","(3,7]","(7,12]","(12,60]"))

Annotations.Therapy$Target <- Annotations.Therapy$Target[order(Annotations.Therapy$Target$Progression.time.months),]

column_ha <- create_column_ha(Annotations.Therapy$Target, c("Progression.time.months.Cat","Progression.event"), colorlist)


draw(Heatmap(t(scale(t(Expr.Table[target.sign$Accession,row.names(Annotations.Therapy$Target)]))), 
             name="Z-score",
             col=colorRamp2(c(-2, 0, 2), c("blue","white", "red")),
             cluster_rows = T, top_annotation = column_ha, 
             cluster_columns = F,
             row_split = target.sign$Prediction_TargetedTh,
             cluster_column_slices = F,
             clustering_distance_rows = "euclidean",
             #column_split = progression.time,
             #column_labels = columnnames,
             show_row_names = F, 
             show_heatmap_legend = T,
             width = unit(15, "cm"),
             #row_title = NULL,
             #row_title_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 10)),
     #row_split = rowsplits, gap = unit(5, "mm"),
     #height  = unit(20, "cm")),
     heatmap_legend_side="right",
     annotation_legend_side="bottom")

pdf(file="Targetedtherapy_response_predicting_proteins.pdf",width = 15, height = 10)
draw(Heatmap(t(scale(t(Expr.Table[target.sign$Accession,row.names(Annotations.Therapy$Target)]))), 
             name="Z-score",
             col=colorRamp2(c(-2, 0, 2), c("blue","white", "red")),
             cluster_rows = T, top_annotation = column_ha, 
             cluster_columns = F,
             row_split = target.sign$Prediction_TargetedTh,
             cluster_column_slices = F,
             clustering_distance_rows = "euclidean",
             #column_split = progression.time,
             #column_labels = columnnames,
             show_row_names = F, 
             show_heatmap_legend = T,
             width = unit(15, "cm"),
             #row_title = NULL,
             #row_title_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 10)),
     #row_split = rowsplits, gap = unit(5, "mm"),
     #height  = unit(20, "cm")),
     heatmap_legend_side="right",
     annotation_legend_side="bottom")
dev.off()




###########################
# Export results
###########################

Protein.Info.ANOVA.Cox <- merge(Protein.Info, Cox.Results.MetStrata.Both, by="Accession", all.x=T)

#write.table(Protein.Info.ANOVA.Cox, "Protein_info_table_SampleClusterANOVA_TherapyCox.txt", 
#            sep="\t", quote = F, row.names = F, na="")



