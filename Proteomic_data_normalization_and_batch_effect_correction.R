##############################
#
# Proteomic data normalization and batch effect correction
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

Expr.Table <- read.delim("Protein_intensity_table.txt", check.names = F)
Annotation <- read.delim("Annotation_table.txt", check.names = F)
row.names(Expr.Table) <- Expr.Table$Accession
Expr.Table$Accession <- NULL
row.names(Annotation) <- Annotation$UniqueID
Annotation <- Annotation[colnames(Expr.Table),]

###########################
# Log2 transformation
###########################

Expr.Table.log2 <- apply(Expr.Table, c(1,2), log2)

###########################
# Median normalization
###########################

Expr.Table.log2.norm <- normalize_median(Expr.Table.log2)

Expr.Table.log2.numbered <- Expr.Table.log2
colnames(Expr.Table.log2.numbered) <- seq(1, ncol(Expr.Table.log2.numbered), 1)
histo.beforenorm <- create_sample_histogram(Expr.Table.log2.numbered, "Before normalization")

Expr.Table.log2.norm.numbered <- Expr.Table.log2.norm
colnames(Expr.Table.log2.norm.numbered) <- seq(1, ncol(Expr.Table.log2.norm.numbered), 1)
histo.afternorm <- create_sample_histogram(Expr.Table.log2.norm.numbered, "After median-normalization")


###########################
# Spike-in protein intensity vs injection order
###########################


Lys <- as.data.frame(t(rbind(Expr.Table.log2["P00698",], Expr.Table.log2.norm["P00698",])))
colnames(Lys) <- c("Before Norm", "After Norm")
Lys <- Lys[row.names(Lys) !="Accession",]
Lys <- as.data.frame(apply(Lys,c(1,2), as.numeric))
Lys$Sample <- factor(row.names(Lys), levels=row.names(Lys), labels=seq(1,nrow(Lys),1))
Lys <- reshape2::melt(Lys)
Lys$Order <- varhandle::unfactor(Lys$Sample)


plot_Lys <- ggplot(Lys, aes(x=Order, y=value, color=variable, group=variable))+geom_point()+geom_smooth(method = "loess", se = F)+#geom_hline(aes(yintercept = 33))+
  ggtitle("Lysozyme C intensity vs injection order")+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color=guide_legend(title="Color"))+ scale_color_manual(values=c("darkgrey", "black"))+ ylab("Lysozyme C Protein Intensity") + xlab("Injection order")+ scale_x_continuous(breaks = seq(0, 90, by = 25))
plot_Lys



###########################
# PCA plot after median normalization
###########################


Annotation$Organ <- ifelse(Annotation$`Organ of the Samples`=="Cut", "C", "L")

plot_PCA_before <- ggbiplot(prcomp(t(na.omit(Expr.Table.log2.norm))),
                            circle=F, scale = T, 
                            labels=Annotation$Organ, 
                            groups = as.factor(Annotation$MSOrder.Cat),
                            ellipse = T,
                            var.axes	=F) + 
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("Median-normalized before batch correction") + guides(color=guide_legend(title="Injection order"))+ scale_color_manual(values=c("red", "#d4af37", "blue", "black"))+
  xlim(-4,6) + ylim(-2.5,2.5)
plot_PCA_before



###########################
# Continuous drift correction
###########################

# Prepare data for proBatch
Annotation$MSbatch <- ifelse(is.na(Annotation$MSbatch), 90, Annotation$MSbatch)
Annotation$MSbatch <- cut(Annotation$MSOrder, breaks=c(0,25,50,75,90), labels = c("1", "2", "3", "4"))

sample_annotation <- Annotation[,c("UniqueID", "MSbatch", "MSOrder", "Organ of the Samples")]
colnames(sample_annotation) <- c("FullRunName", "MS_batch", "order", "Organ")
peptide_annotation <- data.frame(peptide_group_label = row.names(Expr.Table.log2.norm),
                                 Gene =NA,
                                 ProteinName = row.names(Expr.Table.log2.norm))
color_list <- sample_annotation_to_colors(sample_annotation,
                                          factor_columns = c('MS_batch', 'Organ'),
                                          numeric_columns = c('order'))

median_normalized_long <- matrix_to_long(Expr.Table.log2.norm)
plot_spike_in(median_normalized_long, sample_annotation,
              peptide_annotation = peptide_annotation,
              protein_col = 'ProteinName', spike_ins = 'P00698',
              plot_title = 'Lys C',
              color_by_batch = TRUE, color_scheme = color_list[["MS_batch"]])


# Perform continuous drift correction
#loess_fit_90 <- adjust_batch_trend_df(median_normalized_long, sample_annotation,
#                                      span = 0.9)
#save(loess_fit_90, file="loess_fit_90.RData")
load("loess_fit_90.RData")


plot_loess <- plot_with_fitting_curve(feature_name = 'P00698',
                                      fit_df = loess_fit_90, fit_value_col = 'fit',
                                      df_long = median_normalized_long,
                                      sample_annotation = sample_annotation,
                                      color_by_batch = FALSE, color_scheme = color_list[[batch_col]],
                                      plot_title = 'LOESS fit, Span = 90%')+ theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("Lysozyme C Protein Intensity")+xlab("Injection order")
plot_loess


median_normalized_long_loess <- center_feature_batch_medians_df(loess_fit_90, sample_annotation)
plot_single_feature(feature_name = 'P00698', df_long = median_normalized_long_loess,
                    sample_annotation = sample_annotation, measure_col = 'Intensity',
                    plot_title = 'Feature-level Median Centered')

Epxr_batchcorr <- long_to_matrix(median_normalized_long_loess)

plot_sample_mean(Epxr_batchcorr, sample_annotation, order_col = 'order',
                 batch_col = "MS_batch", color_by_batch = TRUE, #ylimits = c(12, 16.5),
                 color_scheme = color_list[["MS_batch"]])



###########################
# iRT intensity vs injection order
###########################


iRT <- as.data.frame(t(rbind(Expr.Table.log2["Biognosys|iRT-Kit_WR_fusion",], 
                             Expr.Table.log2.norm["Biognosys|iRT-Kit_WR_fusion",],
                             Epxr_batchcorr["Biognosys|iRT-Kit_WR_fusion",])))
colnames(iRT) <- c("Before Norm", "After Norm","After Batchcorr")
iRT <- iRT[row.names(iRT) !="Accession",]
iRT <- as.data.frame(apply(iRT,c(1,2), as.numeric))

iRT$Sample <- factor(row.names(iRT), levels=row.names(iRT), labels=seq(1,nrow(iRT),1))
iRT <- reshape2::melt(iRT)
iRT$Order <- varhandle::unfactor(iRT$Sample)

plot_iRT <- ggplot(iRT, aes(x=Order, y=value, color=variable, group=variable))+geom_point()+geom_smooth(method = "loess", se = F)+#geom_hline(aes(yintercept = 33))+
  ggtitle("iRT intensity vs injection order")+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color=guide_legend(title="Color"))+ scale_color_manual(values=c("darkgrey", "black", "red"))+ ylab("iRT Protein Intensity")+ xlab("Injection order")+ scale_x_continuous(breaks = seq(0, 90, by = 25))
plot_iRT




###########################
# PCA plot after batch correction
###########################


Epxr_batchcorr.numbered <- Epxr_batchcorr
colnames(Epxr_batchcorr.numbered) <- seq(1, ncol(Epxr_batchcorr.numbered), 1)
histo.afterbatch <- create_sample_histogram(Epxr_batchcorr.numbered, "After median-normalization and batch correction")


plot_PCA_after <- ggbiplot(prcomp(t(na.omit(Epxr_batchcorr))),
                           circle=F, scale = T, 
                           labels=Annotation$Organ, 
                           groups = as.factor(Annotation$MSOrder.Cat),
                           ellipse = T,
                           var.axes	=F) + 
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("Median-normalized and batch corrected") + guides(color=guide_legend(title="Injection order"))+ scale_color_manual(values=c("red", "#d4af37", "blue", "black"))+
  xlim(-4,6) + ylim(-2.5,2.5)
plot_PCA_after




###########################
# Save suppl figure
###########################


top_row <- plot_grid(plot_Lys, plot_loess,
                     ncol = 2,
                     labels = c("A", "B"),
                     label_fontfamily = '',
                     label_fontface = 'bold',
                     label_size = 16,
                     align = 'h',
                     rel_widths = c(1.25,1))

middle_row <- plot_grid(plot_iRT,
                        ncol = 1,
                        labels = c("C"),
                        label_fontfamily = '',
                        label_fontface = 'bold',
                        label_size = 16,
                        align = 'h',
                        rel_widths = c(0.5))

bottom_row1 <- plot_grid(plot_PCA_before, plot_PCA_after,
                         ncol = 2,
                         labels = c("D"),
                         label_fontfamily = '',
                         label_fontface = 'bold',
                         label_size = 16,
                         axis="rlbt",
                         align = 'v',
                         rel_widths = c(1))


bottom_row2 <- plot_grid(histo.beforenorm,histo.afternorm, histo.afterbatch,
                         ncol = 1,
                         labels = c("E"),
                         label_fontfamily = '',
                         label_fontface = 'bold',
                         label_size = 16,
                         axis="rlbt",
                         align = 'v',
                         rel_widths = c(0.5,0.5))

plot_grid(top_row, middle_row, bottom_row1, bottom_row2, ncol = 1,
          rel_heights = c(0.2,0.2,0.5,0.7))


#ggsave(file="Suppl_Figure1.svg", width =12, height = 20)
#ggsave(file="Suppl_Figure1.jpg", width=12, height = 20)
#ggsave(file="Suppl_Figure1.pdf", width=12, height = 20)



###########################
# Save batch corrected expression table
###########################

Expr.Table.export <- Epxr_batchcorr
Expr.Table.export <- cbind(row.names(Expr.Table.export), Expr.Table.export)
colnames(Expr.Table.export)[1] <- "Accession" 

#write.table(Expr.Table.export, "Protein_intensity_table_batchCorr.txt", 
#            sep="\t", quote = F, row.names = F)

