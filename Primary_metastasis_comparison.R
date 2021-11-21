##############################
#
# Primary vs metastasis comparisons
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
              "ggpubr", "survival", "reshape2", "lmerTest", "clusterProfiler",
              "ggrepel")
lapply(.packages, library, character.only=TRUE)
source("Utility_functions_vs2.R")
source("Color_list.R")


###########################
# Load tables
###########################

Expr.Table <- read.delim("Protein_intensity_table_batchCorr.txt", check.names = F)
Annotation <- read.delim("Annotation_table_SampleClusters.txt", check.names = F)
Protein.Info <- read.delim("Protein_info_table_SampleClusterANOVA_TherapyCox.txt", check.names = F)
row.names(Expr.Table) <- Expr.Table$Accession
Expr.Table$Accession <- NULL
row.names(Annotation) <- Annotation$UniqueID
Annotation <- Annotation[colnames(Expr.Table),]
Protein.Info$GeneName <- unlist(lapply(Protein.Info$`Gene Symbol`, function(x){strsplit(x, split=";")[[1]][1]}))
row.names(Protein.Info) <- Protein.Info$Accession

###########################
# Perform paired t-test on the paired samples
###########################

duplicated.patient <- c("Pilot-1","Pilot-2","Pilot-22","Pilot-25","Pilot-27","Pilot-28","Pilot-3","Pilot-52")
Annotation.Paired <- Annotation[Annotation$ID %in% duplicated.patient,]
Expr.Table.Paired <- as.matrix(Expr.Table[,row.names(Annotation.Paired)])

#Paired.test.Type <- paired_Ttest(m = Expr.Table.Paired, 
#                                 annotation = Annotation.Paired, 
#                                 patientID.column = "ID",
#                                 colname.for.factor = "Type of Samples",
#                                 analysisName = "Paired",
#                                 levelorder = c("Prim","Met"))
#save(Paired.test.Type, file="Paired_test_Type.RData")
load("Paired_test_Type.RData")

###########################
# Perform t-test on all samples
###########################

#Independent.test.Type <- ANOVAandTK_padj_noMVfilter(m = as.data.frame(Expr.Table), 
#                                                    annotation = Annotation, 
#                                                    colname_for_factor = "Type of Samples", 
#                                                    levelOrder = c("Prim", "Met"), 
#                                                    analysisName="Full", 
#                                                    filter=0.8, 
#                                                    stars = F)
#save(Independent.test.Type, file="Independent_test_Type.RData")
load("Independent_test_Type.RData")


###########################
# Create ranking of proteins
###########################

Paired.rank <- na.omit(data.frame(Accession = row.names(Paired.test.Type),
                                  Rank = -log10(Paired.test.Type$`p.v  Type of Samples_Paired`) * 
                                          Paired.test.Type$`Log2FC(Met-Prim)_Paired`,
                                  Gene = Protein.Info[row.names(Paired.test.Type),"GeneName"]))
Paired.rank <- Paired.rank[order(Paired.rank$Rank, decreasing = T),]
Paired.ranking <- Paired.rank$Rank
names(Paired.ranking) <- Paired.rank$Gene

Full.rank <- na.omit(data.frame(Accession = row.names(Independent.test.Type),
                                Rank = -log10(Independent.test.Type$`p.v Type of Samples_Full`) * 
                                        Independent.test.Type$`Log2FC(Met-Prim)_Full`,
                                Gene = Protein.Info[row.names(Independent.test.Type),"GeneName"]))
Full.rank <- Full.rank[order(Full.rank$Rank, decreasing = T),]
Full.ranking <- Full.rank$Rank
names(Full.ranking) <- Full.rank$Gene


###########################
# Load gmt file, select only KEGG genesets (msigdb version 7.4),
# and perform pre-ranked Gene Set Enrichment Analysis
###########################

gmtfile <- read.gmt("msigdb.v7.4.symbols.gmt")
gmtfile.sub <- subset(gmtfile, grepl("KEGG_", gmtfile$term) | 
                              grepl("REACTOME_", gmtfile$term) |
                              grepl("GOBP_", gmtfile$term) |
                              grepl("WP_", gmtfile$term),)
GSEA.results.paired <- GSEA(Paired.ranking, TERM2GENE=gmtfile.sub, verbose=FALSE, pvalueCutoff = 1.0, eps = 0)
GSEA.results.full <- GSEA(Full.ranking, TERM2GENE=gmtfile.sub, verbose=FALSE, pvalueCutoff = 1.0, eps = 0)


###########################
# Merge pGSEA results in one table
###########################

GSEA.results.sub.P <- GSEA.results.paired@result[,-c(2,12)]
colnames(GSEA.results.sub.P)[-1] <- paste(colnames(GSEA.results.sub.P)[-1], "Paired", sep="_")
GSEA.results.sub.F <- GSEA.results.full@result[,-c(2,12)]
colnames(GSEA.results.sub.F)[-1] <- paste(colnames(GSEA.results.sub.F)[-1], "Full", sep="_")
GSEA.results.export <- merge(GSEA.results.sub.P, GSEA.results.sub.F, by="ID", all.x=T, all.y = T)

for (j in 1:nrow(GSEA.results.export)){
        GSEA.results.export$core_enrichment_both[j] <- paste(intersect(unlist(lapply(GSEA.results.export$core_enrichment_Paired[j], 
                                                                                       function(x){strsplit(x, split="/", fixed = T)[[1]]})),
                                                                         unlist(lapply(GSEA.results.export$core_enrichment_Full[j], 
                                                                                       function(x){strsplit(x, split="/", fixed = T)[[1]]}))), 
                                                               collapse = "/")
}

GSEA.results.export <- GSEA.results.export[!is.na(GSEA.results.export$ID),]


GSEA.results.export <- GSEA.results.export[,c("ID", "NES_Paired", "NES_Full", "p.adjust_Paired", "p.adjust_Full",
                                                  "core_enrichment_Paired", "core_enrichment_Full","core_enrichment_both",
                                                  "pvalue_Paired", "pvalue_Full", "setSize_Paired", "setSize_Full", "enrichmentScore_Paired", "enrichmentScore_Full",
                                                  "qvalues_Paired", "qvalues_Full", "rank_Paired", "rank_Full", "leading_edge_Paired", "leading_edge_Full")]

GSEA.results.export$Sign.p_Paired <- ifelse(GSEA.results.export$pvalue_Paired < 0.05 & !is.na(GSEA.results.export$pvalue_Paired), "*","")
GSEA.results.export$Sign.p_Full <- ifelse(GSEA.results.export$pvalue_Full < 0.05 & !is.na(GSEA.results.export$pvalue_Full), "*","")
GSEA.results.export$Sign.p_both <- ifelse(GSEA.results.export$Sign.p_Paired =="*" & GSEA.results.export$Sign.p_Full =="*", "*","")
GSEA.results.export$Sign.padj_Paired <- ifelse(GSEA.results.export$p.adjust_Paired < 0.10 & !is.na(GSEA.results.export$p.adjust_Paired), "*","")
GSEA.results.export$Sign.padj_Full <- ifelse(GSEA.results.export$p.adjust_Full < 0.10 & !is.na(GSEA.results.export$p.adjust_Full), "*","")
GSEA.results.export$Sign.padj_both <- ifelse(GSEA.results.export$Sign.padj_Paired =="*" & GSEA.results.export$Sign.padj_Full =="*", "*","")
GSEA.results.export$Summarized.p <- GSEA.results.export$pvalue_Paired + GSEA.results.export$pvalue_Full

GSEA.results.export <- subset(GSEA.results.export, !is.na(GSEA.results.export$ID))


###########################
# Visualize pGSEA results
###########################

GSEA.results.export$Significance <- ifelse(GSEA.results.export$Sign.padj_both =="*", "Both", "NS")
GSEA.results.export$Significance <- ifelse(GSEA.results.export$Sign.padj_both !="*" & 
                                                   GSEA.results.export$Sign.padj_Paired =="*", "Only.Paired", GSEA.results.export$Significance)
GSEA.results.export$Significance <- ifelse(GSEA.results.export$Sign.padj_both !="*" & 
                                                   GSEA.results.export$Sign.padj_Full =="*", "Only.Full", GSEA.results.export$Significance)

colornames <- c("purple", "orange", "blue","grey")
names(colornames) <- c("Both", "Only.Paired","Only.Full","NS")

GSEA.results.export.fdr <- GSEA.results.export
GSEA.results.export.fdr <- GSEA.results.export.fdr[!is.na(GSEA.results.export.fdr$ID),]
GSEA.results.export.fdr$NES_Paired <- ifelse(is.na(GSEA.results.export.fdr$NES_Paired), 0, GSEA.results.export.fdr$NES_Paired)

#GSEA.results.export.fdr$GenesSet <- lapply(GSEA.results.export.fdr$ID, function(x){strsplit(x, split="KEGG_")[[1]][2]})
GSEA.results.export.fdr$GenesSet <- ifelse(GSEA.results.export$Significance =="NS", "", GSEA.results.export.fdr$ID)

#GSEA.results.export.fdr

GSEA.results.export.fdr <- GSEA.results.export.fdr[GSEA.results.export.fdr$p.adjust_Full < 0.05 | 
                                                           GSEA.results.export.fdr$p.adjust_Paired < 0.05,]
GSEA.results.export.fdr <- GSEA.results.export.fdr[!is.na(GSEA.results.export.fdr$ID),]

GSEA.results.export.fdr$NES_Full <- ifelse(is.na(GSEA.results.export.fdr$NES_Full),0,GSEA.results.export.fdr$NES_Full)

ggplot(GSEA.results.export.fdr, aes(x=NES_Paired, y=NES_Full, label=GenesSet, color = Significance)) + 
        geom_point() +geom_label_repel(fill = "white", max.overlaps = 1000,size = 3)+#geom_text_repel(size = 2,max.overlaps = 1000) + 
        ggtitle("Prim vs Met comparison - on Paired vs on Full dataset")+
        xlim(c(-2.5,2.5)) + ylim (-2.5,2.5)+
        geom_hline(yintercept=0, linetype="dashed", color = "darkgrey")+geom_vline(xintercept=0, linetype="dashed", color = "darkgrey")+
        theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
        xlab("NES - Met vs Prim (Paired dataset)")+ ylab("NES - Met vs Prim (Full dataset)")+scale_color_manual(values=colornames) 

#ggsave("PrimvsMet_GSEA_KEGG.svg", width = 8, height = 5)


###########################
# Export results
###########################

Protein.Info.ANOVA.Cox.Type <- merge(Protein.Info, Paired.test.Type, by="Accession", all.x=T)
Protein.Info.ANOVA.Cox.Type <- merge(Protein.Info.ANOVA.Cox.Type, Independent.test.Type[,c(1,3,4,6)], by="Accession", all.x=T)

#write.table(Protein.Info.ANOVA.Cox.Type, "Protein_info_table_SampleClusterANOVA_TherapyCox_PrimVsMet.txt", 
#            sep="\t", quote = F, row.names = F, na="")
#write.table(GSEA.results.export, "GSEA_results_PrimVsMet.txt", 
#            sep="\t", quote = F, row.names = F, na="")





