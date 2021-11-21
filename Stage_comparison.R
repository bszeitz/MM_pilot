##############################
#
# Stage II - IV comparisons
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
Protein.Info <- read.delim("Protein_info_table_SampleClusterANOVA_TherapyCox_PrimVsMet.txt", check.names = F)
row.names(Expr.Table) <- Expr.Table$Accession
Expr.Table$Accession <- NULL
row.names(Annotation) <- Annotation$UniqueID
Annotation <- Annotation[colnames(Expr.Table),]
Protein.Info$GeneName <- unlist(lapply(Protein.Info$`Gene Symbol`, function(x){strsplit(x, split=";")[[1]][1]}))
row.names(Protein.Info) <- Protein.Info$Accession

###########################
# Extract stage II-IV patients
###########################

Annotation <- Annotation[Annotation$Stage !=1 & !is.na(Annotation$Stage),]
Expr.Table <- Expr.Table[,row.names(Annotation)]

###########################
# Perform ANCOVA (including stage and sample type in the model)
###########################


Annotation$`Type of Samples` <- factor(Annotation$`Type of Samples`, levels = c("Prim","Met"))
Annotation$Sample <- row.names(Annotation)

run = F
if (run){
        min.sample.sizes <- vector(length=3)
        names(min.sample.sizes) <- c("2","3","4")
        for (i in 1:length(min.sample.sizes)){
                min.sample.sizes[i] <- nrow(Annotation[Annotation$Stage==names(min.sample.sizes)[i],]) * 0.8
        }
        
        Stage.ANCOVA <- matrix(nrow=nrow(Expr.Table), ncol=4)
        colnames(Stage.ANCOVA) <- c("coeff stage","p.v stage","coeff type","p.v type")
        
        for (i in 1:nrow(Expr.Table)) {
                prot.data <- as.data.frame(t(Expr.Table[i,]))
                prot.data$'Sample' <- row.names(prot.data)
                data <- merge(prot.data,Annotation[,c("Sample", "Stage", "Type of Samples")],by="Sample")
                row.names(data) <- data$Sample
                data <- data[,which(colnames(data) !="Sample")] # remove the "sample" column
                data <- na.omit(data)
                groups.to.include <- vector()
                
                for (N in 1:length(min.sample.sizes)){
                        groups.to.include <- c(groups.to.include, 
                                               ifelse(nrow(data[data[,2]==names(min.sample.sizes)[N],]) >= min.sample.sizes[N], #this  "=" is new!!!
                                                      names(min.sample.sizes)[N], NA))
                }
                groups.to.include <- na.omit(groups.to.include)
                
                if (length(groups.to.include)!=3){ next}
                
                lm_res <- lm(data[,1] ~ data[,2]+data[,3],data=data, na.action=na.exclude)
                lm_res <- summary(lm_res)
                Stage.ANCOVA[i,1] <- lm_res[["coefficients"]][2]
                Stage.ANCOVA[i,3] <- lm_res[["coefficients"]][3]
                Stage.ANCOVA[i,2] <- lm_res[["coefficients"]][11]
                Stage.ANCOVA[i,4] <- lm_res[["coefficients"]][12]
        }
        row.names(Stage.ANCOVA) <- row.names(Expr.Table)
        Stage.ANCOVA <- as.data.frame(Stage.ANCOVA)
        Stage.ANCOVA$Accession <- row.names(Stage.ANCOVA)
        save(Stage.ANCOVA, file="Stage_ANCOVA.RData")
} else {
        load("Stage_ANCOVA.RData")    
}


Stage.ANCOVA$`p.adj stage` <- p.adjust(Stage.ANCOVA$`p.v stage`, method="fdr")
Stage.ANCOVA$`p.adj type` <- p.adjust(Stage.ANCOVA$`p.v type`, method="fdr")

nrow(Stage.ANCOVA[!is.na(Stage.ANCOVA$`p.v stage`),])
nrow(Stage.ANCOVA[!is.na(Stage.ANCOVA$`p.v stage`) & 
                                 Stage.ANCOVA$`p.v stage` < 0.05 &
                                 Stage.ANCOVA$`coeff stage` > 0,])
nrow(Stage.ANCOVA[!is.na(Stage.ANCOVA$`p.v stage`) & 
                          Stage.ANCOVA$`p.v stage` < 0.05 &
                          Stage.ANCOVA$`coeff stage` < 0,])


###########################
# Create ranking of proteins
###########################

ANCOVA.rank <- na.omit(data.frame(Accession = row.names(Stage.ANCOVA),
                                  Rank = -log10(Stage.ANCOVA$`p.v stage`) * Stage.ANCOVA$`coeff stage`,
                                  Gene = Protein.Info[row.names(Stage.ANCOVA),"GeneName"]))
ANCOVA.rank <- ANCOVA.rank[order(ANCOVA.rank$Rank, decreasing = T),]
ANCOVA.ranking <- ANCOVA.rank$Rank
names(ANCOVA.ranking) <- ANCOVA.rank$Gene


###########################
# Load gmt file, select GOBP, KEGG, Reactome and Wikipathways genesets (msigdb version 7.4),
# and perform pre-ranked Gene Set Enrichment Analysis
###########################

gmtfile <- read.gmt("msigdb.v7.4.symbols.gmt")
gmtfile.sub <- subset(gmtfile, grepl("KEGG_", gmtfile$term) | 
                              grepl("REACTOME_", gmtfile$term) |
                              grepl("GOBP_", gmtfile$term) |
                              grepl("WP_", gmtfile$term),)

GSEA.results.stage <- GSEA(ANCOVA.ranking, TERM2GENE=gmtfile.sub, verbose=FALSE, pvalueCutoff = 1.0, eps = 0)

GSEA.results.stage.matrix <- GSEA.results.stage@result

###########################
# Visualize pGSEA results
###########################


important.genesets <- c("GOBP_HUMORAL_IMMUNE_RESPONSE", "GOBP_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS", 
                        "GOBP_REGULATION_OF_HUMORAL_IMMUNE_RESPONSE", "GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION",
                        "GOBP_INTERSTRAND_CROSS_LINK_REPAIR", "GOBP_EPITHELIAL_CELL_DIFFERENTIATION",
                        "GOBP_NEGATIVE_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY")

de <- GSEA.results.stage.matrix[,c("ID","NES", "pvalue")]
colnames(de) <- c("GeneSet","NES", "pvalue")
de$sign <- ifelse(de$pvalue < 0.05, "yes", "no")
de$name <- ifelse(de$GeneSet %in% c(important.genesets), de$GeneSet, "")
de$color <- "NS"
de$color <- ifelse(de$sign =="yes" & de$NES > 0, "activated", de$color)
de$color <- ifelse(de$sign =="yes" & de$NES < 0, "suppressed", de$color)


de <- de[de$color !="NS",]

ggplot(data=de, aes(x=NES, y=-log10(pvalue), col=color, label=name)) + 
        geom_point() + #ylim(2,4.5)+
        theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        geom_label_repel(fill = "white", max.overlaps = 10,size = 3)+
        scale_color_manual(values=c("red", "blue"))+xlab("NES (Stage IV vs III vs II)")

#ggsave("Stage_GSEA.svg", width = 8, height = 4)


###########################
# Export results
###########################

Protein.Info.ANOVA.Cox.Type.Stage <- merge(Protein.Info, Stage.ANCOVA, by="Accession", all.x=T)

#write.table(Protein.Info.ANOVA.Cox.Type.Stage, "Protein_info_table_SampleClusterANOVA_TherapyCox_PrimVsMet_Stage.txt", 
#            sep="\t", quote = F, row.names = F, na="")
#write.table(GSEA.results.stage.matrix, "GSEA_results_Stage.txt", 
#            sep="\t", quote = F, row.names = F)





