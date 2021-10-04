##############################
#
# Utility functions
#
# Written by Beáta Szeitz
#
##############################


##############################
# Perform median normalization (centering around the global median)
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows

# This function returns the a new protein expression table where 
# samples are median-normalized (centered around the global median)
#
# Packages required: none
#
normalize_median <- function (df) {
  df2 <- df
  med <- median(as.matrix(df), na.rm=T)
  for (i in 1:ncol(df)){
    df2[,i] <- (df[,i] - median(df[,i], na.rm = T) + med)
  }
  return(df2)
}


##############################
# Prepare table in long format for ggplot
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows
#
# This function returns the protein expression table in a long
# format with columns Sample, Protein and Intensity.
#
# Packages required:
# - reshape2
# - varhandle
#
prepare_forggplot <- function (df) {
  df <- t(as.matrix(df))
  df <- cbind(row.names(df),df)
  df <- reshape2::melt(df, id.vars="V1")
  df <- varhandle::unfactor(df)
  df$value <- suppressWarnings(as.numeric(as.character(df$value)))
  df <- df[df[,2] !="",]
  colnames(df) <- c("Sample", "Protein","Intensity")
  return(df)
}


##############################
# Visualize protein intensity distributions across samples
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows
# - name: the title of the plot
#
# This function returns a series of boxplots to visualize the
# protein intensity distributions in each sample.
#
# Packages required:
# - ggplot2
# 
create_sample_histogram <- function(df, name){
  sampleorder <- colnames(df)
  df <- prepare_forggplot(df)
  df.expression.gg_noNA <- na.omit(df)
  df.expression.gg_noNA$Sample <- factor(df.expression.gg_noNA$Sample, levels = sampleorder)
  median.value <- median(df.expression.gg_noNA$Intensity)
  g_cell_raw <- ggplot(df.expression.gg_noNA, aes(x=Sample, y=Intensity))
  g_cell_raw + geom_violin() + geom_hline(yintercept = median.value)+ 
    theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(text= element_text(size=12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.0),
          legend.position = "none")+
    geom_boxplot(width=0.1) +
    stat_summary(fun.y=median, geom="point", size=2, color="red") +
    labs(title=name, 
         y="Protein Intensities", x="Samples")
  
}


##############################
# Perform ANOVA and pairwise Tukey HSD test
# 
# Arguments for this function:
# - m: a protein expression table where samples are in columns and
#   proteins are in rows
# - annotation: the annotation table of the samples
# - colname_for_factor: the column name in the annotation table
#   that should be used as the independent variable
# - levelOrder: how should the levels of the independent variable
#   be ordered
# - analysisName: if the columns in the result table should be 
#   appended by a certain string
# - filter: what is the minimum ratio of valid values in each
#   level of the independent variable to be considered for ANOVA
# - stars: whether the result table should contain columns that
#   indicate significance with stars (*** p<0.005, ** p<0.01,
#   * p<0.05, . p<0.10)
#
# This function returns a result table that contains:
# - Protein accession numbers
# - ANOVA test F, p-value, Benjamini-Hochberg adjusted p-value
# - Pairwise Tukey HSD test p-values
# - Pairwise Log2FC values
# - 95% confidence intervals for the pairwise log2FC values
#
# Packages required:
# - varhandle
# - DTK
# - reshape2
# - car
#
ANOVAandTK_padj_noMVfilter <- function(m, annotation, colname_for_factor, levelOrder, analysisName="", filter, stars) {
  annotation = as.matrix(cbind(row.names(annotation), annotation[,colname_for_factor]))
  colnames(annotation) <- c("Sample", colname_for_factor)
  annotation <- as.data.frame(annotation)
  if (is.factor(annotation[,1])){
    annotation[,1] <- varhandle::unfactor(annotation[,1])
  } 
  
  if (is.factor(annotation[,2])){
    annotation[,2] <- factor(varhandle::unfactor(annotation[,2]), levels=levelOrder)
  } else {
    annotation[,2] <- factor(annotation[,2], levels=levelOrder)
  }
  
  annotation <- na.omit(annotation)
  m <- m[,annotation$Sample]
  
  min.sample.sizes <- vector(length=length(levels(annotation[,2])))
  names(min.sample.sizes) <- levels(annotation[,2])
  for (i in 1:length(min.sample.sizes)){
    min.sample.sizes[i] <- nrow(annotation[annotation[,2]==names(min.sample.sizes)[i],]) * filter
  }
  
  p.v <- matrix(nrow=nrow(m), ncol=1)
  Fvalues <- matrix(nrow=nrow(m), ncol=1)
  n <- as.numeric(length(levels(annotation[,2])))
  k <- 2 
  comb <- factorial(n)/(factorial(k)*factorial(n-k))
  
  
  # Get pairwise comparison names
  m.noNA <- na.omit(m)
  prot.data <- as.data.frame(t(m.noNA[1,]))
  data <- cbind(prot.data,annotation)
  data$Sample <- NULL
  tukey <- DTK::TK.test(x=data[,1],f=data[,2], a= 0.05)
  tukey <- reshape2::melt(tukey[["f"]])
  tukey.pv0 <- subset(tukey, Var2 =="p adj")
  comparisons <- tukey.pv0[,1]
  
  
  tukey.pv <- matrix(nrow=nrow(m), ncol=comb)
  colnames(tukey.pv) <- paste("tukey(",comparisons,")",sep="")
  tukey.diff <- matrix(nrow=nrow(m), ncol=comb)
  colnames(tukey.diff) <- paste("Log2FC(",comparisons,")",sep="")
  tukey.confint <- matrix(nrow=nrow(m), ncol=comb)
  colnames(tukey.confint) <- paste("Log2FC.95%.CI(",comparisons,")",sep="")
  
  for (i in 1:nrow(m)) {
    prot.data <- as.data.frame(t(m[i,]))
    prot.data$Sample <- row.names(prot.data)
    data <- merge(prot.data,annotation,by="Sample")
    row.names(data) <- data$Sample
    data <- data[,which(colnames(data) !="Sample")]
    data <- na.omit(data)
    groups.to.include <- vector()
    
    for (N in 1:length(min.sample.sizes)){
      groups.to.include <- c(groups.to.include, ifelse(nrow(data[data[,2]==names(min.sample.sizes)[N],]) >= min.sample.sizes[N],
                                                       names(min.sample.sizes)[N], NA))
    }
    groups.to.include <- na.omit(groups.to.include)
    
    if (length(groups.to.include)==0 | length(groups.to.include)==1){ next}
    data <- data[which(data[,2] %in% groups.to.include),]
    
    lm_res <- lm(data[,1] ~ data[,2],data=data, na.action=na.exclude)
    ANOVAres <- car::Anova(lm_res, type=2)
    Fvalues[i,1] <- paste0(round(ANOVAres$`F value`[1], 4), " (df=",ANOVAres$Df[1],")" )
    p.v[i,1] <- ANOVAres$`Pr(>F)`[1]
    tukey <- DTK::TK.test(x=data[,1],f=data[,2], a= 0.05)
    tukey <- reshape2::melt(tukey[["f"]])
    
    # tukey.PV
    tukey.pv0 <- subset(tukey, Var2 =="p adj")
    tukey.pv1 <- as.data.frame(t(tukey.pv0[,3]))
    colnames(tukey.pv1) <- paste("tukey(",tukey.pv0[,1],")",sep="")
    for (K in 1:ncol(tukey.pv1)){
      tukey.pv[i,colnames(tukey.pv1)[K]] <- tukey.pv1[1,K]
    }
    
    # tukey.diff
    tukey.diff0 <- subset(tukey, Var2 =="diff")
    tukey.diff1 <- as.data.frame(t(tukey.diff0[,3]))
    colnames(tukey.diff1) <- paste("Log2FC(",tukey.diff0[,1],")",sep="")
    for (K in 1:ncol(tukey.diff1)){
      tukey.diff[i,colnames(tukey.diff1)[K]] <- tukey.diff1[1,K]
    }
    
    # tukey.CI
    tukey.CI0 <- subset(tukey, Var2 =="lwr" | Var2 =="upr")
    tukey.CI0$Var1 <- varhandle::unfactor(tukey.CI0$Var1)
    comps <- unique(tukey.CI0$Var1)
    tukey.CI <- tukey.CI0[1:length(comps),]
    
    for (K in 1:length(comps)){
      uprCI <- round(tukey.CI0[tukey.CI0$Var1 == comps[K] & tukey.CI0$Var2 =="upr","value"],4)
      lwrCI <- round(tukey.CI0[tukey.CI0$Var1 == comps[K] & tukey.CI0$Var2 =="lwr","value"],4)
      tukey.CI[K,3] <- paste0("[",lwrCI," , ", uprCI,"]")
    }
    tukey.CI2 <- as.data.frame(t(tukey.CI[,3]))
    colnames(tukey.CI2) <- paste("Log2FC.95%.CI(",tukey.CI[,1],")",sep="")
    for (K in 1:ncol(tukey.CI2)){
      tukey.confint[i,colnames(tukey.CI2)[K]] <- tukey.CI2[1,K]
    }
  }
  row.names(p.v) <- row.names(m)
  p.v <- as.data.frame(p.v)
  colnames(p.v) <- c(paste("p.v", colnames(annotation)[2]))
  row.names(tukey.pv) <- row.names(m)
  row.names(tukey.diff) <- row.names(m)
  adjp <- p.adjust(as.vector(p.v[,1]), "fdr")
  adjp <- as.data.frame(as.numeric(as.character(adjp)))
  p.v_ANOVA <- as.data.frame(cbind(p.v,adjp,tukey.pv))
  colnames(p.v_ANOVA)[2] <- c(paste("p.adj", colnames(annotation)[2])) 
  if (stars){
    p.v_star <- p.v_ANOVA
    for (C in 1:ncol(p.v_star)){
      p.v_star[,C] <- cut(as.numeric(p.v_ANOVA[,C]), breaks=c(-Inf,0.005, 0.01, 0.05, 0.1, Inf), 
                          label=c("***", "**", "*", ".",""))
    }
    colnames(p.v_star) <- paste("*",colnames(p.v_star) , sep="_")
    table_results <- cbind(Fvalues, p.v_ANOVA,tukey.diff,tukey.confint, p.v_star)
  } else {
    table_results <- cbind(Fvalues, p.v_ANOVA,tukey.diff,tukey.confint)
  }
  
  if (analysisName !=""){
    colnames(table_results) <- paste(colnames(table_results) ,analysisName, sep="_")
  }
  
  table_results <- as.data.frame(cbind(row.names(table_results),table_results))
  colnames(table_results)[1] <- "Accession"
  return(table_results)
}

##############################
# Filter for valid values based on percentage
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows
# - percentage: the minimum percentage of valid values required for 
#   each protein
#
# This function returns a protein expression table that is filtered
# for valid values.
#
# Packages required: none
# 
filter_missingvalues <- function (df, percentage) {
  perc <- round(ncol(df)*(percentage/100))
  vv <- rowSums(!is.na(df))
  to.be.filtered <- unlist(lapply(vv, function(x){ 
    if (x < perc) { y <- F}
    else y <- T}))
  return(as.data.frame(df[to.be.filtered,]))
}

##############################
# Extract Fisher's exact test results
# 
# Arguments for this function:
# - d: the contingency table
#
# This function returns a vector of Fisher test results:
# - confidence interval, odds ratio, p-value
#
# Packages required: none
# 
perform_fisher <- function(d){
  fis <- fisher.test(d,conf.int = TRUE, conf.level = 0.95, alternative = "less")
  CI <- paste0(round(fis$conf.int[1],4)," - ",round(fis$conf.int[2],4))
  OR <- round(fis[["estimate"]][["odds ratio"]],4)
  p <- round(fis$p.value,4)
  return(c(CI, OR, p))
}

##############################
# Summarize Fisher's exact test results for all levels of a categorical
# variable for the elements of one sample cluster
# 
# Arguments for this function:
# - Annot.Table: the sample annotation table
# - samplecluster.column: the column name that indicates the sample clusters
# - attribute: the clinical or histopathological characteristic to be tested
# - cluster: the name of the sample cluster to be tested
#
# This function returns a table of Fisher test results:
# - the sample ratios: the ratio of samples in category within the sample 
#   cluster and the ratio outside of the sample cluster
# - Fisher test confidence interval, odds ratio, p-value
#
# Packages required: none
# 
enrichment_test <- function(Annot.Table,samplecluster.column,attribute,cluster){
  
  Annot.Table <- subset(Annot.Table, !is.na(Annot.Table[,attribute]))
  Annot.Table[,attribute] <- as.factor(Annot.Table[,attribute])
  attr.levels <- levels(Annot.Table[,attribute])
  Annot.Table.in.cluster <- subset(Annot.Table, Annot.Table[,samplecluster.column] ==cluster)
  Annot.Table.not.in.cluster <- subset(Annot.Table, Annot.Table[,samplecluster.column] !=cluster)
  
  ResultTable <- matrix(nrow=0, ncol=4)
  colnames(ResultTable) <- c("PatientRatio","ConfInt","OddsRatio","P.value")
  
  for (i in 1:length(attr.levels)){
    
    in.cluster <- nrow(Annot.Table.in.cluster[Annot.Table.in.cluster[,attribute] ==attr.levels[i],])
    out.cluster <- nrow(Annot.Table.not.in.cluster[Annot.Table.not.in.cluster[,attribute] ==attr.levels[i],])
    
    d <- data.frame(patients.not.interest=c(out.cluster,
                                            nrow(Annot.Table.not.in.cluster) - out.cluster), 
                    patients.in.interest=c(in.cluster, nrow(Annot.Table.in.cluster)- in.cluster), 
                    row.names = c("in.cluster","not.in.cluster"))
    
    ResultTable.sub <- data.frame(PatientRatio = paste0("InCluster: ",in.cluster,"/",nrow(Annot.Table.in.cluster),
                                                        ", OutCluster: ",out.cluster,"/", nrow(Annot.Table.not.in.cluster)),
                                  ConfInt = perform_fisher(d)[1],
                                  OddsRatio = perform_fisher(d)[2],
                                  p.value = perform_fisher(d)[3], row.names = paste0(attribute, " - ",attr.levels[i]))
    
    ResultTable <- rbind(ResultTable, ResultTable.sub)
    
  }
  return(ResultTable)
}

##############################
# Short function to create the column annotations for heatmaps
# 
# Arguments for this function:
# - Annotation: the sample annotation table
# - columnnames: the annotation table's columns to be included
# - colorlist: the list of colors
#
# This function returns the column heatmap annotation.
#
# Packages required:
# - ComplexHeatmap
# 
create_column_ha <- function(Annotation, columnnames, colorlist){
  column_ha = HeatmapAnnotation(df =Annotation[,columnnames],
                                which="col",
                                col=colorlist,
                                annotation_name_side = "left",
                                gp = gpar(col = "grey"),
                                show_legend = TRUE,
                                show_annotation_name = TRUE)
}

##############################
# Perform multiple Cox regression analysis with strata on primary/metastasis
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows
# - VVfilter: the minimum number of valid values required for each protein
#   to be considered for the analysis
# - annot: the annotation table of the samples, containing the columns
#   "time", "status","Met"
# - analysisName: if the columns in the result table should be 
#   appended by a certain string
#
# This function returns a result table that contains:
# - Protein accession numbers
# - Hazard ratios, p-values, 95% confidence intervals for HRs
# - p-value of testing for the proportional hazards assumption
#
# Packages required:
# - survival
# - varhandle
#
Cox_withStrata <- function(df, VVfilter, annot,analysisName){
  df.filt <- apply(df, 1 ,function(x) {
    if (sum(!is.na(x)) >= (ncol(df)*VVfilter/100) ) {
      y <- "YES"
    } else { y <- "NO"} 
  }
  )
  df <- df[which(df.filt == "YES"),]
  df <- df[,row.names(annot)]
  results <- matrix(nrow=nrow(df), ncol=4)
  annot <- annot[,c("time","status","Met")]
  colnames(results) <- c("Cox_HR", "Cox_p.v", "Cox_HR.95%.CI",
                         "Proport.Haz_p.v")
  row.names(results) <- row.names(df)
  for (i in 1:nrow(df)) {
    Covariate <- as.data.frame(t(df[i,]))
    colnames(Covariate) <- "Protein"
    Covariate <- merge(Covariate, annot, by="row.names")
    row.names(Covariate) <- Covariate$Row.names
    Covariate$Row.names <- NULL
    time <- Covariate$time
    status <- Covariate$status
    Covariate$time <- NULL
    Covariate$status <- NULL
    multi_res <- coxph(Surv(time = time, event = status) ~ Protein+ strata(Met), data = Covariate)
    PH.pv <- cox.zph(multi_res)$table[1,3]
    summ <- summary(coxph(Surv(time,status) ~ Protein+ strata(Met), data=Covariate))
    #results[i,1] <- row.names(df)[i]
    melt.summ <- varhandle::unfactor(melt(summ[["coefficients"]]))
    results[i,"Cox_HR"] <- melt.summ[melt.summ$Var2 =="exp(coef)",3]
    results[i,"Cox_p.v"] <- melt.summ[grep("Pr(", melt.summ$Var2, fixed=T),3]
    results[i,"Proport.Haz_p.v"] <- PH.pv
    CIs <- varhandle::unfactor(melt(summ[["conf.int"]]))
    lwrCI <- round(CIs[CIs$Var2 == "lower .95","value"],4)
    uprCI <- round(CIs[CIs$Var2 == "upper .95","value"],4)
    results[i,"Cox_HR.95%.CI"] <- paste0("[",lwrCI," , ", uprCI,"]")
    
  }
  results.p.adj <- p.adjust(results[,"Cox_p.v"], method="fdr")
  results <- as.data.frame(cbind(results, results.p.adj))
  colnames(results)[ncol(results)] <- "Cox_p.adj"
  if (analysisName !=""){
    colnames(results) <- paste(colnames(results) ,analysisName, sep="_")
  }
  
  results <- as.data.frame(cbind(row.names(results),results))
  colnames(results)[1] <- "Accession"
  return(results)
}
