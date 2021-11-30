# MM_pilot scripts

This is the custom R code used for the proteomic data analysis in:

Szadai, L., Velasquez, E. et al. (2021). Deep Proteomic Analysis on Paraffine-Archived Melanoma with Prognostic/Predictive Biomarker read-out. Manuscript accepted for publication.

**This repository consists of 7 scripts used for analysis and plotting of proteomic data:**

| Script name  | Description  |
|---|---|
| Utility_functions.R  | This contains custom functions used in the scripts. |
|  Color_list.R | This contains the colors specified for heatmap annotations.  |
| Proteomic_data_normalization_and_batch_effect_correction.R  | This contains the steps of proteomic data normalization and batch effect correction.  |
| Unsupervised_clustering_of_samples_and_ANOVA.R | This contains the steps of the unsupervised clustering of samples, followed by differential expression analysis with ANOVA. |
|Therapy_response_prediction.R | This contains the steps of multiple Cox regression analysis to identify proteins that are potential predictors of therapy response.  |
|Primary_metastasis_comparison.R | This contains the primary vs metastasis comparisons and subsequent pre-ranked GSEA.  |
|Stage_comparison.R | This contains the clinical stage comparisons and subsequent pre-ranked GSEA.  |

