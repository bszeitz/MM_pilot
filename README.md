# MM_pilot scripts

This is the custom R code used for the proteomic data analysis in:

Szadai, L.; Velasquez, E.; Szeitz, B.; Almeida, N.P.d.; Domont, G.; Betancourt, L.H.; Gil, J.; Marko-Varga, M.; Oskolas, H.; Jánosi, Á.J.; Boyano-Adánez, M.d.C.; Kemény, L.; Baldetorp, B.; Malm, J.; Horvatovich, P.; Szász, A.M.; Németh, I.B.; Marko-Varga, G. Deep Proteomic Analysis on Biobanked Paraffine-Archived Melanoma with Prognostic/Predictive Biomarker Read-Out. Cancers 2021, 13, 6105. https://doi.org/10.3390/cancers13236105

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

