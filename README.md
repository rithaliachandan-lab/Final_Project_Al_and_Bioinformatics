Title : Expression dynamics of NKX2-1, EGFR and TP53 define molecular heterogeneity in lung adenocarcinoma.

Overview
This project performs a comprehensive RNA-seq–based biomarker analysis of TCGA Lung Adenocarcinoma (LUAD) data. The expression of key genes NKX2-1, EGFR, and TP53 is studied to compare Primary Tumor vs Solid Tissue Normal samples and to investigate tumor heterogeneity using correlation, heatmap, and PCA analyses. The processed data are further prepared for interactive Shiny visualization.

Objectives
- Extract TPM expression values from TCGA LUAD RNA-seq data
- Compare NKX2-1 expression between Tumor and Normal samples
- Perform multi-gene analysis for NKX2-1, EGFR, TP53
- Conduct gene–gene correlation analysis
- Visualize tumor heterogeneity using Heatmap and PCA
- Prepare the final dataset for Shiny web dashboard

Dataset used :
At this URL, you will find the following datasets:  https://figshare.com/s/fd7276e3583b457bd61d  
1. Sample sheet file (gdc_sample_sheet.tsv) 
2. Manifest file (gdc_manifest.txt) 
3. Gene expression data (tcga_data.tar.gz)


Outputs
- Tumor vs Normal boxplots
- Multi-gene expression comparison
- Gene–gene correlation scatter plots
- Heatmap of tumor biomarker expression
- PCA for tumor clustering
- Final processed dataset for Shiny: ss_rna_for_shiny.rds

Tools
- R Programming Language
- Shiny (for optional interactive dashboard)


Conclusion
This project establishes a complete bioinformatics workflow for TCGA LUAD biomarker analysis, demonstrating the expression patterns, correlations, and heterogeneity of NKX2-1, EGFR, and TP53 and enabling both static and interactive visualization.# Final_Project_Al_and_Bioinformatics
