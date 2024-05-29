# Chronic stress deficits in reward behaviour co-occur with low nucleus accumbens dopamine activity during reward anticipation specifically

This repo contains the data, code and R Markdown files that accompany the paper by Zhang C., Dulinskas R. *et al.*. These scripts were used to process and analyse the transcriptomics dataset produced for the study.

This work was conducted by the [Preclinical Lab for Translational Research into Affective Disorders (PLATRAD)](https://www.dppp.uzh.ch/en/researchgroups/researchgroups/preclinicallab.html) led by Prof. Christopher R. Pryce in collaboration with scientists from the University of Zurich, the University Hospital of Psychiatry Zurich, ETH Zurich, the Peking University School of Life Sciences and Boehringer Ingelheim Pharma GmbH & Co. KG.

The repo includes the following files:

- `data/1191_DA_neurons_VEH_TPM.xlsx`: Matrix of TPM-normalised gene expression values across samples.
- `data/SummExp_1191_CSS_CON_VEH_VTA_DA_neurons_Pryce.rds`: R data storage file (rds) with a `SummarizedExperiment` object containing all sample metadata, gene metadata, raw and normalised gene expression matrices.
- `data_analysis_DA_neurons_VEH.Rmd`: R Markdown document with the steps followed to analyse the transcriptomics dataset produced for this study.
- `data_analysis_DA_neurons_VEH.html`: The result of knitting the above R Markdown document to produce the full bioinformatics report for this study.
- `volcano_plot_DA_neurons_VEH.R`: R script to generate a volcano plot for the main differential expression contrast in this study.

## R version

Please note that the above analyses were performed with `R version 4.1.2`.