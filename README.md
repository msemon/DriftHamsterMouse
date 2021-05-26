# DriftHamsterMouse

This repository contains code and basic data to reproduce the analyses presented in the manuscript "XX".

# Cusp patterning analysis 

`cuspid_Erlang_2021_GamBoxCoxModel.r`

This script first fits a relationship between embryo weight and age in days post coitum (dpc) by GAM+boxcox or linear+boxcox tranformation. This is using the data of hundreds of embryos in `weightTime_hamstermouse.txt`. Then age is predicted for embryos used for cusp patterning models (`cusp_fgf4_hamster.xlsx` and `cusp_fgf4_mouse.xlsx`) and RNAseq samples (`metadata-samplesRNAseq.txt`). 

In a second time, the models for cusp patterning are fitted with continuous time Marlov chains.
CountTot.txt

# Deconvolutions 

`DeconvoEPIMES_and_BL.md`

This script performs deconvolutions from markers extracted from mesenchyme/epithelium samples and bucco/lingual samples.  This is using the file annotating the samples `metadataTot_FINAL.txt`, the age predictions for plotting `predictions_rnaseq_().txt` made in previous section, and count data `CountTot.txt`. 
Control plots to check the accuracy of prediction in pure tissues and boostraps made by resampling the marker lists are made in the code.


# Spline DE analyses and GO enrichment

`ComputeByPairs_GAM.R`

This code fits splines to temporal expression profiles, and then tests for differential expression with DESeq2. The proportion of differentially expressed genes in several lists of genes (`liste_bite-it.csv`, `keystone_genes.csv`, `dispensable_genes.csv`, `pathways-Margaux-corMS.xlsx`) is measured and plotted, and the general gene ontology enrichment is computed in :

`FigUpLow.R` (for upper and lower comparisons)
`FigHamSou.R` (for mouse and hamster comparisons)


# ROMA pathway activation analyses

Pathway activation is measured by rROMA package, for three pathways with their target genes(`BMP-avril2019.csv`, `WNT-avril2019.csv`, `genesetSHH.txt`). 

# Multivariate analyses



# Crop genes sequences to orthologous parts in hamster and mouse 

Code in `CropGenome.tar.xz`
