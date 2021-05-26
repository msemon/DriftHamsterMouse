# DriftHamsterMouse

This repository contains code and basic data to reproduce the analyses presented in the manuscript "XX".

# Cusp patterning analysis 

`cuspid_Erlang_2021_GamBoxCoxModel.r`

This script first fits a relationship between embryo weight and age in days post coitum (dpc) by GAM+boxcox or linear+boxcox tranformation. This is using the data of hundreds of embryos in `weightTime_hamstermouse.txt`. Then age is predicted for embryos used for cusp patterning models (`cusp_fgf4_hamster.xlsx` and `cusp_fgf4_mouse.xlsx`) and RNAseq samples (`metadata-samplesRNAseq.txt`). 

In a second time, the models for cusp patterning are fitted with continuous time Marlov chains.

# 


# Crop genes sequences to orthologous parts in hamster and mouse 

Code in `CropGenome.tar.xz`
