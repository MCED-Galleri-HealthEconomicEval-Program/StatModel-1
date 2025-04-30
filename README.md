# BayesianModelTestSens-1
Files to support "A Bayesian approach to sharing information on sensitivity of a Multi-Cancer Early Detection test across and within tumour types and stages".

This repository contains data, JAGS and R code, and some additional results supporting work done as part of a comprehensive program of NHS England funded research 
assessing the cost-effectiveness of the Multi-Cancer Early Detection (MCED) Galleri test within a national screening framework to support informed policy decisions.

The files included here support the publication ***"A Bayesian approach to sharing information on sensitivity of a Multi-Cancer Early Detection test across and within tumour types and stages"***
which sough to model currently available test sensitivity data.

Files are included to implement different modelling assumptions, as follows:

1. Text files starting with _'SensModel...'_ contain JAGS code to implement the models used (6 files in total). Basic description of each model is given within each file.
2. _'CCGAdatav1.txt'_ contains the MCED sensitivity data used, to be loaded into R.
3. Two R scripts are provided, one to generate the JAGS datasets for each model (_Galleri JAGS datasets.R_) and one to run the JAGS models (_Galleri binomial models v1.R_).
4. Text files starting with _'data...'_ contain data for the JAGS models.
5. Text files including '_...Groups..._' contain data for the class models.

Supplementary results displaying residual deviances for each datapoint in the models are included in _'HeatMaps Pub v1.xlsx'_ and _'HeatMaps Pub v1.pdf'_ as Microsof Excel and PDF files respectively 
(the information in both files is the same).
