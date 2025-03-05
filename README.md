# Transcriptomic analysis
This repository contains the codes used for the article "Clonal memory of cell division in human healthy hematopoiesis and Acute Myeloid Leukemia".\
Article's authors: Donada A, Hermange G, Tocci T, Midoun A, Prevedello G, Hadj Abed L, Dupré D, Sun W, Milo I, Tenreira Bento S, Pospori C, Willekens C, Vargaftig J, Michonneau D, Lo Celso C, Servant N, Duffy K, Isambert H, Cournède PH, Laplane L and Perié L

The folder "Code" contains the code used for the analysis.\
The folder "Datasets" contains the original data, which must be downloaded from *insert link*.\
The folder "Output" contains the files obtained during the analysis, such as the MIIC input.\
The folder "MIIC_summary" contains the file downloaded from the MIIC network to obtain the list of genes significantly associated with properties of interest.\
The folder "List genes" contains the list of genes reported in the article, significantly associated with properties of interest.\
The folder "References" contains the file needed to apply Azimuth, which must be downloaded from *insert link*.\
The folder "Figures" contains the figures that the code will generate. These are the figures used in the original article.

The analysis is divided in five steps, for each of them a code is available to replicate the results of the article:
1. Azimuth annotation 
2. Cell fate annotation and obtain figures for the paper concerning annotation, gene expression 
3. Cosine distance computation and obtain figures for the paper concerning the cosine distance per family 
4. Run feature selection for MIIC 
5. Obtain list of genes from MIIC summary 

Between point 4 and 5, the MIIC network is obtained with the MIIC server: https://miic.curie.fr/ \
Each code returns intermediate files, stored in the "Output" folder. The files are already present in the "Output" folder in case you want to replicate only one part of the analysis.

R Version used: 4.3.2\
Seurat Version: 4\
Dependencies:\
install.packages("Azimuth")\
install.packages("Seurat")\
install.packages("tidyverse")\
install.packages("matrixStats")\
install.packages("miic")
… Work in progress 
