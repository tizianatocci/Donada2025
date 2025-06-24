# Transcriptomic analysis
This repository contains the codes used for the article "Clonal memory of cell division diverges between human healthy haematopoiesis and acute myeloid leukaemia".\
Article's authors: Alessandro Donada; Gurvan Hermange; Tiziana Tocci; Adil Midoun; Giulio Prevedello; Louisa Hadj Abed; Delia Dupré; Wenjie Sun; Idan Milo; Sabrina Tenreira Bento; Costandina Pospori; Andrew Innes; Christophe Willekens; Jacques Vargaftig; David Michonneau; Cristina Lo Celso; Nicolas Servant; Ken R Duffy; Hervé Isambert; Paul-Henri Cournede; Lucie Laplane; Leila Perié

The folder "Code" contains the code used for the analysis.\
The folder "Datasets" contains the original data.\
The folder "Output" contains the files obtained during the analysis, such as the MIIC input.\

The analysis is divided in four steps, for each of them a code is available to replicate the results of the article:
1. Azimuth annotation 
2. Cell fate annotation and obtain figures for the paper concerning annotation, gene expression 
3. Cosine distance computation and obtain figures for the paper concerning the cosine distance per family 
4. Run feature selection for MIIC 

After point 4, the MIIC network is obtained with the MIIC server: https://miic.curie.fr/ \
Each code returns intermediate files, stored in the "Output" folder. The files are already present in the "Output" folder in case you want to replicate only one part of the analysis.

R Version used: 4.3.2\
Seurat Version: 4\
Dependencies:\
install.packages("Seurat")\
install.packages("tidyverse")\
install.packages("matrixStats")\
install.packages("miic")
