# Transcriptomic analysis
This repository contains the codes used for the article "Clonal memory of cell division diverges between human healthy haematopoiesis and acute myeloid leukaemia".\
Article's authors: Alessandro Donada; Gurvan Hermange; Tiziana Tocci; Adil Midoun; Giulio Prevedello; Louisa Hadj Abed; Delia Dupré; Wenjie Sun; Idan Milo; Sabrina Tenreira Bento; Costandina Pospori; Andrew Innes; Christophe Willekens; Jacques Vargaftig; David Michonneau; Cristina Lo Celso; Nicolas Servant; Ken R Duffy; Hervé Isambert; Paul-Henri Cournede; Lucie Laplane; Leila Perié\
Click here to read the article: ["Clonal memory of cell division diverges between human healthy haematopoiesis and acute myeloid leukaemia"](https://www.biorxiv.org/content/10.1101/2025.06.24.660535v1)

The folder "Code" contains the code used for the analysis.\
The folder "Datasets" contains the original data.\
The folder "Output" contains the files obtained during the analysis, such as the MIIC input.\
The folder "Figures" is the folder where the figure will be saved by running the code. Only the figures concerning the family UMAP have been uploaded.

The analysis is divided in three steps, for each of them a code is available to replicate the results of the article:
1. Cell fate annotation and obtain figures concerning annotation, gene expression, division times
2. Cosine distance computation and obtain figures concerning the cosine distance analysis
3. Run feature selection for MIIC 

After point 3, the MIIC network is obtained with the [MIIC server](https://miic.curie.fr/)\
Each code returns intermediate files, stored in the "Output" folder. The files are already present in the "Output" folder in case you want to replicate only one part of the analysis.

R Version used: 4.3.2

📦 R Packages Used \
This project relies on the following R packages. The specific versions listed were used during the analysis to ensure reproducibility.

| Package       | Version        |
|---------------|----------------|
| matrixStats   | 1.1.0          |
| miic          | 1.9.0          |
| patchwork     | 1.1.3          |
| pals          | 1.8            |
| plotly        | 4.10.4.9000    |
| Seurat        | 4.4.0          |
| tidyverse     | 2.0.0          |

🔧 Installation\
To install all required packages (if not already installed), you can use the following R code:

```r
required_packages <- c(
"matrixStats", "miic", "patchwork",
  "pals", "plotly", "Seurat", "tidyverse"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
```


