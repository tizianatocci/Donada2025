# File to store all routines for the project



#===============================================================================
# General routines, used in all/ most of the files. 
#===============================================================================

# Load a dataset from a .csv file and create a Seurat object
# without metadata.
load_data <- function(path_data){
  # load matrix genes-cells
  # file_name <- path_data
  dataset <- read.csv(path_data, header = TRUE)
  cts <- subset(dataset, select = -c(X))
  name_cells <- colnames(cts)
  name_genes <- dataset$X
  output_list <- list("obj"=dataset,
                      "nm_cells" = name_cells)
  return(output_list)
}

# Load metadata from a .csv file.
load_metadata <- function(path_metadata, name_cells){
  #metadata_name <- sprintf(path_metadata)
  m_csv <- read.csv(path_metadata, header = TRUE)
  names_metadata <- gsub("-",".", m_csv$X)
  m_csv$X <- names_metadata
  correct_order <- match(names_metadata, name_cells)
  new_ordered <- m_csv[order(correct_order),]
  rownames(new_ordered) <- new_ordered$X;
  dataset_metadata<- subset(new_ordered, select = -c(X))
  return(dataset_metadata)
}

# Create average column for dT2, dT3
create_avg_dT <- function(dataset){
  if (!'avg_dT2' %in% colnames(dataset)){
    dataset$avg_dT2 <- rowMeans(cbind(dataset$dT2.1, dataset$dT2.2))}
  if (!'avg_dT3' %in% colnames(dataset)){
  dataset$avg_dT3 <- rowMeans(cbind(dataset$dT3.1, dataset$dT3.2,
                                    dataset$dT3.3, dataset$dT3.4))}
  return(dataset)
}

# Create variance column for dT2, dT3
create_var_dT <- function(dataset){
  if (!'var_dT2' %in% colnames(dataset)){
    dataset$var_dT2 <- rowVars(cbind(dataset$dT2.1, dataset$dT2.2))}
  if (!'var_dT3' %in% colnames(dataset)){
    dataset$var_dT3 <- rowVars(cbind(dataset$dT3.1, dataset$dT3.2,
                                     dataset$dT3.3, dataset$dT3.4))}
  return(dataset)
}

# Create max_var and min_var columns
create_maxminvar <- function(dataset){
  if (!'max_var' %in% colnames(dataset)){
    dataset$max_var <- rowMaxs(cbind(dataset$var_dT2, dataset$var_dT3), na.rm=TRUE)
    dataset[is.infinite(dataset$max_var), c("max_var")] <- NA}
  
  if (!'min_var' %in% colnames(dataset)){
    dataset$min_var <- rowMins(cbind(dataset$var_dT2, dataset$var_dT3), na.rm=TRUE)
    dataset[is.infinite(dataset$min_var), c("min_var")] <- NA}
  
  return(dataset)
}

# Summary function n2: load data, create Seurat object, 
# load metadata, 
# create avg and var columns, create max_var, min_var column
summary_load_v2 <- function(path_data, path_metadata){
  output_load<-load_data(path_data)
  dataset_metadata <- load_metadata(path_metadata,output_load$nm_cells)
  
  #create avg and var columns for dT2 and dT3
  dataset_metadata <- create_avg_dT(dataset_metadata)
  dataset_metadata <- create_var_dT(dataset_metadata)
  dataset_metadata <- create_maxminvar(dataset_metadata)
  
  output_list <- list("obj" = output_load$obj,
                      "nm_cells" = dataset_metadata)
  return(output_list)
}

# Create Seurat object with metadata
seurat_object_wt_metadata <- function(dataset1, dataset1_metadata)
{
  seurat_realdata <- CreateSeuratObject(counts = subset(dataset1, select = -c(X)), 
                                        row.names = dataset1$X)
  seurat_realdata <- AddMetaData(object = seurat_realdata,
                                 metadata = dataset1_metadata)
  return(seurat_realdata)
}



#===============================================================================
# Routines specific for the file: "cellfate_annotation.R"
#===============================================================================

# Clustering for a dataset, running the SCT normalization
run_clustering <- function(dataset, ndims, 
                           n_neigh, k_param, min_dist, res_1, res_2, 
                           reduction_name_pca, reduction_name_umap){
  dataset <- SCTransform(dataset, verbose = FALSE)
  dataset <- RunPCA(dataset, verbose = FALSE, npcs=50, reduction.name=reduction_name_pca, approx=FALSE) #prima 70
  ep <- ElbowPlot(dataset, reduction=reduction_name_pca, ndims=50)
  print(ep)
  dataset <- FindNeighbors(dataset, reduction=reduction_name_pca,
                           dims = 1:ndims, k.param=k_param, verbose = FALSE)
  dataset <- FindClusters(dataset, verbose = FALSE, resolution=res_1)
  dataset <- FindClusters(dataset, verbose = FALSE, resolution=res_2)
  dataset <- RunUMAP(dataset, dims = 1:ndims, reduction=reduction_name_pca,
                     n.neighbors=n_neigh, min.dist=min_dist, 
                     reduction.name=reduction_name_umap, verbose = FALSE)
  return(dataset)}

# Label cells as Neutrophiles if they express ELANE
find_elane <- function(dataset, set_assay, thresh=20){
  DefaultAssay(dataset) <- set_assay
  if (set_assay == 'RNA'){
    dataset <- NormalizeData(dataset)}
  elane_dataset <-WhichCells(object = dataset, expression = ELANE > 0, slot = 'data')
  if (set_assay=='RNA'){
    elane_expr <- dataset@assays$RNA@data["ELANE",colnames(dataset) %in% elane_dataset]}
  else if (set_assay=='SCT'){
    elane_expr <- dataset@assays$SCT@data["ELANE",colnames(dataset) %in% elane_dataset]}
  
  elane_cells <- WhichCells(object = dataset, expression = ELANE > thresh/100*max(elane_expr), slot = 'data')
  output <- list("elane_dataset"=elane_dataset, "elane_expr"=elane_expr, "elane_cells"=elane_cells)
  return(output)}

# Harmonize annotation between the one used in the article and CITESeq
harmonize_citeseq_ann <- function(exp.realdata){
  exp.realdata$citeseq.ann.adj <- NA
  exp.realdata$citeseq.ann.adj <- ifelse(exp.realdata$citeseq.ann == "Granulocytes_Neutrophiles", G_name, exp.realdata$citeseq.ann.adj)
  exp.realdata$citeseq.ann.adj <- ifelse(exp.realdata$citeseq.ann == "Erythrocytes" | exp.realdata$citeseq.ann =="Megakaryocyte", EM_name, exp.realdata$citeseq.ann.adj)
  exp.realdata$citeseq.ann.adj <- ifelse(exp.realdata$citeseq.ann == "CD34+" | exp.realdata$citeseq.ann =="CD34+_CD38+", U_name, exp.realdata$citeseq.ann.adj)
  exp.realdata$citeseq.ann.adj <- ifelse(exp.realdata$citeseq.ann == "Monocytes", "Mo/DC-P", exp.realdata$citeseq.ann.adj)
  return(exp.realdata)
}


#===============================================================================
# Routines specific for the file: "cosinedistance_computation.R"
#===============================================================================

# Find cell types for each family
find_celltype_family <- function(unique_fam, stem.combined){
df_types <- data.frame()
k <- 1
for (single_fam in unique_fam){
  cell_types<-stem.combined@meta.data[stem.combined$family %in% single_fam,]$predicted.celltype
  label_types <- paste(cell_types, collapse=", ")
  df_types[k, "family"] <- single_fam
  df_types[k, "cell.fates"] <- label_types
  k <- k+1
}
return(df_types)}

# Function to compute combinations
find_comb <- function(unique_fam, stem.combined){
  all_comb <- list()
  for (idx_f in 1:length(unique_fam)){
    single_fam <- unique_fam[idx_f]
    idx_fam <- which(stem.combined$family==single_fam)
    #idx_fam <- idx_fam - idx_trasl
    #names_fam <- gsub(".*_","",names(idx_fam))
    names_fam <- names(idx_fam)
    combinations <- combn(names_fam, 2)
    if (length(combinations)==2){
      combinations <- matrix(combinations)
    }
    all_comb[[single_fam]] <- (combinations)
  }
  return(all_comb)
}

# Function to calculate cosine distance between two vectors
cosine_distance <- function(first_vector, second_vector) {
  # compute dot product
  dot_product <- sum(first_vector * second_vector)
  # compute norms
  first_norm <- sqrt(sum(first_vector^2))
  second_norm <- sqrt(sum(second_vector^2))
  # compute cosine distance
  output_cosinedistance <- 1 - (dot_product / (first_norm * second_norm))
  return(output_cosinedistance)
}

# Fuction to compute median distance family
compute_md_cosine <- function(combinations, which_slot='counts') {
  dist_family <- c()
  for (num_cells in 1:dim(combinations)[2]) {
    #print(num_cells)
    if (which_slot == 'data') {
      vectorA <-
        stem.combined@assays$RNA@data[, combinations[1, num_cells]]
      vectorB <-
        stem.combined@assays$RNA@data[, combinations[2, num_cells]]
    }
    else{
      vectorA <-
        stem.combined@assays$RNA@counts[, combinations[1, num_cells]]
      vectorB <-
        stem.combined@assays$RNA@counts[, combinations[2, num_cells]]
    }
    dist_family[num_cells] <- cosine_distance(vectorA, vectorB)
  }
  md <- median(dist_family)
  return(md)
}

# Add distance to the metadata
add_distance_metadata <- function(name_file, file_tosave, md_df_cosine){
  dataset <- read.csv(name_file, header=T)
  dataset$family[is.na(dataset$family)]= 'missing_family'
  dataset_tosave <- merge(dataset, md_df_cosine, by='family') %>% select(union(colnames(dataset), "distance"))
  new_order <- match(dataset_tosave$X, dataset$X)
  dataset_tosave <- dataset_tosave[order(new_order),]
  dataset_tosave$family[(dataset_tosave$family)=='missing_family']= NA
  write.csv(dataset_tosave, file_tosave, row.names=FALSE)
}


# Choose random families, same size and eventually same flag type
create_random_families<-function(unique_fam, stem.combined, 
                                 flag_type, which_slot='counts'){
  # count number of cells per family
  cts <- dplyr::count(stem.combined@meta.data, family)
  random_all_md <- c()
  i <- 1
  for (name_f in unique_fam){
    print(name_f)
    # find size of family
    size_f <- cts[cts$family==name_f, "n"]
    # find names of cells inside family
    names_all <- rownames(stem.combined@meta.data[stem.combined$family %in% name_f, ])
    # choose n=size_f random cell to replace and 
    # replace the with n random cells 
    # if flag_type=T, 
    # replace with a cell of the same cell type 
    # according to the annotation, without repetitions
    # (i.e.: I cannot choose the same cell twice)
    names_replacer <- c()
    if (flag_type) {
      j <- 1
      for (one_name in names_all){
        # type to replace for each cell
        type_toreplace <- stem.combined@meta.data[rownames(stem.combined@meta.data) %in% one_name,]$predicted.celltype
        # look inside the same cell type
        name_choice_replacer <- rownames(stem.combined@meta.data[(stem.combined$predicted.celltype == type_toreplace),])
        name_choice_replacer <- setdiff(name_choice_replacer, names_replacer)
        names_replacer[j] <- sample(name_choice_replacer, 1)
        j <- j+1
      }}
    else {
      names_replacer <- sample(rownames(stem.combined@meta.data), size_f)
    }
    # check I do not pick the same cell twice
    stopifnot(length(unique(names_replacer))==length(names_replacer))
    # find combinations
    combinations <- combn(names_replacer, 2)
    if (length(combinations)==2){
      combinations <- matrix(combinations)
    }
    random_md <- compute_md_cosine(combinations, which_slot)
    random_all_md[i] <- random_md
    i <- i+1}
  return(random_all_md)
}

# repetitions
obtain_random_families <- function(num_it, 
                                   unique_fam, stem.combined, 
                                   flag_type, which_slot, m_seed){
  set.seed(m_seed)
  random_mixed_all <- data.frame()
  for (n_it in 1:num_it){
    random_mixed <- create_random_families(unique_fam, stem.combined, 
                                           flag_type, which_slot)
    random_mixed_all[1:30,n_it] <- random_mixed
  }
  return(random_mixed_all)
}

# obtain plots
obtain_random_plot <- function(random_mixed_wtrep, 
                               r_weights, x_limit, y_limit,
                               h_color, p_title){
  plot_mixed_wtrep <- ggplot(random_mixed_wtrep,
                             aes(x=distance, weights=r_weights)) +
    geom_histogram(binwidth=0.05, color='white', fill=h_color) +
    labs(title=p_title) +
    scale_x_continuous(limits=c(0,x_limit)/10, breaks=seq(0,x_limit)/10)+
    scale_y_continuous(limits=c(0,y_limit), breaks=seq(0,y_limit))
  return(plot_mixed_wtrep)}



#===============================================================================
# Routines specific for the file: "featureselection.R"
#===============================================================================

# Compute homogeneity
compute_homogeneity <- function(df_family)
{
  # all names of families; number of families; number of cells per family
  true_name_families <- unique(df_family$family)
  num_families <- length(true_name_families)
  cts <- dplyr::count(df_family, family)
  # plot joint distribution
  #p <- ggplot(df_family, aes(family))+geom_bar(aes(fill = predicted.celltype))
  #p+scale_y_continuous(breaks = seq(0, max(cts$n)+1, by=1), limits=c(0,max(cts$n)+1))
  #browser()
  # find metric of homogeneity
  celltypes_per_family <- df_family %>% 
    dplyr::count(family, predicted.celltype)
  columns <- c("family","homogeneity", "n_celltypes") 
  df <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
  colnames(df) <- columns
  for (i  in(1: length(true_name_families))){
    name_fam<- true_name_families[i]
    #print(name_fam)
    df[i, "family"] <- name_fam
    # exclude NAs?
    #celltypes_per_family <- celltypes_per_family[!is.na(celltypes_per_family$predicted.celltype),]
    smatrix <- celltypes_per_family[celltypes_per_family$family %in% name_fam, c("n", "predicted.celltype")]
    # numerator: number of cells for most frequent predicted.celltype, excluding NA
    num_hom <- smatrix[!is.na(smatrix$predicted.celltype), ]$n
    #browser()
    # denominator: number of cells for the whole family
    df[i, "homogeneity"] <- max(num_hom)/sum(smatrix$n)
    #print(df[i, "homogeneity"])
    df[i, "n_celltypes"] <- length(num_hom)
    #print(df[i, "n_celltypes"])
  }
  #browser()
  df[is.infinite(df$homogeneity), c("homogeneity")] <- NA
  df[is.infinite(df$n_celltypes), c("n_celltypes")] <- NA
  #if (sum(df$family %in% "missing_family") !=0)
  #{
  #df[df$family %in% "missing_family",]$homogeneity<- NA
  #df[df$family %in% "missing_family",]$n_celltypes<- NA}
  df_merge <- merge(df_family,df,by.x="family") 
  #browser()
  #df_merge <- select(df_merge, -c(family))
  #browser()
  return(df_merge)
}


# Compute average of group of genes (MT, RPL, RPS, MRPL, MRPS)
compute_avg_mitoribo <- function(stem.combined, all_groups){
  names_genes <- rownames(stem.combined)
  #stem.obj <- stem.combined
  avg_genes <- c()
  # for each group of genes, average
  for (single_group in all_groups){
    print(single_group)
    names_group <- names_genes[grepl(single_group, names_genes)]
    print(names_group)
    submatrix_group <- data.frame(stem.combined@assays$RNA$counts[names_group,])
    # do mean for each cell, it means per column and save it in list
    avg_genes[[sub('\\^', '', single_group)]] <- colMeans(submatrix_group)
  }
  # check if colnames are the same for each averaged group
  stopifnot(identical(colnames(avg_genes[[1]]), colnames(avg_genes[[2]])))
  stopifnot(identical(colnames(avg_genes[[1]]), colnames(avg_genes[[3]])))
  stopifnot(identical(colnames(avg_genes[[1]]), colnames(avg_genes[[4]])))
  stopifnot(identical(colnames(avg_genes[[1]]), colnames(avg_genes[[5]])))
  # delete genes who are mithocondrial or ribosomial
  names_genes_todelete <- names_genes[grepl(paste(all_groups, collapse="|"), names_genes)]
  new_df<- stem.combined@assays$RNA@counts[!(row.names(stem.combined@assays$RNA@counts)) %in% names_genes_todelete,]
  new_df <- data.frame(new_df)
  # add the means obtained for mithocondrial and robisomial genes
  rows_avg <- do.call(rbind, avg_genes)
  df_output <- rbind(new_df, rows_avg)
  # check
  print(setdiff(rownames(df_output), rownames(new_df)))
  
  return(df_output)
}

# From the specified genes, go back to the raw matrix
go_back_raw <- function(dataset, name_top_genes_i)
{
  index_top_i <- rownames(dataset@assays$RNA@counts) %in% name_top_genes_i
  truncated_cts_i<- dataset@assays$RNA@counts[index_top_i== TRUE,]
  identical(rownames(truncated_cts_i), sort(name_top_genes_i))
  truncated_cts_df_i <- as.data.frame(t(as.data.frame(truncated_cts_i)))
  return(truncated_cts_df_i)
}

# Save the input of TF-MIIC
save_input <- function(my_TF, stem.combined, df.merged, 
                       vars, vars2, flag_fam, n_try){
  # go back to raw matrix
  sub_matrix <- go_back_raw(stem.combined, vars)
  #browser()
  # add metadata to row matrix
  stopifnot(identical(rownames(df.merged), rownames(sub_matrix)))
  data_frame_i <- merge(sub_matrix, df.merged, 
                        by = 'row.names', all = TRUE)
  rownames(data_frame_i) <- data_frame_i$Row.names
  data_frame_i <- data_frame_i %>% select(-c("Row.names"))
  df_input <- data_frame_i
  
  # for layout
  df_layout <- df_input %>% select(all_of(union(vars2, my_TF)))
  
  # create additional file
  # supplementary file
  is_consequence <- rep(1, length(df_input))
  is_consequence[colnames(df_input) %in% my_TF] <- 0
  is_consequence[colnames(df_input) %in% vars2] <- 0
  #is_consequence[colnames(df_input[,-1]) %in% name_nc] <- 0
  group <- ifelse(is_consequence==0,"TF","Other")
  color <- ifelse(is_consequence==0,"F0C69C","A9EBA9")
  supp_file_i <- data.frame(colnames(df_input),
                            #as.numeric(lens_i), 
                            as.numeric(is_consequence), group, color)
  colnames(supp_file_i) <- c("var_names", 
                             #"var_type", 
                             "is_consequence", "group", "group_color")
  
  # check
  identical(sort(supp_file_i[supp_file_i$group=='TF', ]$var_names), sort(union(vars2, my_TF)))
  
  # save file
  write.csv(df_input, sprintf("filtering\\input_%s_%s.csv", n_try, flag_fam), row.names=FALSE)
  write.csv(supp_file_i, sprintf("filtering\\input_%s_%s_supp.csv", n_try, flag_fam), row.names=FALSE)
  # layout
  write.csv(df_layout, sprintf("filtering\\input_%s_%s_LAYOUT.csv", n_try, flag_fam), row.names=FALSE)
  
  return(df_input)
}


#-------------------------------------------------------------------------------
# seurat_load_tfs
#-------------------------------------------------------------------------------
seurat_load_tfs = function (tfs=NULL)
{
  if ( is.null (tfs) )
    tfs = read.table ("./Data/TFs.tsv", sep="\t", header=T)
  return (unlist (tfs[,1]) )
}

#-------------------------------------------------------------------------------
# seurat_dir
#-------------------------------------------------------------------------------
seurat_dir = function (seurat_step, sub_dir=NULL, pipeline=T, create_dir=T)
{
  dir = "."
  if (pipeline)
  {
    dir = paste0 ("./", SEURAT_STEPS[[seurat_step]])
    if (create_dir)
      dir.create (dir, showWarnings=F)
    if (!is.null (sub_dir))
    {
      dir = paste0 ("./", SEURAT_STEPS[[seurat_step]], "/", sub_dir)
      if (create_dir)
        dir.create (dir, showWarnings=F)
    }
  }
  return (dir)
}

#-----------------------------------------------------------------------------
# catn0
#-----------------------------------------------------------------------------
catn0 <- function  (...)
{
  cat (paste0 (..., "\n") )
}

#-------------------------------------------------------------------------------
# seurat_mi_compute
#-------------------------------------------------------------------------------
# Compute MI with a condition or gene of interest
# Inputs:
# - mat: matrix with UMI
# - value_name: can be a gene (in matrix) or a condition (not in matrix)
# - values: must be supplied for a condition (otherwise can be extracted from mat)
# - cell_types: used from the filename to store computed MI
# - tfs: TF list
# - recompute: if F, load the MI file and recompute only missing MI
# - display: verbose
# - pipeline: used for location of the MI file
# Returns:
# - a dataframe: MI with genes as rows and condition/marker genes as column
#-------------------------------------------------------------------------------
seurat_mi_compute = function (mat, value_name, values=NULL, cell_types=NULL,
                              tfs=NULL, recompute=F, display=T, pipeline=F)
{
  if (is.null(tfs))
    tfs = seurat_load_tfs (tfs)
  normal_gene = F
  if ( is.null (values) )
  {
    if ( ! (value_name %in% rownames(mat) ) )
      stop (value_name, " not found in data.")
    values = mat[value_name,]
    if ( ! (value_name %in% tfs) )
      normal_gene = T
  }
  
  df_mi = data.frame ("to_be_renamed" = rep (NA, nrow(mat)), stringsAsFactors=F)
  colnames (df_mi) = value_name
  rownames(df_mi) = rownames (mat)
  #
  # Prepare file load
  #
  dir = seurat_dir ("mi", "Files", pipeline=pipeline, create_dir=F)
  if ( is.null (cell_types) ){
    file_name = paste0 ( dir, "/MI.tsv")
  } else
  {
    cell_types_str = paste (cell_types, collapse="_")
    file_name = sanityze_filename (cell_types_str, partial_filename=T)
    file_name = paste0 ( dir, "/MI_", file_name, ".tsv")
  }
  #
  # If Some MI has been computed and saved, reload it
  #
  if ( (!is.null (file_name)) & file.exists(file_name) )
  {
    if (display)
      catn0 ("    +-- Loading MI file")
    df_mi = read.table (file_name, sep="\t", header=T)
    if (nrow (df_mi) == 0)
    {
      df_mi = data.frame ("to_be_renamed" = rep (NA, nrow(mat)), stringsAsFactors=F)
      colnames (df_mi) = value_name
      rownames(df_mi) = rownames (mat)
    }
    else
    {
      rownames(df_mi) = df_mi$gene
      df_mi$gene = NULL
      if ( ! (value_name %in% colnames (df_mi) ) )
        df_mi[,value_name] = NA
      if (recompute)
        df_mi[,value_name] = NA
      
      if (normal_gene)
        missing_genes =  rownames(mat) [ (rownames(mat) %in% tfs)
                                         & ( ! (rownames(mat) %in% rownames(df_mi) ) ) ]
      else
        missing_genes =  rownames(mat) [ ! (rownames(mat) %in% rownames(df_mi) ) ]
      if (length (missing_genes) > 0)
      {
        df_add = df_mi[1,]
        df_add[1,] = rep (NA, ncol(df_add))
        for (one_gene in missing_genes) # lent
        {
          rownames(df_add) = one_gene
          df_mi = rbind (df_mi, df_add)
        }
      }
    }
  }
  #
  # If MI already computed for all genes, nothing to do
  #
  if (normal_gene)
    genes_no_mi = rownames(df_mi) [ (rownames(df_mi) %in% tfs)
                                    & is.na (df_mi[,value_name])
                                    & (rownames (df_mi) != value_name)]
  else
    genes_no_mi = rownames(df_mi) [ is.na (df_mi[,value_name])
                                    & (rownames (df_mi) != value_name)]
  if (length (genes_no_mi) == 0)
  {
    if (display)
      catn0 ("    +-- All MI have already been computed")
    return (df_mi)
  }
  if (display)
    catn0 ("    +-- Computing MI with ", value_name, " for ", length(genes_no_mi), " genes")
  #
  # Run miic with 1000 genes and state_order with genes as consequence
  #
  mat = mat[genes_no_mi,]
  row_low = 1
  if (normal_gene)
    loop_range = length (tfs)
  else
    loop_range = 1000
  while (row_low <= nrow(mat))
  {
    row_max = min (nrow (mat), row_low + loop_range - 1)
    if (display)
      catn0 ("    +-- Computing MI with genes from ", row_low, " to ", row_max)
    df_loop = as.data.frame (t (as.matrix (mat[ c (row_low:row_max), ]) ) )
    df_st = data.frame ("var_names" = c (colnames(df_loop), "Metadata"),
                        "is_contextual" = c (rep(0, ncol(df_loop)), 1),
                        "is_consequence" = c (rep(1, ncol(df_loop)), 0),
                        stringsAsFactors = F)
    if (class(values) == 'data.frame'){
      values=as.vector(unlist(values))}
    df_loop$Metadata = values
    miic_res = suppressMessages (miic (input_data=df_loop, state_order=df_st,
                                       orientation=F, latent="no", n_threads=8) )
    miic_res = miic_res$all.edges.summary
    # miic_res = miic_res[,c("x","y","info_shifted")]
    rownames (miic_res) = NULL
    rownames (miic_res)[miic_res$x != "Metadata"] = miic_res[miic_res$x != "Metadata", "x"]
    rownames (miic_res)[miic_res$y != "Metadata"] = miic_res[miic_res$y != "Metadata", "y"]
    if (row_low == 1)
      df_miic = miic_res[, "info_shifted", F]
    else
      df_miic = rbind (df_miic, miic_res[, "info_shifted", F])
    row_low = row_low + loop_range
  }
  df_mi[rownames(df_miic), value_name] = round (df_miic$info_shifted / ncol(mat), 6)
  if (normal_gene)
    df_mi[ (rownames(df_mi) %in% tfs)
           & is.na (df_mi[,value_name]), value_name] = 0
  else
    df_mi[is.na (df_mi[,value_name]), value_name] = 0
  
  return (df_mi)
}


# feature selection based on info
filter_info <- function(df, th1, to_delete=NULL)
{
  idx <- df > th1
  #View(idx)
  idx_keep <- Reduce('|', data.frame(idx))
  #print('reduce per column')
  #print(idx_keep)
  # variables selected
  vars <- rownames(df)[idx_keep]
  if (!is.null(to_delete)){
    vars <- vars[!grepl(to_delete, vars)]}
  
  #print(vars)
  return(vars)
}

# check if the gene expression did not change from the original file
prova_check <- function(df_hom, dataset1, num_exp)
{
  prova1 <- df_hom[df_hom$exp %in% num_exp, ]
  common_names1 <- intersect(colnames(prova1), dataset1$X)
  a <- prova1[, colnames(prova1) %in% common_names1]
  b <- as.data.frame(t(as.data.frame(dataset1)))
  colnames(b) <- b[1,]
  b <- b[-1,]
  b <- b[, colnames(b) %in% common_names1]
  #setequal(a,b)
  # check identical, change order
  new_order <- match(rownames(b), gsub(".*_", "", rownames(df_hom)))
  ordered <- b[order(new_order),]
  setequal(colnames(a), colnames(ordered))
  for (name_col in colnames(a)){
    print(name_col)
    stopifnot(identical(as.numeric(a[, name_col]), 
                        as.numeric(ordered[, name_col])))
  }
  print("All gene expression is equal to the original file.")
}


#===============================================================================
# Routines specific for the file: "scrnaseq_code5_listgenes.R"
#===============================================================================

# Select genes based on on info or info shifted
select_variable <- function(df, name_var, filtering_var, thresh, 
                            top_number, name_metadata)
{
  #browser()
  # select variable based on name
  submatrix_variable <- df[(df$x %in% name_var| df$y %in% name_var),]
  # information must be > thresh
  submatrix_filtered <- submatrix_variable[submatrix_variable[,filtering_var] > thresh, ]
  # take top_n values
  #submatrix_output <- NA
  if (filtering_var=="info"){
    submatrix_output <- top_n(submatrix_filtered, n=top_number, wt=info)}
  else if (filtering_var=="info_shifted"){
    submatrix_output <- top_n(submatrix_filtered, n=top_number, wt=info_shifted)}
  output_var <- NA
  for (i in 1:nrow(submatrix_output)){
    if (submatrix_output[i,"x"] != name_var){
      output_var[i] <- submatrix_output[i,"x"]} 
    else if (submatrix_output[i,"y"] != name_var){
      output_var[i] <- submatrix_output[i,"y"]
    }
  }
  output_var <- setdiff(output_var, name_metadata)
  # names of variables
  #output_var <- union(submatrix_output$x[submatrix_output$x != name_var], 
  #                    submatrix_output$y[submatrix_output$y != name_var])
  #output_var <- setdiff(output_var, c("Colony_size", "T1", "avg_dT2", "avg_dT3",
  #"var_dT2", "var_dT3", "predicted.celltype", "exp"))
  output <- list("submatrix"= submatrix_output,
                 "name_vars"=output_var)
  return(output)
}

# Select genes based on info or info shifted and save file 
filter_and_write<-function(file_name, filtering_var, thresh, n_variables, name_metadata, features, correct_folder)
{
  # load summary
  summary_MIIC <- read.table(
    sprintf("MIIC_summary\\%s\\edgesList.miic.summary.txt", file_name),
    sep='\t',
    header=TRUE)
  # order df based on increasing mutual information
  summary <- summary_MIIC[order(summary_MIIC[,filtering_var], decreasing = TRUE),]
  # initialize list
  N <- length(features)
  total_list <- vector("list", N)
  # save submatrix and name of genes in list
  for (j in 1:N)
  {
    name_feat <- features[j]
    total_list[[j]] <-select_variable(summary, name_feat, 
                                      filtering_var, thresh, n_variables, name_metadata)
  }
  names(total_list) <- features
  # print everything in folder
  i <- 1
  for (element in total_list){
    print(i)
    current_name <- features[i]
    print(current_name)
    names_write <- setdiff(element$name_vars, name_metadata)
    tofilter_submatrix <- element$submatrix
    filt_submatrix <- tofilter_submatrix %>% select(c("x", "y", "type", "ai", "info", "Nxy_ai", "info_shifted", "partial_correlation"))
    write(names_write, 
          file = file.path("list genes", correct_folder,sprintf("names_%s.txt", current_name)), sep = "")
    write.csv(filt_submatrix, 
              file = file.path("list genes", correct_folder, sprintf("matrix_%s.csv", current_name)))
    i <- i+1
  }
}


#===============================================================================
# Routines specific for the file: "azimuth_annotation.R"
#===============================================================================

# Azimuth function to obtain Azimuth annotation.
run_azimuth <- function(
    # seurat object
  query, 
  # parameters: integration; findNeighbors
  mapping_score_k, k_weight, k_param){
  
  annotation.levels=NULL
  assay <- DefaultAssay(query)
  reference <- LoadReference("../References/Azimuth_bonemarrow")$map
  dims <- as.double(slot(reference, "neighbors")$refdr.annoy.neighbors@alg.info$ndim)
  if (is.null(annotation.levels)) {
    annotation.levels <- names(slot(object = reference, name = "meta.data"))
    annotation.levels <- annotation.levels[!grepl(pattern = "^nCount", x = annotation.levels)]
    annotation.levels <- annotation.levels[!grepl(pattern = "^nFeature", x = annotation.levels)]
    annotation.levels <- annotation.levels[!grepl(pattern = "^ori", x = annotation.levels)]
  }
  
  query <- ConvertGeneNames(
    object = query,
    reference.names = rownames(x = reference),
    homolog.table = 'https://seurat.nygenome.org/azimuth/references/homologs.rds'
  )
  # Calculate nCount_RNA and nFeature_RNA if the query does not
  # contain them already
  if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
    calcn <- as.data.frame(x = Seurat:::CalcN(object = query[[assay]]))
    colnames(x = calcn) <- paste(
      colnames(x = calcn),
      assay,
      sep = '_'
    )
    query <- AddMetaData(
      object = query,
      metadata = calcn
    )
    rm(calcn)
  }
  
  # Calculate percent mitochondrial genes if the query contains genes
  # matching the regular expression "^MT-"
  if (any(grepl(pattern = '^MT-', x = rownames(x = query)))) {
    query <- PercentageFeatureSet(
      object = query,
      pattern = '^MT-',
      col.name = 'percent.mt',
      assay = assay
    )
  }
  # Find anchors between query and reference
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = assay,
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = rownames(Loadings(reference[["refDR"]])),
    dims = 1:dims,
    n.trees = 20,
    mapping.score.k = mapping_score_k,
    verbose = T
  )
  
  refdata <- lapply(X = annotation.levels, function(x) {
    reference[[x, drop = TRUE]]
  })
  names(x = refdata) <- annotation.levels
  
  query <- TransferData(
    reference = reference,
    query = query,
    #query.assay = assay,
    dims = 1:dims,
    anchorset = anchors,
    refdata = refdata,
    n.trees = 20,
    store.weights = TRUE,
    k.weight = k_weight,
    verbose = T
  )
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = reference,
    query = query,
    query.assay = assay,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE,
    verbose = T
  )
  # Calculate the query neighbors in the reference
  # with respect to the integrated embeddings
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(reference[["refDR"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE,
    k.param=k_param,
    verbose = T
  )
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  NNTransform <- function(
    object,
    meta.data,
    neighbor.slot = "query_ref.nn",
    key = 'ori.index'
  ) {
    on.exit(expr = gc(verbose = FALSE))
    ind <- Indices(object[[neighbor.slot]])
    ori.index <- t(x = sapply(
      X = 1:nrow(x = ind),
      FUN = function(i) {
        return(meta.data[ind[i, ], key])
      }
    ))
    rownames(x = ori.index) <- rownames(x = ind)
    slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
    return(object)
  }
  query <- NNTransform(
    object = query,
    meta.data = reference[[]]
  )
  # Project the query to the reference UMAP.
  query[["myumap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = reference[["refUMAP"]],
    reduction.key = 'UMAP_',
    verbose = T
  )
  # Calculate mapping score and add to metadata
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors, 
                            ndim = dims,
                            kanchors = 5,
                            ksmooth = 20,
                            ksnn = 10),
    col.name = "mapping.score"
  )
  return(query)
}

# Obtain metadata with citeseq annotation and Azimuth annotation
metadata_citeseq_azimuth <-
  function(path_metadata, path_output, query, compare_ann) {
    # load metadata dataset
    query$name_well <- rownames(query)
    dt1_firstfile <- load_metadata(path_metadata, query$name_well)
    dt1_metadata_azimuth <-
      merge(dt1_firstfile, compare_ann, by = "row.names")
    names(dt1_metadata_azimuth)[names(dt1_metadata_azimuth) == 'Row.names'] <-
      'X'
    names(dt1_metadata_azimuth)[names(dt1_metadata_azimuth) == 'predicted.celltype'] <-
      'citeseq.ann'
    write.csv(dt1_metadata_azimuth,
              path_output,
              row.names = FALSE)
    return(dt1_metadata_azimuth)
  }





