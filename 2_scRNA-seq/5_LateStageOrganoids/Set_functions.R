

#Useful functions for sc analysis

kk2foldchange <- function(kk){
  GO <- as.data.frame(kk)
  GeneRatio <- strsplit(GO$GeneRatio,"/")
  GRdf <- data.frame()
  for (i in seq_along(GeneRatio)){
    GRdf[i,1] <- GeneRatio[[i]][1]
    GRdf[i,2] <- GeneRatio[[i]][2]
  }
  GRdf$V1 <- as.numeric(GRdf$V1)
  GRdf$V2 <- as.numeric(GRdf$V2)
  GO$GeneRatio2 <- GRdf$V1/GRdf$V2
  BgRatio <- strsplit(GO$BgRatio,"/")
  BgRdf <- data.frame()
  for (i in seq_along(BgRatio)){
    BgRdf[i,1] <- BgRatio[[i]][1]
    BgRdf[i,2] <- BgRatio[[i]][2]
  }
  BgRdf$V1 <- as.numeric(BgRdf$V1)
  BgRdf$V2 <- as.numeric(BgRdf$V2)
  GO$BgRatio2 <- BgRdf$V1/BgRdf$V2
  GO$foldchange <- GO$GeneRatio2/GO$BgRatio2
  return(GO)
}

samplePseudoBulkNN_v2 <- function(grouping, knn.object, mat, percSub = 0.1, nSample = 200, aggr = "mean") {
  nPerGroup <- table(grouping)
  nSamplePerGroup <- ceiling(nPerGroup * percSub)
  
  sampleIdxL <- lapply(
    X = names(nSamplePerGroup),
    FUN = function(gn) {
      idx <- which(grouping == gn)
      cell_names <- rownames(knn.object$nn.idx)[idx]  # Use cell names
      
      if (length(cell_names) == 0) {
        warning(paste("No valid cells for group:", gn))
        return(NULL)
      }
      
      N <- min(c(length(cell_names), nSample))
      
      rr <- lapply(
        1:nSamplePerGroup[gn],
        FUN = function(i) {
          random.cell <- sample(cell_names, 1)
          
          if (!random.cell %in% rownames(knn.object$nn.idx)) {
            stop(paste("Invalid cell:", random.cell))
          }
          
          nn <- sampleNearestNeighbors(
            cell = random.cell,
            knn.object = knn.object,
            N = N,
            output = "indices"
          )
          return(nn)
        }
      )
      return(rr)
    }
  )
  
  sampleIdxL <- unlist(sampleIdxL, recursive = FALSE)
  
  aggrFun <- switch(
    aggr,
    "mean" = function(X) rowMeans(X, na.rm = TRUE),
    "sum" = function(X) rowSums(X, na.rm = TRUE),
    stop("Unknown aggregation function.")
  )
  
  resM <- do.call(
    "cbind",
    lapply(sampleIdxL, function(x) {
      aggrFun(mat[, x, drop = FALSE])
    })
  )
  
  return(list(matrix = resM, sampleIdxL = sampleIdxL, grouping = grouping))
}


sampleNearestNeighbors <- function(cell, knn.object, N = 200, output = "indices") {
  # Convert cell name to index
  cell_idx <- match(cell, rownames(knn.object$nn.idx))
  
  if (is.na(cell_idx)) {
    stop(paste("Cell not found in KNN object:", cell))
  }
  
  if (output == "indices") {
    cells <- knn.object$nn.idx[cell_idx, ]
  } else {
    cells <- knn.object$nn.cells[knn.object$nn.idx[cell_idx, ]]
  }
  
  # Debugging the sampled cells
  #print(paste("Number of neighbors:", length(cells)))
  sample(cells, size = N, replace = FALSE)
}

plot_feature_genes <- function(seurat_obj, gene_list, reduction = "umap", pt.size = 0.5, palette_func) {
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Get available genes in the Seurat object
  available_genes <- rownames(seurat_obj)
  
  # Filter out genes that are not present
  valid_genes <- gene_list[gene_list %in% available_genes]
  missing_genes <- setdiff(gene_list, valid_genes)
  
  # Print warning if some genes are missing
  if (length(missing_genes) > 0) {
    message("The following genes were not found in the Seurat object and will be skipped: ", 
            paste(missing_genes, collapse = ", "))
  }
  
  # Generate FeaturePlots only for valid genes
  plots <- lapply(valid_genes, function(gene) {
    p <- FeaturePlot(seurat_obj, features = gene, reduction = reduction, pt.size = pt.size) &
      NoAxes() & NoLegend()
    p <- p + scale_color_gradientn(colors = palette_func(10))
    return(p)
  })
  
  # If no valid genes remain, return NULL
  if (length(plots) == 0) {
    message("No valid genes were found in the Seurat object. Returning NULL.")
    return(NULL)
  }
  
  # Arrange plots in a 4-column grid
  patchwork_obj <- wrap_plots(plots, ncol = 4)
  return(patchwork_obj)
}

plot_GO_barplot <- function(go_results, comparison, ontology, color_palette, n_terms = 5, xlim_range = NULL) {
  # Convert enrichment results to dataframe
  go_df <- as.data.frame(go_results)
  
  # Check if there are any enriched terms
  if (nrow(go_df) == 0) {
    print(paste("No enriched GO terms for", ontology, "in", comparison))
    return(NULL)
  }
  
  # Select top `n_terms` based on Fold Enrichment
  go_top <- go_df %>%
    arrange(desc(FoldEnrichment)) %>%
    slice_head(n = n_terms)
  
  # Compute color scale range only on the selected terms
  color_range <- range(-log10(go_top$p.adjust), na.rm = TRUE)
  
  # Ensure correct ordering for the plot
  go_top$Description <- factor(go_top$Description, levels = rev(go_top$Description))
  
  # Generate colors from the custom palette
  colors <- color_palette(10)
  
  # Create the barplot
  p <- ggplot(go_top, aes(x = FoldEnrichment, y = Description, fill = -log10(p.adjust))) +
    geom_col() +
    scale_fill_gradientn(colors = colors, name = "-log10(p.adjust)", limits = color_range) +
    theme_minimal(base_size = 10) +
    xlab("Fold Enrichment") +
    ylab("") +
    ggtitle(paste(ontology, "-", comparison))# +
  #theme(legend.position = "none")  # Remove legend
  
  # Apply uniform x-axis limits if provided
  if (!is.null(xlim_range)) {
    p <- p + xlim(xlim_range)
  }
  
  return(p)
}

plot_KEGG_barplot <- function(kegg_results, comparison, color_palette, n_terms = 5, xlim_range = NULL) {
  # Convert enrichment results to dataframe
  kegg_df <- as.data.frame(kegg_results)
  
  # Check if there are any enriched pathways
  if (nrow(kegg_df) == 0) {
    print(paste("No enriched KEGG pathways for", comparison))
    return(NULL)
  }
  
  # Select top `n_terms` based on Fold Enrichment (or GeneRatio if FoldEnrichment is missing)
  kegg_top <- kegg_df %>%
    arrange(desc(FoldEnrichment)) %>%  # If FoldEnrichment exists
    slice_head(n = n_terms)
  
  # Compute color scale range only on the selected terms
  color_range <- range(-log10(kegg_top$p.adjust), na.rm = TRUE)
  
  # Ensure correct ordering for the plot
  kegg_top$Description <- factor(kegg_top$Description, levels = rev(kegg_top$Description))
  
  # Generate colors from the custom palette
  colors <- color_palette(10)
  
  # Create the barplot
  p <- ggplot(kegg_top, aes(x = FoldEnrichment, y = Description, fill = -log10(p.adjust))) +
    geom_col() +
    scale_fill_gradientn(colors = colors, name = "-log10(p.adjust)", limits = color_range) +
    theme_minimal(base_size = 10) +
    ylab("") +
    ggtitle(paste("KEGG ", "-", comparison))
  #theme(legend.position = "none")  # Remove legend
  
  # Apply uniform x-axis limits if provided
  if (!is.null(xlim_range)) {
    p <- p + xlim(xlim_range)
  }
  
  return(p)
}

#Useful functions for sc analysis

kk2foldchange <- function(kk){
  GO <- as.data.frame(kk)
  GeneRatio <- strsplit(GO$GeneRatio,"/")
  GRdf <- data.frame()
  for (i in seq_along(GeneRatio)){
    GRdf[i,1] <- GeneRatio[[i]][1]
    GRdf[i,2] <- GeneRatio[[i]][2]
  }
  GRdf$V1 <- as.numeric(GRdf$V1)
  GRdf$V2 <- as.numeric(GRdf$V2)
  GO$GeneRatio2 <- GRdf$V1/GRdf$V2
  BgRatio <- strsplit(GO$BgRatio,"/")
  BgRdf <- data.frame()
  for (i in seq_along(BgRatio)){
    BgRdf[i,1] <- BgRatio[[i]][1]
    BgRdf[i,2] <- BgRatio[[i]][2]
  }
  BgRdf$V1 <- as.numeric(BgRdf$V1)
  BgRdf$V2 <- as.numeric(BgRdf$V2)
  GO$BgRatio2 <- BgRdf$V1/BgRdf$V2
  GO$foldchange <- GO$GeneRatio2/GO$BgRatio2
  return(GO)
}

samplePseudoBulkNN_v2 <- function(grouping, knn.object, mat, percSub = 0.1, nSample = 200, aggr = "mean") {
  nPerGroup <- table(grouping)
  nSamplePerGroup <- ceiling(nPerGroup * percSub)
  
  sampleIdxL <- lapply(
    X = names(nSamplePerGroup),
    FUN = function(gn) {
      idx <- which(grouping == gn)
      cell_names <- rownames(knn.object$nn.idx)[idx]  # Use cell names
      
      if (length(cell_names) == 0) {
        warning(paste("No valid cells for group:", gn))
        return(NULL)
      }
      
      N <- min(c(length(cell_names), nSample))
      
      rr <- lapply(
        1:nSamplePerGroup[gn],
        FUN = function(i) {
          random.cell <- sample(cell_names, 1)
          
          if (!random.cell %in% rownames(knn.object$nn.idx)) {
            stop(paste("Invalid cell:", random.cell))
          }
          
          nn <- sampleNearestNeighbors(
            cell = random.cell,
            knn.object = knn.object,
            N = N,
            output = "indices"
          )
          return(nn)
        }
      )
      return(rr)
    }
  )
  
  sampleIdxL <- unlist(sampleIdxL, recursive = FALSE)
  
  aggrFun <- switch(
    aggr,
    "mean" = function(X) rowMeans(X, na.rm = TRUE),
    "sum" = function(X) rowSums(X, na.rm = TRUE),
    stop("Unknown aggregation function.")
  )
  
  resM <- do.call(
    "cbind",
    lapply(sampleIdxL, function(x) {
      aggrFun(mat[, x, drop = FALSE])
    })
  )
  
  return(list(matrix = resM, sampleIdxL = sampleIdxL, grouping = grouping))
}


sampleNearestNeighbors <- function(cell, knn.object, N = 200, output = "indices") {
  # Convert cell name to index
  cell_idx <- match(cell, rownames(knn.object$nn.idx))
  
  if (is.na(cell_idx)) {
    stop(paste("Cell not found in KNN object:", cell))
  }
  
  if (output == "indices") {
    cells <- knn.object$nn.idx[cell_idx, ]
  } else {
    cells <- knn.object$nn.cells[knn.object$nn.idx[cell_idx, ]]
  }
  
  # Debugging the sampled cells
  #print(paste("Number of neighbors:", length(cells)))
  sample(cells, size = N, replace = FALSE)
}

plot_feature_genes <- function(seurat_obj, gene_list, reduction = "umap", pt.size = 0.5, palette_func) {
  #DefaultAssay(seurat_obj) <- "RNA"
  
  # Get available genes in the Seurat object
  available_genes <- rownames(seurat_obj)
  
  # Filter out genes that are not present
  valid_genes <- gene_list[gene_list %in% available_genes]
  missing_genes <- setdiff(gene_list, valid_genes)
  
  # Print warning if some genes are missing
  if (length(missing_genes) > 0) {
    message("The following genes were not found in the Seurat object and will be skipped: ", 
            paste(missing_genes, collapse = ", "))
  }
  
  # Generate FeaturePlots only for valid genes
  plots <- lapply(valid_genes, function(gene) {
    p <- FeaturePlot(seurat_obj, features = gene, reduction = reduction, pt.size = pt.size) &
      NoAxes() & NoLegend()
    p <- p + scale_color_gradientn(colors = palette_func(10))
    return(p)
  })
  
  # If no valid genes remain, return NULL
  if (length(plots) == 0) {
    message("No valid genes were found in the Seurat object. Returning NULL.")
    return(NULL)
  }
  
  # Arrange plots in a 4-column grid
  patchwork_obj <- wrap_plots(plots, ncol = 4)
  return(patchwork_obj)
}

plot_GO_barplot <- function(go_results, comparison, ontology, color_palette, n_terms = 5, xlim_range = NULL) {
  # Convert enrichment results to dataframe
  go_df <- as.data.frame(go_results)
  
  # Check if there are any enriched terms
  if (nrow(go_df) == 0) {
    print(paste("No enriched GO terms for", ontology, "in", comparison))
    return(NULL)
  }
  
  # Select top `n_terms` based on Fold Enrichment
  go_top <- go_df %>%
    arrange(desc(FoldEnrichment)) %>%
    slice_head(n = n_terms)
  
  # Compute color scale range only on the selected terms
  color_range <- range(-log10(go_top$p.adjust), na.rm = TRUE)
  
  # Ensure correct ordering for the plot
  go_top$Description <- factor(go_top$Description, levels = rev(go_top$Description))
  
  # Generate colors from the custom palette
  colors <- color_palette(10)
  
  # Create the barplot
  p <- ggplot(go_top, aes(x = FoldEnrichment, y = Description, fill = -log10(p.adjust))) +
    geom_col() +
    scale_fill_gradientn(colors = colors, name = "-log10(p.adjust)", limits = color_range) +
    theme_minimal(base_size = 10) +
    xlab("Fold Enrichment") +
    ylab("") +
    ggtitle(paste(ontology, "-", comparison))# +
  #theme(legend.position = "none")  # Remove legend
  
  # Apply uniform x-axis limits if provided
  if (!is.null(xlim_range)) {
    p <- p + xlim(xlim_range)
  }
  
  return(p)
}

plot_KEGG_barplot <- function(kegg_results, comparison, color_palette, n_terms = 5, xlim_range = NULL) {
  # Convert enrichment results to dataframe
  kegg_df <- as.data.frame(kegg_results)
  
  # Check if there are any enriched pathways
  if (nrow(kegg_df) == 0) {
    print(paste("No enriched KEGG pathways for", comparison))
    return(NULL)
  }
  
  # Select top `n_terms` based on Fold Enrichment (or GeneRatio if FoldEnrichment is missing)
  kegg_top <- kegg_df %>%
    arrange(desc(FoldEnrichment)) %>%  # If FoldEnrichment exists
    slice_head(n = n_terms)
  
  # Compute color scale range only on the selected terms
  color_range <- range(-log10(kegg_top$p.adjust), na.rm = TRUE)
  
  # Ensure correct ordering for the plot
  kegg_top$Description <- factor(kegg_top$Description, levels = rev(kegg_top$Description))
  
  # Generate colors from the custom palette
  colors <- color_palette(10)
  
  # Create the barplot
  p <- ggplot(kegg_top, aes(x = FoldEnrichment, y = Description, fill = -log10(p.adjust))) +
    geom_col() +
    scale_fill_gradientn(colors = colors, name = "-log10(p.adjust)", limits = color_range) +
    theme_minimal(base_size = 10) +
    ylab("") +
    ggtitle(paste("KEGG ", "-", comparison))
  #theme(legend.position = "none")  # Remove legend
  
  # Apply uniform x-axis limits if provided
  if (!is.null(xlim_range)) {
    p <- p + xlim(xlim_range)
  }
  
  return(p)
}


#Function based on 006_Integratino&UMAPsweep_2024-12-11
umap_worker <- function(params, df, genes) {
  dims <- params$dims
  neighbors <- params$neighbors
  dist <- params$dist
  
  cat("Running UMAP with dims =", dims, "neighbors =", neighbors, "min_dist =", dist, "\n")
  
  tmp <- RunUMAP(df, dims = 1:dims, n.neighbors = neighbors, min.dist = dist, seed.use = 42)
  
  # UMAP DimPlot
  p0 <- DimPlot(tmp, reduction = "umap") + 
    ggtitle(paste("Dims:", dims, "Neighbors:", neighbors, "Min Dist:", dist))
  print(p0)
  
  # FeaturePlot for specific markers
  DefaultAssay(tmp) <- "integrated"
  p3 <- FeaturePlot(
    tmp, 
    features = c("DCX", "SOX2", "MECP2", "FOXG1", "EMX1", "DLX2", "SATB2", "TBR1", 
                 "BCL11B", "GFAP", "S100B", "EOMES", "OTX2", "LHX9", "TFAP2A"), 
    reduction = "umap", ncol = 4)
  print(p3)
  
  # FeaturePlot for custom genes
  p <- FeaturePlot(tmp, features = genes, ncol = 5)
  print(p)
  
  # Clear tmp to free memory
  rm(tmp)
  gc()
}

# Main function
run_umap_sweep <- function(df, dims_list, neighbors_list, dist_list, genes) {
  # Create a list of parameter combinations
  param_grid <- expand.grid(dims = dims_list, neighbors = neighbors_list, dist = dist_list)
  
  # Run UMAP sweeps in parallel using pblapply
  pblapply(seq_len(nrow(param_grid)), function(i) {
    umap_worker(param_grid[i, ], df, genes)
  }, cl = parallel::detectCores() - 2) # Use all cores but one
}


markers_short <- c(
  "NES", "SOX2",          # NPC marker
  "NHLH1",                # neuroblast (new-born neuron) marker
  "DCX", "MAP2", "MAPT",  # neuron marker
  "FOXG1",                # telencephalon marker
  "EMX1", "EMX2",         # dorsal telencephalon (cortical) marker
  "EOMES",                # cortical intermediate progenitor (IP, proliferating neuroblast) marker
  "NEUROD6", "SLC17A7",   # dorsal telencephalic (cortical) glutamatergic neuron marker
  "BCL11B", "TBR1",              # deeper layer cortical neuron marker
  "SATB2",                # upper layer cortical neuron marker
  "RELN",                 # Cajal-Retzius cell marker
  "DLX2", "DLX5",         # ganglionic eminence (GE) marker
  "ISL1",                 # lateral ganglionic eminence (LGE) inhibitory neuron marker
  "NKX2-1",               # medial ganglionic eminence (MGE) inhibitory neuron marker
  "RSPO3", "TCF7L2", "LHX5", "LHX9", # diencephalon marker (for different neuron subtypes)
  "OTX2", "LMX1A", "EN1", # midbrain marker (for different neuron subtypes)
  "CYP26A1",              # Purkinje cell (cerebellar) progenitor marker
  "TFAP2A", "CA8",        # Purkinje cell marker
  "HOXB2", "HOXB5",       # hindbrain (medulla/pons) and spinal cord marker
  "SLC17A6",              # glutamatergic neuron marker
  "SLC32A1", "GAD1", "GAD2", # GABAergic neuron marker
  "TH",                   # dopaminergic neuron marker
  "CHAT", "ACHE",         # cholinergic neuron marker
  "TTR",                  # choroid plexus marker
  "GFAP", "AQP4", "S100B", # astrocyte marker
  "OLIG1",                # oligodendrocyte precursor cell marker
  "MBP", "SOX10",         # oligodendrocyte marker
  "SOX10",                # neural crest derivative marker
  "AIF1",                 # microglia marker
  "CLDN5",                # endothelial cell marker
  "DCN",                  # mesenchymal cell marker
  "MKI67"                 # cell cycle G2M phase marker (proliferative cells)
)


optimize_css_umap <- function(
    seurat_obj,
    dims_use_range = list(2:30, c(2:6, 8:30)),
    cluster_resolution_range = seq(0.2, 0.8, 0.2),
    n_neighbors_range = seq(10, 50, 10),
    min_dist_range = seq(0.1, 0.5, 0.2),
    label_tag = "orig.ident",
    output_file = "grid_search_results.csv"
) {
  results_list <- list()
  counter <- 1
  
  # Loop over each parameter combination explicitly
  for (dims in dims_use_range) {
    dims_numeric <- as.numeric(dims)
    for (cluster_res in cluster_resolution_range) {
      for (n_neighbors in n_neighbors_range) {
        for (min_dist in min_dist_range) {
          cat(sprintf("Running: dims_use=%s, cluster_res=%.2f, n_neighbors=%d, min_dist=%.2f\n",
                      toString(dims_numeric), cluster_res, n_neighbors, min_dist))
          res <- tryCatch({
            seurat_copy <- seurat_obj  # Make a copy to avoid modifying the original
            DefaultAssay(seurat_copy) <- "ATAC"
            # Apply SimSpec clustering
            seurat_copy <- cluster_sim_spectrum(
              seurat_copy,
              label_tag = label_tag,
              use_dr = "lsi",
              dims_use = dims_numeric,
              cluster_resolution = cluster_res,
              reduction.name = "css_atac_temp",
              reduction.key = "CSSATAC_"
            )
            
            # Run UMAP using as many dimensions as provided
            seurat_copy <- RunUMAP(
              seurat_copy,
              reduction = "css_atac_temp",
              dims = 1:length(dims_numeric),
              n.neighbors = n_neighbors,
              min.dist = min_dist,
              reduction.name = "umap_css_atac_temp",
              reduction.key = "UMAPCSSATAC_"
            )
            
            # Extract UMAP embedding and clusters
            umap_embedding <- Embeddings(seurat_copy, "umap_css_atac_temp")
            clusters <- seurat_copy$seurat_clusters
            
            # Debugging: print the dimensions of the embedding and the number of clusters
            cat("UMAP embedding dimensions:", dim(umap_embedding), 
                "Number of cluster assignments:", length(clusters), "\n")
            
            # Compute silhouette score
            silhouette_scores <- cluster::silhouette(as.numeric(factor(clusters)), dist(umap_embedding))
            avg_silhouette <- mean(silhouette_scores[, 3])
            
            rm(seurat_copy)
            gc(reset = TRUE)
            
            data.table(
              avg_silhouette = avg_silhouette,
              dims_use = toString(dims_numeric),
              cluster_res = cluster_res,
              n_neighbors = n_neighbors,
              min_dist = min_dist
            )
          }, error = function(e) {
            cat("Error for dims_use:", toString(dims_numeric),
                "cluster_res:", cluster_res,
                "n_neighbors:", n_neighbors,
                "min_dist:", min_dist, "\n")
            cat("Error message:", e$message, "\n")
            data.table(
              avg_silhouette = NA,
              dims_use = toString(dims_numeric),
              cluster_res = cluster_res,
              n_neighbors = n_neighbors,
              min_dist = min_dist
            )
          })
          
          results_list[[counter]] <- res
          counter <- counter + 1
        }
      }
    }
  }
  
  # Combine all results and save to CSV
  results <- rbindlist(results_list)
  fwrite(results, file = output_file)
  
  # Find best parameters (if any valid silhouette scores exist)
  best_params <- results[!is.na(avg_silhouette)][order(-avg_silhouette)][1]
  
  gc(reset = TRUE)
  # Optional Linux-specific memory clearing:
  # system("sync; echo 3 > /proc/sys/vm/drop_caches")
  
  return(list(results = results, best_params = best_params))
}



library(simspec)
library(Seurat)
library(cluster)
library(data.table)

optimize_integration_umap <- function(
    seurat_obj,
    dims_use = 2:30,
    cluster_resolution_range = seq(0.4, 1.0, 0.2),
    n_neighbors_range = seq(10, 50, 10),
    min_dist_range = seq(0.1, 0.5, 0.1)
) {
  param_grid <- expand.grid(
    cluster_resolution = cluster_resolution_range,
    n_neighbors = n_neighbors_range,
    min_dist = min_dist_range
  )
  
  results_list <- list()  # Store results before combining
  
  for (i in seq_len(nrow(param_grid))) {
    params <- param_grid[i, ]
    cluster_res <- params$cluster_resolution
    n_neighbors <- params$n_neighbors
    min_dist <- params$min_dist
    
    cat(sprintf("Testing: cluster_res=%.2f, n_neighbors=%d, min_dist=%.2f\n", 
                cluster_res, n_neighbors, min_dist))
    
    result <- tryCatch({
      seurat_copy <- seurat_obj  # Prevent modifying the original object
      
      # Apply SimSpec clustering
      seurat_copy <- cluster_sim_spectrum(
        seurat_copy,
        label_tag = "orig.ident",
        use_dr = "lsi",
        dims_use = dims_use,
        cluster_resolution = cluster_res,
        reduction.name = "css_atac",
        reduction.key = "CSSATAC_"
      )
      
      # Run UMAP
      seurat_copy <- RunUMAP(
        seurat_copy,
        reduction = "css_atac",
        dims = 1:ncol(Embeddings(seurat_copy, "css_atac")),
        n.neighbors = n_neighbors,
        min.dist = min_dist,
        reduction.name = "umap_css_atac",
        reduction.key = "UMAPCSSATAC_"
      )
      
      # Compute silhouette score
      umap_embedding <- Embeddings(seurat_copy, "umap_css_atac")
      silhouette <- silhouette(
        x = as.numeric(factor(seurat_copy$orig.ident)),
        dist = dist(umap_embedding)
      )
      avg_silhouette <- mean(silhouette[, 3])
      
      rm(seurat_copy)
      gc(reset = TRUE)
      
      data.table(cluster_res, n_neighbors, min_dist, avg_silhouette)
    }, error = function(e) {
      data.table(cluster_res, n_neighbors, min_dist, avg_silhouette = NA)
    })
    
    results_list[[i]] <- result
  }
  
  # Combine results
  results <- rbindlist(results_list)
  
  # Find best parameters
  best_params <- results[!is.na(avg_silhouette)][order(-avg_silhouette)][1]
  
  gc(reset = TRUE)
  
  
  return(list(results = results, best_params = best_params))
}


#library(lisi)

grid_search_umap_lisi <- function(
    seurat_obj,
    dims_use_range = list(2:30, c(2:6, 8:30)),
    cluster_resolution_range = seq(0.2, 0.8, 0.2),
    n_neighbors_range = seq(10, 50, 10),
    min_dist_range = seq(0.1, 0.5, 0.2),
    label_tag = "orig.ident",
    output_file = "umap_lisi_grid_search.csv"
) {
  results_list <- list()
  counter <- 1
  
  for (dims in dims_use_range) {
    dims_numeric <- as.numeric(dims)
    for (cluster_res in cluster_resolution_range) {
      for (n_neighbors in n_neighbors_range) {
        for (min_dist in min_dist_range) {
          cat(sprintf("Running: dims_use=%s, cluster_res=%.2f, n_neighbors=%d, min_dist=%.2f\n",
                      toString(dims_numeric), cluster_res, n_neighbors, min_dist))
          
          res <- tryCatch({
            # Make a copy so the original object remains unchanged
            seurat_copy <- seurat_obj
            DefaultAssay(seurat_copy) <- "ATAC"
            
            # Run the clustering step (using your custom cluster_sim_spectrum function)
            seurat_copy <- cluster_sim_spectrum(
              seurat_copy,
              label_tag = label_tag,
              use_dr = "lsi",
              dims_use = dims_numeric,
              cluster_resolution = cluster_res,
              reduction.name = "css_atac_temp",
              reduction.key = "CSSATAC_"
            )
            
            # Run UMAP on the new reduction
            seurat_copy <- RunUMAP(
              seurat_copy,
              reduction = "css_atac_temp",
              dims = 1:length(dims_numeric),
              n.neighbors = n_neighbors,
              min.dist = min_dist,
              reduction.name = "umap_css_atac_temp",
              reduction.key = "UMAPCSSATAC_"
            )
            
            # Extract UMAP embeddings and create metadata
            embedding <- Embeddings(seurat_copy, "umap_css_atac_temp")
            meta_data <- data.frame(orig.ident = seurat_copy[[label_tag]], 
                                    row.names = colnames(seurat_copy))
            
            # Compute LISI scores. Note: the computed column name will match the label_colnames.
            lisi_scores <- compute_lisi(embedding, meta_data, label_colnames = label_tag)
            # Rename for clarity
            names(lisi_scores) <- "LISI"
            
            # Compute summary statistics (using median as the primary metric)
            median_lisi <- median(lisi_scores$LISI, na.rm = TRUE)
            mean_lisi   <- mean(lisi_scores$LISI, na.rm = TRUE)
            
            data.table(
              dims_use = toString(dims_numeric),
              cluster_res = cluster_res,
              n_neighbors = n_neighbors,
              min_dist = min_dist,
              median_lisi = median_lisi,
              mean_lisi = mean_lisi
            )
          }, error = function(e) {
            cat("Error for dims_use:", toString(dims_numeric), 
                "cluster_res:", cluster_res, 
                "n_neighbors:", n_neighbors, 
                "min_dist:", min_dist, "\n")
            cat("Error message:", e$message, "\n")
            data.table(
              dims_use = toString(dims_numeric),
              cluster_res = cluster_res,
              n_neighbors = n_neighbors,
              min_dist = min_dist,
              median_lisi = NA,
              mean_lisi = NA
            )
          })
          
          results_list[[counter]] <- res
          counter <- counter + 1
        }
      }
    }
  }
  
  # Combine results into one data.table and write to CSV
  results <- rbindlist(results_list)
  fwrite(results, file = output_file)
  
  # Choose the best combination based on the highest median LISI (if available)
  best_params <- results[!is.na(median_lisi)][order(-median_lisi)][1]
  
  gc(reset = TRUE)
  
  return(list(results = results, best_params = best_params))
}


QC_VlnPlot <- function(seurat_obj, hline_values = c(1000, 2000, 10)) {
  features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  
  # Generate violin plots without combining them into one plot
  vln_plots <- VlnPlot(seurat_obj, features = features, ncol = 3, 
                       cols = "#025e8d", pt.size = 0.1, combine = FALSE)
  
  # Add jitter, horizontal lines, and annotations to each plot
  vln_plots <- lapply(seq_along(vln_plots), function(i) {
    vln_plots[[i]] + 
      geom_jitter(width = 0.2, size = 0.5, alpha = 0.2, color = "#363633") +
      geom_hline(yintercept = hline_values[i], linetype = "dashed", 
                 color = "#be1818", size = 0.7) +
      annotate("text", x = 0.5, 
               y = hline_values[i] + (0.25 * hline_values[i]), 
               label = paste0(hline_values[i]), color = "#be1818", size = 4)
  })
  
  # Arrange plots in a single row using patchwork and return the result
  wrap_plots(vln_plots, ncol = 3)
}

QC_VlnPlot_atac <- function(seurat_obj, hline_values = c(3000, 2, 0.05, 4, 15)) {
  features <- c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks')
  
  # Generate violin plots without combining them into one plot
  vln_plots <- VlnPlot(seurat_obj, features = features, ncol = 5, 
                       cols = "#00a69d", pt.size = 0.1, combine = FALSE)
  
  # Adjust blacklist_ratio and pct_reads_in_peaks scaling (log10 transformation for better visualization if needed)
  vln_plots <- lapply(seq_along(vln_plots), function(i) {
    vln_plots[[i]] + 
      geom_jitter(width = 0.2, size = 0.5, alpha = 0.2, color = "#363633") +
      geom_hline(yintercept = hline_values[i], linetype = "dashed", 
                 color = "#D65A00", size = 0.7) +
      annotate("text", x = 0.6, 
               y = hline_values[i] + (0.05 * max(layer_data(vln_plots[[i]])$y, na.rm = TRUE)), 
               label = paste0(hline_values[i]), color = "#D65A00", size = 4) +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  })
  
  # Arrange plots in a single row using patchwork and return the result
  wrap_plots(vln_plots, ncol = 5)
}


library(ggplot2)
library(patchwork)

QC_HistPlot <- function(seurat_obj) {
  # List of features to plot
  features <- c("nFeature_RNA", "nCount_RNA")
  
  # Function to generate the histogram plot for a single feature
  make_hist_plot <- function(feature) {
    # Extract the data from the metadata
    dat <- seurat_obj@meta.data[[feature]]
    
    # Define breaks in steps of 100
    brks <- seq(0, max(dat) + 100, by = 100)
    
    # Compute the histogram (without plotting)
    hist_data <- hist(dat, breaks = brks, plot = FALSE)
    
    # Log-transform the counts (adding 1 to avoid log10(0))
    log_counts <- log10(hist_data$counts + 1)
    
    # Smooth the log-transformed histogram using a spline
    fit <- smooth.spline(x = hist_data$mids, y = log_counts, spar = 0.5)
    smoothed <- fit$y
    
    # Identify local minima in the smoothed curve
    # Local minima are identified as positions where the sign of the second difference is 2.
    local_min_idx <- which(diff(sign(diff(smoothed))) == 2) + 1
    # Convert indices to x-values using histogram midpoints
    local_min_vals <- hist_data$mids[local_min_idx]
    
    # Create a data frame for plotting
    df_hist <- data.frame(mid = hist_data$mids,
                          log_counts = log_counts)
    
    # Determine a y position for annotations (90% of maximum log-count)
    label_y <- max(df_hist$log_counts) * 0.9
    
    # Build the plot
    p <- ggplot(df_hist, aes(x = mid, y = log_counts)) +
      geom_col(fill = "#025e8d") +
      labs(x = feature, y = "Log10(Frequency + 1)", 
           title = paste("Log-Transformed Histogram of", feature)) +
      theme_minimal()
    
    # Add vertical lines and annotations for each local minimum
    for(i in seq_along(local_min_vals)){
      p <- p +
        geom_vline(xintercept = local_min_vals[i], linetype = "dashed", 
                   color = "#be1818", size = 1) +
        annotate("text", x = local_min_vals[i], y = label_y, 
                 label = paste("min:", round(local_min_vals[i], 1)), 
                 color = "#be1818", vjust = -0.5, angle = 90, size = 3)
    }
    return(p)
  }
  
  # Generate a list of plots for the two features
  plot_list <- lapply(features, make_hist_plot)
  
  # Combine the plots into a patchwork object (side-by-side)
  combined_plot <- wrap_plots(plot_list, ncol = 2)
  
  return(combined_plot)
}


library(Seurat)
library(ggplot2)
library(patchwork)

QC_ScatterPlots <- function(seurat_obj, 
                            low_nFeature = 1500, 
                            high_nFeature = 9000, 
                            mt_threshold = 10) {
  # Extract metadata for convenience
  md <- seurat_obj@meta.data
  
  # Determine extremes for annotations
  max_nFeature <- max(md$nFeature_RNA, na.rm = TRUE)
  max_percent  <- max(md$percent.mt, na.rm = TRUE)
  min_percent  <- min(md$percent.mt, na.rm = TRUE)
  
  # Calculate counts for each region based on thresholds
  # Lower row: percent.mt < mt_threshold
  reg1 <- sum(md$nFeature_RNA < low_nFeature & md$percent.mt < mt_threshold, na.rm = TRUE)
  reg2 <- sum(md$nFeature_RNA >= low_nFeature & md$nFeature_RNA < high_nFeature & md$percent.mt < mt_threshold, na.rm = TRUE)
  reg3 <- sum(md$nFeature_RNA >= high_nFeature & md$percent.mt < mt_threshold, na.rm = TRUE)
  # Upper row: percent.mt >= mt_threshold
  reg4 <- sum(md$nFeature_RNA < low_nFeature & md$percent.mt >= mt_threshold, na.rm = TRUE)
  reg5 <- sum(md$nFeature_RNA >= low_nFeature & md$nFeature_RNA < high_nFeature & md$percent.mt >= mt_threshold, na.rm = TRUE)
  reg6 <- sum(md$nFeature_RNA >= high_nFeature & md$percent.mt >= mt_threshold, na.rm = TRUE)
  
  # Create p1: FeatureScatter of nFeature_RNA vs. percent.mt
  p1 <- FeatureScatter(seurat_obj, 
                       feature1 = "nFeature_RNA", 
                       feature2 = "percent.mt", 
                       cols = "#025e8d") +
    # Add horizontal and vertical dashed lines
    geom_hline(yintercept = mt_threshold, linetype = "dashed", color = "#be1818", size = 1) +
    geom_vline(xintercept = low_nFeature, linetype = "dashed", color = "#be1818", size = 1) +
    geom_vline(xintercept = high_nFeature, linetype = "dashed", color = "#be1818", size = 1) +
    # Annotations for threshold lines
    annotate("text", 
             x = max_nFeature * 0.95, 
             y = mt_threshold + 3, 
             label = paste("y =", mt_threshold), 
             color = "#be1818", size = 4, hjust = 1) +
    annotate("text", 
             x = low_nFeature + low_nFeature, 
             y = min_percent + 50, 
             label = paste("x =", low_nFeature), 
             color = "#be1818", size = 4, vjust = -1) +
    annotate("text", 
             x = high_nFeature - low_nFeature, 
             y = min_percent + 50, 
             label = paste("x =", high_nFeature), 
             color = "#be1818", size = 4, vjust = -1) +
    # Annotations for cell counts in each region
    # Lower row (percent.mt < mt_threshold)
    annotate("text", 
             x = low_nFeature / 2, 
             y = (min_percent + mt_threshold) / 2, 
             label = reg1, 
             color = "black", size = 4) +
    annotate("text", 
             x = (low_nFeature + high_nFeature) / 2, 
             y = (min_percent + mt_threshold) / 2, 
             label = reg2, 
             color = "#ffcc00", size = 4) +
    annotate("text", 
             x = (high_nFeature + max_nFeature) / 2, 
             y = (min_percent + mt_threshold) / 2, 
             label = reg3, 
             color = "black", size = 4) +
    # Upper row (percent.mt >= mt_threshold)
    annotate("text", 
             x = low_nFeature / 2, 
             y = (mt_threshold + max_percent) / 2, 
             label = reg4, 
             color = "black", size = 4) +
    annotate("text", 
             x = (low_nFeature + high_nFeature) / 2, 
             y = (mt_threshold + max_percent) / 2, 
             label = reg5, 
             color = "black", size = 4) +
    annotate("text", 
             x = (high_nFeature + max_nFeature) / 2, 
             y = (mt_threshold + max_percent) / 2, 
             label = reg6, 
             color = "black", size = 4)
  
  # Create p2: Scatter plot of log10(nCount_RNA) vs. log10(nFeature_RNA)
  # Extract metadata and compute log10-transformed values (adding 1 to avoid log10(0))
  md <- seurat_obj@meta.data
  md$log_nCount_RNA <- log10(md$nCount_RNA + 1)
  md$log_nFeature_RNA <- log10(md$nFeature_RNA + 1)
  
  p2 <- ggplot(md, aes(x = log_nCount_RNA, y = log_nFeature_RNA)) +
    geom_point(color = "#025e8d", alpha = 0.7) +
    theme_minimal() +
    labs(x = "log10(nCount_RNA + 1)", y = "log10(nFeature_RNA + 1)",
         title = "Scatter: log10 Transformed") +
    geom_smooth(method = "lm", se = FALSE, color = "#be1818")
  
  # Combine the two plots side-by-side using patchwork and return the patchwork object
  return(p1 + p2)
}

library(dplyr)
library(ggplot2)
library(ggrepel)
library(rlang)

make_volcano <- function(
    df,
    motif_col       = "motif",        # name of the motif/feature column
    tf_col          = "TF",           # name of the TF column
    fc_col          = "avg_log2FC",   # name of the log2 fold‐change column
    pval_col        = "p_val_adj",
    col_up          = "#D00000",
    col_down        = "#025E8D",
    top_n           = 10,             # how many up/down to label if no highlights
    highlight_genes = NULL,           ### NEW: vector of genes (labels) to force‐highlight
    fc_cutoff       = 0.5,            # vertical cutoff lines at ±fc_cutoff
    p_cutoff        = 0.05,           # horizontal cutoff line at p_cutoff
    title           = NULL,           # optional plot title
    save_path       = NULL,           # if non‐NULL, will ggsave() here
    width           = 5,
    height          = 4,
    dpi             = 300
) {
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  
  # 1. Tidy & compute new columns
  dat <- df %>%
    rename(
      motif      = !!sym(motif_col),
      TF         = !!sym(tf_col),
      avg_log2FC = !!sym(fc_col),
      p_val_adj  = !!sym(pval_col)
    ) %>%
    mutate(
      negLog10padj = -log10(p_val_adj),
      label        = ifelse(is.na(TF) | TF == "", motif, TF)
    )
  
  # 2. Decide which to highlight
  if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
    # use only the user‐supplied list
    to_up   <- dat %>% filter(label %in% highlight_genes, avg_log2FC >  0)
    to_down <- dat %>% filter(label %in% highlight_genes, avg_log2FC <  0)
    dat <- dat %>%
      mutate(selectLab = ifelse(label %in% highlight_genes, label, ""))
  } else {
    # fallback to top_n logic
    to_up   <- dat %>% filter(avg_log2FC >  0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = top_n)
    to_down <- dat %>% filter(avg_log2FC <  0) %>% arrange( avg_log2FC)        %>% slice_head(n = top_n)
    dat <- dat %>%
      mutate(selectLab = ifelse(label %in% c(to_up$label, to_down$label), label, ""))
  }
  
  # 3. Symmetric x‐axis limit
  max_fc <- max(abs(dat$avg_log2FC), na.rm = TRUE) + 1
  
  # 4. Build the plot
  p <- ggplot(dat, aes(x = avg_log2FC, y = negLog10padj)) +
    geom_point(color = "grey60", size = 1) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dotted") +
    geom_hline(yintercept = -log10(p_cutoff),           linetype = "dotted") +
    # highlighted points
    geom_point(data = to_up,   aes(x = avg_log2FC, y = negLog10padj), color = col_up,   size = 1) +
    geom_point(data = to_down, aes(x = avg_log2FC, y = negLog10padj), color = col_down, size = 1) +
    # labels
    geom_text_repel(
      data          = subset(dat, selectLab != ""),
      aes(label      = selectLab),
      size          = 4,
      box.padding   = 0.2,
      point.padding = 0.2,
      segment.size  = 0.5,
      segment.color = "#363633",
      max.overlaps  = Inf
    ) +
    scale_x_continuous(limits = c(-max_fc, max_fc)) +
    labs(
      title = title,
      x     = "log2 Fold Change",
      y     = "-log10 Adjusted P-value"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      legend.position = "none",
      axis.title      = element_text(size = 10),
      axis.text       = element_text(size = 8),
      panel.grid      = element_blank()
    )
  
  # 5. Save if requested
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
  }
  
  return(p)
}

make_volcano_v2 <- function(
    data,
    # for Seurat objects: specify ident.1, ident.2, assay, min.pct, log2.threshold
    ident.1         = NULL,
    ident.2         = NULL,
    assay           = NULL,
    min.pct         = 0.1,
    log2.threshold  = 0.25,
    # for both data.frame and Seurat outputs:
    motif_col       = "gene",        # feature column; will pull from rownames if missing
    fc_col          = "avg_log2FC",   # log2 fold‐change column
    pval_col        = "p_val_adj",    # adjusted p‐value column
    col_up          = "#D00000",
    col_down        = "#025E8D",
    top_n           = 10,              # top_n labels if no highlight_genes
    highlight_genes = NULL,            # vector of genes to highlight instead of top_n
    fc_cutoff       = 0.5,
    p_cutoff        = 0.05,
    title           = NULL,
    save_path       = NULL,
    width           = 5,
    height          = 4,
    dpi             = 300
) {
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(rlang)
  library(tibble)
  
  # 1. Get marker dataframe
  if (inherits(data, "Seurat")) {
    if (is.null(ident.1) || is.null(ident.2)) {
      stop("ident.1 and ident.2 must be provided for Seurat objects")
    }
    df <- FindMarkers(
      object         = data,
      ident.1        = ident.1,
      ident.2        = ident.2,
      assay          = assay,
      min.pct        = min.pct,
      log2.threshold = log2.threshold
    )
  } else {
    df <- data
  }
  
  # 2. Ensure gene column exists
  if (!motif_col %in% colnames(df)) {
    df <- df %>% rownames_to_column(var = motif_col)
  }
  
  # 3. Rename and compute stats
  dat <- df %>%
    rename(
      motif      = !!sym(motif_col),
      avg_log2FC = !!sym(fc_col),
      p_val_adj  = !!sym(pval_col)
    ) %>%
    mutate(
      negLog10padj = -log10(p_val_adj),
      label        = motif,
      # flag significant points
      signif = (p_val_adj < p_cutoff) & (abs(avg_log2FC) > fc_cutoff)
    )
  
  # 4. Determine labels (only among significant)
  if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
    dat <- dat %>%
      mutate(selectLab = ifelse(label %in% highlight_genes & signif, label, ""))
  } else {
    top_up <- dat %>% filter(signif & avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC)) %>% slice_head(n = top_n)
    top_down <- dat %>% filter(signif & avg_log2FC < 0) %>%
      arrange(avg_log2FC) %>% slice_head(n = top_n)
    dat <- dat %>%
      mutate(selectLab = ifelse(label %in% c(top_up$label, top_down$label), label, ""))
  }
  
  # 5. Split data for plotting
  to_up_lab   <- dat %>% filter(selectLab != "" & avg_log2FC >  0)
  to_down_lab <- dat %>% filter(selectLab != "" & avg_log2FC <  0)
  to_sig_uns  <- dat %>% filter(signif & selectLab == "")
  to_ns       <- dat %>% filter(!signif)
  
  # 6. Axis limits
  max_fc <- max(abs(dat$avg_log2FC), na.rm = TRUE) + 1
  
  # 7. Plot
  p <- ggplot() +
    # non-significant points (below threshold) in light gray
    geom_point(data = to_ns, aes(x = avg_log2FC, y = negLog10padj),
               color = "grey90", size = 1) +
    # significant but unlabeled in dark gray
    geom_point(data = to_sig_uns, aes(x = avg_log2FC, y = negLog10padj),
               color = "grey70", size = 1) +
    # labeled points in input colors
    geom_point(data = to_up_lab,   aes(x = avg_log2FC, y = negLog10padj),
               color = col_up,   size = 1) +
    geom_point(data = to_down_lab, aes(x = avg_log2FC, y = negLog10padj),
               color = col_down, size = 1) +
    # threshold lines
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dotted") +
    geom_hline(yintercept = -log10(p_cutoff),           linetype = "dotted") +
    # labels for selected genes only
    geom_text_repel(
      data          = subset(dat, selectLab != ""),
      aes(x = avg_log2FC, y = negLog10padj, label = selectLab),
      size          = 6,
      box.padding   = 0.2,
      point.padding = 0.2,
      segment.size  = 0.5,
      segment.color = "#363633",
      max.overlaps  = Inf
    ) +
    scale_x_continuous(limits = c(-max_fc, max_fc)) +
    labs(
      title = title,
      x     = "log2 Fold Change",
      y     = "-log10 Adjusted P-value"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      legend.position = "none",
      axis.title      = element_text(size = 10),
      axis.text       = element_text(size = 8),
      panel.grid      = element_blank()
    )
  
  # 8. Save if requested
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = p, width = width, height = height, dpi = dpi)
  }
  
  return(p)
}



