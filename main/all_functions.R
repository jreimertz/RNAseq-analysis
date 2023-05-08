#!/usr/bin/Rscript
## Author: Justin Reimertz

# Load in libraries
library(shiny)
library(bslib)
library(DT)
library(BiocManager)
library(biomaRt)
library(data.table)
library(R.utils)
library(glue)
library(colourpicker)
library(pheatmap)
library(ggbeeswarm)
library(DESeq2)
library(limma)
library(edgeR)
library(pheatmap)
library(igraph)
library(tidyverse)
theme_set(theme_bw())

#### Loading and processing data ####
#' Load Normalized Count Data
#'
#' @param filepath A text string of the full filepath to the file to load.
#'
#' @return A dataframe containing the data loaded from the CSV in `filepath`. 
#'
#' @examples 
#' `data <- load_meta_data('data/norm_count_data.csv')`
load_meta_data <- function(input) {
  reactive({
    req(input$metaFile)
    return(as_tibble(fread(input$metaFile$datapath)))
  })
}

##### Making a summary table #####
#' Make a summary table from the tidy data frame
#' 
#' @param meta_data A tidy dataframe with one column of gene names (or ids) and
#'              columns with sample names containing the expression or
#'              count data for each sample
#' @return A dataframe with summary information including each column name
#'          and the mean and standard deviation for that column
#'          
#' @examples 
#'  `summary_table <- make_summary_table(data)`
make_summary_table <- function(meta_data) {
  out <- tryCatch(
    {
      # Save the column names
      data_cols <- meta_data %>%
        sapply(typeof) %>%
        enframe() %>%
        rename(`Column Name` = name, Type = value)
      # Calculate the mean and standard deviation for all numeric columns
      num_means <- meta_data %>% 
        select_if(is.numeric) %>% 
        # Calculate mean for each column
        colMeans(na.rm=T) %>% 
        enframe() %>%
        rename(`Column Name` = name, mean = value)
      num_sd <- meta_data %>%
        select_if(is.numeric) %>% 
        # Calculate standard deviation for each column
        apply(2, function(x) sd(x, na.rm=T)) %>% 
        enframe() %>% 
        rename(`Column Name` = name, std = value)
      # Combine numeric values into one column
      num_cols <- merge(num_means, num_sd) %>%
        mutate(`Mean (sd) or Distinct Values` = as.list(glue(
          "{format(mean, digits=3)} (+/-{format(std, digits=3)})"))) %>%
        select(`Column Name`, `Mean (sd) or Distinct Values`)
      # Get unique values for non-numeric values
      unique_vals <- meta_data %>%
        select_if(is.character) %>%
        apply(2, unique) %>% enframe() %>%
        rename(`Column Name` = name, `Mean (sd) or Distinct Values` = value)
      # Combine everything together and return the assembled summary table
      summary_table <- bind_rows(unique_vals, num_cols) %>%
        merge(data_cols, by="Column Name") %>%
        # Get the columns in the right order
        select(`Column Name`, Type, `Mean (sd) or Distinct Values`)
      return(summary_table)
    },
    error=function(cond) {
      return(NULL)
    },
    warning=function(cond) {
      return(NULL)
    }
  )
  return(out)
}

##### Choose a plot #####
#' Helper function to call the appropriate plot function to based on the user
#' inputted parameters
#' 
#' @param data data frame of meta data to build a plot with
#' @param column User inputted column to make a plot of
#' @param group User inputted column to group to the data by
#' @param plot_type Type of plot to make using the inputted data and parameters
#' 
#' @examples 
#'   `make_meta_plot(data, age_of_death, Diagnosis, "density")`
make_meta_plot <- function(meta_data, column, group, plot_type) {
  if (plot_type == "Histogram") {
    out <- meta_data %>% ggplot(aes(
      x=meta_data[[as.character(column)]],
      fill=meta_data[[as.character(group)]])) + geom_histogram(alpha=0.6) +
      labs(x = glue("{column}"),
           fill = glue("{group}")) +
      theme(legend.position = "bottom", text = element_text(size = 18))
  }
  else if (plot_type == "Density") {
    out <- meta_data %>% ggplot(aes(
      x=meta_data[[as.character(column)]],
      fill=meta_data[[as.character(group)]])) + geom_density(alpha=0.6) +
      labs(x = glue("{column}"),
           fill = glue("{group}")) + 
      theme(legend.position = "bottom", text = element_text(size = 18))
  }
  else {
    out <- meta_data %>% ggplot(aes(
      x=meta_data[[as.character(group)]], 
      y=meta_data[[as.character(column)]], 
      fill=meta_data[[as.character(group)]])) + geom_violin() +
      labs(x = glue("{group}"),
           y = glue("{column}"),
           fill = glue("{group}")) +
      theme(legend.position = "bottom", text = element_text(size = 18))
  }
  return(out)
}

##### Load in Normalized Counts data #####
#' Function to read in a normalized counts matrix file
#' and store it in a tibble
#' 
#' @param filepath A text string to the full filepath of the file to load
#' 
#' @return A tibble containing a normalized counts matrix 
#'
#'@examples
#' `load_count_data(input)`
load_count_data <- function(input) {
  reactive({
    req(input$countsFile)
    out <- tryCatch(
      {
        # Sometimes the first column name isn't read in correctly so rename that
        # column if it doesn't get read in correctly
        count_data <- as_tibble(fread(input$countsFile$datapath)) %>%
          dplyr::rename(gene = V1)
        return(count_data)
      },
      error=function(cond) {
        # If everything everything works appropriately just return in the dataframe
        count_data <- as_tibble(fread(input$countsFile$datapath)) 
        return(count_data)
      }
    )
    return(out)
  })
}

##### Apply filters to count data #####
#' Helper function for load_count_data() to apply filters based on user
#' specifications and return the filtered counts matrix
#' 
#' @param count_data (tibble) normalized counts data to apply filters to
#' @param variance (double) percentile of variance to filter genes by
#' @param sample_threshold (double) number of non-zero samples to filter genes by
#' 
#' @return A tibble containing the statistics of the original counts matrix and
#'          the counts matrix after filters have been applied
#'
#'@examples
#' `apply_filters(count_data, 1000, 5)`
apply_filters <- function(counts_data, var_threshold, sample_threshold) {
  # Store original data counts for samples and genes before filters are
  # applied to the data
  num_samples <- ncol(counts_data) - 1
  num_genes <- nrow(counts_data)
  # Apply filters to the count data
  filtered_counts <- get_filtered_counts(counts_data, var_threshold, 
                                         sample_threshold)
  # Get the number of genes after filtering the data and the number removed
  num_genes_new <- nrow(filtered_counts)
  num_genes_filtered <- num_genes - num_genes_new
  # Return a named list of all the information
  return(tibble(
    `Number of Samples` = format(num_samples, digits=2), 
    `Total Genes Before Filtering` = num_genes,
    `Total Genes Passing Filters` = num_genes_new,
    `Percent Genes Passing Filters` = glue(
      "{format((num_genes_new/num_genes)*100, digits=4)}%"),
    `Total Genes Removed` = num_genes_filtered,
    `Percent Genes Removed` = glue(
      "{format((num_genes_filtered/num_genes)*100, digits=4)}%")
  ))
}

##### Get filtered counts matrix #####
#' Function to get filtered counts matrix based on the specified filters
#' 
#' @param counts_data (tibble) normalized count data to apply filters to
#' @param var_threshold (double) percentile of variance to filter genes by
#' @param sample_threshold (double) number of non-zero samples to filter genes by
#' 
#' @return A tibble containing a processed normalized counts matrix filtered
#'          based on user specified parameters
#'
#' @examples
#' `get_filtered_counts(counts_data, 0.3, 20)`
get_filtered_counts <- function(counts_data, var_threshold, sample_threshold) {
  # Apply filters to the count data
  filtered_counts <- counts_data %>%
    # Pivot data to long format so that it can be filtered appropriately
    pivot_longer(cols = -c(gene), names_to = "sample", values_to = "counts") %>% 
    # Group by gene and calculate the variance and total samples with zero
    # for each gene
    group_by(gene) %>% 
    reframe(sample=sample, 
            counts = counts, 
            variance = var(counts), 
            zero_counts = ifelse(counts==0, 0,1), 
            total_counts = sum(zero_counts)) %>%
    # Calculate the percentile of variance
    mutate(percent_var = variance/median(variance)) %>%
    # Filter genes based on the user specified parameters
    filter(percent_var >= var_threshold, total_counts > sample_threshold) %>% 
    # Select only the important rows and return the counts data to wide format
    select(gene, sample, counts) %>%
    pivot_wider(names_from = "sample", values_from = "counts")
  return(filtered_counts)
}

##### Perform PCA #####
#' Function to perform PCA on the filtered counts data
#' 
#' @param filtered_count_data (tibble) Normalized counts matrix after filters 
#'                              have been applied
#'                            
#' @return The named list object returned from prcomp
#' 
#' @examples
#' `do_PCA(filtered_count_data)`
do_PCA <- function(filtered_count_data) {
  # Convert the filtered count data to a dataframe and store gene ids as the
  # row names
  filtered_count_df <- filtered_count_data %>% 
    column_to_rownames("gene") %>% as.data.frame()
  # Perform PCA on the filtered counts data
  pca_results <- prcomp(scale(t(filtered_count_df)), center=FALSE, scale=FALSE)
  
  return(pca_results)
}

##### PCA beeswarm plot #####
#' Function to generate beeswarm plots of inputted top N principal components
#' 
#' @param pca_results (list) output from the prcomp function
#' @param N (double) Number of principal components to include in plot
#' 
#' @return A ggplot object containing the PCA beeswarm plot
#' 
#' @examples
#' `PCA_plot(pca_results, 3)`
PCA_plot <- function(pca_results, N) {
  # Calculate the variance explained by each PC and store as a tibble
  pca_var_explained <- as_tibble(pca_results$sdev^2/sum(pca_results$sdev^2)) %>%
    rename(variance_explained = value)
  # Store the PCs as a tibble
  principal_components <- as_tibble(colnames(pca_results$x)) %>% 
    rename(PCs = value)
  # Create a tibble of PCs and respective variance explained by each component
  pca_var_tibble <- bind_cols(principal_components, pca_var_explained)
  # Merge the scores for the user-specified PCs and the variance explained for
  # each PC
  p <- as.data.frame(pca_results$x[,1:N]) %>% 
    pivot_longer(cols = everything(), names_to = "PCs", values_to = "scores") %>%
    merge(pca_var_tibble, by="PCs") %>%
    # Calculate the percentage of variance explained for each PC to include in
    # the generated plot
    mutate(PCs = glue(
      "{PCs} ({format(variance_explained*100, digits=4)}% of variance)")) %>%
    # Make a beeswarm plot of PCA results for the requested PCs
    ggplot(aes(x=PCs, y=scores, color=PCs)) + geom_beeswarm(cex=3) +
    theme(legend.position = "bottom", text = element_text(size = 18),
          axis.text.x = element_text(angle=90, hjust=1))
  
  return(p)
}

##### Diagnostic count data scatter plot #####
#' Function to generate diagnostic scatter plots for exploration of the counts data
#' 
#' @param count_data (tibble) Unfiltered normalized counts matrix
#' @param filtered_count_data (tibble) Normalized counts matrix after filters 
#'                            have been applied
#' @param plot_type (character) Specifies which diagnostic plot to output 
#'                    between `varianceScatter` and `sampleScatter`
#' 
#' @return A ggplot object
#'
#' @examples
#' `make_diag_plot(count_data, filtered_count_data, "varianceScatter")`
make_diag_plot <- function(count_data, filtered_count_data, filtered_genes, plot_type) {
  # Make a scatter plot of all genes in the uploaded counts matrix color
  # points based on whether or not the gene passed the user defined filters
  # or not
  full_counts <- count_data %>% 
    # pivot dataframe to long format
    pivot_longer(-c(gene), names_to = "sample", values_to = "counts") %>%
    # Add column to keep track of genes with 0 counts in each sample
    mutate(zero_counts = ifelse(counts==0, 1,0)) %>%
    # Calculate the overall statistics for each gene across all samples
    group_by(gene) %>% 
    summarize( 
      med_counts = median(counts),
      var_counts = var(counts), 
      total_counts = sum(zero_counts)) %>%
    # Add a column to distinguish genes which passed all the user-defined 
    # filters and those that did not
    mutate(filtered = ifelse(gene %in% filtered_genes,
                             "Pass","Not Pass"))
  # Output a plot based on the specified plot type
  if (plot_type == "varianceScatter") {
    p <- full_counts %>% 
      ggplot(aes(x=med_counts, y=var_counts, color=filtered)) +
      geom_point(alpha=0.8) + scale_x_continuous(trans = "log10") + 
      scale_y_continuous(trans = "log10") +
      scale_color_manual(values = c("#07a6ab","#5536d3"))
  }
  else {
    p <- full_counts %>% 
      ggplot(aes(x=med_counts, y=total_counts, color=filtered)) +
      geom_point(alpha=0.8) + scale_x_continuous(trans = "log10") + 
      scale_color_manual(values = c("#07a6ab","#5536d3"))
  }
  
  return(p)
}

##### Clustered Heatmap #####
#' Function to generate a clustered heatmap of genes and samples
#' 
#' @param filtered_count_data (tibble or matrix) Normalized counts matrix 
#'                              after filters have been applied
#'
#' @return A pheatmap object
#' 
#' @examples
#' `make_heatmap(filtered_count_data)`
make_heatmap <- function(filtered_count_data, gene_names=FALSE) {
  # Check if the counts data passed is a matrix or not
  if (!is.matrix(filtered_count_data)) {
    # Transform filtered count data to a matrix setting gene ids to rownames
    filtered_count_mat <- filtered_count_data %>% 
      column_to_rownames("gene") %>%
      as.matrix() %>% 
      # Perform a log10 transformation on all pseudo count values
      apply(2, function(x) log10(x+1))
  }
  else {
    filtered_count_mat <- filtered_count_data %>%
      # Perform a log10 transformation on all pseudo count values
      apply(2, function(x) log10(x+1))
  }
  
  # Make the heatmap scaling by gene counts
  p <- pheatmap(
    filtered_count_mat,
    scale="row",
    cluster_rows=TRUE,
    cluster_cols=TRUE,
    # cluster genes using pearson correlation
    clustering_distance_rows = "correlation",
    clustering_method_rows = "ward.D",
    # Cluster samples using euclidean distances
    clustering_distance_cols = "euclidean",
    clustering_method_cols = "ward.D",
    show_rownames=gene_names
  )
  return(p)
}

##### Load in differential expression results data #####
#' Function to read in a differntial expression results file
#' and store it in a tibble
#' 
#' @param filepath A text string to the full filepath of the file to load
#' 
#' @return A tibble containing results of differential expression
#'
#'@examples
#' `load_de_results(input)`
load_de_results <- function(input) {
  reactive({
    req(input$deResultsFile)
    out <- tryCatch(
      {
        # Sometimes the first column name isn't read in correctly so rename that
        # column if it doesn't get read in correctly
        de_results <- as_tibble(fread(input$deResultsFile$datapath)) %>%
          rename(gene = V1)
        return(de_results)
      },
      error=function(cond) {
        # If everything everything works appropriately just return in the dataframe
        re_results <- as_tibble(fread(input$deResultsFile$datapath)) 
        return(de_results)
      }
    )
    return(out)
  })
}

##### Differential Gene Expression w/DeSeq2 #####
#' Function to perform a DESeq2 analysis of rna seq data
#'
#' @param filtered_count_data (tibble ) Normalized counts matrix after filters 
#'                            have been applied.
#' @param meta_data (dataframe) The meta data for the counts matrix
#' @param count_filter (double) An arbitrary number of genes each row should 
#'                      contain or be excluded. DESeq2 suggests 10, but this 
#'                      could be customized while running. 
#'
#' @return A dataframe of DESeq results. It has a header describing the 
#'          condition, and 6 columns with genes as row names. 
#'          
#' @examples 
#' `run_deseq(filtered_count_data, coldata, 10)`
run_deseq <- function(filtered_count_data, meta_data, 
                      pheno_vector , count_filter) {
  # Find the column in the meta data that contains the two specified phenotypes
  pheno_column <- names(purrr::keep(meta_data, ~ any(. %in% pheno_vector)))
  # Subset the meta data by the GEO accession numbers and the column with the
  # phenotypes
  meta_data_subset <- meta_data %>% 
    select(`GEO_Accession (exp)`, all_of(pheno_column)) %>%
    rownames_to_column()
  # Pivot the counts matrix to long format to aid in later steps
  counts_long <- filtered_count_data %>% select(!gene) %>% 
    pivot_longer(everything(), names_to = "samples", values_to = "counts")
  # Create the design matrix using the specified phenotype column and the type
  # of library that was used 
  coldata <- counts_long %>%
    group_by(samples) %>% 
    rownames_to_column() %>% 
    ungroup() %>%
    merge(meta_data_subset, by="rowname") %>% 
    select(-c(rowname, `GEO_Accession (exp)`, counts)) %>% 
    # Ensure rows are in the same order as the columns of the count matrix
    arrange(samples, colnames(count_filter)[-1]) %>%
    column_to_rownames("samples") %>% 
    mutate(type=meta_data$LibraryLayout[1])
  # Set the column names for easier reference in DESeq2
  colnames(coldata) <- c("condition", "type")
  # Convert the tibble to a dataframe
  count_df <- filtered_count_data %>% column_to_rownames("gene")
  # Create a DESeq data set using the given data
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = as.matrix(round(count_df)),
    colData = coldata,
    rowData = row.names(count_df),
    design = ~ condition)
  # Filter reads before running DESeq2
  keep <- rowSums(BiocGenerics::counts(dds)) >= count_filter
  dds <- dds[keep,]
  
  # Run DESeq2
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  return(as.data.frame(res))
}

##### Volcano plot of Differential Expression Results #####
#' Function to make a volcano plot of differential expression results using
#' two columns specified by the user for the axes
#' 
#' @param de_results (Dataframe) of differential expression results
#' @param x_name The column name to plot on the x-axis
#' @param y_name The column name to plot on the y-axis
#' @param slider A negative integer value representing the magnitude of
#'                p-adjusted values to color
#' @param color1 One of the colors for the points.
#' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
#'
#' @return A ggplot object of a volcano plot.
#'
#' @examples 
#' `volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")`
volcano_plot <-
  function(de_results, x_name, y_name, slider, color1, color2) {
    de_results <- de_results %>% drop_na()
    group <- ifelse(de_results[[y_name]]<1*(10^slider),TRUE,FALSE)
    p <- ggplot(de_results, aes(x=!!sym(x_name), y=-log10(!!sym(y_name)), 
                           group=group, color=group)) + 
      geom_point() + 
      scale_color_manual(name=glue("padj < 1e{slider}"), 
                         values=c(color1, color2)) + 
      labs(x = x_name, y = glue("-log10({y_name})")) +
      theme(plot.title = element_blank(), legend.position = "bottom")
    return(p)
  }

#### Gene name conversion ####

#' Convert ensembl gene IDs into hgnc_symbol IDs using biomaRt. Inputs and 
#' outputs will likely not be the same size.
#'
#' @param gene_vector (vector) A vector of ensembl gene IDs or HGNC symbols
#' @param id_type (character) Denotes if gene_vector contains HGNC symbols or
#'                ensembl IDs
#'
#' @return A 2 column tibble that contains ensembl IDs in the first column,
#' and their corresponding HGNC symbol in the second column
#'
#' @examples 
#' `> ensembl_to_hgnc(tibble(c('202860_at', '1553551_s_at')))`
#' `affy_hg_u133_plus_2 hgnc_symbol`
#' `1        1553551_s_at      MT-ND1`
#' `2        1553551_s_at       MT-TI`
#' `3        1553551_s_at       MT-TM`
#' `4        1553551_s_at      MT-ND2`
#' `5           202860_at     DENND4B`
ensembl_to_hgnc <- function(gene_vector, id_type) {
  # Begin by creating a reference genome to to search through
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  if (id_type == "ensembl_id") {
    # Use BiomaRt to search for the corresponding hgnc gene symbols
    out <- getBM(
      attributes = c("hgnc_symbol", "ensembl_gene_id"),
      filters = ("ensembl_gene_id"),
      values = gene_vector,
      mart = ensembl
    )}
  else {
    out <- getBM(
      attributes = c("hgnc_symbol", "ensembl_gene_id"),
      filters = ("hgnc_symbol"),
      values = gene_vector,
      mart = ensembl
    )}
  return(as_tibble(out))
}

##### Match up gene symbols and ensembl IDs #####
#' Function to match up HGNC gene symbols and ensembl gene IDs based on a given
#' list of genes containing one or both and a counts matrix of ensembl IDs
#' 
#' @param filtered_count_data (tibble) Normalized counts matrix after filters
#'                            have been applied
#' @param gene_list (tibble) Single column tibble of user inputted genes
#' 
#' @return the normalized counts matrix with both ensembl ids and hgnc symbols 
#'          filtered to include only genes in the inputted gene list
#' 
#' @examples
#' `gene_match(filtered_count_data), c("A1CF", "APOE", "ENSG00000019995.6")`
gene_match <- function(filtered_count_data, gene_list) {
  # Remove version numbers from filtered_count_data if counts matrix uses
  # ensembl ids w/ version numbers
  filtered_count_data <- filtered_count_data %>% 
    separate(gene, c("gene"), sep = "\\.", extra = "drop")
  # Split the genes into two lists
  ensembl_list <- gene_list %>% 
    filter(str_detect(gene, "ENSG")) %>% 
    separate(gene, c("gene"), sep = "\\.", extra = "drop") %>% 
    pull(gene)
  hgnc_list <- gene_list %>%
    filter(!str_detect(gene, "ENSG")) %>% pull(gene)
  # If the gene list contains both ensembl IDs and HGNC symbols
  if (length(ensembl_list) != 0 & length(hgnc_list) != 0) {
    hgnc_tibble <- ensembl_to_hgnc(ensembl_list, "ensembl_id")
    ensembl_tibble <- ensembl_to_hgnc(hgnc_list, "hgnc_symbol") %>% 
      select(ensembl_gene_id, hgnc_symbol)
    full_gene_tibble <- bind_rows(hgnc_tibble, ensembl_tibble)
  }
  # If the list only contains ensembl IDs
  else if (length(ensembl_list) != 0) {
    full_gene_tibble <- ensembl_to_hgnc(ensembl_list, "ensembl_id")
  }
  # If the list only contains HGNC symbols
  else if (length(hgnc_list) != 0) {
    full_gene_tibble <- ensembl_to_hgnc(hgnc_list, "hgnc_symbol")
  } 
  # If the list is empty
  else {
    stop("Please input a gene")
  }
  # Merge the complete gene tibble with the counts matrix
  fil_gene_counts <- merge(filtered_count_data, full_gene_tibble, 
                           by.x="gene", by.y="ensembl_gene_id")
  if (length(ensembl_list) != 0 & !all(ensembl_list %in% fil_gene_counts$gene)) {
    stop(glue("Not all inputted genes were found {gene_list}"))
  }
  if (length(hgnc_list) != 0 & !all(hgnc_list %in% fil_gene_counts$hgnc_symbol)) {
    stop("Not all inputted genes were found")
  }
  return(fil_gene_counts)
}

##### Pairwise Gene Expression Correlation #####
#' Function to calculate the pairwise gene expression correlation for a set of
#' user inputted genes
#' 
#' @param fil_gene_counts (tibble) Output from gene_match()
#' @param threshold (double) Minimum correlation for drawing an edge between
#'                  two genes
#' @param output (character) specifies what output type should be expected
#' 
#' @return tibble with pairwise correlation stats for each gene
#' 
#' @examples
#' `pairwise_correlation(filtered_count_data, c("A1CF", "APOE"), 0.9)`
pairwise_correlation <- function(fil_gene_counts, threshold, output) {
  # Set hgnc symbols of the filtered genes to the rownames
  corr_dataframe <- fil_gene_counts %>% select(!gene) %>% 
    column_to_rownames("hgnc_symbol")
  # Make a clustered heatmap using the subsetted counts matrix
  if (output == "heatmap") {
    counts_mat <- corr_dataframe %>% as.matrix()
    p <- make_heatmap(counts_mat, TRUE)
    return(p)
  }
  else {
    # Create a graph adjacency based on correlation distances between genes in  
    # pairwise fashion using pearson distances
    g <- graph.adjacency(
      as.matrix(as.dist(cor(t(corr_dataframe), method="pearson"))),
      mode="undirected",
      weighted=TRUE,
      diag=FALSE
    )
    # Color negative correlation edges as blue
    E(g)[which(E(g)$weight<0)]$color <- "#00058e"
    # Color positive correlation edges as red
    E(g)[which(E(g)$weight>0)]$color <- "#b90b2c"
    # Convert edge weights to absolute values
    E(g)$weight <- abs(E(g)$weight)
    # Remove edges below absolute Pearson correlation specified by user
    g <- delete_edges(g, E(g)[which(E(g)$weight<threshold)])
    # Assign gene names to the graph vertices
    V(g)$name <- V(g)$name
    # Change shape of graph vertices
    #V(g)$shape <- "sphere"
    # Change color of graph vertices
    V(g)$color <- "#a283d6"
    # Change color of vertex frames
    V(g)$vertex.frame.color <- "white"
    # Scale the size of the vertices to be proportional to the counts of each gene 
    # represented by each vertexul and multiply scaled vales by a factor of 10
    scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
    vSizes <- (scale01(apply(corr_dataframe, 1, mean)) + 1.0) * 10
    # Amplify or decrease the width of the edges
    edgeweights <- E(g)$weight * 2.0
    # Convert the graph adjacency object into a minimum spanning tree based on 
    # Prim's algorithm
    mst <- mst(g, algorithm="prim")
    # Output the correlation graph
    if (output == "graph") {
      # Plot the correlation network
      p <- plot(
        mst,
        layout=layout.fruchterman.reingold,
        edge.curved=TRUE,
        vertex.size=vSizes,
        vertex.label.dist=-0.5,
        vertex.label.color="black",
        asp=FALSE,
        vertex.label.cex=2,
        edge.width=edgeweights,
        edge.arrow.mode=0,
        main="Pairwise Correlation Network")
      return(p)
    }
    # Output a table of the correlation metrics
    else {
      return(tibble(
        Gene = V(g)$name,
        Degree = degree(g),
        Closeness = closeness(g),
        Betweenness = betweenness(g)
      )) 
    }
  }
}
