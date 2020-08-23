#' Compute FCs and log2-FCs.
#'
#' \code{log2_FC} extends the results obtained via \code{\link{distinct_test}}, by computing fold changes (FC) and log2-FC between conditions.
#' 
#' @param res a \code{\linkS4class{data.frame}} with results as returned from \code{\link{distinct_test}}.
#' @param x a \code{linkS4class{SummarizedExperiment}} or a \code{linkS4class{SingleCellExperiment}} object.
#' @param name_assays_expression a character ("cpm" by default), 
#' indicating the name of the assays(x) element which stores the expression data (i.e., assays(x)$name_assays_expression).
#' We strongly encourage using normalized data, such as counts per million (CPM).
#' Do not use logarithm transformed data to compute FCs.
#' @param name_group a character ("group_id" by default), 
#' indicating the name of the colData(x) element which stores the group id of each cell (i.e., colData(x)$name_group).
#' @param name_cluster a character ("cluster_id" by default), 
#' indicating the name of the colData(x) element which stores the cluster id of each cell (i.e., colData(x)$name_cluster).
#' @return A \code{\linkS4class{data.frame}} object, extending the results in `res`.
#' Two additional columns are added: FC_group1/group2 and log2FC_group1/group2, inicating the FC and log2-FC of group1/group2.
#' A FC > 1 (or log2FC > 0) indicates up-regulation of group1 (compared to group2); while a FC < 1 (or log2FC < 0) indicates down-regulation of group1 (compared to group2).
#' @examples
<<<<<<< HEAD
#' # load pre-computed results (obtaines via `distinct_test`)
#' data("res", package = "distinct")
#' # load the input data:
#' data("Kang_subset", package = "distinct")
#' 
#' # We can optionally add the fold change (FC) and log2-FC between groups:
#' res = log2_FC(res = res,
#'   x = Kang_subset, 
#'   name_assays_expression = "cpm",
#'   name_group = "stim",
#'   name_cluster = "cell")
#' 
#' # Visualize significant results:
#' head(top_results(res))
=======
#' # Check the examples of 'distinct_test'
>>>>>>> upstream/RELEASE_3_11
#' 
#' @author Simone Tiberi \email{simone.tiberi@uzh.ch}
#' 
#' @seealso \code{\link{distinct_test}}, \code{\link{top_results}}
#' 
#' @export
log2_FC = function(res, 
<<<<<<< HEAD
                   x, 
                   name_assays_expression = "cpm",
                   name_group = "group_id",
                   name_cluster = "cluster_id"
){
=======
                       x, 
                       name_assays_expression = "cpm",
                       name_group = "group_id",
                       name_cluster = "cluster_id"
                       ){
>>>>>>> upstream/RELEASE_3_11
  stopifnot(
    is.data.frame(res),
    ( is(x, "SummarizedExperiment") | is(x, "SingleCellExperiment") ),
    is.character(name_assays_expression), length(name_assays_expression) == 1L,
    is.character(name_group), length(name_group) == 1L,
    is.character(name_cluster), length(name_cluster) == 1L
  )
  
  # check if FC/log2_FC are already present in res:
  sel_FC_cols = grep("FC", colnames(res))
  if( length(sel_FC_cols) > 0 ){
    message("'res' already contains columns 'FC' and/or 'log2FC': they will be overwritten")
    res = res[,-sel_FC_cols]
  }
  sel_mean_cols = grep("mean", colnames(res))
  if( length(sel_mean_cols) > 0 ){
    res = res[,-sel_mean_cols]
  }
  
<<<<<<< HEAD
  
  # count matrix:
=======
  # load count matrix:
>>>>>>> upstream/RELEASE_3_11
  sel = which(names(assays(x)) == name_assays_expression)
  if( length(sel) == 0 ){
    message("'name_assays_expression' not found in names(assays(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'name_assays_expression' found multiple times in names(assays(x))")
    return(NULL)
  }
  counts = assays(x)[[sel]]
<<<<<<< HEAD
  
  # cluster ids:
  sel = which(names(colData(x)) == name_cluster)
  if( length(sel) == 0 ){
    message("'name_cluster' not found in names(colData(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'name_cluster' found multiple times in names(colData(x))")
    return(NULL)
  }
  cluster_ids = colData(x)[[sel]]
  
  # group ids (1 per cell)
=======
  gene_ids = rownames(counts)
  
  # load group ids:
>>>>>>> upstream/RELEASE_3_11
  sel = which(names(colData(x)) == name_group)
  if( length(sel) == 0 ){
    message("'name_group' not found in names(colData(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'name_group' found multiple times in names(colData(x))")
    return(NULL)
  }
  group_ids = factor(colData(x)[[sel]])
<<<<<<< HEAD
  
  # TODO: GET GROUP
  group_levels = levels(group_ids)
  
  if(length(group_levels) != 2){
=======
  groups = unique(group_ids)
  
  # only for 2 groups:
  if(length(groups) != 2){
>>>>>>> upstream/RELEASE_3_11
    message("2 groups have to be provided in colData(x)$name_group")
    return(NULL)
  }
  
<<<<<<< HEAD
  group_1 = group_ids == group_levels[1]
  group_2 = group_ids == group_levels[2]
  
  pb_1 = assays(sumCountsAcrossCells(x = counts[,group_1],
                                     ids = cluster_ids[group_1],
                                     average = TRUE))$average
  # average computes the MEAN expression across cells (it divides by the number of cells in each cluster)
  
  pb_2 = assays(sumCountsAcrossCells(x = counts[,group_2],
                                     ids = cluster_ids[group_2],
                                     average = TRUE))$average
  
  # store cluster and gene names:
  n_genes = nrow(pb_1)
  
  cluster_levels = colnames(pb_1)
  
  # in each cluster, compute the mean for groups 1 and 2 (and FC-log2FC)
  FC_by_cluster = lapply(seq_along(cluster_levels), function(i){
    exp_1 = pb_1[,i]
    exp_2 = pb_2[,i]
    FC = exp_1/exp_2
    log2_FC = log2(FC)
    
    cbind( cluster_num = i, gene_num = seq_len(n_genes), exp_1, exp_2, FC, log2_FC)
  })
  FC_by_cluster = do.call(rbind,FC_by_cluster)
  
  # numeric transformations of res and cluster from "res":
  cluster_res_num = as.numeric(factor(res$cluster_id, levels =  cluster_levels ))
  gene_res_num = as.numeric(factor(res$gene, levels = rownames(pb_1)))
  
  # paste cluster and gene in res:
  cluster_gene_res = paste(cluster_res_num, gene_res_num, sep=".")
  
  # paste cluster and gene in FC_by_cluster
  cluster_gene_FC = paste(FC_by_cluster[,1], FC_by_cluster[,2], sep=".")
  
  matching = match(cluster_gene_res, cluster_gene_FC)
  tmp = FC_by_cluster[matching, 3:6 ]
  
  colnames(tmp) = c( paste0("mean_", group_levels[1]),
                     paste0("mean_", group_levels[2]),
                     paste0("FC_", group_levels[1], "/", group_levels[2]),
                     paste0("log2FC_", group_levels[1], "/", group_levels[2]) )
  rownames(tmp) = NULL
  
  res_final = cbind(res, tmp)
  
  message("FC and log2_FC computed, returning results")
  
=======
  group_1 = group_ids == groups[1]
  group_2 = group_ids == groups[2]
  
  # cluster ids:
  sel = which(names(colData(x)) == name_cluster)
  if( length(sel) == 0 ){
    message("'name_cluster' not found in names(colData(x))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'name_cluster' found multiple times in names(colData(x))")
    return(NULL)
  }
  cluster_ids = colData(x)[[sel]]
  
  # check that all clusters in res are present in cluster_id
  tmp = t(apply(res, 1, function(res_one_row){
    sel_cells = cluster_ids == res_one_row[[2]]
    sel_gene = gene_ids == res_one_row[[1]]
    
    exp_1 = mean(counts[sel_gene, group_1 & sel_cells])
    exp_2 = mean(counts[sel_gene, group_2 & sel_cells])
    FC = exp_1/exp_2
    log2_FC = log2(FC)
    
    c(exp_1, exp_2, FC, log2_FC)
  }))
  
  colnames(tmp) = c( paste0("mean_", groups[1]),
                     paste0("mean_", groups[2]),
                     paste0("FC_", groups[1], "/", groups[2]),
                     paste0("log2FC_", groups[1], "/", groups[2]) )
                     
  res_final = cbind(res, tmp)
  
>>>>>>> upstream/RELEASE_3_11
  res_final
}
