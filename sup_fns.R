##plotting functions for Sleuth
## G. Stein-O'Brien 
## 18 July 2021

#####################################

plot_bootstrap_pool<-function (obj, target_id, units = "est_counts", color_by = setdiff(colnames(obj$sample_to_covariates), 
                      "sample"), x_axis_angle = 50, divide_groups = TRUE) 
{
  stopifnot(sleuth:::check_norm_status(obj))
  units <- sleuth:::check_quant_mode(obj, units)
  df <- sleuth:::get_bootstrap_summary(obj, target_id, units)
  p <- ggplot(df, aes(x = color_by, y=mid))
  p <- p + geom_boxplot(aes_string(fill = color_by))
  p <- p + theme(legend.position = "none", axis.text.x = element_blank()) + xlab(color_by)
  p <- p + ggtitle(paste(sleuth_significant[sleuth_significant$target_id==target_id,"ext_gene"],":", target_id))
  p <- p + ylab(units)
  if (divide_groups) {
    p <- p + facet_wrap(color_by, strip.position = "bottom", scales = "free_x")
  }
  p
}


#####################################

plot_bootstrap_pool<-function (obj, target_id, units = "est_counts", color_by = setdiff(colnames(obj$sample_to_covariates), 
                                   "sample"), x_axis_angle = 50, divide_groups = TRUE, group_by=NA) 
{
  stopifnot(sleuth:::check_norm_status(obj))
  units <- sleuth:::check_quant_mode(obj, units)
  df <- sleuth:::get_bootstrap_summary(obj, target_id, units)
  p <- ggplot(df, aes(x = color_by, y=mid))
  p <- p + geom_boxplot(aes_string(fill = group_by))
  p <- p + theme(legend.position = "none", axis.text.x = element_blank()) + xlab(color_by)
  p <- p + ggtitle(paste(sleuth_significant[sleuth_significant$target_id==target_id,"ext_gene"],":", target_id))
  p <- p + ylab(units)
  if (divide_groups) {
    p <- p + facet_wrap(group_by, strip.position = "bottom", scales = "free_x")
  }
  p
}

#####################################




#####################################

plot_bootstrap_mod<-function (obj,target_id, target_id_source, units = "est_counts", color_by = setdiff(colnames(obj$sample_to_covariates), 
                              "sample"), x_axis_angle = 50, divide_groups = TRUE) 
{
  stopifnot(sleuth:::check_norm_status(obj))
  units <- sleuth:::check_quant_mode(obj, units)
  df <- sleuth:::get_bootstrap_summary(obj, target_id, units)
  p <- ggplot(df, aes(x = sample, ymin = min, lower = lower,  middle = mid, upper = upper, ymax = max))
  p <- p + geom_boxplot(stat = "identity", aes_string(fill = color_by))
  p <- p + theme(axis.text.x = element_text(angle = x_axis_angle, hjust = 1))
  p <- p + ggtitle(paste(target_id_source[target_id_source$target_id==target_id,"ext_gene"],":", target_id))
  p <- p + ylab(units)
  if (divide_groups) {
    p <- p + facet_wrap(color_by, strip.position = "bottom",  scales = "free_x")
  }
  p
}


#####################################

plot_transcript_heatmap_mod <-function (obj, transcripts_source, units = "tpm", trans = "log", cluster_transcripts = FALSE, 
          offset = 1, color_high = "#581845", color_mid = "#FFC300", color_low = "#DAF7A6", x_axis_angle = 50, 
          annotation_cols = setdiff(colnames(obj$sample_to_covariates), "sample"), ...) 
{
  transcripts<-transcripts_source$target_id
  stopifnot(sleuth:::check_norm_status(obj))
  units <- sleuth:::check_quant_mode(obj, units)
  if (!all(transcripts %in% obj$obs_norm$target_id)) {
    stop("Couldn't find the following transcripts: ", paste(transcripts[!(transcripts %in% 
                                                                            obj$obs_norm$target_id)], collapse = ", "), "\n\tIt is highly likely that some of them were filtered out.")
  }
  tabd_df <- obj$obs_norm[obj$obs_norm$target_id %in% transcripts, 
  ]
  if (units == "tpm") {
    tabd_df <- dplyr::select(tabd_df, target_id, sample, 
                             tpm)
    tabd_df <- reshape2::dcast(tabd_df, target_id ~ sample, 
                               value.var = "tpm")
  }
  else if (units == "est_counts") {
    tabd_df <- dplyr::select(tabd_df, target_id, sample, 
                             est_counts)
    tabd_df <- reshape2::dcast(tabd_df, target_id ~ sample, 
                               value.var = "est_counts")
  }
  else if (units == "scaled_reads_per_base") {
    tabd_df <- dplyr::select(tabd_df, target_id, sample, 
                             scaled_reads_per_base)
    tabd_df <- reshape2::dcast(tabd_df, target_id ~ sample, 
                               value.var = "scaled_reads_per_base")
  }
  else {
    stop("Didn't recognize the following unit: ", units)
  }
  rownames(tabd_df) <- tabd_df$target_id
  tabd_df$target_id <- NULL
  if (nchar(trans) > 0 && !is.null(trans)) {
    tFunc <- eval(parse(text = trans))
    trans_mat <- as.matrix(tFunc(tabd_df + offset))
  }
  else if (is.function(trans)) {
    trans_mat <- as.matrix(trans(tabd_df + offset))
  }
  else {
    trans_mat <- as.matrix(tabd_df)
  }
  s2c <- obj$sample_to_covariates
  if (is.null(annotation_cols)) {
    s2c <- NA
  }
  else if (!all(annotation_cols %in% colnames(s2c))) {
    bad_cols <- which(!(annotation_cols %in% colnames(s2c)))
    formatted_cols <- paste(annotation_cols[bad_cols], collapse = ", ")
    stop("At least one covariate selected in 'annotation_cols' does not exist.", 
         "\nHere are the covariates that do not exist: ", 
         formatted_cols)
  }
  else {
    rownames(s2c) <- s2c$sample
    s2c <- s2c[, annotation_cols, drop = FALSE]
  }
  colors <- colorRampPalette(c(color_low, color_mid, color_high))(100)
  pdf(file = NULL)
  genes<-sapply(rownames(trans_mat), function(x) transcripts_source[transcripts_source$target_id==x,"ext_gene"])
  if (cluster_transcripts) {
    p <- pheatmap::pheatmap(trans_mat, annotation_col = s2c, 
                            color = colors, cluster_cols = TRUE, cluster_rows = cluster_transcripts, 
                            labels_row = genes, ...)
  }
  else {
    p <- pheatmap::pheatmap(trans_mat, annotation_col = s2c, labels_row = genes,
                            color = colors, cluster_cols = TRUE, ...)
  }
  invisible(dev.off())
  p$gtable$grobs[[3]]$rot <- 360 - x_axis_angle
  gridExtra::grid.arrange(p$gtable)
  return(p)
}

#####################################
plot_sample_heatmap<- function (obj, use_filtered = TRUE, color_high = "white", color_low = "dodgerblue", 
          x_axis_angle = 50, annotation_cols = setdiff(colnames(obj$sample_to_covariates), 
                                                       "sample"), cluster_bool = TRUE, ...) 
{
  stopifnot(check_norm_status(obj))
  abund <- NULL
  if (use_filtered) {
    abund <- spread_abundance_by(obj$obs_norm_filt, "tpm", 
                                 obj$sample_to_covariates$sample)
  }
  else {
    abund <- spread_abundance_by(obj$obs_norm, "tpm", obj$sample_to_covariates$sample)
  }
  all_pairs <- apply_all_pairs(abund, jsd)
  s2c <- obj$sample_to_covariates
  if (is.null(annotation_cols)) {
    s2c <- NA
  }
  else if (!all(annotation_cols %in% colnames(s2c))) {
    bad_cols <- which(!(annotation_cols %in% colnames(s2c)))
    formatted_cols <- paste(annotation_cols[bad_cols], collapse = ", ")
    stop("At least one covariate selected in 'annotation_cols' does not exist.", 
         "\nHere are the covariates that do not exist: ", 
         formatted_cols)
  }
  else {
    rownames(s2c) <- s2c$sample
    s2c <- s2c[, annotation_cols, drop = FALSE]
  }
  colors <- colorRampPalette(c(color_high, color_low))(100)
  pdf(file = NULL)
  p <- pheatmap::pheatmap(all_pairs, annotation_col = s2c, 
                          color = colors, cluster_rows = cluster_bool, cluster_cols = cluster_bool, 
                          clustering_distance_cols = dist(all_pairs), clustering_distance_rows = dist(all_pairs), 
                          treeheight_row = 0, ...)
  invisible(dev.off())
  p$gtable$grobs[[3]]$rot <- 360 - x_axis_angle
  p$gtable$grobs[[4]]$label <- NULL
  gridExtra::grid.arrange(p$gtable)
}


#####################################


