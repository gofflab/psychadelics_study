##plotting functions for Slueth

plot_bootstrap_mod<-function (obj, target_id, units = "est_counts", color_by = setdiff(colnames(obj$sample_to_covariates), 
                      "sample"), x_axis_angle = 50, divide_groups = TRUE) 
{
  stopifnot(sleuth:::check_norm_status(obj))
  units <- sleuth:::check_quant_mode(obj, units)
  df <- sleuth:::get_bootstrap_summary(obj, target_id, units)
  p <- ggplot(df, aes(x = sample, ymin = min, lower = lower, 
                      middle = mid, upper = upper, ymax = max))
  p <- p + geom_boxplot(stat = "identity", aes_string(fill = color_by))
  p <- p + theme(axis.text.x = element_text(angle = x_axis_angle, 
                                            hjust = 1))
  p <- p + ggtitle(paste(sleuth_significant[sleuth_significant$target_id==target_id,"ext_gene"],":", target_id))
  p <- p + ylab(units)
  if (divide_groups) {
    p <- p + facet_wrap(color_by, strip.position = "bottom", 
                        scales = "free_x")
  }
  p
}