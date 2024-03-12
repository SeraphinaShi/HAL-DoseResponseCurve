rename_df <- function(df){
  columns_to_rename <- grep("^oracal", names(df), value = TRUE)
  new_names <- sub("^oracal", "oracle", columns_to_rename)
  names(df)[names(df) %in% columns_to_rename] <- new_names
  
  return(df)
}

add_bound <- function(summary_df){
  
  bounds = c(0,1)
  
  summary_df$y_hat <- pmax(bounds[1], summary_df$y_hat)
  summary_df$y_hat <- pmin(summary_df$y_hat, bounds[2])
  
  summary_df$ci_lwr <- pmax(bounds[1], summary_df$ci_lwr)
  summary_df$ci_lwr <- pmin(summary_df$ci_lwr, bounds[2])
  
  summary_df$ci_upr = pmax(bounds[1], summary_df$ci_upr)
  summary_df$ci_upr <- pmin(summary_df$ci_upr, bounds[2])
  
  if("oracle_ci_lwr" %in% names(summary_df)){
    summary_df$oracle_ci_lwr <- pmax(bounds[1], summary_df$oracle_ci_lwr)
    summary_df$oracle_ci_lwr <- pmin(summary_df$oracle_ci_lwr, bounds[2])
    
    summary_df$oracle_ci_upr = pmax(bounds[1], summary_df$oracle_ci_upr)
    summary_df$oracle_ci_upr <- pmin(summary_df$oracle_ci_upr, bounds[2])
  }
  
  return(summary_df)
}



run_simu_rep_npcausal_summary <- function(results_list, method){
  
  result_all <-  do.call("rbind", results_list$all_results) %>% as.data.frame()
  result_all <- merge(as.data.frame(psi0_pnt), result_all, by=c("a"))
  
  result_summary <- result_all %>% 
    filter(SE != 0, ! is.na(y_hat)) %>% 
    mutate(SE = SE/sqrt(nn)) %>%
  mutate(bias = abs(y_hat - psi0),
         bias_se_ratio = bias / SE,
         cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr)) %>% 
    group_by(a) %>% 
    mutate(oracal_SE = sqrt(var(y_hat)),
           oracal_bias_se_ratio = bias / oracal_SE,
           oracal_ci_lwr = y_hat - 1.96 * oracal_SE,
           oracal_ci_upr = y_hat + 1.96 * oracal_SE,
           oracal_cover_rate = as.numeric(oracal_ci_lwr <= psi0 & psi0 <= oracal_ci_upr)) %>%
    summarise(across(where(is.numeric), mean)) %>% 
    ungroup() %>%
    mutate(method = method)
  
  results <- list(result_summary = result_summary,
                  all_results = results_list$all_results)
  
  return(results)
}

run_simu_rep_GAM_poly_summary <- function(results_list, method){
  
  result_all <-  do.call("rbind", results_list$all_results) %>% as.data.frame()
  result_all <- merge(as.data.frame(psi0_pnt), result_all, by=c("a"))
  
  result_summary <- result_all %>% 
    filter(SE != 0) %>% 
    mutate(bias = abs(y_hat - psi0),
           bias_se_ratio = bias / SE,
           cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr)) %>% 
    group_by(a) %>% 
    mutate(oracal_SE = sqrt(var(y_hat)),
           oracal_bias_se_ratio = bias / oracal_SE,
           oracal_ci_lwr = y_hat - 1.96 * oracal_SE,
           oracal_ci_upr = y_hat + 1.96 * oracal_SE,
           oracal_cover_rate = as.numeric(oracal_ci_lwr <= psi0 & psi0 <= oracal_ci_upr)) %>%
    summarise(across(where(is.numeric), mean)) %>% 
    ungroup() %>%
    mutate(method = method)
  
  
  results <- list(result_summary = result_summary,
                  all_results = results_list$all_results)
  
  
  return(results)
}


plot_performences_cv_ug_alla <- function(df, save_plot=NA, return_plot = "all"){
  
  color_cv =  "#2b6a99" 
  color_u_g = "#f16c23"

  a_max <- max(df$a)
  
  p_lambda <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=lambda_scaler, color=method), alpha = 0.7) +
    geom_point(aes(y=lambda_scaler, color=method), shape=17, size=2, alpha= 0.7) +
    ylim(0, 1.2) +
    labs(x="a", y="lambda_scaler", 
         title = "Lambda scaler",
         subtitle = paste0("upon CV_lambda ", round(mean(df$lambda[df$method == 'CV']), 6)))+
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Selector',
                       breaks=c('CV', 'Undersmooth'),
                       values=c('CV'=color_cv, 'Undersmooth'=color_u_g)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw()+
    theme(legend.position='none') 
  
  p_est_avg <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=psi0), alpha = 0.5, color="darkgrey") +
    geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr, color=method, fill=method, linetype='Delta'),  alpha=0.1) +
    geom_ribbon(aes(ymin=oracle_ci_lwr, ymax=oracle_ci_upr, color=method, fill=method, linetype = "Oracle"),  width=0.7, alpha=0.1) +
    geom_point(aes(y=psi0), color = "black") +
    geom_point(aes(y=y_hat, color=method), shape=17, size=2, alpha= 0.7) +
    labs(x="Treatment", y="Outcome", title = "(a) Estimations & 95% CIs") +
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Selector',
                       breaks=c('CV', 'Undersmooth'),
                       values=c('CV'=color_cv, 'Undersmooth'=color_u_g)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw() +
    theme(legend.box = "horizontal",
          legend.position='none')
  
  if (return_plot == "p_est_avg"){
    return(p_est_avg)
  }
  
  ymin_cr = max(0.95, min(df$cover_rate, df$oracle_cover_rate))
  p_cr <- ggplot(df, aes(x = a)) +  
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=ymin_cr,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
    geom_line(aes(y = cover_rate, color=method, linetype='Delta'), alpha=0.7) +
    geom_point(aes(y = cover_rate, color=method, linetype='Delta'), alpha=0.7) + 
    geom_line(aes(y = oracle_cover_rate, color=method, linetype='Oracle'), alpha=0.7) +
    geom_point(aes(y = oracle_cover_rate, color=method, linetype='Oracle'), alpha=0.7) + 
    labs(x="Treatment", y = "Coverage Rate", title="(b) 95% CI Coverage Rate") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(name='Selector',
                       breaks=c('CV', 'Undersmooth'),
                       values=c('CV'=color_cv, 'Undersmooth'=color_u_g)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw() +
    theme(legend.position='none') 
  
  if (return_plot == "p_cr"){
    return(p_cr)
  }
  
  df_bias_se <- df[! df$a %in% c(0,5), ]
  p_bias_se <- ggplot(df_bias_se, aes(x = a)) +  
    xlim(0,5) +
    geom_line(aes(y = bias_se_ratio, color=method, linetype='Delta'), alpha=0.7) +
    geom_point(aes(y = bias_se_ratio, color=method, linetype='Delta'), alpha=0.7) + 
    geom_line(aes(y = oracle_bias_se_ratio, color=method, linetype='Oracle'), alpha=0.7) +
    geom_point(aes(y = oracle_bias_se_ratio, color=method, linetype='Oracle'), alpha=0.7) + 
    labs(x="Treatment", y = "|Bias| / Standard Error", title="(c) Bias-SE Ratio") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    geom_hline(aes(yintercept=1/log(nn)), linetype = "dashed") +
    scale_color_manual(name='Selector',
                       breaks=c('CV', 'Undersmooth'),
                       values=c('CV'=color_cv, 'Undersmooth'=color_u_g)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw() +
    theme(legend.position='none') 
  if (return_plot == "p_bias_se"){
    return(p_bias_se)
  }
  
  p_mse <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = MSE, color=method, linetype='Oracle'),alpha=0.7) +
    geom_point(aes(y = MSE, color=method, linetype='Oracle'),alpha=0.7) + 
    labs(x="Treatment", y = "MSE", title="(f) Mean Squared Error") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Selector',
                       breaks=c('CV', 'Undersmooth'),
                       values=c('CV'=color_cv, 'Undersmooth'=color_u_g)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw() + 
    theme(legend.position='none')
  
  if (return_plot == "p_mse"){
    return(p_mse)
  }
  
  
  p_bias <- ggplot(df, aes(x = a, y = bias)) +  
    geom_line(aes(color=method)) +
    geom_point(aes(color=method)) + 
    labs(x="Treatment", y = "|Bias|", title="(d) Absolute Bias") +
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Selector',
                       breaks=c('CV', 'Undersmooth'),
                       values=c('CV'=color_cv, 'Undersmooth'=color_u_g)) +
    theme_bw() +
    theme(legend.position='none') 

  
  if (return_plot == "legend"){
    p_se <- ggplot(df, aes(x = a)) +  
      geom_line(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) +
      geom_point(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) + 
      geom_line(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) +
      geom_point(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) + 
      labs(x="Treatment", y = "SE", title="(e) Standard Error") +
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      scale_color_manual(name='Selector',
                         breaks=c('CV', 'Undersmooth'),
                         values=c('CV'=color_cv, 'Undersmooth'=color_u_g)) +
      scale_linetype_manual(name='Method',
                            breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() 
    
    legend <- get_legend(p_se)
    
    return(legend)
    }
  
  p_se <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) +
    geom_point(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) + 
    geom_line(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) +
    geom_point(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) + 
    labs(x="Treatment", y = "SE", title="(e) Standard Error") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Selector',
                       breaks=c('CV', 'Undersmooth'),
                       values=c('CV'=color_cv, 'Undersmooth'=color_u_g)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw() 
    # theme(legend.box = "horizontal")
  
  legend <- get_legend(p_se)
  p_se <- p_se + theme(legend.position='none')
  
  if (return_plot == "p_se"){
    return(p_se)
  }

  
  
  p <- grid.arrange(p_est_avg, p_cr, p_bias_se,  p_bias, p_se, p_mse, legend,
                    layout_matrix = rbind(c(1,1,2,2,3,3,NA),
                                          c(1,1,2,2,3,3,7),
                                          c(4,4,5,5,6,6,7),
                                          c(4,4,5,5,6,6,NA))# ,
                    # top = textGrob(paste0("HAL-based plug-in estimator performences"), gp=gpar(fontsize=17))
                    )

  if(!any(is.na(save_plot))){
    for (i in 1:length(save_plot)) {
      save_loc <- save_plot[i]
      ggsave(save_loc, plot=p, width = 10, height = 5, dpi = 800)
    }
  }
  
  return(p)
  
}


estimation_qqplot_cv_ug_alla <- function(results_list, save_plot=NA){
  df <- data.frame()
  for (method in c("CV", "U_G")){
    result_all <-  do.call("rbind", results_list[[method]]$all_results) %>% as.data.frame()
    result_all$method = method
    
    df <- rbind(df, result_all)
  }
  
  df <- df[!is.na(df$a), ]
  df <- df[df$a %in%  seq(from = 0.25, to = 4.75, by = 0.5),]
  
  p <- ggplot(df, aes(sample = y_hat)) + 
    stat_qq() + stat_qq_line() +
    labs(x = "Theoretical Quantiles",
         y = "Sample Quantiles",
         title = paste0("Q-Q Plot for Estimated Dose-Response")) +
    facet_grid(method ~ a, scales = 'free') +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          axis.text = element_text(size=10),
          axis.title = element_text(size=14),
          strip.text = element_text(size = 12))
  

  if(!any(is.na(save_plot))){
    for (i in 1:length(save_plot)) {
      save_loc <- save_plot[i]
      ggsave(save_loc, plot=p, width = 12, height = 4, dpi = 800)
    }
  }
  
  return(p)
}



plot_performences_adapt <- function(df, save_plot=NA){
  
  df$smooth_order = round(df$smooth_order, 4)
  df$smooth_order = factor(df$smooth_order)
  
  so_colors = c( "#2b6a99", "#A3A500", "#00B0F6", "#E76BF3", "#f16c23")
  if(sum(! unique(df$smooth_order) %in% 0:3) > 0) {
    mean_sl_pick_SO = unique(df$smooth_order)[! unique(df$smooth_order) %in% 0:3]
  } else {
    mean_sl_pick_SO = unique(df$smooth_order[df$sl_pick==1])
    so_colors[which(0:3 == mean_sl_pick_SO) ] = "#f16c23"
  }
  
  
  df$if_n_knots_default = round(df$if_n_knots_default, 4)
  df$if_n_knots_default = factor(df$if_n_knots_default)
  
  line_types = c('solid', 'dotted', 'twodash')
  if(sum(! unique(df$if_n_knots_default) %in% c(0,1)) > 0) {
    mean_sl_pick_if_n_knots_default = unique(df$if_n_knots_default)[! unique(df$if_n_knots_default) %in% c(0,1)]
  } else {
    mean_sl_pick_if_n_knots_default = unique(df$if_n_knots_default[df$sl_pick==1])
    line_types[which(c(0,1) == mean_sl_pick_SO) ] = "twodash"
  }
  
  
  a_max <- max(df$a)
  
  p_est_avg_e <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=psi0), alpha = 0.5, color="darkgrey") +
    geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr, color=smooth_order, fill=smooth_order, linetype=if_n_knots_default),  alpha=0.1) +
    geom_point(aes(y=psi0), color = "black") +
    geom_point(aes(y=y_hat, color=smooth_order), shape=17, size=2, alpha= 0.7) +
    labs(x="Treatment", y="Outcome", title = "(a) Estimations & 95% CIs") +
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    theme_bw() + 
    scale_color_manual(name='smooth order',
                       breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                       values = so_colors) +
    scale_fill_manual(name='smooth order',
                      breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                      values=so_colors) +
    scale_linetype_manual(name='default number of knots',
                          breaks=c("0", "1", as.character(mean_sl_pick_if_n_knots_default)),
                          values=line_types) +
    theme(legend.position='none')
  
  p_est_avg_o <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=psi0), alpha = 0.5, color="darkgrey") +
    geom_ribbon(aes(ymin=oracle_ci_lwr, ymax=oracle_ci_upr, color=smooth_order, fill=smooth_order, linetype=if_n_knots_default),  alpha=0.1) +
    geom_point(aes(y=psi0), color = "black") +
    geom_point(aes(y=y_hat, color=smooth_order), shape=17, size=2, alpha= 0.7) +
    labs(x="Treatment", y="Outcome", title = "(a) Estimations & 95% CIs") +                
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    theme_bw() +
    theme(legend.box = "horizontal",
          legend.position='none') + 
    scale_color_manual(name='smooth order',
                       breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                       values=so_colors) +
    scale_fill_manual(name='smooth order',
                      breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                      values=so_colors) +
    scale_linetype_manual(name='default number of knots',
                          breaks=c("0", "1", as.character(mean_sl_pick_if_n_knots_default)),
                          values=line_types) 
  ymin_cr_e = max(0.95, min(df$cover_rate))
  
  p_cr_e <- ggplot(df, aes(x = a)) +  
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=ymin_cr_e,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
    geom_line(aes(y = cover_rate, color=smooth_order, linetype=if_n_knots_default), alpha=0.7) +
    geom_point(aes(y = cover_rate, color=smooth_order, linetype=if_n_knots_default), alpha=0.7) + 
    labs(x="Treatment", y = "Coverage Rate", title="(b) 95% CI Coverage Rate") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(name='smooth order',
                       breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                       values=so_colors) +
    scale_linetype_manual(name='default number of knots',
                          breaks=c("0", "1", as.character(mean_sl_pick_if_n_knots_default)),
                          values=line_types)  +
    theme_bw() +
    theme(legend.position='none') 
  
  ymin_cr_o = max(0.95, min(df$oracle_cover_rate))
  p_cr_o <- ggplot(df, aes(x = a)) +  
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=ymin_cr_o,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
    geom_line(aes(y = oracle_cover_rate, color=smooth_order, linetype=if_n_knots_default), alpha=0.7) +
    geom_point(aes(y = oracle_cover_rate, color=smooth_order, linetype=if_n_knots_default), alpha=0.7) + 
    labs(x="Treatment", y = "Coverage Rate", title="(b) 95% CI Coverage Rate") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(name='smooth order',
                       breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                       values=so_colors) +
    scale_linetype_manual(name='default number of knots',
                          breaks=c("0", "1", as.character(mean_sl_pick_if_n_knots_default)),
                          values=line_types)  +
    theme_bw() +
    theme(legend.position='none') 
  
  p_mse_e <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = MSE, color=smooth_order, linetype=if_n_knots_default),alpha=0.7) +
    geom_point(aes(y = MSE, color=smooth_order, linetype=if_n_knots_default),alpha=0.7) + 
    # labs(title="Standard Error, Delta-method") +
    labs(x="Treatment", y = "MSE", title="(f) Mean Squared Error") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    theme_bw() + 
    guides(color = guide_legend(nrow = 3, byrow = F)) +
    scale_color_manual(name='smooth order',
                       breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                       values=so_colors) +
    scale_linetype_manual(name='default number of knots',
                          breaks=c("0", "1", as.character(mean_sl_pick_if_n_knots_default)),
                          values=line_types) + 
    theme(legend.position='none')
  
  p_mse_o <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = MSE, color=smooth_order, linetype=if_n_knots_default),alpha=0.7) +
    geom_point(aes(y = MSE, color=smooth_order, linetype=if_n_knots_default),alpha=0.7) +
    # labs(title="Standard Error, Oracle") +
    labs(x="Treatment", y = "MSE", title="(f) Mean Squared Error") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='smooth order',
                       breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                       values=so_colors) +
    scale_linetype_manual(name='default number of knots',
                          breaks=c("0", "1", as.character(mean_sl_pick_if_n_knots_default)),
                          values=line_types)  +
    theme_bw() + 
    theme(legend.position='none')
  
  p_bias <- ggplot(df, aes(x = a, y = bias, color=smooth_order, linetype=if_n_knots_default)) +  
    geom_line() +
    geom_point() + 
    labs(x="Treatment", y = "|Bias|", title="(d) Absolute Bias") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    theme_bw() +
    theme(legend.position='none') + 
    scale_color_manual(name='smooth order',
                       breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                       values=so_colors) +
    scale_linetype_manual(name='default number of knots',
                          breaks=c("0", "1", as.character(mean_sl_pick_if_n_knots_default)),
                          values=line_types)
  
  p_se_e <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = SE, color=smooth_order, linetype=if_n_knots_default),alpha=0.7) +
    geom_point(aes(y = SE, color=smooth_order, linetype=if_n_knots_default),alpha=0.7) + 
    # labs(title="Standard Error, Delta-method") +
    labs(x="Treatment", y = "SE", title="(e) Standard Error") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    theme_bw() + 
    theme(legend.box = "horizontal") +
    guides(color = guide_legend(nrow = 3, byrow = F)) +
    scale_color_manual(name='smooth order',
                       breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                       values=so_colors) +
    scale_linetype_manual(name='default number of knots',
                          breaks=c("0", "1", as.character(mean_sl_pick_if_n_knots_default)),
                          values=line_types)
  
  
  legend <- get_legend(p_se_e)
  p_se_e <- p_se_e + theme(legend.position='none')
  
  p_se_o <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = oracle_SE, color=smooth_order, linetype=if_n_knots_default),alpha=0.7) +
    geom_point(aes(y = oracle_SE, color=smooth_order, linetype=if_n_knots_default),alpha=0.7) +
    # 0labs(title="Standard Error, Oracle") +
    labs(x="Treatment", y = "SE", title="(e) Standard Error") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='smooth order',
                       breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                       values=so_colors) +
    scale_linetype_manual(name='default number of knots',
                          breaks=c("0", "1", as.character(mean_sl_pick_if_n_knots_default)),
                          values=line_types)  +
    theme_bw() + 
    theme(legend.position='none')
  
  df_bias_se <- df[! df$a %in% c(0,5), ]
  p_bias_se_e <- ggplot(df_bias_se, aes(x = a)) +  
    xlim(0,5) +
    geom_line(aes(y = bias_se_ratio, color=smooth_order, linetype=if_n_knots_default), alpha=0.7) +
    geom_point(aes(y = bias_se_ratio, color=smooth_order, linetype=if_n_knots_default), alpha=0.7) + 
    geom_hline(aes(yintercept=1/log(nn)), linetype = "dashed") +
    labs(x="Treatment", y = "|Bias| / Standard Error", title="(c) Bias-SE Ratio") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='smooth order',
                       breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                       values=so_colors) +
    scale_linetype_manual(name='default number of knots',
                          breaks=c("0", "1", as.character(mean_sl_pick_if_n_knots_default)),
                          values=line_types)  +
    theme_bw() +
    theme(legend.position='none') 

  p_bias_se_o <- ggplot(df_bias_se, aes(x = a)) +  
    xlim(0,5) +
    geom_line(aes(y = oracle_bias_se_ratio, color=smooth_order, linetype=if_n_knots_default), alpha=0.7) +
    geom_point(aes(y = oracle_bias_se_ratio, color=smooth_order, linetype=if_n_knots_default), alpha=0.7) + 
    geom_hline(aes(yintercept=1/log(nn)), linetype = "dashed") +
    labs(x="Treatment", y = "|Bias| / Standard Error", title="(c) Bias-SE Ratio") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='smooth order',
                       breaks=c("0", "1", "2", "3", as.character(mean_sl_pick_SO)),
                       values=so_colors) +
    scale_linetype_manual(name='default number of knots',
                          breaks=c("0", "1", as.character(mean_sl_pick_if_n_knots_default)),
                          values=line_types)  +
    theme_bw() +
    theme(legend.position='none') 
  
  p <- grid.arrange(p_est_avg_e, p_est_avg_o, 
                    p_cr_e, p_cr_o,
                    p_bias_se_e, p_bias_se_o,
                    p_bias, legend,
                    p_se_e, p_se_o, 
                    p_mse_o,
                    layout_matrix = rbind(c(1,1,2,2),
                                          c(1,1,2,2),
                                          c(3,3,4,4),
                                          c(3,3,4,4),
                                          c(5,5,6,6),
                                          c(5,5,6,6),
                                          c(7,7,8,8),
                                          c(7,7,8,8),
                                          c(9,9,10,10),
                                          c(9,9,10,10),
                                          c(NA,NA,11,11),
                                          c(NA,NA,11,11)),
                    top = textGrob(paste0("   Delta-Method                               Oracle    \n"), 
                                   gp=gpar(fontsize=17)))
  
  if(!any(is.na(save_plot))){
    for (i in 1:length(save_plot)) {
      save_loc <- save_plot[i]
      ggsave(save_loc, plot=p, width = 7, height = 14, dpi = 800)
    }
  }
  return(p)
  
}




plot_perforences_grid <- function(df, u_g_scaler=NA, save_plot=NA, max_bias_sd=NA){
  
  curve_pnts <- sort(c(seq(from = 0.25, to = 4.75, by = 0.5), 2, 4)) 
  
  df <- df[df$a %in%  curve_pnts,]
  
  legend_undersmoothing = ggplot(df) +  
    geom_line(aes(x = lambda_scaler, y = y_hat, color = "Undersmooth"), lty=2) + 
    geom_line(aes(x = lambda_scaler, y = y_hat, color = "CV"), lty=2) + 
    theme_bw() +
    scale_colour_manual(name="Selector",
                        values=c(Undersmooth= "#f16c23", CV="#2b6a99" )) #, Local="#619CFF",
  legend_undersmoothing <- get_legend(legend_undersmoothing)
  
  p_est_avg_list = list()
  p_bias_list = list()
  p_se_list = list()
  p_mse_list = list()
  p_bias_se_list = list()
  p_cr_list = list()
  # p_n_basis_list = list()
  legend = NA
  
  for (i in 1:length(curve_pnts)) {
    df_a <- df %>% filter(a == curve_pnts[i])
    
    p_est_avg = ggplot(df_a) +  
      geom_ribbon(aes(x = lambda_scaler, ymin=ci_lwr, ymax=ci_upr, color='Delta', fill = 'Delta'),  alpha=0.5) +
      geom_ribbon(aes(x = lambda_scaler, ymin=oracle_ci_lwr, ymax=oracle_ci_upr,  color='Oracle', fill = 'Oracle'),  alpha=0.5) +
      geom_line(aes(x = lambda_scaler, y = y_hat), color = "grey") + 
      geom_point(aes(x = lambda_scaler, y = y_hat)) + 
      geom_hline(aes(yintercept=psi0)) + 
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      scale_color_manual(breaks=c('Oracle', 'Delta'),
                         values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      scale_fill_manual(breaks=c('Oracle', 'Delta'),
                        values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      theme_bw() +
      labs(x = "", title = paste0('a = ', curve_pnts[i])) +
      theme(axis.title=element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.position='none')
    
    
    p_bias <- ggplot(df_a, aes(x = lambda_scaler, y = bias)) +  
      geom_line(color = "grey") + 
      geom_point() + 
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      theme_bw()+
      theme(axis.title=element_blank())
    
    p_se <- ggplot(df_a, aes(x = lambda_scaler)) +  
      geom_line(aes(y = SE, color='Delta')) + 
      geom_point(aes(y = SE, color='Delta')) + 
      geom_line(aes(y = oracle_SE, color='Oracle')) + 
      geom_point(aes(y = oracle_SE, color='Oracle')) +
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      scale_color_manual(name='Method',
                         breaks=c('Oracle', 'Delta'),
                         values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      theme_bw()+
      theme(axis.title=element_blank())
    
    legend <- get_legend(p_se)
    p_se <- p_se + theme(legend.position='none')
    
    p_mse <- ggplot(df_a, aes(x = lambda_scaler)) +  
      # geom_line(aes(y = MSE, color='Delta')) + 
      # geom_point(aes(y = MSE, color='Delta')) + 
      geom_line(aes(y = MSE, color='Oracle')) + 
      geom_point(aes(y = MSE, color='Oracle')) +
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      scale_color_manual(name='Method',
                         breaks=c('Oracle', 'Delta'),
                         values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      theme_bw()+
      theme(axis.title=element_blank()) +
      theme(legend.position='none')
      
    p_bias_se <- ggplot(df_a, aes(x = lambda_scaler)) + 
      geom_line(aes(y = bias_se_ratio, color = "Delta")) + 
      geom_point(aes(y = bias_se_ratio, color = "Delta")) +
      geom_line(aes(y = oracle_bias_se_ratio, color = "Oracle")) + 
      geom_point(aes(y = oracle_bias_se_ratio, color = "Oracle")) +
      geom_hline(aes(yintercept=1/log(nn)), linetype = "dashed") +
      scale_color_manual(name='Method',
                         breaks=c('Oracle', 'Delta'),
                         values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      theme_bw() +
      theme(axis.title=element_blank(),
            legend.position='none')
    
    if(!is.na(max_bias_sd)){
      p_bias_se <- p_bias_se + ylim(0,max_bias_sd)
    }
    
    p_cr <- ggplot(df_a, aes(x = lambda_scaler)) +  
      geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=0.95,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
      geom_line(aes(y = cover_rate, color = "Delta")) + 
      geom_point(aes(y = cover_rate, color = "Delta")) + 
      geom_line(aes(y = oracle_cover_rate, color = "Oracle")) + 
      geom_point(aes(y = oracle_cover_rate, color = "Oracle")) +
      scale_color_manual(name='Method',
                         breaks=c('Oracle', 'Delta'),
                         values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_bw()  + 
      theme(axis.title=element_blank(),
            legend.position='none')
    
    if(! is.na(u_g_scaler) ){
      p_est_avg <- p_est_avg +      
        geom_vline(xintercept = u_g_scaler, lty=2, col = "#f16c23") +
        geom_vline(xintercept = 1, lty=2, col = "#2b6a99") 
      
      p_bias <- p_bias +      
        geom_vline(xintercept = u_g_scaler, lty=2, col = "#f16c23") +
        geom_vline(xintercept = 1, lty=2, col = "#2b6a99") 
      
      p_se <- p_se +      
        geom_vline(xintercept = u_g_scaler, lty=2, col = "#f16c23") +
        geom_vline(xintercept = 1, lty=2, col = "#2b6a99") 
      
      p_mse <- p_mse +      
        geom_vline(xintercept = u_g_scaler, lty=2, col = "#f16c23") +
        geom_vline(xintercept = 1, lty=2, col = "#2b6a99") 
      
      p_bias_se <- p_bias_se +      
        geom_vline(xintercept = u_g_scaler, lty=2, col = "#f16c23") +
        geom_vline(xintercept = 1, lty=2, col = "#2b6a99") 
      
      p_cr <- p_cr +      
        geom_vline(xintercept = u_g_scaler, lty=2, col = "#f16c23") +
        geom_vline(xintercept = 1, lty=2, col = "#2b6a99") 
      
    }
    
    p_est_avg_list[[i]] = p_est_avg
    p_bias_list[[i]] = p_bias
    p_se_list[[i]] = p_se
    p_mse_list[[i]] = p_mse
    p_bias_se_list[[i]] = p_bias_se
    p_cr_list[[i]] = p_cr
  }
  
  
  g1 <- arrangeGrob(grobs = p_est_avg_list, nrow=1, left = grid::textGrob("Estimation, 95% CI", rot=90, gp=gpar(fontsize=12)))
  g2 <- arrangeGrob(grobs = p_bias_list, nrow=1, left = grid::textGrob("|Bias|", rot=90, gp=gpar(fontsize=12)))
  g3 <- arrangeGrob(grobs = p_se_list, nrow=1, left = grid::textGrob("Standard Error", rot=90, gp=gpar(fontsize=12)))
  g4 <- arrangeGrob(grobs = p_mse_list, nrow=1, left = grid::textGrob("MSE", rot=90, gp=gpar(fontsize=12)))
  g5 <- arrangeGrob(grobs = p_bias_se_list, nrow=1, left = grid::textGrob("Bias-SE Ratio", rot=90, gp=gpar(fontsize=12)))
  g6 <- arrangeGrob(grobs = p_cr_list, nrow=1, left = grid::textGrob("Coverage rate", rot=90, gp=gpar(fontsize=12)) )#,
  # bottom = grid::textGrob("lambda scalers", gp=gpar(fontsize=15)))
  
  
  
  p <- grid.arrange(g1, g6,  g5, g2, g3, g4,legend, legend_undersmoothing,
                    layout_matrix = rbind(c(1,NA),
                                          c(2,NA),
                                          c(3,8),
                                          c(4,7),
                                          c(5,NA),
                                          c(6,NA)),
                    widths=c(13, 1), 
                    top = textGrob("HAL-based plug-in estimator performances \n", 
                                   gp=gpar(fontsize=18)))  
  
  if(!any(is.na(save_plot))){
    for (i in 1:length(save_plot)) {
      save_loc <- save_plot[i]
      ggsave(save_loc, plot=p, width = 25, height = 9, dpi = 800)
    }
  }
  
  return(p)
}


results_grid_summary <- function(results_grid_in){
  
  lambda_scalers = c(1.2, 1.1, 10^seq(from=0, to=-3, length=20))
  
  results <- list()
  # no_empirical_CI_proportion <- c()
  
  for (i in 1:length(lambda_scalers)){
    
    lambda_scaler = lambda_scalers[i]
    
    all_results = results_grid_in[[i]]$all_results
    result_all <-  do.call("rbind", all_results) %>% as.data.frame()
    result_all <- merge(as.data.frame(psi0_pnt), result_all, by=c("a"))
    
    result_summary <- result_all %>%
      filter(SE != 0) %>%
      mutate(bias = abs(y_hat - psi0),
             bias_se_ratio = bias / SE,
             # bias_se_ratio_bt = bias / SE_bt,
             # cover_rate_bt = as.numeric(ci_lwr_bt <= psi0 & psi0 <= ci_upr_bt) ,
             cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr)) %>%
      group_by(a) %>%
      mutate(oracle_SE = sqrt(var(y_hat)),
             oracle_bias_se_ratio = bias / oracle_SE,
             oracle_ci_lwr = y_hat - 1.96 * oracle_SE,
             oracle_ci_upr = y_hat + 1.96 * oracle_SE,
             oracle_cover_rate = as.numeric(oracle_ci_lwr <= psi0 & psi0 <= oracle_ci_upr)) %>%
      summarise(across(where(is.numeric), mean)) %>%
      ungroup() %>%
      mutate(hal_fit_time_unit = 'secs',
             method = 'scale')
    
    
    results[[paste0("scale=", round(lambda_scaler, 4))]] <- list(result_summary = result_summary)
    
  }
  
  result_summary <- results[[1]]$result_summary
  for (i in 2:length(lambda_scalers)) {
    result_summary <- rbind(result_summary, results[[i]]$result_summary)
  }
  
  return(result_summary)
}


plot_compare_methods_performances <- function(df, save_plot=NA){
  
  color_0_hal = '#619CFF'
  color_u_adapt = "#85c876"
  color_gam = '#f7b722'
  color_poly = '#ef9db6'
  color_npcausal =  "blueviolet"
  
  df <- df[df$method %in% c("0_HAL", "U_SOadapt_HAL", "GAM", "Poly", "npcausal"), ]
  
  a_max <- max(df$a)
  
  #-------------------------------------------------
  ci_min = min(df$ci_lwr, df$oracle_ci_lwr)
  ci_max = max(df$ci_upr, df$oracle_ci_upr)
  
  p_est_avg <- list()
  for(i in 1:2){
    p_est_avg[[i]] <- ggplot(data=df, aes(x=a)) +
      geom_line(aes(y=psi0), alpha = 0.5, color="darkgrey") +
      geom_point(aes(y=psi0), color = "black") + 
      geom_point(aes(y=y_hat, color=method), shape=17, size=2, alpha= 0.7) +
      labs(x="a") +
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      scale_y_continuous(limits = c(ci_min, ci_max)) + 
      
      scale_color_manual(name='Method',
                         breaks=c('U_SOadapt_HAL', "0_HAL", 'GAM', 'Poly', 'npcausal'),
                         values=c('U_SOadapt_HAL'=color_u_adapt, "0_HAL"=color_0_hal, 'npcausal'=color_npcausal, 'GAM'=color_gam, 'Poly'=color_poly)) +
      scale_fill_manual(name='Method',
                        breaks=c('U_SOadapt_HAL', "0_HAL", 'GAM', 'Poly', 'npcausal'),
                        values=c('U_SOadapt_HAL'=color_u_adapt, "0_HAL"=color_0_hal, 'npcausal'=color_npcausal, 'GAM'=color_gam, 'Poly'=color_poly)) +
      scale_linetype_manual(breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() +
      theme(legend.box = "horizontal", legend.position='none') 
  }
  
  p_est_avg[[1]] <- p_est_avg[[1]] + geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr, color=method, fill=method, linetype='Delta'),  alpha=0.1) +
    labs(title = "(a.1) Estimations & 95% CIs [Delta]" , y="Outcome", x = "Treatment")
  p_est_avg[[2]] <- p_est_avg[[2]] + geom_ribbon(aes(ymin=oracle_ci_lwr, ymax=oracle_ci_upr, color=method, fill=method, linetype = "Oracle"),  width=0.7, alpha=0.1) +
    labs(title = "(a.2) Estimations & 95% CIs [Oracle]" , y="Outcome", x = "Treatment")
  
  #-------------------------------------------------
  p_bias <- ggplot(df, aes(x = a, y = bias)) +  
    geom_line(aes(color=method)) +
    geom_point(aes(color=method)) + 
    labs(title="(d) Absolute Bias", x="Treatment", y="|Bias|") +
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Method',
                       breaks=c('U_SOadapt_HAL', "0_HAL", 'GAM', 'Poly', 'npcausal'),
                       values=c('U_SOadapt_HAL'=color_u_adapt, "0_HAL"=color_0_hal, 'npcausal'=color_npcausal, 'GAM'=color_gam, 'Poly'=color_poly)) +
    theme_bw()
  
  legend <- get_legend(p_bias)
  p_bias <- p_bias + theme(legend.position='none')
  
  
  #-------------------------------------------------
  se_min = min(df$SE, df$oracle_SE)
  se_max = max(df$SE, df$oracle_SE)
  
  p_se <- list()
  
  p_se[[1]] <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) +
    geom_point(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) + 
    labs(title = "(e.1) SE [Delta]", y = "SE")
  p_se[[2]] <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) +
    geom_point(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) + 
    labs(title = "(e.2) SE [Oracle]", y = "")
  
  for(i in 1:2){
    p_se[[i]] <- p_se[[i]] +
      labs(x="Treatment") +
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      # scale_y_continuous(limits = c(se_min, se_max)) +
      scale_color_manual(name='Method',
                         breaks=c('U_SOadapt_HAL', "0_HAL", 'GAM', 'Poly', 'npcausal'),
                         values=c('U_SOadapt_HAL'=color_u_adapt, "0_HAL"=color_0_hal, 'npcausal'=color_npcausal, 'GAM'=color_gam, 'Poly'=color_poly)) +
      scale_linetype_manual(breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() + 
      theme(legend.box = "horizontal", legend.position='none')
  }
  
  #-------------------------------------------------
  mse_min = min(df$MSE, df$MSE)
  mse_max = max(df$MSE, df$MSE)
  
  p_mse <- list()
  
  p_mse[[1]] <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = MSE, color=method, linetype='Delta'),alpha=0.7) +
    geom_point(aes(y = MSE, color=method, linetype='Delta'),alpha=0.7) + 
    labs(title = "(f.1) MSE [Delta]", y = "MSE")
  p_mse[[2]] <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = MSE, color=method, linetype='Oracle'),alpha=0.7) +
    geom_point(aes(y = MSE, color=method, linetype='Oracle'),alpha=0.7) + 
    labs(title = "(f.2) MSE [Oracle]", y = "")
  
  for(i in 1:2){
    p_mse[[i]] <- p_mse[[i]] +
      labs(x="Treatment") +
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      # scale_y_continuous(limits = c(se_min, se_max)) +
      scale_color_manual(name='Method',
                         breaks=c('U_SOadapt_HAL', "0_HAL", 'GAM', 'Poly', 'npcausal'),
                         values=c('U_SOadapt_HAL'=color_u_adapt, "0_HAL"=color_0_hal, 'npcausal'=color_npcausal, 'GAM'=color_gam, 'Poly'=color_poly)) +
      scale_linetype_manual(breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() + 
      theme(legend.box = "horizontal", legend.position='none')
  }
  
  #-------------------------------------------------
  bias_se_min = min(df$bias_se_ratio, df$oracle_bias_se_ratio)
  bias_se_max = max(df$bias_se_ratio, df$oracle_bias_se_ratio)
  
  p_bias_sd <- list()
  
  df_bias_se <- df[! df$a %in% c(0,5), ]
  p_bias_sd[[1]] <- ggplot(df_bias_se, aes(x = a)) +  
    xlim(0,5) +
    geom_line(aes(y = bias_se_ratio, color=method, linetype='Delta'), alpha=0.7) +
    geom_point(aes(y = bias_se_ratio, color=method, linetype='Delta'), alpha=0.7) +
    geom_hline(aes(yintercept=1/log(nn)), linetype = "dashed") +
    labs(y = "|Bias| / Standard Error", title = "(c.1) Bias-SE Ratio [Delta]") 
  
  p_bias_sd[[2]] <- ggplot(df_bias_se, aes(x = a)) +  
    xlim(0,5) +
    geom_line(aes(y = oracle_bias_se_ratio, color=method, linetype='Oracle'), alpha=0.7) +
    geom_point(aes(y = oracle_bias_se_ratio, color=method, linetype='Oracle'), alpha=0.7) + 
    geom_hline(aes(yintercept=1/log(nn)), linetype = "dashed") +
    labs(y="", title = "(c.1) Bias-SE Ratio [Oracle]") 
  
  for(i in 1:2){
    p_bias_sd[[i]] <- p_bias_sd[[i]] +  
      labs(x='Treatment') +
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      scale_color_manual(name='Method',
                         breaks=c('U_SOadapt_HAL', "0_HAL", 'GAM', 'Poly', 'npcausal'),
                         values=c('U_SOadapt_HAL'=color_u_adapt, "0_HAL"=color_0_hal, 'npcausal'=color_npcausal, 'GAM'=color_gam, 'Poly'=color_poly)) +
      scale_linetype_manual(breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() +
      theme(legend.position='none') 
  }
  
  #-------------------------------------------------
  p_cr <- list()
  
  ymin_cr_e = max(0.95, min(df$cover_rate))
  p_cr[[1]] <- ggplot(df, aes(x = a)) +  
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=ymin_cr_e,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
    geom_line(aes(y = cover_rate, color=method, linetype='Delta'), alpha=0.7) +
    geom_point(aes(y = cover_rate, color=method, linetype='Delta'), alpha=0.7) + 
    labs(y="95% CI Coverage Rate", title = "(b.1) CI Coverage Rate [Delta]")
  
  ymin_cr_o = max(0.95, min(df$oracle_cover_rate))
  p_cr[[2]] <- ggplot(df, aes(x = a)) +  
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=ymin_cr_o,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
    geom_line(aes(y = oracle_cover_rate, color=method, linetype='Oracle'), alpha=0.7) +
    geom_point(aes(y = oracle_cover_rate, color=method, linetype='Oracle'), alpha=0.7) + 
    labs(y = "", title = "(b.1) CI Coverage Rate [Oracle]")
  
  for (i in 1:2) {
    p_cr[[i]] <- p_cr[[i]] +
      labs(x='Treatment')+
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_manual(name='Method',
                         breaks=c('U_SOadapt_HAL', "0_HAL", 'GAM', 'Poly', 'npcausal'),
                         values=c('U_SOadapt_HAL'=color_u_adapt, "0_HAL"=color_0_hal, 'npcausal'=color_npcausal, 'GAM'=color_gam, 'Poly'=color_poly)) +
      scale_linetype_manual(breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() +
      theme(legend.position='none') 
  }
  
  
  #-------------------------------------------------
  p <- grid.arrange(p_est_avg[[1]], p_est_avg[[2]],
                    p_cr[[1]], p_cr[[2]], 
                    p_bias_sd[[1]], p_bias_sd[[2]], 
                    p_bias, legend, 
                    p_se[[1]], p_se[[2]], 
                     p_mse[[2]],
                    layout_matrix = rbind(c(1,1,2,2),
                                          c(1,1,2,2),
                                          c(3,3,4,4),
                                          c(3,3,4,4),
                                          c(5,5,6,6),
                                          c(5,5,6,6),
                                          c(7,7,8,8),
                                          c(7,7,8,8),
                                          c(9,9,10,10),
                                          c(9,9,10,10),
                                          c(NA,NA,11,11),
                                          c(NA,NA,11,11)), 
                    top = textGrob(paste0("Compare Methods"), 
                                   gp=gpar(fontsize=17)))
  

  if(!any(is.na(save_plot))){
    for (i in 1:length(save_plot)) {
      save_loc <- save_plot[i]
      ggsave(save_loc, plot=p, width = 8, height = 15, dpi = 800)
      
      p_1 = p_est_avg[[1]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_est_delta.png"), plot=p_1, width = 4, height = 3, dpi = 800)
      
      p_2 = p_est_avg[[2]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_est_oracle.png"), plot=p_2, width = 4, height = 3, dpi = 800)
      
      p_bias = p_bias + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_bias.png"), plot=p_bias, width = 4, height = 3, dpi = 800)
      
      p_1 = p_se[[1]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_se_delta.png"), plot=p_1, width = 4, height = 3, dpi = 800)
      
      p_2 = p_se[[2]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_se_oracle.png"), plot=p_2, width = 4, height = 3, dpi = 800)
      
      p_1 = p_mse[[1]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_mse_delta.png"), plot=p_1, width = 4, height = 3, dpi = 800)
      
      p_2 = p_mse[[2]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_mse_oracle.png"), plot=p_2, width = 4, height = 3, dpi = 800)
      
      p_1 = p_bias_sd[[1]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_bias_sd_delta.png"), plot=p_1, width = 4, height = 3, dpi = 800)
      
      p_2 = p_bias_sd[[2]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_bias_sd_oracle.png"), plot=p_2, width = 4, height = 3, dpi = 800)
      
      p_1 = p_cr[[1]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_cr_delta.png"), plot=p_1, width = 4, height = 3, dpi = 800)
      
      p_2 = p_cr[[2]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_cr_oracle.png"), plot=p_2, width = 4, height = 3, dpi = 800)
      
    }
  }
  
  
  return(p)
  
}


plot_compare_methods_performances_npcausal <- function(df, save_plot=NA){
  
  color_u_adapt = "#85c876"
  color_npcausal =  "blueviolet"
  
  df <- df[df$method %in% c("U_SOadapt_HAL", "npcausal"), ]
  
  a_max <- max(df$a)
  
  #-------------------------------------------------
  ci_min = min(df$ci_lwr, df$oracle_ci_lwr)
  ci_max = max(df$ci_upr, df$oracle_ci_upr)
  
  p_est_avg <- list()
  for(i in 1:2){
    p_est_avg[[i]] <- ggplot(data=df, aes(x=a)) +
      geom_line(aes(y=psi0), alpha = 0.5, color="darkgrey") +
      geom_point(aes(y=psi0), color = "black") + 
      geom_point(aes(y=y_hat, color=method), shape=17, size=2, alpha= 0.7) +
      labs(x="a") +
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      scale_y_continuous(limits = c(ci_min, ci_max)) + 
      scale_color_manual(name='Method',
                         breaks=c('U_SOadapt_HAL', 'npcausal'),
                         values=c('U_SOadapt_HAL'=color_u_adapt, 'npcausal'=color_npcausal)) +
      scale_fill_manual(name='Method',
                        breaks=c('U_SOadapt_HAL', 'npcausal'),
                        values=c('U_SOadapt_HAL'=color_u_adapt,  'npcausal'=color_npcausal )) +
      scale_linetype_manual(breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() +
      theme(legend.box = "horizontal", legend.position='none') 
  }
  
  p_est_avg[[1]] <- p_est_avg[[1]] + geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr, color=method, fill=method, linetype='Delta'),  alpha=0.1) +
    labs(title = "(a.1) Estimations & 95% CIs [Delta]" , y="Outcome", x = "Treatment")
  p_est_avg[[2]] <- p_est_avg[[2]] + geom_ribbon(aes(ymin=oracle_ci_lwr, ymax=oracle_ci_upr, color=method, fill=method, linetype = "Oracle"),  width=0.7, alpha=0.1) +
    labs(title = "(a.2) Estimations & 95% CIs [Oracle]" , y="Outcome", x = "Treatment")
  
  #-------------------------------------------------
  p_bias <- ggplot(df, aes(x = a, y = bias)) +  
    geom_line(aes(color=method)) +
    geom_point(aes(color=method)) + 
    labs(title="(d) Absolute Bias", x="Treatment", y="|Bias|") +
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Method',
                       breaks=c('U_SOadapt_HAL', 'npcausal'),
                       values=c('U_SOadapt_HAL'=color_u_adapt, 'npcausal'=color_npcausal)) +
    theme_bw()
  
  legend <- get_legend(p_bias)
  p_bias <- p_bias + theme(legend.position='none')
  
  
  #-------------------------------------------------
  se_min = min(df$SE, df$oracle_SE)
  se_max = max(df$SE, df$oracle_SE)
  
  p_se <- list()
  
  p_se[[1]] <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) +
    geom_point(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) + 
    labs(title = "(e.1) SE [Delta]", y = "SE")
  p_se[[2]] <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) +
    geom_point(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) + 
    labs(title = "(e.2) SE [Oracle]", y = "")
  
  for(i in 1:2){
    p_se[[i]] <- p_se[[i]] +
      labs(x="Treatment") +
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      # scale_y_continuous(limits = c(se_min, se_max)) +
      scale_color_manual(name='Method',
                         breaks=c('U_SOadapt_HAL', 'npcausal'),
                         values=c('U_SOadapt_HAL'=color_u_adapt, 'npcausal'=color_npcausal)) +
      scale_linetype_manual(breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() + 
      theme(legend.box = "horizontal", legend.position='none')
  }
  
  #-------------------------------------------------
  mse_min = min(df$MSE, df$MSE)
  mse_max = max(df$MSE, df$MSE)
  
  p_mse <- list()
  
  p_mse[[1]] <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = MSE, color=method, linetype='Delta'),alpha=0.7) +
    geom_point(aes(y = MSE, color=method, linetype='Delta'),alpha=0.7) + 
    labs(title = "(f.1) MSE [Delta]", y = "MSE")
  p_mse[[2]] <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = MSE, color=method, linetype='Oracle'),alpha=0.7) +
    geom_point(aes(y = MSE, color=method, linetype='Oracle'),alpha=0.7) + 
    labs(title = "(f.2) MSE [Oracle]", y = "")
  
  for(i in 1:2){
    p_mse[[i]] <- p_mse[[i]] +
      labs(x="Treatment") +
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      # scale_y_continuous(limits = c(se_min, se_max)) +
      scale_color_manual(name='Method',
                         breaks=c('U_SOadapt_HAL', 'npcausal'),
                         values=c('U_SOadapt_HAL'=color_u_adapt, 'npcausal'=color_npcausal)) +
      scale_linetype_manual(breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() + 
      theme(legend.box = "horizontal", legend.position='none')
  }
  
  #-------------------------------------------------
  bias_se_min = min(df$bias_se_ratio, df$oracle_bias_se_ratio)
  bias_se_max = max(df$bias_se_ratio, df$oracle_bias_se_ratio)
  
  p_bias_sd <- list()
  
  df_bias_se <- df[! df$a %in% c(0,5), ]
    
  p_bias_sd[[1]] <- ggplot(df_bias_se, aes(x = a)) +  
    xlim(0,5) +
    geom_line(aes(y = bias_se_ratio, color=method, linetype='Delta'), alpha=0.7) +
    geom_point(aes(y = bias_se_ratio, color=method, linetype='Delta'), alpha=0.7) +
    geom_hline(aes(yintercept=1/log(nn)), linetype = "dashed") +
    labs(y = "|Bias| / Standard Error", title = "(c.1) Bias-SE Ratio [Delta]") 
  
  p_bias_sd[[2]] <- ggplot(df_bias_se, aes(x = a)) +  
    xlim(0,5) +
    geom_line(aes(y = oracle_bias_se_ratio, color=method, linetype='Oracle'), alpha=0.7) +
    geom_point(aes(y = oracle_bias_se_ratio, color=method, linetype='Oracle'), alpha=0.7) + 
    geom_hline(aes(yintercept=1/log(nn)), linetype = "dashed") +
    labs(y="", title = "(c.1) Bias-SE Ratio [Oracle]") 
  
  for(i in 1:2){
    p_bias_sd[[i]] <- p_bias_sd[[i]] +  
      labs(x='Treatment') +
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      scale_color_manual(name='Method',
                         breaks=c('U_SOadapt_HAL', 'npcausal'),
                         values=c('U_SOadapt_HAL'=color_u_adapt, 'npcausal'=color_npcausal)) +
      scale_linetype_manual(breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() +
      theme(legend.position='none') 
  }
  
  #-------------------------------------------------
  p_cr <- list()
  
  ymin_cr_e = max(0.95, min(df$cover_rate))
  p_cr[[1]] <- ggplot(df, aes(x = a)) +  
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=ymin_cr_e,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
    geom_line(aes(y = cover_rate, color=method, linetype='Delta'), alpha=0.7) +
    geom_point(aes(y = cover_rate, color=method, linetype='Delta'), alpha=0.7) + 
    labs(y="95% CI Coverage Rate", title = "(b.1) CI Coverage Rate [Delta]")
  
  ymin_cr_o = max(0.95, min(df$oracle_cover_rate))
  p_cr[[2]] <- ggplot(df, aes(x = a)) +  
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=ymin_cr_o,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
    geom_line(aes(y = oracle_cover_rate, color=method, linetype='Oracle'), alpha=0.7) +
    geom_point(aes(y = oracle_cover_rate, color=method, linetype='Oracle'), alpha=0.7) + 
    labs(y = "", title = "(b.1) CI Coverage Rate [Oracle]")
  
  for (i in 1:2) {
    p_cr[[i]] <- p_cr[[i]] +
      labs(x='Treatment')+
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_manual(name='Method',
                         breaks=c('U_SOadapt_HAL', 'npcausal'),
                         values=c('U_SOadapt_HAL'=color_u_adapt, 'npcausal'=color_npcausal)) +
      scale_linetype_manual(breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() +
      theme(legend.position='none') 
  }
  
  
  #-------------------------------------------------
  p <- grid.arrange(p_est_avg[[1]], p_est_avg[[2]],
                    p_cr[[1]], p_cr[[2]], 
                    p_bias_sd[[1]], p_bias_sd[[2]], 
                    p_bias, legend, 
                    p_se[[1]], p_se[[2]], 
                    p_mse[[2]],
                    layout_matrix = rbind(c(1,1,2,2),
                                          c(1,1,2,2),
                                          c(3,3,4,4),
                                          c(3,3,4,4),
                                          c(5,5,6,6),
                                          c(5,5,6,6),
                                          c(7,7,8,8),
                                          c(7,7,8,8),
                                          c(9,9,10,10),
                                          c(9,9,10,10),
                                          c(NA,NA,11,11),
                                          c(NA,NA,11,11)),
                    top = textGrob(paste0("Compare Methods"), 
                                   gp=gpar(fontsize=17)))
  

  if(!any(is.na(save_plot))){
    for (i in 1:length(save_plot)) {
      save_loc <- save_plot[i]
      ggsave(save_loc, plot=p, width = 8, height = 15, dpi = 800)
      
      p_1 = p_est_avg[[1]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_est_delta.png"), plot=p_1, width = 4, height = 3, dpi = 800)
      
      p_2 = p_est_avg[[2]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_est_oracle.png"), plot=p_2, width = 4, height = 3, dpi = 800)
      
      p_bias = p_bias + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_bias.png"), plot=p_bias, width = 4, height = 3, dpi = 800)
      
      p_1 = p_se[[1]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_se_delta.png"), plot=p_1, width = 4, height = 3, dpi = 800)
      
      p_2 = p_se[[2]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_se_oracle.png"), plot=p_2, width = 4, height = 3, dpi = 800)
      
      p_1 = p_mse[[1]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_mse_delta.png"), plot=p_1, width = 4, height = 3, dpi = 800)
      
      p_2 = p_mse[[2]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_mse_oracle.png"), plot=p_2, width = 4, height = 3, dpi = 800)
      
      p_1 = p_bias_sd[[1]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_bias_sd_delta.png"), plot=p_1, width = 4, height = 3, dpi = 800)
      
      p_2 = p_bias_sd[[2]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_bias_sd_oracle.png"), plot=p_2, width = 4, height = 3, dpi = 800)
      
      p_1 = p_cr[[1]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_cr_delta.png"), plot=p_1, width = 4, height = 3, dpi = 800)
      
      p_2 = p_cr[[2]] + theme(axis.title.y=element_blank(), title = element_blank())
      ggsave( paste0(gsub(".png", "", save_loc), "_cr_oracle.png"), plot=p_2, width = 4, height = 3, dpi = 800)
      
    }
  }
  
  return(p)
  
}


plot_perforences_grid_lambda <- function(df, u_g_lambda=NA, cv_lambda=NA, save_plot=NA, max_bias_sd=NA){
  
  
  curve_pnts <- sort(c(seq(from = 0.25, to = 4.75, by = 0.5), 2, 4)) 
  
  df <- df[df$a %in%  curve_pnts,]
  
  legend_undersmoothing = ggplot(df) +  
    geom_line(aes(x = lambda, y = y_hat, color = "Undersmooth"), lty=2) + 
    geom_line(aes(x = lambda, y = y_hat, color = "CV"), lty=2) + 
    theme_bw() +
    scale_colour_manual(name="Selector",
                        values=c(Undersmooth="#f16c23", CV="#2b6a99")) #, Local="#619CFF",
  legend_undersmoothing <- get_legend(legend_undersmoothing)
  
  p_est_avg_list = list()
  p_bias_list = list()
  p_se_list = list()
  p_mse_list = list()
  p_bias_se_list = list()
  p_cr_list = list()
  # p_n_basis_list = list()
  legend = NA
  
  for (i in 1:length(curve_pnts)) {
    df_a <- df %>% filter(a == curve_pnts[i])
    
    p_est_avg = ggplot(df_a) +  
      geom_ribbon(aes(x = lambda, ymin=ci_lwr, ymax=ci_upr, color='Delta', fill = 'Delta'),  alpha=0.5) +
      geom_ribbon(aes(x = lambda, ymin=oracle_ci_lwr, ymax=oracle_ci_upr,  color='Oracle', fill = 'Oracle'),  alpha=0.5) +
      geom_line(aes(x = lambda, y = y_hat), color = "grey") + 
      geom_point(aes(x = lambda, y = y_hat)) + 
      geom_hline(aes(yintercept=psi0)) + 
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      scale_color_manual(breaks=c('Oracle', 'Delta'),
                         values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      scale_fill_manual(breaks=c('Oracle', 'Delta'),
                        values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      theme_bw() +
      labs(x = "", title = paste0('a = ', curve_pnts[i])) +
      theme(axis.title=element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.position='none')
    
    
    p_bias <- ggplot(df_a, aes(x = lambda, y = bias)) +  
      geom_line(color = "grey") + 
      geom_point() + 
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      theme_bw()+
      theme(axis.title=element_blank())
    
    p_se <- ggplot(df_a, aes(x = lambda)) +  
      geom_line(aes(y = SE, color='Delta')) + 
      geom_point(aes(y = SE, color='Delta')) + 
      geom_line(aes(y = oracle_SE, color='Oracle')) + 
      geom_point(aes(y = oracle_SE, color='Oracle')) +
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      scale_color_manual(name='Method',
                         breaks=c('Oracle', 'Delta'),
                         values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      theme_bw()+
      theme(axis.title=element_blank())
    
    legend <- get_legend(p_se)
    p_se <- p_se + theme(legend.position='none')
    
    p_mse <- ggplot(df_a, aes(x = lambda)) +  
      geom_line(aes(y = MSE, color='Delta')) + 
      geom_point(aes(y = MSE, color='Delta')) + 
      # geom_line(aes(y = MSE, color='Oracle')) + 
      # geom_point(aes(y = MSE, color='Oracle')) +
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      scale_color_manual(name='Method',
                         breaks=c('Oracle', 'Delta'),
                         values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      theme_bw()+
      theme(axis.title=element_blank()) +
      theme(legend.position='none')
    

    p_bias_se <- ggplot(df_a, aes(x = lambda)) + 
      geom_line(aes(y = bias_se_ratio, color = "Delta")) + 
      geom_point(aes(y = bias_se_ratio, color = "Delta")) +
      geom_line(aes(y = oracle_bias_se_ratio, color = "Oracle")) + 
      geom_point(aes(y = oracle_bias_se_ratio, color = "Oracle")) +
      geom_hline(aes(yintercept=1/log(nn)), linetype = "dashed") +
      scale_color_manual(name='Method',
                         breaks=c('Oracle', 'Delta'),
                         values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      theme_bw() +
      theme(axis.title=element_blank(),
            legend.position='none')
    
    if(!is.na(max_bias_sd)){
      p_bias_se <- p_bias_se + ylim(0,max_bias_sd)
    }
    
    p_cr <- ggplot(df_a, aes(x = lambda)) +  
      geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=0.95,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
      geom_line(aes(y = cover_rate, color = "Delta")) + 
      geom_point(aes(y = cover_rate, color = "Delta")) + 
      geom_line(aes(y = oracle_cover_rate, color = "Oracle")) + 
      geom_point(aes(y = oracle_cover_rate, color = "Oracle")) +
      scale_color_manual(name='Method',
                         breaks=c('Oracle', 'Delta'),
                         values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
      # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_bw()  + 
      theme(axis.title=element_blank(),
            legend.position='none')
    
    if(! is.na(u_g_scaler) ){
      p_est_avg <- p_est_avg +      
        geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
        geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
      
      p_bias <- p_bias +      
        geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
        geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
      
      p_se <- p_se +      
        geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
        geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
      
      p_mse <- p_mse +      
        geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
        geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
      
      p_bias_se <- p_bias_se +      
        geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
        geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
      
      p_cr <- p_cr +      
        geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
        geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
      
    }
    
    p_est_avg_list[[i]] = p_est_avg
    p_bias_list[[i]] = p_bias
    p_se_list[[i]] = p_se
    p_mse_list[[i]] = p_mse
    p_bias_se_list[[i]] = p_bias_se
    p_cr_list[[i]] = p_cr
  }
  
  
  g1 <- arrangeGrob(grobs = p_est_avg_list, nrow=1, left = grid::textGrob("Estimation, 95% CI", rot=90, gp=gpar(fontsize=12)))
  g2 <- arrangeGrob(grobs = p_bias_list, nrow=1, left = grid::textGrob("|Bias|", rot=90, gp=gpar(fontsize=12)))
  g3 <- arrangeGrob(grobs = p_se_list, nrow=1, left = grid::textGrob("Standard Error", rot=90, gp=gpar(fontsize=12)))
  g4 <- arrangeGrob(grobs = p_mse_list, nrow=1, left = grid::textGrob("MSE", rot=90, gp=gpar(fontsize=12)))
  g5 <- arrangeGrob(grobs = p_bias_se_list, nrow=1, left = grid::textGrob("Bias-SE Ratio", rot=90, gp=gpar(fontsize=12)))
  g6 <- arrangeGrob(grobs = p_cr_list, nrow=1, left = grid::textGrob("Coverage rate", rot=90, gp=gpar(fontsize=12)) )#,
  # bottom = grid::textGrob("lambda scalers", gp=gpar(fontsize=15)))
  
  
  
  p <- grid.arrange(g1, g6, g5, g2, g3, g4, legend, legend_undersmoothing,
                    layout_matrix = rbind(c(1,NA),
                                          c(2,NA),
                                          c(3,8),
                                          c(4,7),
                                          c(5,NA),
                                          c(6,NA)),
                    widths=c(13, 1), 
                    top = textGrob("HAL-based plug-in estimator performances \n", 
                                   gp=gpar(fontsize=18)))  
  

  if(!any(is.na(save_plot))){
    for (i in 1:length(save_plot)) {
      save_loc <- save_plot[i]
      ggsave(save_loc, plot=p, width = 25, height = 9, dpi = 500)
    }
  }
  
  return(p)
}




plot_perforences_grid_lambda_a <- function(df, a=1, u_g_lambda=NA, cv_lambda=NA, save_plot=NA, max_bias_sd=NA){
  
  df_a = df[df$a == a, ]
  
  legend_undersmoothing = ggplot(df_a) +  
    geom_line(aes(x = lambda_scaler, y = y_hat, color = "Undersmooth"), lty=2) + 
    geom_line(aes(x = lambda_scaler, y = y_hat, color = "CV"), lty=2) + 
    theme_bw() +
    scale_colour_manual(name="Selector",
                        values=c(Undersmooth="#f16c23", CV="#2b6a99")) #, Local="#619CFF",
  legend_undersmoothing <- get_legend(legend_undersmoothing)
  
  
  p_est_avg = ggplot(df_a) +  
    geom_ribbon(aes(x = lambda, ymin=ci_lwr, ymax=ci_upr, color='Delta', fill = 'Delta'),  alpha=0.5) +
    geom_ribbon(aes(x = lambda, ymin=oracle_ci_lwr, ymax=oracle_ci_upr,  color='Oracle', fill = 'Oracle'),  alpha=0.5) +
    geom_line(aes(x = lambda, y = y_hat), color = "grey") + 
    geom_point(aes(x = lambda, y = y_hat)) + 
    geom_hline(aes(yintercept=psi0)) + 
    # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
    scale_color_manual(breaks=c('Oracle', 'Delta'),
                       values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
    scale_fill_manual(breaks=c('Oracle', 'Delta'),
                      values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
    theme_bw() +
    labs(x = "Penalty", y = "Outcome", title = "(a) Estimations & 95% CIs") +
    theme(axis.title=element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position='none')
  
  
  p_bias <- ggplot(df_a, aes(x = lambda, y = bias)) +  
    geom_line(color = "grey") + 
    geom_point() + 
    # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
    theme_bw()+
    theme(axis.title=element_blank()) + 
    labs(x="Penalty", y = "|Bias|", title="(d) Absolute Bias") 
  
  p_se <- ggplot(df_a, aes(x = lambda)) +  
    geom_line(aes(y = SE, color='Delta')) + 
    geom_point(aes(y = SE, color='Delta')) + 
    geom_line(aes(y = oracle_SE, color='Oracle')) + 
    geom_point(aes(y = oracle_SE, color='Oracle')) +
    # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
    scale_color_manual(name='Method',
                       breaks=c('Oracle', 'Delta'),
                       values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
    theme_bw()+
    theme(axis.title=element_blank()) + 
    labs(x="Penalty", y = "SE", title="(e) Standard Error") 
  
  legend <- get_legend(p_se)
  p_se <- p_se + theme(legend.position='none')
  
  p_mse <- ggplot(df_a, aes(x = lambda)) +  
    geom_line(aes(y = MSE, color='Delta')) + 
    geom_point(aes(y = MSE, color='Delta')) + 
    # geom_line(aes(y = MSE, color='Oracle')) + 
    # geom_point(aes(y = MSE, color='Oracle')) +
    # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
    scale_color_manual(name='Method',
                       breaks=c('Oracle', 'Delta'),
                       values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
    theme_bw()+
    theme(axis.title=element_blank()) +
    theme(legend.position='none') + 
    labs(x="Penalty", y = "MSE", title="(f) Mean Squared Error") 
  
  p_bias_se <- ggplot(df_a, aes(x = lambda)) + 
    geom_line(aes(y = bias_se_ratio, color = "Delta")) + 
    geom_point(aes(y = bias_se_ratio, color = "Delta")) +
    geom_line(aes(y = oracle_bias_se_ratio, color = "Oracle")) + 
    geom_point(aes(y = oracle_bias_se_ratio, color = "Oracle")) +
    geom_hline(aes(yintercept=1/log(nn)), linetype = "dashed") +
    scale_color_manual(name='Method',
                       breaks=c('Oracle', 'Delta'),
                       values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
    # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
    theme_bw() +
    theme(axis.title=element_blank(),
          legend.position='none') + 
    labs(x="Penalty", y = "|Bias| / Standard Error", title="(c) Bias-SE Ratio") 
  
  if(!is.na(max_bias_sd)){
    p_bias_se <- p_bias_se + ylim(0,max_bias_sd)
  }
  
  p_cr <- ggplot(df_a, aes(x = lambda)) +  
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=0.95,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
    geom_line(aes(y = cover_rate, color = "Delta")) + 
    geom_point(aes(y = cover_rate, color = "Delta")) + 
    geom_line(aes(y = oracle_cover_rate, color = "Oracle")) + 
    geom_point(aes(y = oracle_cover_rate, color = "Oracle")) +
    scale_color_manual(name='Method',
                       breaks=c('Oracle', 'Delta'),
                       values=c('Oracle'='#27B2AF', 'Delta'='#F4A7C1')) +
    # scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw()  + 
    theme(axis.title=element_blank(),
          legend.position='none') +
    labs(x="Penalty", y = "Coverage Rate", title="(b) 95% CI Coverage Rate") 
  
  if(! is.na(u_g_scaler) ){
    p_est_avg <- p_est_avg +      
      geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
      geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
    
    p_bias <- p_bias +      
      geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
      geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
    
    p_se <- p_se +      
      geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
      geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
    
    p_mse <- p_mse +      
      geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
      geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
    
    p_bias_se <- p_bias_se +      
      geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
      geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
    
    p_cr <- p_cr +      
      geom_vline(xintercept = u_g_lambda, lty=2, col = "#f16c23") +
      geom_vline(xintercept = cv_lambda, lty=2, col = "#2b6a99") 
    
  }
  
  
  p <- grid.arrange(p_est_avg, p_cr,  p_bias_se,   
                    p_bias, p_se, p_mse, 
                    legend_undersmoothing, legend,
                    layout_matrix = rbind(c(1,1,2,2,3,3,NA),
                                          c(1,1,2,2,3,3,7),
                                          c(4,4,5,5,6,6,8),
                                          c(4,4,5,5,6,6,NA))
  )
  

  if(!any(is.na(save_plot))){
    for (i in 1:length(save_plot)) {
      save_loc <- save_plot[i]
      ggsave(save_loc, plot=p, width = 10, height = 5, dpi = 580)
    }
  }
  
  return(p)
}



plot_performences_knots <- function(df, save_plot=NA, return_plot = "all"){
  
  color_cv =  "#2b6a99" 
  color_u_g = "#f16c23"
  color_cv_10 =  "purple2"
  color_u_g_10 = "#1b7c3d"
  
  a_max <- max(df$a)
  
  p_lambda <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=lambda_scaler, color=method), alpha = 0.7) +
    geom_point(aes(y=lambda_scaler, color=method), shape=17, size=2, alpha= 0.7) +
    ylim(0, 1.2) +
    labs(x="a", y="lambda_scaler", 
         title = "Lambda scaler",
         subtitle = paste0("upon CV_lambda ", round(mean(df$lambda[df$method == 'CV']), 6)))+
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Selector',
                       breaks=c('CV_knots20', 'Undersmooth_knots20','CV_knots10', 'Undersmooth_knots10'),
                       values=c('CV_knots20'=color_cv, 'Undersmooth_knots20'=color_u_g,'CV_knots10'=color_cv_10, 'Undersmooth_knots10'=color_u_g_10)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw()+
    theme(legend.position='none') 
  
  p_est_avg <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=psi0), alpha = 0.5, color="darkgrey") +
    geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr, color=method, fill=method, linetype='Delta'),  alpha=0.1) +
    geom_ribbon(aes(ymin=oracle_ci_lwr, ymax=oracle_ci_upr, color=method, fill=method, linetype = "Oracle"),  width=0.7, alpha=0.1) +
    geom_point(aes(y=psi0), color = "black") +
    geom_point(aes(y=y_hat, color=method), shape=17, size=2, alpha= 0.7) +
    labs(x="Treatment", y="Outcome", title = "(a) Estimations & 95% CIs") +
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Selector',
                       breaks=c('CV_knots20', 'Undersmooth_knots20','CV_knots10', 'Undersmooth_knots10'),
                       values=c('CV_knots20'=color_cv, 'Undersmooth_knots20'=color_u_g,'CV_knots10'=color_cv_10, 'Undersmooth_knots10'=color_u_g_10)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw() +
    theme(legend.box = "horizontal",
          legend.position='none')
  
  if (return_plot == "p_est_avg"){
    return(p_est_avg)
  }
  
  ymin_cr = max(0.95, min(df$cover_rate, df$oracle_cover_rate))
  p_cr <- ggplot(df, aes(x = a)) +  
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=ymin_cr,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
    geom_line(aes(y = cover_rate, color=method, linetype='Delta'), alpha=0.7) +
    geom_point(aes(y = cover_rate, color=method, linetype='Delta'), alpha=0.7) + 
    geom_line(aes(y = oracle_cover_rate, color=method, linetype='Oracle'), alpha=0.7) +
    geom_point(aes(y = oracle_cover_rate, color=method, linetype='Oracle'), alpha=0.7) + 
    labs(x="Treatment", y = "Coverage Rate", title="(b) 95% CI Coverage Rate") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(name='Selector',
                       breaks=c('CV_knots20', 'Undersmooth_knots20','CV_knots10', 'Undersmooth_knots10'),
                       values=c('CV_knots20'=color_cv, 'Undersmooth_knots20'=color_u_g,'CV_knots10'=color_cv_10, 'Undersmooth_knots10'=color_u_g_10)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw() +
    theme(legend.position='none') 
  
  if (return_plot == "p_cr"){
    return(p_cr)
  }
  
  df_bias_se <- df[! df$a %in% c(0,5), ]
  p_bias_se <- ggplot(df_bias_se, aes(x = a)) +  
    xlim(0,5) +
    geom_line(aes(y = bias_se_ratio, color=method, linetype='Delta'), alpha=0.7) +
    geom_point(aes(y = bias_se_ratio, color=method, linetype='Delta'), alpha=0.7) + 
    geom_line(aes(y = oracle_bias_se_ratio, color=method, linetype='Oracle'), alpha=0.7) +
    geom_point(aes(y = oracle_bias_se_ratio, color=method, linetype='Oracle'), alpha=0.7) + 
    labs(x="Treatment", y = "|Bias| / Standard Error", title="(c) Bias-SE Ratio") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    geom_hline(aes(yintercept=1/log(nn)), linetype = "dashed") +
    scale_color_manual(name='Selector',
                       breaks=c('CV_knots20', 'Undersmooth_knots20','CV_knots10', 'Undersmooth_knots10'),
                       values=c('CV_knots20'=color_cv, 'Undersmooth_knots20'=color_u_g,'CV_knots10'=color_cv_10, 'Undersmooth_knots10'=color_u_g_10)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw() +
    theme(legend.position='none') 
  if (return_plot == "p_bias_se"){
    return(p_bias_se)
  }
  
  p_mse <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = MSE, color=method, linetype='Delta'),alpha=0.7) +
    geom_point(aes(y = MSE, color=method, linetype='Delta'),alpha=0.7) + 
    # geom_line(aes(y = MSE, color=method, linetype='Oracle'),alpha=0.7) +
    # geom_point(aes(y = MSE, color=method, linetype='Oracle'),alpha=0.7) + 
    labs(x="Treatment", y = "MSE", title="(f) Mean Squared Error") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Selector',
                       breaks=c('CV_knots20', 'Undersmooth_knots20','CV_knots10', 'Undersmooth_knots10'),
                       values=c('CV_knots20'=color_cv, 'Undersmooth_knots20'=color_u_g,'CV_knots10'=color_cv_10, 'Undersmooth_knots10'=color_u_g_10)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw() + 
    theme(legend.position='none')
  
  if (return_plot == "p_mse"){
    return(p_mse)
  }
  
  
  p_bias <- ggplot(df, aes(x = a, y = bias)) +  
    geom_line(aes(color=method)) +
    geom_point(aes(color=method)) + 
    labs(x="Treatment", y = "|Bias|", title="(d) Absolute Bias") +
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Selector',
                       breaks=c('CV_knots20', 'Undersmooth_knots20','CV_knots10', 'Undersmooth_knots10'),
                       values=c('CV_knots20'=color_cv, 'Undersmooth_knots20'=color_u_g,'CV_knots10'=color_cv_10, 'Undersmooth_knots10'=color_u_g_10)) +
    theme_bw() +
    theme(legend.position='none') 
  
  
  if (return_plot == "legend"){
    p_se <- ggplot(df, aes(x = a)) +  
      geom_line(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) +
      geom_point(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) + 
      geom_line(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) +
      geom_point(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) + 
      labs(x="Treatment", y = "SE", title="(e) Standard Error") +
      scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
      scale_color_manual(name='Selector',
                         breaks=c('CV_knots20', 'Undersmooth_knots20','CV_knots10', 'Undersmooth_knots10'),
                         values=c('CV_knots20'=color_cv, 'Undersmooth_knots20'=color_u_g,'CV_knots10'=color_cv_10, 'Undersmooth_knots10'=color_u_g_10)) +
      scale_linetype_manual(name='Method',
                            breaks=c('Oracle', 'Delta'),
                            values=c('Oracle'=1, 'Delta'=5)) +
      theme_bw() 
    
    legend <- get_legend(p_se)
    
    return(legend)
  }
  
  p_se <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) +
    geom_point(aes(y = SE, color=method, linetype='Delta'),alpha=0.7) + 
    geom_line(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) +
    geom_point(aes(y = oracle_SE, color=method, linetype='Oracle'),alpha=0.7) + 
    labs(x="Treatment", y = "SE", title="(e) Standard Error") + 
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Selector',
                       breaks=c('CV_knots20', 'Undersmooth_knots20','CV_knots10', 'Undersmooth_knots10'),
                       values=c('CV_knots20'=color_cv, 'Undersmooth_knots20'=color_u_g,'CV_knots10'=color_cv_10, 'Undersmooth_knots10'=color_u_g_10)) +
    scale_linetype_manual(name='Method',
                          breaks=c('Oracle', 'Delta'),
                          values=c('Oracle'=1, 'Delta'=5)) +
    theme_bw() 
  # theme(legend.box = "horizontal")
  
  legend <- get_legend(p_se)
  p_se <- p_se + theme(legend.position='none')
  
  if (return_plot == "p_se"){
    return(p_se)
  }
  
  
  
  p <- grid.arrange(p_est_avg, p_cr, p_bias_se,  p_bias, p_se, p_mse, legend,
                    layout_matrix = rbind(c(1,1,2,2,3,3,NA),
                                          c(1,1,2,2,3,3,7),
                                          c(4,4,5,5,6,6,7),
                                          c(4,4,5,5,6,6,NA))# ,
                    # top = textGrob(paste0("HAL-based plug-in estimator performences"), gp=gpar(fontsize=17))
  )
  
  if(!any(is.na(save_plot))){
    for (i in 1:length(save_plot)) {
      save_loc <- save_plot[i]
      ggsave(save_loc, plot=p, width = 10, height = 5, dpi = 800)
    }
  }
  
  return(p)
  
}

