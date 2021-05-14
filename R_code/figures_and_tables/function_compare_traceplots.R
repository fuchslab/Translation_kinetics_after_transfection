#' Traceplots for two stanfit objects with same parameters next to each othe
#'
#' @param stanfit1 stanfit object for which plots will be shown in the lower and the diagonal panels
#' @param stanfit2 stanfit object for which plots will be shown in the upper and the diagonal panels
#' @param pars character vector of parameter names that should be contained in both stanfit objects stanfit1 and stanfit2
#' @param inc_warmup
#' @param true_values vector of the true parameter values with which the data was generated. default = NA
#' @param pars_names character vector or vector of expressions giving the parameter names as to be printed in the plot
#'
#' @return 
#'
compare_traceplots <- function(stanfit1, stanfit2, pars, pars_names,
                               inc_warmup = FALSE, common_scale = FALSE){
  if(common_scale){ # in case the  traceplots are to have the same y axis for 
    # both samples, determine the range of values for each parameter
    sample1 <- as.matrix(stanfit1, pars = pars)
    sample2 <- as.matrix(stanfit2, pars = pars)
    y_lims <- apply(rbind(sample1, sample2), 2, range)
  }
  
  plots <- list()
  n_pars <- length(pars)
  for(i in 1:n_pars){
    p1 <- rstan::traceplot(stanfit1, pars = pars[i], inc_warmup = inc_warmup) +
      theme(legend.position = "none") +
      labs(y = pars_names[i]) + 
      theme(
        axis.title.y = element_text(size=12, face="bold", angle = 90, vjust = 0.5),
        plot.margin = unit(c(0,0,0,0), "lines")
      ) 
    p2 <- rstan::traceplot(stanfit2, pars = pars[i], inc_warmup = inc_warmup) +
      theme(legend.position = "none") +
      labs(y = pars_names[i]) + 
      theme(
        axis.title.y = element_text(size=12, face="bold", angle = 90, vjust = 0.5),
        plot.margin = unit(c(0,0.6,0,0.5), "lines")
      )  
    if(common_scale){
      p1 <- p1 + ylim(y_lims[1, pars[i]], y_lims[2, pars[i]])
      p2 <- p2 + ylim(y_lims[1, pars[i]], y_lims[2, pars[i]])
    }
    if(i != n_pars){
      p1 <- p1 + theme(axis.text.x=element_blank())
      p2 <- p2 + theme(axis.text.x=element_blank())
    }
    
    plots[[2 * (i-1) + 1]] <- p1
    plots[[2 * i]] <- p2
  }
  
  title1 <- cowplot::ggdraw() + 
    cowplot::draw_label("SDE model", fontface = 'bold', x = .57, size = 14) +
    theme(# add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  title2 <- cowplot::ggdraw() + 
    cowplot::draw_label("ODE model", fontface = 'bold', x = .57, size = 14) +
    theme(# add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0)
    )
  
  titles <- cowplot::plot_grid(title1, title2, ncol = 2)
  
  # arrange plots in grd with two columns
  plotgrid <- cowplot::plot_grid(plotlist = plots, ncol = 2, align = 'v', 
                                 rel_heights = c(rep(1, times = n_pars - 1), 1.18))
  # add common legend
  legend <- cowplot::get_legend(
    # create some space to the left of the legend
    rstan::traceplot(stanfit1, pars = pars[1], inc_warmup = inc_warmup) + 
      theme(legend.box.margin = margin(0, 0, 0, 10))
  )
  
  plotgrid <- cowplot::plot_grid(titles, plotgrid, ncol = 1, align = 'v', rel_heights = c(0.07,1))
  cowplot::plot_grid(plotgrid, legend, rel_widths = c(3, .3))
}


compare_traceplots_quick_fix <- function(stanfit1, stanfit2, pars, inc_warmup = FALSE){
  plot1 <- rstan::traceplot(stanfit1, pars = pars, inc_warmup = inc_warmup,
                            ncol = 1) +
    theme(legend.position = "none")
  plot2 <- rstan::traceplot(stanfit2, pars = pars, inc_warmup = inc_warmup,
                            ncol = 1)
  
  gridExtra::grid.arrange(plot1,plot2, ncol=2)
}
