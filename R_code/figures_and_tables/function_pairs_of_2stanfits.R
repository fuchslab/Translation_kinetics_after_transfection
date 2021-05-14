#' Pairs plot of two stanfit objects
#'
#' @param stanfit1 stanfit object for which plots will be shown in the lower and the diagonal panels
#' @param stanfit2 stanfit object for which plots will be shown in the upper and the diagonal panels
#' @param name1 character string with heading for the stanfit1
#' @param name2 character string with heading for the stanfit2
#' @param pars character vector of parameter names that should be contained in both stanfit objects stanfit1 and stanfit2
#' @param true_values vector of the true parameter values with which the data was generated. default = NA
#' @param pars_names character vector or vector of expressions giving the parameter names as to be printed in the plot
#' @param inc_boxplot Logical scalar indicating whether boxplots should be included
#' @param col1 color for 1d-densities and boxplots for stanfit1
#' @param col2 color for 1d-densities and boxplots for stanfit2
#' @param v_col1 color(s) for 2d-densites for stanfit1
#' @param v_col2 color(s) for 2d-densites for stanfit2
#' @param cex_title scaling factor for the titles given by name1 and name2
#'
#' @return
#' @export
#'
#' @examples
#'

##  Define default colors ------------------------------------------------------------
HMGU_red <- '#E50030'
TUM_blue <- '#3070B3'
Helmholtz_green <- '#8CB423'
Helmholtz_blue <- '#005aa0'
Helmholtz_dark_blue <- '#0a2D6e'
Helmholtz_green_rgb <- c(140,180,35)
Helmholtz_blue_rgb <- c(0,90,160)
Helmholtz_dark_blue_rgb <- c(10,45,110)

##  main function --------------------------------------------------------------------
pairs_of_2stanfits <- function(stanfit1, stanfit2 = NULL, 
                               name1 = "stanfit1", name2=NULL,
                               pars = c("theta[1]", "theta[3]"), 
                               true_values = rep(NA, length(pars)),
                               pars_names = pars,
                               inc_boxplot = FALSE, incl_priors = FALSE, 
                               col1 = Helmholtz_blue, col2 = Helmholtz_green, 
                               v_col1 = col1, v_col2 = col2,
                               cex_title = 1.2){
  # check whether samples for specified pars exist
  if(! all(pars %in% dimnames(stanfit1)[[3]])){
    stop(paste0('There are no samples for parameter "',
                pars[! pars %in% dimnames(stanfit1)[[3]]], '" in stanfit1.'))
  }
  if(!is.null(stanfit2) && ! all(pars %in% dimnames(stanfit2)[[3]])){
    stop(paste0('There are no samples for parameter "',
                pars[! pars %in% dimnames(stanfit1)[[3]]], '" in stanfit2.'))
  }
  # extract the MCMC samples from the stanfit objects and convert them to a matrix
  sample1 <- do.call(cbind, rstan::extract(stanfit1, pars = pars))
  if(!is.null(stanfit2)){
    sample2 <- do.call(cbind, rstan::extract(stanfit2, pars = pars))
  }else{
    sample2 <- NULL
  }
  
  if(incl_priors){
    priors <-  readRDS("intermediate_output_files/prior_dens_fct_or_sample.rds")
  }
  
  # determine the plot range for each parameter
  plot_range <- apply(rbind(sample1, sample2, true_values), MARGIN = 2, 
                      FUN = range, na.rm = TRUE)
  
  num_par <- length(pars)
  layout(matrix(1:(num_par^2), ncol = num_par))

  # determine margins and title lines
  # margins needed for axis labels
  log10_of_max_of_range_max <- log(max(plot_range), base=10)
  log10_of_min_of_range_max <- log(min(plot_range[2,]), base=10)
  num_digits_max <- ceiling(log10_of_max_of_range_max)
  num_digits_min <- ceiling(abs(log10_of_min_of_range_max))
  
  if(log10_of_min_of_range_max > 0 || 
     num_digits_max >  num_digits_min + 2){
    ref_string <- paste0(rep("0", num_digits_max), collapse = "")
  }else{
    ref_string <- paste("0.", paste0(rep("0", num_digits_min + 1), collapse = ""))
  }
  label_margin_inches <- strwidth(ref_string, units="inches")
  line_height_inches <- par("csi")
  
  # margins for titles
  cex_title <- cex_title
  title1_line <- ceiling((label_margin_inches / line_height_inches +
                          par("tcl") + par("mgp")[2]) * 10) / 10
  
  name1_height_inches <- strheight(name1, units = "inches") * cex_title
  name2_height_inches <- strheight(name2, units = "inches") * cex_title
  
  # set outer margins
  oma_2 <- title1_line + name1_height_inches / line_height_inches + 0.8 * cex_title
  oma_3 <- name2_height_inches / line_height_inches + 1 * cex_title
  oma_4 <- title1_line + 0.1
  
  par(oma = c(1.5, oma_2, 1.5, oma_4), 
      mar= .4 * c(1, 1, 1, 1), 
      tcl = -0.4)
  
  for(k in 1:(num_par^2)){
    j <- ceiling(k / num_par)
    i <- k - num_par * (j-1)
    incl_x_lab <- ifelse(i == num_par, TRUE, FALSE)
    incl_y_lab <- ifelse(j == 1 || j == num_par, TRUE,  FALSE)
    
    if(i == j){
      true_value <- ifelse(is.null(true_values), NA, true_values[i])
      if(is.null(sample2)){
        s2 <- NULL
      }else{
        s2 <- sample2[ ,i]
      } 
      if(incl_priors){
        prior <-  priors[[param[i]]]
      }else{
        prior <- NULL
      }
      density_1d_plots(sample1[ ,i], s2, p = pars[i], 
                       pars_name = pars_names[i], true_value = true_value,
                       xlim = plot_range[,i], inc_boxplot = inc_boxplot,
                       col1 = col1, col2 = col2, incl_x_lab = incl_x_lab,
                       prior = prior)
    }else{
      true_val <- ifelse(rep(is.null(true_values),2), NA, true_values[c(i,j)])
      
      if(i > j){# lower triangle
        density_2d_plots(sample1[ ,i], sample1[ ,j], pars[c(i,j)],
                         xlim = plot_range[,j], ylim = plot_range[,i],
                         col_vec = v_col1, true_values = true_val,
                         pos_y_lab = 2, incl_y_lab = incl_y_lab, 
                         incl_x_lab = incl_x_lab)
        
      }else if(i < j){# upper triangle
        if(!is.null(stanfit2)){
          density_2d_plots( sample2[ ,i], sample2[ ,j],pars[c(i,j)],
                            xlim = plot_range[,j], ylim = plot_range[,i],
                            col_vec = v_col2, true_values = true_val,
                            pos_y_lab = 4, incl_y_lab = incl_y_lab, 
                            incl_x_lab = incl_x_lab)
        }else{
          plot(0,type='n',axes=FALSE,ann=FALSE)
        }
      }
    }
  }
  if(!is.null(as.character(name1))) 
    mtext(name1, side = 2, line = title1_line, outer = TRUE, cex = cex_title)
  if(!is.null(as.character(name2))) 
    mtext(name2, side = 3, line = -0.2, outer = TRUE, cex = cex_title)

}

density_1d_plots <- function(s1, s2, p, pars_name = p, xlim, true_value = NA, 
                             inc_boxplot = FALSE,
                             col1 = TUM_blue, col2 = Helmholtz_green,
                             incl_x_lab = TRUE, prior = NULL){
  dens1 <- density(s1)
  if(!is.null(s2)){
    dens2 <- density(s2)
    bw <- max(dens1$bw, dens2$bw)
    dens1 <- density(s1, bw = bw)
    dens2 <- density(s2, bw = bw)
    max_y <- max(dens1$y, dens2$y) 
  }else{
    bw <- dens1$bw
    max_y <- max(dens1$y) 
  }
  if(names(dev.cur()) == "png"){
    lwd <- 2
    lwd_bp_median <- 3
  }else{
    lwd <- 1.5
    lwd_bp_median <- 2
  }
  if(inc_boxplot){
    min_y <- - 0.3 * max_y
    scaling_factor_max_y  <- 1.5
    scaling_factor_true_val_line <- 1.2
  }else{
    min_y <- 0
    scaling_factor_max_y <- 1 + 0.5/1.8
    scaling_factor_true_val_line <- 1.05
  }
  ylim <- c(min_y, max_y * scaling_factor_max_y)
  plot(dens1, xlim = xlim, ylim = ylim, main = "", xlab = "", ylab = "",
       col = col1, las = 1,  yaxt='n',  xaxt='n', lwd = lwd) # bty='l',
  axis(side = 1, label = incl_x_lab)
  if(!is.null(s2)) lines(dens2, col = col2, lwd = lwd)
  text(mean(xlim), ylim[2] * 0.9, pars_name, cex = 1.5, font = 1) # font = 2 for bold
  
  if(inc_boxplot){
    if(names(dev.cur()) == "png"){
      abline(h = 0, col= "gray", lwd = .7)
    }
    draw_box_plot <- function(sample, positions, col){
      col_rgb <- col2rgb(col)
      pale_col <- rgb(col_rgb[1], col_rgb[2], col_rgb[3], maxColorValue = 255,
                       alpha = 90)
      quant_sample <- quantile(sample, probs = c(0,0.025,0.5,0.975,1))
      # draw box
      rect(xleft = quant_sample[2], ybottom = positions[1],
           xright = quant_sample[4], ytop = positions[3],
           lwd = lwd, col = pale_col, border = col)
      # draw whiskers
      segments(x0 = c(quant_sample[1], quant_sample[4]),
               y0 = rep(positions[2], 2),
               x1 = c(quant_sample[2], quant_sample[5]),
               y1 = rep(positions[2], 2),
               lwd = lwd, col = col)
      # draw median
      segments(x0 = quant_sample[3],
               y0 = positions[1],
               x1 = quant_sample[3],
               y1 = positions[3],
               lwd = lwd_bp_median)
    }
    draw_box_plot(s1, min_y * c(1/7, 2/7, 3/7), col1)
    if(!is.null(s2)) draw_box_plot(s2, min_y * c(4/7, 5/7, 6/7), col2)
  }
  if(!is.na(true_value)){
    window_coord <- par("usr")
    segments(x0 = true_value,
             y0 = window_coord[3],
             x1 = true_value,
             y1 = max_y * scaling_factor_true_val_line,
             lty = 3,
             lwd = lwd)
  }
  
  if(!is.null(prior)){
    if(prior$dens_fct_exists){
      x <- seq(from = xlim[1], to = xlim[2], length = 100)
      lines(x, prior$densfct(x), lty = 2, lwd = lwd * 0.9)
    }else{
      ind_in_range <- prior$sample > xlim[1] & prior$sample < xlim[2] 
      proportion_include <- sum(ind_in_range) / length(prior$sample)
      prior_dens <- density(prior$sample[ind_in_range], bw = bw)
      lines(prior_dens$x, prior_dens$y * proportion_include, lty = 2, 
            lwd = lwd * 0.9)
    }
  }
}

density_2d_plots <- function(s_p1, s_p2, p, ylim, xlim, col_vec = "black",
                             true_values = NA, pos_y_lab = 2, incl_y_lab =TRUE,
                             incl_x_lab = TRUE){
  colramp = colorRampPalette(c("white",  col_vec))
  smoothScatter(s_p2, s_p1, nrpoints = 0, xlab = "", ylab = "", colramp = colramp,
                ylim = ylim, xlim = xlim, las = 1, yaxt='n', xaxt='n')
  axis(side = 1, label = incl_x_lab)
  axis(side = pos_y_lab, las = 1, label = incl_y_lab)
  if(!all(is.na(true_values))){
    if(names(dev.cur()) == "png"){
      lwd <- 2
    }else{
      lwd <- 1.5
    }
    abline(v = true_values[2], lty = 3, lwd = lwd)
    abline(h = true_values[1], lty = 3, lwd = lwd)
  }
}

