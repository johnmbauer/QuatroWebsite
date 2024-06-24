# Script to analyze E. coli microplate growth curves as measured through OD600
# This version used the updated data and the growthrates package to analyze the
# growth curve data, instead of doing by hand

# Import the libraries
library(tidyverse)
library(growthrates)
library(growthcurver)
library(ggpubr)

# Import the data and summarize
Ecoli_growth_data <- read_csv('Ecoli_plate_growth_data_v2.csv') %>% 
  pivot_longer(cols = -c(Strain, CRISPRi_system, Replicate, Time),
               names_to = 'Inducers',
               values_to = 'OD') %>% 
  mutate(name = ifelse(Inducers == 'WT', 'WT', 
                       ifelse(CRISPRi_system == 'Sp', paste0('pSR2017, ', Inducers, ' ng/mL aTc'),
                              ifelse(CRISPRi_system == 'Fn', paste0('pSR2022, ', Inducers, ' ng/mL aTc'),
                                     paste0('pSR2032, ', Inducers, ' ng/mL aTc')))))
Ecoli_growth_data[Ecoli_growth_data$Strain == 'Nissle', ]$Strain <- 'EcN'

Ecoli_growth_summary <- Ecoli_growth_data %>% 
  group_by(Strain, CRISPRi_system, Inducers, name, Time) %>% 
  summarize(OD_mean = mean(OD),
            OD_sd = sd(OD))

# Define some base themes for different types of plots (same as all others)
text_size = 11
base_cat_theme <- theme_classic() +
  theme(axis.title = element_text(face = 'plain', size = text_size + 2),
        axis.text = element_text(size = text_size, color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        axis.text.y = element_text(vjust = 0.2),
        axis.ticks.length = unit(0.2, 'cm'),
        legend.title = element_text(face = 'plain', size = text_size),
        legend.position = 'bottom',
        plot.title = element_text(face = 'plain', size = text_size + 3, hjust = 0.5)
  )
base_cont_theme <- theme_classic() +
  theme(axis.title = element_text(face = 'plain', size = text_size + 2),
        axis.text = element_text(size = text_size, color = 'black'),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.3),
        axis.text.y = element_text(vjust = 0.2),
        axis.ticks.length = unit(0.2, 'cm'),
        legend.title = element_text(face = 'plain', size = text_size),
        legend.position = 'bottom',
        plot.title = element_text(face = 'plain', size = text_size + 3, hjust = 0.5)
  )
facet_theme <- theme(strip.text = element_text(hjust = 0, vjust = 0, size = text_size + 2, face = 'plain'),
                     strip.background = element_rect(color = NA),
)

mg_color <- '#D55E00' # '#fbc46a'
ecn_color <- '#CC79A7' # '#aa73ab'
cft_color <- '#0072B2' #'#7f8fc4'
umn_color <- '#009E73' # '#62ae83'
Sp_color <- '#E69F00' # '#2e3092'
Fn_color <- '#F0E442' # '#e57525'
Lb_color <- '#56B4E9' # '#e01a8e'

strain_clrs <- c(mg_color, ecn_color, cft_color, umn_color)
system_clrs <- c(Fn_color, Lb_color, Sp_color)


# Create a function to add significant marks to a plot with 

calc_sig_pos_dodge <- function(data, group_var, dodge_var, y_var, 
                               label_col, combos = 'all', scale = 'linear',
                               width = 0.9, buffer = 0,
                               line_y_add = 1, line_y_mult = 1, 
                               text_y_add = 0.05, text_y_mult = 1,
                               tip_len = 0) {
  if (!(combos %in% c('all', 'left', 'right'))) {
    stop('combos must be one of the following: all, left, or right')
  }
  if (!(scale %in% c('linear', 'log'))) {
    stop('scale must be one of the following: linear, log')
  }
  line_tibble <- tibble(set = as.character(),
                        x = as.numeric(),
                        y = as.numeric())
  label_tibble <- tibble(set = as.character(),
                         x_text = as.numeric(),
                         y_text = as.numeric(),
                         label = as.character())
  # Determine order of x-axis variables for plotting
  lvls <- levels(pull(data, group_var)) 
  if (is.null(lvls)) {lvls <- sort(unique(pull(data, group_var)))}
  # Determine the ordering for the dodging for each x-axis variable
  dodge_order <- levels(pull(filter(data, !!sym(group_var) == lvls[1]), dodge_var))
  if (is.null(dodge_order)) {
    dodge_order <- sort(unique(pull(filter(data, !!sym(group_var) == lvls[1]), dodge_var)))
  }
  dodge_tibble <- tibble(!!dodge_var := dodge_order,
                         order = seq(length(dodge_order)))
  # Determine number of lines that must be calculated
  n_vars <- length(pull(filter(data, !!sym(group_var) == lvls[1]), y_var))
  if (combos == 'all') {n_cond = choose(n_vars, 2)} else {n_cond = n_vars - 1}
  spacing <- width/n_vars 
  i <- 1 # Track position of the grouping variable
  k <- 1 # Track the set for plotting the lines
  for (lvl in lvls) {
    j <- 1 # Track position of line relative to the maximum y
    df <- filter(data, !!sym(group_var) == lvl) %>% ungroup()
    max_y <- max(pull(df, y_var))
    if (combos == 'all') {
      start <- c()
      end <- c()
      for (z in seq(1,n_cond-1)) {
        start <- c(start, rep(z, n_cond-z))
        end <- c(end, seq(z+1, n_cond))
      }
      # Map the labels to the positions
      
    } else {
      start <- seq(n_cond)
      end <- rep(n_vars, times = n_cond)
      if (combos == 'left') {
        start <- n_vars - start +1
        end <- rep(1, times = n_cond)
      }
      # Map the labels to the positions
      temp_text <- select(df, c(!!sym(dodge_var), !!sym(label_col))) %>% 
        left_join(dodge_tibble, by = names(dodge_tibble)[names(dodge_tibble) %in% names(df)])
    }
    temp_tibble <- tibble(x1 = start, x2 = end,
                          diff = abs(end - start)) %>% 
      arrange(diff)
    if (combos %in% c('left', 'right')) {
      temp_tibble <- temp_tibble %>% 
        left_join(rename(temp_text, x1 = order), by = 'x1')
    } else {
      temp_tibble <- temp_tibble %>% 
        left_join(rename(temp_text, x2 = order), by = 'x2')
    }
    
    # Need to map the positions (via temp tibble) to the text labels via the 
    # filtered df. However, this depends on the type of combinations
    if (tip_len != 0) {
      for (row in seq(length(temp_tibble$x1))) {
        temp <- temp_tibble[row,]
        a <- c(rep(temp$x1, times = 2), rep(temp$x2, times = 2))
        if (scale == 'linear') {
          b <- max_y + buffer + (c(-tip_len, 0, 0, -tip_len) + j*line_y_add)*line_y_mult
        } else {
          b <- max_y + buffer + 10**c(c(-tip_len, 0, 0, -tip_len) + log10(j*line_y_add+line_y_mult))
        }
        x_ <- i - width/2 + spacing*(a-0.5)
        line_tibble <- rbind(line_tibble,
                             tibble(set = k,
                                    x = x_,
                                    y = b))
        label_tibble <- rbind(label_tibble,
                              tibble(set = k, 
                                     x_text = mean(c(x_[2], x_[3])),
                                     y_text = b[2]*text_y_mult + text_y_add,
                                     label = temp$label))
        j <- j + 1
        k <- k + 1
      } 
    } else {
      for (row in seq(length(temp_tibble$x1))) {
        temp <- temp_tibble[row,]
        a <- i - width/2 + spacing*(c(temp$x1, temp$x2)-0.5)
        if (scale == 'linear') {
          #b <- rep(c(max_y + buffer + j*line_y_add)*line_y_mult, times = 2)
          b <- max_y + buffer + rep(j*line_y_add*line_y_mult, times = 2)
        } else {
          b <- max_y + buffer + rep(10**j, times = 2)
        }
        line_tibble <- rbind(line_tibble,
                             tibble(set = k,
                                    x = a,
                                    y = b))
        label_tibble <- rbind(label_tibble,
                              tibble(set = k, 
                                     x_text = mean(a),
                                     y_text = b[1]*text_y_mult + text_y_add,
                                     label = temp$label))
        j <- j + 1
        k <- k + 1
      }
    }
    i <- i + 1
  }
  return_tibble <- full_join(line_tibble, label_tibble, by = 'set')
  return_tibble
}



# Fit the data to the exponential growth equation from the growthrates package
# Testing with single sample first 
test_data <- filter(Ecoli_growth_data, Strain == 'MG1655', 
                    Inducers == 'WT', Replicate == 1, CRISPRi_system == 'Sp')

test_exp_fit <- fit_easylinear(test_data$Time, test_data$OD)
summary(test_exp_fit)

# Fit looks relatively well, so fit to all the datapoints using the all_easylinear function
all_fits <- all_easylinear(OD ~ Time | Strain + CRISPRi_system + Replicate + Inducers, 
                           data = Ecoli_growth_data)

# Loop through the data and fit data to extract the results while also 
# calculating model points for each dataset for plotting
exp_fit_params <- tibble(Strain = as.character(),
                         CRISPRi_system = as.character(),
                         Replicate = as.character(),
                         Inducers = as.character(),
                         name = as.character(),
                         pts = as.numeric(),
                         x0 = as.numeric(),
                         mu = as.numeric(),
                         lag = as.numeric(),
                         r2 = as.numeric())
model_fits <- tibble(Strain = as.character(),
                     CRISPRi_system = as.character(),
                     Replicate = as.character(),
                     Inducers = as.character(),
                     name = as.character(),
                     x = as.numeric(),
                     y = as.numeric())

for (strain in unique(Ecoli_growth_data$Strain)) {
  for (system in unique(Ecoli_growth_data$CRISPRi_system)) {
    for (rep in unique(Ecoli_growth_data$Replicate)) {
      for (set in unique(Ecoli_growth_data$Inducers)) {
        nm <- filter(Ecoli_growth_data, Inducers == set,
                     CRISPRi_system == system)$name[1]
        id <- paste(strain, system, rep, set, sep = ':')
        fit <- all_fits@fits[[id]]
        pts <- (fit@ndx-1)/4
        x0 <- coef(fit)[2]
        mu <- coef(fit)[3]
        lag <- coef(fit)[4]
        r2 <- fit@rsquared
        x <- seq(min(pts), max(pts), length.out = 10)
        y <- x0 * exp(mu * x)
        exp_fit_params <- rbind(exp_fit_params,
                                tibble(Strain = strain,
                                       CRISPRi_system = system,
                                       Replicate = rep,
                                       Inducers = set,
                                       name = nm,
                                       pts = paste0(pts, collapse = '-'),
                                       x0 = x0,
                                       mu = mu,
                                       lag = lag,
                                       r2 = r2))
        model_fits <- rbind(model_fits,
                            tibble(
                              Strain = strain,
                              CRISPRi_system = system,
                              Replicate = rep,
                              Inducers = set,
                              name = nm,
                              x = x,
                              y = y))
      }
    }
  }
}

#write_csv(exp_fit_params, 'Ecoli_growth_fits_v2.csv')

# Iterate through each growth dataset (Strain + system) and make plots
# of each replicate with the fits to check that they look fine
Sp_clrs <- c(RColorBrewer::brewer.pal(9, 'OrRd')[2:4], Sp_color, 'gray50')
Fn_clrs <- c(RColorBrewer::brewer.pal(9, 'YlOrRd')[1:3], Fn_color, 'gray50')
Lb_clrs <- c(RColorBrewer::brewer.pal(9, 'Blues')[2:5], 'gray50')

# for (strain in unique(Ecoli_growth_data$Strain)) {
#   for (system in unique(Ecoli_growth_data$CRISPRi_system)) {
#     rep_data <- filter(Ecoli_growth_data, Strain == strain,
#                        CRISPRi_system == system)
#     model_data <- filter(model_fits, Strain == strain,
#                          CRISPRi_system == system)
#     clrs <- if (system == 'Fn') Fn_clrs else if (system == 'Lb') Lb_clrs else Sp_clrs
#     plt <- ggplot(rep_data,
#                   aes(x = Time, y = OD, fill = Inducers)) +
#       geom_point(size = 0.75, shape = 21, color = 'black',
#                  stroke = 0.25) +
#       geom_line(data = model_data,
#                 mapping = aes(x = x, y = y, color = Inducers),
#                 color = 'black', linewidth = 0.65, alpha = 0.75) +
#       scale_y_log10(limits = c(0.01, 1),
#                     labels = scales::trans_format('log10', scales::math_format(10^.x))) +
#       scale_fill_manual(values = clrs) +
#       annotation_logticks(sides = 'l', outside = T) +
#       coord_cartesian(clip = 'off') +
#       base_cont_theme +
#       theme(legend.position = 'none') +
#       facet_wrap(~name + Replicate, scales = 'free', nrow = 5) +
#       facet_theme + 
#       theme(strip.text = element_text(size = text_size),
#             axis.text = element_text(size = text_size - 2),
#             axis.title = element_text(size = text_size))
#     filename <- paste0('images/', paste(strain, system, 'growth_rep_plots.pdf', sep = '_'))
#     ggsave(filename, plt, units = 'in', height = 7, width = 6, dpi = 800)
#   }
# }

# Fits look good, so continue with summarizing the fit data and performing 
# Welch's one-sided t-test relative to WT for both growth rate and lag phase
exp_fit_params_summary  <- exp_fit_params %>% 
  group_by(Strain, CRISPRi_system, Inducers, name) %>% 
  summarize(mu_avg = mean(mu),
            mu_sd = sd(mu),
            lag_avg = mean(lag),
            lag_sd = sd(lag))

#write_csv(exp_fit_params_summary, 'Ecoli_fits_summary_v2.csv')

t_tests <- tibble(Strain = as.character(),
                  CRISPRi_system = as.character(),
                  Set = as.character(),
                  name = as.character(),
                  mu_pval = as.numeric(),
                  lag_pval = as.numeric())

for (strain in unique(exp_fit_params$Strain)) {
  for (system in unique(exp_fit_params$CRISPRi_system)) {
    for (set in unique(exp_fit_params$Inducers)) {
      data <- filter(exp_fit_params, Strain == strain, 
                        CRISPRi_system == system, Inducers == set)
      wt_data <- filter(exp_fit_params, Strain == strain, 
                        CRISPRi_system == system, Inducers == 'WT')
      nm <- data$name[1]
      test_mu <- data$mu
      test_lag <- data$lag
      wt_mu <- wt_data$mu
      wt_lag <- wt_data$lag
      mu_test <- t.test(test_mu, wt_mu, paired = F, alternative = 'less')
      lag_test <- t.test(test_lag, wt_lag, paired = F, alternative = 'less')
      mu_pval <- mu_test$p.value
      lag_pval <- lag_test$p.value
      t_tests <- rbind(t_tests,
                       tibble(Strain = strain,
                              CRISPRi_system = system,
                              Set = set,
                              name = nm,
                              mu_pval = mu_pval,
                              lag_pval = lag_pval))
    }
  }
}

#write_csv(t_tests, 'Ecoli_growth_ttests.csv')

# Add asterisks for significance to the summarized exponential fit parameters 
# dataframe for plotting
exp_fits_summary <- exp_fit_params_summary %>% 
  full_join(rename(t_tests, Inducers = Set)) %>% 
  mutate(label = ifelse(mu_pval >= 0.05, 'ns',
                        ifelse(mu_pval > 0.01, '*',
                               ifelse(mu_pval > 0.001, '**', '***'))))
exp_fits_summary[exp_fits_summary$Inducers == 'WT',]$label <- ''
exp_fits_summary$Strain <- factor(exp_fits_summary$Strain, 
                                  levels = c('MG1655', 'EcN',
                                             'CFT073', 'UMN026')) 

# Plots of specific growth rate across strains for each system ---------------
# Make plots of the data fits for mu and the lag phase to visualize
#Sp_clrs <- c(RColorBrewer::brewer.pal(9, 'Blues')[5:8], 'gray20')
#Fn_clrs <- c(RColorBrewer::brewer.pal(9, 'OrRd')[3:6], 'gray20')
#Lb_clrs <- c(RColorBrewer::brewer.pal(9, 'RdPu')[3:6], 'gray20')
#Sp_clrs <- c(RColorBrewer::brewer.pal(9, 'OrRd')[3:6], 'gray20')
Sp_clrs <- c(RColorBrewer::brewer.pal(9, 'OrRd')[3:6], 'gray30')
Fn_clrs <- c(RColorBrewer::brewer.pal(9, 'YlOrRd')[1:3], Fn_color, 'gray30')
Lb_clrs <- c(RColorBrewer::brewer.pal(9, 'Blues')[3:6], 'gray30')
#Lb_clrs <- c(RColorBrewer::brewer.pal(9, 'Blues')[2:4], Lb_color, 'gray20')



for (system in unique(exp_fits_summary$CRISPRi_system)) {
  data_summary <- filter(exp_fits_summary, CRISPRi_system == system)
  data <- filter(exp_fit_params, CRISPRi_system == system)
  clrs <- if (system == 'Fn') Fn_clrs else if (system == 'Lb') Lb_clrs else Sp_clrs
  plot <- ggplot(data_summary, 
                 aes(x = Strain, y = mu_avg, color = Inducers)) +
    geom_point(size = 5, shape = '\u2013', 
               position = position_dodge(width = 0.75)) +
    geom_errorbar(aes(ymin = (mu_avg - mu_sd), ymax = (mu_avg + mu_sd)),
                  width = 0.2, position = position_dodge(width = 0.75)) +
    geom_point(data = data,
               mapping = aes(x = Strain, y = mu, color = Inducers),
               shape = 16, size = 1.2,
               position = position_jitterdodge(dodge.width = 0.75,
                                               jitter.width = 0.08)) +
    geom_text(aes(x = Strain, y = mu_avg+0.1, color = Inducers, label = label), 
              position = position_dodge(width  = 0.75), size = 4) +
    scale_color_manual(values = clrs, 
                       guide = guide_legend(override.aes = list(label = '',
                                                                shape = '\u2014',
                                                                size = 6))) +
    scale_y_continuous(limits = c(0, 1.0),
                       expand = ) +
    labs(x = 'Strain',
        y = 'Specific growth rate (1/hr)',
        color = 'ng/mL aTc') +
    base_cat_theme
  plotname <- paste0('images/', paste(sym(system), 'mu_plots_v2.svg', sep = '_'))
  #ggsave(plotname, plot, units = 'in', height = 4.5, width = 5)
}

line_tibble <- tibble(x = as.numeric(),
                      xend = as.numeric(),
                      y = as.numeric(),
                      yend = as.numeric())
text_tibble <- tibble(x = as.numeric(), 
                      y = as.numeric(),
                      Strain = as.character(),
                      CRISPRi_system = as.character(),
                      label = as.character())
i <- 0
for (strain in c('MG1655', 'EcN',
                 'CFT073', 'UMN026')) {
  data <- filter(exp_fits_summary, Strain == strain)
  n_conditions <- length(filter(data,
                                CRISPRi_system == 'Sp')$Inducers)
  labels <- filter(data, Inducers != 'WT')
  spacing <- (1.4 - 0.6)/n_conditions
  a <- i + seq(1-spacing*2, 1+spacing*1, by = spacing)
  b <- i + rep(1.4-spacing*0.5, times = 4)
  text_x <- a + (b - a)/2
  min_y <- unique(filter(exp_fit_params_summary, Strain == strain,
                  Inducers == 'WT')$mu_avg) + 0.15
  y <- c(0, 0.075, 0.15, 0.225) + rep(min_y, times = 4)
  text_y <- y + 0.035
  line_tibble <- rbind(line_tibble,
                       tibble(x = a, xend = b, y = y, yend = y))
  text_tibble <- rbind(text_tibble, 
                       tibble(x = rep(text_x, times = 3), 
                              y = rep(text_y, times = 3), 
                              Strain = strain,
                              CRISPRi_system = labels$CRISPRi_system,
                              label = labels$label))
  i <-  i + 1
}


# Make the same plots as bar charts
for (system in unique(exp_fits_summary$CRISPRi_system)) {
  data_summary <- filter(exp_fits_summary, CRISPRi_system == system)
  data <- filter(exp_fit_params, CRISPRi_system == system)
  clrs <- if (system == 'Fn') Fn_clrs else if (system == 'Lb') Lb_clrs else Sp_clrs
  # Calculate the positions for the lines and scales for the stat results
  buf <- max(filter(data_summary, Inducers == 'WT')$mu_sd)
  scales_tibble <- calc_sig_pos_dodge(data_summary, 'Strain', 'name', 'mu_avg', 'label',
                                      combos = 'right', line_y_add = 0.05, text_y_add = 0.02,
                                      buffer = buf)
  plot <- ggplot(data_summary, 
                 aes(x = Strain, y = mu_avg, fill = name)) +
    geom_bar(stat = 'identity', color = 'black',
               position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = (mu_avg - mu_sd), ymax = (mu_avg + mu_sd)),
                  width = 0.2, position = position_dodge(width = 0.9)) +
    geom_point(data = data,
               mapping = aes(x = Strain, y = mu, color = name),
               shape = 16, size = 1.2, color = 'black',
               position = position_jitterdodge(dodge.width = 0.9,
                                               jitter.width = 0.08),
               show.legend = F) +
    # geom_text(aes(x = Strain, y = mu_avg+0.14, color = Inducers, label = label), 
    #           position = position_dodge(width  = 0.9), size = 3.5, color = 'black') +
   scale_fill_manual(values = clrs,
                      guide = guide_legend(override.aes = list(label = ''))) +
    scale_color_manual(values = clrs, 
                      guide = 'none') +
    scale_y_continuous(limits = c(0, 1.1),
                       expand = c(0,0)) +
    labs(x = 'Strain',
         y = 'Specific growth rate (1/hr)',
         fill = '') +
    base_cat_theme +
    theme(legend.position = 'right')
  # Add stat results to plot 
  for (st in unique(scales_tibble$set)) {
    stat_data <- filter(scales_tibble, set == st)
    plot <- plot + 
      geom_line(data = stat_data, 
                mapping = aes(x = x, y = y), linewidth = 0.75, inherit.aes = F) +
      geom_text(data = stat_data[1,],
                mapping = aes(x = x_text, y = y_text, label = label),
                size = 3, inherit.aes = F)
  }
  # for (row in seq(length(line_tibble$x))) {
  #   line_data <- line_tibble[row,]
  #   plot <- plot +
  #     annotate('line', x = c(line_data$x, line_data$xend),
  #              y = c(line_data$y, line_data$yend), linewidth = 0.8)
  # }
  # text_labels <- filter(text_tibble, CRISPRi_system == system)
  # for (row in seq(length(text_labels$x))) {
  #   text_data <- text_labels[row,]
  #   plot <- plot +
  #     annotate('text', x = text_data$x, y = text_data$y,
  #              label = text_data$label, size = 3.5)
  # }
  plot
  plotname <- paste0('images/', paste(sym(system), 'mu_barplots_v2.3.svg', sep = '_'))
  ggsave(plotname, plot, units = 'in', height = 4, width = 6)
}


# OD at specific time point -----------------------------------------
# Extract data at a specific time point (4 hrs)
growth_data_4hr <- filter(Ecoli_growth_data, Time == 4)

growth_summary_4hr <- growth_data_4hr %>% 
  group_by(Strain, CRISPRi_system, Time, Inducers, name) %>% 
  summarize(OD_mean = mean(OD),
            OD_sd = sd(OD),
            n_reps = n())

# Run t-tests compared to the WT (like the growth curve data) 

# Comparing ANOVA with Dunnett's post-hoc test to t-tests for this dataset
# Dunnett's feels more appropriate because multiple comparisons, but not 
# sure in Lauren will want this
t_tests_OD <- tibble(Strain = as.character(),
                     CRISPRi_system = as.character(),
                     Set = as.character(),
                     od_pval = as.numeric())

dunnett_OD <- tibble(Strain = as.character(),
                     CRISPRi_system = as.character(),
                     inducer_set = as.character(),
                     pval = as.numeric())


for (strain in unique(growth_data_4hr$Strain)) {
  for (system in unique(growth_data_4hr$CRISPRi_system)) {
    for (set in unique(growth_data_4hr$Inducers)) {
      # t tests
      data <- filter(growth_data_4hr, Strain == strain, 
                     CRISPRi_system == system, Inducers == set)
      wt_data <- filter(growth_data_4hr, Strain == strain, 
                        CRISPRi_system == system, Inducers == 'WT')
      test_od <- data$OD
      wt_od <- wt_data$OD
      od_test <- t.test(test_od, wt_od, paired = F, alternative = 'less')
      od_pval <- od_test$p.value
      t_tests_OD <- rbind(t_tests_OD,
                       tibble(Strain = strain,
                              CRISPRi_system = system,
                              Set = set,
                              od_pval = od_pval))
    }
    # ANOVA and Dunnett's post-hoc test, if needed
    aov_data <- filter(exp_fit_params, Strain == strain, 
                       CRISPRi_system == system)
    anova <- aov(mu ~ Inducers, data = aov_data)
    if (summary(anova)[[1]][['Pr(>F)']][1] < 0.05) {
      dun_test <- DescTools::DunnettTest(mu ~ Inducers, data = aov_data, control = 'WT')
      df <- data.frame(dun_test[[1]])
      dunnett_OD <- rbind(dunnett_OD, 
                          tibble(Strain = strain,
                                 CRISPRi_system = system,
                                 inducer_set = rownames(df),
                                 pval = df[,'pval']))
    }
  }
}

# Add the p-value marker to the OD tibble
growth_summary_4hr <- growth_summary_4hr %>% 
  left_join(rename(t_tests_OD, Inducers = Set), 
            by = c('Strain','CRISPRi_system', 'Inducers')) %>% 
  mutate(label = ifelse(od_pval >= 0.05, 'ns',
                        ifelse(od_pval > 0.01, '*',
                               ifelse(od_pval > 0.001, '**', '***'))))
growth_summary_4hr$Strain <- factor(growth_summary_4hr$Strain, 
                                    levels = c('MG1655', 'EcN', 'CFT073', 'UMN026'))

# Make a plot of the OD data for each CRISPRi system
for (system in unique(growth_summary_4hr$CRISPRi_system)) {
  data_summary <- filter(growth_summary_4hr, CRISPRi_system == system)
  data <- filter(growth_data_4hr, CRISPRi_system == system)
  clrs <- if (system == 'Fn') Fn_clrs else if (system == 'Lb') Lb_clrs else Sp_clrs
  # Calculate stat results line and text data 
  buf <- max(filter(data_summary, Inducers == 'WT')$OD_sd)
  scales_tibble <- calc_sig_pos_dodge(data_summary, 'Strain', 'name', 'OD_mean', 'label',
                                      combos = 'right', line_y_add = 0.03, text_y_add = 0.016,
                                      buffer = buf)
  plot <- ggplot(data_summary, 
                 aes(x = Strain, y = OD_mean, fill = name)) +
    geom_bar(stat = 'identity', color = 'black',
             position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = (OD_mean - OD_sd), ymax = (OD_mean + OD_sd)),
                  width = 0.2, position = position_dodge(width = 0.9)) +
    geom_point(data = data,
               mapping = aes(x = Strain, y = OD, color = name),
               shape = 16, size = 1.2, color = 'black',
               position = position_jitterdodge(dodge.width = 0.9,
                                               jitter.width = 0.08),
               show.legend = F) +
    # geom_text(aes(x = Strain, y = OD_mean+0.035, color = name, label = label),
    #           position = position_dodge(width  = 0.9), size = 3.5, color = 'black') +
    scale_fill_manual(values = clrs,
                      guide = guide_legend(override.aes = list(label = ''))) +
    scale_color_manual(values = clrs, 
                       guide = 'none') +
    scale_y_continuous(limits = c(0, 0.64),
                       expand = c(0,0)) +
    labs(x = 'Strain',
         y = 'Optical density (600 nm) at 4 hr',
         fill = '') +
    base_cat_theme +
    theme(legend.position = 'right')
  # for (row in seq(length(line_tibble$x))) {
  #   line_data <- line_tibble[row,]
  #   plot <- plot +
  #     annotate('line', x = c(line_data$x, line_data$xend),
  #              y = c(line_data$y, line_data$yend), linewidth = 0.8)
  # }
  # text_labels <- filter(text_tibble, CRISPRi_system == system)
  # for (row in seq(length(text_labels$x))) {
  #   text_data <- text_labels[row,]
  #   plot <- plot +
  #     annotate('text', x = text_data$x, y = text_data$y,
  #              label = text_data$label, size = 3.5)
  # }
  for (st in unique(scales_tibble$set)) {
    stat_data <- filter(scales_tibble, set == st)
    plot <- plot + 
      geom_line(data = stat_data, 
                mapping = aes(x = x, y = y), linewidth = 0.75, inherit.aes = F) +
      geom_text(data = stat_data[1,],
                mapping = aes(x = x_text, y = y_text, label = label),
                size = 3, inherit.aes = F)
  }
  plot
  plotname <- paste0('images/', paste(sym(system), 'od_4hr_barplots_v2.3.svg', sep = '_'))
  ggsave(plotname, plot, units = 'in', height = 4, width = 6)
}

# Determine stats for OD and specific growth rate ----------------
# Stats relative to the WT of each strain

growth_sum_stats <- select(ungroup(growth_summary_4hr), -Time) %>% 
  full_join(exp_fits_summary, 
            by = c('Strain', 'CRISPRi_system', 'Inducers', 'name')) %>%
  group_by(Strain, CRISPRi_system) %>% 
  mutate(wt_od = OD_mean[Inducers == 'WT'],
         wt_mu = mu_avg[Inducers == 'WT']) %>% 
  group_by(Strain, CRISPRi_system, Inducers, name) %>% 
  summarize(od_rel_WT = OD_mean / wt_od, 
            mu_rel_WT = mu_avg / wt_mu) 

write_csv(growth_sum_stats, 'stat_results/growth_rate_stats.csv')

#
# ANOVA across CRISPRi systems -------------------------------------------
# ANOVA for the highest concentration of aTc for each strain to determine 
# if there is one with the greatest toxicity
toxicity_tukey_results <- tibble(Strain = as.character(),
                                 Inducer_conc = as.character(),
                                 system_combo = as.character(),
                                 p.adj = as.numeric())

for (strain in unique(exp_fit_params$Strain)) {
  for (conc in c('0','0.25','2','4')) {
    data <- filter(exp_fit_params, Strain == strain, Inducers == conc)
    anova_res <- aov(mu ~ CRISPRi_system, data = data)
    if (summary(anova_res)[[1]][["Pr(>F)"]][1] < 0.05) {
      tukey <- TukeyHSD(anova_res, 'CRISPRi_system')
      df <- data.frame(tukey$CRISPRi_system)
      toxicity_tukey_results <- rbind(toxicity_tukey_results,
                                      tibble(Strain = strain, 
                                             Inducer_conc = conc,
                                             system_combo = rownames(df),
                                             p.adj = df$p.adj))
    }
  }
}

write_csv(toxicity_tukey_results, 'tukey_results_systems_toxicity.csv')

# Plot of these data to display results
high_atc_data <- filter(exp_fits_summary, Inducers == '0')
high_atc_data[high_atc_data$CRISPRi_system == 'Sp',]$CRISPRi_system <- 'Sp-dCas9'
high_atc_data[high_atc_data$CRISPRi_system == 'Fn',]$CRISPRi_system <- 'Fn-dCas12a'
high_atc_data[high_atc_data$CRISPRi_system == 'Lb',]$CRISPRi_system <- 'Lb-dCas12a'
high_atc_data$CRISPRi_system <- factor(high_atc_data$CRISPRi_system,
                                       levels = c('Sp-dCas9', 'Fn-dCas12a', 
                                                  'Lb-dCas12a'))
high_atc_pts <- filter(exp_fit_params, Inducers == '0')
high_atc_pts[high_atc_pts$CRISPRi_system == 'Sp',]$CRISPRi_system <- 'Sp-dCas9'
high_atc_pts[high_atc_pts$CRISPRi_system == 'Fn',]$CRISPRi_system <- 'Fn-dCas12a'
high_atc_pts[high_atc_pts$CRISPRi_system == 'Lb',]$CRISPRi_system <- 'Lb-dCas12a'
high_atc_pts$CRISPRi_system <- factor(high_atc_pts$CRISPRi_system,
                                       levels = c('Sp-dCas9', 'Fn-dCas12a', 
                                                  'Lb-dCas12a'))

toxicity_across_systems_plot <- ggplot(high_atc_data,
                                       aes(x = Strain, y = mu_avg, fill = CRISPRi_system)) +
  geom_bar(stat = 'identity', color = 'black',
           position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = (mu_avg - mu_sd), ymax = (mu_avg + mu_sd)),
                width = 0.2, position = position_dodge(width = 0.9)) +
  geom_point(data = high_atc_pts,
             mapping = aes(x = Strain, y = mu, color = CRISPRi_system),
             shape = 16, size = 1.2, color = 'black', alpha = 0.75,
             position = position_jitterdodge(dodge.width = 0.9,
                                             jitter.width = 0.2),
             show.legend = F) +
  scale_fill_manual(values = c(Sp_color, Fn_color, Lb_color),
                    guide = guide_legend(override.aes = list(label = ''))) +
  scale_y_continuous(limits = c(0, 1.1),
                     expand = c(0,0)) +
  labs(x = 'Strain',
       y = 'Specific growth rate (1/hr)',
       fill = 'CRISPRi system') +
  base_cat_theme +
  theme(legend.position = 'right',
        legend.spacing.y = unit(1, 'mm'))
toxicity_across_systems_plot

# add significant marks to plot
toxicity_sig_labs <- toxicity_tukey_results %>% 
  filter(Inducer_conc == '0') %>% 
  mutate(label = ifelse(p.adj >= 0.05, 'ns',
                        ifelse(p.adj > 0.01, '*',
                               ifelse(p.adj > 0.001, '**', '***'))))

line_tibble <- tibble(x = as.numeric(),
                      xend = as.numeric(),
                      y = as.numeric(),
                      yend = as.numeric())
text_tibble <- tibble(x = as.numeric(), 
                      y = as.numeric(),
                      Strain = as.character(),
                      label = as.character())
i <- 0
for (strain in c('MG1655', 'EcN',
                 'CFT073', 'UMN026')) {
  data <- filter(high_atc_data, Strain == strain)
  n_conditions <- length(data$CRISPRi_system)
  temp <- filter(toxicity_sig_labs, Strain == strain)
  ordering <- c('Sp-Fn', 'Lb-Fn', 'Sp-Lb')
  labels <- temp[match(ordering, temp$system_combo),]
  spacing <- (1.45 - 0.55)/n_conditions
  a <- i + c(1-spacing, 1, 1-spacing)
  b <- i + c(1, 1+spacing, 1+spacing)
  text_x <- a + (b - a)/2
  min_y <- max(filter(high_atc_data, Strain == strain)$mu_avg) + 0.1
  y <- c(0, 0.06, 0.12) + rep(min_y, times = 3)
  text_y <- y + 0.035
  line_tibble <- rbind(line_tibble,
                       tibble(x = a, xend = b, y = y, yend = y))
  text_tibble <- rbind(text_tibble, 
                       tibble(x = rep(text_x), 
                              y = rep(text_y), 
                              Strain = strain,
                              label = labels$label))
  i <-  i + 1
}

text_tibble[is.na(text_tibble$label), ]$label <- 'ns'

for (row in seq(length(line_tibble$x))) {
  line_data <- line_tibble[row,]
  toxicity_across_systems_plot <- toxicity_across_systems_plot +
    annotate('line', x = c(line_data$x,line_data$x, line_data$xend,line_data$xend),
             y = c(line_data$y-0.01,line_data$y, line_data$yend,line_data$yend-0.01), linewidth = 0.8)
}
for (row in seq(length(text_tibble$x))) {
  text_data <- text_tibble[row,]
  toxicity_across_systems_plot <- toxicity_across_systems_plot +
    annotate('text', x = text_data$x, y = text_data$y,
             label = text_data$label, size = 3.5)
}
toxicity_across_systems_plot
ggsave('images/growth_rates/toxicity_across_systems_0ngml_plot.pdf', toxicity_across_systems_plot,
       units = 'cm', height = 12, width = 12)


# Growth correlations across strains ----------------------
growth_corr_results <- tibble(CRISPRi_system = as.character(),
                              Correlation = as.character(),
                              Strains = as.character(),
                              corr_val = as.numeric(),
                              p_val = as.numeric())

rel_mu <- ungroup(exp_fits_summary) %>% 
  select(Strain, CRISPRi_system, Inducers, mu_avg) %>%
  group_by(Strain, CRISPRi_system) %>% 
  mutate(wt_mu = mu_avg[Inducers == 'WT'],
         rel_mu = mu_avg / wt_mu) %>% 
  ungroup() %>% 
  filter(Inducers != 'WT')

strains <- c('MG1655', 'EcN', 'CFT073', 'UMN026')
i <- 2
for (system in c('Sp', 'Fn', 'Lb')) {
  df <- filter(rel_mu, CRISPRi_system == system) %>% 
    select(Inducers, Strain, mu_avg) %>%
    pivot_wider(names_from = Strain, values_from = mu_avg)
  for (j in seq(1, 3)) {
    strain1 <- strains[j]
    x_vals <- (pull(df, sym(strain1))) # Extract the column and convert to vector
    for (k in seq(i, 4)) {
      strain2 <- strains[k]
      y_vals <- (pull(df, sym(strain2)))
      strain_combo <- paste(strain1, strain2, sep = '-')
      # Correlations
      for (method in c('spearman', 'pearson')) {
        test <- cor.test(x_vals,
                         y_vals,
                         alternative = 'two.sided',
                         method = method)
        value <- test$estimate[[1]]
        p_val <- test$p.value[[1]]
        growth_corr_results <- rbind(growth_corr_results,
                                     tibble(CRISPRi_system = system,
                                            Correlation = method,
                                            Strains = strain_combo,
                                            corr_val = value,
                                            p_val = p_val))
      }
      # Linear fits
    #   frmla <- paste0(sym(strain2), ' ~ ', sym(strain1), ' - 1')
    #   lin_fit <- lm(frmla, data = df)
    #   coeff <- lin_fit$coefficients[[1]]
    #   r2 <- summary(lin_fit)$r.squared
    #   pval <- summary(lin_fit)$coefficients[,4]
    #   active_lm_fits <- rbind(active_lm_fits,
    #                           tibble(CRISPRi_system = system,
    #                                  Strains = strain_combo,
    #                                  b = coeff,
    #                                  r2 = r2,
    #                                  p_val = pval))
    }
    i <- i + 1
  }
  i <- 2
} 

#write_csv(growth_corr_results, './stat_results/growth_rate_correlations.csv')
Sp_clrs <- c(RColorBrewer::brewer.pal(9, 'OrRd')[3:6], 'gray30')
Fn_clrs <- c(RColorBrewer::brewer.pal(9, 'YlOrRd')[1:3], Fn_color, 'gray30')
Lb_clrs <- c(RColorBrewer::brewer.pal(9, 'Blues')[3:6], 'gray30')

for (system in unique(rel_mu$CRISPRi_system)) {
  df <- filter(rel_mu, )
}
ggplot(rel_mu,
       aes(x = Strain, y = rel_mu, fill = Inducers)) +
  geom_bar(stat = 'identity', color = 'black', position = position_dodge()) +
  base_cat_theme +
  facet_wrap(~CRISPRi_system, scales = 'free') +
  facet_theme
