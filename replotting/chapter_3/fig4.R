#### Set working dir and source some functions ####

setwd("/Users/derekvantilborg/Dropbox/PhD/Thesis/replotting/chapter_3")
source('theme.R')
library(patchwork)

detach("package:plyr", unload = TRUE)

# GGplot default theme I use
default_theme = theme(
  panel.border = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "plain", size=8, margin = margin(b = 0)),
  axis.text.y = element_text(size=7, face="plain", colour = "#101e25"),
  axis.text.x = element_text(size=7, face="plain", colour = "#101e25"),
  axis.title.x = element_text(size=8, face="plain", colour = "#101e25"),
  axis.title.y = element_text(size=8, face="plain", colour = "#101e25"),
  axis.ticks = element_line(color="#101e25", size=0.35),
  axis.line.x.bottom=element_line(color="#101e25", size=0.35),
  axis.line.y.left=element_line(color="#101e25", size=0.35),
  legend.key = element_blank(),
  legend.position = 'right',
  legend.title = element_text(size=8),
  legend.background = element_blank(),
  legend.text = element_text(size=8),
  legend.spacing.y = unit(0., 'cm'),
  legend.key.size = unit(0.25, 'cm'),
  legend.key.width = unit(0.5, 'cm'),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())


#### Data prep ####

# Import data
benchmark = read_csv('data/MoleculeACE_results.csv')



# Rename some stuff
benchmark$descriptor[benchmark$descriptor == 'TOKENS'] = 'SMILES'
benchmark$descriptor[benchmark$descriptor == 'PHYSCHEM'] = 'Physchem'
benchmark$descriptor[benchmark$descriptor == 'MACCS'] = 'MACCs'
benchmark$descriptor[benchmark$descriptor == 'GRAPH'] = 'Graph'

#### Figures ####
# df = benchmark

figure5 = function(df){
  
  # add dataset abbreviation
  df$label = dataset_abbrv$abbrv[match(df$dataset, dataset_abbrv$dataset)]
  df$label = as.factor(df$label)
  # Calculate RMSE delta
  df$rmse_delta = df$cliff_rmse - df$rmse
  
  # Order datasets
  mean_dataset <- df %>% group_by(label) %>%  dplyr::summarise(MinDataset = cor(cliff_rmse, rmse)) %>% ungroup
  
  dataset_order = mean_dataset$label[order(-mean_dataset$MinDataset)]
  df$label = factor(df$label, levels = dataset_order)
  
  df$corr = mean_dataset$MinDataset[match(df$label, mean_dataset$label)]
  df$corr_label = as.character(round(df$corr,2))
  
  df$descriptor = factor(df$descriptor, levels = c("Graph","ECFP","MACCs","Physchem","SMILES","WHIM"))
  colours = descr_cols$cols[match(levels(df$descriptor), descr_cols$descr)]
  
  # Plot RMSE RMSEcliff scatters
  p_good = ggplot(subset(df, dataset == "CHEMBL214_Ki"), aes(x = rmse, y = cliff_rmse))+
    # geom_point(size=1, shape=1, alpha=0.75, aes(color = descriptor))+
    geom_point(size=1, shape=21, alpha=0.80, aes(fill = descriptor), stroke = 0.1)+
    geom_abline(slope=1, intercept = 0, linetype='dashed', alpha=0.75, size=0.35)+
    geom_abline(slope=1, intercept = 0.25, linetype='dashed', alpha=0.25, size=0.35)+
    geom_abline(slope=1, intercept = -0.25, linetype='dashed', alpha=0.25, size=0.35)+
    labs(x='RMSE', y=bquote("RMSE"[cliff]), fill = 'Descriptor')+
    scale_color_manual(values = colours)+
    scale_fill_manual(values = colours)+
    guides(fill = 'none', color = 'none')+
    scale_x_continuous(breaks = seq(0.25,1.75,0.5), limits = c(0.01,1.75), expand = expansion(mult = c(0.01, 0.01)))+
    scale_y_continuous(breaks = seq(0.25,1.75,0.5), limits = c(0.01,1.75), expand = expansion(mult = c(0.01, 0.01)))+
    coord_cartesian(ylim=c(0.25, 1.75), xlim=c(0.25, 1.75))+
    default_theme + theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
    # theme(
    #   panel.border = element_rect(colour = "#101e25", size = 1, fill = NA),
    #   panel.background = element_blank(),
    #   plot.title = element_text(hjust = 0.5, face = "plain"),
    #   axis.ticks.y = element_line(colour = "#101e25"),
    #   axis.ticks.x = element_line(colour = "#101e25"),
    #   axis.text.y = element_text(size=6, face="plain", colour = "#101e25"),
    #   axis.text.x = element_text(size=6, face="plain", colour = "#101e25"),
    #   axis.title.x = element_text(size=6, face="plain", colour = "#101e25"),
    #   axis.title.y = element_text(size=6, face="plain", colour = "#101e25"),
    #   legend.key = element_blank(),
    #   legend.text = element_text(colour = "#101e25"),
    #   legend.position = 'right', legend.box = "vertical",
    #   legend.title = element_blank(),
    #   legend.background = element_blank(),
    #   panel.grid.major = element_blank(),
    #   panel.grid.minor = element_blank())
  
  p_bad = ggplot(subset(df, dataset == "CHEMBL4203_Ki"), aes(x = rmse, y = cliff_rmse))+
    # geom_point(size=1, shape=1, alpha=0.75, aes(color = descriptor))+
    geom_point(size=1, shape=21, alpha=0.80, aes(fill = descriptor), stroke = 0.1)+
    geom_abline(slope=1, intercept = 0, linetype='dashed', alpha=0.75, size=0.35)+
    geom_abline(slope=1, intercept = 0.25, linetype='dashed', alpha=0.25, size=0.35)+
    geom_abline(slope=1, intercept = -0.25, linetype='dashed', alpha=0.25, size=0.35)+
    labs(x='RMSE', y="", fill = 'Descriptor')+
    scale_color_manual(values = colours)+
    scale_fill_manual(values = colours)+
    guides(fill = 'none', color = 'none')+
    scale_x_continuous(breaks = seq(0.25,1.75,0.5), limits = c(0.01,1.75), expand = expansion(mult = c(0.01, 0.01)))+
    scale_y_continuous(breaks = seq(0.25,1.75,0.5), limits = c(0.01,1.75), expand = expansion(mult = c(0.01, 0.01)))+
    coord_cartesian(ylim=c(0.25, 1.75), xlim=c(0.25, 1.75))+
    default_theme + theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
    # theme(
    #   panel.border = element_rect(colour = "#101e25", size = 1, fill = NA),
    #   panel.background = element_blank(),
    #   plot.title = element_text(hjust = 0.5, face = "plain"),
    #   axis.ticks.y = element_line(colour = "#101e25"),
    #   axis.ticks.x = element_line(colour = "#101e25"),
    #   axis.text.y = element_text(size=6, face="plain", colour = "#101e25"),
    #   axis.text.x = element_text(size=6, face="plain", colour = "#101e25"),
    #   axis.title.x = element_text(size=6, face="plain", colour = "#101e25"),
    #   axis.title.y = element_text(size=6, face="plain", colour = "#101e25"),
    #   legend.key = element_blank(),
    #   legend.text = element_text(colour = "#101e25"),
    #   legend.position = 'right', legend.box = "vertical",
    #   legend.title = element_blank(),
    #   legend.background = element_blank(),
    #   panel.grid.major = element_blank(),
    #   panel.grid.minor = element_blank())
    # 
  
  # Make plot
  diff_per_ds = ggplot(df, aes(x=rmse_delta, y=label, fill = descriptor))+
    geom_vline(xintercept = 0, linetype='dashed', alpha=0.5, colour='#101e25', size=0.35)+
    # stat_boxplot(geom ='errorbar', coef=NULL, aes(group=dataset), lwd=0.25, width = 0.5, size=0.5) +
    geom_boxplot(alpha=0.5, outlier.size = 0, position = position_dodge(0.75), width = 0.5,
                 outlier.shape=NA, varwidth = FALSE, lwd=0.3, fatten=0.75, aes(group=dataset), fill='#97a4ab', colour='#192930') +
    scale_y_discrete(position = "left") +
    scale_x_continuous(labels = c('-0.2', '0', '0.2', '0.4', '0.6')) +
    labs(y='Dataset', x=bquote("RMSE"[cliff]~"- RMSE"), fill = 'Descriptor')+
    guides(fill = 'none', color = 'none')+
    coord_cartesian(xlim=c(-0.2, 0.6))+
    default_theme + theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
  

  corr_per_label = df %>% group_by(label) %>% dplyr::summarize(cor=mean(corr))
  corr_per_label$label = factor(corr_per_label$label, levels = dataset_order)
  corr_per_label$cor_label = as.character(round(corr_per_label$cor,2))
  
  corr_bar = ggplot(corr_per_label, aes(y=label, x=1, fill=cor))+
    geom_tile()+
    geom_text_repel(aes(x = 1, label = cor_label), size = 1, 
                    segment.size = 0.25, force = 0, fontface="plain", max.iter = 10, 
                    max.overlaps = 1000, show.legend = FALSE, color='#101e25')+
    scale_fill_gradient(low = "#577788", high = "#bfbab8") +
    guides(fill = 'none', color = 'none')+
    theme(
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0, 0, 0.75, 0), "cm"))
  
  
  ## Train size effects
  datasets <- read_csv("data/datasets.csv")
  df$train_size = datasets$`Train compounds`[match(df$dataset, datasets$Dataset)]
  corr_per_dataset = df %>% group_by(dataset) %>% dplyr::summarize(cor=cor(rmse, cliff_rmse))
  corr_per_dataset$train_size = datasets$`Train compounds`[match(corr_per_dataset$dataset, datasets$Dataset)]
  
  p_corr = ggplot(corr_per_dataset, aes(y = cor, x = train_size))+
    geom_point(size=1, shape=21, alpha=0.80, fill = '#97a4ab', stroke = 0.1)+
    # geom_text_repel(aes(label = dataset))+
    labs(x="Train molecules", y=bquote("r(RMSE, RMSE"[cliff]~")")) +
    scale_x_continuous(breaks = seq(500,3000,1000), limits = c(450,3000), expand = expansion(mult = c(0.01, 0.01))) +
    default_theme + theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
  
  
  rmse_per_dataset <- group_by(df, dataset)
  rmse_per_dataset = dplyr::summarise(rmse_per_dataset, mean_rmse_delta=mean(rmse_delta), min_rmse_delta=min(rmse_delta),
                               max_rmse_delta=max(rmse_delta), train_size=mean(train_size) )
  
  p_datasize = ggplot(rmse_per_dataset, aes(y = mean_rmse_delta, x = train_size))+
    geom_errorbar(aes(ymin=min_rmse_delta, ymax=max_rmse_delta), 
                  colour='#97a4ab', width=0, size=0.35, alpha=0.5)+
    geom_point(size=1, shape=21, alpha=0.80, fill = '#97a4ab', stroke = 0.1)+
    labs(x="Train molecules", y=bquote("RMSE"[cliff]~"- RMSE")) +
    scale_x_continuous(breaks = seq(500,3000,1000), limits = c(450,3000), expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(ylim=c(-0.2, 0.6))+
    scale_color_manual(values = colours)+
    guides(fill = 'none', color = 'none')+
    default_theme + theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
  
  
  # p_data =plot_grid(p_datasize, p_corr,
  #                   label_size = 8, labels = c('d', 'e'),
  #                   ncol=2, nrow =1, scale = 1)
  
  scatters = plot_grid(p_good, p_bad, p_datasize, p_corr,
                       label_size = 10, labels = c('b', 'c', 'd', 'e'),
                       ncol=2, nrow =2, scale = 1)
  
  plot_corbar = plot_grid(diff_per_ds, corr_bar, ncol=2, nrow =1, scale = 1,
                          rel_widths = c(0.85, 0.15))
  
  scatters = plot_grid(scatters, plot_spacer(), ncol=1, rel_heights = c(0.75, 0.25))
    
  
  p = plot_grid(plot_corbar, scatters,
                       label_size = 10, labels = c('a', ''),
                       ncol=2, nrow =1, scale = 1, rel_widths = c(0.38, 0.62))
  
  p = plot_grid(p, plot_spacer(), ncol=2, rel_widths = c(0.99, 0.01))
  

  return(p)
  
}


new_height = 2.35 + 0.787 

# df = benchmark
# 130 mm = 5.118 inch
5.118 * 5/7.205

# 0.787 inch = 2cm

0.586 * 4

print(figure5(benchmark))
dev.print(pdf, 'Fig_4.pdf', width = 5.118*0.8, height = 2.35 + 0.787) # width = 7.205, height = 4




