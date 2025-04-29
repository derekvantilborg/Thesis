#### Set working dir and source some functions ####

setwd("/Users/derekvantilborg/Dropbox/PhD/Thesis/replotting/chapter_3")
source('theme.R')
library(patchwork)

cols = c('#aa3535','#e54230', '#f76c3c', '#ffb361','#fde0b0', '#72af9c', '#098e91', '#0d677c','#113959', '#16263a')


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
figure3 = function(df){

  # Select only the classical methods
  df = subset(df, algorithm %in% c('RF', 'SVM', 'GBM', 'KNN'))
  
  # order factors of algorithm and descriptor
  mean_algo = df %>% group_by(algorithm) %>% summarise_all("mean")
  mean_descr = df %>% group_by(descriptor) %>% summarise_all("mean")
  
  algo_order = mean_algo[rev(order(mean_algo$cliff_rmse)),]$algorithm
  descr_order = mean_descr[rev(order(mean_descr$cliff_rmse)),]$descriptor
  
  # re-oder factors based on mean performance
  df$algorithm = factor(df$algorithm, levels = algo_order)
  df$descriptor = factor(df$descriptor, levels = descr_order)
  
  
  ### Boxplot ###
  colours = descr_cols$cols[match(c('WHIM', 'Physchem', 'MACCs', 'ECFP'), descr_cols$descr)]
  
  box_plot = ggplot(df, aes(x=algorithm, y=cliff_rmse, fill = descriptor))+
    # geom_jitter(aes(color=descriptor), position=position_jitterdodge(0), size=1, shape=1, alpha=0.5) +
    geom_jitter(aes(fill=descriptor), position=position_jitterdodge(0), size=1, shape=21, alpha=0.80, color = "black", stroke = 0.1) +
    # geom_point(alpha = ifelse(df$descriptor == descriptor, 0.80, 0), shape=21, color = "black", size = 1.5, stroke = 0.1)+
    geom_boxplot(alpha=0.5, outlier.size = 0, position = position_dodge(0.75), width = 0.55,
                 outlier.shape=NA, varwidth = FALSE, lwd=0.3, fatten=0.75) +
    scale_y_continuous(breaks = seq(0.25,1.75,0.5), limits = c(0.25,1.75), expand = expansion(mult = c(0.01, 0.01)))+
    scale_color_manual(values = colours)+
    scale_fill_manual(values = colours)+
    labs(x='', y=bquote("RMSE"[cliff]), fill = 'Descriptor')+
    guides(fill = 'none', color = 'none')+
    default_theme + # theme(  #axis.line.x.bottom=element_blank(),
                            #axis.ticks.x = element_blank(),
                            #axis.text.x=element_text(vjust=15)) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))


  ### PCA ###
  pca_all = data_to_biplot(df, val_var="cliff_rmse" )
  bi_all = pca_all$bi
  
  # Fix some names
  bi_all$algorithm = unlist(strsplit(as.character(bi_all$name),  ' - '))[2*(1:length(bi_all$name))-1 ]
  bi_all$descriptor = unlist(strsplit(as.character(bi_all$name),  ' - '))[2*(1:length(bi_all$name)) ]
  bi_all$algorithm[grepl('CHEMBL', bi_all$name)] = ''
  bi_all$descriptor[grepl('CHEMBL', bi_all$name)] = ''
  
  # rename some stuff
  bi_all$algorithm = gsub('Svm', 'SVM', bi_all$algorithm)
  bi_all$algorithm = gsub('Knn', 'KNN', bi_all$algorithm)
  bi_all$algorithm = gsub('Gbm', 'GBM', bi_all$algorithm)
  bi_all$algorithm = gsub('Rf', 'RF', bi_all$algorithm)
  bi_all$descriptor =  gsub('maccs', 'MACCs', bi_all$descriptor)
  bi_all$descriptor =  gsub('physchem', 'Physchem', bi_all$descriptor)
  bi_all$descriptor =  gsub('whim', 'WHIM', bi_all$descriptor)
  bi_all$descriptor =  gsub('ecfp', 'ECFP', bi_all$descriptor)
  
  # Give colours to the data
  bi_all$col = descr_cols$cols[match(bi_all$descriptor, descr_cols$descr)]
  
  # Get xy the coordinates for the best and worst points
  best = unlist(subset(bi_all, algorithm == 'Best')[c(2,3)])
  worst = unlist(subset(bi_all, algorithm == 'Worst')[c(2,3)])
  
  x_axis_label = paste0('PC1 (',round(pca_all$scree$data$eig[1],1),'%)')
  y_axis_label = paste0('PC2 (',round(pca_all$scree$data$eig[2],1),'%)')
  
  # Make the actual plot
  pca_plot = ggplot(bi_all, aes(x = x, y =y)) +
    geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.25, size=0.35) +
    geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.25, size=0.35) +
    geom_segment(aes(x = worst[1], y = worst[2], xend = best[1], yend = best[2]),
                 linetype='solid',  alpha = 0.8, colour='#101e25', size=0.25)+
    # geom_point(aes(x, y ), colour = bi_all$col, shape = 1,  size = 1, alpha = ifelse(bi_all$type == 'Score', 0.5, 0)) +
    geom_point(aes(x, y ), fill = bi_all$col, shape = 21,  size = 1.5, alpha = ifelse(bi_all$type == 'Score', 0.8, 0), size = 1.25, stroke = 0.1) +
    # geom_point(alpha = ifelse(df$descriptor == descriptor, 0.80, 0), shape=21, color = "black", size = 1.5, stroke = 0.1)+
    geom_text_repel(aes(label = algorithm), colour = bi_all$col, alpha = ifelse(bi_all$type == 'Score', 1, 0), 
                    size = 2.75, segment.size = 0.25, force = 20, fontface="bold", max.iter = 1505, 
                    max.overlaps = 30, show.legend = FALSE)+
    scale_y_continuous(limits = c(-1.35, 1.35), expand = expansion(mult = c(0.01, 0.01)), breaks = seq(-2,2,1))+
    scale_x_continuous(limits = c(-1.35, 1.35), expand = expansion(mult = c(0.01, 0.01)), breaks = seq(-2,2,1))+
    coord_cartesian(ylim=c(-1.1, 1.1))+
    labs(x = x_axis_label, y = y_axis_label, shape = 'Algorithm',  color = 'Descriptor') +
    guides(shape = guide_legend(override.aes = list(color = "#101e25"), order = 2),
           color = guide_legend(override.aes = list(shape = 16), order = 1 )) +
    default_theme + theme(  axis.line.x.bottom=element_blank(),
                            axis.line.y.left=element_blank(),
                            panel.border = element_rect(colour = "#101e25", size = 0.75, fill = NA)) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))

  ### Scatter plots ###
  relative_scatter = function(descriptor, colours){
    p = ggplot(df, aes(x=rmse, y=cliff_rmse, fill=descriptor))+
      # geom_point(alpha = ifelse(df$descriptor == descriptor, 0, 0.80), shape=21, size=1.5, fill='#efeae8' )+
      geom_point(alpha = ifelse(df$descriptor == descriptor, 0, 0.80), shape=16, size = 1, color='#efeae8') +
      # geom_point(alpha = ifelse(df$descriptor == descriptor, 0.5, 0), shape=1, size=1)+
      # geom_point(alpha = ifelse(df$descriptor == descriptor, 0.80, 0), shape=16, size=1.5)+
      geom_point(alpha = ifelse(df$descriptor == descriptor, 0.80, 0), shape=21, color = "black", size = 1, stroke = 0.1)+
      geom_abline(slope=1, intercept = 0, linetype='dashed', alpha=0.75, size=0.35)+
      geom_abline(slope=1, intercept = 0.25, linetype='dashed', alpha=0.25, size=0.35)+
      geom_abline(slope=1, intercept = -0.25, linetype='dashed', alpha=0.25, size=0.35)+
      scale_x_continuous(breaks = seq(0.25,1.75,0.5), limits = c(0.25,1.75), expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_continuous(breaks = seq(0.25,1.75,0.5), limits = c(0.25,1.75), expand = expansion(mult = c(0.01, 0.01))) +
      labs(x = 'RMSE',  y = bquote("RMSE"[cliff])) +
      scale_fill_manual(values = colours)+
      default_theme +
      theme(legend.position = 'none') +
      theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
    
    return(p)
  }
  
  ecfp4_relative_plot = relative_scatter('ECFP', c('#efeae8', '#efeae8', '#efeae8', descr_cols$cols[which(descr_cols$descr == 'ECFP')]))
  maccs_relative_plot = relative_scatter('MACCs', c('#efeae8', '#efeae8', descr_cols$cols[which(descr_cols$descr == 'MACCs')], '#b6b6b6'))
  physchem_relative_plot = relative_scatter('Physchem', c('#efeae8', descr_cols$cols[which(descr_cols$descr == 'Physchem')], '#b6b6b6', '#b6b6b6'))
  whim_relative_plot = relative_scatter('WHIM', c(descr_cols$cols[which(descr_cols$descr == 'WHIM')], '#b6b6b6', '#b6b6b6', '#b6b6b6'))

  ### combine subplots ###
  box_pca = plot_grid(box_plot, pca_plot, labels = c('a', 'b'), label_size = 10, ncol=2, nrow =1, scale = 1)
  scatters = plot_grid(whim_relative_plot, physchem_relative_plot, maccs_relative_plot, ecfp4_relative_plot, ncol=4, nrow =1, scale = 1)
  fig = plot_grid(box_pca, scatters, ncol=1, nrow =2, scale = 1, rel_heights=c(0.6, 0.4), labels = c('', 'c'), label_size = 10)
  
  fig = plot_grid(fig, plot_spacer(), ncol=2, rel_widths = c(0.99, 0.01))

  return(fig)
}


# 130 mm = 5.118 inch
5.118 * 5/7.205

print(figure3(benchmark))
dev.print(pdf, 'Fig_2.pdf', width = 5.118, height = 2.84136) # width = 7.205, height = 4


