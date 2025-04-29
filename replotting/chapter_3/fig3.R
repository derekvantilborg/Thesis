#### Set working dir and source some functions ####

setwd("/Users/derekvantilborg/Dropbox/PhD/Thesis/replotting/chapter_3")
source('theme.R')

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

figure4 = function(df){
  
  # Select the DL methods + a good and bad classical method
  df = rbind(subset(df, algorithm %in% c('GCN', 'MPNN', 'AFP', "GAT", 'MLP', 'Transformer')),
             subset(df, algorithm == 'SVM' & descriptor == 'ECFP'),
             subset(df, algorithm == 'LSTM' & augmentation == 10),
             subset(df, algorithm == 'CNN' & augmentation == 10))
  
  # order factors of algorithm and descriptor
  algo_order = c('GAT', 'GCN', 'AFP', 'MPNN', 'CNN', 'Transformer', 'LSTM', 'MLP', "SVM")
  df$algorithm = factor(df$algorithm, levels = algo_order)
  df$descriptor = factor(df$descriptor, levels = c('ECFP', 'SMILES', 'Graph'))
  
  colours = descr_cols$cols[match(c('ECFP', 'SMILES', 'Graph'), descr_cols$descr)]
  
  
  ### Box plot ###
  
  box_plot = ggplot(df, aes(x=algorithm, y=cliff_rmse, fill=descriptor))+
    # geom_jitter(aes(color=descriptor), position=position_jitterdodge(0), size=1, shape=1, alpha=0.75) +
    geom_jitter(aes(fill=descriptor), position=position_jitterdodge(0), size=1, shape=21, alpha=0.80, color = "black", stroke = 0.1) +
    # geom_point(alpha = ifelse(df$descriptor == descriptor, 0.80, 0), shape=21, color = "black", size = 1.5, stroke = 0.1)+
    geom_boxplot(alpha=0.5, outlier.size = 0, position = position_dodge(0.75), width = 0.30,
                 outlier.shape=NA, varwidth = FALSE, lwd=0.3, fatten=0.75) +
    
    # geom_jitter(aes(color=descriptor), position=position_jitterdodge(0), size=1, shape=16, alpha=0.75) +
    # geom_boxplot(alpha=0.1, outlier.size = 0, position = position_dodge(0.75), width = 0.25,
    #              outlier.shape=NA, varwidth = FALSE, lwd=0.4, fatten=0.75) +
    scale_y_continuous(breaks = seq(0.25,1.75,0.5), expand = expansion(mult = c(0.01, 0.01)))+
    coord_cartesian(ylim=c(0.25, 1.75))+
    scale_color_manual(values = colours)+
    scale_fill_manual(values = colours)+
    labs(x='Algorithm', y=bquote("RMSE"[cliff]))+
    guides(fill = 'none', 
           color = 'none')+
    default_theme +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
  
  
  ### PCA  ###
  
  # Compute pca coordinates
  pca_all = data_to_biplot(df, val_var="cliff_rmse")
  bi_all = pca_all$bi
  
  bi_all$algorithm = unlist(strsplit(as.character(bi_all$name), ' - '))[2*(1:length(bi_all$name))-1 ]
  bi_all$descriptor = unlist(strsplit(as.character(bi_all$name), ' - '))[2*(1:length(bi_all$name)) ]
  bi_all$algorithm[grepl('CHEMBL', bi_all$name)] = ''
  bi_all$descriptor[grepl('CHEMBL', bi_all$name)] = ''
  
  # rename some stuff
  bi_all$algorithm = gsub('Svm', 'SVM', bi_all$algorithm)
  bi_all$algorithm = gsub('Knn', 'KNN', bi_all$algorithm)
  bi_all$algorithm = gsub('Cnn', 'CNN', bi_all$algorithm)
  bi_all$algorithm = gsub('Mlp', 'MLP', bi_all$algorithm)
  bi_all$algorithm = gsub('Lstm', 'LSTM', bi_all$algorithm)
  bi_all$algorithm = gsub('Afp', 'AFP', bi_all$algorithm)
  bi_all$algorithm = gsub('Gat', 'GAT', bi_all$algorithm)
  bi_all$algorithm =  gsub('Gcn', 'GCN', bi_all$algorithm)
  bi_all$algorithm =  gsub('Mpnn', 'MPNN', bi_all$algorithm)
  bi_all$descriptor =  gsub('ecfp', 'ECFP', bi_all$descriptor)
  bi_all$descriptor =  gsub('graph', 'Graph', bi_all$descriptor)
  bi_all$descriptor =  gsub('smiles', 'SMILES', bi_all$descriptor)
  
  # Define the colours for the descriptors
  bi_all$col = descr_cols$cols[match(bi_all$descriptor, descr_cols$descr)]
  
  # If 'Best' is on the left side of the plot, mirror everything. Makes it easer to compare with the previous plots
  if (subset(bi_all, algorithm == 'Best')$x < 0){
    bi_all$x = bi_all$x*-1
  }
  
  # Get the coordinates of 'best' and 'worst' and find the best axis limit
  best = unlist(subset(bi_all, algorithm == 'Best')[c(2,3)])
  worst = unlist(subset(bi_all, algorithm == 'Worst')[c(2,3)])
  axis_limit = ceiling(max(abs(subset(bi_all, type == 'Score')[c(2,3)]))) 
  
  # Make the actual plot
  pca_plot = ggplot(bi_all, aes(x = x, y =y)) +
    geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.25, size=0.35) +
    geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.25, size=0.35) +
    geom_segment(aes(x = worst[1], y = worst[2], xend = best[1], yend = best[2]), linetype='solid',  alpha = 0.8, size=0.25)+
    # geom_point(aes(x, y ), colour = bi_all$col,  size = 1, shape=1, alpha = ifelse(bi_all$type == 'Score', 0.75, 0)) +
    geom_point(aes(x, y ), fill = bi_all$col, shape = 21,  size = 1.5, alpha = ifelse(bi_all$type == 'Score', 0.8, 0), size = 1.25, stroke = 0.1) +
    geom_text_repel(aes(label = algorithm), colour = bi_all$col,   alpha = ifelse(bi_all$type == 'Score', 1, 0), 
                    size = 2.5, segment.size = 0.25, force = 10, fontface="bold", max.iter = 1505, 
                    max.overlaps = 30, show.legend = FALSE)+
    scale_y_continuous(limits = c(-2.5, 2.5), expand = expansion(mult = c(0.01, 0.01)), breaks = seq(-3,3,1))+
    scale_x_continuous(limits = c(-2.5, 2.5), expand = expansion(mult = c(0.01, 0.01)), breaks = seq(-3,3,1))+
    coord_cartesian(ylim=c(-1.1, 1.1), xlim=c(-1.65, 1.65))+
    labs(x = paste0('PC ',1, ' (',round(pca_all$scree$data$eig[1],1),'%)'),
         y = paste0('PC ',2, ' (',round(pca_all$scree$data$eig[2],1),'%)'),
         shape = 'Algorithm',  color = 'Descriptor') +
    guides(color = guide_legend(override.aes = list(shape = 16), order = 1 )) +
    default_theme + theme(  axis.line.x.bottom=element_blank(),
                            axis.line.y.left=element_blank(),
                            panel.border = element_rect(colour = "#101e25", size = 0.75, fill = NA)) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
  
  
  ### Scatter plot s###
  
  df$descriptor = factor(df$descriptor, levels = c('ECFP', 'SMILES', 'Graph'))
  df$algorithm = factor(df$algorithm, levels = c('GCN', 'LSTM', 'MLP', 'SVM', 'AFP', 'CNN', 'GAT', 'MPNN', 'Transformer')) 
  # algo = 'GAT'
  # descriptor = 'Graph'. 
  # "#578d88"
  scatter_plot = function(algo, descriptor){
    
    # colours = c(rep('#efeae8', length(levels(df$algorithm))))
    # colours[which(levels(df$algorithm) == algo)] = descr_cols$cols[which(descr_cols$descr == descriptor)]
    collor = descr_cols$cols[which(descr_cols$descr == descriptor)]
    p = ggplot(df, aes(x=rmse, y=cliff_rmse, colour=algorithm))+
      geom_point(alpha = ifelse(df$algorithm == algo, 0, 0.80), shape=16, size=1, color='#efeae8' )+
      # geom_point(alpha = ifelse(df$algorithm == algo, 0.5, 0), shape=1, size=1)+
      # geom_point(alpha = ifelse(df$algorithm == algo, 0.75, 0), shape=16, size=1, color=collor)+
      geom_point(alpha = ifelse(df$algorithm == algo, 0.80, 0), shape=21, color = "black", size = 1, stroke = 0.1, fill=collor)+
      geom_abline(slope=1, intercept = 0, linetype='dashed', alpha=0.75, size=0.35)+
      geom_abline(slope=1, intercept = 0.25, linetype='dashed', alpha=0.25, size=0.35)+
      geom_abline(slope=1, intercept = -0.25, linetype='dashed', alpha=0.25, size=0.35)+
      geom_text(x=0.35, y=1.6, label=algo, color = '#101e25', size=2, fontface="bold", hjust=0) +
      scale_x_continuous(breaks = seq(0.25,1.75,0.5), limits = c(0.01,1.75), expand = expansion(mult = c(0.01, 0.01)))+
      scale_y_continuous(breaks = seq(0.25,1.75,0.5), limits = c(0.01,1.75), expand = expansion(mult = c(0.01, 0.01)))+
      labs(x = '',  y = bquote("RMSE"[cliff])) +
      scale_color_manual(values = colours)+
      coord_cartesian(ylim=c(0.25, 1.75), xlim=c(0.25, 1.75))+
      default_theme +
      theme(legend.position = 'none') +
      theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
    
    return(p)
  }
  
  
  gat_relative_plot = scatter_plot('GAT', 'Graph')
  mpnn_relative_plot = scatter_plot('MPNN', 'Graph')
  afp_relative_plot = scatter_plot('AFP', 'Graph')
  gcn_relative_plot = scatter_plot('GCN', 'Graph')
  cnn_relative_plot = scatter_plot('CNN', 'SMILES')
  trans_relative_plot = scatter_plot('Transformer', 'SMILES')
  lstm_relative_plot = scatter_plot('LSTM', 'SMILES')
  mlp_relative_plot= scatter_plot('MLP', 'ECFP')
  
  
  ### combine all plots ###
  box_pca = plot_grid(box_plot, pca_plot, labels = c('a', 'b'), label_size=10, ncol=2, nrow=1, scale=1)
  scatters = plot_grid(gat_relative_plot, gcn_relative_plot, afp_relative_plot, mpnn_relative_plot, 
                       cnn_relative_plot, trans_relative_plot, lstm_relative_plot, mlp_relative_plot, 
                       ncol=4, nrow=2, scale=1)
  fig = plot_grid(box_pca, scatters, ncol=1, nrow =2, scale = 1, rel_heights=c(0.43, 0.58), labels = c('', 'c'), label_size = 10)
  fig = plot_grid(fig, plot_spacer(), ncol=2, rel_widths = c(0.99, 0.01))
  
  return(fig)
  
}



# 130 mm = 5.118 inch
5.118 * 5/7.205

print(figure4(benchmark))
dev.print(pdf, 'Fig_3.pdf', width = 5.118, height = 4) # width = 7.205, height = 4


