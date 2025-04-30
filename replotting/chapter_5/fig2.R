
library(readr)
library(dplyr)
library(data.table)
library(readr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(data.table)
library(patchwork)
library(scico)


se <- function(x, na.rm = FALSE) {sd(x, na.rm=na.rm) / sqrt(sum(1*(!is.na(x))))}
zip <- function(...) {mapply(list, ..., SIMPLIFY = FALSE)}


'#0d677c'
'#83b0b7'
'#9dc6c6'

cols = c('#aa3535','#e54230','#f76c3c','#ffb361','#fde0b0', '#72af9c', '#098e91', '#0d677c','#113959', '#16263a')

bias_colours = c("#0d677c", "#efc57b", "#aa3535")


acq_funcs = c("BALD (least mutual information)", 
              "Exploit (best predictions)", 
              "Exploit, no retraining", 
              "Explore (most uncertain)", 
              "Random",
              "Similarity")

acq_cols = list(custom_colours, acq_funcs)


# GGplot default theme I use
default_theme = theme(
  panel.border = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, face = "plain", size=8, margin = margin(b = 0)),
  axis.text.y = element_text(size=7, face="plain", colour = "#101e25"),
  axis.text.x = element_text(size=7, face="plain", colour = "#101e25"),
  axis.title.x = element_text(size=8, face="plain", colour = "#101e25"),
  axis.title.y = element_text(size=8, face="plain", colour = "#101e25"),
  axis.ticks.x = element_line(color="#101e25", size=0.35),
  axis.ticks.y = element_line(color="#101e25", size=0.35),
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

setwd("/Users/derekvantilborg/Dropbox/PhD/Thesis/replotting/chapter_5")

df2 <- read_csv('data/fig2.csv')


df2 = df2 %>%
  group_by(acquisition_method, bias, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "mean_total_sims", 'mean_total_hits_sims', 
                     'n_unique_scaffolds', 'enrichment', 'mean_tani_per_batch',
                     'mean_tani_batch_to_start_batch', 'mean_tani_all_mols_to_start_batch'), 
                   list(mean = mean, se = se), na.rm = TRUE)) %>% ungroup()

df2$architecture = factor(df2$architecture, levels=c('mlp', 'gcn'))
df2$bias = factor(df2$bias, levels=c('Innate', 'Small', 'Large'))


fig2a = ggplot(subset(df2, acquisition_method %in% c("Exploit (best predictions)")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = 'Similarity of acquired\nbatch to start data', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.5, by=0.1)) +
  default_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0), "cm"),
                        legend.position = 'None',)

fig2b = ggplot(subset(df2, acquisition_method %in% c("BALD (least mutual information)")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = '', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.5, by=0.1)) +
  default_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.1), "cm"),
                        legend.position = 'None',)

fig2c = ggplot(subset(df2, acquisition_method %in% c("Exploit, no retraining")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = '', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.5, by=0.1)) +
  default_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.1), "cm"),
                        legend.position = 'None',)

fig2d = ggplot(subset(df2, acquisition_method %in% c("Similarity")), 
               aes(x = train_cycle, y=mean_tani_batch_to_start_batch_mean, 
                   color=bias, fill=bias, group_by=bias, linetype=architecture))+
  geom_ribbon(aes(ymin = mean_tani_batch_to_start_batch_mean - mean_tani_batch_to_start_batch_se, ymax = mean_tani_batch_to_start_batch_mean + mean_tani_batch_to_start_batch_se), color=NA, alpha=0.1) +
  labs(y = '', x='Active learning cycle', color = 'Method', fill = 'Method')+
  scale_color_manual(values=bias_colours) +
  scale_fill_manual(values=bias_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(0, 15.5), ylim = c(0.10, 0.5), expand=F) +
  scale_y_continuous(breaks = seq(0.1, 0.5, by=0.1)) +
  default_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.1), "cm"),
                        legend.position = 'None',)



fig2 = plot_grid(fig2a, fig2b, fig2c, fig2d,
                 ncol=4, 
                 label_size = 10, labels= c('a', 'b', 'c', 'd'))

pdf("fig2.pdf", width = 5.118, height = 1.15)
print(fig2)
dev.off()
fig2
