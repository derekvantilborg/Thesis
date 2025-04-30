
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
umap_cols = c("#dddddd", "#0d677c", "#ffb361")
custom_colours = c("#0d677c", "#72af9c", "#9dc6c6", "#ffb361", "#bbbbbb", "#efc57b")

col_gradient = c("#7b5435", "#c28a5e", "#efc57b", "#85b8d6", "#5c8095", "#3e5a6a")
col_gradient = c("#aa3535", "#ffb361", "#fde0b0", "#72af9c", "#098e91", "#0d677c")

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

df5 <- read_csv("data/exfig1.csv")

df5$n_start = factor(df5$n_start, levels=sort(unique(df5$n_start)))


df5 = df5 %>%
  group_by(acquisition_method, n_start, dataset, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()


fig5a = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) +
  labs(y = 'Increase in enrichment\nover similarity search', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_theme + theme(plot.margin = unit(c(0.25, 0, 0.5, 0.2), "cm"), legend.position = 'None')

fig5b = ggplot(subset(df5, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) + 
  labs(y = '', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.5, 0.25), "cm"), legend.position = 'None')


fig5c = ggplot(subset(df5, dataset == 'VDR' & architecture == 'mlp'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) + 
  labs(y = '', x='', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  default_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.5, 0.25), "cm"), legend.position = 'None')




fig5d = ggplot(subset(df5, dataset == 'PKM2' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) +
  labs(y = 'Increase in enrichment\nover similarity search', x='active learning cycle', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_theme + theme(plot.margin = unit(c(0.25, 0, 0.5, 0.2), "cm"), legend.position = 'None')


fig5e = ggplot(subset(df5, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) + 
  labs(y = '', x='active learning cycle', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.5, 0.25), "cm"), legend.position = 'None')


fig5f = ggplot(subset(df5, dataset == 'VDR' & architecture == 'gcn'), aes(x = train_cycle, y=enrichment_mean, color=n_start, linetype=architecture, fill=n_start))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se), 
              color=NA, alpha=0.1) +
  geom_line(size=0.45) + 
  labs(y = '', x='active learning cycle', color = 'Start size', fill = 'Method')+
  scale_color_manual(values=col_gradient) +
  scale_fill_manual(values=col_gradient) +
  coord_cartesian(xlim = c(0, 15.1), ylim = c(-2, 4), expand=F) +
  scale_y_continuous(breaks = seq(-2, 4, by=1)) +
  scale_linetype_manual(values=c("dashed"))+
  default_theme + theme(plot.margin = unit(c(0.25, 0.25, 0.5, 0.25), "cm"), legend.position = 'None')

fig5 = plot_grid(fig5a, fig5b, fig5c, plot_spacer(),
                 fig5d, fig5e, fig5f, plot_spacer(),
                 labels = c('a', 'b', 'c', '', 
                            'd', 'e', 'f', ''), 
                 ncol=4, label_size=10)



pdf('fig5.pdf', width = 5.118, height = 2.7)
print(fig5)
dev.off()
fig5

