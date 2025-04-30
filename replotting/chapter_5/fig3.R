
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

custom_colours = c("#0d677c", "#72af9c", "#9dc6c6", "#ffb361", "#bbbbbb", "#efc57b")


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

df3 <- read_csv('data/fig3.csv')


###### plots ######

df3 = df3 %>%
  group_by(acquisition_method, bias, dataset, total_mols_screened, train_cycle, architecture) %>%
  summarise(across(c("hits_discovered", "test_tpr", 'test_roc_auc', 
                     'test_balanced_accuracy', 'enrichment'), 
                   list(mean = mean, sd = sd, se = se))) %>% ungroup()

fig3a = ggplot(subset(df3, dataset == 'PKM2' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  default_theme + theme(plot.margin = unit(c(0.25, 0, 0, 0.25), "cm"),
                        legend.position = 'None')
# u r b l
fig3b = ggplot(subset(df3, dataset == 'PKM2' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed"))+
  default_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"),
                        legend.position = 'None')

fig3c = ggplot(subset(df3, dataset == 'PKM2' & architecture == 'rf'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dotted"))+
  default_theme + theme(plot.margin = unit(c(0.25, 0.25, 0, 0.25), "cm"),
                        legend.position = 'None')

fig3d = ggplot(subset(df3, dataset == 'ALDH1' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  default_theme + theme(plot.margin = unit(c(0, 0, 0.25, 0.25), "cm"),
                        legend.position = 'None')

fig3e = ggplot(subset(df3, dataset == 'ALDH1' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed"))+
  default_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"),
                        legend.position = 'None')

fig3f = ggplot(subset(df3, dataset == 'ALDH1' & architecture == 'rf'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dotted"))+
  default_theme + theme(plot.margin = unit(c(0, 0.25, 0.25, 0.25), "cm"),
                        legend.position = 'None')

fig3g = ggplot(subset(df3, dataset == 'VDR' & architecture == 'mlp'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = 'Mean enrichment in\nacquired molecules', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  default_theme + theme(plot.margin = unit(c(-0.25, 0, 0.5, 0.25), "cm"),
                        legend.position = 'None')

fig3h = ggplot(subset(df3, dataset == 'VDR' & architecture == 'gcn'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dashed"))+
  default_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"),
                        legend.position = 'None')

fig3i = ggplot(subset(df3, dataset == 'VDR' & architecture == 'rf'), aes(x = total_mols_screened, y=enrichment_mean, color=acquisition_method, fill=acquisition_method, linetype=architecture))+
  geom_ribbon(aes(ymin = enrichment_mean - enrichment_se, ymax = enrichment_mean + enrichment_se),
              color=NA, alpha=0.1) +
  labs(y = '', x='n molecules screened', color = 'Method', fill = 'Method')+
  scale_color_manual(values=custom_colours) +
  scale_fill_manual(values=custom_colours) +
  geom_line(size=0.45) + 
  coord_cartesian(xlim = c(64, 1010), ylim = c(0, 7), expand=F) + 
  scale_y_continuous(breaks = seq(0,7, by=1)) +
  scale_x_continuous(breaks = c(64, 250, 500, 750, 1000)) +
  scale_linetype_manual(values=c("dotted"))+
  default_theme + theme(plot.margin = unit(c(-0.25, 0.25, 0.5, 0.25), "cm"),
                        legend.position = 'None')

# "#005f73" "#6f9cbc" "#6d7889" "#ee9b00" "#bbbbbb" "#f8cd48"


fig3 = plot_grid(fig3a, fig3b, fig3c, plot_spacer(),
                 fig3d, fig3e, fig3f, plot_spacer(),
                 fig3g, fig3h, fig3i, plot_spacer(),
                 labels = c('a', 'b', 'c', '', 
                            'd', 'e', 'f', '', 
                            'g', 'h' ,'i', ''), 
                 ncol=4, label_size=10)


pdf("fig3.pdf", width = 5.118, height = 3.6)
print(fig3)
dev.off()
fig3



