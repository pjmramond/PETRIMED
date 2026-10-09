################################################################################################
# Script to analyze the results of simka analysis in SOLA/MOLA metagenomics data
################################################################################################

# Libraries
library(ggplot2)
library(stringr)
library(reshape2)
library(dplyr)
library(vegan)
library(ggtree)
library(ggtreeExtra)
library(dendextend)
library(patchwork)

# 1) Import data
# ---------------

# SOLA: PICO 
pico.mat<-read.csv("~/Desktop/PETRIMED/DATA/METAG/simka/kmer70_all/all_mat_abundance_braycurtis.csv", sep = ";", row.names = 1)
pico.k<-read.table("~/Desktop/PETRIMED/DATA/METAG/simka/clusters/sample_cluster");colnames(pico.k)<-c("sample", "cluster")

# SOLA: NANO
nano.mat<-read.csv("~/Desktop/PETRIMED/DATA/METAG/simka/kmer70_n/mat_abundance_braycurtis.csv", sep = ";", row.names = 1)
nano.k<-read.table("~/Desktop/PETRIMED/DATA/METAG/simka/clusters/nano_sample_cluster");colnames(nano.k)<-c("sample", "cluster")

# MOLA: PICO
mola.mat<-read.csv("~/Desktop/PETRIMED/DATA/METAG/simka/kmer70_mola/mat_abundance_braycurtis.csv", sep = ";", row.names = 1)
mola.k<-read.table("~/Desktop/PETRIMED/DATA/METAG/simka/clusters/mola_sample_cluster");colnames(mola.k)<-c("sample", "cluster")

# Environmental Data
clusters<-rbind(pico.k, nano.k, mola.k)
clusters$cluster<-factor(clusters$cluster, levels = as.character(1:13), ordered = TRUE)
info<-read.csv("~/Desktop/PETRIMED/DATA/METAG/metadata/metag_env.csv")
info<-info[grep("NANO|PICO",info$ID),]
info$ID<-gsub("PICO_SO","SO",info$ID)
info<-merge(info, clusters,by.x = "ID", by.y = "sample")
#info<-unique(info)
info$Year<-format(as.Date(info$Date), "%Y")
info$Season<-factor(info$Season, levels = c("Winter","Spring", "Summer", "Autumn"), ordered = TRUE)
info$Depth<-factor(info$Depth, levels = c("3","5","40", "150", "500"), ordered = TRUE)

# 2) Multivariate analysis: NMDS and Dendrograms
# -------------------------------------------

# Random color palette for clusters
pk=c("#ebbcb0", "green4", "black", "coral3", "coral","#d741a7","#3a1772","#5398be","#f2cd5d","#dea54b", "#a7ba4e", "red", "steelblue1")

## a) SOLA: PICO
# -------------------------------------------

# NMDS
nm.ps<-metaMDS(pico.mat, k = 4)
st<-nm.ps$stress
nm.ps<-merge(scores(nm.ps)$sites, info, by.x = "row.names", by.y = "ID")

# Seasons
ggplot(data = nm.ps, aes(x =NMDS1, y = NMDS2, colour = Season, shape =))+
  geom_point(size = 3, alpha = 0.8)+
  scale_color_manual(values = c("Spring" = "green", "Summer" = "gold3", "Autumn" = "coral3", "Winter" = "steelblue4") )+
  annotate("label",  x=Inf, y = Inf, label = paste0("Stress:", round(st,2)), vjust=1.5, hjust=1.5)+
  theme(panel.background = element_rect(fill = "#edf2f7"),aspect.ratio = 1)

# clusters
gc<-ggplot(data = nm.ps, aes(x =NMDS1, y = NMDS2, colour = cluster))+
  geom_point(size = 3, alpha = 0.8)+
  scale_color_manual(values = pk, name = "Cluster - k")+
  annotate("label",  x=Inf, y = Inf, label = paste0("Stress:", round(st,2)), vjust=1.5, hjust=1.5)+
  theme(panel.background = element_rect(fill = "#edf2f7"),aspect.ratio = 1);gc

# Hierarchical clustering
hc <- hclust(as.dist(pico.mat), method = "ward.D")
cutree(hc, h = 0.91)

gd<-ggtree(hc, layout = "dendrogram") %<+% 
  pico.k +
  geom_vline(xintercept =  -0.91, linetype = "dotted") +
  geom_tippoint(aes(color = factor(cluster)), size = 2) +
  scale_color_manual(values = pk, guide = "none")+
  labs(title = "SOLA - PICO / simka (kmer=70)\nBray-Curtis Distance - Ward's Clustering Method")+
  theme(plot.title = element_text(size = 12, hjust = 0.5));gd
gd+gc    

gs<-ggplot(info[grep("PICO|NANO",info$ID, invert = TRUE),], aes(x = as.POSIXct(Date), y = 1, fill = cluster) )+
  geom_area(stat = "identity", position = "stack")+
  geom_text(aes(y = 0.5, label = ".") )+
  scale_fill_manual(values = pk)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_datetime()+
  theme_void()+
  theme(axis.text.x = element_text(size = 12) );gs

gd+gs + plot_annotation(tag_levels = "A")

## b) SOLA: NANO
# -------------------------------------------

# NMDS
nm.ps<-metaMDS(nano.mat, k = 2)
st<-nm.ps$stress
nm.ps<-merge(scores(nm.ps), info, by.x = "row.names", by.y = "ID")

# Seasons
gs<-ggplot(data = nm.ps, aes(x =NMDS1, y = NMDS2, colour = Season, shape = Year))+
  geom_point(size = 3, alpha = 0.8)+
  scale_color_manual(values = c("Spring" = "green", "Summer" = "gold3", "Autumn" = "coral3", "Winter" = "steelblue4") )+
  annotate("label",  x=Inf, y = Inf, label = paste0("Stress: ", round(st,2)), vjust=1.5, hjust=1.5)+
  theme(panel.background = element_rect(fill = "#edf2f7"),aspect.ratio = 1);gs

# clusters
gc<-ggplot(data = nm.ps, aes(x =NMDS1, y = NMDS2, colour = cluster))+
  geom_point(size = 3, alpha = 0.8)+
  scale_color_manual(values = pk, name = "Cluster - k")+
  annotate("label",  x=Inf, y = Inf, label = paste0("Stress: ", round(st,2)), vjust=1.5, hjust=1.5)+
  theme(panel.background = element_rect(fill = "#edf2f7"),aspect.ratio = 1);gc

library(patchwork)
# Hierarchical clustering
hc <- hclust(as.dist(nano.mat), method = "ward.D")
gd<-ggtree(hc, layout = "dendrogram") %<+% 
  nano.k +
  geom_vline(xintercept =  -0.91, linetype = "dotted") +
  geom_tippoint(aes(color = factor(cluster)), size = 2) +
  scale_color_manual(values = pk, guide = "none")+
  labs(title = "SOLA - NANO / simka (kmer=70)\nBray-Curtis Distance - Ward's Clustering Method")+
  theme(plot.title = element_text(size = 12, hjust = 0.5));gd
gd+gc+gs+ plot_annotation(tag_levels = "A")

ggplot(info[grep("NANO",info$ID, invert = FALSE),], aes(x = as.POSIXct(Date), y = 1, fill = cluster) )+
  geom_area(stat = "identity", position = "stack")+
  geom_text(aes(y = 0.5, label = ".") )+
  scale_fill_manual(values = pk)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_datetime()+
  theme_void()+
  theme(axis.text.x = element_text(size = 12) )

## c) MOLA: PICO
# -------------------------------------------

# NMDS
nm.ps<-metaMDS(mola.mat, k = 2)
st<-nm.ps$stress
nm.ps<-merge(scores(nm.ps), info, by.x = "row.names", by.y = "ID")

# Seasons
gs<-ggplot(data = nm.ps, aes(x =NMDS1, y = NMDS2, colour = Season, shape = Year))+
  geom_point(size = 3, alpha = 0.8)+
  scale_color_manual(values = c("Spring" = "green", "Summer" = "gold3", "Autumn" = "coral3", "Winter" = "steelblue4") )+
  scale_shape_manual(values = c(3,4,8,15,16,17,18))+
  annotate("label",  x=Inf, y = Inf, label = paste0("Stress:", round(st,2)), vjust=1.5, hjust=1.5)+
  theme(panel.background = element_rect(fill = "#edf2f7"),aspect.ratio = 1);gs

# Depth
gdp<-ggplot(data = nm.ps, aes(x =NMDS1, y = NMDS2, colour = Depth))+
  geom_point(size = 3, alpha = 0.8)+
  scale_color_manual(values = c("5" = "#a7ba4e", "40" = "#8db8e0", "150" = "#145591", "500"="#010f1c"))+
  annotate("label",  x=Inf, y = Inf, label = paste0("Stress:", round(st,2)), vjust=1.5, hjust=1.5)+
  theme(panel.background = element_rect(fill = "#edf2f7"),aspect.ratio = 1);gdp

# clusters
gc<-ggplot(data = nm.ps, aes(x =NMDS1, y = NMDS2, colour = cluster))+
  geom_point(size = 3, alpha = 0.8)+
  scale_color_manual(values = pk, name = "Cluster - k")+
  annotate("label",  x=Inf, y = Inf, label = paste0("Stress:", round(st,2)), vjust=1.5, hjust=1.5)+
  theme(panel.background = element_rect(fill = "#edf2f7"),aspect.ratio = 1);gc

# Hierarchical clustering
hc <- hclust(as.dist(mola.mat), method = "ward.D")
gd<-ggtree(hc, layout = "dendrogram") %<+% 
  mola.k +
  geom_vline(xintercept =  -0.9, linetype = "dotted") +
  geom_tippoint(aes(color = factor(cluster)), size = 2) +
  scale_color_manual(values = pk, guide = "none")+
  labs(title = "MOLA - PICO / simka (kmer=70)\nBray-Curtis Distance - Ward's Clustering Method")+
  theme(plot.title = element_text(size = 12, hjust = 0.5));gd

wrap_plots(gd, gc, gdp, gs)+ plot_annotation(tag_levels = "A", )

ggplot(unique(info[grep("MO",info$ID, invert = FALSE),-match("ID", colnames(info))]), aes(x = as.POSIXct(Date), y = 1, fill = cluster) )+
  facet_grid(Depth~., switch = "y")+
  geom_area(stat = "identity")+
  geom_text(aes(y = 0.5, label = ".") )+
  scale_fill_manual(values = pk)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_datetime()+
  theme_void()+
  theme(axis.text.x = element_text(size = 12), strip.text.x = element_text(size = 12))



