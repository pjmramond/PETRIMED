##############################################################################################################
# Run machine-learning training (NNET or XGBoost) based on psbo data (without heterotrophs)
##############################################################################################################

# Load required libraries
library(adiv); library(boot); library(DescTools); library(dendextend)
library(dplyr); library(doFuture); library(doParallel); library(doRNG)
library(foreach); library(forecast); library(future); library(ggdendro)
library(ggplot2); library(ggpubr); library(ggtree); library(ggtreeDendro)
library(ggtreeExtra); library(hydroTSM); library(iNEXT); library(lubridate)
library(maps); library(Metrics); library(microseq); library(ModelMetrics)
library(neuralnet); library(NeuralNetTools); library(nnet); library(parallelDist)
library(patchwork); library(phyloseq); library(phyloseqCompanion); library(phytoclass)
library(progressr); library(readxl); library(reshape2); library(scales)
library(stringr); library(tidyr); library(tidyverse); library(vegan)
library(xgboost); library(tidymodels); library(bonsai); library(purrr); library(caret)
library(ggh4x)

month_labeller <- function(x) {
  m <- format(x, "%b")
  m
}
pal.pig <- c("Chl_a"="seagreen4","Chl_b"="seagreen3","Chl_C2"="#bfe0d0","X19BF"="#cf9e65","X19HF"="#ffd700","Alloxanthin"="#ff4500","Diadinoxanthin"="#ffff00","Fucoxanthin"="#ff8c00","Peridinine"="gold4","Prasinoxanthin"="#fffacd","Zeaxanthin"="#ff0000","bb_Carotene"="#800000")

replace_zero_halfmin <- function(x){
  min_nonzero <- min(x[x > 0], na.rm = TRUE)
  x[x == 0] <- min_nonzero / 2
  return(x)
}

# 1) Metabarcoding
#---------------

### Actual Data
## DATA Import
# PFT relative abundance estimated on 16S and 18$
data.pft.pig<-read.csv("~/Desktop/PETRIMED/DATA/PIGMENTS/sola_pigments_pst_ptot.csv")

# colors for PFTs
pal.pft<-readRDS("~/Desktop/PETRIMED/DATA/PIGMENTS/pst_palette.RData")
colnames(data.pft.pig)[!colnames(data.pft.pig) %in% names(pal.pft)]
pal.pft=pal.pft[names(pal.pft) %in% colnames(data.pft.pig)]

# determine two sets of variables
info=colnames(data.pft.pig)[1:7]
pigments<-colnames(data.pft.pig)[30:ncol(data.pft.pig)]
pfts=names(pal.pft)[-c(1,19)]
proks=pfts[18:length(pfts)]

# Transform datasets
data.pft.pig <- do.call(rbind, lapply(split(data.pft.pig, interaction(data.pft.pig$SF)), function(sf_df) {
  sf_df[pigments] <- lapply(sf_df[pigments], replace_zero_halfmin)
  sf_df[pigments] <- decostand(sf_df[pigments], MARGIN = 1, method = "total")
  sf_df
}))

## Plot Actual Data
mds<-reshape2::melt(data.pft.pig[, c(info, pfts)], id.vars = info)
mds$Date<-as.POSIXct(mds$Date)
mds$Domain<-ifelse(mds$variable %in% proks, "Prokaryotes", "Eukaryotes")

custom_y <- list(
  scale_y_continuous(limits = c(0, 90), breaks = c(0,30,60,90)),
  scale_y_continuous(limits = c(0, 80), breaks = c(0,25,50,75)),
  scale_y_continuous(limits = c(0, 20), breaks = c(0,10,20))
)

gmba <- ggplot(mds[!(mds$SF == "NANO" & mds$Domain == "Prokaryotes"), ], 
               aes(x = Date, y = value, fill = variable)) +
  facet_grid(SF + Domain ~ Year, scale = "free", space = "free") +
  facetted_pos_scales(y = custom_y)+
  geom_area() +
  guides(fill = guide_legend(ncol = 1))+
  #scale_y_continuous(limits = c(0,126), breaks = c(0,25,50,75,100))+
  #geom_text(aes(y = 20, label = "."), size = 4, data = . %>% group_by(ID, SF, Domain, Year, Date) %>% slice(1)) +
  scale_fill_manual(values = pal.pft) +
  scale_x_datetime(date_breaks = "1 month", labels = month_labeller, expand = c(0.025,0)) +
  labs(y = "Observed relative abundance (%)", x="", fill = "PST") +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold", angle = 90),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        axis.ticks.length.x= unit(0.25, "cm"),
        axis.text.x = element_text(size = 14, angle = 90, color = "gray25", hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 14, color = "gray25"));gmba

### Modeled Data
## Model Import
# NANO
mb_nano_models<-readRDS("Desktop/PETRIMED/ANALYSES/ML/FIN/metab_nano_models.RData")
mb_nano_models_info<-readRDS("Desktop/PETRIMED/ANALYSES/ML/FIN/metab_nano_models_info.RData")
mb_nano_models_impo<-readRDS("Desktop/PETRIMED/ANALYSES/ML/FIN/metab_nano_models_importance.RData")
mb_nano_ds<-readRDS("Desktop/PETRIMED/ANALYSES/ML/FIN/metab_nano_ds.RData")

# PICO
mb_pico_models<-readRDS("Desktop/PETRIMED/ANALYSES/ML/FIN/metab_pico_models.RData")
mb_pico_models_info<-readRDS("Desktop/PETRIMED/ANALYSES/ML/FIN/metab_pico_models_info.RData")
mb_pico_models_impo<-readRDS("Desktop/PETRIMED/ANALYSES/ML/FIN/metab_pico_models_importance.RData")
mb_pico_ds<-readRDS("Desktop/PETRIMED/ANALYSES/ML/FIN/metab_pico_ds.RData")
nrow(mb_nano_ds$Micromonas$valid_df)
nrow(mb_pico_ds$Micromonas$valid_df)

# Check number of qualifying models per PFT across a gradient of cor_valid thresholds
thresholds <- seq(0.3, 0.9, by = 0.1)
count_models <- function(models_info, label) {
  do.call(rbind, lapply(thresholds, function(thr) {
    do.call(rbind, lapply(names(models_info), function(pft) {
      x <- models_info[[pft]]
      n <- if (is.data.frame(x)) sum(x$cor_valid >= thr, na.rm = TRUE) else 0
      data.frame(SF = label, PFT = pft, threshold = thr, n_models = n)
    }))
  }))
}
df_counts <- rbind(count_models(mb_nano_models_info, "NANO"),count_models(mb_pico_models_info, "PICO"))

# Quick visual overview
ggplot(df_counts, aes(x = factor(threshold), y = PFT, fill = n_models)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_models), size = 3.5) +
  scale_fill_gradient(low = "white", high = "steelblue4", name = "N models") +
  facet_wrap(~ SF) +
  labs(x = "cor_valid threshold", y = "", title = "Models passing threshold per PFT") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        strip.text = element_text(face = "bold"))

# FILTER NANO
for(i in names(mb_nano_models_info)){mb_nano_models_info[[i]]<-as.data.frame(mb_nano_models_info[[i]][mb_nano_models_info[[i]]$cor_valid >=0.6  ,]) }
for(i in names(mb_nano_models_impo)){mb_nano_models_impo[[i]]<-mb_nano_models_impo[[i]][mb_nano_models_impo[[i]]$Model %in% mb_nano_models_info[[i]]$model_id,] }
for(i in names(mb_nano_models)){mb_nano_models[[i]]<-mb_nano_models[[i]][names(mb_nano_models[[i]]) %in% mb_nano_models_info[[i]]$model_id] }

# FILTER PICO
for(i in names(mb_pico_models_info)){mb_pico_models_info[[i]]<-as.data.frame(mb_pico_models_info[[i]][mb_pico_models_info[[i]]$cor_valid >=0.6  ,]) }
for(i in names(mb_pico_models_impo)){mb_pico_models_impo[[i]]<-mb_pico_models_impo[[i]][mb_pico_models_impo[[i]]$Model %in% mb_pico_models_info[[i]]$model_id,] }
for(i in names(mb_pico_models)){mb_pico_models[[i]]<-mb_pico_models[[i]][names(mb_pico_models[[i]]) %in% mb_pico_models_info[[i]]$model_id] }

# Prediction from nano models
mb_nano_list_pred <- data.pft.pig[,info]
for (n in names(mb_nano_models)) {
  # 2) Predict with HQ models
  p <- NULL;pp<-NULL;m<-length(mb_nano_models[[n]])
  if(m == 0){mb_nano_list_pred[,n]<-NA;mb_nano_ds[[n]]$valid_df[,paste0(n,"_pred")]<-NA;next}
  for (mod in 1:m) {
    pred <- predict(mb_nano_models[[n]][[mod]], newdata=as.matrix(data.pft.pig[,pigments]))
    #pred <- 10^predict(mb_nano_models[[n]][[mod]], newdata=as.matrix(data.pft.pig[,pigments]))
    p <- cbind(p, ifelse(pred <0,0, pred ))
    vpred <- predict(mb_nano_models[[n]][[mod]], newdata=as.matrix(mb_nano_ds[[n]]$valid_df[,pigments]))
    #vpred <- 10^predict(mb_nano_models[[n]][[mod]], newdata=as.matrix(mb_nano_ds[[n]]$valid_df[,pigments]))
    pp <- cbind(pp, ifelse(vpred <0,0, vpred ))
  }
  mb_nano_list_pred[,n] <- rowMeans(p, na.rm=TRUE)
  mb_nano_ds[[n]]$valid_df[,paste0(n,"_pred")]<-rowMeans(pp, na.rm=TRUE)
}
# Prediction from pico models
mb_pico_list_pred <- data.pft.pig[,info]
for (n in names(mb_pico_models)) {
  # 2) Predict with HQ models
  p <- NULL;pp<-NULL;m<-length(mb_pico_models[[n]])
  if(m == 0){mb_pico_list_pred[,n]<-NA;mb_pico_ds[[n]]$valid_df[,paste0(n,"_pred")]<-NA;next}
  for (mod in 1:m) {
    pred <- predict(mb_pico_models[[n]][[mod]], newdata=as.matrix(data.pft.pig[,pigments]))
    #pred <- 10^predict(mb_pico_models[[n]][[mod]], newdata=as.matrix(data.pft.pig[,pigments]))
    p <- cbind(p, ifelse(pred <0,0, pred ))
    vpred <- predict(mb_pico_models[[n]][[mod]], newdata=as.matrix(mb_pico_ds[[n]]$valid_df[,pigments]))
    #vpred <- 10^predict(mb_pico_models[[n]][[mod]], newdata=as.matrix(mb_pico_ds[[n]]$valid_df[,pigments]))
    pp <- cbind(pp, ifelse(vpred <0,0, vpred ))
  }
  mb_pico_list_pred[,n] <- rowMeans(p, na.rm=TRUE)
  mb_pico_ds[[n]]$valid_df[,paste0(n,"_pred")]<-rowMeans(pp, na.rm=TRUE)
}
# merge prediction
mdsp<-as.data.frame(rbind(mb_nano_list_pred[mb_nano_list_pred$SF %in% "NANO",],mb_pico_list_pred[mb_pico_list_pred$SF %in% "PICO",]))
mdsp<-reshape2::melt(mdsp, id.vars = info)
mdsp$Date<-as.POSIXct(mdsp$Date)
mdsp$Domain<-ifelse(mdsp$variable %in% proks, "Prokaryotes", "Eukaryotes")

# a) Make prediction of the patterns of the full dataset
mdsp[mdsp$variable %in% c("Coscinodiscophyceae"),]
gmbp <- ggplot(mdsp[!(mdsp$SF == "NANO" & mdsp$Domain == "Prokaryotes") & mdsp$variable %in% pfts,], 
               aes(x = Date, y = value, fill = variable)) +
  facet_grid(SF + Domain ~ Year, scale = "free", space = "free") +
  facetted_pos_scales(y = custom_y)+
  geom_area() +
  #geom_text(aes(y = 20, label = "."), size = 4, data = . %>% group_by(ID, SF, Domain, Year, Date) %>% slice(1)) +
  guides(fill = guide_legend(ncol = 1))+
  scale_fill_manual(values = pal.pft) +
  scale_x_datetime(date_breaks = "1 month", labels = month_labeller, expand = c(0.025,0)) +
  labs(y = "Predicted relative abundance (%)", x="", fill = "PST") +
  scale_y_continuous(limits = c(0,126), breaks = c(0,25,50,75,100))+
  theme_void(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold", angle = 90),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"),
        axis.ticks.length.x= unit(0.25, "cm"),
        axis.text.x = element_text(size = 14, angle = 90, color = "gray25", hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 14, color = "gray25"));gmbp

# assemble patterns and predictions
first<-gmba+theme(strip.text.y = element_blank()) + gmbp + plot_layout(guides = "collect" )

# b) Make correlation plots 
gg_cor_nano <- list()
metrics_nano <- data.frame(PFT = pfts, cor = NA, MAE = NA, RMSE = NA, avg.ab = NA, nb.q.model = NA)
for (i in seq_along(pfts)) {
  pft <- pfts[i]
  df  <- mb_nano_ds[[pft]]$valid_df
  metrics_nano$avg.ab[i] <- mean(bind_rows(mb_nano_ds[[pft]])[,pft])
  metrics_nano$nb.q.model[i] <- length(mb_nano_models[[pft]])
  if(nrow(na.omit(df)) == 0){next}
  obs  <- df[[pft]]
  pred <- df[[paste0(pft, "_pred")]]
  ## metrics
  metrics_nano$cor[i]  <- cor(obs, pred, use = "complete.obs")
  metrics_nano$MAE[i]  <- Metrics::mae(obs, pred)
  metrics_nano$RMSE[i] <- Metrics::rmse(obs, pred)
  ## plot
  gg_cor_nano[[pft]] <-ggplot(df, aes_string(x = pft, y = paste0(pft, "_pred"))) +
    geom_point(col = "gray1", alpha = 0.75, size = 3) +
    geom_abline(size = 0.5, col = "gray5" ,linetype = "dashed")+
    geom_smooth(method = "lm", size = 1, col = "coral3", se = FALSE) +
    stat_cor(label.x.npc = "left", label.y.npc = "top") +
    labs(title = pft,x = "Observed relative abundance (%)",y = "Predicted relative abundance (%)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text  = element_text(size = 14))
}
gg_cor_pico <- list()
metrics_pico <- data.frame(PFT = pfts, cor = NA, MAE = NA, RMSE = NA, avg.ab = NA, nb.q.model = NA)
for (i in seq_along(pfts)) {
  pft <- pfts[i]
  df  <- mb_pico_ds[[pft]]$valid_df
  metrics_pico$avg.ab[i] <- mean(bind_rows(mb_pico_ds[[pft]])[,pft])
  metrics_pico$nb.q.model[i] <- length(mb_pico_models[[pft]])
  if(nrow(na.omit(df)) == 0){next}
  obs  <- df[[pft]]
  pred <- df[[paste0(pft, "_pred")]]
  ## metrics
  metrics_pico$cor[i]  <- cor(obs, pred, use = "complete.obs")
  metrics_pico$MAE[i]  <- Metrics::mae(obs, pred)
  metrics_pico$RMSE[i] <- Metrics::rmse(obs, pred)
  ## plot
  gg_cor_pico[[pft]] <-ggplot(df, aes_string(x = pft, y = paste0(pft, "_pred"))) +
    geom_point(col = "gray1", alpha = 0.75, size = 3) +
    geom_abline(size = 0.5, col = "gray5" ,linetype = "dashed")+
    geom_smooth(method = "lm", size = 1, col = "coral3", se = FALSE) +
    stat_cor(label.x.npc = "left", label.y.npc = "top") +
    labs(title = pft,x = "Observed relative abundance (%)",y = "Predicted relative abundance (%)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0, size = 14, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text  = element_text(size = 14))
}
metrics_nano
metrics_pico

# c) Make Importance plots 
ggolden_nano<-list()
for (p in names(mb_nano_models_impo)){
  # plot importance per model for all quality models
  ggolden_nano[[p]]<-ggplot(mb_nano_models_impo[[p]], aes(x = Feature, y= Model, fill = Gain ))+
    geom_tile()+
    scale_fill_gradientn(colours = c("gray85", "#011945"), limits = c(0,1), name = "Relative\nimportance", na.value = "transparent")+
    ggtitle(p)+
    ylab(paste("Models (",length(unique(mb_nano_models_impo[[p]]$Model)), ")", sep = "" ))+
    xlab("")+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    theme(axis.title = element_text(size = 14,hjust = 0.5),
          plot.title = element_text(size = 14, face = "bold",hjust = 0),
          plot.subtitle  = element_text(size = 14,hjust = 0),
          panel.spacing = unit(1, "lines"),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 14,  angle = 60, hjust = 1, vjust = 1),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          panel.grid = element_blank(),
          legend.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = "transparent", linewidth = 2),
          axis.ticks = element_line(colour = "gray5"),
          plot.background = element_rect(fill = "transparent",colour = NA))
}
ggolden_pico<-list()
for (p in names(mb_pico_models_impo)){
  # plot importance per model for all quality models
  ggolden_pico[[p]]<-ggplot(mb_pico_models_impo[[p]], aes(x = Feature, y= Model, fill = Gain ))+
    geom_tile()+
    scale_fill_gradientn(colours = c("gray85", "#011945"), limits = c(0,1), name = "Relative\nimportance", na.value = "transparent")+
    ggtitle(p)+
    ylab(paste("Models (",length(unique(mb_pico_models_impo[[p]]$Model)), ")", sep = "" ))+
    xlab("")+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    theme(axis.title = element_text(size = 14,hjust = 0.5),
          plot.title = element_text(size = 14, face = "bold",hjust = 0),
          panel.spacing = unit(1, "lines"),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 14,  angle = 60, hjust = 1, vjust = 1),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          panel.grid = element_blank(),
          legend.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = "transparent", linewidth = 2),
          axis.ticks = element_line(colour = "gray5"),
          plot.background = element_rect(fill = "transparent",colour = NA))
}

# plot the olden importance of pigments for each model predicting a pft
ggolden_nano$Micromonas
ggolden_pico$Micromonas

#gain<-ggpubr::ggarrange(plotlist = ggolden_nano, nrow = ceiling(length(ggolden_nano)/5), ncol = 5,align = "hv", common.legend = TRUE,legend = "right") #15x10
#annotate_figure(gain,bottom = text_grob("Pigments",hjust = 0.5, face = "bold", size = 16))
sec<-gg_cor_nano$Dinophyceae + 
  gg_cor_nano$Bacillariophyceae +
  gg_cor_nano$Haptophyta +
  gg_cor_pico$Bathycoccus +
  gg_cor_pico$Synechococcus + plot_layout(ncol = 2, byrow = FALSE, axis_titles = "collect")

third<-ggolden_nano$Dinophyceae + 
  ggolden_nano$Bacillariophyceae +
  ggolden_nano$Haptophyta +
  ggolden_pico$Bathycoccus +
  ggolden_pico$Synechococcus + plot_layout(ncol = 2, byrow = FALSE,guides = "collect")

des<-"AA
      AA
      BC"

wrap_plots(first, sec, third, design = des)


#ggpubr::ggarrange(plotlist = ggolden_nano, nrow = ceiling(length(ggolden_nano)/5), ncol = 5,align = "hv", common.legend = TRUE,legend = "right") #15x10
#wrap_plots(sec, third, ncol = 2)

metrics_nano
metrics_pico

#plot(metrics_nano$MAE, metrics_nano$avg.ab)
cor.test(metrics_nano$RMSE, metrics_nano$avg.ab, use = "complete.obs", method = "pearson")

metrics_nano[order(metrics_nano$cor),]
metrics_pico[order(metrics_pico$cor),]
summary(rbind(metrics_nano,metrics_pico))

left<-names(ggolden_nano)[!names(ggolden_nano) %in% names(gg_cor_nano)]

on<-ggarrange(plotlist = ggolden_nano[-match(left, names(ggolden_nano))], ncol = 5, nrow = 4)
op<-ggarrange(plotlist = ggolden_pico[-match(left[c(1,3)], names(ggolden_pico))], ncol = 5, nrow = 4)

cn<-ggarrange(plotlist = gg_cor_nano, ncol = 5, nrow = 4)
cp<-ggarrange(plotlist = gg_cor_pico, ncol = 5, nrow = 4)


wrap_plots(on,cn, ncol = 1) +plot_annotation(tag_levels = "A")
wrap_plots(op,cp, ncol = 1) +plot_annotation(tag_levels = "A")

