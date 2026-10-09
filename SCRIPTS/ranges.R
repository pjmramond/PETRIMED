# Script to make additional figures for the PETRIMED Pigment paper
# ------------------------------------------------------------------

# Load packages
library(patchwork)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(vegan)
library(lubridate)
library(tidyverse)
library(scales)

# create objects
pal.pig <- c(
  "Chl_a"="seagreen4",
  "Chl_b"="seagreen3",
  "Chl_C2"="#bfe0d0",
  "X19BF"="#cf9e65",
  "X19HF"="#ffd700",
  "Alloxanthin"="#ff4500",
  "Diadinoxanthin"="#ffff00",
  "Fucoxanthin"="#ff8c00",
  "Peridinine"="gold4",
  "Prasinoxanthin"="#fffacd",
  "Zeaxanthin"="#ff0000",
  "bb_Carotene"="#800000"
)

# 1) Study ranges of pigments values across SOLA and BBMO
# -------------------------------------------------------

data.pft.pig <- read.csv("Desktop/PETRIMED/DATA/METAG/psbo/mola_bbmo_sola_pigments_psbo_pst_rpkm.csv")
info=colnames(data.pft.pig)[1:8]
mpig<-reshape2::melt(data.pft.pig[,c(info,names(pal.pig) )], id.vars = info)
mpig$pigment<-factor(mpig$variable, levels = names(pal.pig), ordered = TRUE)
mpig$ST<-factor(mpig$ST, levels = c("SOLA", "MOLA", "BBMO"), ordered = TRUE)
mpig$Month<-factor(mpig$Month, levels = month.abb, ordered = TRUE)

# Summary table: KW test per pigment + Dunn's post-hoc
resc<-data.frame()
for (i in names(pal.pig)){
  ds<-mpig[mpig$SF %in% "PICO" & mpig$pigment %in% i,]
  resc[i,"kw"]<-kruskal.test(value ~ ST, ds)$statistic
  resc[i,"pv"]<-kruskal.test(value ~ ST, ds)$p.value
  dun<-dunn.test::dunn.test(ds$value, ds$ST)
  resc[i,"dunn"]<-paste0(dun$comparisons[dun$P<0.05], collapse = "; ")
};resc
pos<-rownames(resc[resc$pv < 0.05,])
pos.tab<-aggregate(value ~ pigment, data = mpig[mpig$SF %in% "PICO" & mpig$pigment %in% pos,], FUN = max)

gc<-ggplot(mpig[mpig$SF %in% "PICO",],aes(x = pigment, y = value, fill = pigment)) +
  geom_jitter(aes(shape = ST),position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2),size = 2.5, alpha = 0.6) +
  geom_boxplot(aes(color = ST),position = position_dodge(1),width = 0.6,alpha = 0.6,outlier.shape = NA,notch=FALSE) +
  geom_text(data = pos.tab,aes(x = pigment, y = value) , label = "*", fontface = "bold", size = 10 )+
  #scale_y_log10(limits = c(0.0001,5), expand = c(0,0)) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10, sigma = 0.0001),breaks = c(0, 10^(-3:2)),labels = label_number(accuracy = 0.001))+
  scale_color_manual(values = c("SOLA"="black","BBMO"="black", "MOLA" ="black"),guide = "none") +
  scale_shape_manual(values = c(21,22, 24)) +
  scale_fill_manual(values = pal.pig, guide = "none") +
  labs(x = "", y = "Concentration (µg/L)") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1));gc

mpig.per<-reshape2::melt(cbind(data.pft.pig[,info], decostand(data.pft.pig[,names(pal.pig)], "total") ), id.vars = info)
mpig.per$pigment<-factor(mpig.per$variable, levels = names(pal.pig), ordered = TRUE)
mpig.per$ST<-factor(mpig.per$ST, levels = c("SOLA", "MOLA", "BBMO"), ordered = TRUE)
mpig.per$Month<-factor(mpig.per$Month, levels = month.abb, ordered = TRUE)

# Summary table: KW test per pigment + Dunn's post-hoc
resp<-data.frame()
for (i in names(pal.pig)){
  ds<-mpig.per[mpig.per$SF %in% "PICO" & mpig.per$pigment %in% i,]
  resp[i,"kw"]<-kruskal.test(value ~ ST, ds)$statistic
  resp[i,"pv"]<-kruskal.test(value ~ ST, ds)$p.value
  dun<-dunn.test::dunn.test(ds$value, ds$ST)
  resp[i,"dunn"]<-paste0(dun$comparisons[dun$P<0.05], collapse = "; ")
};resp
pos<-rownames(resp[resp$pv < 0.05,])
pos.tab<-aggregate(value ~ pigment, data = mpig.per[mpig.per$SF %in% "PICO" & mpig.per$pigment %in% pos,], FUN = max)

gp<-ggplot(mpig.per[mpig.per$SF %in% "PICO",],aes(x = pigment, y = value, fill = pigment)) +
  geom_jitter(aes(shape = ST),position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2),size = 2.5, alpha = 0.6) +
  geom_boxplot(aes(color = ST),position = position_dodge(1),width = 0.6,alpha = 0.6,outlier.shape = NA,notch=FALSE) +
  geom_text(data = pos.tab,aes(x = pigment, y = value) , label = "*", fontface = "bold", size = 10 )+
  #scale_y_log10(limits = c(0.0001,5), expand = c(0,0)) +
  #scale_y_continuous(trans = pseudo_log_trans(base = 10, sigma = 0.0001),breaks = c(0, 10^(-3:2)),labels = label_number(accuracy = 0.001))+
  scale_color_manual(values = c("SOLA"="black","BBMO"="black", "MOLA" ="black"),guide = "none") +
  scale_shape_manual(values = c(21,22, 24)) +
  scale_fill_manual(values = pal.pig, guide = "none") +
  labs(x = "", y = "Proportion (0-1)") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1));gp

library(vegan)
# Concentration PERMANOVA
(perm_conc <- adonis2(data.pft.pig[data.pft.pig$SF=="PICO", names(pal.pig)] ~ data.pft.pig[data.pft.pig$SF=="PICO", "ST"], by="margin"))
# Proportion PERMANOVA  
(perm_prop <- adonis2(decostand(data.pft.pig[data.pft.pig$SF=="PICO", names(pal.pig)], method = "total") ~ data.pft.pig[data.pft.pig$SF=="PICO", "ST"], by="margin"))
# Dispersion (within-site homogeneity)
(disp_conc <- betadisper(vegdist(data.pft.pig[data.pft.pig$SF=="PICO", names(pal.pig)]), data.pft.pig[data.pft.pig$SF=="PICO", "ST"]))
(disp_prop <- betadisper(vegdist(decostand(data.pft.pig[data.pft.pig$SF=="PICO", names(pal.pig)], method = "total")), data.pft.pig[data.pft.pig$SF=="PICO", "ST"]))

# 2) Study ranges of PST psbO values across SOLA and BBMO
# -------------------------------------------------------

pal.pft <- readRDS("~/Desktop/PETRIMED/DATA/PIGMENTS/pst_palette.RData")
pal.pft <- pal.pft[names(pal.pft) %in% colnames(data.pft.pig)]
pfts<-names(pal.pft)
proks<-pfts[17:19]

mpft<-reshape2::melt(data.pft.pig[,c(info,names(pal.pft) )], id.vars = info)
mpft$PFT<-factor(mpft$variable, levels = names(pal.pft), ordered = TRUE)
mpft$ST<-factor(mpft$ST, levels = c("SOLA", "MOLA","BBMO"), ordered = TRUE)
mpft$Month<-factor(mpft$Month, levels = month.abb, ordered = TRUE)

# Summary table: KW test per PST + Dunn's post-hoc
resps<-data.frame()
for (i in pfts){
  ds<-mpft[mpft$SF %in% "PICO" & mpft$PFT %in% i,]
  resps[i,"kw"]<-kruskal.test(value ~ ST, ds)$statistic
  resps[i,"pv"]<-kruskal.test(value ~ ST, ds)$p.value
  dun<-dunn.test::dunn.test(ds$value, ds$ST)
  resps[i,"dunn"]<-paste0(dun$comparisons[dun$P<0.05], collapse = "; ")
};resps
pos<-rownames(resps[resps$pv < 0.05,])
pos.tab<-aggregate(value ~ PFT, data = mpft[mpft$SF %in% "PICO" & mpft$PFT %in% pos,], FUN = max)

aggregate(value ~ PFT+ST, data = mpft[mpft$SF %in% "PICO" & mpft$PFT %in% pos,], FUN = mean)

gpk<-ggplot(mpft[mpft$SF %in% "PICO",],aes(x = PFT, y = value, fill = PFT)) +
  geom_jitter(aes(shape = ST),position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),size = 2.5, alpha = 0.6) +
  geom_boxplot(aes(color = ST),position = position_dodge(0.8),width = 0.6,alpha = 0.6,outlier.shape = NA,notch=FALSE) +
  geom_text(data = pos.tab,aes(x = PFT, y = value) , label = "*", fontface = "bold", size = 10 )+
  scale_y_continuous(trans = pseudo_log_trans(base = 10, sigma = 1),breaks = c(0, 10^(1:6)),labels = label_number(accuracy = 1))+
  scale_color_manual(values = c("SOLA"="black","BBMO"="black", "MOLA" ="black"),guide = "none") +
  scale_shape_manual(values = c(21,22, 24)) +
  scale_fill_manual(values = pal.pft, guide = "none") +
  labs(x = "", y = "psbO RPKM") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1));gpk

ggplot(mpft[mpft$SF %in% "PICO",],aes(x = Month, y = value, fill = PFT)) +
  facet_grid(PFT~., scale = "free")+
  geom_jitter(aes(shape = ST),position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),size = 2.5, alpha = 0.6) +
  geom_boxplot(aes(color = ST),position = position_dodge(0.8),width = 0.6,alpha = 0.6,outlier.shape = NA,notch=FALSE) +
  scale_color_manual(values = c("SOLA"="black","BBMO"="black", "MOLA" ="black"),guide = "none") +
  scale_shape_manual(values = c(21,22, 24)) +
  scale_fill_manual(values = pal.pft, guide = "none") +
  labs(x = "", y = "psbO RPKM") +
  theme_minimal(base_size = 14) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank())

# 3) Envrionmental Dfferences
# -------------------------------------------------

env<-read.csv("/Users/pierreramond/Desktop/PETRIMED/DATA/METAG/metadata/metag_env.csv", header = TRUE)
env<-env[env$SF %in% c("0.2-3 um"),]
env<-env[env$Depth < 10,]

var<-c("Temperature","Salinity","NH4","NO3","NO2","PO4","SIOH4","COP","CHLA")
pal.env<-c("Temperature"="#224282","Salinity" = "gray25","Oxygen"="#913e33", "MES" = "gray50",  "NH4"="gold4", "NO3" = "coral4", "NO2" = "coral2", "PO4" = "steelblue","SIOH4","#f0624f","COP"="#71bf9b","NOP" = "purple" ,"CHLA"= "#46705d")  

menv<-na.omit(reshape2::melt(env[,c("ST", var)] ))
menv$ST<-factor(menv$ST, levels = c("SOLA", "MOLA", "BBMO"))

# Summary table: KW test per pigment + Dunn's post-hoc
rese<-data.frame()
for (i in unique(menv$variable)){
  ds<-menv[menv$variable %in% i,]
  rese[i,"kw"]<-kruskal.test(value ~ ST, ds)$statistic
  rese[i,"pv"]<-kruskal.test(value ~ ST, ds)$p.value
  dun<-dunn.test::dunn.test(ds$value, ds$ST)
  rese[i,"dunn"]<-paste0(dun$comparisons[dun$P<0.05], collapse = "; ")
};rese
pos<-rownames(rese[rese$pv < 0.05,])
pos.tab<-aggregate(value ~ variable, data = menv[menv$variable %in% pos,], FUN = max)

ge<-ggplot(menv,aes(x = variable, y = value, fill = variable)) +
  geom_jitter(aes(shape = ST),position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2),size = 2.5, alpha = 0.6) +
  geom_boxplot(aes(color = ST),position = position_dodge(1),width = 0.6,alpha = 0.6,outlier.shape = NA,notch=FALSE) +
  #scale_y_log10(limits = c(0.0001,10), expand = c(0,0)) +
  geom_text(data = pos.tab,aes(x = variable, y = value) , label = "*", fontface = "bold", size = 10 )+
  scale_y_continuous(trans = pseudo_log_trans(base = 10, sigma = 0.0001),breaks = c(0, 10^(-3:2)),labels = label_number(accuracy = 0.001))+
  scale_color_manual(values = c("SOLA"="black","BBMO"="black", "MOLA" ="black"),guide = "none") +
  scale_shape_manual(values = c(21,22, 24)) +
  scale_fill_manual(values = pal.env, guide = "none") +
  labs(x = "", y = "Values") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1));ge

gpk/gc/gp/ge + plot_annotation(tag_levels = "A")


# 3) Study Correlations between the patterns of the same group across methods
# ----------------------------------------------------------------------------

psbo <- read.csv("~/Desktop/PETRIMED/DATA/METAG/psbo/bbmo_sola_pigments_psbo_pst_rpkm.csv")
metab <-read.csv("~/Desktop/PETRIMED/DATA/PIGMENTS/sola_pigments_pst_ptot.csv")
mags <- read.csv("~/Desktop/PETRIMED/DATA/PIGMENTS/sola_pigments_mags_pst_rpkm.csv")

# a. Harmonise + filter SOLA
# -------------------------------------------------
clean_sf <- function(x){
  case_when(
    x %in% c("0.2-3 um","PICO") ~ "PICO",
    x %in% c("3-20 um","NANO") ~ "NANO",
    TRUE ~ x)
}

prep <- function(df){
  df |>
    filter(ST=="SOLA") |>
    mutate(Date = as.Date(Date),
           SF   = clean_sf(SF)) |>
    arrange(SF,Date)
}

psbo  <- prep(psbo)
metab <- prep(metab)
mags  <- prep(mags)


# b. Inter-comparison of the datasets and relative abundances
# ------------------------------------------------------------

# Helper function
extract_pft_block <- function(df, sf_filter, pfts_to_use, method_name, sf_name) {
  vals <- (colSums(df[df$SF %in% sf_filter, pfts_to_use]) /
             sum(df[df$SF %in% sf_filter, pfts_to_use])) * 100
  tibble(method = method_name,SF     = sf_name,PFT    = names(vals),pct    = as.numeric(vals))}

# create dataset
pfts_no_prok <- pfts[-match(proks, pfts)]
result <- bind_rows(
  extract_pft_block(metab, "PICO", pfts_no_prok, "metab", "PICO"),
  extract_pft_block(metab, "NANO", pfts_no_prok, "metab", "NANO"),
  extract_pft_block(metab, "PICO", proks, "metab", "PICO"),
  extract_pft_block(metab, "NANO", proks, "metab", "NANO"),
  extract_pft_block(psbo,  "PICO", pfts,          "psbo",  "PICO"),
  extract_pft_block(psbo,  "NANO", pfts,          "psbo",  "NANO"),
  extract_pft_block(mags,  "PICO", pfts,          "mags",  "PICO"),
  extract_pft_block(mags,  "NANO", pfts,          "mags",  "NANO")
)

# Define log10 breaks and labels
log_breaks <- c(0, 0.01, 0.1, 1, 10, 100)
log_labels <- c("0", "0.01", "0.1", "1", "10", "100")
bin_colors <- c("0"= "grey92","< 0.1%"= "#d1e5f0","0.1–1%"= "#4393c3","1–10%"= "#f4a582","> 10%"= "#b2182b")

result$PFT<-factor(result$PFT, levels = rev(pfts), ordered = TRUE)
result$SF<-factor(result$SF, levels = c("PICO", "NANO"), ordered = TRUE)
result$method<-factor(result$method, levels = c("metab", "psbo","mags") ,labels =  c("Metabarcoding", "psbO", "MAGs"), ordered = TRUE)
result$Domain<-factor(ifelse(result$PFT %in% proks, "Prokaryotes",  "Eukaryotes"),levels = c("Eukaryotes", "Prokaryotes"), ordered = TRUE)
result$pct_bin<-cut(result$pct, breaks = c(0,0.01, 0.1, 1, 10,100), labels = names(bin_colors),ordered_result = TRUE, right = FALSE)
  
grab<-ggplot(result, aes(x = method, y = PFT, fill = pct_bin)) +
  facet_grid(Domain~SF, scales = "free", space = "free")+
  geom_tile(color = "transparent", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct > 0.05, sprintf("%.1f", pct), "")),size = 2.8, color = "grey20") +
  scale_fill_manual(values = bin_colors,name   = "Relative\nabundance",drop   = FALSE,guide  = guide_legend(reverse = FALSE)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.y =  element_blank(),
        panel.grid      = element_blank(),
        legend.position = "right") +
  labs(x = "", y = "PST");grab

# c. Z-score per SF (safe)
# -------------------------------------------------
safe_scale <- function(x){
  if(sd(x,na.rm=TRUE)==0) return(rep(0,length(x)))
  as.numeric(scale(x))
}

zscore <- function(df){
  df |>
    group_by(SF) |>
    mutate(across(all_of(pfts), safe_scale)) |>
    ungroup()
}

psbo_z  <- zscore(psbo)
metab_z <- zscore(metab)
mags_z  <- zscore(mags)

# d. Merge ALL THREE (shared Date + SF only)
# -------------------------------------------------

merged_all <- metab_z |>
  select(Date,SF,all_of(pfts)) |>
  rename_with(~paste0(.x,"_meta"),all_of(pfts)) |>
  inner_join(
    psbo_z |>
      select(Date,SF,all_of(pfts)) |>
      rename_with(~paste0(.x,"_psbo"),all_of(pfts)),
    by=c("Date","SF")
  ) |>
  inner_join(
    mags_z |>
      select(Date,SF,all_of(pfts)) |>
      rename_with(~paste0(.x,"_mags"),all_of(pfts)),
    by=c("Date","SF")
  ) |>
  arrange(SF,Date)

# e. Per-PFT Spearman trend concordance
# -------------------------------------------------

trend_cor <- function(df, tag1, tag2){
  df %>%
    group_by(SF) %>%                  # compute separately for each size-fraction
    group_modify(~{
      map_dfr(pfts, function(p){
        x <- .x[[paste0(p,"_",tag1)]]
        y <- .x[[paste0(p,"_",tag2)]]
        if(sd(x,na.rm=TRUE)==0 | sd(y,na.rm=TRUE)==0){
          tibble(PFT=p, rho=NA, p.value=NA)
        } else {
          t <- cor.test(x, y, method="pearson")
          tibble(PFT=p, rho=t$estimate, p.value=t$p.value)
        }
      })
    }) %>% ungroup()
}

cor_meta_psbo <- trend_cor(merged_all,"meta","psbo") %>% mutate(comp="Metabarcoding-psbO")
cor_meta_mags <- trend_cor(merged_all,"meta","mags") %>% mutate(comp="Metabarcoding-MAGs")
cor_psbo_mags <- trend_cor(merged_all,"mags","psbo") %>% mutate(comp="MAGs-psbO")
cor_results <- bind_rows(cor_meta_psbo, cor_meta_mags, cor_psbo_mags)
cor_results$PFT<-factor(cor_results$PFT, levels = rev(pfts), ordered = TRUE)
cor_results$SF<-factor(cor_results$SF, levels = c("PICO", "NANO"), ordered = TRUE)
cor_results$Domain<-factor(ifelse(cor_results$PFT %in% proks, "Prokaryotes",  "Eukaryotes"),levels = c("Eukaryotes", "Prokaryotes"), ordered = TRUE)

gcor<-ggplot(cor_results,aes(comp,PFT,fill=rho))+
  facet_grid(Domain~SF, scale = "free", space = "free")+
  geom_tile()+
  geom_text(aes(label = ifelse(p.value > 0.05, "NS", "")),size = 2.8, color = "grey20") +
  scale_fill_gradient2(limits=c(-1,1), low = "coral4", mid = "white", high = "steelblue4")+
  theme_minimal(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        panel.grid      = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0))+
  labs(fill="Correlation\nPearson p",x="",y="");gcor

grab + gcor + plot_layout(axis_titles = "collect", guides = "collect") +  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 14, face = "bold"))

aggregate(rho~comp+SF, cor_results, FUN = mean)

# f. Δ (first derivative) concordance
# -------------------------------------------------

delta_cor <- function(df, tag1, tag2){
  df %>%
    group_by(SF) %>%
    group_modify(~{
      map_dfr(pfts, function(p){
        x <- diff(.x[[paste0(p,"_",tag1)]])
        y <- diff(.x[[paste0(p,"_",tag2)]])
        if(sd(x,na.rm=TRUE)==0 | sd(y,na.rm=TRUE)==0){
          tibble(PFT=p, rho_delta=NA, p.value_delta=NA)
        } else {
          t <- cor.test(x, y, method="spearman")
          tibble(PFT=p, rho_delta=t$estimate, p.value_delta=t$p.value)
        }
      })
    }) %>% ungroup()
}

delta_meta_psbo <- delta_cor(merged_all,"meta","psbo") %>% mutate(comp="meta-psbo")
delta_meta_mags <- delta_cor(merged_all,"meta","mags") %>% mutate(comp="meta-mags")
delta_mags_psbo <- delta_cor(merged_all,"mags","psbo") %>% mutate(comp="mags-psbo")
delta_results <- bind_rows(delta_meta_psbo, delta_meta_mags, delta_mags_psbo)
delta_results$PFT<-factor(delta_results$PFT, levels = rev(pfts), ordered = TRUE)
delta_results$SF<-factor(delta_results$SF, levels = c("PICO", "NANO"), ordered = TRUE)

gds<-ggplot(delta_results,aes(comp,PFT,fill=rho_delta))+
  facet_grid(.~SF)+
  geom_tile()+
  scale_fill_gradient2(limits=c(-1,1), low = "coral4", mid = "white", high = "seagreen4")+
  theme_minimal(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  labs(fill="Directional\nSynchrony\nSpearman ρ",x="Method comparison",y="");gds

gcor + gds

# g. Procrustes per SF
# -------------------------------------------------

get_pair_matrices <- function(df, sf, tag1, tag2){
  sub <- df |> filter(SF==sf)
  mat1 <- sub |> select(ends_with(tag1)) |> setNames(pfts) |> as.matrix()
  mat2 <- sub |> select(ends_with(tag2)) |> setNames(pfts) |> as.matrix()
  keep <- apply(mat1,2,sd)!=0 & apply(mat2,2,sd)!=0
  list(mat1=mat1[,keep,drop=FALSE],mat2=mat2[,keep,drop=FALSE]) }

run_procrustes <- function(df, sf, tag1, tag2){
  mats <- get_pair_matrices(df,sf,tag1,tag2)
  protest(rda(mats$mat1),rda(mats$mat2),permutations=999)}

proc_meta_psbo_pico <- run_procrustes(merged_all,"PICO","meta","psbo")
proc_meta_mags_pico <- run_procrustes(merged_all,"PICO","meta","mags")
proc_mags_psbo_pico <- run_procrustes(merged_all,"PICO","mags","psbo")

proc_meta_psbo_nano <- run_procrustes(merged_all,"NANO","meta","psbo")
proc_meta_mags_nano <- run_procrustes(merged_all,"NANO","meta","mags")
proc_mags_psbo_nano <- run_procrustes(merged_all,"NANO","mags","psbo")

proc_results<-data.frame(ss = c(proc_meta_psbo_pico$ss,proc_meta_mags_pico$ss, proc_mags_psbo_pico$ss, proc_meta_psbo_nano$ss, proc_meta_mags_nano$ss, proc_mags_psbo_nano$ss),
                  #scale = c(proc_meta_psbo_pico$scale,proc_meta_mags_pico$scale, proc_mags_psbo_pico$scale, proc_meta_psbo_nano$scale, proc_meta_mags_nano$scale, proc_mags_psbo_nano$scale),
                  SF = c(rep("PICO",3),rep("NANO",3)),comp = rep(c("meta-psbo", "meta-mags", "mags-psbo"), 2))

# h. plot zscores
# -------------------------------------------------

# Combine datasets for plotting
plot_df <- bind_rows(
  metab_z  |> select(Date,SF,all_of(pfts)) |> pivot_longer(all_of(pfts), names_to="PFT", values_to="z") |> mutate(Method="metab"),
  psbo_z   |> select(Date,SF,all_of(pfts)) |> pivot_longer(all_of(pfts), names_to="PFT", values_to="z") |> mutate(Method="psbo"),
  mags_z   |> select(Date,SF,all_of(pfts)) |> pivot_longer(all_of(pfts), names_to="PFT", values_to="z") |> mutate(Method="mags")
) |> arrange(SF, PFT, Method, Date)

# For better separation in Y-axis, create a combined PFT-Method label
plot_df <- plot_df |> mutate(PFT_method = paste(PFT, Method, sep="_"))
plot_df$PFT<-factor(plot_df$PFT, levels = pfts, ordered = TRUE)
plot_df$SF<-factor(plot_df$SF, levels = c("PICO", "NANO"), ordered = TRUE)
plot_df$Method<-factor(plot_df$Method, levels = rev(c("metab", "psbo", "mags")), labels = rev(c("Metabarcoding", "psbO", "MAGs")), ordered = TRUE)

ggplot(plot_df,aes(x=Date, y= Method, fill=z)) +
  facet_grid(PFT~SF)+
  geom_tile() +
  scale_fill_gradient2(low="steelblue4", mid="white", high="coral3", midpoint=0, limits = c(-2, 2),oob = scales::oob_squish) +
  labs(y="PFT + Method", x="Date", fill="Z-score") +
  scale_x_datetime(breaks = "3 months", expand = c(0.025, 0), date_labels = "%Y-%b")+
  theme_minimal(base_size=14) +
  theme(axis.text.y = element_text(size=10),
        strip.text.y = element_text(angle=0, hjust=0),panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))


# Modeling results
model<-readxl::read_xlsx("Desktop/PETRIMED/ANALYSES/ML/results_per_method_and_PFTs.xlsx", sheet = 2)
model<-model[!model$PFT %in% c("MOCH"),]
model$PFT<-factor(model$PFT, levels = rev(pfts), ordered = TRUE)
model$SF<-factor(model$SF, levels = c("PICO", "NANO"), ordered = TRUE)
model$Method<-factor(model$Method, levels = c("Metabarcoding", "psbO", "MAGs"), ordered = TRUE)
model[,2:6]<-apply(model[,2:6], 2, as.numeric)
model<-rbind(model,cbind(PFT = "Average",aggregate(.~SF+Method,model[,-1], FUN = mean)[,c(3:7,1:2)]))
model$Domain<-factor(ifelse(model$PFT %in% proks, "Prokaryotes", ifelse(model$PFT %in% "Average", "", "Eukaryotes")),levels = c("Eukaryotes", "Prokaryotes", ""), ordered = TRUE)

ggplot(model,aes(x= Method, y= PFT, fill=cor)) +
  facet_grid(Domain~SF, scales = "free", space = "free")+
  geom_tile() +
  scale_fill_gradient2(low="coral3", mid="white", high="steelblue4", midpoint=0, limits = c(0, 1),na.value = "gray95") +
  theme_minimal(base_size=14) +
  labs(y="PST", x="Method", fill= bquote('Prediction'~(R^2))) +
  theme(strip.text.y = element_text(angle=0, hjust=0),panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))

# colors for PFTs
pal.pft<-readRDS("~/Desktop/PETRIMED/DATA/PIGMENTS/pst_palette.RData")
pal.pft=pal.pft[names(pal.pft) %in% model$PFT]

ggplot(model[!model$Domain %in% "",],aes(x= avg.ab, y= cor, color=PFT, group = SF)) +
  facet_grid(SF~Method, scale = "free")+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)+
  scale_color_manual(values = pal.pft)+
  scale_x_log10()+
  theme_minimal(base_size=14) +
  labs(y="Prediction", x="Average Abundance", color= "PST") +
  theme(strip.text.y = element_text(angle=0, hjust=0),panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))







