library(readxl)
library(ggplot2)
library(dplyr)
library("astronomyengine")
library(patchwork)

# Season function
get_season <- function(Date, lat) {
  Date <- as.Date(Date)
  years <- as.integer(format(Date, "%Y"))
  out <- character(length(Date))
  for (y in unique(years)) {
    i <- which(years == y)
    s <- astro_seasons(y)
    dates <- as.Date(as.POSIXct(unlist(s), tz = "UTC"))
    if (mean(lat[i], na.rm = TRUE) >= 0) {
      out[i] <- ifelse(Date[i] >= dates[1] & Date[i] < dates[2], "Spring",
                       ifelse(Date[i] >= dates[2] & Date[i] < dates[3], "Summer",
                              ifelse(Date[i] >= dates[3] & Date[i] < dates[4], "Autumn", "Winter")))
    } else {
      out[i] <- ifelse(Date[i] >= dates[1] & Date[i] < dates[2], "Autumn",
                       ifelse(Date[i] >= dates[2] & Date[i] < dates[3], "Winter",
                              ifelse(Date[i] >= dates[3] & Date[i] < dates[4], "Spring", "Summer")))
    }
  }
  factor(out, levels = c("Winter", "Spring", "Summer", "Autumn"))
}

# Abundance from coverM
ab<-readRDS("~/Desktop/PETRIMED/DATA/METAG/mags/MAGs_BBMO_SOLA_MOLA_coverm_euk01_prok20.RData")

# Distinguish Prokaryotes and Eukaryotes
ab$Domain<-ifelse(grepl("EUK", ab$Genome), "Eukaryotes", "Prokaryotes" )

# Sites
ab$ST<-ifelse(grepl("SO",ab$samples),"SOLA", 
              ifelse(grepl("MO",ab$samples),"MOLA",
                     ifelse(grepl("BL",ab$samples),"BBMO","VIDA")))

# coordinates
coord<-data.frame(lat = c(41.665, 45.5, 42.488702, 42.5055556), long = c(2.805, 13.6, 3.142906, 3.6999999999999997), ST = c("BBMO", "VIDA", "SOLA", "MOLA"))
ab<-merge(ab, coord, by = "ST")

# Date
dates<-gsub("_.*|SO|PICO_MO|NANO_SO|PICO_SO|BL|BF|_2_1_QC|_1_QC|.1_QC","",ab$samples)
dates<-as.Date(ifelse(nchar(dates) < 6 ,paste0("0", dates),dates), format = "%y%m%d")
ab$Date<-c(dates)
ab$Season<-get_season(Date = ab$Date, lat = ab$lat)

# Depth
ab$Depth<-as.numeric(ifelse(stringr::str_split_fixed(ab$samples, pattern = "_", n = 4)[,3] %in% c("", "QC", "1"), 1, stringr::str_split_fixed(ab$samples, pattern = "_", n = 4)[,3]))

# Size-fraction
ab$SF<-ifelse(grepl("NANO",ab$samples), "NANO", "PICO" )
length(grep("EUK",unique(ab$Genome)))
length(grep("PROK",unique(ab$Genome)))

ab.dom<-aggregate(Relative.Abundance ~ Season+Domain+ST+samples+Date+Depth+lat+long+SF,data = ab, FUN = sum)
ab.dom<-ab.dom[ab.dom$Depth < 10,]
ab.dom$STSF<-factor(paste(ab.dom$ST,ab.dom$SF, sep = " "), levels = c("SOLA PICO", "SOLA NANO","BBMO PICO","MOLA PICO" ), ordered = TRUE)

summary(ab.dom[ab.dom$Domain == "Eukaryotes","Relative.Abundance"])
summary(ab.dom[ab.dom$Domain == "Prokaryotes","Relative.Abundance"])

gcomp<-ggplot(data =ab.dom, aes( x = samples, y= Relative.Abundance, fill = Domain))+
  facet_grid(.~STSF, space = "free", scales = "free")+
  geom_bar(position = "stack", stat = "identity")+
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"), expand = c(0, 0))+
  scale_fill_manual(values = c("steelblue4", "seagreen4") )+
  labs(x = "", y = "") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_blank());gcomp

# Summary table: KW test per pigment + Dunn's post-hoc
ab.tot<-aggregate(Relative.Abundance ~ Season+ST+samples+Date+Depth+lat+long+SF,data = ab, FUN = sum)
ab.tot<-ab.tot[ab.tot$Depth < 10,]
summary(ab.tot$Relative.Abundance)

ds<-ab.tot[ab.tot$SF %in% "PICO",]
kruskal.test(Relative.Abundance ~ ST, ds)$statistic
kruskal.test(Relative.Abundance ~ ST, ds)$p.value
dun<-dunn.test::dunn.test(ds$Relative.Abundance, ds$ST)
dun$comparisons[dun$P<0.05]

summary(ds[ds$ST == "MOLA","Relative.Abundance"])
summary(ds[ds$ST == "SOLA","Relative.Abundance"])
summary(ds[ds$ST == "BBMO","Relative.Abundance"])

gbst<-ggplot(ds, aes(x = ST, y = Relative.Abundance)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2),size = 2.5, alpha = 0.6) +
  geom_boxplot(position = position_dodge(1),width = 0.6,alpha = 0.6,outlier.shape = NA,notch=TRUE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"), expand = c(0, 0))+
  labs(x = "", y = "") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1));gbst

kruskal.test(Relative.Abundance ~ Season, ds)$statistic
kruskal.test(Relative.Abundance ~ Season, ds)$p.value
dun<-dunn.test::dunn.test(ds$Relative.Abundance, ds$Season)
dun$comparisons[dun$P<0.05]

summary(ds[ds$Season == "Autumn","Relative.Abundance"])
summary(ds[ds$Season == "Winter","Relative.Abundance"])
summary(ds[ds$Season == "Spring","Relative.Abundance"])
summary(ds[ds$Season == "Summer","Relative.Abundance"])

gbseas<-ggplot(ds, aes(x = Season, y = Relative.Abundance)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2),size = 2.5, alpha = 0.6) +
  geom_boxplot(position = position_dodge(1),width = 0.6,alpha = 0.6,outlier.shape = NA,notch=TRUE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"), expand = c(0, 0))+
  labs(x = "", y = "") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1));gbseas

# Summary table: KW test per pigment + Dunn's post-hoc
ds<-ab.tot
table(ds$ST, ds$SF)
kruskal.test(Relative.Abundance ~ SF, ds)$statistic
kruskal.test(Relative.Abundance ~ SF, ds)$p.value
dun<-dunn.test::dunn.test(ds$Relative.Abundance, ds$SF)
dun$comparisons[dun$P<0.05]

summary(ab.tot[ab.tot$SF == "NANO","Relative.Abundance"])
summary(ab.tot[ab.tot$SF == "PICO","Relative.Abundance"])

gbsf<-ggplot(ds, aes(x = SF, y = Relative.Abundance)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 1, jitter.width = 0.2),size = 2.5, alpha = 0.6) +
  geom_boxplot(position = position_dodge(1),width = 0.6,alpha = 0.6,outlier.shape = NA,notch=TRUE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"), expand = c(0, 0))+
  labs(x = "", y = "") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1));gbsf

gcomp/(gbsf+ gbst+gbseas) & 
  plot_annotation(tag_levels = "A","Percentage of metagenomics reads recruited by our set of MAGs",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))


#####################
# 
library(ggplot2)
library(cowplot)

euk <- read.csv("Desktop/PETRIMED/DATA/METAG/mags/euk_mags_metadata.csv")
prok <- read.csv("Desktop/PETRIMED/DATA/METAG/mags/prok_mags_metadata.csv")

gi <- as.data.frame(rbind(euk[,1:9], prok[,1:9]))
gi$Domain <- factor(
  ifelse(grepl("EUK", gi$genome), "Eukaryote", "Prokaryote"),
  levels = c("Eukaryote", "Prokaryote"),
  ordered = TRUE
)

domain_cols <- c("Prokaryote" = "steelblue4", "Eukaryote" = "seagreen4")

# Main plot
aggregate(contamination~ Domain,data = gi, FUN = mean)
aggregate(completeness~ Domain,data = gi, FUN = mean)
aggregate(sum_len~ Domain,data = gi, FUN = mean)

pmain <- ggplot() +
  geom_point(data = gi[gi$Domain == "Prokaryote",], aes(y = completeness, x = contamination ), color = "steelblue4",size = 2.2, alpha = 0.65 )+
  geom_point(data = gi[gi$Domain == "Eukaryote",], aes(y = completeness, x = contamination), color = "seagreen4" ,size = 2.2, alpha = 0.65)+
  #scale_color_manual(values = domain_cols) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Contamination (%)", y = "Completeness (%)", color = "") +
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(size = 13),
    axis.text = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey85"),
    legend.position = c(0.82, 0.18),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

# Marginal boxplot: contamination
xdens <- axis_canvas(pmain, axis = "x") +
  geom_boxplot(
    data = gi,
    aes(x = contamination, y = 0, fill = Domain),
    width = 0.8,
    alpha = 0.55,
    linewidth = 0.3,
    outlier.shape = NA
  ) +
  scale_fill_manual(values = domain_cols) +
  theme_void() +
  theme(legend.position = "none")

# Marginal boxplot: completeness
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
  geom_boxplot(
    data = gi,
    aes(x = completeness, y = 0, fill = Domain),
    width = 0.8,
    alpha = 0.55,
    linewidth = 0.3,
    outlier.shape = NA
  ) +
  coord_flip() +
  scale_fill_manual(values = domain_cols) +
  theme_void() +
  theme(legend.position = "none")

# Assemble
p1 <- insert_xaxis_grob(
  pmain, xdens,
  grid::unit(0.20, "null"),
  position = "top"
)

p2 <- insert_yaxis_grob(
  p1, ydens,
  grid::unit(0.20, "null"),
  position = "right"
)

ggdraw(p2)
