## This program aims to analyze parameters obatined using the RK4 method of solving ODEs written for RNAi in the python program '2022_6_29_Celegans_RNAi_ODEs_RK4_method_d6_param_sweep.py'.

experiment <- "2022_7_14_RNAiDynamics"

root <- "/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/RNA_silencing/"
setwd(paste0(root, "code/R"))
dir.create(paste0(root, "analyses/R/Figures/", experiment))
dir.create(paste0(root, "analyses/R/Tables/", experiment))
path_to_figure_outputs <- paste0(root, "analyses/R/Figures/", experiment)
path_to_table_outputs <- paste0(root, "analyses/R/Tables/", experiment)
path_to_parameters <- paste0(root, "analyses/Python/Tables/", experiment)

# Import libraries for analysis and plotting.
library(tidyverse)
library(grid)
library(reshape2)
library(ggtext)
library(readxl)

## Get the file with all parameters.
C_elegans_RNADynamics <- data.frame(read.csv(paste0(path_to_parameters, "/RNAiDynamics_param_sweep.csv"), sep = ",", header = TRUE))
Silenced <- C_elegans_RNADynamics %>% filter(kd_time != 0)

### Set plot themes
my_theme <- theme(panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_line(linetype = 1, size = 0.5, lineend="round"),
                  axis.title = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10),
                  axis.ticks = element_line(size = 0.5), axis.ticks.length = unit(.1, "cm"),
                  axis.text = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10))

## PLOT MANY
## kd_time vs. time_to_kd
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary lm
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = lm), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "lm") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_lm.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary lds
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = lds), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "lds") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_lds.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary Tm
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = Tm), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "Tm") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_Tm.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary Tp
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = Tp), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "Tp") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_Tp.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary Tug
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = Tug), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "Tug") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_Tug.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary Tpri
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = Tpri), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "Tpri") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_Tpri.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary Tsec
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = Tsec), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "Tsec") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_Tsec.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary k1
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = k1), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "k1") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_k1.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary k2
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = k2), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "k2") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_k2.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary k3
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = k3), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "k3") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_k3.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary k4
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = k4), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "k4") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_k4.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary k5
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = k5), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "k5") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_k5.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary k6
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = k6), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "k6") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_k6.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## kd_time vs. time_to_kd vary k8
max_kd_time = round(max(Silenced$`kd_time`))
max_time_to_kd = round(max(Silenced$`time_to_kd`))
kd_subset <- ggplot(Silenced, aes(x = `kd_time`, y = `time_to_kd`, color = k8), alpha = 0.5) +
  my_theme + ylab("time to kd") + xlab("kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "k8") +
  scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 500), breaks= c(0, 250, 500), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_with_kd_highlight_k8.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## max_sec vs. min_m for low values of secondary - inset
max_of_max_sec = round(max(C_elegans_RNADynamics$`max_sec`))
max_of_min_m = round(max(C_elegans_RNADynamics$`min_m`))
kd_subset <- ggplot(C_elegans_RNADynamics, aes(x = `min_m`, y = `max_sec`), alpha = 0.5) +
  my_theme + ylab("max 2º siRNA") + xlab("min mRNA") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_of_min_m + 1), breaks = c(0, 0.5*(max_of_min_m + 1), (max_of_min_m + 1)), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_of_max_sec + 1), breaks= c(0, 0.5*(max_of_max_sec+ 1), (max_of_max_sec + 1)), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_sec_vs_mRNA.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## max_sec vs. min_m for low values of secondary - inset
kd_subset <- ggplot(C_elegans_RNADynamics, aes(x = `min_m`, y = `max_sec`), alpha = 0.5) +
  my_theme + ylab("max 2º siRNA") + xlab("min mRNA") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, 20), breaks = c(0, 10, 20), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 200), breaks= c(0, 100, 200), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_sec_vs_mRNA_inset.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## max_sec vs. min_m for low values of secondary - inset expanded.
kd_subset <- ggplot(C_elegans_RNADynamics, aes(x = `min_m`, y = `max_sec`), alpha = 0.5) +
  my_theme + ylab("max 2º siRNA") + xlab("min mRNA") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, 50), breaks = c(0, 25, 50), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 50), breaks= c(0, 25, 50), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_sec_vs_mRNA_inset_low_siRNA.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## max_sec vs. min_m for low values of secondary - inset expanded colored by k5
kd_subset <- ggplot(C_elegans_RNADynamics, aes(x = `min_m`, y = `max_sec`, color = k5), alpha = 0.5) +
  my_theme + ylab("max 2º siRNA") + xlab("min mRNA") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "k5") +
  scale_x_continuous(limits = c(0, 50), breaks = c(0, 25, 50), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 50), breaks= c(0, 25, 50), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_sec_vs_mRNA_inset_low_siRNA.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)
