## This program aims to analyze parameters obatined using the RK4 method of solving ODEs written for RNAi in the python program '2022_6_29_Celegans_RNAi_ODEs_RK4_method_d6_param_sweep_non_steady_state.py'.

experiment <- "2022_7_27_RNAiDynamics"

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
C_elegans_RNADynamics <- data.frame(read.csv(paste0(path_to_parameters, "/RNAiDynamics_non_SS_param_sweep.csv"), sep = ",", header = TRUE))
Trancription_change <- C_elegans_RNADynamics %>% group_by(lm,lds,Tm,Tp,Tug,Tpri,Tsec,k1,k2,k3,k4,k5,k6,k8) %>% mutate(kd_time_change = (max(kd_time)-min(kd_time))/max(kd_time), time_to_kd_change = (max(time_to_kd)-min(time_to_kd))/max(time_to_kd)) %>%
                                filter(kd_time_change != 0 && time_to_kd_change != 0)

### Set plot themes
my_theme <- theme(panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_line(linetype = 1, size = 0.5, lineend="round"),
                  axis.title = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10),
                  axis.ticks = element_line(size = 0.5), axis.ticks.length = unit(.1, "cm"),
                  axis.text = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10))

## PLOT
## change in kd_time vs. change in time_to_kd with doubling of transcription rate.
max_kd_time = round(max(Trancription_change$`kd_time_change`))
max_time_to_kd = round(max(Trancription_change$`time_to_kd_change`))
kd_subset <- ggplot(Trancription_change, aes(x = `kd_time_change`, y = `time_to_kd_change`, color = `kd_time`), alpha = 0.5) +
  my_theme + ylab("% change in time to kd") + xlab("% change in kd time") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "kd time") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 1), breaks= c(0, 0.5, 1), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAiDynamics_on_doubling_transcription.eps"), kd_subset, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)
