## This program explores the parameter subspace of interest that were found after running '2022_6_13_RNAi_in_Celegans_linear_modified.R'.
## Assumptions:
# 1. All reactions reach equilibrium between one instance of calculation time and the next. Specifically, in this model, 1ºsiRNAs are made, then polyUG RNAs are made, then 2ºsiRNAs are made. This means that different amounts of mRNA are available for 1ºsiRNA and 2ºsiRNA to act on. Also, effect of 2ºsiRNA on pre-mRNA are simulated before its effects on mRNA.
# 2. No tertiary small RNAs are made for silencing of somatic targets as supported by Pak et al., Cell, 2012 ("Protection from feed-forward amplification in an amplified RNAi mechanism"). Note that Sapetschnig et al., PLoS Genetics, 2015 provide evidence for the formation of tertiary small RNAs during paramutation-like silencing of a gene expressed in the germline ("Tertiary siRNAs mediate paramutation in C. elegans").
# 3. There is no recycling of full-length mRNA or full-length pre-mRNA after small RNA binding and there are no other mechanisms for the turnover of species in the timescale considered.

## Question: What are the mutual constraints on different RNA species that result from the RNA interference pathway?

## Definitions:
# pri = amount of 1ºsiRNA
# pri_m = amount of 1ºsiRNA-mRNA complex. Formation of this complex leads to processing of mRNA into poly-UG RNAs.
# ug = amount of poly UG RNA
# sec = amount of 2ºsiRNA
# sec_m = amount of 2ºsiRNA complex. Formation of this complex leads to the degradation of mRNA
# sec_pm = amount of 2ºsiRNA-pre-mRNA complex.
# pm = amount of pre-mRNA.
# m = amount of mRNA
# p = amount of protein.
# ds = amount of dsRNA.
# ds_l = length of dsRNA.
# m_l = length of mRNA.
# DNA = number of transcribed loci, typically 2.
# pm_i, m_i = pre-mRNA and mRNA after RNA interference

## Equations for each species:
# pri = a*ds_l/22*ds
# pri_m = A*pri*m
# ug = b*pri_m
# sec = c*ug*m_l
# sec_m = B*sec*(m - A*pri*m) # This simulates loss of mRNA to 1ºsiRNA-mRNA and subsequent pUG RNA formation.
# sec_pm = C*sec*pm
# pm_i = pm + e*DNA(1 - d*sec_pm)
# m_i = m - A*pri*m + f*pm - sec_m ## mRNA after 1º + newly spliced mRNA - mRNA bound by 2º

## constants:
# A - 1ºsiRNA binding to mRNA
# B - 2ºsiRNA binding to mRNA
# C - 2ºsiRNA binding to pre-mRNA
# a - rate of 1ºsiRNA production from dsRNA
# b - rate of polyUG-RNA production from 1ºsiRNA-bound mRNA*
# c - rate of 2ºsiRNA production from polyUG-RNA
# d - rate of inhibition of transcription
# e - rate of transcription
# f - rate of splicing and mRNA export
# g - rate of protein synthesis

## Equations for the measured species:

# ug = ds_l/22*A*a*b*ds*m

# pm_i = pm + e*DNA(1-C*c*d*pm*ug*m_l)
# pm_i = pm + e*DNA - e*DNA*d*C*c*ds_l/22*A*a*b*ds*m*m_l*pm
# Simplifying above:
# pm_i = pm + k1 - k1*k2*m*m_l*pm, where k1 = e*DNA, k2 = A*a*ds_l/22*ds and k3 = b*c*d*C

# m_i = (m - A*pri*m) + f*pm - (B*sec*(m - A*pri*m)) ## (mRNA after 1º) + newly spliced mRNA - mRNA lost to secondary.
# Simplifying above:
# m_i = f*pm + m*(1 - A*pri - B*sec + B*sec*A*pri)
# m_i = f*pm + m*(A*pri - 1)(B*sec - 1)
# m_i = f*pm + m*(1 - A*a*ds_l/22*ds - B*c*ds_l/22*A*a*b*ds*m*m_l - B*c*m*m_l*(A*a*ds_l/22*ds)^2)
# m_i = k5*pm_i + m*(1 - k2 - k2*k4*m*m_l - k2^2*k4*m*m_l), where k4 = B*c and k5 = f

# pm_i and m_i are values after RNA interference.

# These consolidated constants and their approx. meaning are:
# k1 = e*DNA          : a measure of transcription
# k2 = A*a*ds_l/22*ds : a measure of dsRNA processing, 1ºsiRNA binding mRNA, and subsequent effect on mRNA
# k3 = b*c*d*C        : a measure of 2ºsiRNA production and effect on pre-mRNA
# k4 = B*c            : a measure of 2ºsiRNA production and effect on mRNA
# k5 = f              : a measure of mRNA splicing and export
# Simulating these minimal equations provides an idea of how mRNA and pre-mRNA vary in response to RNAi
# Constraining conditions for the constants are:
# 1. pm_i and m_i cannot be < 0, and
# 2. m_i < m to signify knockdown.

experiment <- "2022_6_7_RNAi_in_Celegans_linear"

root <- "/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/RNA_silencing/"
setwd(paste0(root, "code/R"))
dir.create(paste0(root, "analyses/R/Figures/", experiment))
dir.create(paste0(root, "analyses/R/Tables/", experiment))
path_to_figure_outputs <- paste0(root, "analyses/R/Figures/", experiment)
path_to_table_outputs <- paste0(root, "analyses/R/Tables/", experiment)

## run conditions for varying number of experiments to be considered. Comment out different conditions below and reassign variables accordingly.
run_condition <- "vary_pm_m_k1_to_k5_subspeces_premRNA_increase_good_kd"

library(tidyverse)
library(reshape2)

parameters <- data.frame()
parameters_varying_length <- data.frame()
r_low <- 0 # lower bound for random number range.
r_high <- 10 # upper bound for random number range.

for (s in 1:10) {
  set.seed(s) ## Set a number of different random seeds to get reproducible code and yet explore a range of random numbers.
  for (i in 1:10000000) {
    m = runif(1, r_low, r_high)
    pm = runif(1, r_low, r_high)
    m_p_ratio = m / pm ## Starting ratio of mRNA to pre-mRNA
    k1 = runif(1, r_low, r_high)
    k2 = runif(1, r_low, r_high)
    k3 = runif(1, r_low, r_high)
    k4 = runif(1, r_low, r_high)
    k5 = runif(1, r_low, r_high)
    ## Simulating for gene of a given length.
    m_l = 1
    pm_i = pm + k1 - k1*k2*k3*m*m_l*pm
    m_i = k5*pm_i + m*(1 - k2 - k2*k4*m*m_l - k2^2*k4*m*m_l)
    kd_i = m_i / m
    m_i_pm_i_ratio = m_i / pm_i ## Ending ratio of mRNA to pre-mRNA
    if ( m_i < m &&
       pm_i > 0 &&
       m_i > 0 &&
       kd_i < 0.1 &&
       m_i_pm_i_ratio > m_p_ratio){ ## These conditions are minimal constraints on the constants k1 through k4
      parameters <- rbind(parameters, cbind(s, m, pm, k1, k2, k3, k4, k5, m_i, pm_i, m_p_ratio, m_i_pm_i_ratio, kd_i))
    }
    ## Simulating genes of different lengths. Valid k1 through k5 will now take on different values based on the values of m_l.
    m_l = runif(1, 0, 1) # this simulates differences in mRNA lengths
    pm_il = pm + k1 - k1*k2*k3*m*m_l*pm
    m_il = k5*pm_il + m*(1 - k2 - k2*k4*m*m_l - k2^2*k4*m*m_l)
    kd_il = m_il / m
    m_il_pm_il_ratio = m_il / pm_il ## Ending ratio of mRNA to pre-mRNA
    if ( m_il < m &&
        pm_il > 0 &&
        m_il > 0 &&
        kd_il < 0.1 &&
        m_i_pm_i_ratio > m_p_ratio){ ## These conditions are minimal constraints on the constants k1 through k4
      parameters_varying_length <- rbind(parameters_varying_length, cbind(s, m, pm, k1, k2, k3, k4, k5, m_il, pm_il, m_l, m_p_ratio, m_il_pm_il_ratio, kd_il))
    }
## Provide column names and write out tables with results.
colnames(parameters) <- c("seed", "mRNA", "pre-mRNA", "k1","k2","k3","k4", "k5", "mRNA_after_RNAi", "pre-mRNA_after_RNAi", "mRNA/pre-mRNA_initial", "mRNA/pre-mRNA_after_RNAi", "knockdown")
write.table(parameters, file = paste0(path_to_table_outputs, "/RNAi_plausible_quantitative_values_", run_condition, ".csv"), sep = ",", row.names = FALSE)

colnames(parameters_varying_length) <- c("seed", "mRNA", "pre-mRNA", "k1","k2","k3","k4", "k5", "mRNA_after_RNAi", "pre-mRNA_after_RNAi", "length", "mRNA/pre-mRNA_initial", "mRNA/pre-mRNA_after_RNAi", "knockdown")
write.table(parameters_varying_length, file = paste0(path_to_table_outputs, "/RNAi_plausible_quantitative_values_vs_length_", run_condition, ".csv"), sep = ",", row.names = FALSE)

### Set plot themes
my_theme <- theme(panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_line(linetype = 1, size = 0.5, lineend="round"),
                  axis.title = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10),
                  axis.ticks = element_line(size = 0.5), axis.ticks.length = unit(.1, "cm"),
                  axis.text = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10))

## PLOT ALL PARAMETERS for (1) fixed length of mRNA, (2) varying length of mRNA, (3) good knockdown with fixed mRNA length, and (4) poor knockdown with fixed mRNA length.
######### fixed mRNA lengths
## plot pre-mRNA vs mRNA after RNAi as a function of initial mRNA concentration.
max_m_i = round(max(parameters$`mRNA_after_RNAi`))
max_pm_i = round(max(parameters$`pre-mRNA_after_RNAi`))
m_vs_m_i_vs_pm_i <- ggplot(parameters, aes(x = `mRNA_after_RNAi`, y = `pre-mRNA_after_RNAi`, color = mRNA), alpha = 0.5) +
  my_theme + ylab("[pre-mRNA] after RNAi") + xlab("[mRNA] after RNAi") + geom_point(size = 2) +
  scale_color_gradient(low = "darkolivegreen1", high = "darkgreen") +
  labs(color = "initial [mRNA]") +
  scale_x_continuous(limits = c(0, max_m_i), breaks = c(0, 0.5*max_m_i, max_m_i), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_pm_i), breaks= c(0, 0.5*max_pm_i, max_pm_i), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_initial_mRNA_vs_mRNA_and_premRNA_after_RNAi_fixed_mRNA_length_", run_condition, ".eps"), m_vs_m_i_vs_pm_i, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot mRNA to pre-mRNA ratio before and after RNAi as it relates to levels of knockdown.
max_m_p_ratio = max(parameters$`mRNA/pre-mRNA_initial`)
max_m_i_pm_i_ratio = max(parameters$`mRNA/pre-mRNA_after_RNAi`)
m_vs_pm_vs_m_i_vs_pm_i <- ggplot(parameters, aes(x = `mRNA/pre-mRNA_initial`, y = `mRNA/pre-mRNA_after_RNAi`, color = knockdown), alpha = 0.5) +
  my_theme + ylab("[mRNA]/[pre-mRNA] after RNAi") + xlab("[mRNA]/[pre-mRNA] initial") + geom_point(size = 2) +
  scale_color_gradient(low = "cyan", high = "blue") +
  labs(color = "knock down") +
  scale_x_continuous(limits = c(0, 10), breaks = c(0, 0.5*10, 10), expand = c(0.05,0.05)) + ## Note the axis limits here have been hard coded to see trends better.
  scale_y_continuous(limits = c(0, 10), breaks= c(0, 0.5*10, 10), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_initial_mRNA_premRNA_vs_mRNA_premRNA_after_RNAi_fixed_mRNA_length_", run_condition, ".eps"), m_vs_pm_vs_m_i_vs_pm_i, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot constants k1 through k5 in plausible parameter sets for mRNAs of fixed length.
plausible_constants <- parameters[,4:8]
to_plot_constants <- melt(plausible_constants)
constants_jitter <- ggplot(to_plot_constants, aes(variable, value)) +
  my_theme + ylab("value") + xlab("constants") + geom_point(size = 2) + geom_jitter(position = "jitter")

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_plausible_constants_fixed_mRNA_length_", run_condition, ".eps"), constants_jitter, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot knockdown efficiency vs. pre-mRNA/mRNA ratio with fixed mRNA lengths.
max_m_i_pm_i_ratio = max(parameters$`mRNA/pre-mRNA_after_RNAi`)
max_kd_i = max(parameters$knockdown)
KD_vs_m_i_pm_i_ratio <- ggplot(parameters, aes(x = `mRNA/pre-mRNA_after_RNAi`, y = knockdown), alpha = 0.5) +
  my_theme + ylab("(initial [mRNA])/([mRNA] after RNAi)") + xlab("[mRNA]/[pre-mRNA] after RNAi") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_m_i_pm_i_ratio), breaks = c(0, 0.5*max_m_i_pm_i_ratio, max_m_i_pm_i_ratio), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 1), breaks= c(0, 0.5, 1), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_mRNA_pre-mRNA_ratio_vs_KD_fixed_mRNA_length_", run_condition, ".eps"), KD_vs_m_i_pm_i_ratio, height = 3, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot knockdown efficiency vs. amount of initial mRNA for fixed lengths of mRNA.
max_length = max(parameters$mRNA)
max_kd_il = max(parameters$knockdown)
KD_vs_m_i_pm_i_ratio <- ggplot(parameters, aes(x = mRNA, y = knockdown), alpha = 0.5) +
  my_theme + ylab("(initial [mRNA])/([mRNA] after RNAi)") + xlab("[mRNA] before RNAi") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_length), breaks = c(0, 0.5*max_length, max_length), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 1), breaks= c(0, 0.5, 1), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_mRNA_conc_vs_KD_fixed_mRNA_length_", run_condition, ".eps"), KD_vs_m_i_pm_i_ratio, height = 3, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

######### varying mRNA lengths
## plot pre-mRNA vs mRNA after RNAi for mRNAs of varying lengths.
max_m_il = round(max(parameters_varying_length$`mRNA_after_RNAi`))
max_pm_il = round(max(parameters_varying_length$`pre-mRNA_after_RNAi`))
m_vs_m_il_vs_pm_il <- ggplot(parameters_varying_length, aes(x = `mRNA_after_RNAi`, y = `pre-mRNA_after_RNAi`, color = length), alpha = 0.5) +
  my_theme + ylab("[pre-mRNA] after RNAi") + xlab("[mRNA] after RNAi") + geom_point(size = 2) + labs(color = "mRNA length") +
  scale_color_gradient(low = "darkolivegreen1", high = "darkgreen") +
  guides(guide_legend(expression(`gene length`))) +
  scale_x_continuous(limits = c(0, max_m_il), breaks = c(0, 0.5*max_m_il, max_m_il), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_pm_il), breaks= c(0, 0.5*max_pm_il, max_pm_il), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_initial_mRNA_vs_mRNA_and_premRNA_after_RNAi_varying_mRNA_length_", run_condition, ".eps"), m_vs_m_il_vs_pm_il, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot relative lengths of mRNA vs pre-mRNA to mRNA ratios after RNAi.
max_m_il_pm_il_ratio = max(parameters_varying_length$`mRNA/pre-mRNA_after_RNAi`)
max_m_l = max(parameters_varying_length$length)
length_vs_m_il_pm_il_ratio <- ggplot(parameters_varying_length, aes(x = `mRNA/pre-mRNA_after_RNAi`, y = length), alpha = 0.5) +
  my_theme + ylab("relative mRNA length") + xlab("[mRNA]/[pre-mRNA] after RNAi") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_m_il_pm_il_ratio), breaks = c(0, 0.5*max_m_il_pm_il_ratio, max_m_il_pm_il_ratio), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_m_l), breaks= c(0, 0.5*max_m_l, max_m_l), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_mRNA_pre-mRNA_ratio_vs_length_varying_mRNA_length_", run_condition, ".eps"), length_vs_m_il_pm_il_ratio, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot constants k1 through k5 in plausible parameter sets for mRNAs of varying lengths.
plausible_constants <- parameters_varying_length[,4:8]
to_plot_constants <- melt(plausible_constants)
constants_jitter <- ggplot(to_plot_constants, aes(variable, value)) +
  my_theme + ylab("value") + xlab("constants") + geom_point(size = 2) + geom_jitter(position = "jitter")

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_plausible_constants_varying_mRNA_length_", run_condition, ".eps"), constants_jitter, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot knockdown efficiency vs. pre-mRNA/mRNA ratio for mRNAs of varying lengths.
max_m_il_pm_il_ratio = max(parameters_varying_length$`mRNA/pre-mRNA_after_RNAi`)
max_kd_il = max(parameters_varying_length$knockdown)
KD_vs_m_il_pm_il_ratio <- ggplot(parameters_varying_length, aes(x = `mRNA/pre-mRNA_after_RNAi`, y = knockdown), alpha = 0.5) +
  my_theme + ylab("(initial [mRNA])/([mRNA] after RNAi)") + xlab("[mRNA]/[pre-mRNA] after RNAi") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_m_il_pm_il_ratio), breaks = c(0, 0.5*max_m_il_pm_il_ratio, max_m_il_pm_il_ratio), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 1), breaks= c(0, 0.5, 1), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_mRNA_pre-mRNA_ratio_vs_KD_varying_mRNA_length_", run_condition, ".eps"), KD_vs_m_il_pm_il_ratio, height = 3, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot knockdown efficiency vs. length of mRNA.
max_length = max(parameters_varying_length$length)
max_kd_il = max(parameters_varying_length$knockdown)
KD_vs_m_i_pm_i_ratio <- ggplot(parameters_varying_length, aes(x = length, y = knockdown), alpha = 0.5) +
  my_theme + ylab("(initial [mRNA])/([mRNA] after RNAi)") + xlab("relative length of mRNA") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_length), breaks = c(0, 0.5*max_length, max_length), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 1), breaks= c(0, 0.5, 1), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_mRNA_length_vs_KD_", run_condition, ".eps"), KD_vs_m_i_pm_i_ratio, height = 3, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot knockdown efficiency vs. amount of initial mRNA for varying lengths of mRNA.
max_length = max(parameters_varying_length$mRNA)
max_kd_il = max(parameters_varying_length$knockdown)
KD_vs_m_i_pm_i_ratio <- ggplot(parameters_varying_length, aes(x = mRNA, y = knockdown), alpha = 0.5) +
  my_theme + ylab("(initial [mRNA])/([mRNA] after RNAi)") + xlab("[mRNA] before RNAi") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_length), breaks = c(0, 0.5*max_length, max_length), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 1), breaks= c(0, 0.5, 1), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_mRNA_conc_vs_KD_varying_mRNA_length_", run_condition, ".eps"), KD_vs_m_i_pm_i_ratio, height = 3, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot mRNA to pre-mRNA ratio before and after RNAi as it relates to levels of knockdown and varying lengths of mRNA.
max_m_p_ratio = max(parameters_varying_length$`mRNA/pre-mRNA_initial`)
max_m_il_pm_il_ratio = max(parameters_varying_length$`mRNA/pre-mRNA_after_RNAi`)
m_vs_pm_vs_m_i_vs_pm_i <- ggplot(parameters_varying_length, aes(x = `mRNA/pre-mRNA_initial`, y = `mRNA/pre-mRNA_after_RNAi`, color = knockdown), alpha = 0.5) +
  my_theme + ylab("[mRNA]/[pre-mRNA] after RNAi") + xlab("[mRNA]/[pre-mRNA] initial") + geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(color = "knock down") +
  scale_x_continuous(limits = c(0, 10), breaks = c(0, 0.5*10, 10), expand = c(0.05,0.05)) + ## Note the axis limits here have been hard coded to see trends better.
  scale_y_continuous(limits = c(0, 10), breaks= c(0, 0.5*10, 10), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_initial_mRNA_premRNA_vs_mRNA_premRNA_after_RNAi_varying_mRNA_length_", run_condition, ".eps"), m_vs_pm_vs_m_i_vs_pm_i, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

### Plot each constant against all other constants with fixed length.
### k1 vs k2, k1 vs k3, k1 vs k4, k1 vs k5, k2 vs k3, k2 vs k4, k2 vs k5, k3 vs k4, k3 vs k5, and k4 vs k5
# k1 vs k2
max_k1 = round(max(parameters$k1), 1)
max_k2 = round(max(parameters$k2), 1)
constants_vs_constants_fixed <- ggplot(parameters, aes(x = k1, y = k2), alpha = 0.5) +
  my_theme + ylab("k2") + xlab("k1") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k1), breaks = c(0, 0.5*max_k1, max_k1), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k2), breaks= c(0, 0.5*max_k2, max_k2), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_fixed_gene_length_k1_vs_k2_", run_condition, ".eps"), constants_vs_constants_fixed, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k1 vs k3
max_k1 = round(max(parameters$k1), 1)
max_k3 = round(max(parameters$k3), 1)
constants_vs_constants_fixed <- ggplot(parameters, aes(x = k1, y = k3), alpha = 0.5) +
  my_theme + ylab("k3") + xlab("k1") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k1), breaks = c(0, 0.5*max_k1, max_k1), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k3), breaks= c(0, 0.5*max_k3, max_k3), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_fixed_gene_length_k1_vs_k3_", run_condition, ".eps"), constants_vs_constants_fixed, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k1 vs k4
max_k1 = round(max(parameters$k1), 1)
max_k4 = round(max(parameters$k4), 1)
constants_vs_constants_fixed <- ggplot(parameters, aes(x = k1, y = k4), alpha = 0.5) +
  my_theme + ylab("k4") + xlab("k1") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k1), breaks = c(0, 0.5*max_k1, max_k1), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k4), breaks= c(0, 0.5*max_k4, max_k4), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_fixed_gene_length_k1_vs_k4_", run_condition, ".eps"), constants_vs_constants_fixed, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k1 vs k5
max_k1 = round(max(parameters$k1), 1)
max_k5 = round(max(parameters$k5), 1)
constants_vs_constants_fixed <- ggplot(parameters, aes(x = k1, y = k5), alpha = 0.5) +
  my_theme + ylab("k5") + xlab("k1") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k1), breaks = c(0, 0.5*max_k1, max_k1), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k5), breaks= c(0, 0.5*max_k5, max_k5), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_fixed_gene_length_k1_vs_k5_", run_condition, ".eps"), constants_vs_constants_fixed, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k2 vs k3
max_k2 = round(max(parameters$k2), 1)
max_k3 = round(max(parameters$k3), 1)
constants_vs_constants_fixed <- ggplot(parameters, aes(x = k2, y = k3), alpha = 0.5) +
  my_theme + ylab("k3") + xlab("k2") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k2), breaks = c(0, 0.5*max_k2, max_k2), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k3), breaks= c(0, 0.5*max_k3, max_k3), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_fixed_gene_length_k2_vs_k3_", run_condition, ".eps"), constants_vs_constants_fixed, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k2 vs k4
max_k2 = round(max(parameters$k2), 1)
max_k4 = round(max(parameters$k4), 1)
constants_vs_constants_fixed <- ggplot(parameters, aes(x = k2, y = k4), alpha = 0.5) +
  my_theme + ylab("k4") + xlab("k2") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k2), breaks = c(0, 0.5*max_k2, max_k2), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k4), breaks= c(0, 0.5*max_k4, max_k4), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_fixed_gene_length_k2_vs_k4_", run_condition, ".eps"), constants_vs_constants_fixed, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k2 vs k5
max_k2 = round(max(parameters$k2), 1)
max_k5 = round(max(parameters$k5), 1)
constants_vs_constants_fixed <- ggplot(parameters, aes(x = k2, y = k5), alpha = 0.5) +
  my_theme + ylab("k5") + xlab("k2") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k2), breaks = c(0, 0.5*max_k2, max_k2), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k5), breaks= c(0, 0.5*max_k5, max_k5), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_fixed_gene_length_k2_vs_k5_", run_condition, ".eps"), constants_vs_constants_fixed, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k3 vs k4
max_k3 = round(max(parameters$k3), 1)
max_k4 = round(max(parameters$k4), 1)
constants_vs_constants_fixed <- ggplot(parameters, aes(x = k3, y = k4), alpha = 0.5) +
  my_theme + ylab("k4") + xlab("k3") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k3), breaks = c(0, 0.5*max_k3, max_k3), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k4), breaks= c(0, 0.5*max_k4, max_k4), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_fixed_gene_length_k3_vs_k4_", run_condition, ".eps"), constants_vs_constants_fixed, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k3 vs k5
max_k3 = round(max(parameters$k3), 1)
max_k5 = round(max(parameters$k5), 1)
constants_vs_constants_fixed <- ggplot(parameters, aes(x = k3, y = k5), alpha = 0.5) +
  my_theme + ylab("k5") + xlab("k3") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k3), breaks = c(0, 0.5*max_k3, max_k3), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k5), breaks= c(0, 0.5*max_k5, max_k5), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_fixed_gene_length_k3_vs_k5_", run_condition, ".eps"), constants_vs_constants_fixed, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k4 vs k5
max_k4 = round(max(parameters$k4), 1)
max_k5 = round(max(parameters$k5), 1)
constants_vs_constants_fixed <- ggplot(parameters, aes(x = k4, y = k5), alpha = 0.5) +
  my_theme + ylab("k5") + xlab("k4") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k4), breaks = c(0, 0.5*max_k4, max_k4), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k5), breaks= c(0, 0.5*max_k5, max_k5), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_fixed_gene_length_k4_vs_k5_", run_condition, ".eps"), constants_vs_constants_fixed, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

### Plot each constant against all other constants with varying length.
### k1 vs k2, k1 vs k3, k1 vs k4, k1 vs k5, k2 vs k3, k2 vs k4, k2 vs k5, k3 vs k4, k3 vs k5, and k4 vs k5
# k1 vs k2
max_k1 = round(max(parameters_varying_length$k1), 1)
max_k2 = round(max(parameters_varying_length$k2), 1)
constants_vs_constants_varying <- ggplot(parameters_varying_length, aes(x = k1, y = k2), alpha = 0.5) +
  my_theme + ylab("k2") + xlab("k1") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k1), breaks = c(0, 0.5*max_k1, max_k1), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k2), breaks= c(0, 0.5*max_k2, max_k2), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_varying_gene_length_k1_vs_k2_", run_condition, ".eps"), constants_vs_constants_varying, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k1 vs k3
max_k1 = round(max(parameters_varying_length$k1), 1)
max_k3 = round(max(parameters_varying_length$k3), 1)
constants_vs_constants_varying <- ggplot(parameters_varying_length, aes(x = k1, y = k3), alpha = 0.5) +
  my_theme + ylab("k3") + xlab("k1") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k1), breaks = c(0, 0.5*max_k1, max_k1), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k3), breaks= c(0, 0.5*max_k3, max_k3), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_varying_gene_length_k1_vs_k3_", run_condition, ".eps"), constants_vs_constants_varying, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k1 vs k4
max_k1 = round(max(parameters_varying_length$k1), 1)
max_k4 = round(max(parameters_varying_length$k4), 1)
constants_vs_constants_varying <- ggplot(parameters_varying_length, aes(x = k1, y = k4), alpha = 0.5) +
  my_theme + ylab("k4") + xlab("k1") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k1), breaks = c(0, 0.5*max_k1, max_k1), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k4), breaks= c(0, 0.5*max_k4, max_k4), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_varying_gene_length_k1_vs_k4_", run_condition, ".eps"), constants_vs_constants_varying, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k1 vs k5
max_k1 = round(max(parameters_varying_length$k1), 1)
max_k5 = round(max(parameters_varying_length$k5), 1)
constants_vs_constants_varying <- ggplot(parameters_varying_length, aes(x = k1, y = k5), alpha = 0.5) +
  my_theme + ylab("k5") + xlab("k1") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k1), breaks = c(0, 0.5*max_k1, max_k1), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k5), breaks= c(0, 0.5*max_k5, max_k5), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_varying_gene_length_k1_vs_k5_", run_condition, ".eps"), constants_vs_constants_varying, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k2 vs k3
max_k2 = round(max(parameters_varying_length$k2), 1)
max_k3 = round(max(parameters_varying_length$k3), 1)
constants_vs_constants_varying <- ggplot(parameters_varying_length, aes(x = k2, y = k3), alpha = 0.5) +
  my_theme + ylab("k3") + xlab("k2") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k2), breaks = c(0, 0.5*max_k2, max_k2), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k3), breaks= c(0, 0.5*max_k3, max_k3), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_varying_gene_length_k2_vs_k3_", run_condition, ".eps"), constants_vs_constants_varying, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k2 vs k4
max_k2 = round(max(parameters_varying_length$k2), 1)
max_k4 = round(max(parameters_varying_length$k4), 1)
constants_vs_constants_varying <- ggplot(parameters_varying_length, aes(x = k2, y = k4), alpha = 0.5) +
  my_theme + ylab("k4") + xlab("k2") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k2), breaks = c(0, 0.5*max_k2, max_k2), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k4), breaks= c(0, 0.5*max_k4, max_k4), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_varying_gene_length_k2_vs_k4_", run_condition, ".eps"), constants_vs_constants_varying, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k2 vs k5
max_k2 = round(max(parameters_varying_length$k2), 1)
max_k5 = round(max(parameters_varying_length$k5), 1)
constants_vs_constants_varying <- ggplot(parameters_varying_length, aes(x = k2, y = k5), alpha = 0.5) +
  my_theme + ylab("k5") + xlab("k2") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k2), breaks = c(0, 0.5*max_k2, max_k2), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k5), breaks= c(0, 0.5*max_k5, max_k5), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_varying_gene_length_k2_vs_k5_", run_condition, ".eps"), constants_vs_constants_varying, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k3 vs k4
max_k3 = round(max(parameters_varying_length$k3), 1)
max_k4 = round(max(parameters_varying_length$k4), 1)
constants_vs_constants_varying <- ggplot(parameters_varying_length, aes(x = k3, y = k4), alpha = 0.5) +
  my_theme + ylab("k4") + xlab("k3") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k3), breaks = c(0, 0.5*max_k3, max_k3), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k4), breaks= c(0, 0.5*max_k4, max_k4), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_varying_gene_length_k3_vs_k4_", run_condition, ".eps"), constants_vs_constants_varying, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k3 vs k5
max_k3 = round(max(parameters_varying_length$k3), 1)
max_k5 = round(max(parameters_varying_length$k5), 1)
constants_vs_constants_varying <- ggplot(parameters_varying_length, aes(x = k3, y = k5), alpha = 0.5) +
  my_theme + ylab("k5") + xlab("k3") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k3), breaks = c(0, 0.5*max_k3, max_k3), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k5), breaks= c(0, 0.5*max_k5, max_k5), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_varying_gene_length_k3_vs_k5_", run_condition, ".eps"), constants_vs_constants_varying, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

# k4 vs k5
max_k4 = round(max(parameters_varying_length$k4), 1)
max_k5 = round(max(parameters_varying_length$k5), 1)
constants_vs_constants_varying <- ggplot(parameters_varying_length, aes(x = k4, y = k5), alpha = 0.5) +
  my_theme + ylab("k5") + xlab("k4") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_k4), breaks = c(0, 0.5*max_k4, max_k4), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, max_k5), breaks= c(0, 0.5*max_k5, max_k5), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_varying_gene_length_k4_vs_k5_", run_condition, ".eps"), constants_vs_constants_varying, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## Subscaled plotting of k2 vs k4 for fixed mRNA lengths to show exclusion zone of parameters.
constants_vs_constants_fixed <- ggplot(parameters, aes(x = k2, y = k4), alpha = 0.5) +
  my_theme + ylab("k4") + xlab("k2") + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.05,0.05)) +
  scale_y_continuous(limits = c(0, 1), breaks= c(0, 0.5, 1), expand = c(0.05,0.05))

ggsave(file = paste0(path_to_figure_outputs, "/RNAi_constants_vs_constants_fixed_gene_length_k2_vs_k4_LOW_VALUES_", run_condition, ".eps"), constants_vs_constants_fixed, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)
