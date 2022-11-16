## This program aims to simulate various amounts of 'silencing signals' requiring mut-16, rde-10, nrde-3, or other Argonautes.
## Sets of values that satisfy the experimentally observed dependence for silencing suggest plausible models.

## Definitions:
# Nm = amount of silencing signals requiring MUT-16 and NRDE-3
# Om = amount of silencing signals requiring MUT-16 and Other Argonautes
# Nr = amount of silencing signals requiring RDE-10 and NRDE-3
# Or = amount of silencing signals requiring RDE-10 and Other Argonautes
# T_unc22 = threshold of silencing signals required for observing loss of unc-22 function, which accounts for target-specific factors (length, sequence bias, accessibility, etc.)
# T_bli1 = threshold of silencing signals required for observing loss of bli-1 function, which accounts for target-specific factors (length, sequence bias, accessibility, etc.)

## Experimental Conditions:
# (Nm + Om) > T_unc22 ; loss of rde-10 alone does not prevent unc-22 silencing.
# (Nr + Or) > T_unc22 ; loss of mut-16 alone does not prevent unc-22 silencing.
# (Or + Om) > T_unc22 ; loss of nrde-3 alone does not prevent unc-22 silencing.
# Or < T_unc22 ; loss of both mut-16 and nrde-3 prevents unc-22 silencing.
# Om < T_unc22 ; loss of both rde-10 and nrde-3 prevents unc-22 silencing.
# (Nm + Om) < T_bli1 ; loss of rde-10 prevents bli-1 silencing.
# (Or + Om) < T_bli1 ; loss of nrde-3 prevents bli-1 silencing.
# (Nr + Or) < T_bli1 ; loss of mut-16 prevents bli-1 silencing.
# (Nm + Nr + Or + Om) > T_bli ; wild-type animals can silence bli-1.

## Question: Are there values for all 6 variables that satisfy all 9 conditions?
## If we find possible values, then a quantitative model based upon fixed relative amounts of silencing signals from different routes can account for the observations.
## If not, the silencing of different genes could rely on different RNA response networks.

## Approach: Set T_unc22 to 1 and determine the values of all other variables with respect to that: 5 variables and 9 conditions.

experiment <- "2022_2_9_RNAi_network_thresholds_simpler"

root <- "/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/RNA_silencing/"
setwd(paste0(root, "code/R"))
dir.create(paste0(root, "analyses/R/Figures/", experiment))
dir.create(paste0(root, "analyses/R/Tables/", experiment))
path_to_figure_outputs <- paste0(root, "analyses/R/Figures/", experiment)
path_to_table_outputs <- paste0(root, "analyses/R/Tables/", experiment)

## run conditions for varying number of experiments to be considered. Comment out different conditions below and reassign variables accordingly to explore how each experimental result constrains values.
run_condition <- "unc22_and_bli1"

library(tidyverse)

parameters <- data.frame()
reduced_parameters <- data.frame()

for (s in 1:100) {
  set.seed(s) ## Set a number of different random seeds to get reproducible code and yet explore a range of random numbers.
  T_unc22 = 1 ## This makes T_bli = T_bli / T_unc22.
  for (i in 1:1000) {
  Nm = runif (1, 0, 2) ## Can be arbitrarily high or low and contributes to setting the value for T_bli1. [i.e., fr. of nAgo-dep. process is indeterminate]
  Om = runif (1, 0, 2) ## Has to be less than 1 (i.e., T_unc22)
  Nr = runif (1, 0, 2) ## Can be arbitrarily high or low and contributes to setting the value for T_bli1. [i.e., fr. of nAgo-dep. process is indeterminate]
  Or = runif (1, 0, 2) ## Has to be less than 1 (i.e., T_unc22).
  T_bli1 = runif (1, 0, 100)  ## T_bli1 > T_unc22, given our experiments. Consistently, no systems with T_bli1 < 1 satisfies all conditions.
    if ((Nm + Om) > T_unc22 &&
       (Nr + Or) > T_unc22 &&
       (Or + Om) > T_unc22 &&
       Om < T_unc22 &&
       Or < T_unc22 &&
       (Nm + Om) < T_bli1 &&
       (Or + Om) < T_bli1 &&
       (Nr + Or) < T_bli1 &&
       (Nm + Nr + Om + Or) > T_bli1
        ) {
      parameters <- rbind(parameters, cbind(s, i, T_unc22, T_bli1, Nm, Om, Nr, Or))
      mut16_to_nrde3 <- Nm / (Nm + Om) ## fraction MUT-16-dependent silencing signals that require NRDE-3
      rde10_to_nrde3 <- Nr / (Nr + Or) ## fraction RDE-10-dependent silencing signals that require NRDE-3
      rel_nrde3 <- (Nm + Nr) / (Nm + Nr + Om + Or) ## relative contribution of NRDE-3
      rel_mut16 <- (Nm + Om) / (Nm + Nr + Om + Or) ## relative contribution of MUT-16
      rel_rde10 <- (Or + Nr) / (Nm + Nr + Om + Or) ## relative contribution of RDE-10
      rel_oAgo <- (Or + Om) / (Nm + Nr + Om + Or) ## relative contribution of other Agos
      reduced_parameters <- rbind(reduced_parameters, cbind(T_bli1, mut16_to_nrde3, rde10_to_nrde3, rel_nrde3, rel_mut16, rel_rde10, rel_oAgo))

colnames(parameters) <- c("seed", "test", "T_unc22", "T_bli1","Nm","Om","Nr","Or")
colnames(reduced_parameters) <- c("T_bli1", "mut16_to_nrde3", "rde10_to_nrde3", "rel_nrde3", "rel_mut16", "rel_rde10", "rel_oAgo")
write.table(parameters, file = paste0(path_to_table_outputs, "/RNAi_plausible_quantitative_values_", run_condition, ".csv"), sep = ",", row.names = FALSE)
write.table(reduced_parameters, file = paste0(path_to_table_outputs, "/RNAi_plausible_quantitative_values-relative_", run_condition, ".csv"), sep = ",", row.names = FALSE)

### Set plot themes
my_theme <- theme(panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_line(linetype = 1, size = 0.5, lineend="round"),
                  axis.title = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10),
                  axis.ticks = element_line(size = 0.5), axis.ticks.length = unit(.1, "cm"),
                  axis.text = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10))

## plot T_bli1 vs mut16_to_nrde3 vs rde10_vs_nrde3
rde10_vs_mut16_to_nrde3 <- ggplot(reduced_parameters, aes(x = mut16_to_nrde3, y = rde10_to_nrde3, color = T_bli1), alpha = 0.5) +
  my_theme + ylab("RDE-10 to NRDE-3") + xlab("MUT-16 to NRDE-3") + geom_point(size = 2) +
  scale_color_gradient(low = "cyan", high = "red") +
  guides(guide_legend(expression(T[italic(bli-1)]))) +
  scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.05), breaks= c(0, 0.5, 1), expand = c(0,0))

ggsave(file = paste0(path_to_figure_outputs, "/RDE-10_vs_MUT-16_to_NRDE-3_", run_condition, ".eps"), rde10_vs_mut16_to_nrde3, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

max_T_bli1 <- round(max(reduced_parameters$T_bli1), 1)

## plot T_bli1 vs rel_nrde3
relative_nrde3 <- ggplot(reduced_parameters, aes(x = T_bli1, y = rel_nrde3), alpha = 0.5) +
  my_theme + ylab("relative NRDE-3") + xlab(expression(T[italic(bli-1)])) + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_T_bli1*1.05), breaks = c(0, max_T_bli1/2, max_T_bli1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.05), breaks= c(0, 0.5, 1), expand = c(0,0))

ggsave(file = paste0(path_to_figure_outputs, "/Relative_NRDE-3_", run_condition, ".eps"), relative_nrde3, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot T_bli1 vs rel_mut16
relative_mut16 <- ggplot(reduced_parameters, aes(x = T_bli1, y = rel_mut16), alpha = 0.5) +
  my_theme + ylab("relative MUT-16") + xlab(expression(T[italic(bli-1)])) + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_T_bli1*1.05), breaks = c(0, max_T_bli1/2, max_T_bli1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.05), breaks= c(0, 0.5, 1), expand = c(0,0))

ggsave(file = paste0(path_to_figure_outputs, "/Relative_MUT-16_", run_condition, ".eps"), relative_mut16, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot T_bli1 vs rel_rde10
relative_rde10 <- ggplot(reduced_parameters, aes(x = T_bli1, y = rel_rde10), alpha = 0.5) +
  my_theme + ylab("relative RDE-10") + xlab(expression(T[italic(bli-1)])) + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_T_bli1*1.05), breaks = c(0, max_T_bli1/2, max_T_bli1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.05), breaks= c(0, 0.5, 1), expand = c(0,0))

ggsave(file = paste0(path_to_figure_outputs, "/Relative_RDE-10_", run_condition, ".eps"), relative_rde10, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)

## plot T_bli1 vs rel_oAgo
relative_oAgo <- ggplot(reduced_parameters, aes(x = T_bli1, y = rel_oAgo), alpha = 0.5) +
  my_theme + ylab("relative other Ago") + xlab(expression(T[italic(bli-1)])) + geom_point(size = 2) +
  scale_x_continuous(limits = c(0, max_T_bli1*1.05), breaks = c(0, max_T_bli1/2, max_T_bli1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.05), breaks= c(0, 0.5, 1), expand = c(0,0))

ggsave(file = paste0(path_to_figure_outputs, "/Relative_other_Ago_", run_condition, ".eps"), relative_oAgo, height = 2.5, width = 3.5, device = "eps", units = "in", dpi = 300, limitsize = FALSE)
