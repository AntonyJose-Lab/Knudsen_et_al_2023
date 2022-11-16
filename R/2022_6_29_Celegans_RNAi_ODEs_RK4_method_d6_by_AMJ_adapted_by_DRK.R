### This is a program for using Runge Kutta order 4 method to solve RNAi pathways written as Ordinary Differential Equations.
### This program considers all binding events as irreversible.
### Arbitrary parameter values have been given for all constants. Adapted to R by Daphne Knudsen from program in Python by Antony Jose
###

#Load required packages
library(ggplot2)
library(dplyr)


#Characteristics of gene
lm = 10 # length of mRNA in kb

#time step for the simulation and duration
Dt = 0.01       # Time step delta t
Dt2 = Dt/2
t0 = 0          # start time
tn = 500       # end time
n_steps = as.integer(round(((tn-t0)/Dt), digits = 0))   # number of time steps

# Set starting conditions of the simulation
ds0 = 10        # Double-stranded RNA avvailble for silencing
pri0 = 0        # Initial amount of primary siRNA
ug0 = 0         # Initial amount of poly-UG RNA
sec0 = 0        # Initial amount of secondary siRNA
p0 = 50         # Initial amount of pre-mRNA
m0 = 100        # Initial amount of mRNA


#Set values for all parameters
lds = 220       #length of dsRNA in base pairs
Tm = 0.05       # faction of mRNA degraded by other processes per time step
Tp = 0.05       # fraction of pre-mRNA degraded by other processes per time step
Tug = 0.05      # fraction of pUG RNA degraded per time step
Tpri = 0.05     # fraction of 1º siRNA degraded per time step
Tsec = 0.05     # fraction of 2º siRNA degraded per time step

k1 = 1     # 1 dsRNA imported and processed into siRNA / time steps.
k2 = 0.01   # binding constant of 1º siRNAs with mRNA.
k3 = 1       # 1 pUG RNA made per 1º siRNA-bound mRNAs.
k4 = 0.05*lm  # 1 2º siRNA made using every 10 pUG RNA and more made from longer RNA.
k5 = 0.01   # binding constant of 2º siRNA with mRNA.
k6 = 0.01   # binding constant of 2º siRNA with pre-mRNA.
k7 = Tm*m0/p0 # proportion of spliced and exported transcripts per time step.
k8 = 0.05    # repression of transcription by secondary siRNAs in the nucleus.
k9 = Tp*p0 + Tm*m0 # Number of transcripts transcribed per time step.

## base values
#k1 = 1     # 1 dsRNA imported and processed into siRNA / time steps.
#k2 = 0.01   # binding constant of 1º siRNAs with mRNA.
#k3 = 1       # 1 pUG RNA made per 1º siRNA-bound mRNAs.
#k4 = 0.05*lm  # 1 2º siRNA made using every 10 pUG RNA and more made from longer RNA.
#k5 = 0.01   # binding constant of 2º siRNA with mRNA.
#k6 = 0.01   # binding constant of 2º siRNA with pre-mRNA.
#k7 = Tm*m0/p0 # proportion of spliced and exported transcripts per time step.
#k8 = 0.05    # repression of transcription by secondary siRNAs in the nucleus.
#k9 = Tp*p0 + Tm*m0

#Create arrays with zeros for storing the values of all intermediates throughout the simulation.
ds_exact_arr = rep(0, (n_steps +1))
ds_arr = rep(0, (n_steps +1))
pri_arr = rep(0, (n_steps +1))
ug_arr = rep(0, (n_steps +1))
sec_arr = rep(0, (n_steps +1))
p_arr = rep(0, (n_steps +1))
m_arr = rep(0, (n_steps +1))
t_arr = rep(0, (n_steps +1))

## Add initial values of all variables to 1st element of array.
ds_arr[1] = ds0
pri_arr[1] = pri0
ug_arr[1] = ug0
sec_arr[1] = sec0
p_arr[1] = p0
m_arr[1] = m0
t_arr[1] = t0



# initialize all variables.
##DK: since R index values start at 1 and not 0, the range for the loop starts at index 2 and proceeds to x_arr[2-1]
for(i in 2:(n_steps + 1)) {
  ds = ds_arr[i-1]
  pri = pri_arr[i-1]
  ug = ug_arr[i-1]
  sec = sec_arr[i-1]
  p = p_arr[i-1]
  m = m_arr[i-1]
  t = t_arr[i-1]

# Calculate del for all
  ds_del = -k1*ds
  pri_del = (-ds_del)*lds/22 - k2*pri*m - Tpri*pri
  ug_del = k3*k2*pri*m - Tug*ug
  sec_del = k4*ug - k5*m*sec - k6*p*sec - Tsec*sec
  p_del = k9*(1 - k8*k6*p*sec) - k7*p - Tp*p
  m_del = k7*p - k2*m*pri - k5*m*sec - Tm*m

# Calculate c1, the first slope approximation
  ds_c1 = ds_del
  pri_c1 = pri_del
  ug_c1 = ug_del
  sec_c1 = sec_del
  p_c1 = p_del
  m_c1 = m_del

# Calculate c2, the second slope approximation
  ds_c2 = Dt2*ds_c1
  pri_c2 = Dt2*pri_c1
  ug_c2 = Dt2*ug_c1
  sec_c2 = Dt2*sec_c1
  p_c2 = Dt2*p_c1
  m_c2 = Dt2*m_c1

# Calculate c3, the third slope approximation
  ds_c3 = Dt2*ds_c2
  pri_c3 = Dt2*pri_c2
  ug_c3 = Dt2*ug_c2
  sec_c3 = Dt2*sec_c2
  p_c3 = Dt2*p_c2
  m_c3 = Dt2*m_c2

# Calculate c4, the fourth slope approximation
  ds_c4 = Dt*ds_c3
  pri_c4 = Dt*pri_c3
  ug_c4 = Dt*ug_c3
  sec_c4 = Dt*sec_c3
  p_c4 = Dt*p_c3
  m_c4 = Dt*m_c3

# Update all arrays after RK4 calculations.
  ds_arr[i] = ds + (Dt/6)*(ds_c1 + 2*ds_c2 + 2*ds_c3 + ds_c4)
  pri_arr[i] = pri + (Dt/6)*(pri_c1 + 2*pri_c2 + 2*pri_c3 + pri_c4)
  ug_arr[i] = ug + (Dt/6)*(ug_c1 + 2*ug_c2 + 2*ug_c3 + ug_c4)
  sec_arr[i] = sec + (Dt/6)*(sec_c1 + 2*sec_c2 + 2*sec_c3 + sec_c4)
  p_arr[i] = p + (Dt/6)*(p_c1 + 2*p_c2 + 2*p_c3 + p_c4)
  m_arr[i] = m + (Dt/6)*(m_c1 + 2*m_c2 + 2*m_c3 + m_c4)
  t_arr[i] = t + Dt
}



#Formatting theme
my_theme <- theme(panel.background = element_blank(),

                  ##set blank background for facet label region
                  strip.background = element_blank(),

                  ##set blank grid
                  panel.grid = element_blank(),

                  ##specify line format
                  axis.line = element_line(linetype = 1, size = 1, lineend = "round"),

                  ##specify axis title font
                  axis.title = element_text(family = "Helvetica", face = "plain", color = "grey0", size = 15),

                  ##set axis tick thickness and length
                  axis.ticks = element_line(size = 0.9, color = "grey0"), axis.ticks.length = unit(.2, "cm"),

                  ##set axis text font
                  axis.text = element_text(family = "Helvetica", face = "plain", color = "grey0", size = 15),

                  ##set legend font size
                  legend.text = element_text(size = 15),

                  ##set facet label font size
                  strip.text = element_text(size = 15))

#plotting
ggplot()+
  my_theme+
 # ylim(0,250)+
  geom_line(aes(x=t_arr, y=sec_arr), color="lightsalmon3", size=2)+
  geom_line(aes(x=t_arr, y=ug_arr), color="green3", size=2)+
  geom_line(aes(x=t_arr, y=ds_arr), color="red", size=2)+
  geom_line(aes(x=t_arr, y=m_arr), color="dodgerblue2", size=2)+
  geom_line(aes(x=t_arr, y=pri_arr), color="mediumpurple", size=2)+
  geom_line(aes(x=t_arr, y=p_arr), color="goldenrod2", size=2) +
  xlab("time (a.u.)")+
  ylab("[RNA]t")+
  scale_color_continuous(labels = c("2º siRNA", "poly-UG RNA", "dsRNA", "mRNA", "1º RNA", "pre-mRNA"))
