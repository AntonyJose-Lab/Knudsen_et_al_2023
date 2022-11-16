# -*- coding: utf-8 -*-
"""
This is a program for using Runge Kutta order 4 method to solve RNAi pathways written as Ordinary Differential Equations.
All binding events are considered as irreversible and arbitrary parameter values have been given for all constants. This program generates an illustrative plot.
Author: Antony Jose
"""
## Import required libraries.
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import os

# Directory for saving files, changing into, and setting working directory.
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/RNA_silencing/analyses/Python/Figures', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/RNA_silencing/analyses/Python/Tables', exist_ok=True)
os.chdir('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/RNA_silencing/analyses/Python/')
curr_dir = os.getcwd()

## time step for the simulation and duration.
Dt = 0.01           # Time step delta t
Dt2 = Dt/2
t0 = 0              # start time
tn = 500           # end time
n_steps = int(round((tn-t0)/Dt))    # number of time steps

## Set starting conditions of the simulation
ds0 = 10            # Double-stranded RNA available for silencing
pri0 = 0            # Initial amount of primary siRNA
ug0 = 0             # Initial amount of poly UG RNA
sec0 = 0            # Initial amount of secondary siRNA
p0 = 50             # Initial amount of pre-mRNA
m0 = 100            # initial amount of mRNA

## Set values for all parameters.
lds = 220   # length of dsRNA in base pairs
lm = 1      # length of mRNA in kb

Tm = 0.05   # fraction of mRNA degraded by other processes per time step
Tp = 0.05   # fraction of pre-mRNA degraded by other processes per time step.
Tug = 0.05  # fraction of pUG RNA degraded per time step.
Tpri = 0.05 # fraction of 1º siRNA degraded per time step.
Tsec = 0.05 # fraction of 2º siRNA degraded per time step.

k1 = 1     # 1 dsRNA imported and processed into siRNA / time steps.
k2 = 0.01   # binding constant of 1º siRNAs with mRNA.
k3 = 1       # 1 pUG RNA made per 1º siRNA-bound mRNAs.
k4 = 0.5*lm  # 1 2º siRNA made using every 10 pUG RNA and more made from longer RNAs.
k5 = 0.01   # binding constant of 2º siRNA with mRNA.
k6 = 0.01   # binding constant of 2º siRNA with pre-mRNA.
k7 = Tm*m0/p0 # proportion of spliced and exported transcripts per time step.
k8 = 0.05    # repression of transcription by secondary siRNAs in the nucleus.
k9 = Tp*p0 + Tm*m0 # Number of transcripts transcribed per time step.

## Create arrays with zeros for storing the values of all intermediates throughout the simulation.
ds_exact_arr = np.zeros(n_steps + 1)
ds_arr = np.zeros(n_steps + 1)
pri_arr = np.zeros(n_steps + 1)
ug_arr = np.zeros(n_steps + 1)
sec_arr = np.zeros(n_steps + 1)
p_arr = np.zeros(n_steps + 1)
m_arr = np.zeros(n_steps + 1)
t_arr = np.zeros(n_steps + 1)

## Add initial values of all variables to 1st element of array.
ds_arr[0] = ds0
pri_arr[0] = pri0
ug_arr[0] = ug0
sec_arr[0] = sec0
p_arr[0] = p0
m_arr[0] = m0
t_arr[0] = t0

for i in range (1, n_steps + 1):
   # initialize all variables.
   ds = ds_arr[i-1]
   pri = pri_arr[i-1]
   ug = ug_arr[i-1]
   sec = sec_arr[i-1]
   p = p_arr[i-1]
   m = m_arr[i-1]
   t = t_arr[i-1]

   # Calculate increment per time-step for all
   ds_del = -k1*ds
   pri_del = (-ds_del)*lds/22 - k2*pri*m - Tpri*pri
   ug_del = k3*k2*pri*m - Tug*ug
   sec_del = k4*lm*ug - k5*m*sec - k6*p*sec - Tsec*sec
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


# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(t_arr, m_arr, linewidth = 1, label = 'mRNA')
plt.plot(t_arr, p_arr, linewidth = 1, label = 'pre-mRNA')
plt.plot(t_arr, ug_arr, linewidth = 1, label = 'pUG RNA')
plt.plot(t_arr, ds_arr, linewidth = 1, label = 'dsRNA')
plt.plot(t_arr, pri_arr, linewidth = 1, label = '1º siRNA')
plt.plot(t_arr, sec_arr, linewidth = 1, label = '2º siRNA')
plt.title('The dynamics of RNA interference in C. elegans', fontsize = 12)
plt.xlabel('t (hours)', fontsize = 8)
plt.ylabel('molecules(t)', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)                 # do not show grid
ymax = max(max(m_arr), max(p_arr), max(ug_arr), max(ds_arr), max(pri_arr), max(sec_arr)) + 10
plt.axis([t0, tn, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_7_14_RNAiDynamics_arbitary_parameters'
path = os.path.join(curr_dir, 'Figures/2022_7_14_RNAiDynamics/', file_name + '.svg')
fig.savefig(path, dpi=300)
