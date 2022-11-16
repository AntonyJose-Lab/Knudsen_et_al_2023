# -*- coding: utf-8 -*-
"""
This is a program for using Runge Kutta order 4 method to solve RNAi pathways written as Ordinary Differential Equations.
All binding events are considered as irreversible and arbitrary parameter values have been given for all constants.
The steady state assumption for mRNA levels has been relaxed and the change in time to knockdown and duration of knockdown in response to doubling the rate of transcription is calculated.
Author: Antony Jose
"""
## Import required libraries.
import numpy as np
import math
import matplotlib.pyplot as plt
import os
import csv

# Directory for saving files, changing into, and setting working directory.
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/RNA_silencing/analyses/Python/Figures', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/RNA_silencing/analyses/Python/Tables', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/RNA_silencing/analyses/Python/Tables/2022_7_27_RNAiDynamics', exist_ok=True)
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

## Writing out the parameters from the simulation as a .csv file.
f = open('Tables/2022_7_27_RNAiDynamics/RNAiDynamics_non_SS_param_sweep.csv', 'w', encoding='UTF8', newline='')
writer = csv.writer(f)
header = ['lm','lds','Tm','Tp','Tug','Tpri','Tsec','k1','k2','k3','k4','k5','k6','k8','k9','kd_time','time_to_kd','max_sec', 'min_m']
writer.writerow(header)

for lm in (1, 10) :
    for lds in (220, 440):
        for Tm in (0.05, 0.005):
            for Tp in (0.05, 0.005):
                for Tug in (0.05, 0.005):
                    for Tpri in (0.05, 0.005):
                        for Tsec in (0.05, 0.005):
                            for k1 in (0.1, 1):
                                for k2 in (0.01, 0.1):
                                    for k3 in (0.1, 1):
                                        for k4 in (0.1, 1):
                                            for k5 in (0.01, 0.1):
                                                for k6 in (0.01, 0.1):
                                                    for k8 in (0.01, 0.1):
                                                        for k9 in (7.5, 15):
                                                            k7 = Tm*m0/p0 # (5/50 = 0.1) Proportion of spliced and exported transcripts per time step.
                                                            # k9 = Tp*p0 + Tm*m0 # (2.5 + 5 = 7.5) Number of transcripts transcribed per time step.
                                                            ## Measurements and criteria for each parameter set.
                                                            time_to_kd = tn # initial value assuming no detectable knockdown
                                                            kd_time = 0 # duration of knockdown before recovery
                                                            kd_m = 0.1*m0 # threshold for detectable knockdown
                                                            max_sec = 0 # for calculating maximum 2ºsiRNA accumulation
                                                            min_m = m0
                                                            ## Create arrays with zeros for storing the values of all intermediates throughout the simulation.
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

                                                               # force negative values to equal 0.
                                                               if (ds_arr[i] < 0):
                                                                   ds_arr[i] = 0
                                                               if (pri_arr[i] < 0):
                                                                   pri_arr[i] = 0
                                                               if (ug_arr[i] < 0):
                                                                   ug_arr[i] = 0
                                                               if (sec_arr[i] < 0):
                                                                   sec_arr[i] = 0
                                                               if (p_arr[i] < 0):
                                                                   p_arr[i] = 0
                                                               if (m_arr[i] < 0):
                                                                   m_arr[i] = 0

                                                               ## update measurements for each parameter set.
                                                               if (m_arr[i] < kd_m):
                                                                   kd_time = kd_time + Dt
                                                               if ((m_arr[i] < kd_m) and (time_to_kd == tn)):
                                                                   time_to_kd = t_arr[i]
                                                               if (max_sec < sec_arr[i]):
                                                                   max_sec = sec_arr[i]
                                                               if (min_m > m_arr[i]):
                                                                   min_m = m_arr[i]
                                                            param_val = [lm, lds, Tm, Tp, Tug, Tpri, Tsec, k1, k2, k3, k4, k5, k6, k8, k9, kd_time, time_to_kd, max_sec, min_m]
                                                            writer.writerow(param_val)
f.close()
