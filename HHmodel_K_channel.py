import matplotlib.pyplot as plt
import numpy as np
import math

# INTRO
#n is probability that a “particle” in potassium channel is in the open position. 
#In order for the channel to conduct, you need all four of these particles to be in the open position.


#Functions to calculate the alpha and beta rates (voltage converted to mV)

def alpha_n(v):
    v *= 1000
    return 0.01 * (-v-55)/( math.exp((-v-55)/10.0) -1) * 1000

def beta_n(v):
    v *= 1000
    return 0.125*math.exp((-v-65)/80.0) * 1000

dt = 10E-6              #10 us timestep
Cm = 100E-12            #membrane capacitance 100pF
v_init = -70E-3         #initial membrane potential

#HH Potassium channel
n = alpha_n(v_init)/(alpha_n(v_init)+beta_n(v_init))

Gbar_k = 1E-6              #potassium conductance
G_leak = 5E-9           #leak
E_k = -80E-3            #reversal for HH potassium current
E_leak = -70E-3         #reversal potential of -70 mV

current = 200E-12       #200pA current

#Injected current 0.2 seconds of no current, 0.3 seconds of current,and 0.5 seconds of no current
i_inj = np.concatenate((np.zeros([round(0.2/dt), 1]), current*np.ones([round(0.3/dt), 1]), np.zeros([round(0.5/dt), 1])))

#Spikes calculation
v_out = np.zeros(np.size(i_inj))

for i in range(np.size(v_out)):
    if i == 1:
        v_out[i] = v_init
    else:
        dn = (alpha_n(v_out[i-1]) * (1 - n) - beta_n(v_out[i-1]) * n) * dt      #particles
        n = n + dn
    
    G_k = Gbar_k*n**4                       #potassium channel current
    i_k = G_k*(v_out[i-1] - E_k)            #HH potassium channel current
    i_leak = G_leak * (v_out[i-1] - E_leak)    #current through ion channels
    i_cap = i_inj[i]  - i_leak - i_k        #real current
    dv = (i_cap/Cm)*dt                      #calculation of dv
    v_out[i]=v_out[i-1] + dv                #last known voltage increased by vd
    

#Make the graph
t_vec = np.linspace(0, dt*np.size(v_out), np.size(v_out))
plt.plot(t_vec, v_out)
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.show()