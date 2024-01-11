import matplotlib.pyplot as plt
import numpy as np
import math

#INTRO
#the sodium channel has 2 gating particles, m and h. The m particle is called the activation particle and opens when the cell is depolarized. 
#The h particle is called the inactivation particle, and “closes” (it’s value approaches 0) when the cell is depolarized. 
#However, when the cell is depolarized, the m particle opens faster than the h particle closes, giving you the transient nature of the sodium current. 


#Functions to calculate the alpha and beta rates (voltage converted to mV)

def alpha_n(v):
    v = v*1000
    return 0.01 * (-v-55)/( math.exp((-v-55)/10.0) -1) * 1000
 
def beta_n(v):
    v = v*1000
    return 0.125*math.exp((-v-65)/80.0) * 1000
 
def alpha_m(v):
    v = v*1000
    return 0.1 * (-v-40)/( math.exp((-v-40)/10.0) -1) * 1000
 
def beta_m(v):
    v = v*1000
    return 4*math.exp((-v-65)/18.0) * 1000
 
def alpha_h(v):
    v = v*1000
    return 0.07*math.exp((-v-65)/20.0) * 1000
 
def beta_h(v):
    v = v*1000
    return 1/( math.exp((-v-35)/10.0) +1) * 1000


dt = 10E-6              #10 us timestep
Cm = 100E-12            #membrane capacitance 100pF
v_init = -80E-3         #initial membrane potential


#HH Potassium channel
n = alpha_n(v_init)/(alpha_n(v_init)+beta_n(v_init))
m = alpha_m(v_init)/(alpha_m(v_init)+beta_m(v_init))
h = alpha_h(v_init)/(alpha_h(v_init)+beta_h(v_init))

Gbar_k = 1E-6           #potassium conductance
Gbar_na = 7E-6          #sodium conductance
G_leak = 5E-9           #leak

E_k = -80E-3            #reversal for HH potassium current
E_na = 40E-3            #reversal for HH potassium current
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
        n += dn

        dm = (alpha_m(v_out[i-1]) * (1 - m) - beta_m(v_out[i-1]) * m) * dt      #particles
        m += dm

        dn = (alpha_h(v_out[i-1]) * (1 - h) - beta_h(v_out[i-1]) * h) * dt      #particles
        h += dn
    
    G_k = Gbar_k*n**4                       #potassium channel current
    i_k = G_k*(v_out[i-1] - E_k)            #HH potassium channel current

    G_na = Gbar_na*m**3*h                   #sodium channel current
    i_na = G_na*(v_out[i-1] - E_na)         #HH sodium channel current

    i_leak = G_leak * (v_out[i-1] - E_leak) #current through ion channels
    i_cap = i_inj[i] - i_leak - i_k - i_na  #real current
    dv = (i_cap/Cm)*dt                      #calculation of dv
    v_out[i]=v_out[i-1] + dv                #last known voltage increased by vd
    

#Make the graph
t_vec = np.linspace(0, dt*np.size(v_out), np.size(v_out))
plt.plot(t_vec, v_out)
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.show()