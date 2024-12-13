import numpy as np
import matplotlib.pyplot as plt

def kramers_kronig_relation(frequencies, absorption_spectrum, ref_index_initial=1):
    n_change = np.zeros_like(frequencies)
    delta_omega = frequencies[1] - frequencies[0]  # Frequency step size

    for i, omega in enumerate(frequencies):
        integral = 0.0
        for j, omega_prime in enumerate(frequencies):
            if j != i:  # Skip the singularity
                integral += (omega_prime * absorption_spectrum[j]) / (omega_prime**2 - omega**2) * delta_omega
        n_change[i] = (2 / np.pi) * integral

    return n_change + ref_index_initial
