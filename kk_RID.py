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

#Absorption spectra simulation 
voltages = np.linspace(0, 10, 10)  
frequencies = np.linspace(3e15, 9e15, 1000)  
absorption_spectra = []
for voltage in voltages:
    peak1_center = 4.5e15 - voltage * 1e14  
    peak2_center = 7.0e15 + voltage * 1e14  
    peak1 = np.exp(-((frequencies - peak1_center) / 5e14)**2) * (1 - voltage / 10)
    peak2 = np.exp(-((frequencies - peak2_center) / 5e14)**2) * (voltage / 10)
    absorption_spectra.append(peak1 + peak2)

# Calculate refractive index changes for each voltage
ref_index_initial = 1.4  
refractive_index_changes = [
    kramers_kronig_relation(frequencies, spectrum, ref_index_initial)
    for spectrum in absorption_spectra
]

# Plot results
plt.figure(figsize=(12, 8))

# Plot absorption spectra
plt.subplot(2, 1, 1)
for i, voltage in enumerate(voltages):
    plt.plot(frequencies, absorption_spectra[i], label=f"{voltage:.1f} V")
plt.title("Absorption Spectra with Shifting Peaks at Different Voltages")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Absorption")
plt.legend()

# Plot refractive index changes
plt.subplot(2, 1, 2)
for i, voltage in enumerate(voltages):
    plt.plot(frequencies, refractive_index_changes[i], label=f"{voltage:.1f} V")
plt.title("Refractive Index Change at Different Voltages")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Refractive Index")
plt.legend()

plt.tight_layout()
plt.show()