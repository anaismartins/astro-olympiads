from astropy import units as u
from astropy import constants as const
import numpy as np
import matplotlib.pyplot as plt

planck_bands = [30, 44, 70, 100, 143, 217, 353, 545, 857]


h = 6.62607015e-34 * 1e9  # J GHz-1
c = 299792458e-9  # m GHz
k = 1.380649e-23  # J K-1

def planck(freq, temp):
    """
    Planck function returning in units of MJy/sr.
    Input frequency in GHz and temperature in K.
    """

    b = 2 * h * freq**3 / c**2 / (np.exp(h * freq / (k * temp)) - 1) * 1e20  # MJy sr-1

    return b

nu = np.linspace(10, 1000, 3000)
Tcmb = 2.73
monopole = planck(nu, Tcmb)

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple']

plt.plot(nu, monopole)
for band in planck_bands:
    plt.vlines(band, min(monopole), max(monopole), label=f'{band} GHz', linestyle='--', color=colors[planck_bands.index(band)])
plt.xlabel('Frequency (GHz)')
plt.ylabel('Intensity (MJy/sr)')
plt.xscale('log')
plt.yscale('log')
plt.title('CMB Monopole Spectrum')
plt.legend()
plt.tight_layout()
plt.savefig("1_monopole_spectrum.png")
plt.clf()

Tanisotropies = 0.0000035

anisotropies = monopole * np.exp(h * nu / (k * Tcmb)) / (np.exp(h * nu / (k * Tcmb)) - 1) * k * nu / (k * Tcmb * Tcmb) * Tanisotropies

plt.plot(nu, anisotropies)
for band in planck_bands:
    plt.vlines(band, min(anisotropies), max(anisotropies), label=f'{band} GHz', linestyle='--', color=colors[planck_bands.index(band)])
plt.xlabel('Frequency (GHz)')
plt.ylabel('Intensity (MJy/sr)')
plt.xscale('log')
plt.yscale('log')
plt.title('CMB Anisotropy Spectrum')
plt.legend()
plt.tight_layout()
plt.savefig("2_anisotropy_spectrum.png")
plt.clf()

A_synch = 30
nu0_synch = 0.408 # GHz
beta_synch = -3.07
synchrotron = A_synch * (nu / nu0_synch)**(beta_synch)

tau_dust = 9.6e-7
Tdust = 19.7 # K
nu0_dust = 353 # GHz
beta_dust = 1.62
dust = tau_dust * planck(nu, Tdust) * (nu / nu0_dust)**(beta_dust)

plt.plot(nu, anisotropies, label='CMB Anisotropies')
plt.plot(nu, synchrotron, label='Synchrotron')
plt.plot(nu, dust, label='Dust')
plt.plot(nu, synchrotron + dust, label='Sum of foregrounds')
for band in planck_bands:
    plt.vlines(band, min(synchrotron), max(dust), linestyle='--', color=colors[planck_bands.index(band)], label=f'{band} GHz')
plt.xlabel('Frequency (GHz)')
plt.ylabel('Intensity (MJy/sr)')
plt.xscale('log')
plt.yscale('log')
plt.title('CMB Anisotropy and Foregrounds Spectrum')
plt.legend()
plt.tight_layout()
plt.savefig("3_foregrounds_spectrum.png")
plt.clf()