# Let's start by importing everything we need.
import os
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import matplotlib.pyplot as plt
from petitRADTRANS import Radtrans
from petitRADTRANS import nat_cst as nc
from petitRADTRANS.retrieval.parameter import Parameter
from petitRADTRANS.retrieval import plot_style as ps

from petitRADTRANS.retrieval.models import emission_model_diseq,\
                                           emission_model_diseq_patchy_clouds,\
                                           emission_model_diseq_simple_patchy_clouds,\
                                           guillot_emission,\
                                           guillot_patchy_emission,\
                                           interpolated_profile_emission,\
                                           gradient_profile_emission



parameters = {}
parameters['D_pl'] = Parameter('D_pl', False, value = 26.84*nc.pc)
parameters['mass'] = Parameter('mass',False,value = 3.5 * nc.m_jup)
parameters['R_pl'] = Parameter('R_pl', False, value = 1.3 * nc.r_jup_mean)

parameters['T_int'] = Parameter('T_int', False, value = 1500.0)
parameters['T1'] = Parameter('T1', False, value = 0.5)
parameters['T2'] = Parameter('T2', False, value = 0.4)
parameters['T3'] = Parameter('T3', False, value = 0.8)
parameters['log_delta'] = Parameter('log_delta', False, value = 0.65)
parameters['alpha'] = Parameter('alpha', False, value = 1.70)

parameters['Fe/H'] = Parameter('[Fe/H]', False, value = 1.0)
parameters['C/O'] = Parameter('C/O', False, value = 0.7)
parameters['log_pquench'] = Parameter('log_pquench', False, value = 2.5)

parameters['sigma_lnorm'] = Parameter('sigma_lnorm', False, value = 2.0)
parameters['fsed_MgSiO3(c)'] = Parameter('fsed_MgSiO3(c)', False, value = 3.)
parameters['fsed_Fe(c)'] = Parameter('fsed_Fe(c)', False, value = 8.)
parameters['log_kzz'] = Parameter('log_kzz', False, value = 10.)

parameters['eq_scaling_Fe(c)'] = Parameter('Fe(c)', False, value = -1.00)
parameters['eq_scaling_MgSiO3(c)'] = Parameter('MgSiO3(c)', False, value = -1.00)


resolution = 300

# line_species = [f'H2O_Exomol_R_{resolution}',
#                 f'CO_all_iso_HITEMP_R_{resolution}',
#                 f'CH4_R_{resolution}',
#                 f'CO2_R_{resolution}',
#                 f'HCN_R_{resolution}',
#                 f'FeH_R_{resolution}',
#                 f'H2S_R_{resolution}',
#                 f'NH3_R_{resolution}',
#                 f'PH3_R_{resolution}',
#                 f'Na_allard_R_{resolution}',
#                 f'K_allard_R_{resolution}',
#                 f'TiO_all_Exomol_R_{resolution}',
#                 f'VO_R_{resolution}',
#                 f'SiO_R_{resolution}']

line_species = [f'H2O_Exomol',
                f'CO_all_iso_HITEMP',
                f'CH4',
                f'CO2',
                f'HCN',
                f'FeH',
                f'H2S',
                f'NH3',
                f'PH3',
                f'Na_allard',
                f'K_allard',
                f'TiO_all_Exomol',
                f'VO',
                f'SiO']

rayleigh_species = ['H2', 'He']
continuum_opacities = ['H2-H2', 'H2-He']
cloud_species = ['MgSiO3(c)_cd', 'Fe(c)_cd']

atmosphere = Radtrans(line_species = line_species,
                      rayleigh_species= rayleigh_species,
                      continuum_opacities = continuum_opacities,
                      cloud_species = cloud_species,
                      mode='c-k',
                      do_scat_emis=True,
                      wlen_bords_micron = [0.8,5.2])
pressures = np.logspace(-6,2,100)
atmosphere.setup_opa_structure(pressures)

wavelength, model = emission_model_diseq(atmosphere, parameters, AMR = False, PT_plot_mode = False)
pressure,temperature = emission_model_diseq(atmosphere, parameters, AMR = False, PT_plot_mode = True)

fig, ax = plt.subplots(figsize = (8,5))
ax.plot(wavelength, model, linewidth = 3)
ax.set_xlabel("Wavelength [micron]")
ax.set_ylabel(r"F$_{\lambda}$ [W/m$^{2}/\mu$m]")
plt.show()

fig, ax = plt.subplots(figsize = (8,5))
ax.plot(temperature, pressure, linewidth = 3)
ax.set_xlabel("Temperature [K]")
ax.set_ylabel(r"Pressure [bar]")
ax.set_ylim(1e2,1e-6)
ax.set_yscale('log')
plt.show()
