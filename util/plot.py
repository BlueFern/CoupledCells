import numpy as np
import matplotlib.pyplot as plt


columns = [
#### SMC
# Fluxes
'smc_J_IP3',
'smc_J_SERCA',
'smc_J_CICR',
'smc_J_Extrusion',
'smc_J_Leak',
'smc_J_IP3_deg',

'smc_J_VOCC',
'smc_J_Na_Ca',
'smc_Na_K',
'smc_Cl',
'smc__K',
'smc_activation', # 11

# Homogeneous coupling_species
'smc_smc_cpl_Ca',
'smc_smc_cpl_Vm',
'smc_smc_cpl_IP3',

# Heterogeneous coupling_species
'smc_ec_cpl_Ca',
'smc_ec_cpl_Vm',
'smc_ec_cpl_IP3',

# Equations
'smc_dcdt', 
'smc_dsdt',
'smc_dvdt',
'smc_dwdt',
'smc_didt',

# Concentrations
'smc_Ca',
'smc_SR',
'smc_Vm',
'smc_w',
'smc_IP3',


### EC
# Fluxes
'ec_J_IP3', # 28
'ec_J_SERCA',
'ec_J_CICR',
'ec_J_Extrusion',
'ec_J_Leak',
'ec_J_IP3_deg',

'ec_J_NSC',
'ec_J_BK_Ca',
'ec_J_SK_Ca',
'ec_J_Ktot',
'ec_J_Residual',
'ec_J_trivial_Ca',

# Homogeneous coupling_species
'ec_ec_cpl_Ca', # 40
'ec_ec_cpl_Vm',
'ec_ec_cpl_IP3',

# Heterogeneous coupling_species
'ec_smc_cpl_Ca',
'ec_smc_cpl_Vm',
'ec_smc_cpl_IP3',

# Equations
'ec_dcdt',  # 46
'ec_dsdt',
'ec_dvdt',
'ec_didt', 

# Concentrations
'ec_Ca',
'ec_SR',
'ec_Vm',
'ec_IP3',

# Time
't'
]

# For saving plots
PATH = 'plots/rk_suite/'

data = np.genfromtxt('t55.csv', delimiter=',', skip_header=0,
                     skip_footer=0, names=columns)
              

# 3 x 3 subplots of equations
fig1 = plt.figure() 
val = 0
for i in range(18, 23):
    tmp = fig1.add_subplot(3,3,val)
    val += 1
    tmp.set_xlabel('Time (s)')
    tmp.set_ylabel(columns[i] + ' value')
    tmp.set_title(columns[i], y=0.7)
    tmp.plot(data[columns[-1]], data[columns[i]], c=np.random.rand(3,))
    
for i in range(46, 50):
    tmp = fig1.add_subplot(3,3,val)
    val += 1
    
    tmp.set_xlabel('Time (s)')
    tmp.set_ylabel(columns[i] + ' value')
    tmp.set_title(columns[i], y=0.7)
    tmp.plot(data[columns[-1]], data[columns[i]], c=np.random.rand(3,))

# 3 x 3 subplots of concentrations
fig2 = plt.figure() 
val = 0
for i in range(23, 28):
    tmp = fig2.add_subplot(3,3,val)
    val += 1
    tmp.set_xlabel('Time (s)')
    tmp.set_ylabel(columns[i] + ' value')
    tmp.set_title(columns[i], y=0.7)
    tmp.plot(data[columns[-1]], data[columns[i]], c=np.random.rand(3,))
    
for i in range(50, 54):
    tmp = fig2.add_subplot(3,3,val)
    val += 1
    
    tmp.set_xlabel('Time (s)')
    tmp.set_ylabel(columns[i] + ' value')
    tmp.set_title(columns[i], y=0.7)
    tmp.plot(data[columns[-1]], data[columns[i]], c=np.random.rand(3,))

plt.ioff()

# All plots into sub directories
for i in range(0, 54):

    fig = plt.figure()
    tmp = fig.add_subplot(111)
    tmp.set_xlabel('Time (s)')
    tmp.set_ylabel(columns[i] + ' value')
    tmp.set_title(columns[i], y=0.7)
    tmp.plot(data[columns[-1]], data[columns[i]], c=np.random.rand(3,))
    fig.savefig(PATH + 'pdf/' + columns[i] + '_plot.pdf', bbox_inches='tight')
    fig.savefig(PATH + 'svg/' + columns[i] + '_plot.svg', bbox_inches='tight')
