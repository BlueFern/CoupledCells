import numpy as np
import matplotlib.pyplot as plt


#TODO fix orderings... 54 things now D: + time
columns = [
'smc_J_IP3',
'smc_J_SERCA',
'smc_J_CICR',
'smc_J_Extrusion',
'smc_J_Leak',
'smc_J_VOCC',
'smc_J_Na_Ca',

'smc_gamma_J_Na_K',
'smc_gamma_J_Cl',
'smc_gamma_J_VOCC',
'smc_gamma_J_Na_Ca',
'smc_gamma_J_K',

'smc_J_IP3_deg',

'smc_lambda_K_activation',
'smc_lambda_smc_w',

'smc_c',
'smc_s',
'smc_v',
'smc_i',
'smc_w',

'ec_J_IP3',
'ec_J_SERCA',
'ec_J_CICR',
'ec_J_Extrusion',
'ec_J_Leak',
'ec_J_NSC',
'ec_J_trivial_Ca',

'ec_J_Ktot',
'ec_J_Residual',
'ec_JPLC',
'ec_J_IP3_deg',

'ec_c',
'ec_s',
'ec_v',
'ec_i', 
't'
]

data = np.genfromtxt('everything_no_atp_300s.csv', delimiter=',', skip_header=0,
                     skip_footer=0, names=columns)
              

# TODO: write out files of sections (eg smc_fluxes.txt, smc_equations.txt etc)
fig = plt.figure()         
val = 0
for i in range(15, 20):
    tmp = fig.add_subplot(3,3,val)
    val += 1
    tmp.set_xlabel('Time (s)')
    tmp.set_ylabel(columns[i] + ' value')
    tmp.set_title(columns[i], y=0.7)
    tmp.plot(data[columns[-1]], data[columns[i]], c=np.random.rand(3,))
    
for i in range(31, 35):
    tmp = fig.add_subplot(3,3,val)
    val += 1
    
    tmp.set_xlabel('Time (s)')
    tmp.set_ylabel(columns[i] + ' value')
    tmp.set_title(columns[i], y=0.7)
    tmp.plot(data[columns[-1]], data[columns[i]], c=np.random.rand(3,))
'''

for i in range(0, 35):
    if i in range(15,20) or i in range(31, 35):
        continue
    fig = plt.figure()
    tmp = fig.add_subplot(111)
    tmp.set_xlabel('Time (s)')
    tmp.set_ylabel(columns[i] + ' value')
    tmp.set_title(columns[i], y=0.7)
    tmp.plot(data[columns[-1]], data[columns[i]], c=np.random.rand(3,))
    fig.savefig('plots/' + columns[i] + '_plot.png', bbox_inches='tight')
    
    
    '''
#    max_x = max(data[filenames[i]])
#    min_x = min(data[filenames[i]])
#    plt.plot([(x - min_x) / (max_x - min_x) for x in data[filenames[i]]])
    

#plt.legend(legend)
#plt.show()

