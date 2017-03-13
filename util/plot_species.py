import numpy as np
import matplotlib.pyplot as plt
import argparse
from pylab import rcParams
from itertools import cycle

rcParams['figure.figsize'] = 9, 5
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18


import colorsys
N = 13
HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
color_cycle = cycle(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))


columns = [
    'smc_vm', 
    'smc_ca',
    'smc_sr',
    'smc_w',
    'smc_i',
    
    'ec_vm',
    'ec_ca',
    'ec_sr',
    'ec_i',
    
    't'
]
labels = [
    'SMC membrane potential (mV)',
    'SMC cytosolic $Ca^{2+}$ $(\mu M)$',
    'SMC SR $Ca^{2+}$ $(\mu M)$',
    'Open state probability of \nthe calcium-activated potassium channels ',
    'SMC $IP_3$ $(\mu M)$',

    'EC membrane potential (mV)',
    'EC cytosolic $Ca^{2+}$ $(\mu M)$',
    'EC SR $Ca^{2+}$ $(\mu M)$',
    'EC $IP_3$ ($\mu M)$',
]


def plot_singles(args):
    
        data = np.genfromtxt(args.data, delimiter=',', skip_header=0,
                         skip_footer=0, names=columns)
                         
        for i in range(0, len(columns) - 1):
   
            fig = plt.figure() 

            plt.xlabel('Time (s)',fontsize=20)
            plt.ylabel(labels[i],fontsize=20,labelpad=18)
            plt.plot(data[columns[-1]][1:], data[columns[i]][1:], c=color_cycle.next(),linewidth=2.0)

            plt.tight_layout()

            #fig.savefig('path'+columns[i]+'.png', dpi=150)
                         
                         
    
def plot(args):

    data = np.genfromtxt(args.data, delimiter=',', skip_header=0,
                         skip_footer=0, names=columns)
    
    # 3 x 3 subplots of equations
    fig1 = plt.figure() 
    val = 1
    
    # If the colour is not provided, use a defualt.
    if args.colour:
        colour = args.colour
    else:
        print "No colour has been set (or detected). Using default."
        colour = (0.7,0.2,0.1)
    
    for i in range(0, len(columns) - 1):
        
        tmp = fig1.add_subplot(3, 3, val)
        val += 1
        tmp.set_xlabel('Time (s)')
        tmp.set_ylabel(columns[i] + ' value')
        tmp.set_title(columns[i], y=0.5, x=0.1)
        tmp.plot(data[columns[-1]], data[columns[i]], c=colour)

    plt.show()


if __name__ == "__main__":

    argParser = argparse.ArgumentParser()
    argParser.add_argument('data', type=str, help='Input data name/path.')
    argParser.add_argument('-c', '--colour',type=str, help='Hexidecimal STRING specifying the colour of the output.')
    args = argParser.parse_args()
    
    plot(args)