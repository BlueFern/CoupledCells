import numpy as np
import matplotlib.pyplot as plt
import argparse

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
    'ec_Gprot', #active g proteins
    'ec_P_P2Y', #ratio of bound to total p2y
    'ec_R_PIP2_H', #rate of pip2 hydrolysis
    'ec_J_IP3_deg',
    'ec_J_ind_I',
    't'
]
    
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
        
        tmp = fig1.add_subplot(5, 3, val)
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