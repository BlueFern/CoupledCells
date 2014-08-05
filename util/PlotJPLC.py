# This will plot the Sigmoid function controlling JPLC contentration in the
# simulated vessel bifurcation according to the JPLC min, max, gradient, number
# of ECs per quad, number of quads, length of a branch in this bifurcation,
# and angle of the branches. These parameters are to be written in the
# JPLC_params.py file.

import numpy as np
import matplotlib.pyplot as plt
import math

# Read parameters from file.
import JPLC_Params

EC_length = JPLC_Params.length / (JPLC_Params.num_quads * JPLC_Params.num_ECs)

trunk_points = np.arange(-JPLC_Params.length, 0, EC_length)
branch_points = np.arange(0, (math.sin(JPLC_Params.angle) * JPLC_Params.length), (math.sin(JPLC_Params.angle) * EC_length))

#print trunk_points.size
#print branch_points.size

all_points = np.concatenate((trunk_points,branch_points))

def JPLC(x):
    return JPLC_Params.min_jplc + (JPLC_Params.max_jplc / (1.0 + np.exp(-JPLC_Params.gradient * x)))

plt.plot(all_points, JPLC(all_points), 'b')
plt.show()
