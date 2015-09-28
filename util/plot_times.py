# -*- coding: utf-8 -*-

"Provide a glob expression selecting a group of files to plot the values."

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt

arg = sys.argv[1]

base_name = os.path.splitext(os.path.basename(arg))[0].replace('*', 'ALL')

matching_files = glob.glob(arg)
print 'Matching files:', matching_files

plt.title(base_name)
plt.ylabel("seconds")
plt.xlabel("core #")

for file in matching_files:
    data = np.loadtxt(file)
    plt.plot(data, 'o')
    
plt.savefig(base_name + ".pdf", bbox_inches='tight')
plt.show()