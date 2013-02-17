#!/bin/sh
# @ output = $(job_name).$(Cluster).out
# @ error= $(job_name).$(Cluster).err
# @ job_type = bluegene
# @ wall_clock_limit = 1:00:00,1:00:00
# @ job_name =bifurcation_test
# @ class=bg512_100
# @ bg_size=180
# @ notification = error
# @ group=bg01
# @ queue
mpirun -mode CO  -np 180 -cwd `pwd` -connect MESH -exe `pwd`/test -verbose 2
