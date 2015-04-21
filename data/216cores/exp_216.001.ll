#!/bin/ksh
#
# MPI LoadLeveler Job file for BG/L
#
# @ group                = bg01
# @ job_type             = bluegene
# @ bg_connection        = prefer_torus
# @ output               = $(job_name).$(jobid).out
# @ error                = $(job_name).$(jobid).err
# @ notification         = always
# @ notify_user          = you@wherever.com

# @ bg_size              = 128
# @ class                = bg128_400
# @ wall_clock_limit     = 00:30:00
# @ job_name             = 216.001

mpirun -mode VN -np 216 -verbose 2 -env BG_COREDUMP_BINARY='*' -cwd `pwd` -exe `pwd`/coupledCellsModel -args "-f config.txt -S solution -T profiling -t 20.00 -w 1.0 -i 1e-2"

# @ queue

