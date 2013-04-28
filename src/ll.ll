# MPI LoadLeveler Job file for BGP
# @ shell = /bin/ksh
#
# @ job_name = tree
#
# @ job_type = bluegene
#
# @ wall_clock_limit     = 10:00:00
#
# @ group = NZ_merit
# @ account_no = nesi00052
# @ output               = $(job_name).$(jobid).out
# @ error                = $(job_name).$(jobid).err
# @ notification         = never
# @ class                = bgp
# @ bg_connection = prefer_torus
# @ bg_size = 256
# @ queue
/bgsys/drivers/ppcfloor/bin/mpirun -mode VN -np 1024 -cwd `pwd` -exe `pwd`/model -verbose 2

