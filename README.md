CoupledCells
============

This project uses massively parallel simulations to investigate the relationship
between arterial geometry, the transport of information via ionic
species through the gap junctions along arterial wall and the onset of
atherosclerotic plaques.

Code Dependencies
-----------------

The project depends on MPI and SUNDIALS libraries.

In addition, the Python scripts for converting the simulations output from HDF5 format
to VTK (VTU) format require VTK with Python bindings.


How to Compile
--------------

The project files can be generated with CMake. The project has been previously compiled and tested
on Linux machines with CMake and Eclipse. Configure the project with CMake and specify the *out-of-source* build
directory. CMake can be run in GUI or CLI modes.

The desired ODE solver should be chosen in the CMake interface. There are three options available:

ODEOPTION=RK_SUITE, ODEOPTION=ARK_ODE, or ODEOPTION=BOOST_ODEINT.

It is recommended to use the ARK_ODE solver, which is a part of the SUNDIALS library.

After the the project has been configured, it can be opened and compiled with
the target IDE, or in the case of Unix Makefile configuration simply run make in
the build directory. The generated executable is called *coupledCellsModel*.

There are two makefiles for compiling the project on BlueGene/L and BlueGene/P. For the BlueGene/L build, the
Makefile.bgp can be used as follows:

    make -f Makefile.bgp ODEOPTION=ARK_ODE

In the current BlueGene/P build environment CMake should be launched in the following manner in order to 
specify the MPI compiler and the location of the Sundials library:

    CXX=mpixlcxx CC=mpixlc ccmake ../CoupledCells/. -DSUNDIALS_DIR=/bgp/local/pkg/sundials/2.5.0

How to Run
----------

    coupledCellsModel -args "-f <configFile> -S <solutionDirectory> -T <profilingDirectory> -t <duration> -w <checkpointFrequency> -i <delta>"

 where the command-line have the following meaning:

* **f** - Configuration file. Each file must start with a line stating the total
  number of domains. Each line in the config file must be terminated with a ';'.
  For each subdomain the file must contain a line in the following format:
    * Element 0: *Key_val* or serial number of the subdomain.
    * Element 1: Subdomain Type. There are two possibilities with their corresponding values:
        + Straight Segment (STRSEG) (0)
		+ Bifurcation (BIF) (1)
    * Element 2: Number of quads in the axial extent of current *Key_val* subdomain.
    * Element 3: Number of quads in the circumferential extent of the current
      *Key_val* subdomain.
    * Element 4: Parent subdomain *Key_val* of current *Key_val*.
    * Element 5: Left child subdomain *Key_val* of the current *Key_val*.
    * Element 6: Right child subdomain *Key_val* of the current *Key_val*.
    * Element 7: Required number of ECs axially per core.
    * Element 8: Required number of SMCs circumferentially per core.

For example, for a fist subdomain of type BIF, with 12 cores along the axial
direction, 112 cores in the cirfumferrential direction, with no parent or child
subdomains, and with 32 ECs and 3 SMCs in the axial and cirfumferrential
directions respectively, the config file will look like
this:

    1;
    0,1,12,112,-1,-1,-1,32,3;

* **S** - Location of the generated output files.
* **T** - Location of the profiling output files.
* **t** - Duration of the simulation in seconds.
* **w** - How often to write checkpont data in seconds.
* **i** - Delta in milliseconds specifying the frequency of exchanging data betwen
  UV quads/MIP processes.

For example, to run a simulation for 500 seconds, with checkpoints written every
second, the arguments would look like this:

    coupledCellsModel -f config.txt -S solution -T profiling -t 500.00 -w 1000 -i 1e-2

The following in an example of a load-leveller script to run a 100 physiological second
simulation with a 1008 quad/core bifurcation mesh.

	#!/bin/ksh
	#
	# MPI LoadLeveler Job file for BG/P
	#
	# @ group                = UC
	# @ account_no           = bfcs00321
	# @ job_type             = bluegene
	# @ bg_connection        = prefer_torus
	# @ output               = $(job_name).$(jobid).out
	# @ error                = $(job_name).$(jobid).err

	# @ bg_size              = 256
	# @ class                = bgp
	# @ wall_clock_limit     = 04:00:00
	# @ job_name             = CC_1008
	
	# @ queue
	
	mpirun -mode VN -np 1008 -verbose 2 -env BG_COREDUMP_BINARY='*' -cwd `pwd` -exe `pwd`/coupledCellsModel -args "-f config.txt -S solution -T profiling -t 100.00 -w 1.0 -i 1e-2"



Input Files
-----------

The *files* subdirectory of the current working directory must contain the input
ATP files in TXT format: `parent_atp.txt`, `left_daughter_atp.txt`, and `right_daughter_atp.txt`
for a simulation involving a bifurcation. Only 'parent_atp.txt' is required for a simulation
involving only a trunk. Details on creating these files can be found in the
[DBiharMesher](https://github.com/BlueFern/DBiharMesher) repository.

Output Files and Conversion
---------------------------

Once the simulation completes the solution directory will contain a number of HDF5-format files for
both ECs and SMCs. To visualise the results these files should be converted to the VTK's VTU format.
The scripts to convert EC and SMC output for both trunk and bifurcation simulations are provided in
the util directory. These sripts write to the solution directory and expect start and stop timestep
parameters (between which timesteps to convert). An example of such a conversion for a bifurcation
simulation run for 100 physiological seconds for both cell types would be the two commands:

	python /path/to/util/bifurcation_smc_hdf5ToVTU.py 0 100
	python /path/to/util/bifurcation_ec_hdf5ToVTU.py 0 100

A `jplc_*_.h5` file is written out for each branch of a simulation, intended for verifying the
input ATP files where read in correctly. These output files also allo to verify the mesh
dimensions and the corresponding ATP file are in agreement. The output `jplc_*_.h5` are
converted to VTU format using either `bifurcation_jplc_hdf5ToVTU.py` or `trunk_jplc_hdf5ToVTU.py`
scripts, depending on the type of the input mesh.

For these conversions the subdirectory `vtk` of the current working directory must exist and 
contain `*.vtp' files for each branch of the mesh configuration used in the simulation. Again,
for details on creating these files see the [DBiharMesher](https://github.com/BlueFern/DBiharMesher)
repository.

The script `genCommands.py` is provided in the util directory for parallel conversions --- it simply prints 
out a number of calls to the conversion scripts, where these calls can be redirected to a file. 
It should be copied into the working directory and modified according to the geometry of the mesh, (and so whether `trunk_*_hdf5ToVTU.py` or `bifurcation_*_hdf5ToVTU.py` is called) and the path to these scripts. It expects the
number of tasks being used and the total number of timesteps to process as arguments. `genCommands.py` is 
currently being used in conjunction with a load-leveller script, an example of which follows:

	# Example multiple Serial LoadLeveler Jobs file on multiple nodes.
	
	# @ shell              = /bin/bash
	# @ job_name           = hdf52vtu
	
	# This type of job is parallel even though your executable is serial, batcher is a MPI parallel program.
	# @ job_type           = parallel
	# @ wall_clock_limit   = 00:30:00
	# @ class              = p7linux_dev
	# @ group              = UC
	# @ account_no         = bfcs00321
	
	# Use the unlimited options for blocking to run tasks accross nodes.
	# @ blocking           = unlimited
	
	# Total_tasks needs to be more than 1.
	# @ total_tasks        = 8
	
	# Affinity options to improve performance.
	# @ task_affinity      = core(1)
	# @ rset               = rset_mcm_affinity
	# @ mcm_affinity_options = mcm_mem_pref
	
	# @ output             = $(job_name).$(schedd_host).$(jobid).out
	# @ error              = $(job_name).$(schedd_host).$(jobid).err
	# @ environment = COPY_ALL
	# @ queue
	
	# All commands that follow will be run as part of the job.
	# Display name of host running serial job.
	hostname
	
	# Display current time.
	date
	
	# Filename with commands to parse to batcher.
	export INPUT_FILE="batcher_input_cmds"-$LOADL_STEP_ID
	echo $INPUT_FILE
	
	module purge
	module load vtk
	
	# Create the Input Command file for batcher if using the create_input_command.sh script.
	python ./genCommands-p7.py 8 100 >> $INPUT_FILE
	
	# Do not modify the following command, it will launch your serial jobs with batcher.
	poe batcher $INPUT_FILE


Project Documentation
---------------------

The documentation for the project is very much *work in progress*. Currently two types of
documentation can be generated:

 1. Doxygen API documentation (if Doxygen is installed on the system).
 1. Sphinx high-level project documentation (if Sphinx and Doxylink are
 installed on the system). See [Sphinx documentation pages](http://sphinx-doc.org)
 and [Doxylink documentation pages](https://pypi.python.org/pypi/sphinxcontrib-doxylink)
 for details.
