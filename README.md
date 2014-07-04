CoupledCells
============

This project uses massively parallel simulations to investigate the relationship
between arterial geometry, the transport of information via ionic
species through the gap junctions along arterial wall and the onset of
atherosclerotic plaques.

How to Compile
--------------

Run `make` in the `src` directory. The executable file is called `model`.

How to Run
----------

    model -args "-f <configFile> -S <solutionDirectory> -T <profilingDirectory> -t <duration> -w <checkpointFrequency> -i <delta>"

 where the command-line have the following meaning:

* `f` - Configuration file. For each subdomain the file must contain a line in the following format: 
    * Element 0: `Key_val` or serial number of the subdomain.
    * Element 1: Subdomain Type. There are two possibilities with their corresponding values:
        + Straight Segment (STRSEG) (0)
		+ Bifurcation (BIF) (1)
    * Element 2: Number of quads in the axial extent of current `Key_val` subdomain.
    * Element 3: Number of quads in the circumferential extent of the current `Key_val` subdomain.
    * Element 4: Parent subdomain `Key_val` of current `Key_val`.
    * Element 5: Left child subdomain `Key_val` of the current `Key_val`.
    * Element 6: Right child subdomain `Key_val` of the current `Key_val`.
    * Element 7: Required number of ECs axially per core.
    * Element 8: Required number of SMCs circumferentially per core.
 
For example for a fist subdomain of type BIF, with 12 cores along the axial
direction, 112 cores in the cirfumferrential direction, with no parent or child
subdomains, and with 32 ECs and 3 SMCs in the axial and cirfumferrential
directions respectively, the single line in the config file will look like
this:

    0,1,12,112,-1,-1,-1,32,3

* `S` - Location of the generated output files.
* `T` - Location of the pofiling output files.
* `w` - Duration of the simulation in seconds.
* `i` - Delta in milliseconds.

    model -f config.txt -S solution -T profiling -t 500.00 -w 1.0 -i 1e-2

Input Files
-----------

The input files which contain the subdomain geometry mapped to the UV quads/MPI cores shoud be in
the `files` subdirectory of the current working directory.
    

*WORK IN PROGRESS*




