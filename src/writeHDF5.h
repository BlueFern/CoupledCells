#ifndef _WRITE_HDF5_
#define _WRITE_HDF5_

#include <hdf5.h>
// #include <hdf5_hl.h>

#include "computelib.h"
#include "gather.h"

#define FAIL -1

/* Use %ld to print the value because long should cover most cases. */
#define CHECK(ret, val, where) do { \
    if ((ret) == (val)) { \
	printf("*** UNEXPECTED RETURN from %s is %ld at line %4d in %s\n", where, (long)(ret), (int)__LINE__, __FILE__); \
    } \
} while(0)

void write_dataset_HDF5(grid_parms*, double *, char *path, char*);

void write_EC_data_HDF5(grid_parms*, ec_data_buffer*, int, char*);

void write_SMC_data_HDF5(grid_parms*, smc_data_buffer*, int, char*);

#endif

