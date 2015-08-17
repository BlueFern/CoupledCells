#ifndef _WRITE_HDF5_
#define _WRITE_HDF5_

#include <hdf5.h>
#include <hdf5_hl.h>

#include "computelib.h"
#include "gather.h"

void gather_JPLC(grid_parms* grid, double *jplc_buffer, EC_cell** ec);

void write_HDF5_JPLC(grid_parms* grid, double *jplc_buffer, char *path);

void write_EC_data_HDF5(grid_parms* grid, ec_data_buffer *ec_buffer, int write_count, char* path);

#endif

