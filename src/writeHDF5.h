#ifndef _WRITE_HDF5_
#define _WRITE_HDF5_

#include <hdf5.h>
#include <hdf5_hl.h>

#include "computelib.h"

void gather_JPLC(grid_parms* grid, double *jplc_buffer, EC_cell** ec);

void write_HDF5_JPLC(grid_parms* grid, double *jplc_buffer, char *path);

#endif

