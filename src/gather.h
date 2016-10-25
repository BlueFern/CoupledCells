#ifndef __GATHER_H
#define __GATHER_H

#include <mpi.h>

#include "computelib.h"

/* Collect state variables from the EC cells in one branch. */
void gather_EC_data(grid_parms *grid, double *ec_buffer, EC_cell **ec);

/* Collect state variables from the SMC cells in one branch. */
void gather_SMC_data(grid_parms *grid, double *smc_buffer, SMC_cell **smc);

/* Collect JPLC values from the EC cells in one branch. */
void gather_JPLC(grid_parms* grid, double *jplc_buffer, EC_cell** ec, int);

#endif
