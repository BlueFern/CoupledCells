#ifndef __GATHER_H
#define __GATHER_H

#include <mpi.h>

#include "computelib.h"

/* The ec_data_buffer provides storage required for collecting all EC state variables
 * in one branch as required for writing the state variables to disk. */
typedef struct
{
	int chunk_size;
	int num_chunks;
	double *_ec_Ca;
	double *_ec_cpl_Ca;
	double *_ec_IP3;
	double *_ec_cpl_IP3;
	double *_ec_Vm;
	double *_ec_cpl_Vm;
	double *_ec_SR;

} ec_data_buffer;

/* The smc_data_buffer provides storage required for collecting all SMC state variables
 * in one branch as required for writing the state variables to disk. */
typedef struct
{
	int chunk_size;
	int num_chunks;
	double *_smc_Ca;
	double *_smc_cpl_Ca;
	double *_smc_IP3;
	double *_smc_cpl_IP3;
	double *_smc_Vm;
	double *_smc_cpl_Vm;
	double *_smc_SR;
	double *_smc_W;

} smc_data_buffer;

/* Allocate memory for the ec_data_buffer. */
ec_data_buffer *allocate_EC_data_buffer(int tasks_per_branch, int elements_per_task, int deep);

/* Allocate memory for the smc_data_buffer. */
smc_data_buffer *allocate_SMC_data_buffer(int tasks_per_branch, int elements_per_task, int deep);

/* Free the memory in the ec_data_buffer. */
void free_EC_data_buffer(ec_data_buffer *smc_buffer, int deep);

/* Free the memory in the smc_data_buffer. */
void free_SMC_data_buffer(smc_data_buffer *smc_buffer, int deep);

/* Collect state variables from the EC cells in one branch. */
void gather_EC_data(grid_parms *grid, ec_data_buffer *ec_buffer, EC_cell **ec);

/* Collect state variables from the SMC cells in one branch. */
void gather_SMC_data(grid_parms *grid, smc_data_buffer *smc_buffer, SMC_cell **smc);

/* Collect JPLC values from the EC cells in one branch. */
void gather_JPLC(grid_parms* grid, double *jplc_buffer, EC_cell** ec);

/* Collect WSS values from the EC cells in one branch. */
void gather_WSS(grid_parms* grid, double *jplc_buffer, EC_cell** ec);

#endif
