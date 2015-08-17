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
	double *_ec_SR;
	double *_ec_Vm;
	double *_ec_cpl_Vm;


} ec_data_buffer;

/* Allocate memory for the ec_data_buffer. */
ec_data_buffer *allocate_EC_data_buffer(int tasks_per_branch, int elements_per_task, int deep);

/* Free the memory in the ec_data_buffer. */
void free_EC_data_buffer(ec_data_buffer *ec_buffer);

/* Collect state variables from the EC cells in one branch. */
void gather_EC_data(grid_parms *grid, ec_data_buffer *ec_buffer, EC_cell **ec);

#endif
