#include <malloc.h>

#include "computelib.h"
#include "koenigsberger_model.h"
#include "gather.h"

void gather_EC_data(grid_parms *grid, double *ec_buffer, EC_cell **ec)
{
	// printf("[%d] ++++++ Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	int chunk_size = grid->num_ec_axially * grid->num_ec_circumferentially * (grid->num_coupling_species_ec + grid->neq_ec);

	// Allocate displacements for MPI comms.
	int *disp = (int*)checked_malloc((int)(grid->num_ranks_write_group * sizeof(int)), SRC_LOC);

	// Allocate counts for MPI comms.
	int *recv_count = (int*)checked_malloc(grid->num_ranks_write_group * sizeof(int), SRC_LOC);

	// printf("[%d] Checkpoint %s:%d\n", grid->universal_rank, __FILE__, __LINE__);

	if(grid->rank_write_group == 0)
	{
		// Compute displacement values.
		for (int i = 0; i < grid->num_ranks_write_group; i++)
		{
			disp[i] = i * chunk_size;
			recv_count[i] = chunk_size;
		}
	}

	// Allocate send buffer for all attributes for each cell.
	double* send_buffer = (double *)checked_malloc(sizeof(double) * chunk_size, SRC_LOC);

	int p = 0;

	for(int i = 0; i < grid->num_ec_axially; i++)
	{
		for(int j = 0; j < grid->num_ec_circumferentially; j++)
		{
			send_buffer[p++] = ec[j + 1][i + 1].vars[ec_Ca];
			send_buffer[p++] = ec[j + 1][i + 1].homo_fluxes[cpl_Ca];
			send_buffer[p++] = ec[j + 1][i + 1].vars[ec_IP3];
			send_buffer[p++] = ec[j + 1][i + 1].homo_fluxes[cpl_IP3];
			send_buffer[p++] = ec[j + 1][i + 1].vars[ec_SR];
			send_buffer[p++] = ec[j + 1][i + 1].vars[ec_Vm];
			send_buffer[p++] = ec[j + 1][i + 1].homo_fluxes[cpl_Vm];
			send_buffer[p++] = ec[j + 1][i + 1].vars[ec_Gprot];

		}
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, chunk_size, MPI_DOUBLE, ec_buffer, recv_count, disp, MPI_DOUBLE, 0, grid->write_group));

	free(disp);
	free(recv_count);
	free(send_buffer);

	// printf("[%d] Leaving %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);
}

void gather_SMC_data(grid_parms *grid, double *smc_buffer, SMC_cell **smc)
{
	// printf("[%d] ++++++ Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	int chunk_size = grid->num_smc_axially * grid->num_smc_circumferentially * (grid->num_coupling_species_smc + grid->neq_smc);

	// Allocate displacements for MPI comms.
	int *disp = (int*)checked_malloc((int)(grid->num_ranks_write_group * sizeof(int)), SRC_LOC);

	// Allocate counts for MPI comms.
	int *recv_count = (int*)checked_malloc(grid->num_ranks_write_group * sizeof(int), SRC_LOC);

	// printf("[%d] Checkpoint %s:%d\n", grid->universal_rank, __FILE__, __LINE__);

	if(grid->rank_write_group == 0)
	{
		// Compute displacement values.
		for (int i = 0; i < grid->num_ranks_write_group; i++)
		{
			disp[i] = i * chunk_size;
			recv_count[i] = chunk_size;
		}
	}

	// Allocate send buffer for all attributes for each cell.
	double* send_buffer = (double *)checked_malloc(sizeof(double) * chunk_size, SRC_LOC);

	int p = 0;

	for(int i = 0; i < grid->num_smc_axially; i++)
	{
		for(int j = 0; j < grid->num_smc_circumferentially; j++)
		{
			send_buffer[p++] = smc[j + 1][i + 1].vars[smc_Ca];
			send_buffer[p++] = smc[j + 1][i + 1].homo_fluxes[cpl_Ca];
			send_buffer[p++] = smc[j + 1][i + 1].vars[smc_IP3];
			send_buffer[p++] = smc[j + 1][i + 1].homo_fluxes[cpl_IP3];
			send_buffer[p++] = smc[j + 1][i + 1].vars[smc_Vm];
			send_buffer[p++] = smc[j + 1][i + 1].homo_fluxes[cpl_Vm];
			send_buffer[p++] = smc[j + 1][i + 1].vars[smc_SR];
			send_buffer[p++] = smc[j + 1][i + 1].vars[smc_w];
		}
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, chunk_size, MPI_DOUBLE, smc_buffer, recv_count, disp, MPI_DOUBLE, 0, grid->write_group));

	free(disp);
	free(recv_count);
	free(send_buffer);
}
/*
 * Collect JPLC for each cart grid into an array on the writing cores.
 */
void gather_JPLC(grid_parms* grid, double *jplc_buffer, EC_cell** ec, int timestep)
{

	int local_jplc_buffer_size = grid->num_ec_axially * grid->num_ec_circumferentially;

	// Allocate local buffer for jplc values.
	double* local_jplc_buffer = (double *)checked_malloc(sizeof(double) * local_jplc_buffer_size, SRC_LOC);

	// Collect jplc values.
	int seq_count = 0;
	for(int i = 0; i < grid->num_ec_axially; i++)
	{
		for(int j = 0; j < grid->num_ec_circumferentially; j++, seq_count++)
		{
			local_jplc_buffer[seq_count] = ec[j + 1][i + 1].JPLC[timestep];
		}
	}

	// Allocate displacements.
	int *disp = (int*)checked_malloc(grid->num_ranks_write_group * sizeof(int), SRC_LOC);
	// Allocate count values.
	int *recv_count = (int*)checked_malloc(grid->num_ranks_write_group * sizeof(int), SRC_LOC);

	if(grid->rank_write_group == 0)
	{
		// Compute displacement values.
		for (int i = 0; i < grid->num_ranks_write_group; i++)
		{
			disp[i] = i * local_jplc_buffer_size;
			recv_count[i] = local_jplc_buffer_size;
		}
	}

	// Gather all local jplc buffers into the output jplc buffer.
	CHECK_MPI_ERROR(MPI_Gatherv(local_jplc_buffer, local_jplc_buffer_size, MPI_DOUBLE, jplc_buffer, recv_count, disp, MPI_DOUBLE, 0, grid->write_group));

	free(disp);
	free(recv_count);
	free(local_jplc_buffer);

	// printf("[%d] Leaving %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);
}
