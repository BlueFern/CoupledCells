#include <malloc.h>

#include "computelib.h"
#include "koenigsberger_model.h"
#include "gather.h"

ec_data_buffer *allocate_EC_data_buffer(int tasks_per_branch, int elements_per_task, int deep)
{
	ec_data_buffer *ec_buffer = (ec_data_buffer *)checked_malloc(sizeof(ec_data_buffer), SRC_LOC);

	ec_buffer->num_chunks = tasks_per_branch;
	ec_buffer->chunk_size = elements_per_task;

	if(deep)
	{
		int buffer_size = ec_buffer->num_chunks * ec_buffer->chunk_size;

		// This is not the neatest way of allocating (and, consequently, releasing memory).
		// Perhaps we can use a double ** array here, where the elements are accessed with enumerated constants.
		ec_buffer->_ec_Ca = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		ec_buffer->_ec_SR = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		ec_buffer->_ec_Vm = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		ec_buffer->_ec_IP3 = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		ec_buffer->_ec_cpl_Ca = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		ec_buffer->_ec_cpl_Vm = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		ec_buffer->_ec_cpl_IP3 = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
	}

	return ec_buffer;
}

void free_EC_data_buffer(ec_data_buffer *ec_buffer, int deep)
{
	if(deep)
	{
		free(ec_buffer->_ec_Ca);
		free(ec_buffer->_ec_SR);
		free(ec_buffer->_ec_Vm);
		free(ec_buffer->_ec_IP3);
		free(ec_buffer->_ec_cpl_Ca);
		free(ec_buffer->_ec_cpl_Vm);
		free(ec_buffer->_ec_cpl_IP3);
	}

	free(ec_buffer);

	return;
}

void gather_EC_data(grid_parms *grid, ec_data_buffer *ec_buffer, EC_cell **ec)
{
	// printf("[%d] ++++++ Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	// Allocate displacements for MPI comms.
	int *disp = (int*)checked_malloc((int)(grid->num_ranks_branch * sizeof(int)), SRC_LOC);

	int chunk_size = grid->num_ec_axially * grid->num_ec_circumferentially * (grid->neq_ec + grid->num_coupling_species_ec);
	double *recv_buffer = (double*)checked_malloc(grid->num_ranks_branch * chunk_size * sizeof(double), SRC_LOC);

	// Allocate counts for MPI comms.
	int *recv_count = (int*)checked_malloc(grid->num_ranks_branch * sizeof(int), SRC_LOC);

	// printf("[%d] Checkpoint %s:%d\n", grid->universal_rank, __FILE__, __LINE__);

	if(grid->rank_branch == 0)
	{
		// Compute displacement values.
		for (int i = 0; i < grid->num_ranks_branch; i++)
		{
			disp[i] = i * chunk_size;
			recv_count[i] = chunk_size;
		}
	}

	// Allocate send buffer.
	double* send_buffer = (double *)checked_malloc(sizeof(double) * chunk_size, SRC_LOC);

	for(int i = 0, p = 0; i < grid->num_ec_axially; i++)
	{
		for(int j = 0; j < grid->num_ec_circumferentially; j++, p += (grid->neq_ec + grid->num_coupling_species_ec))
		{
			send_buffer[p] = ec[j + 1][i + 1].vars[ec_Ca];
			send_buffer[p + 1] = ec[j + 1][i + 1].vars[ec_IP3];
			send_buffer[p + 2] = ec[j + 1][i + 1].vars[ec_SR];
			send_buffer[p + 3] = ec[j + 1][i + 1].vars[ec_Vm];
			send_buffer[p + 4] = ec[j + 1][i + 1].homo_fluxes[cpl_Ca];
			send_buffer[p + 5] = ec[j + 1][i + 1].homo_fluxes[cpl_IP3];
			send_buffer[p + 6] = ec[j + 1][i + 1].homo_fluxes[cpl_Vm];
		}
	}



	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, chunk_size, MPI_DOUBLE, recv_buffer, recv_count, disp, MPI_DOUBLE, 0, grid->cart_comm));

	if(grid->rank_branch == 0)
	{
		// write recv buffer elements to their appropriate locations in ec_buffer
		for (int i = 0, j = 0; i < grid->num_ec_axially * grid->num_ec_circumferentially; i++)
		{
			ec_buffer->_ec_Ca[i] = recv_buffer[j++];
			ec_buffer->_ec_IP3[i] = recv_buffer[j++];
			ec_buffer->_ec_SR[i] = recv_buffer[j++];
			ec_buffer->_ec_Vm[i] = recv_buffer[j++];

			ec_buffer->_ec_cpl_Ca[i] = recv_buffer[j++];
			ec_buffer->_ec_cpl_IP3[i] = recv_buffer[j++];
			ec_buffer->_ec_cpl_Vm[i] = recv_buffer[j++];

		}
	}

	free(recv_buffer);
	free(disp);
	free(recv_count);
	free(send_buffer);

	// printf("[%d] Leaving %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);
}


/* Allocate memory for the smc_data_buffer. */
smc_data_buffer *allocate_SMC_data_buffer(int tasks_per_branch, int elements_per_task, int deep)
{
	smc_data_buffer *smc_buffer = (smc_data_buffer *)checked_malloc(sizeof(smc_data_buffer), SRC_LOC);

	smc_buffer->num_chunks = tasks_per_branch;
	smc_buffer->chunk_size = elements_per_task;

	if(deep)
	{
		int buffer_size = smc_buffer->num_chunks * smc_buffer->chunk_size;

		// This is not the neatest way of allocating (and, consequently, releasing memory).
		// Perhaps we can use a double ** array here, where the elements are accessed with enumerated constants.
		smc_buffer->_smc_Ca = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		smc_buffer->_smc_SR = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		smc_buffer->_smc_W = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		smc_buffer->_smc_Vm = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		smc_buffer->_smc_IP3 = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		smc_buffer->_smc_cpl_Ca = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		smc_buffer->_smc_cpl_Vm = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
		smc_buffer->_smc_cpl_IP3 = (double *)checked_malloc(sizeof(double) * buffer_size, SRC_LOC);
	}

	return smc_buffer;
}

/* Free the memory in the smc_data_buffer. */
void free_SMC_data_buffer(smc_data_buffer *smc_buffer, int deep)
{
	if(deep)
	{
		free(smc_buffer->_smc_Ca);
		free(smc_buffer->_smc_SR);
		free(smc_buffer->_smc_W);
		free(smc_buffer->_smc_Vm);
		free(smc_buffer->_smc_IP3);
		free(smc_buffer->_smc_cpl_Ca);
		free(smc_buffer->_smc_cpl_Vm);
		free(smc_buffer->_smc_cpl_IP3);
	}

	free(smc_buffer);

	return;
}

/* Collect state variables from the SMC cells in one branch. */
void gather_SMC_data(grid_parms *grid, smc_data_buffer *smc_buffer, SMC_cell **smc)
{
	// printf("[%d] ++++++ Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	int chunk_size = grid->num_smc_axially * grid->num_smc_circumferentially * (grid->neq_smc + grid->num_coupling_species_smc);
	double *recv_buffer = (double*)checked_malloc(grid->num_ranks_branch * chunk_size * sizeof(double), SRC_LOC);

	// Allocate displacements for MPI comms.
	int *disp = (int*)checked_malloc((int)(grid->num_ranks_branch * sizeof(int)), SRC_LOC);

	// Allocate counts for MPI comms.
	int *recv_count = (int*)checked_malloc(grid->num_ranks_branch * sizeof(int), SRC_LOC);

	// printf("[%d] Checkpoint %s:%d\n", grid->universal_rank, __FILE__, __LINE__);

	if(grid->rank_branch == 0)
	{
		// Compute displacement values.
		for (int i = 0; i < grid->num_ranks_branch; i++)
		{
			disp[i] = i * chunk_size;
			recv_count[i] = chunk_size;
		}
	}

	// Allocate send buffer for values.
	double* send_buffer = (double *)checked_malloc(sizeof(double) * chunk_size, SRC_LOC);

	for(int i = 0, p = 0; i < grid->num_smc_axially; i++)
	{
		for(int j = 0; j < grid->num_smc_circumferentially; j++, p += (grid->neq_smc + grid->num_coupling_species_smc))
		{
			send_buffer[p] = smc[j + 1][i + 1].vars[smc_Ca];
			send_buffer[p + 1] = smc[j + 1][i + 1].vars[smc_IP3];
			send_buffer[p + 2] = smc[j + 1][i + 1].vars[smc_SR];
			send_buffer[p + 3] = smc[j + 1][i + 1].vars[smc_Vm];
			send_buffer[p + 4] = smc[j + 1][i + 1].vars[smc_w];
			send_buffer[p + 5] = smc[j + 1][i + 1].homo_fluxes[cpl_Ca];
			send_buffer[p + 6] = smc[j + 1][i + 1].homo_fluxes[cpl_IP3];
			send_buffer[p + 7] = smc[j + 1][i + 1].homo_fluxes[cpl_Vm];
		}
	}

	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, chunk_size, MPI_DOUBLE, recv_buffer, recv_count, disp, MPI_DOUBLE, 0, grid->cart_comm));

	if(grid->rank_branch == 0)
	{
		// write recv buffer elements to their appropriate locations in ec_buffer
		for (int i = 0, j = 0; i < grid->num_smc_axially * grid->num_smc_circumferentially; i++)
		{
			smc_buffer->_smc_Ca[i] = recv_buffer[j++];
			smc_buffer->_smc_IP3[i] = recv_buffer[j++];
			smc_buffer->_smc_SR[i] = recv_buffer[j++];
			smc_buffer->_smc_Vm[i] = recv_buffer[j++];
			smc_buffer->_smc_W[i] = recv_buffer[j++];
			smc_buffer->_smc_cpl_Ca[i] = recv_buffer[j++];
			smc_buffer->_smc_cpl_IP3[i] = recv_buffer[j++];
			smc_buffer->_smc_cpl_Vm[i] = recv_buffer[j++];
		}
	}

	free(recv_buffer);
	free(disp);
	free(recv_count);
	free(send_buffer);
}
/*
 * Collect JPLC for each cart grid into an array on the writing cores.
 */
void gather_JPLC(grid_parms* grid, double *jplc_buffer, EC_cell** ec)
{
	// printf("[%d] ++++++ Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	// printf("%s, %d, [%d]\n", __FILE__, __LINE__, grid->rank / grid->num_ranks_branch);

	int local_jplc_buffer_size = grid->num_ec_axially * grid->num_ec_circumferentially;

	// Allocate local buffer for jplc values.
	double* local_jplc_buffer = (double *)checked_malloc(sizeof(double) * local_jplc_buffer_size, SRC_LOC);

	// Collect jplc values.
	int seq_count = 0;
	for(int i = 0; i < grid->num_ec_axially; i++)
	{
		for(int j = 0; j < grid->num_ec_circumferentially; j++, seq_count++)
		{
			local_jplc_buffer[seq_count] = ec[j + 1][i + 1].JPLC;
		}
	}

	// Allocate displacements.
	int *disp = (int*)checked_malloc(grid->num_ranks_branch * sizeof(int), SRC_LOC);
	// Allocate count values.
	int *recv_count = (int*)checked_malloc(grid->num_ranks_branch * sizeof(int), SRC_LOC);

	if(grid->rank_branch == 0)
	{
		// Compute displacement values.
		for (int i = 0; i < grid->num_ranks_branch; i++)
		{
			disp[i] = i * local_jplc_buffer_size;
			recv_count[i] = local_jplc_buffer_size;
		}
	}

	// Gather all local jplc buffers into the output jplc buffer.
	CHECK_MPI_ERROR(MPI_Gatherv(local_jplc_buffer, local_jplc_buffer_size, MPI_DOUBLE, jplc_buffer, recv_count, disp, MPI_DOUBLE, 0, grid->cart_comm));

	free(disp);
	free(recv_count);
	free(local_jplc_buffer);

	// printf("[%d] Leaving %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);
}

