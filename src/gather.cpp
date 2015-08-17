#include "macros.h"
#include "computelib.h"
#include "gather.h"


ec_data_buffer *allocate_EC_data_buffer(int tasks_per_branch, int elements_per_task, int deep)
{
	ec_data_buffer *ec_buffer = (ec_data_buffer *)malloc(sizeof(ec_data_buffer));

	ec_buffer->num_chunks = tasks_per_branch;
	ec_buffer->chunk_size = elements_per_task;

	if(!deep)
	{
		return ec_buffer;
	}

	int buffer_size = ec_buffer->num_chunks * ec_buffer->chunk_size;

	// This is not the neatest way of allocating (and, consequently, releasing memory).
	// Perhaps we can use a double ** array here, where the elements are accessed with enumerated constants.
	ec_buffer->_ec_Ca = (double *)malloc(sizeof(double) * buffer_size);
	ec_buffer->_ec_SR = (double *)malloc(sizeof(double) * buffer_size);
	ec_buffer->_ec_Vm = (double *)malloc(sizeof(double) * buffer_size);
	ec_buffer->_ec_IP3 = (double *)malloc(sizeof(double) * buffer_size);
	ec_buffer->_ec_cpl_Ca = (double *)malloc(sizeof(double) * buffer_size);
	ec_buffer->_ec_cpl_Vm = (double *)malloc(sizeof(double) * buffer_size);
	ec_buffer->_ec_cpl_IP3 = (double *)malloc(sizeof(double) * buffer_size);

	return ec_buffer;
}

void free_EC_data_buffer(ec_data_buffer *ec_buffer)
{
	free(ec_buffer->_ec_Ca);
	free(ec_buffer->_ec_SR);
	free(ec_buffer->_ec_Vm);
	free(ec_buffer->_ec_IP3);
	free(ec_buffer->_ec_cpl_Ca);
	free(ec_buffer->_ec_cpl_Vm);
	free(ec_buffer->_ec_cpl_IP3);

	free(ec_buffer);

	return;
}

void gather_EC_data(grid_parms *grid, ec_data_buffer *ec_buffer, EC_cell **ec)
{

	printf("[%d] Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	// Allocate displacements for MPI comms.
	int *disp = (int*)checked_malloc((int)(grid->tasks * sizeof(int)), SRC_LOC);

	// Allocate counts for MPI comms.
	int *recv_count = (int*)checked_malloc(grid->tasks * sizeof(int), SRC_LOC);

	// printf("[%d] Checkpoint %s:%d\n", grid->universal_rank, __FILE__, __LINE__);

	// Allocate send buffer for values.
	double* send_buffer = (double *)checked_malloc(sizeof(double) * ec_buffer->chunk_size, SRC_LOC);

	// printf("[%d] Checkpoint %s:%d\n", grid->universal_rank, __FILE__, __LINE__);

	if(grid->rank / grid->tasks == 0)
	{
		// Compute displacement values.
		for (int i = 0; i < grid->tasks; i++)
		{
			disp[i] = i * ec_buffer->chunk_size;
			recv_count[i] = ec_buffer->chunk_size;
		}
	}

	// printf("[%d] Checkpoint %s:%d\n", grid->universal_rank, __FILE__, __LINE__);

	/************* ec_Ca field data *************/
	for(int i = 0, p = 0; i < grid->num_ec_axially; i++)
	{
		for(int j = 0; j < grid->num_ec_circumferentially; j++, p++)
		{
			send_buffer[p] = ec[j + 1][i + 1].q[ec_Ca];
		}
	}
	check_flag(MPI_Gatherv(send_buffer, ec_buffer->chunk_size, MPI_DOUBLE, ec_buffer->_ec_Ca, recv_count, disp, MPI_DOUBLE, grid->rank / grid->tasks, grid->split_comm), SRC_LOC);

	/************* ec_cpl_Ca field data *************/
	for(int i = 0, p = 0; i < grid->num_ec_axially; i++)
	{
		for(int j = 0; j < grid->num_ec_circumferentially; j++, p++)
		{
			send_buffer[p] = ec[j + 1][i + 1].B[cpl_Ca];
		}
	}
	check_flag(MPI_Gatherv(send_buffer, ec_buffer->chunk_size, MPI_DOUBLE, ec_buffer->_ec_cpl_Ca, recv_count, disp, MPI_DOUBLE, grid->rank / grid->tasks, grid->split_comm), SRC_LOC);

	/************* ec_IP3 field data *************/
	for(int i = 0, p = 0; i < grid->num_ec_axially; i++)
	{
		for(int j = 0; j < grid->num_ec_circumferentially; j++, p++)
		{
			send_buffer[p] = ec[j + 1][i + 1].q[ec_IP3];
		}
	}
	check_flag(MPI_Gatherv(send_buffer, ec_buffer->chunk_size, MPI_DOUBLE, ec_buffer->_ec_IP3, recv_count, disp, MPI_DOUBLE, grid->rank / grid->tasks, grid->split_comm), SRC_LOC);

	/************* ec_cpl_IP3 field data *************/
	for(int i = 0, p = 0; i < grid->num_ec_axially; i++)
	{
		for(int j = 0; j < grid->num_ec_circumferentially; j++, p++)
		{
			send_buffer[p] = ec[j + 1][i + 1].B[cpl_IP3];
		}
	}
	check_flag(MPI_Gatherv(send_buffer, ec_buffer->chunk_size, MPI_DOUBLE, ec_buffer->_ec_cpl_IP3, recv_count, disp, MPI_DOUBLE, grid->rank / grid->tasks, grid->split_comm), SRC_LOC);

	/************* ec_SR field data *************/
	for(int i = 0, p = 0; i < grid->num_ec_axially; i++)
	{
		for(int j = 0; j < grid->num_ec_circumferentially; j++, p++)
		{
			send_buffer[p] = ec[j + 1][i + 1].q[ec_SR];
		}
	}
	check_flag(MPI_Gatherv(send_buffer, ec_buffer->chunk_size, MPI_DOUBLE, ec_buffer->_ec_SR, recv_count, disp, MPI_DOUBLE, grid->rank / grid->tasks, grid->split_comm), SRC_LOC);

	/************* ec_Vm field data *************/
	for(int i = 0, p = 0; i < grid->num_ec_axially; i++)
	{
		for(int j = 0; j < grid->num_ec_circumferentially; j++, p++)
		{
			send_buffer[p] = ec[j + 1][i + 1].q[ec_Vm];
		}
	}
	check_flag(MPI_Gatherv(send_buffer, ec_buffer->chunk_size, MPI_DOUBLE, ec_buffer->_ec_Vm, recv_count, disp, MPI_DOUBLE, grid->rank / grid->tasks, grid->split_comm), SRC_LOC);

	/************* ec_cpl_Vm field data *************/
	for(int i = 0, p = 0; i < grid->num_ec_axially; i++)
	{
		for(int j = 0; j < grid->num_ec_circumferentially; j++, p++)
		{
			send_buffer[p] = ec[j + 1][i + 1].B[cpl_Vm];
		}
	}
	check_flag(MPI_Gatherv(send_buffer, ec_buffer->chunk_size, MPI_DOUBLE, ec_buffer->_ec_cpl_Vm, recv_count, disp, MPI_DOUBLE, grid->rank / grid->tasks, grid->split_comm), SRC_LOC);

	free(disp);
	free(recv_count);
	free(send_buffer);

	// printf("[%d] Leaving %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);
}


