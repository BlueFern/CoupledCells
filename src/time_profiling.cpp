#include <malloc.h>
#include <stdlib.h> // free

#include "computelib.h"

void initialize_t_stamp(time_stamps* t_stamp)
{
	t_stamp->aggregate_compute = 0;
	t_stamp->aggregate_comm = 0;
	t_stamp->aggregate_ec_gather = 0;
	t_stamp->aggregate_smc_gather = 0;
	t_stamp->aggregate_ec_write = 0;
	t_stamp->aggregate_smc_write = 0;
}

void dump_time_profiling(grid_parms grid, time_stamps* t_stamp)
{
	dump_time_field((char *)"aggregated_compute", grid, t_stamp->aggregate_compute);
	dump_time_field((char *)"aggregated_comm", grid, t_stamp->aggregate_comm);
	dump_time_field((char *)"aggregated_ec_gather", grid, t_stamp->aggregate_ec_gather);
	dump_time_field((char *)"aggregated_smc_gather", grid, t_stamp->aggregate_smc_gather);
	dump_time_field((char *)"aggregated_ec_write", grid, t_stamp->aggregate_ec_write);
	dump_time_field((char *)"aggregated_smc_write", grid, t_stamp->aggregate_smc_write);
}

#define NUM_DBL_TO_CHAR_BYTES 64

void dump_time_field(char* file_prefix, grid_parms grid, double field)
{
	MPI_Status status;
	MPI_Offset displacement = 0;
	MPI_File fw;
	char* buffer = (char*)checked_malloc(NUM_DBL_TO_CHAR_BYTES * sizeof(char), SRC_LOC);
	char* write_buffer;
	int root = 0;
	char filename[50];

	int length = sprintf(buffer, "%2.12lf\n", field);

	int *recv_count = (int*) checked_malloc(grid.num_ranks_branch * sizeof(int), SRC_LOC);
	int *disp = (int*) checked_malloc(grid.num_ranks_branch * sizeof(int), SRC_LOC);

	/// Gathering the lengths of buffer from each MPI process.
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid.cart_comm));

	int total_buffer_length = 0;
	for (int i = 0; i < grid.num_ranks_branch; i++)
	{
		disp[i] = total_buffer_length;
		total_buffer_length += recv_count[i];
	}

	if (grid.rank_branch == 0)
	{
		write_buffer = (char*) checked_malloc(total_buffer_length * sizeof(char), SRC_LOC);
	}

	// Gathering the buffers from all MPI processes.
	CHECK_MPI_ERROR(MPI_Gatherv(buffer, length, MPI_CHAR, write_buffer, recv_count, disp, MPI_CHAR, root, grid.cart_comm));

	if (grid.rank_branch == 0)
	{
		sprintf(filename, "%s/%s_%d_%d.txt", grid.time_profiling_dir, file_prefix, grid.domain_index, grid.branch_tag);

		CHECK_MPI_ERROR(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fw));
		CHECK_MPI_ERROR(MPI_File_write_at(fw, 0, write_buffer, total_buffer_length, MPI_CHAR, &status));
		MPI_File_close(&fw);
	}

	if (grid.rank_branch == 0)
	{
		free(write_buffer);
	}

	free(recv_count);
	free(buffer);
	free(disp);
}
