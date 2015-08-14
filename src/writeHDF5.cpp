#include "writeHDF5.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define SRC_LOC __FILE__ ":" TOSTRING(__LINE__)

#define _2D 2

/*
 * Collect JPLC into an array on each core.
 */
void gather_JPLC(grid_parms* grid, double *jplc_buffer, EC_cell** ec)
{
	//printf("[%d] Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	int root = 0;
	int local_jplc_buffer_size = grid->num_ec_axially * grid->num_ec_circumferentially;

	// Allocate local buffer for jplc values.
	double* local_jplc_buffer = (double *)checked_malloc(sizeof(double) * local_jplc_buffer_size, SRC_LOC);

	// Collect jplc values.
	for(int i = 0; i < grid->num_ec_circumferentially; i++)
	{
		for(int j = 0; j < grid->num_ec_axially; j++)
		{
			local_jplc_buffer[i * grid->num_ec_axially + j] = ec[i + 1][j + 1].JPLC;
		}
	}

	// Allocate displacements.
	int *disp = 0;
	// Allocate count values.
	int *recv_count = (int*)checked_malloc(grid->numtasks * sizeof(int), SRC_LOC);

	if(grid->universal_rank == 0)
	{
		disp = (int*)checked_malloc(grid->numtasks * sizeof(int), SRC_LOC);
		recv_count = (int*)checked_malloc(grid->numtasks * sizeof(int), SRC_LOC);

		// Compute displacement values.
		for (int i = 0; i < grid->numtasks; i++)
		{
			disp[i] = i * local_jplc_buffer_size;
			recv_count[i] = local_jplc_buffer_size;
		}
	}

	// Gather all local jplc buffers into the output jplc buffer.
	check_flag(MPI_Gatherv(local_jplc_buffer, local_jplc_buffer_size, MPI_DOUBLE, jplc_buffer, recv_count, disp, MPI_DOUBLE, root, grid->universe), SRC_LOC);

	free(disp);
	free(recv_count);
	free(local_jplc_buffer);

	//printf("[%d] Leaving %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);
}

void write_HDF5_JPLC(grid_parms* grid, double *jplc_buffer, char *path)
{
	printf("[%d] >>>>>> Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	char filename[256];
	int err = sprintf(filename, "%s/jplc_%d.h5", path, grid->branch_tag);

	printf("JPLC file: %s\n", filename);

	hid_t file_id;
	hid_t space_id;
	hid_t dset_id;

#if 0

	// Does this number exist somewhere in the grid?
	int num_branches = (grid->numtasks / (grid->m * grid->n));
	printf("num_branches %d, m: %d, n: %d\n", num_branches, grid->m, grid->n);

	hsize_t dims[_2D] = {grid->num_ec_axially * grid->m * num_branches, grid->num_ec_circumferentially * grid->n};

	printf("dimensions: %d x %d\n", grid->num_ec_axially * grid->m * num_branches, grid->num_ec_circumferentially * grid->n);
	herr_t status;

	// Create a HDF5 file.
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// Create dataspace.
	space_id = H5Screate_simple(_2D, dims, NULL);

	dset_id = H5Dcreate(file_id, "/jplc", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Create and write an integer type dataset named "jplc".
	//status = H5LTmake_dataset(file_id,"/jplc", _2D, dims, H5T_NATIVE_DOUBLE, jplc_buffer);

	status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, jplc_buffer);

	// Close everything.
	status = H5Dclose(dset_id);
	status = H5Sclose(space_id);
	status = H5Fclose(file_id);
#endif
}
