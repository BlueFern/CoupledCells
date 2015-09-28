#include <assert.h>
#include "computelib.h"

/**
 * Dump info/debug output to a log file.
 */
void dump_rank_info(conductance cpl_cef, grid_parms grid) //, IO_domain_info* my_IO_domain_info)
{
	MPI_Status status;
	MPI_Offset displacement = 0;
	char* buffer = (char*) checked_malloc(2 * 1024 * sizeof(char), SRC_LOC);
	int root = 0;
	char filename[50];
	int length =
			sprintf(buffer,
					"BRANCH_TAG	= %d\n[Universal_Rank, Cart_Rank= (%d,%d)] \tcoords= %d,%d\t nbrs: local (u,d,l,r)=(%d %d %d %d)\t "
					"remote: (up,down)=(%d %d)\nflip_array: (%d,%d,%d,%d)\n"
					"Boundary_tag = %c\n(T = Top\t B= Bottom\t I=Interior edges of the bifurcation segments, parent or children\t N=Interior of the subdomain)\n"
					"COUPLING COEFFICIENTS\n"
					"Vm_hm_smc=%2.5lf\nVm_hm_ec=%2.5lf\nCa_hm_smc=%2.5lf\nCa_hm_ec=%2.5lf\nIP3_hm_smc=%2.5lf\n"
					"IP3_hm_ec=%2.5lf\nVm_ht_smc=%2.5lf\nVm_ht_ec=%2.5lf\nCa_ht_smc=%2.5lf\nCa_ht_ec=%2.5lf\n"
					"IP3_ht_smc=%2.5lf\nIP3_ht_ec=%2.5lf\n\n"
					"Spatial Gradient info:\nUniform JPLC\t=%2.5lf\nMinimum JPLC\t=%2.5lf\nMaximum JPLC\t=%2.5lf\nGradient\t=%2.5lf\n"
					"Total Tasks=%d\n"
					"Number of grid points in axial direction =%d\n"
					"Number of grid points in circumferential direction =%d\n"
					"Number of ECs per node (axially) =%d\n"
					"Number of SMCs per node (circumferentially) =%d\n"
					"Total ECs on this node =%d\n"
					"Total SMCs on this node =%d\n"
					"Total number of cells on this node =%d\n"
					"\nTotal ECs in the full computational domain =%d\n"
					"Total SMCs in the full computational domain =%d\n"
					"Total number of cells in the full computational domain =%d\n"
					"Total number of equations in the full computational domain =%d\n"
					"------------------------------------------------------------------\n\n",
					grid.branch_tag, grid.universal_rank, grid.rank_branch, grid.coords[0], grid.coords[1], grid.nbrs[local][UP], grid.nbrs[local][DOWN],
					grid.nbrs[local][LEFT], grid.nbrs[local][RIGHT], grid.nbrs[remote][UP], grid.nbrs[remote][DOWN],
					grid.flip_array[0], grid.flip_array[1], grid.flip_array[2], grid.flip_array[3],
					grid.my_domain.internal_info.boundary_tag, cpl_cef.Vm_hm_smc, cpl_cef.Vm_hm_ec, cpl_cef.Ca_hm_smc, cpl_cef.Ca_hm_ec,
					cpl_cef.IP3_hm_smc, cpl_cef.IP3_hm_ec, cpl_cef.Vm_ht_smc, cpl_cef.Vm_ht_ec, cpl_cef.Ca_ht_smc, cpl_cef.Ca_ht_ec,
					cpl_cef.IP3_ht_smc, cpl_cef.IP3_ht_ec, grid.uniform_jplc, grid.min_jplc, grid.max_jplc, grid.gradient, grid.num_ranks, grid.m,
					grid.n, grid.num_ec_axially, grid.num_smc_circumferentially, grid.num_ec_axially * grid.num_ec_circumferentially,
					grid.num_smc_axially * grid.num_smc_circumferentially,
					(grid.num_ec_axially * grid.num_ec_circumferentially) + (grid.num_smc_axially * grid.num_smc_circumferentially),
					(grid.num_ec_circumferentially * grid.num_ec_axially * grid.num_ranks),
					(grid.num_smc_circumferentially * grid.num_smc_axially * grid.num_ranks),
					((grid.num_ec_axially * grid.num_ec_circumferentially) + (grid.num_smc_axially * grid.num_smc_circumferentially)) * grid.num_ranks,
					grid.NEQ * grid.num_ranks);

	int *recv_count = (int*) checked_malloc(grid.num_ranks_branch * sizeof(int), SRC_LOC);
	int *disp = (int*) checked_malloc(grid.num_ranks_branch * sizeof(int), SRC_LOC);

	// Gathering and summing the length of all the CHARs contained in every send_buffer containing coordinates from each MPI process.
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid.cart_comm));

	grid.logfile_displacements = 0;
	for (int i = 0; i < grid.num_ranks_branch; i++)
	{
		disp[i] = grid.logfile_displacements;
		grid.logfile_displacements += recv_count[i];
	}

	if (grid.rank_branch == 0)
	{
		grid.logfile_write_buffer = (char*) checked_malloc(grid.logfile_displacements * sizeof(char), SRC_LOC);
	}
	CHECK_MPI_ERROR(MPI_Gatherv(buffer, length, MPI_CHAR, grid.logfile_write_buffer, recv_count, disp, MPI_CHAR, root, grid.cart_comm));

	if (grid.rank_branch == 0)
	{
		sprintf(filename, "Logfile_%s.txt", grid.suffix);
		MPI_File rank_info_file;
		CHECK_MPI_ERROR(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &rank_info_file));
		CHECK_MPI_ERROR(MPI_File_write_at(rank_info_file, displacement, grid.logfile_write_buffer, grid.logfile_displacements, MPI_CHAR, &status));
		MPI_File_close(&rank_info_file);
		free(grid.logfile_write_buffer);
	}

	free(recv_count);
	free(disp);
}

/// Read data from the config.txt file to retrieve the information related to how the domain is set up.
/// All tasks open the same file to read.
/// Every task has the same displacement so each start to read data from the same position.
/// Each task decides the read buffer size to be allocated by looking at the file size of the file opened for read operation.
/// The data is read and sorted into delimited strings; numbers and stored into corresponding place holders in the array
/// domains[][] in the structure grid_parms grid.
///
///	For a bifurcation as well as a straight segment there's only one domain.
///
/// The arrays in the domains[][] store the following information:
///
/// Element 0: 	Key_val or serial number of the subdomain
/// Element 1:	Subdomain Type (2 possibilities and their values)
/// 				1. Straight Segment (STRSEG): (0)
/// 				2. Bifurcation	(BIF): (1)
/// Element 2: Number of quads/tasks in the axial extent for the current key_val.
/// Element 3: Number of quads/tasks in the circumferential for the current key_val.
/// Element 4: Parent subdomain key_val of current key_val.
/// Element 5: Left child subdomain key_val of the current key_val.
/// Element 6: Right child subdomain key_val of the current key_val.
/// Element 7: Required number of endothelial cells axially per quad/task.
/// Element 8: Required number of smooth muscle cells circumferentially per quad/task.
///
/// In the case of elements 5 and 6, if subdomain type of current key_val is a straight segment,
/// left child is positive or zero, and right child is negative. ???
/// If subdomain type of current key_val is a bifurcation, then both right and left child subdomains are non-negative. ???

void read_config_file(int rank, char* filename, grid_parms* grid)
{
	int err;
	MPI_File input_file;
	MPI_Offset file_size;
	MPI_Status status;
	const char *delimiter = ";,\n";
	char *buffer;
	char *token;
	int *values;

	CHECK_MPI_ERROR(MPI_File_open(grid->universe, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &input_file));

	err = MPI_File_get_size(input_file, &file_size);

	buffer = (char*)checked_malloc((int)file_size * sizeof(char) + 1, SRC_LOC); // Extra char is for the null character.

	err = MPI_File_read_all(input_file, buffer, file_size, MPI_CHAR, &status);
	buffer[file_size] = '\0'; // Null-terminate the damn string!

	// Parse the first value indicating the number of domains.
	token = strtok(buffer, delimiter);

	if(token != NULL)
	{
		grid->num_domains = atoi(token);
	}
	else
	{
		printf("[%d] Unable to read the number of domains value from file %s.\n", grid->universal_rank, filename);
		MPI_Abort(MPI_COMM_WORLD, 911);
	}

	// Allocate the array for the rest of the values in the config file.
	values = (int*)checked_malloc(NUM_CONFIG_ELEMENTS * grid->num_domains * sizeof(int), SRC_LOC);

	// Parse the rest of the values in the config file.
	int index = 0;
	token = strtok(NULL, delimiter);
	while(token != NULL)
	{
		values[index++] = atoi(token);
		token = strtok(NULL, delimiter);
	}

	// Error checking.
	if(index != (NUM_CONFIG_ELEMENTS * grid->num_domains))
	{
		printf("[%d] Insufficient number of values in the config file %s.\n", grid->universal_rank, filename);
		MPI_Abort(MPI_COMM_WORLD, 911);
	}

	// TODO: Check the allocated memory is released when appropriate.

	// Allocate first dimension array.
	grid->domain_params = (int**) checked_malloc(grid->num_domains * sizeof(int*), SRC_LOC);

	// Allocate second dimension arrays.
	for (int i = 0; i < grid->num_domains; i++) {
		grid->domain_params[i] = (int*) checked_malloc(NUM_CONFIG_ELEMENTS * sizeof(int), SRC_LOC);
	}

	// Copy the data into the domains array from the values array.
	for (int i = 0; i < grid->num_domains; i++) {
		for (int j = 0; j < NUM_CONFIG_ELEMENTS; j++) {
			grid->domain_params[i][j] = values[(i * NUM_CONFIG_ELEMENTS) + j];
		}
	}

	MPI_File_close(&input_file);

	free(buffer);
	free(values);
}

// Prepare the suffix which indicates our subdomain information, bifurcation or tube segment suffix, and the containing branch info.
void set_file_naming_strings(grid_parms* grid)
{
	int subdomain, branch, err;

	if(grid->my_domain.internal_info.domain_type == STRSEG)
	{
		subdomain = grid->my_domain.internal_info.domain_index;
		err = sprintf(grid->suffix, "%d", subdomain);
	}
	else if(grid->my_domain.internal_info.domain_type == BIF)
	{
		subdomain = grid->my_domain.internal_info.domain_index;
		branch = grid->branch_tag;
		err = sprintf(grid->suffix, "%d_%d", subdomain, branch);
	}
}

void read_init_ATP(grid_parms *grid, EC_cell **ECs)
{
	int branch;
	if (grid->my_domain.internal_info.domain_type == STRSEG)
	{
		branch = P;
	}
	else if (grid->my_domain.internal_info.domain_type == BIF)
	{
		branch = grid->branch_tag;
	}

	int jplc_per_task_count = grid->num_ec_circumferentially * grid->num_ec_axially;

	int jplc_in_size = jplc_per_task_count * grid->num_ranks_branch;

	// This can be allocated only on the IO nodes.
	double *send_jplc = (double *)checked_malloc(jplc_in_size * sizeof(double), SRC_LOC);

	// Only the IO nodes read the input files.
	if (grid->rank_branch == 0)
	{
		char jplc_file_name[64];
		switch(branch)
		{
			case P:
				sprintf(jplc_file_name, "files/parent_atp.txt");
				break;
			case L:
				sprintf(jplc_file_name, "files/left_daughter_atp.txt");
				break;
			case R:
				sprintf(jplc_file_name, "files/right_daughter_atp.txt");
				break;
			default:
				; // Do something sensible here otherwise all hell breaks loose...
		}

		FILE *fr = fopen(jplc_file_name, "r+");
		printf("Reading ATP from %s, FILE is %s...\n", jplc_file_name, fr == NULL ? "NULL" : "OK");

		int count_in = 0;
		double jplc_val = 0;
		while(fscanf(fr, "%lf", &jplc_val) == 1)
		{
			send_jplc[count_in] = jplc_val;
			count_in++;
		}
		if(feof(fr))
		{
			// Check the number of values read.
			assert(jplc_in_size == count_in);
		}
		fclose(fr);
	}

	int *send_jplc_counts = (int *)checked_malloc(grid->num_ranks_branch * sizeof(int), SRC_LOC);
	int *send_jplc_offsets = (int *)checked_malloc(grid->num_ranks_branch * sizeof(int), SRC_LOC);

	for(int task = 0; task < grid->num_ranks_branch; task++)
	{
		send_jplc_counts[task] = jplc_per_task_count;
		send_jplc_offsets[task] = task * jplc_per_task_count;
	}

	int recv_jplc_count = jplc_per_task_count;
	double *recv_jplc = (double *)checked_malloc(recv_jplc_count * sizeof(double), SRC_LOC);

	int root = 0;

	// printf("%s, grid->cart_comm: %p\n", __FUNCTION__, (void *)grid->cart_comm);

	// Scatter JPLC values to the nodes in this Cartesian grid.
	CHECK_MPI_ERROR(MPI_Scatterv(send_jplc, send_jplc_counts, send_jplc_offsets, MPI_DOUBLE, recv_jplc, recv_jplc_count, MPI_DOUBLE, root, grid->cart_comm));

	// Assign received JPLC values to the cells.
	for(int m = 1; m <= grid->num_ec_circumferentially; m++)
	{
		for(int n = 1; n <= grid->num_ec_axially; n++)
		{
			// Fortran array referencing!
			ECs[m][n].JPLC = recv_jplc[(n - 1) * grid->num_ec_circumferentially + m - 1];
		}
	}

	free(send_jplc);
	free(send_jplc_counts);
	free(send_jplc_offsets);
	free(recv_jplc);
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

void dump_time_field(char* file_prefix, grid_parms grid, double field) //, IO_domain_info* my_IO_domain_info)
{
	MPI_Status status;
	MPI_Offset displacement = 0;
	MPI_File fw;
	char* buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * sizeof(char), SRC_LOC);
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
		sprintf(filename, "%s/%s_%s.txt", grid.time_profiling_dir, file_prefix, grid.suffix);
		CHECK_MPI_ERROR(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fw));
		CHECK_MPI_ERROR(MPI_File_write_at(fw, 0, write_buffer, total_buffer_length, MPI_CHAR, &status));
		MPI_File_close(&fw);
		free(write_buffer);
	}
	free(recv_count);
	free(buffer);
	free(disp);
}
