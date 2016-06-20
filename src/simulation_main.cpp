#include <stdlib.h>
#include <cstring>

#include "computelib.h"
#include "koenigsberger_model.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

void read_config_file(grid_parms* grid);

conductance cpl_cef;
double **sendbuf, **recvbuf;
grid_parms grid;

FILE* var_file;
double* plotttingBuffer;
int bufferPos;

/**
 * The following steps in the ::main function are necessary for setting up
 * the simulation.
 */
int main(int argc, char* argv[])
{
	/// - Global declaration of request and status update place-holders.
	/// Request and status handles for nonblocking send and receive operations,
	/// for communicating with each of the four neighbours.

	MPI_Request reqs[8];
	MPI_Status stats[8];

	/// - Initialise MPI.
	MPI_Init(&argc, &argv);

	grid.universe = MPI_COMM_WORLD;

	/// - Reveal information of myself and size of MPI_COMM_WORLD
	CHECK_MPI_ERROR(MPI_Comm_rank(grid.universe, &grid.universal_rank));
	CHECK_MPI_ERROR(MPI_Comm_size(grid.universe, &grid.num_ranks));

#if PLOTTING && EXPLICIT_ONLY
	if (grid.universal_rank == RANK)
	{
		printf("\nWriting to %s\n\n", FILENAME);
		var_file = fopen(FILENAME, "w");
		double buf[OUTPUT_PLOTTING_SIZE];
		plotttingBuffer = buf;
		bufferPos = 0;
	}
#endif

#ifdef HAVE_OPENMP
#pragma omp parallel
 {
        int thread_id = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        if (thread_id == 0)
	        printf("Running with %d thread(s)\n", num_threads);
 }
#endif

	char filename[50];
	int coupling_case = 1;

	/// Time variables
	double tfinal = 1e-2;
	double interval = 1e-2;
	double data_writing_frequency = 10.00;
	grid.random = 0;

	// Read command line input
	// t - T_END for the simulation
	// w - A number deciding how frequent the data should be recorded, default is every 10 seconds
	// i - Time interval between two steps. Default is 1e-2
	if (argc > 1) {
		for (int i = 0; i < argc; i++) {
			if (argv[i][0] == '-') {
				if (argv[i][1] == 'f') {
					sprintf(grid.config_file, "%s", argv[i + 1]);
				} else if (argv[i][1] == 'C') {
					coupling_case = atoi(argv[i + 1]);
				} else if (argv[i][1] == 'R') {
					grid.random = atoi(argv[i + 1]);
				} else if (argv[i][1] == 'S') {
					sprintf(grid.solution_dir, "%s", argv[i + 1]);
				} else if (argv[i][1] == 'T') {
					sprintf(grid.time_profiling_dir, "%s", argv[i + 1]);
				} else if (argv[i][1] == 't') {
					tfinal = atof(argv[i + 1]);
				} else if (argv[i][1] == 'w') {
					data_writing_frequency = atof(argv[i + 1]);
				} else if (argv[i][1] == 'i') {
					interval = atof(argv[i + 1]);
				}
			}
		}
	}

	if(grid.universal_rank == 0)
	{
		printf("Running with coupling case %d.\n", coupling_case);

#ifdef RK_SUITE
		printf("Running with RK_SUITE solver.\n");
#elif defined ARK_ODE
		printf("Running with ARK_ODE solver.\n");
#elif defined BOOST_ODEINT
		printf("Running with BOOST_ODEINT solver.\n");
#else
		fprintf(stderr, "No solver? How is this possible? Does not compute.\n");
		MPI_Abort(grid.universe, 911);
#endif
	}

	//grid.NO_path = 0;
	//grid.cGMP_path = 0;
	grid.smc_model = KNBGR;
	grid.ec_model = KNBGR;
	grid.uniform_jplc = 0.3;
	grid.stimulus_onset_time = 5.00;

	grid.num_smc_fundblk_circumferentially = 1;
	grid.num_ec_fundblk_circumferentially = 5;
	grid.num_smc_fundblk_axially = 13;
	grid.num_ec_fundblk_axially = 1;
	grid.num_ghost_cells = 2;

	grid.num_fluxes_smc = NUM_FLUXES_SMC;
	grid.num_fluxes_ec = NUM_FLUXES_EC;

	grid.num_coupling_species_smc = NUM_COUPLING_SPECIES_SMC;
	grid.num_coupling_species_ec = NUM_COUPLING_SPECIES_EC;

	grid.neq_smc = NUM_VARS_SMC;
	grid.neq_ec = NUM_VARS_EC;


	// File written every 1 second.
	int file_write_per_unit_time = (int) (data_writing_frequency * int(1 / interval));

	/// - Read domain configuration from input config file.
	read_config_file(&grid);

	/// - Calculate the number of cells per task.
	set_task_parameters(&grid);

	if (grid.domain_type == STRSEG)
	{
		make_straight_cart_grid(&grid);
	}
	else if (grid.domain_type == BIF)
	{
		make_bifucation_cart_grids(&grid);
	}

	/// Now allocate memory for the structures representing the cells and the various members of those structures.
	/// Each of the two cell grids have two additional rows and two additional columns as ghost cells.
	/// Following is an example of a 5x7 grid with added ghost cells on all four sides. the '0's are the actual
	/// members of the grid whereas the '+'s are the ghost cells.

	/**
\verbatim
+ + + + + + + + +
+ 0 0 0 0 0 0 0 +
+ 0 0 0 0 0 0 0 +
+ 0 0 0 0 0 0 0 +
+ 0 0 0 0 0 0 0 +
+ 0 0 0 0 0 0 0 +
+ + + + + + + + +
\endverbatim
*/

	int nc2 = grid.num_smc_circumferentially + grid.num_ghost_cells;
        int na2 = grid.num_smc_axially + grid.num_ghost_cells;
	SMC_cell** smc = (SMC_cell**) checked_malloc(nc2 * sizeof(SMC_cell*), SRC_LOC);
	for (int i = 0; i < nc2; i++)
	{
	        smc[i] = (SMC_cell*) checked_malloc(na2* sizeof(SMC_cell), SRC_LOC);
	}

        nc2 = grid.num_ec_circumferentially + grid.num_ghost_cells;
        na2 = grid.num_ec_axially + grid.num_ghost_cells;
        EC_cell** ec =  (EC_cell**) checked_malloc(nc2 * sizeof(EC_cell*), SRC_LOC);
	for (int i = 0; i < nc2; i++)
	{
		ec[i] = (EC_cell*) checked_malloc(na2* sizeof(EC_cell), SRC_LOC);
	}

        /// To pass the fields to the solvers
        All_cell all_cell;
        all_cell.smc = smc;
        all_cell.ec = ec;

	/// Allocating memory for coupling data to be sent and received through MPI.
	/// sendbuf and recvbuf are 2D arrays with up, down, left and right directions as their first dimension.

	sendbuf = (double**) checked_malloc(4 * sizeof(double*), SRC_LOC);
	recvbuf = (double**) checked_malloc(4 * sizeof(double*), SRC_LOC);

	/// Each processor now allocates the memory for send and recv buffers those will hold the coupling information.
	/// Since sendbuf must contain the information of number of SMCs and ECs being sent in all directions,
	/// the first two elements contain the total count of SMCs located on the rank in the relevant dimension (circumferential or axial)
	/// and the count of SMCs for which information is being sent, respectively.
	/// The next two elements contain the same information for ECs.

	/// Variables to calculate the length of the prospective buffer based on number of cell in either orientations (circumferential or axial).
	int extent_s, extent_e;

	/// Number of elements containing additional information at the beginning of the send buffer.
	grid.added_info_in_send_buf = 4;

	extent_s = grid.num_smc_circumferentially;
	extent_e = grid.num_ec_circumferentially;
	//printf("num_smc_circumferentially: %d\nnum_ec_circumferentially: %d\n\n", extent_s, extent_e);

	/// Data to send to the neighbour in UP direction.
	/// Recording the number of elements in Send buffer in Up direction for use in update routine as count of elements.
	grid.num_elements_send_up = grid.added_info_in_send_buf + (grid.num_coupling_species_smc * extent_s + grid.num_coupling_species_ec * extent_e);

	sendbuf[UP] = (double*) checked_malloc(grid.num_elements_send_up * sizeof(double), SRC_LOC);

	sendbuf[UP][0] = 1.0; //Start of the 1st segment of SMC array in  in UP direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[UP][1] = (double) extent_s; //End of the 1st segment of SMC array in UP direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[UP][2] = 1.0; //Start of the 1st segment of EC array in UP direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[UP][3] = (double) extent_e; ///End of the 1st segment of EC array in UP direction (circumferential direction) to be sent to neighbouring processor

	/// Data to send to the neighbour in  direction.
	/// Recording the number of elements in Send buffer in DOWN direction for use in update routine as count of elements.
	grid.num_elements_send_down = grid.added_info_in_send_buf + (grid.num_coupling_species_smc * extent_s + grid.num_coupling_species_ec * extent_e);

	sendbuf[DOWN] = (double*) checked_malloc(grid.num_elements_send_down * sizeof(double), SRC_LOC);

	sendbuf[DOWN][0] = 1.0; //Start of the 1st segment of SMC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[DOWN][1] = (double) extent_s; //End of the 1st segment of SMC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[DOWN][2] = 1.0; //Start of the 1st segment of EC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[DOWN][3] = (double) extent_e; ///End of the 1st segment of EC array in DOWN direction (circumferential direction) to be sent to neighbouring processor

	extent_s = grid.num_smc_axially;
	extent_e = grid.num_ec_axially;

	/// Data to send to the neighbour in LEFT direction.
	/// Recording the number of elements in Send buffer in Left direction for use in update routine as count of elements.
	grid.num_elements_send_left = grid.added_info_in_send_buf + (grid.num_coupling_species_smc * extent_s + grid.num_coupling_species_ec * extent_e);
	sendbuf[LEFT] = (double*) checked_malloc(grid.num_elements_send_left * sizeof(double), SRC_LOC);

	sendbuf[LEFT][0] = 1.0; //Start of the 1st segment of SMC array in  in LEFT direction (axial direction) to be sent to neighbouring processor
	sendbuf[LEFT][1] = (double) extent_s; //END of the 1st segment of SMC array in  in LEFT direction (axial direction) to be sent to neighbouring processor
	sendbuf[LEFT][2] = 1.0; //Start of the 1st segment of EC array in  in LEFT direction (axial direction) to be sent to neighbouring processor
	sendbuf[LEFT][3] = (double) extent_e; //END of the 1st segment of EC array in  in LEFT direction (axial direction) to be sent to neighbouring processor

	/// Data to send to the neighbour in RIGHT direction.
	/// Recording the number of elements in Send buffer in RIGHT direction for use in update routine as count of elements.
	grid.num_elements_send_right = grid.added_info_in_send_buf + (grid.num_coupling_species_smc * extent_s + grid.num_coupling_species_ec * extent_e);

	sendbuf[RIGHT] = (double*) checked_malloc(grid.num_elements_send_right * sizeof(double), SRC_LOC);

	sendbuf[RIGHT][0] = 1.0; //Start of the 1st segment of SMC array in  in RIGHT direction (axial direction) to be sent to neighbouring processor
	sendbuf[RIGHT][1] = (double) extent_s; //END of the 1st segment of SMC array in  in RIGHT direction (axial direction) to be sent to neighbouring processor
	sendbuf[RIGHT][2] = 1.0; //Start of the 1st segment of EC array in  in RIGHT direction (axial direction) to be sent to neighbouring processor
	sendbuf[RIGHT][3] = (double) extent_e; //END of the 1st segment of EC array in  in RIGHT direction (axial direction) to be sent to neighbouring processor

	/// Tell neighbours the number of elements to be sent to them and allocate memory of recvbuf for each direction accordingly.
	communication_update_recv_size(&grid);

	/// data to receive from the neighbour in UP direction
	recvbuf[UP] = (double*) checked_malloc(grid.num_elements_recv_up * sizeof(double), SRC_LOC);

	/// data to recv from the neighbour in DOWN direction
	recvbuf[DOWN] = (double*) checked_malloc(grid.num_elements_recv_down * sizeof(double), SRC_LOC);

	/// data to receive from the neighbour in LEFT direction
	recvbuf[LEFT] = (double*) checked_malloc(grid.num_elements_recv_left * sizeof(double), SRC_LOC);

	/// data to receive from the neighbour in RIGHT direction
	recvbuf[RIGHT] = (double*) checked_malloc(grid.num_elements_recv_right * sizeof(double), SRC_LOC);

	/// Setting up the solver.
	double tnow = 0.0;

	// Error control variables.
	double TOL = 1e-6, absTOL = 1e-7;
	double *thres = (double*)checked_malloc(grid.NEQ * sizeof(double), SRC_LOC);
	for (int i = 0; i < grid.NEQ; i++)
	{
		thres[i] = absTOL;
	}

	// State variables.
	double* y = (double*) checked_malloc(grid.NEQ * sizeof(double), SRC_LOC);

	// Internal work space.
	double* yp = (double*) checked_malloc(grid.NEQ * sizeof(double), SRC_LOC);

	/// Initialise state variables and coupling data values.

	initialize_koenigsberger_smc(grid, y, smc);
	initialize_koenigsberger_ec(grid, y, ec);

	// Reverse mapping from state vector to cells.
	// Essential for restarts, when data is loaded from a checkpoint.
	int err = map_solver_output_to_cells(grid, y, smc, ec);

	if (err != 0) {
		printf("[%d] error in mapping y to cells\n", grid.universal_rank);
		MPI_Abort(grid.universe, 1000);
	}

	// Initialising the coupling coefficients to be used in the ODEs.
	set_coupling_parms(coupling_case, &cpl_cef);

	// Debug output.
	dump_rank_info(cpl_cef, grid);

	// This is read in here for validation purposes in the output.
	// the solver will reset JPLC and read later it when the time is right.
	read_init_ATP(&grid, ec);

	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("[%d] %s:%d\n", grid.universal_rank, __FILE__, __LINE__);

// TODO: Move yp to solver_rk.

#ifdef RK_SUITE

	rksuite_solver_CT(tnow, tfinal, interval, y, yp, grid.NEQ, TOL, thres, file_write_per_unit_time, grid.solution_dir);

#elif defined ARK_ODE

	arkode_solver(tnow, tfinal, interval, y, grid.NEQ, TOL, absTOL, file_write_per_unit_time, grid.solution_dir, all_cell);

#elif defined BOOST_ODEINT

	odeint_solver(tnow, tfinal, interval, y, grid.NEQ, TOL, absTOL, file_write_per_unit_time, grid.solution_dir);

#else
#error ODE solver not selected. Use -DRK_SUITE | -DARK_ODE | -DBOOST_ODEINT during compilation
#endif

	// Free grid->domain_params.
	for(int domain_num = 0; domain_num < grid.num_domains; domain_num++)
	{
		free(grid.domain_params[domain_num]);
	}
	free(grid.domain_params);


	// Free SMCs.
        nc2 = grid.num_smc_circumferentially + grid.num_ghost_cells;
        for (int i = 0; i < nc2; i++) 
	{
	        free(smc[i]);
	}
	free(smc);

	// Free ECs.
        nc2 = grid.num_ec_circumferentially + grid.num_ghost_cells;
        for (int i = 0; i < nc2; i++) 
	{
	        free(ec[i]);
	}
	free(ec);

	// Free send/receive buffers.
	for (int i = 0; i < 4; i++)
	{
		free(sendbuf[i]);
		free(recvbuf[i]);
	}
	free(sendbuf);
	free(recvbuf);

	// TODO: Do all solvers have the same arrays or do we need to have some ifdefs here?
	// Free the solver-related arrays.
	free(thres);
	free(y);
	free(yp);

#if PLOTTING && EXPLICIT_ONLY
	if (grid.universal_rank == RANK)
	{
		fclose(var_file);
	}
#endif

	MPI_Finalize();

	return EXIT_SUCCESS;
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

void read_config_file(grid_parms* grid)
{
	char line[1024];
	char *token;
	FILE *input_file;
	const char *delims = ";,\n";

	// Open config file for reading.
	input_file = fopen(grid->config_file, "r");
	if(input_file == NULL)
	{
		fprintf(stderr, "[%d] Unable to open file %s.\n", grid->universal_rank, grid->config_file);
		MPI_Abort(MPI_COMM_WORLD, 911);
	}

	// Read the number of domains value from the first line.
	if(fgets(line, 1024, input_file) != NULL)
	{
		sscanf(line, "%d", &grid->num_domains);
	}
	else
	{
		fprintf(stderr, "[%d] Unable to read the number of domains value from file %s.\n", grid->universal_rank, grid->config_file);
		MPI_Abort(MPI_COMM_WORLD, 911);
	}

	// Allocate first dimension array.
	grid->domain_params = (int**)checked_malloc(grid->num_domains * sizeof(int*), SRC_LOC);

	// Read domain parameters line-by-line.
	for(int domain_num = 0; domain_num < grid->num_domains; domain_num++)
	{

		if(fgets(line, 1024, input_file) != NULL)
		{
			// Allocate second dimension array.
			grid->domain_params[domain_num] = (int*) checked_malloc(NUM_CONFIG_ELEMENTS * sizeof(int), SRC_LOC);

			// Get first token.
			token = strtok(line, delims);

			for(int elem_num = 0; elem_num < NUM_CONFIG_ELEMENTS; elem_num++)
			{
				// Parse the token and remember the value.
				grid->domain_params[domain_num][elem_num] = atoi(token);

				// Get next token.
				token = strtok(NULL, delims);
			}
		}
		else
		{
			fprintf(stderr, "[%d] Unable to read domains parameters for domain num %d from file %s.\n", grid->universal_rank, domain_num, grid->config_file);
			MPI_Abort(MPI_COMM_WORLD, 911);
		}
	}

	fclose(input_file);
}

