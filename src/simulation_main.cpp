#include "computelib.h"

using namespace std;

conductance cpl_cef;
SMC_cell **smc;
EC_cell **ec;
double **sendbuf, **recvbuf;
grid_parms grid;
time_keeper elps_t;

int CASE = 1;

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

	elps_t.t_old = MPI_Wtime();

	grid.universe = MPI_COMM_WORLD;

	/// - Reveal information of myself and size of MPI_COMM_WORLD
	CHECK_MPI_ERROR(MPI_Comm_rank(grid.universe, &grid.universal_rank));
	CHECK_MPI_ERROR(MPI_Comm_size(grid.universe, &grid.numtasks));

	char filename[50];
	int error;

	/// Time variables
	double tfinal = 1e-2;
	double interval = 1e-2;
	double data_writing_frequency = 10.00;

	// Read command line input
	// t - T_END for the simulation
	// w - A number deciding how frequent the data should be recorded, default is every 10 seconds
	// i - Time interval between two steps. Default is 1e-2
	if (argc > 1) {
		for (int i = 0; i < argc; i++) {
			if (argv[i][0] == '-') {
				if (argv[i][1] == 'f') {
					sprintf(grid.config_file, "%s", argv[i + 1]);
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

	/// - Read domain configuration from input config file (domain_info.txt).
	/// Why does this function need three parameters? They all come from the same struct.
	error = read_config_file(grid.universal_rank, grid.config_file, &grid);

	/// - Make subdomains according to the information read from the config file (domain_info.txt).

	// WARNING: Why does this function need three parameters? They all come from the same struct.
	// WARNING: In this case and many other cases the structs are passed by value. That's why
	// they all return copies of structs overwriting the initial structs.
	grid = configure_subdomains_topology(grid, grid.num_domains, grid.domains);

	// File written every 1 second.
	int file_write_per_unit_time = (int) (data_writing_frequency * int(1 / interval));
	grid.NO_path = 0;
	grid.cGMP_path = 0;
	grid.smc_model = KNBGR;
	grid.ec_model = KNBGR;
	grid.uniform_jplc = 0.3;
	grid.min_jplc = 0.20;
	grid.max_jplc = 2.5;
	grid.gradient = 0.288e3;
	grid.stimulus_onset_time = 10.00;

	/// - Calculate the number of cells per task.
	grid = set_task_parameters(grid);

	if (grid.my_domain.internal_info.domain_type == STRSEG) {
		grid = make_straight_segment_cart_grids(grid);
	} else if (grid.my_domain.internal_info.domain_type == BIF) {
		grid = make_bifucation_cart_grids(grid);
	}

	/// Set prefixes for output files.
	set_file_naming_strings(&grid);

	/// Initialise checkpoint routine which opens files.
	checkpoint_handle *check = initialise_checkpoint(grid);

	/// Initialising IO_domain for creating writers.
	IO_domain_info* my_IO_domain_info = make_io_domains(&grid);

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

	// TODO: Write a flippin' macro for allocating 2D arrays.
	smc = (SMC_cell**) checked_malloc((grid.num_smc_circumferentially + grid.num_ghost_cells) * sizeof(SMC_cell*), SRC_LOC);
	for (int i = 0; i < (grid.num_smc_circumferentially + grid.num_ghost_cells); i++)
	{
		smc[i] = (SMC_cell*) checked_malloc((grid.num_smc_axially + grid.num_ghost_cells) * sizeof(SMC_cell), SRC_LOC);
	}

	// TODO: Write a flippin' macro for allocating 2D arrays.
	ec = (EC_cell**) checked_malloc((grid.num_ec_circumferentially + grid.num_ghost_cells) * sizeof(EC_cell*), SRC_LOC);
	for (int i = 0; i < (grid.num_ec_circumferentially + grid.num_ghost_cells); i++)
	{
		ec[i] = (EC_cell*) checked_malloc((grid.num_ec_axially + grid.num_ghost_cells) * sizeof(EC_cell), SRC_LOC);
	}

	/// Memory allocation for state vector, the single cell evaluation placeholders (the RHS of the ODEs for each cell) and coupling fluxes is implemented in this section.
	/// In ghost cells, only the state vector array for each type of cells exists including all other cells.
	/// The memory is allocated for all the cells except the ghost cells, hence the ranges 1 to grid.num_ec/smc_circumferentially (inclusive).

	/// SMC domain.
	for (int i = 0; i < (grid.num_smc_circumferentially + grid.num_ghost_cells); i++)
	{
		for (int j = 0; j < (grid.num_smc_axially + grid.num_ghost_cells); j++)
		{
			smc[i][j].vars = (double*) checked_malloc(grid.neq_smc * sizeof(double), SRC_LOC);
			smc[i][j].fluxes = (double*) checked_malloc(grid.num_fluxes_smc * sizeof(double), SRC_LOC);
			smc[i][j].homo_fluxes = (double*) checked_malloc(grid.num_coupling_species_smc * sizeof(double), SRC_LOC);
			smc[i][j].hetero_fluxes = (double*) checked_malloc(grid.num_coupling_species_smc * sizeof(double), SRC_LOC);
		}
	}

	/// EC domain.
	for (int i = 0; i < (grid.num_ec_circumferentially + grid.num_ghost_cells); i++)
	{
		for (int j = 0; j < (grid.num_ec_axially + grid.num_ghost_cells); j++)
		{
			ec[i][j].vars = (double*) checked_malloc(grid.neq_ec * sizeof(double), "Allocation of array for state variables failed.");
			ec[i][j].fluxes = (double*) checked_malloc(grid.num_fluxes_ec * sizeof(double), "Matrix A in ec.");
			ec[i][j].homo_fluxes = (double*) checked_malloc(grid.num_coupling_species_ec * sizeof(double), "Matrix B in ec.");
			ec[i][j].hetero_fluxes = (double*) checked_malloc(grid.num_coupling_species_ec * sizeof(double), "Matrix C in ec.");
		}
	}

	/// Allocating memory for coupling data to be sent and received by MPI communication routines.

	/// sendbuf and recvbuf are 2D arrays having up, down, left and right directions as their first dimension.

	sendbuf = (double**) checked_malloc(4 * sizeof(double*), "sendbuf dimension 1");
	recvbuf = (double**) checked_malloc(4 * sizeof(double*), "recvbuf dimension 1");

	/// Each processor now allocates the memory for send and recv buffers those will hold the coupling information.
	/// Since sendbuf must contain the information of number of SMCs and ECs being sent in the directions,
	/// the first two elements contain the total count of SMCs located on the rank in the relevant dimension (circumferential or axial) and the count of SMCs for which
	/// information is being sent, respectively.
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

	sendbuf[UP] = (double*) checked_malloc(grid.num_elements_send_up * sizeof(double), "sendbuf[UP] dimension 2");

	sendbuf[UP][0] = 1.0; //Start of the 1st segment of SMC array in  in UP direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[UP][1] = (double) extent_s; //End of the 1st segment of SMC array in UP direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[UP][2] = 1.0; //Start of the 1st segment of EC array in UP direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[UP][3] = (double) (extent_e); ///End of the 1st segment of EC array in UP direction (circumferential direction) to be sent to neighbouring processor

	/// Data to send to the neighbour in  direction.
	/// Recording the number of elements in Send buffer in DOWN direction for use in update routine as count of elements.
	grid.num_elements_send_down = grid.added_info_in_send_buf + (grid.num_coupling_species_smc * extent_s + grid.num_coupling_species_ec * extent_e);

	sendbuf[DOWN] = (double*) checked_malloc(grid.num_elements_send_down * sizeof(double), "sendbuf[DOWN] dimension 2");

	sendbuf[DOWN][0] = 1.0; //Start of the 1st segment of SMC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[DOWN][1] = (double) extent_s; //End of the 1st segment of SMC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[DOWN][2] = 1.0; //Start of the 1st segment of EC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[DOWN][3] = (double) extent_e; ///End of the 1st segment of EC array in DOWN direction (circumferential direction) to be sent to neighbouring processor

	extent_s = grid.num_smc_axially;
	extent_e = grid.num_ec_axially;

	/// Data to send to the neighbour in LEFT direction.
	/// Recording the number of elements in Send buffer in Left direction for use in update routine as count of elements.
	grid.num_elements_send_left = grid.added_info_in_send_buf + (grid.num_coupling_species_smc * extent_s + grid.num_coupling_species_ec * extent_e);
	sendbuf[LEFT] = (double*) checked_malloc(grid.num_elements_send_left * sizeof(double), "sendbuf[LEFT] dimension 2");

	sendbuf[LEFT][0] = 1.0; //Start of the 1st segment of SMC array in  in LEFT direction (axial direction) to be sent to neighbouring processor
	sendbuf[LEFT][1] = (double) extent_s; //END of the 1st segment of SMC array in  in LEFT direction (axial direction) to be sent to neighbouring processor
	sendbuf[LEFT][2] = 1.0; //Start of the 1st segment of EC array in  in LEFT direction (axial direction) to be sent to neighbouring processor
	sendbuf[LEFT][3] = (double) extent_e; //END of the 1st segment of EC array in  in LEFT direction (axial direction) to be sent to neighbouring processor

	/// Data to send to the neighbour in RIGHT direction.
	/// Recording the number of elements in Send buffer in RIGHT direction for use in update routine as count of elements.
	grid.num_elements_send_right = grid.added_info_in_send_buf + (grid.num_coupling_species_smc * extent_s + grid.num_coupling_species_ec * extent_e);

	sendbuf[RIGHT] = (double*) checked_malloc(grid.num_elements_send_right * sizeof(double), "sendbuf[RIGHT] dimension 2");

	sendbuf[RIGHT][0] = 1.0; //Start of the 1st segment of SMC array in  in RIGHT direction (axial direction) to be sent to neighbouring processor
	sendbuf[RIGHT][1] = (double) extent_s; //END of the 1st segment of SMC array in  in RIGHT direction (axial direction) to be sent to neighbouring processor
	sendbuf[RIGHT][2] = 1.0; //Start of the 1st segment of EC array in  in RIGHT direction (axial direction) to be sent to neighbouring processor
	sendbuf[RIGHT][3] = (double) extent_e; //END of the 1st segment of EC array in  in RIGHT direction (axial direction) to be sent to neighbouring processor

	/// Call communication to the number of elements to be recieved by neighbours and allocate memory of recvbuf for each direction accordingly.
	grid = communicate_num_recv_elements_to_nbrs(grid);
	/// memory allocation

	/// data to receive from the neighbour in UP direction
	recvbuf[UP] = (double*) checked_malloc(grid.num_elements_recv_up * sizeof(double), "recvbuf[UP] dimension 2");

	/// data to recv from the neighbour in DOWN direction
	recvbuf[DOWN] = (double*) checked_malloc(grid.num_elements_recv_down * sizeof(double), "recvbuf[DOWN] dimension 2");

	/// data to receive from the neighbour in LEFT direction
	recvbuf[LEFT] = (double*) checked_malloc(grid.num_elements_recv_left * sizeof(double), "recvbuf[LEFT] dimension 2");

	/// data to receive from the neighbour in RIGHT direction
	recvbuf[RIGHT] = (double*) checked_malloc(grid.num_elements_recv_right * sizeof(double), "recvbuf[RIGHT] dimension 2");

	grid.NEQ = grid.neq_smc * (grid.num_smc_axially * grid.num_smc_circumferentially) + grid.neq_ec * (grid.num_ec_axially * grid.num_ec_circumferentially);

	/// Setup output streams to write data in files. Each node opens an independent set of files and write various state variables into it.
	///	checkpoint_handle *check = initialise_checkpoint(myRank);

	/// Setting up the solver.
	double tnow = 0.0;

	// Error control variables.
	double TOL = 1e-6, absTOL = 1e-7;
	double *thres = (double*)checked_malloc(grid.NEQ * sizeof(double), "Threshold array for RKSUITE");
	for (int i = 0; i < grid.NEQ; i++)
	{
		thres[i] = absTOL;
	}

	// State variables.
	double* y = (double*) checked_malloc(grid.NEQ * sizeof(double), "Solver array y for RKSUITE");

	// Internal work space.
	double* yp = (double*) checked_malloc(grid.NEQ * sizeof(double), "Workspace array y for RKSUITE");

	/// Initialise state variables and coupling data values.

	// checkpoint(check, grid, &tnow, y, smc, ec);
	Initialize_koeingsberger_smc(grid, y, smc);
	Initialize_koeingsberger_ec(grid, y, ec);

	// Reverse mapping from state vector to cells.
	// Essential for restarts, when data is loaded from a checkpoint.
	int err = map_solver_output_to_cells(grid, y, smc, ec);

	if (err != 0) {
		printf("[%d] error in mapping y to cells\n", grid.universal_rank);
		MPI_Abort(grid.universe, 1000);
	}

	// Initialising the coupling coefficients to be used in the ODEs.
	int state = couplingParms(CASE, &cpl_cef);

	// Debug output.
	dump_rank_info(check, cpl_cef, grid, my_IO_domain_info);

	// Debug/validation/reporting.
	Total_cells_in_computational_domain(grid);

	// Reading all points coordinates.
	// int ret = read_topology_info((char *)"files/configuration_info.txt", &grid, smc, ec);

	// This is read in here for validation purposes in the output.
	// the solver will reset JPLC and read later it when the time is right.
	read_init_ATP(&grid, ec);

#ifdef RK_SUITE
	rksuite_solver_CT(tnow, tfinal, interval, y, yp, grid.NEQ, TOL, thres, file_write_per_unit_time, check, grid.solution_dir, my_IO_domain_info);
#elif defined ARK_ODE
	arkode_solver(tnow, tfinal, interval, y, grid.NEQ, TOL, absTOL, file_write_per_unit_time, check, grid.solution_dir, my_IO_domain_info);
#elif defined BOOST_ODEINT
	odeint_solver(tnow, tfinal, interval, y, grid.NEQ, TOL, absTOL, file_write_per_unit_time, check, grid.solution_dir, my_IO_domain_info);
#else
#error ODE solver not selected. Use -DRK_SUITE | -DARK_ODE | -DBOOST_ODEINT during compilation.
#endif

	// Final_checkpoint(check, grid);
	update_elapsed_time(check, grid, &elps_t, my_IO_domain_info);

	// fclose(grid.logptr);

	MPI_Finalize();
	return (0);
}
