#ifndef _COMPUTE_LIB_
#define _COMPUTE_LIB_

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include "macros.h"

#define success 0
#define	none -1

#define local 0
#define remote 1

/**
 * When converting double to char for writing via MPI-IO to write in ASCII format
 * the double is to be truncated to 12 characters including the decimal point.
 * The 13th char is a white space (new line or tab).
 */
#define NUM_DBL_TO_CHAR_BYTES 64

#define WRITER_COLOR 43
#define WRITER_KEY 43
#define COMPUTE_ONLY_COLOR 1
#define COMPUTE_ONLY_KEY 1

#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

#define STRSEG 0 ///< A straight segment.
#define BIF 1 ///< A bifurcation.
#define P 1 ///< Parent.
#define L 2 ///< Left branch.
#define R 3 ///< Right branch.

/**
 * Macros for use in retrieving mesh topological data from grid_parms.info.
 * @{ */
#define 	TOTAL_POINTS		0
#define		POINTS_m			1
#define     POINTS_n			2
#define		TOTAL_CELLS			3
#define		CELLS_m				4
#define		CELLS_n				5
/** @} */

/**
 * Macros representing mesh types.
 * @{ */
#define     PROCESS_MESH 		0
#define 	SMC_MESH 			1
#define		EC_MESH 			2
#define		EC_CENT_MESH 		3

#define 	ProcessCell			4
#define 	smcCell				5
#define		ecCell				6
#define		ecCentroidCell		7

#define 	ProcessCellType		8
#define		smcCellType			9
#define		ecCellType			10
#define		ecCentroidCellType	11
#define 	ec_ATP_Conc			12
#define 	ec_WSS_val			13
/** @} */

#define		smcDataLength		0
#define		ecDataLength		1

/**
 * Helper functions for exponentiation to integer powers.
 * @{ */
#define P2(x) ((x)*(x))
#define P3(x) ((x)*(x)*(x))
#define P4(x) ((x)*(x)*(x)*(x))
/** @} */

/**
 * Conductance / coupling coefficients.
 */
struct conductance {
	double Vm_hm_smc,	///< Homocellular membrane potential coupling between SMCs.
			Vm_hm_ec,	///< Homocellular membrane potential coupling between ECs.
			Ca_hm_smc,	///< Homocellular Ca coupling between SMCs.
			Ca_hm_ec,	///< Homocellular Ca coupling between ECs.
			IP3_hm_smc,	///< Homocellular IP3 coupling between SMCs.
			IP3_hm_ec,	///< Homocellular IP3 coupling between ECs.
			Vm_ht_smc,	///< Heterocellular membrane potential coupling between SMCs.
			Vm_ht_ec,	///< Heterocellular membrane potential coupling between ECs.
			Ca_ht_smc,	///< Heterocellular Ca coupling between SMCs.
			Ca_ht_ec,	///< Heterocellular Ca coupling between ECs.
			IP3_ht_smc,	///< Heterocellular IP3 coupling between SMCs.
			IP3_ht_ec;	///< Heterocellular IP3 coupling between ECs.
};

#if 1
// TODO: Pretty sure this can be chucked out.
struct node {
	int domain_type, ///< Bifurcation or a straight segment.
	    domain_index,
	    domain_start, ///< Universal ranks from MPI_COMM_WORLD.
	    domain_end, ///< Universal ranks from MPI_COMM_WORLD.
		parent_branch_case_bifurcation,	///< If my parent is a bifurcation which branch am I a child of?
		m, n; ///< Row and columns in my MPI_sub_world.
	char boundary_tag; ///< An identifier showing whether I am a rank from top or bottom edge of a subdomain.
	int half_marker; ///< A marker for demarcating the bottom edge of the Left/Right daughter artery.
					///exists which couples not to the parent artery but the other daughter segment. This can have following values:
					/// 1. half coupling to parent
					/// 2. half coupling to other daughter
					/// 3. half splitting in the middle with left portion coupling to parent, and right portion of data coupling
					///    to other daughter segment.

	//double d, l; ///< Diameter and length scales.

};
// TODO: Pretty sure this can be chucked out.
struct my_tree {
	node internal_info;
	node left_child, right_child, parent;
	//double z_offset_start, z_offset_end;				/// These are domain offsets start and end points to demacated distance in
	//double local_z_start, local_z_end;					/// z direction spanned by a processor's own sub-domain that it belongs to.
};
#endif

#if 0
struct glb_domn_inf {
	int num_subdomains,							///number of total subdomains
			*m, *n,	///number of grid points axially, number of grid points circumferentially
			*list_type_subdomains,	///list of types of subdomains, either STRSEG=0 or BIF=1
			*list_num_ec_axially_per_domain;	///list of number of ECs axially in each subdomain in sequence of increasing z_coordinate of distance.
	double **list_domain_z_coord_index;	///stores the start and end of each subdomain in axial direction. This will be used to estimate the agonist
};
#endif

///on that particular z coordinate. First two elements store the coords for any STRSEG or Left/Right child of
///bifurcation where as the last two elements store coords for the parent segment of a bifurcation, if the domain
///type is BIF

#define DOMAIN_NUM 0
#define DOMAIN_TYPE 1
#define AX_QUADS 2
#define CR_QUADS 3
#define PARENT_DOMAIN_NUM 4
#define LEFT_DOMAIN_NUM 5
#define RIGHT_DOMAIN_NUM 6
#define AX_ECS 7
#define CR_SMCS 8
#define NUM_CONFIG_ELEMENTS 9

typedef struct {
	// double tfinal;

	///General information on cell geometry and the geometric primitive constructed.
	// double hx_smc, hx_ec, hy_smc, hy_ec, requested_length, requested_diameter, corrected_length, corrected_diameter, new_circ;

	int num_smc_fundblk_circumferentially;
	int num_ec_fundblk_circumferentially;
	int num_smc_fundblk_axially;
	int num_ec_fundblk_axially;
	int num_ghost_cells;
	int num_fluxes_smc; // Number of SMC ioinic currents to be evaluated for eval of LHS of the d/dt terms of the ODEs.
	int num_fluxes_ec; // Number of EC ioinic currents to be evaluated for eval of LHS of the d/dt terms of the ODEs.
	int num_coupling_species_smc; // Number of SMC coupling species homogenic/heterogenic.
	int num_coupling_species_ec; // Number of SMC coupling species homogenic/heterogenic.
	int neq_smc; // Number of SMC ODEs for a single cell.
	int neq_ec; // Number of EC ODEs for a single cell.

	int
	/// Global domain information storage with these elements:
	// [DOMAIN_NUM, DOMIN_TYPE, AX_QUADS, CR_QUADS, PARENT_DOMAIN_NUM, LEFT_DOMAIN_NUM, RIGHT_DOMAIN_NUM, AX_ECS, CR_SMCS].
	num_domains, **domain_params,
	///total grid points axially
	m,
	///total grid points circumferentially
	n,
	///placeholder to retrieve information from topology files
	// **info,
	/// My coordinates
	coords[2],
	///Coordinates for neighbour tasks
	nbrs[2][4],
	///Node payload information(number of cells laid out on a node)
	num_ec_axially, num_ec_circumferentially, num_smc_axially, num_smc_circumferentially, neq_ec_axially, neq_smc_axially,
	///Total number of state variables in the computational domain
	NEQ,
	///Number of elements added to the Send buffer for sending relevant information on the content of the buffer to receiving task
	added_info_in_send_buf,
	///This is global and local MPI information
	num_ranks, universal_rank, /// numtasks = total CPUs in MPI_COMM_WORLD,
	num_ranks_branch, rank_branch, /// tasks = total CPUs in my-subdomain's comm (branch)
	color, key,

	//Each processor on the edges of each branch contains brach_tag can have one of four values P=parent = 1, L=Left branch = 2, R=Right branch = 3.
	//If branch_tag=0, this implies that the rank is located interior or doesn't  contain a remote neighbour on any other branch.
	branch_tag,
	/// Variables for remote MPI information (P=parent, L & R = Left & Right branch respectively).
	scheme, offset_P, offset_L, offset_R, flip_array[4],
	/// Number of elements being sent and received.
	num_elements_send_up, num_elements_send_down, num_elements_send_left, num_elements_send_right,
	num_elements_recv_up, num_elements_recv_down, num_elements_recv_left, num_elements_recv_right;

	// double **coordinates;

	///Information for spatial variation in agonist
	double min_jplc, max_jplc, gradient, uniform_jplc, stimulus_onset_time;	/// the time when spatially varying agonist kicks in

	my_tree my_domain;

	// glb_domn_inf global_domain_info;

	// Allow three types of communicators to exist, first resulting from subdomain allocation, second resulting from comm_split
	// operation on MPI_COMM_WORLD and the other a Cartesian communicator arising from Cart_create operation.
	MPI_Comm universe, cart_comm;

	//MPI_Comm sub_universe;
	//int sub_universe_numtasks, sub_universe_rank;
	//int my_domain_color, my_domain_key;

	int smc_model, ec_model;	// These are placeholders for the selection of model to be simulated in each cell.
	int NO_path, cGMP_path;	// Specific for Tsoukias model to signal whether to activate NO and cGMP pathways for vasodilation.

	char suffix[16];	// this is for use in the naming convention of the IO files to recognise and record
						// which files are associated with a given task/processor.

						///Temporary array for use in time profiling checkpointing
	// double **time_profile;
	// FILE* logptr;
	// int num_parameters;			///Number of parameters e.g. JPLC, ATP, WSS etc those are to be used to stimulate the discrete cell models.

	int logfile_displacements;
	char *logfile_write_buffer;
	char solution_dir[1024], time_profiling_dir[1024], config_file[1024];
} grid_parms;

#if 0
///Structure to store coupling data received from the neighbouring task.
typedef struct {
	double c, v, I;
} nbrs_data;
#endif

typedef struct {
	double *vars;		///storage for the state variables corresponding to an SMC.
	double NO, NE, I_stim;		///specific to Tsoukias model
	int node_row, node_col;	///stores coordinates of the node on which I am located.
	int my_row, my_col;		///stores my location on the node.
	double* fluxes;			    ///stores single cell fluxes
	double* homo_fluxes;			    ///stores homogeneous coupling fluxes
	double* hetero_fluxes;			    ///stores heterogeneous coupling fluxes
	// double x_coordinate[4], y_coordinate[4], z_coordinate[4];
	int cell_index[4];
	conductance cpl_cef;
} SMC_cell;

typedef struct {
	double *vars;		///storage for the state variables corresponding to an SMC.
	int node_row, node_col;	///stores coordinates of the node on which I am located.
	int my_row, my_col;		///stores my location on the node.
	double* fluxes;			    ///stores single cell fluxes
	double* homo_fluxes;			    ///stores homogeneous coupling fluxes
	double* hetero_fluxes;			    ///stores heterogeneous coupling fluxes
	// double z_coord;
	// double x_coordinate[4], y_coordinate[4], z_coordinate[4];
	int cell_indx[4];
	// double centeroid_point[3];
	int centeroid_cell;
	double JPLC;			    ///local agonsit concentration  on my GPCR receptor (an ith EC)
	conductance cpl_cef;
} EC_cell;


#if 1
// TODO: This is to be thrown away along with all related functions. They do nothing.
typedef struct {
	MPI_File
	/* common handlers */
	logptr,
	jplc,
	coords,
	elapsed_time,
	itter_count,
	//Task topology file
	task_mesh,
	//SMC Data file
	smc_data_file,
	//EC Data file
	ec_data_file;
} checkpoint_handle;
#endif

typedef struct {
	double

	/// Communication timings.
	async_comm_calls_t1, async_comm_calls_t2,
	async_comm_calls_wait_t1, async_comm_calls_wait_t2,
	remote_async_comm_calls_t1, remote_async_comm_calls_t2,
	remote_async_comm_calls_wait_t1, remote_async_comm_calls_wait_t2,
	update_sendbuf_t1, update_sendbuf_t2,
	update_recvbuf_t1, update_recvbuf_t2,
	barrier_in_solver_before_comm_t1, barrier_in_solver_before_comm_t2,
	total_comms_cost_t1, total_comms_cost_t2,
	diff_update_sendbuf,
	diff_update_recvbuf,
	diff_async_comm_calls,
	diff_async_comm_calls_wait,
	diff_remote_async_comm_calls,
	diff_remote_async_comm_calls_wait,
	diff_barrier_in_solver_before_comm,
	diff_total_comms_cost,

	///Solver timings.
	map_function_t1, map_function_t2,
	single_cell_fluxes_t1, single_cell_fluxes_t2,
	coupling_fluxes_t1, coupling_fluxes_t2,
	solver_t1, solver_t2,
	write_t1, write_t2,
	diff_map_function,
	diff_single_cell_fluxes,
	diff_coupling_fluxes,
	diff_solver,
	diff_write;
	int computeDerivatives_call_counter;
	double aggregate_write, aggregate_compute, aggregate_comm;
	double max_compute, max_comm, max_write, min_compute, min_comm, min_write;
	int max_compute_index, max_comm_index, max_write_index, min_compute_index, min_comm_index, min_write_index;
} time_stamps;

typedef struct {
	double t_new, t_old, elapsed_time;
} time_keeper;

#if 0
typedef struct {
	double** points;
	int** cells;
} vtk_info;
#endif

typedef struct {
	///IO_domain related members
	int IO_rank, IO_domain_ID, num_IO_tasks, num_IO_domains, writer_rank, writer_tasks;
	int my_IO_domain_color, my_IO_domain_key;
	int *my_IO_domain_members;
	char ***data_filenames;				/// root rank stores receives data into this array
	char **my_data_filenames;			/// filenames of files each Rank is suppose to have its data in.
	int *my_IO_domain_member_disp;
	MPI_Comm IO_COMM, writer_comm;
} IO_domain_info;

#if 0
typedef struct {
	char *process_mesh_points, *smc_mesh_points, *ec_mesh_points, *ec_centroid_points;
	char *process_mesh_cells, *smc_mesh_cells, *ec_mesh_cells, *ec_centroid_cells;
	char *process_mesh_type, *smc_mesh_type, *ec_mesh_type, *ec_centroid_type;
	char *ci, *si, *vi, *wi, *Ii, *cj, *sj, *vj, *Ij, *cpCi, *cpVi, *cpIi, *cpCj, *cpVj, *cpIj;
	char *jplc,*atp,*wss;
	int	 *buffer_length,*smc_stat_var_buffer_length,*ec_stat_var_buffer_length, *smc_cpl, *ec_cpl,jplc_buffer_length,atp_buffer_length,wss_buffer_length;
} data_buffer;
#endif

int couplingParms(int CASE, conductance* cpl_cef);
void Initialize_koeingsberger_smc(grid_parms, double*, SMC_cell**);
void Initialize_koeingsberger_ec(grid_parms, double*, EC_cell**);
int map_solver_output_to_cells(grid_parms, double*, SMC_cell**, EC_cell**);

grid_parms communicate_num_recv_elements_to_nbrs(grid_parms);
void communication_update_sendbuf(grid_parms, double**, SMC_cell**, EC_cell**);
void communication_update_recvbuf(grid_parms, double**, SMC_cell**, EC_cell**);
void determine_source_destination(grid_parms, int*, int*);
void communication_async_send_recv(grid_parms, double**, double**, SMC_cell**, EC_cell**);

//Cell dynamics evaluation handlers. These contain the ODEs for representative models from different sources.
void single_cell(double, double*, grid_parms, SMC_cell**, EC_cell**);
void coupling(double, double*, grid_parms, SMC_cell**, EC_cell**, conductance);
void tsoukias_smc(grid_parms, SMC_cell**);
void koenigsberger_smc(grid_parms, SMC_cell**);
void tsoukias_smc_derivatives(double*, grid_parms, SMC_cell**);
void koenigsberger_smc_derivatives(double*, grid_parms, SMC_cell**);
void koenigsberger_ec(grid_parms, EC_cell**);
void koenigsberger_ec_derivatives(double, double*, grid_parms, EC_cell**);

///Checkpoint functions.
checkpoint_handle* initialise_checkpoint(grid_parms);
checkpoint_handle* initialise_time_wise_checkpoint(checkpoint_handle*, grid_parms, int, char*, IO_domain_info*);
int checkpoint(checkpoint_handle*, grid_parms, double*, double*, SMC_cell**, EC_cell**);

void open_common_checkpoint(checkpoint_handle*, grid_parms);
void open_tsoukias_smc_checkpoint(checkpoint_handle*, grid_parms, char*);
void open_koenigsberger_smc_checkpoint(checkpoint_handle*, grid_parms, int, char*, IO_domain_info*);
void open_tsoukias_ec_checkpoint(checkpoint_handle*, grid_parms, char*);
void open_koenigsberger_ec_checkpoint(checkpoint_handle*, grid_parms, int, char*, IO_domain_info*);
void open_coupling_data_checkpoint(checkpoint_handle*, grid_parms, int, char*, IO_domain_info*);

//void dump_smc_async(grid_parms, SMC_cell**, checkpoint_handle*, int);
//void dump_ec_async(grid_parms, EC_cell**, checkpoint_handle*, int);
//void write_smc_and_ec_data(checkpoint_handle*, grid_parms*, double, SMC_cell**, EC_cell**, int, IO_domain_info*,data_buffer*);
void final_checkpoint(checkpoint_handle*, grid_parms);
void close_common_checkpoints(checkpoint_handle*);
void close_time_wise_checkpoints(checkpoint_handle*);
void close_time_profiling_checkpoints(checkpoint_handle*);

void dump_rank_info(checkpoint_handle*, conductance, grid_parms, IO_domain_info*);
void dump_smc_with_ghost_cells(grid_parms, SMC_cell**, checkpoint_handle*, int);
void dump_ec_with_ghost_cells(grid_parms, EC_cell**, checkpoint_handle*, int);
void checkpoint_with_ghost_cells(checkpoint_handle*, grid_parms, double, SMC_cell**, EC_cell**, int);

// Solver wrapper functions.
#ifdef RK_SUITE
void rksuite_solver_CT(double, double, double, double*, double*, int, double, double*, int, checkpoint_handle*, char*, IO_domain_info*);
#endif
#ifdef ARK_ODE
void arkode_solver(double, double, double, double*, int, double, double, int, checkpoint_handle*, char*, IO_domain_info*);
#endif
#ifdef BOOST_ODEINT
void odeint_solver(double, double, double, double*, int, double, double, int, checkpoint_handle*, char*, IO_domain_info*);
#endif

int compute_with_time_profiling(time_stamps*, grid_parms, SMC_cell**, EC_cell**, conductance cpl_cef, double, double*, double*);
int compute(grid_parms, SMC_cell**, EC_cell**, conductance cpl_cef, double, double*, double*);

///These are debugging functions, not used in production runs.
void print_domains(FILE*, grid_parms, SMC_cell**, EC_cell**);
void print_send_buffer(FILE*, grid_parms, double**);
void print_recv_buffer(FILE*, grid_parms, double**);
void print_compare(double, double*, grid_parms, SMC_cell**, EC_cell**);

//Topology related functions
grid_parms make_bifucation_cart_grids(grid_parms);
grid_parms make_straight_segment_cart_grids(grid_parms);
void set_task_parameters(grid_parms *);
void configure_subdomains_topology(grid_parms *);

void checkpoint_timing_data(grid_parms, checkpoint_handle*, double, time_stamps, int, int);
double agonist_profile(double, grid_parms, int, int, double);
void initialize_t_stamp(time_stamps*);

void Initialize_tsoukias_smc(grid_parms grid, double y[], SMC_cell** smc);
void read_config_file(int, char*, grid_parms*);
void set_file_naming_strings(grid_parms* grid);
void write_elapsed_time(checkpoint_handle*, grid_parms, time_keeper*); //,IO_domain_info*);
int determine_file_offset_for_timing_data(checkpoint_handle* check, grid_parms grid);

void Total_cells_in_computational_domain(grid_parms gird);

void Record_timing_data_in_arrays(grid_parms, double, time_stamps, int, double**);
void process_time_profiling_data(grid_parms, double**, int);
void minimum(double* table, int size, double *value, int *index);
void maximum(double* table, int size, double *value, int *index);
void average(double* table, int size, double *value);

void rksuite_solver_CT_debug(double tnow, double tfinal, double interval, double *y, double* yp, int total, double TOL, double* thres,
		int file_write_per_unit_time, int line_number, checkpoint_handle *check);

int read_topology_info(char*, grid_parms*, SMC_cell**, EC_cell**);
void read_init_ATP(grid_parms *grid, EC_cell **ECs);
//void read_coordinates(int** info, vtk_info* mesh, int branch, int mesh_type, int points, int cells, int read_counts[2]);
IO_domain_info* make_io_domains(grid_parms* grid);

//void gather_tasks_mesh_point_data_on_writers(grid_parms*, IO_domain_info*, data_buffer*, SMC_cell**, EC_cell**);
//void gather_smc_mesh_data_on_writers(grid_parms*, IO_domain_info*, data_buffer*, SMC_cell**);
//void gather_ec_mesh_data_on_writers(grid_parms*, IO_domain_info*, data_buffer*, EC_cell**);

//void gather_smcData(grid_parms* , IO_domain_info* , data_buffer* , SMC_cell**, int );
//void gather_ecData(grid_parms*, IO_domain_info*, data_buffer*, EC_cell**, int);
//void gather_JPLC_map(grid_parms*, IO_domain_info*, data_buffer*, EC_cell**);

//void write_process_mesh(checkpoint_handle*, grid_parms* , IO_domain_info* , data_buffer*, char*);
//void write_JPLC_map(checkpoint_handle*, grid_parms*, IO_domain_info*, data_buffer*, EC_cell**,char* path);
//void dump_smc_data(checkpoint_handle*, grid_parms* , IO_domain_info* , data_buffer* , SMC_cell**, int);
//void dump_ec_data(checkpoint_handle*, grid_parms*, IO_domain_info*, data_buffer*, EC_cell**,int);

void memory_diagnostics(FILE*);

void write_timing(char* file_prefix, grid_parms grid, double field, IO_domain_info* my_IO_domain_info);
void checkpoint_coarse_time_profiling_data(grid_parms grid, time_stamps* t_stamp, IO_domain_info* my_IO_domain_info);
void write_min_max_timing(char* file_prefix, grid_parms grid, double field, IO_domain_info* my_IO_domain_info);

/**
 * \brief Catch failed memory allocation.
 */
void* checked_malloc(size_t bytes, const char* errmsg);

#endif
