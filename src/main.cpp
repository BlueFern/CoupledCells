#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"

#ifdef CVODE
#include <nvector/nvector_serial.h>
#endif

using namespace std;

    conductance		cpl_cef;
    celltype1** 	smc;
    celltype2** 	ec;
    double		**sendbuf,**recvbuf;
    grid_parms		grid;

int 	CASE=1;
///***************************************************************************************/
///************ checked_malloc(size_t bytes, FILE* errfile, const char* errmsg)*************/
///***************************************************************************************/
void* checked_malloc(size_t bytes, const char* errmsg){
	void *pval = malloc(bytes);

	if (pval == NULL) {
		fprintf(stdout, "%s",errmsg);
		MPI_Abort(MPI_COMM_WORLD, 100);
	}

	return pval;
}

int main(int argc, char* argv[]) {




    ///Global declaration of request and status update place holders.
    ///Request and Status handles for nonblocking Send and Receive operations, for communicating with each of the four neighbours.
    MPI_Request reqs[8];
    MPI_Status stats[8];

    ///Initialize MPI
    MPI_Init(&argc, &argv);

    grid.universe = MPI_COMM_WORLD;


//Reveal information of myself and size of MPI_COMM_WORLD
    check_flag(MPI_Comm_rank(grid.universe, &grid.universal_rank),
	    "error Comm_rank");
    check_flag(MPI_Comm_size(grid.universe, &grid.numtasks),
	    "error Comm_size");

//Test case
//config
    int num_subdomains	=	1;
    int
         m[num_subdomains],
         n[num_subdomains],
         e = 4,	//ECs per processor in axial direction
         s = 4;	//SMCs per processor in circumferential direction

         m[0]	=	512;
        /* m[1]	=	2;
         m[2]	=	4;
         m[3]	=	32;
         m[4]	=	32;
         m[5]	=	4;
         m[6]	=	32;
         m[7]	=	32;*/

    	for (int i = 0; i < num_subdomains; i++) {
    		n[i] = 4;
    	}

   	int **domains;
    domains = (int**) checked_malloc(num_subdomains * sizeof(int*),"Subdomain information allocation");

    /*second coordinate has following information:
     * Element 0: 	Key_val or serial number of the subdomain
     * Element 1:	Subdomain Type (2 possibilities and their values)
     * 				1. Straight Segment (STRSEG)					(0)
     * 				2. Bifurcation	(BIF)							(1)
     * Element 2:	Axial extent of processor of current key_val
     * Element 3: 	circumferential extent of processors of current key_val
     * Element 4:	Parent subdomain key_val of current Key_val.
     * Element 5: 	Left Child subdomain key_val of the current Key_val.
     * Element 6:   Right Child subdomain key_val of the current Key_val.
     *
     * In the case of elements 6 & 7, if subdomain type of current key_val is a straight segment, left Child is positive or zero, and right Child is negative.
     * If subdomain type of current key_val is a bifurcation, then both right and left child subdomains are non-negative.
     */
    for (int i = 0; i < num_subdomains; i++)
	{
	domains[i] = (int*) checked_malloc(7 * sizeof(int),
		"Subdomains array elements allocation");
	}
    domains[0][0] 	= 	0;
    domains[0][1] 	= 	STRSEG;
    domains[0][2] 	= 	m[0];
    domains[0][3] 	= 	n[0];
    domains[0][4] 	= 	none;
    domains[0][5] 	= 	none;
    domains[0][6] 	= 	none;
   /*
    domains[1][0]	=	1;
    domains[1][1]	=	STRSEG;
    domains[1][2]	=	m[1];
    domains[1][3]	=	n[1];
    domains[1][4]	=	0;
    domains[1][5]	=	5;
    domains[1][6]	=	none;
    
    domains[2][0]	=	2;
    domains[2][1]	=	BIF;
    domains[2][2]	=	m[2];
    domains[2][3]	=	n[2];
    domains[2][4]	=	0;
    domains[2][5]	=	4;
    domains[2][6]	=	3;
    
    domains[3][0]	=	3;
    domains[3][1]	=	STRSEG;
    domains[3][2]	=	m[3];
    domains[3][3]	=	n[3];
    domains[3][4]	=	2;
    domains[3][5]	=	none;
    domains[3][6]	=	none;
    
    domains[4][0]	=	4;
    domains[4][1]	=	STRSEG;
    domains[4][2]	=	m[4];
    domains[4][3]	=	n[4];
    domains[4][4]	=	2;
    domains[4][5]	=	none;
    domains[4][6]	=	none;
    
    domains[5][0]	=	5;
    domains[5][1]	=	BIF;
    domains[5][2]	=	m[5];
    domains[5][3]	=	n[5];
    domains[5][4]	=	1;
    domains[5][5]	=	7;
    domains[5][6]	=	6;
    
    domains[6][0]	=	6;
    domains[6][1]	=	STRSEG;
    domains[6][2]	=	m[6];
    domains[6][3]	=	n[6];
    domains[6][4]	=	5;
    domains[6][5]	=	none;
    domains[6][6]	=	none;
    
    domains[7][0]	=	7;
    domains[7][1]	=	STRSEG;
    domains[7][2]	=	m[7];
    domains[7][3]	=	n[7];
    domains[7][4]	=	5;
    domains[7][5]	=	none;
    domains[7][6]	=	none;
    */

    grid =  make_subdomains(grid, num_subdomains, domains);

///Time variables
	double tfinal =1.00;
	double interval = 1e-2;
//File written every 1 second
	int file_write_per_unit_time = int(1/interval);

	grid.uniform_jplc = 0.1, grid.min_jplc = 0.3, grid.max_jplc = 1.195, grid.gradient =
			0.5e-2; grid.stimulus_onset_time	=0.10;


grid = set_geometry_parameters(grid,e,s);

	if (grid.my_domain.internal_info.domain_type == STRSEG) {
		grid = make_straight_segment(grid);
	} else if (grid.my_domain.internal_info.domain_type == BIF) {
		grid = make_bifucation(grid);
	}

	///Initialize checkpoint routine which opens files
	    checkpoint_handle *check = initialise_checkpoint(grid);


///Now allocate memory space for the structures represegird.nting the cells and the various members of those structures.


//Each of the two cell grids have two additional rows and two additional columns as ghost cells.
//Follwing is an example of a 5x7 grid with added ghost cells on all four sides. the 0s are the actual
//members of the grid whereas the + are the ghost cells.

// + + + + + + + + +
// + 0 0 0 0 0 0 0 +
// + 0 0 0 0 0 0 0 +
// + 0 0 0 0 0 0 0 +
// + 0 0 0 0 0 0 0 +
// + 0 0 0 0 0 0 0 +
// + + + + + + + + +

	smc 	= (celltype1**) checked_malloc((grid.num_smc_circumferentially+grid.num_ghost_cells)* sizeof(celltype1*), "smc");
		for (int i=0; i<(grid.num_smc_circumferentially+grid.num_ghost_cells); i++){
			smc[i]	= (celltype1*) checked_malloc((grid.num_smc_axially+grid.num_ghost_cells)* sizeof(celltype1), "smc column dimension");

		}
		ec 	= (celltype2**) checked_malloc((grid.num_ec_circumferentially+grid.num_ghost_cells)* sizeof(celltype2*), "ec");
		for (int i=0; i<(grid.num_ec_circumferentially+grid.num_ghost_cells); i++){
			ec[i]	= (celltype2*) checked_malloc((grid.num_ec_axially+grid.num_ghost_cells)* sizeof(celltype2), "ec column dimension");
		}



	///Memory allocation for state vector, the single cell evaluation placeholders (The RHS of the ODEs for each cell) and coupling fluxes is implemented in this section.
	///In ghost cells, only the state vector array for each type of cells exists including all other cells.
	///The memory is allocated for all the cells except the ghost cells, hence the ranges 1 to grid.num_ec_circumferentially(inclusive).
		///SMC domain
		for (int i = 0; i < (grid.num_smc_circumferentially+grid.num_ghost_cells); i++) {
			for (int j = 0; j < (grid.num_smc_axially+grid.num_ghost_cells); j++) {

				smc[i][j].p	=	(double*)checked_malloc(grid.neq_smc*sizeof(double),"allocation of array for state variables failed");
				smc[i][j].A = (double*) checked_malloc(
						grid.num_fluxes_smc * sizeof(double),
						"matrix A in smc");
				smc[i][j].B = (double*) checked_malloc(
						grid.num_coupling_species_smc * sizeof(double),
						"matrix B in smc");
				smc[i][j].C = (double*) checked_malloc(
						grid.num_coupling_species_smc * sizeof(double),
						"matrix C in smc");
			}
		}

		///EC domain
		for (int i = 0; i < (grid.num_ec_circumferentially+grid.num_ghost_cells); i++) {
			for (int j = 0; j < (grid.num_ec_axially+grid.num_ghost_cells); j++) {
				ec[i][j].q	=	(double*)checked_malloc(grid.neq_ec*sizeof(double),"allocation of array for state variables failed");

				ec[i][j].A = (double*) checked_malloc(
						grid.num_fluxes_ec * sizeof(double),
						"matrix A in ec");
				ec[i][j].B = (double*) checked_malloc(
						grid.num_coupling_species_ec * sizeof(double),
						"matrix B in ec");
				ec[i][j].C = (double*) checked_malloc(
						grid.num_coupling_species_ec * sizeof(double),
						"matrix C in ec");
			}
		}




		///Allocating memory space for coupling data to be sent and received by MPI communication routines.

		///sendbuf and recvbuf are 2D arrays having up,down,left and right directions as their first dimension.
		///Each dimension is broken down into two segments, e.g. up1,up2,down1 & down2,etc..
		///The length of the second dimension is equal to half the number of cells for which the information is to be sent and received.
		///Thus each communicating pair will exchange data twice to get the full lenght.

		sendbuf = (double**) checked_malloc(8 * sizeof(double*),
				"sendbuf dimension 1");
		recvbuf = (double**) checked_malloc(8 * sizeof(double*),
				"recvbuf dimension 1");

		///Each processor now allocates the memory for send and recv buffers those will hold the coupling information.
		///Since sendbuf must contain the information of number of SMCs and ECs being sent in the directions,
		///the first two elements contain the total count of SMCs located on the rank in the relevant dimension (circumferential or axial) and the count of SMCs for which
		///information is being sent, respectively.
		///The next two elements contain the same information for ECs.

		int extent_s, extent_e;	///Variables to calculate the length of the prospective buffer based on number of cell in either orientations (circumferential or axial).

		grid.added_info_in_send_buf = 4;///Number of elements containing additional information at the beginning of the send buffer.
		int seg_config_s, seg_config_e;	///Integers to decided whether the row or column being sent is overlapping or exactly divisible into two halfs.
		/// data to send to the neighbour in UP1 direction
		extent_s = (int) (ceil((double) (grid.num_smc_circumferentially) / 2));
		extent_e = (int) (ceil((double) (grid.num_ec_circumferentially) / 2));

		/// The seg_config variables are to recording the configuration of the split of the buffering in each direction.
		/// If the total number of cell in on a face (UP, DOWN, LEFT or RIGHT) are EVEN (i.e. seg_congif=0), the split will be non-overlapping
		/// (eg. if total SMCs are 26 in UP direction, UP1 buffer will send 13 and UP2 will send the other 13 to corresponding nbrs.
		/// If the total number of cells is ODD (i.e. seg_congif=1), then the split will be over lapping.
		/// (eg. if total SMCs are 13 (or any multiple of 13) UP1 will send elements from 0 - 6 and UP2 will send 6 - 12  to corresponding nbrs.
		/// These variables are used in send and recv buffers update before and after the MPI-communication routine is called.
		seg_config_s = grid.num_smc_circumferentially % 2;
		seg_config_e = grid.num_ec_circumferentially % 2;

		///Recording the number of elements in Send buffer in Up direction (UP1 or UP2) for use in update routine as count of elements.
		grid.num_elements_send_up = grid.added_info_in_send_buf
				+ (grid.num_coupling_species_smc * extent_s
						+ grid.num_coupling_species_ec * extent_e);

		sendbuf[UP1] = (double*) checked_malloc(
				grid.num_elements_send_up * sizeof(double),
				"sendbuf[UP1] dimension 2");

		sendbuf[UP1][0] = (double) (1); //Start of the 1st segment of SMC array in  in UP direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[UP1][1] = (double) (extent_s); //End of the 1st segment of SMC array in UP direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[UP1][2] = (double) (1); //Start of the 1st segment of EC array in UP direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[UP1][3] = (double) (extent_e); ///End of the 1st segment of EC array in UP direction (circumferential direction) to be sent to neighbouring processor

		sendbuf[UP2] = (double*) checked_malloc(
				grid.num_elements_send_up * sizeof(double),
				"sendbuf[UP2] dimension 2");

		if (seg_config_s != 0) {
			sendbuf[UP2][0] = (double) (extent_s); //Start of the 2nd segment of SMC array in UP direction (circumferential direction) to be sent to neighbouring processor
		} else if (seg_config_s == 0) {
			sendbuf[UP2][0] = (double) (extent_s + 1); //Start of the 2nd segment of SMC array in UP direction (circumferential direction) to be sent to neighbouring processor
		}
		sendbuf[UP2][1] = (double) (grid.num_smc_circumferentially); //End of the 2nd segment of SMC array in  in UP direction (circumferential direction) to be sent to neighbouring processor

		if (seg_config_e != 0) {
			sendbuf[UP2][2] = (double) (extent_e); //Start of the 2nd segment of EC array in  in UP direction (circumferential direction) to be sent to neighbouring processor
		} else if (seg_config_e == 0) {
			sendbuf[UP2][2] = (double) (extent_e + 1); //Start of the 2nd segment of EC array in  in UP direction (circumferential direction) to be sent to neighbouring processor
		}
		sendbuf[UP2][3] = (double) (grid.num_ec_circumferentially); //End of the 2nd segment of EC array in  in UP direction (circumferential direction) to be sent to neighbouring processor


		/// data to send to the neighbour in DOWN direction
		///Recording the number of elements in Send buffer in DOWN direction (DOWN1 or DOWN2) for use in update routine as count of elements.
		grid.num_elements_send_down = grid.added_info_in_send_buf
				+ (grid.num_coupling_species_smc * extent_s
						+ grid.num_coupling_species_ec * extent_e);

		sendbuf[DOWN1] = (double*) checked_malloc(
				grid.num_elements_send_down * sizeof(double),
				"sendbuf[DOWN1] dimension 2");

		sendbuf[DOWN1][0] = (double) (1); //Start of the 1st segment of SMC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[DOWN1][1] = (double) (extent_s); //End of the 1st segment of SMC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[DOWN1][2] = (double) (1); //Start of the 1st segment of EC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[DOWN1][3] = (double) (extent_e); ///End of the 1st segment of EC array in DOWN direction (circumferential direction) to be sent to neighbouring processor

		sendbuf[DOWN2] = (double*) checked_malloc(
				grid.num_elements_send_down * sizeof(double),
				"sendbuf[DOWN2] dimension 2");
		if (seg_config_s != 0) {
			sendbuf[DOWN2][0] = (double) (extent_s); //Start of the 2nd segment of SMC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
		} else if (seg_config_s == 0) {
			sendbuf[DOWN2][0] = (double) (extent_s + 1); //Start of the 2nd segment of SMC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
		}

		sendbuf[DOWN2][1] = (double) (grid.num_smc_circumferentially); //End of the 2nd segment of SMC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor
		if (seg_config_e != 0) {
			sendbuf[DOWN2][2] = (double) (extent_e); //Start of the 2nd segment of EC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor
		} else if (seg_config_e == 0) {
			sendbuf[DOWN2][2] = (double) (extent_e + 1); //Start of the 2nd segment of EC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor
		}
		sendbuf[DOWN2][3] = (double) (grid.num_ec_circumferentially); //End of the 2nd segment of EC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor

		/// data to send to the neighbour in LEFT direction
		extent_s = (int) (ceil((double) (grid.num_smc_axially) / 2));
		extent_e = (int) (ceil((double) (grid.num_ec_axially) / 2));
		seg_config_s = grid.num_smc_axially % 2;
		seg_config_e = grid.num_ec_axially % 2;

		///Recording the number of elements in Send buffer in Left direction (LEFT1 or LEFT2) for use in update routine as count of elements.
		grid.num_elements_send_left = grid.added_info_in_send_buf
				+ (grid.num_coupling_species_smc * extent_s
						+ grid.num_coupling_species_ec * extent_e);
		sendbuf[LEFT1] = (double*) checked_malloc(
				grid.num_elements_send_left * sizeof(double),
				"sendbuf[LEFT1] dimension 2");
		sendbuf[LEFT1][0] = (double) (1); //Start of the 1st segment of SMC array in  in LEFT direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[LEFT1][1] = (double) (extent_s); //END of the 1st segment of SMC array in  in LEFT direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[LEFT1][2] = (double) (1); //Start of the 1st segment of EC array in  in LEFT direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[LEFT1][3] = (double) (extent_e); //END of the 1st segment of EC array in  in LEFT direction (circumferential direction) to be sent to neighbouring processor

		sendbuf[LEFT2] = (double*) checked_malloc(
				grid.num_elements_send_left * sizeof(double),
				"sendbuf[LEFT2] dimension 2");
		if (seg_config_s != 0) {
			sendbuf[LEFT2][0] = (double) (extent_s); //Start of the 2nd segment of SMC array in LEFT direction (circumferential direction) to be sent to neighbouring processor
		} else if (seg_config_s == 0) {
			sendbuf[LEFT2][0] = (double) (extent_s + 1); //Start of the 2nd segment of SMC array in LEFT direction (circumferential direction) to be sent to neighbouring processor
		}
		sendbuf[LEFT2][1] = (double) (grid.num_smc_axially); //END of the 2nd segment of SMC array in LEFT direction (circumferential direction) to be sent to neighbouring processor
		if (seg_config_e != 0) {
			sendbuf[LEFT2][2] = (double) (extent_e); //Start of the 2nd segment of EC array in LEFT direction (circumferential direction) to be sent to neighbouring processor
		} else if (seg_config_e == 0) {
			sendbuf[LEFT2][2] = (double) (extent_e + 1); //Start of the 2nd segment of EC array in LEFT direction (circumferential direction) to be sent to neighbouring processor
		}
		sendbuf[LEFT2][3] = (double) (grid.num_ec_axially); //END of the 2nd segment of EC array in LEFT direction (circumferential direction) to be sent to neighbouring processor

		/// data to send to the neighbour in RIGHT direction

		///Recording the number of elements in Send buffer in RIGHT direction (RIGHT1 or RIGHT2) for use in update routine as count of elements.
		grid.num_elements_send_right = grid.added_info_in_send_buf
				+ (grid.num_coupling_species_smc * extent_s
						+ grid.num_coupling_species_ec * extent_e);

		sendbuf[RIGHT1] = (double*) checked_malloc(
				grid.num_elements_send_right * sizeof(double),
				"sendbuf[RIGHT1] dimension 2");
		sendbuf[RIGHT1][0] = (double) (1); //Start of the 1st segment of SMC array in  in RIGHT direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[RIGHT1][1] = (double) (extent_s); //END of the 1st segment of SMC array in  in RIGHT direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[RIGHT1][2] = (double) (1); //Start of the 1st segment of EC array in  in RIGHT direction (circumferential direction) to be sent to neighbouring processor
		sendbuf[RIGHT1][3] = (double) (extent_e); //END of the 1st segment of EC array in  in RIGHT direction (circumferential direction) to be sent to neighbouring processor

		sendbuf[RIGHT2] = (double*) checked_malloc(
				grid.num_elements_send_right * sizeof(double),
				"sendbuf[RIGHT2] dimension 2");
		if (seg_config_s != 0) {
			sendbuf[RIGHT2][0] = (double) (extent_s); //Start of the 2nd segment of SMC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor
		} else if (seg_config_s == 0) {
			sendbuf[RIGHT2][0] = (double) (extent_s + 1); //Start of the 2nd segment of SMC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor
		}
		sendbuf[RIGHT2][1] = (double) (grid.num_smc_axially); //END of the 2nd segment of SMC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor
		if (seg_config_e != 0) {
			sendbuf[RIGHT2][2] = (double) (extent_e); //Start of the 2nd segment of EC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor
		} else if (seg_config_e == 0) {
			sendbuf[RIGHT2][2] = (double) (extent_e + 1); //Start of the 2nd segment of EC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor
		}
		sendbuf[RIGHT2][3] = (double) (grid.num_ec_axially); //END of the 2nd segment of EC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor


		///Call communication to the number of elements to be recieved by neighbours and allocate memory of recvbuf for each direction accordingly.
			grid = communicate_num_recv_elements_to_nbrs(grid);
			///memory allocation

		/// data to receive from the neighbour in UP direction
		recvbuf[UP1] = (double*) checked_malloc(
				grid.num_elements_recv_up * sizeof(double),
				"recvbuf[UP1] dimension 2");
		recvbuf[UP2] = (double*) checked_malloc(
				grid.num_elements_recv_up * sizeof(double),
				"recvbuf[UP2] dimension 2");

		/// data to recv from the neighbour in DOWN direction
		recvbuf[DOWN1] = (double*) checked_malloc(
				grid.num_elements_recv_down * sizeof(double),
				"recvbuf[DOWN1] dimension 2");
		recvbuf[DOWN2] = (double*) checked_malloc(
				grid.num_elements_recv_down * sizeof(double),
				"recvbuf[DOWN2] dimension 2");

		/// data to receive from the neighbour in LEFT direction
		recvbuf[LEFT1] = (double*) checked_malloc(
				grid.num_elements_recv_left * sizeof(double),
				"recvbuf[LEFT1] dimension 2");
		recvbuf[LEFT2] = (double*) checked_malloc(
				grid.num_elements_recv_left * sizeof(double),
				"recvbuf[LEFT2] dimension 2");

		/// data to receive from the neighbour in RIGHT direction
		recvbuf[RIGHT1] = (double*) checked_malloc(
				grid.num_elements_recv_right * sizeof(double),
				"recvbuf[RIGHT1] dimension 2");
		recvbuf[RIGHT2] = (double*) checked_malloc(
				grid.num_elements_recv_right * sizeof(double),
				"recvbuf[RIGHT2] dimension 2");

	grid.NEQ =	grid.neq_smc*(grid.num_smc_axially*grid.num_smc_circumferentially) + grid.neq_ec*(grid.num_ec_axially*grid.num_ec_circumferentially);

	///Setup output streams to write data in files. Each node opens an independent set of files and write various state variables into it.
//	checkpoint_handle *check = initialise_checkpoint(myRank);


	///Setting up the solver

	double 	tnow	= 0.0;

	//Error control variables
	double 	TOL	= 1e-6, absTOL = 1e-7;

	double * thres = (double*)checked_malloc(grid.NEQ*sizeof(double),"Threshod array for RKSUITE");

	for (int i=0; i<grid.NEQ; i++)
		thres[i]	=	absTOL;

	//Variables holding new and old values
	#ifdef CVODE
	N_Vector ny;
	ny = N_VNew_Serial(grid.NEQ);
	double *y = NV_DATA_S(ny);
	#else
	double* y =  (double*)checked_malloc(grid.NEQ*sizeof(double),"Solver array y for RKSUITE");
	#endif

	double* yp =  (double*)checked_malloc(grid.NEQ*sizeof(double),"Solver array y for RKSUITE");
	///Initialize different state variables and coupling data values.
		Initialize_koeingsberger_smc(grid,y,smc);
		Initialize_koeingsberger_ec(grid,y,ec);

		map_solver_to_cells(grid,y,smc,ec);

int state 	=  couplingParms(CASE,&cpl_cef);
dump_rank_info(check,cpl_cef,grid);



double t1	=	MPI_Wtime();
#ifdef CVODE
cvode_solver(tnow, tfinal, interval, ny, grid.NEQ, TOL, absTOL,file_write_per_unit_time,check);
#endif
#ifndef CVODE
rksuite_solver_CT(tnow, tfinal, interval, y, yp, grid.NEQ , TOL, thres, file_write_per_unit_time, check);
//rksuite_solver_UT(tnow, tfinal, interval, y, yp, grid.NEQ,TOL,thres, file_write_per_unit_time,check);
#endif
double t2	=	MPI_Wtime();
	final_checkpoint(grid,check, t1, t2);
MPI_Finalize();
}// end main()
























