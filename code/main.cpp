#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"


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
void* checked_malloc(size_t bytes, FILE* errfile, const char* errmsg){
	void *pval = malloc(bytes);

	if (pval == NULL) {
		fprintf(errfile, "Allocation failed for %s\n", errmsg);
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

    //// These are input parameters for the simulation
     int
     m = 5,
     n = 4,
     e = 4,	///number of ECs per node
     s = 4;	///number of SMCs per node



//Reveal information of myself and size of MPI_COMM_WORLD
    check_flag(MPI_Comm_rank(grid.universe, &grid.universal_rank), stdout,
	    "error Comm_rank");
    check_flag(MPI_Comm_size(grid.universe, &grid.numtasks), stdout,
	    "error Comm_size");

///Initialize checkpoint routine which opens files
    checkpoint_handle *check = initialise_checkpoint(grid.universal_rank);

//Test case
//config
    int num_subdomains = 4;
    int **domains;

    domains = (int**) checked_malloc(num_subdomains * sizeof(int*), stdout,
	    "Subdomain information allocation");

    /*second coordinate has following information:
     * Element 1: 	Key_val or serial number of the subdomain
     * Element 2:	Subdomain Type (4 possibilities and their values)
     * 				1. Straight Segment														(0)
     * 				2. Bifurcation															(1)
     * Element 3:	Axial extent of processor of current key_val
     * Element 4: 	circumferential extent of processors of current key_val
     * Element 5:	Parent subdomain key_val of current Key_val.
     * Element 6: 	In case if the Parent is a bifurcation, which branch owns the current domain as child (i.e. which branch is my parent)
     * Element 7: 	Left Child subdomain key_val of the current Key_val.
     * Element 8:   Right Child subdomain key_val of the current Key_val.
     *
     * In the case of elements 6 & 7, if subdomain type of current key_val is a straight segment, left Child is positive or zero, and right Child is negative.
     * If subdomain type of current key_val is a bifurcation, then both right and left child subdomains are non-negative.
     */
    for (int i = 0; i < num_subdomains; i++)
	{
	domains[i] = (int*) checked_malloc(7 * sizeof(int), stdout,
		"Subdomains array elements allocation");
	}
    domains[0][0] = 0;
    domains[0][1] = STRSEG;
    domains[0][2] = m;
    domains[0][3] = n;
    domains[0][4] = -1;
    domains[0][5] = 1;
    domains[0][6] = -1;

    domains[1][0] = 1;
    domains[1][1] = BIF;
    domains[1][2] = m;
    domains[1][3] = n;
    domains[1][4] = 0;
    domains[1][5] = 2;
    domains[1][6] = 3;

    domains[2][0] = 2;
    domains[2][1] = BIF;
    domains[2][2] = m;
    domains[2][3] = n;
    domains[2][4] = 1;
    domains[2][5] = -1;
    domains[2][6] = -1;

    domains[3][0] = 3;
    domains[3][1] = STRSEG;
    domains[3][2] = m;
    domains[3][3] = n;
    domains[3][4] = 1;
    domains[3][5] = -1;
    domains[3][6] = -1;

    grid =  make_subdomains(grid, num_subdomains, domains, check->logptr);

///Time variables
	double tfinal = 1.0;
	double interval = 1e-2;
//File written every 1 second
	int file_write_per_unit_time = 1;//int(0.5/interval);

	grid.uniform_jplc = 0.1, grid.min_jplc = 0.35, grid.max_jplc = 0.4, grid.gradient =
			2.5e-2; grid.stimulus_onset_time	=100.00;


grid = set_geometry_parameters(grid,check->logptr,e,s);

	if (grid.my_domain.internal_info.domain_type == STRSEG) {
		grid = make_straight_segment(grid, check->logptr);
	} else if (grid.my_domain.internal_info.domain_type == BIF) {
		grid = make_bifucation(grid, check->logptr);
	}


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

int a = allocate_memory(grid,smc,ec,sendbuf,recvbuf,check);
if (a == success){
	fprintf(check->logptr,"[%d] error returned in memory allocation\n",grid.universal_rank);
	MPI_Abort(MPI_COMM_WORLD, 100);
}

	grid.NEQ =	grid.neq_smc*(grid.num_smc_axially*grid.num_smc_circumferentially) + grid.neq_ec*(grid.num_ec_axially*grid.num_ec_circumferentially);

	///Setup output streams to write data in files. Each node opens an independent set of files and write various state variables into it.
//	checkpoint_handle *check = initialise_checkpoint(myRank);


	///Setting up the solver

	double 	tnow	= 0.0;

	//Error control variables
	double 	TOL	= 1e-6;

	double * thres = (double*)checked_malloc(grid.NEQ*sizeof(double),stdout,"Threshod array for RKSUITE");

	for (int i=0; i<grid.NEQ; i++)
		thres[i]	=	1e-6;

	//Variables holding new and old values
	double* y =  (double*)checked_malloc(grid.NEQ*sizeof(double),stdout,"Solver array y for RKSUITE");
	double* yp = (double*) checked_malloc(grid.NEQ * sizeof(double), stdout,"Solver array yp for RKSUITE");



	///Initialize different state variables and coupling data values.
		Initialize_koeingsberger_smc(grid,y,smc);
		Initialize_koeingsberger_ec(grid,y,ec);

		map_solver_to_cells(grid,y,smc,ec);




int state 	=  couplingParms(CASE,&cpl_cef);
dump_rank_info(check,cpl_cef,grid);

/*communication_async_send_recv(check->logptr,grid,sendbuf,recvbuf,smc,ec);
print_domains(check->logptr,grid,smc,ec);*/
double t1	=	MPI_Wtime();

	rksuite_solver_CT(tnow, tfinal, interval, y, yp, grid.NEQ , TOL, thres, file_write_per_unit_time, check);

	//rksuite_solver_UT(tnow, tfinal, interval, y, yp, grid.NEQ,TOL,thres, file_write_per_unit_time,check);

double t2	=	MPI_Wtime();
	final_checkpoint(check, grid,t1, t2);
	
MPI_Finalize();
}// end main()
























