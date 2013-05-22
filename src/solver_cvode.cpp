#ifdef CVODE

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"

#include <cvode/cvode.h>
#include <cvode/cvode_spgmr.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>


extern conductance cpl_cef;
extern celltype1** smc;
extern celltype2** ec;
extern double **sendbuf, **recvbuf;
extern grid_parms grid;
extern time_stamps		t_stamp;



/********************************************************************************************/
/**/        static int check_cvode_flag(void *flagvalue, char *funcname, int opt)               /**/
/********************************************************************************************/
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}

///***************************************************************************************/
///************ComputeDerivates(tnow, state_variables[], first_derivatives[])*************/
///***************************************************************************************/
void computeDerivatives(double t, double y[], double f[]) {


	const double gama = 1970.00;
	const double lambda = 45.00;
	const double Cmj = 25.8;
	int k = 0, talk = 0;

	t_stamp.computeDerivatives_call_counter = t_stamp.computeDerivatives_call_counter +1;

	t_stamp.map_function_t1	=	MPI_Wtime();
	map_solver_to_cells(grid,y,smc,ec);
	t_stamp.map_function_t2	=	MPI_Wtime();
	t_stamp.diff_map_function	=	t_stamp.diff_map_function + (t_stamp.map_function_t2 - t_stamp.map_function_t1);


	t_stamp.single_cell_fluxes_t1	=	MPI_Wtime();
	single_cell( t, y,  grid, smc,  ec);
	t_stamp.single_cell_fluxes_t2	=	MPI_Wtime();
	t_stamp.diff_single_cell_fluxes = t_stamp.diff_single_cell_fluxes + (t_stamp.single_cell_fluxes_t2 - t_stamp.single_cell_fluxes_t1);

	t_stamp.coupling_fluxes_t1	=	MPI_Wtime();
	coupling(t,y, grid, smc, ec,cpl_cef);
	t_stamp.coupling_fluxes_t2	=	MPI_Wtime();
	t_stamp.diff_coupling_fluxes = t_stamp.diff_coupling_fluxes + (t_stamp.coupling_fluxes_t2 - t_stamp.coupling_fluxes_t1);

	///In previous version of code (agonist_variation_main.cpp) the indexing in this section was
	///handelled differently but the end result is same.
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
			k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
			k = 0;
			f[k + ((j - 1) * grid.neq_smc) +smc_Ca]= smc[i][j].A[J_IP3] - smc[i][j].A[J_SERCA]
			+ smc[i][j].A[J_CICR] - smc[i][j].A[J_Extrusion]
			+ smc[i][j].A[J_Leak] - smc[i][j].A[J_VOCC]
			+ smc[i][j].A[J_Na_Ca] + smc[i][j].B[cpl_Ca]
			+ smc[i][j].C[cpl_Ca];

			f[k + ((j - 1) * grid.neq_smc) +smc_SR]= smc[i][j].A[J_SERCA] - smc[i][j].A[J_CICR]
			- smc[i][j].A[J_Leak];

			f[k + ((j - 1) * grid.neq_smc) +smc_Vm]= gama
			* (-smc[i][j].A[J_Na_K] - smc[i][j].A[J_Cl]
					- (2 * smc[i][j].A[J_VOCC]) - smc[i][j].A[J_Na_Ca]
					- smc[i][j].A[J_K]) + smc[i][j].B[cpl_Vm]
			+ smc[i][j].C[cpl_Vm];

			f[k + ((j - 1) * grid.neq_smc) +smc_w] = lambda
			* (smc[i][j].A[K_activation] - smc[i][j].p[smc_w]);

			f[k + ((j - 1) * grid.neq_smc) +smc_IP3] = -smc[i][j].A[J_IP3_deg]
			+ smc[i][j].B[cpl_IP3] + smc[i][j].C[cpl_IP3];
		}
	}

	int offset = (grid.neq_smc * grid.num_smc_circumferentially
			* grid.num_smc_axially);

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			if (i > 1)
			k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
			k = offset + 0;

			ec[i][j].JPLC = agonist_profile(t,grid,i,j,ec[i][j].z_coord);

			f[k + ((j - 1) * grid.neq_ec) +ec_Ca] = ec[i][j].A[J_IP3]
			- ec[i][j].A[J_SERCA] + ec[i][j].A[J_CICR]
			- ec[i][j].A[J_Extrusion] + ec[i][j].A[J_Leak]
			+ ec[i][j].A[J_NSC] + ec[i][j].A[J_trivial_Ca]
			+ ec[i][j].B[cpl_Ca] + ec[i][j].C[cpl_Ca];

			f[k + ((j - 1) * grid.neq_ec) +ec_SR] = ec[i][j].A[J_SERCA]
			- ec[i][j].A[J_CICR] - ec[i][j].A[J_Leak];

			f[k + ((j - 1) * grid.neq_ec) +ec_Vm] = ((-1 / Cmj)
					* (ec[i][j].A[J_Ktot] + ec[i][j].A[J_Residual]))
			+ ec[i][j].B[cpl_Vm] + ec[i][j].C[cpl_Vm];

			f[k + ((j - 1) * grid.neq_ec) +ec_IP3] =	 ec[i][j].JPLC
			- ec[i][j].A[J_IP3_deg] + ec[i][j].B[cpl_IP3]
			+ ec[i][j].C[cpl_IP3];

		}
	}
}    //end of computeDerivatives()

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data) {
	double tnow = (double)(t);
	computeDerivatives(tnow, NV_DATA_S(y), NV_DATA_S(ydot));
	return 0;
}

void cvode_solver(double tnow, double tfinal, double interval, N_Vector y, int total, double TOL, double absTOL,
		int file_write_per_unit_time, checkpoint_handle *check){

void* cvode_mem;
int flag;
realtype t;
int itteration = 0;
int write_count=0;
int write_once=0;
int count=0;
initialize_t_stamp(t_stamp);
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_cvode_flag((void *) cvode_mem, "CVodeCreate", 0)) {
		MPI_Abort(MPI_COMM_WORLD, 321);
	}
	flag = CVodeInit(cvode_mem, f, tnow, y);
	if (check_cvode_flag(&flag, "CVodeInit", 1)) {
		MPI_Abort(MPI_COMM_WORLD, 321);
	}

	flag = CVodeSStolerances(cvode_mem, TOL, absTOL);
	if (check_cvode_flag(&flag, "CVodeSStolerances", 1))
		MPI_Abort(MPI_COMM_WORLD, 321);

	flag = CVSpgmr(cvode_mem, PREC_NONE, 0);
	if (check_cvode_flag(&flag, "CVSpgmr", 1)) {
		MPI_Abort(MPI_COMM_WORLD, 321);
	}

	flag = CVodeSetInitStep(cvode_mem, interval / 2);
	if (check_cvode_flag(&flag, "CVodeSetInitStep", 1)) {
		MPI_Abort(MPI_COMM_WORLD, 321);
	}

	///Iterative  calls to the solver start here.
	/* time stamp this*/
	t_stamp.solver_t1 	=	MPI_Wtime();
	for (double k = tnow; k < tfinal; k += interval) {
			 flag = CVode(cvode_mem, k, y, &t, CV_NORMAL);
			 if(check_cvode_flag(&flag, "CVode", 1)){
			 	MPI_Abort(MPI_COMM_WORLD, 400);
			 }
	t_stamp.solver_t2 	=	MPI_Wtime();
	t_stamp.diff_solver = t_stamp.solver_t2 - t_stamp.solver_t1;
		///Increament the itteration as rksuite has finished solving between bounds tnow<= t <= tend.
		itteration++;

		/// Call for interprocessor communication
		t_stamp.barrier_in_solver_before_comm_t1	=	MPI_Wtime();
		MPI_Barrier(grid.universe); 									/* time stamp this*/
		t_stamp.barrier_in_solver_before_comm_t2	=	MPI_Wtime();

		t_stamp.diff_barrier_in_solver_before_comm = t_stamp.barrier_in_solver_before_comm_t2 - t_stamp.barrier_in_solver_before_comm_t1;

		communication_async_send_recv(grid,sendbuf,recvbuf,smc,ec);

		/*if (itteration == 5) {
			dump_JPLC(grid, ec, check, "Local agonist before t=100s");
		}*/
		tnow = (double) (t);
		if (/*(itteration == int(grid.stimulus_onset_time/interval)) && */(write_once<=1) && (tnow>=grid.stimulus_onset_time)) {
			write_once++;
			if (grid.rank%grid.n == 0){
			dump_JPLC(grid, ec, check, "Local agonist after t=100s");
			}
		}


		/* time stamp this*/
		t_stamp.write_t1	=	MPI_Wtime();
		if ((itteration % file_write_per_unit_time) == 0) {
					checkpoint(check, grid, tnow, smc, ec,write_count);
				write_count++;
				}		//end itteration
		t_stamp.write_t2	=	MPI_Wtime();
		t_stamp.diff_write  =   t_stamp.write_t2-t_stamp.write_t1;

		checkpoint_timing_data(grid,check,tnow,t_stamp,count);
		initialize_t_stamp(t_stamp);
		count++;
		//MPI_Barrier(grid.universe);
	}		//end of for loop on TEND
}
#endif	/* CVODE */
