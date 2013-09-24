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
extern time_stamps t_stamp;

/********************************************************************************************/
/**/static int check_cvode_flag(void *flagvalue, char *funcname, int opt) /**/
/********************************************************************************************/
{
	int *errflag;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL) {
		fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1);}

	/* Check if flag < 0 */
	else if (opt == 1) {
		errflag = (int *) flagvalue;
		if (*errflag < 0) {
			fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
					funcname, *errflag);
			return(1);}}

	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && flagvalue == NULL) {
		fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname);
		return(1);}

	return(0);
}

///***************************************************************************************/
///************ComputeDerivates(tnow, state_variables[], first_derivatives[])*************/
///***************************************************************************************/
void computeDerivatives(double t, double y[], double f[]) {

	compute_with_time_profiling(&t_stamp, grid, smc, ec, cpl_cef, t, y,f);
	//compute(grid, smc, ec, cpl_cef, t, y,f);
	t_stamp.computeDerivatives_call_counter += 1;

}    //end of computeDerivatives()

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data) {
	double tnow = (double)(t);
	computeDerivatives(tnow, NV_DATA_S(y), NV_DATA_S(ydot));
	return 0;
}

void cvode_solver(double tnow, double tfinal, double interval, N_Vector y, int total, double TOL, double absTOL,
		int file_write_per_unit_time,int line_number, checkpoint_handle *check, time_keeper* elps_t) {

	void* cvode_mem;
	int flag;
	realtype t;
	int itteration = 0;
	int write_count=0;
	int write_once=0;
	int count=0;
	double hmin;
	initialize_t_stamp(&t_stamp);
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
	if (line_number > 0) {
		hmin = interval/1e3;

	} else if (line_number ==0) {
		hmin = interval/2;
	}
	flag = CVodeSetInitStep(cvode_mem, hmin);
	if (check_cvode_flag(&flag, "CVodeSetInitStep", 1)) {
		MPI_Abort(MPI_COMM_WORLD, 321);
	}

	int file_offset_for_timing_data = determine_file_offset_for_timing_data(check,grid);

	///Temporary insertion Time profiling min,max,avg
	double **time_profiler;
	time_profiler = (double**)checked_malloc(11*sizeof(double*),"Time profile array allocation failed.");
	for (int i =0;i<11;i++) {
		time_profiler[i] = (double*)checked_malloc(int((tfinal/1e-2)+(2/1e-2))*sizeof(double),"Time profile array column allocation failed.");
	}
	///Iterative  calls to the solver start here.
	for (double k = tnow; k < tfinal; k += interval) {
		t_stamp.solver_t1 = MPI_Wtime();
		flag = CVode(cvode_mem, k, y, &t, CV_NORMAL);
		if(check_cvode_flag(&flag, "CVode", 1)) {
			MPI_Abort(MPI_COMM_WORLD, 400);
		}
		t_stamp.solver_t2 = MPI_Wtime();
		t_stamp.diff_solver = t_stamp.solver_t2 - t_stamp.solver_t1;
		///Increament the itteration as rksuite has finished solving between bounds tnow<= t <= tend.
		itteration++;

		/// Call for interprocessor communication
		t_stamp.barrier_in_solver_before_comm_t1 = MPI_Wtime();
		MPI_Barrier(grid.universe); /* time stamp this*/
		t_stamp.barrier_in_solver_before_comm_t2 = MPI_Wtime();

		t_stamp.diff_barrier_in_solver_before_comm = t_stamp.barrier_in_solver_before_comm_t2 - t_stamp.barrier_in_solver_before_comm_t1;

		communication_async_send_recv(grid,sendbuf,recvbuf,smc,ec);

		/*if (itteration == 5) {
		 dump_JPLC(grid, ec, check, "Local agonist before t=100s");
		 }*/
		tnow = (double) (t);
		if ((write_once<=1) && (tnow>=grid.stimulus_onset_time)) {
			write_once++;
			if(((grid.rank+1)%grid.n)==0) {
				dump_JPLC(grid, ec, check, "Local agonist after t=100s");
			}
		}
		/* time stamp this*/
		t_stamp.write_t1 = MPI_Wtime();
		if ((itteration % file_write_per_unit_time) == 0) { 
			dump_data(check, grid, line_number,tnow, smc, ec,write_count);
			update_line_number(check, grid,write_count);
			write_count++;
		}		//end itteration
		t_stamp.write_t2 = MPI_Wtime();
		t_stamp.diff_write = t_stamp.write_t2-t_stamp.write_t1;


		checkpoint_timing_data(grid,check,tnow,t_stamp,count,file_offset_for_timing_data);
		//Record_timing_data_in_arrays(grid,tnow,t_stamp, itteration,time_profiler);
		initialize_t_stamp(&t_stamp);
		count++;
		//update_elapsed_time(check,grid,elps_t);
		//MPI_Barrier(grid.universe);
	}		//end of for loop on TEND

		process_time_profiling_data(grid,time_profiler,count);
}
#endif	/* CVODE */
