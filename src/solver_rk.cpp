#ifndef CVODE
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"

extern "C" {
#include "rksuite.h"
}

extern conductance cpl_cef;
extern celltype1** smc;
extern celltype2** ec;
extern double **sendbuf, **recvbuf;
extern grid_parms grid;
extern time_stamps t_stamp;

///***************************************************************************************/
///************ComputeDerivates(tnow, state_variables[], first_derivatives[])*************/
///***************************************************************************************/
void computeDerivatives(double t, double y[], double f[]) {

	//compute_with_time_profiling(&t_stamp, grid, smc, ec, cpl_cef, t, y, f);
	compute(grid, smc, ec, cpl_cef, t, y, f);
	t_stamp.computeDerivatives_call_counter += 1;

}    //end of computeDerivatives()
void rksuite_solver_CT(double tnow, double tfinal, double interval, double *y,
		double* yp, int total, double TOL, double* thres,
		int file_write_per_unit_time, int line_number,
		checkpoint_handle *check) {

	RKSUITE rksuite;
	//Solver method
	int method = 2;		//RK(4,5)
	double tend;
	int cflag = 0;
	int itteration = 0;
	int write_count = 0;
	int write_once = 0;
	int count = 0;
	initialize_t_stamp(&t_stamp);
	tend = interval;
	rksuite.setup(total, tnow, y, tend, TOL, thres, method, "CT", false, 0.0,
			false);
	communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);
	if (grid.universal_rank==0){
		grid.logptr = fopen("logfile.txt","w+");
		print_domains(grid.logptr, grid, smc,ec);
		fclose(grid.logptr);
	}
	//computeDerivatives(tnow, y, yp);
	MPI_Barrier(grid.universe);
	int file_offset_for_timing_data = determine_file_offset_for_timing_data(
			check, grid);
	int totf, stpcst, stpsok;
	double waste, hnext;
	while (tnow <= tfinal) {
		// the ct() function does not guarantee to advance all the
		// way to the stop time.  Keep stepping until it does.
		t_stamp.solver_t1 = MPI_Wtime();
		do {
			rksuite.ct(computeDerivatives, tnow, y, yp, cflag);
			if (cflag >= 5) {
				fprintf(stdout,
						"[%d] \t RKSUITE failed with error flag %d at t=%lf\n\n",
						grid.rank, cflag, tnow);
				MPI_Abort(MPI_COMM_WORLD, 300);
			}
		} while (tnow < tend);
		/// rksuite.stat() routine calls the to gather statistic on the performance of the solver.
		/// Amongst other things it also informs about what the next step size should be.
		rksuite.stat(totf, stpcst, waste, stpsok, hnext);
		t_stamp.solver_t2 = MPI_Wtime();
		t_stamp.diff_solver = t_stamp.solver_t2 - t_stamp.solver_t1;
		///Increament the itteration as rksuite has finished solving between bounds tnow<= t <= tend.
		itteration++;
		/// Call for interprocessor communication
		/*t_stamp.barrier_in_solver_before_comm_t1 = MPI_Wtime();
		 MPI_Barrier(grid.universe);
		 t_stamp.barrier_in_solver_before_comm_t2 = MPI_Wtime();

		 t_stamp.diff_barrier_in_solver_before_comm =
		 t_stamp.barrier_in_solver_before_comm_t2
		 - t_stamp.barrier_in_solver_before_comm_t1;*/

		 t_stamp.total_comms_cost_t1 = MPI_Wtime();
		 communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);
		 t_stamp.total_comms_cost_t2 = MPI_Wtime();
		 t_stamp.diff_total_comms_cost = t_stamp.total_comms_cost_t2
		 - t_stamp.total_comms_cost_t1;
		/*if (itteration == 5) {
		 dump_JPLC(grid, ec, check, "Local agonist before t=100s\n");
		 }*/

		if ((write_once <= 1) && (tnow >= grid.stimulus_onset_time)) {
			write_once++;
			if (grid.rank % grid.n == 0) {
				dump_JPLC(grid, ec, check, "Local agonist after t=100s");
			}
		}

		t_stamp.write_t1 = MPI_Wtime();
		if ((itteration % file_write_per_unit_time) == 0) {
			dump_data(check, grid, line_number, tnow, smc, ec, write_count);
			update_line_number(check, grid, write_count);
			write_count++;
		}		//end itteration
		t_stamp.write_t2 = MPI_Wtime();
		t_stamp.diff_write = t_stamp.write_t2 - t_stamp.write_t1;

		//checkpoint_timing_data(grid, check, tnow, t_stamp, count,	file_offset_for_timing_data);
		initialize_t_stamp(&t_stamp);
		count++;
		tend += interval;
		rksuite.reset(tend);
	}			//end while()

}

void rksuite_solver_UT(double tnow, double tfinal, double interval, double *y,
		double* yp, int total, double TOL, double* thres,
		int file_write_per_unit_time, int line_number,
		checkpoint_handle *check) {
	printf("[%d]: I have called RKSUITE\n", grid.universal_rank);

	RKSUITE rksuite;
	//Solver method
	int method = 2;		//RK(4,5)
	//Error Flag
	int uflag = 0;
	int itteration = 0;
	double tend;
	int write_count = 0;
	int write_once = 0;
	double* ymax = (double*) checked_malloc(grid.NEQ * sizeof(double),
			"Solver array ymax for RKSUITE");

	communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);
	computeDerivatives(tnow, y, yp);
	rksuite.setup(grid.NEQ, tnow, y, tfinal, TOL, thres, method, "UT", false,
			0.00, false);

	///Iterative  calls to the solver start here.
	for (tend = interval; tend < tfinal; tend += interval) {

		/// RKSUITE UT Call
		/// prototype (f,twant,tgot,ygot,ypgot,ymax,work,uflag)
		rksuite.ut(computeDerivatives, tend, tnow, y, yp, ymax, uflag);
		if (uflag >= 5) {
			fprintf(stdout,
					"[%d] \t RKSUITE failed with error flag %d at t=%lf\n\n",
					grid.rank, uflag, tnow);
			MPI_Abort(MPI_COMM_WORLD, 300);
		}

		///Increament the itteration as rksuite has finished solving between bounds tnow<= t <= tend.
		itteration++;
		/// Call for interprocessor communication
		MPI_Barrier(grid.universe);
		communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);

		/*if (itteration == 5) {
		 dump_JPLC(grid, ec, check, "Local agonist before t=100s");
		 }

		 if ((itteration == int(grid.stimulus_onset_time/interval)) && (write_once<=1)) {
		 write_once++;
		 dump_JPLC(grid, ec, check, "Local agonist after t=100s");
		 }*/

		if ((itteration % file_write_per_unit_time) == 0) {
			dump_data(check, grid, line_number, tnow, smc, ec, write_count);
			write_count++;
		}		//end itteration
		//MPI_Barrier(grid.universe);
	}		//end of for loop on TEND

}

void rksuite_solver_CT_debug(double tnow, double tfinal, double interval,
		double *y, double* yp, int total, double TOL, double* thres,
		int file_write_per_unit_time, int line_number,
		checkpoint_handle *check) {

	RKSUITE rksuite;
	//Solver method
	int method = 2;		//RK(4,5)
	double tend;
	int cflag = 0;
	int itteration = 0;
	int write_count = 0;
	int write_once = 0;
	int count = 0;
	initialize_t_stamp(&t_stamp);
	tend = interval;
	rksuite.setup(total, tnow, y, tend, TOL, thres, method, "CT", false, 0.0,
			false);

	int totf, stpcst, stpsok;
	double waste, hnext;
	dump_rank_info(check, cpl_cef, grid);

	while (tnow <= tfinal) {
		// the ct() function does not guarantee to advance all the
		// way to the stop time.  Keep stepping until it does.
		t_stamp.solver_t1 = MPI_Wtime();
		do {
			rksuite.ct(computeDerivatives, tnow, y, yp, cflag);
			if (cflag >= 5) {
				fprintf(stdout,
						"[%d] \t RKSUITE failed with error flag %d at t=%lf\n\n",
						grid.rank, cflag, tnow);
				MPI_Abort(MPI_COMM_WORLD, 300);
			}
		} while (tnow < tend);
		/// rksuite.stat() routine calls the to gather statistic on the performance of the solver.
		/// Amongst other things it also informs about what the next step size should be.
		rksuite.stat(totf, stpcst, waste, stpsok, hnext);
		t_stamp.solver_t2 = MPI_Wtime();
		t_stamp.diff_solver = t_stamp.solver_t2 - t_stamp.solver_t1;
		///Increament the itteration as rksuite has finished solving between bounds tnow<= t <= tend.
		itteration++;

		/// Call for interprocessor communication
		/*	t_stamp.barrier_in_solver_before_comm_t1 = MPI_Wtime();
		 MPI_Barrier(grid.universe);
		 t_stamp.barrier_in_solver_before_comm_t2 = MPI_Wtime();

		 t_stamp.diff_barrier_in_solver_before_comm =
		 t_stamp.barrier_in_solver_before_comm_t2
		 - t_stamp.barrier_in_solver_before_comm_t1;

		 t_stamp.total_comms_cost_t1 = MPI_Wtime();
		 communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);
		 t_stamp.total_comms_cost_t2 = MPI_Wtime();
		 t_stamp.diff_total_comms_cost = t_stamp.total_comms_cost_t2
		 - t_stamp.total_comms_cost_t1;*/
		/*if (itteration == 5) {
		 dump_JPLC(grid, ec, check, "Local agonist before t=100s\n");
		 }*/
		t_stamp.write_t1 = MPI_Wtime();
		if ((itteration % file_write_per_unit_time) == 0) {
			dump_data(check, grid, line_number, tnow, smc, ec, write_count);
			//update_line_number(check, grid, write_count);
			write_count++;
		}		//end itteration
		t_stamp.write_t2 = MPI_Wtime();
		t_stamp.diff_write = t_stamp.write_t2 - t_stamp.write_t1;
		initialize_t_stamp(&t_stamp);
		count++;
		tend += interval;
		rksuite.reset(tend);
	}	//end while()

}

#endif
