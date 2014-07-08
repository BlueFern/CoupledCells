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
void rksuite_solver_CT(double tnow, double tfinal, double interval, double *y, double* yp, int total, double TOL, double* thres,
		int file_write_per_unit_time, int line_number, checkpoint_handle *check, char* path, IO_domain_info* my_IO_domain_info) {

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
	rksuite.setup(total, tnow, y, tend, TOL, thres, method, "CT", false, 0.0, false);
	communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);

	//computeDerivatives(tnow, y, yp);
	MPI_Barrier(grid.universe);
	//int file_offset_for_timing_data = determine_file_offset_for_timing_data(check, grid);
	int totf, stpcst, stpsok;
	double waste, hnext;
	int err;
	data_buffer* writer_buffer = (data_buffer*) checked_malloc(sizeof(data_buffer), "allocation failed for data_buffer");
	writer_buffer->buffer_length = (int*) checked_malloc(12 * sizeof(int), "allocation failed for buffer_length array in data_buffer structure.");
	writer_buffer->smc_stat_var_buffer_length = (int*) checked_malloc((grid.neq_smc) * sizeof(int),
			"allocation failed for smc_buffer_length array in data_buffer structure.");
	writer_buffer->ec_stat_var_buffer_length = (int*) checked_malloc((grid.neq_ec) * sizeof(int),
			"allocation failed for ec_buffer_length array in data_buffer structure.");
	writer_buffer->smc_cpl = (int*) checked_malloc((grid.num_coupling_species_smc) * sizeof(int),
			"allocation failed for smc_cpl_buffer_length array in data_buffer structure.");
	writer_buffer->ec_cpl = (int*) checked_malloc((grid.num_coupling_species_ec) * sizeof(int),
			"allocation failed for ec_cpl_buffer_length array in data_buffer structure.");

	/// Dump MPI task mesh representation into vtk file to menifest task map.
	gather_tasks_mesh_point_data_on_writers(&grid, my_IO_domain_info, writer_buffer, smc, ec);
	if (grid.rank == 0) {
		dump_process_data(check, &grid, my_IO_domain_info, writer_buffer, path);
	}
	/// Dump JPLC map on bifurcation into a vtk file.
	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			ec[i][j].JPLC = agonist_profile((grid.stimulus_onset_time + 1), grid, i, j, ec[i][j].centeroid_point[1]);
		}
	}
	gather_ec_mesh_data_on_writers(&grid, my_IO_domain_info, writer_buffer, ec);
	gather_JPLC_map(&grid, my_IO_domain_info, writer_buffer, ec);
	if (grid.rank == 0) {
		dump_agonists_map(check, &grid, my_IO_domain_info, writer_buffer, ec, path);
	}
	double palce_holder_for_timing_max_min[3][int(tfinal / interval)];
/// ITERATION loop to go from INITIAL time to FINAL time.
	while (tnow <= tfinal) {
		// the ct() function does not guarantee to advance all the
		// way to the stop time.  Keep stepping until it does.
		t_stamp.solver_t1 = MPI_Wtime();
		do {
			rksuite.ct(computeDerivatives, tnow, y, yp, cflag);
			if (cflag >= 5) {
				fprintf(stdout, "[%d] \t RKSUITE failed with error flag %d at t=%lf\n\n", grid.rank, cflag, tnow);
				MPI_Abort(MPI_COMM_WORLD, 300);
			}
		} while (tnow < tend);
		/// rksuite.stat() routine calls the to gather statistic on the performance of the solver.
		/// Amongst other things it also informs about what the next step size should be.
		rksuite.stat(totf, stpcst, waste, stpsok, hnext);
		t_stamp.solver_t2 = MPI_Wtime();
		t_stamp.diff_solver = t_stamp.solver_t2 - t_stamp.solver_t1;
		palce_holder_for_timing_max_min[0][itteration] = t_stamp.diff_solver;
		t_stamp.aggregate_compute += t_stamp.diff_solver;

		/// Call for interprocessor communication
		t_stamp.total_comms_cost_t1 = MPI_Wtime();
		communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);
		t_stamp.total_comms_cost_t2 = MPI_Wtime();
		t_stamp.diff_total_comms_cost = t_stamp.total_comms_cost_t2 - t_stamp.total_comms_cost_t1;
		palce_holder_for_timing_max_min[1][itteration] = t_stamp.diff_total_comms_cost;
		t_stamp.aggregate_comm += t_stamp.diff_total_comms_cost;
		/*if (itteration == 5) {
		 dump_JPLC(grid, ec, check, "Local agonist before t=100s\n");
		 }*/

		if ((write_once <= 1) && (tnow >= grid.stimulus_onset_time)) {
			write_once++;
			if (grid.rank % grid.n == 0) {
				//dump_JPLC(grid, ec, check, "Local agonist after t=100s");
			}
		}

		if ((itteration % file_write_per_unit_time) == 0) {
			t_stamp.write_t1 = MPI_Wtime();

			gather_smc_mesh_data_on_writers(&grid, my_IO_domain_info, writer_buffer, smc);
			gather_ec_mesh_data_on_writers(&grid, my_IO_domain_info, writer_buffer, ec);
			gather_smcData(&grid, my_IO_domain_info, writer_buffer, smc, write_count);
			gather_ecData(&grid, my_IO_domain_info, writer_buffer, ec, write_count);

			if (grid.rank == 0) {
				initialise_time_wise_checkpoint(check, grid, write_count, path, my_IO_domain_info);
				dump_data(check, &grid, line_number, tnow, smc, ec, write_count, my_IO_domain_info, writer_buffer);
				close_time_wise_checkpoints(check);
			}

			t_stamp.write_t2 = MPI_Wtime();
			t_stamp.diff_write = t_stamp.write_t2 - t_stamp.write_t1;
			palce_holder_for_timing_max_min[2][write_count] = t_stamp.diff_write;
			t_stamp.aggregate_write += t_stamp.diff_write;
			write_count++;
			if (my_IO_domain_info->writer_rank == 0)
				cout << "tnow = " << tnow << endl;
		}		//end itteration

		//checkpoint_timing_data(grid, check, tnow, t_stamp, itteration, file_offset_for_timing_data);
		initialize_t_stamp(&t_stamp);
		///Increament the itteration as rksuite has finished solving between bounds tnow<= t <= tend.
		itteration++;
		tend += interval;
		rksuite.reset(tend);
	}			//end while()
	//t_stamp.aggregate_compute = t_stamp.aggregate_compute / itteration;
	//t_stamp.aggregate_comm = t_stamp.aggregate_comm / itteration;
	//t_stamp.aggregate_write = t_stamp.aggregate_write / write_count;

	double tmp_array[write_count];
	for (int i = 0; i < write_count; i++) {
		tmp_array[i] = palce_holder_for_timing_max_min[2][i];
	}
	maximum(palce_holder_for_timing_max_min[0], itteration, &t_stamp.max_compute, &t_stamp.max_compute_index);
	maximum(palce_holder_for_timing_max_min[1], itteration, &t_stamp.max_comm, &t_stamp.max_comm_index);
	maximum(tmp_array, write_count, &t_stamp.max_write, &t_stamp.max_write_index);

	minimum(palce_holder_for_timing_max_min[0], itteration, &t_stamp.min_compute, &t_stamp.min_compute_index);
	minimum(palce_holder_for_timing_max_min[1], itteration, &t_stamp.min_comm, &t_stamp.min_comm_index);
	minimum(tmp_array, write_count, &t_stamp.min_write, &t_stamp.min_write_index);

	checkpoint_coarse_time_profiling_data(grid, &t_stamp, my_IO_domain_info);

}

#endif
