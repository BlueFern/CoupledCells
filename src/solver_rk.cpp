#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>

#include "computelib.h"
#include "gather.h"
#include "writeHDF5.h"

#include "rksuite.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define SRC_LOC __FILE__ ":" TOSTRING(__LINE__)

extern conductance cpl_cef;
extern SMC_cell** smc;
extern EC_cell** ec;
extern double **sendbuf, **recvbuf;
extern grid_parms grid;
extern time_stamps t_stamp;

char RK_suite_error[] = "[%d] RK solver return value at tnow = %f: %d = \"%s\".\n";

void report_RK_suite_error(int cflag, double tnow, int rank)
{
	switch(cflag)
	{
	// 1 - success.
	case 1:
		break;
	// 2 - inefficient usage.
	case 2:
		fprintf(stderr, RK_suite_error, rank, tnow, cflag, "inefficient usage");
		break;
	// 3 - work intensive.
	case 3:
		fprintf(stderr, RK_suite_error, rank, tnow, cflag, "work intensive");
		break;
	// 4 - problem probably stiff.
	case 4:
		fprintf(stderr, RK_suite_error, rank, tnow, cflag, "problem probably stiff");
		break;
	// 5 - too much accuracy requested.
	case 5:
		fprintf(stderr, RK_suite_error, rank, tnow, cflag, "too much accuracy requested");
		break;
	// 6 - global error assessment unreliable beyond this point.
	case 6:
		fprintf(stderr, RK_suite_error, rank, tnow, cflag, "global error assessment unreliable beyond this point");
		break;
	// 911 - catastrophic failure reported on stdout.
	case 911:
		fprintf(stderr, RK_suite_error, rank, tnow, cflag, "catastrophic failure");
		MPI_Abort(MPI_COMM_WORLD, 911);
		break;
	default:
		fprintf(stderr, "[%d] Unexpected RK Solver Error at tnow = %f: %d\n", rank, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);
	}
}

void computeDerivatives(double t, double y[], double f[])
{
	// compute_with_time_profiling(&t_stamp, grid, smc, ec, cpl_cef, t, y, f);
	compute(grid, smc, ec, cpl_cef, t, y, f);
	t_stamp.computeDerivatives_call_counter += 1;

}

// line_number: redundant parameter.
void rksuite_solver_CT(double tnow, double tfinal, double interval, double *y, double* yp, int neq, double TOL, double* thres,
		int file_write_per_unit_time, char* path) //, IO_domain_info *my_IO_domain_info)
{
	RKSUITE rksuite;

	// Solver method.
	int method = 2; // RK(4,5): Forth order evaluation and fifth order correction.
	double tend;
	int cflag = 0;
	int iteration = 0;
	int write_count = 0;
	int count = 0;
	initialize_t_stamp(&t_stamp);
	tend = interval;
	char CTstr[] = "CT";
	rksuite.setup(neq, tnow, y, tend, TOL, thres, method, CTstr, false, 0.0, false);

	// Exchange SMC and EC variables in the ghost cells.
	// Essential for restarts when data is loaded from a checkpoint.
	communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);

	// This is the way it is done in RKSUITE examples to ensure realistic values of y.
	// Do not need this for a restart on initialised run.
	// computeDerivatives(tnow, y, yp);
	MPI_Barrier(grid.universe);

	// int file_offset_for_timing_data = determine_file_offset_for_timing_data(check, grid);
	int totf, stpcst, stpsok;
	double waste, hnext;

	// Start HDF5 Output.
	// Buffer for the JPLC values for the current branch.
	double *jplc_buffer = 0;

	// Allocate the jplc_buffer.
	if (grid.rank_branch == 0)
	{
		jplc_buffer = (double *)checked_malloc(grid.num_ranks_branch * grid.num_ec_axially * grid.num_ec_circumferentially * sizeof(double), SRC_LOC);
	}

	// Collect all JPLC values in a single buffer on root node.
	gather_JPLC(&grid, jplc_buffer, ec);

	// MPI_Barrier(MPI_COMM_WORLD);

	// Write jplc values to HDF5.
	if(grid.rank_branch == 0)
	{
		write_HDF5_JPLC(&grid, jplc_buffer, path);
		free(jplc_buffer);
	}
	// End HDf5 Output.

	// Reset JPLC to the uniform map.
	// The input file will have to be read later when the time is right.
	for (int i = 1; i <= grid.num_ec_circumferentially; i++)
	{
		for (int j = 1; j <= grid.num_ec_axially; j++)
		{
			ec[i][j].JPLC = grid.uniform_jplc;
		}
	}

	bool jplc_read_in = false;

	// Profiling.
	double timing_values[3][int(tfinal / interval)];

	// Data buffers for collecting and writing HDF5 output.
	ec_data_buffer *ec_buffer;
	smc_data_buffer *smc_buffer;

	// Allocate the EC and SMC buffers.
	if(grid.rank_branch == 0)
	{
		ec_buffer = allocate_EC_data_buffer(grid.num_ranks_branch, grid.num_ec_axially * grid.num_ec_circumferentially, 1);
		smc_buffer = allocate_SMC_data_buffer(grid.num_ranks_branch, grid.num_smc_axially * grid.num_smc_circumferentially, 1);
	}
	else
	{
		ec_buffer = allocate_EC_data_buffer(grid.num_ranks_branch, grid.num_ec_axially * grid.num_ec_circumferentially, 0);
		smc_buffer = allocate_SMC_data_buffer(grid.num_ranks_branch, grid.num_smc_axially * grid.num_smc_circumferentially, 0);
	}

	// ITERATION loop to go from INITIAL time to FINAL time.
	while(tnow <= tfinal)
	{
		t_stamp.solver_t1 = MPI_Wtime();

		// printf("[%d] line: %d, tnow: %f, %f, %f\n", grid.universal_rank, __LINE__, tnow, tend, tfinal);

		// The ct() function does not guarantee to advance all the way to the stop time. Keep stepping until it does.
		// The solver decides step magnitude depending on the stiffness of the problem.
		do
		{
			// Read JPLC in if it is time to do so.
			if(tnow >= grid.stimulus_onset_time && !jplc_read_in)
			{
				read_init_ATP(&grid, ec);
				jplc_read_in = true;
			}

			rksuite.ct(computeDerivatives, tnow, y, yp, cflag);

			// report_RK_suite_error(cflag, tnow, grid.universal_rank);

			// printf("\t[%d] line: %d, tnow: %f, %f, %f\n", grid.universal_rank, __LINE__, tnow, tend, tfinal);

		}
		while(tnow < tend);

		/// rksuite.stat() routine gathers the statistics on the performance of the solver.
		/// Amongst other things it also informs about what the next step size should be.
		rksuite.stat(totf, stpcst, waste, stpsok, hnext);

		t_stamp.solver_t2 = MPI_Wtime();
		t_stamp.diff_solver = t_stamp.solver_t2 - t_stamp.solver_t1;

		timing_values[0][iteration] = t_stamp.diff_solver;
		t_stamp.aggregate_compute += t_stamp.diff_solver;

		/// Call for interprocessor communication.
		t_stamp.total_comms_cost_t1 = MPI_Wtime();

		communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);

		t_stamp.total_comms_cost_t2 = MPI_Wtime();
		t_stamp.diff_total_comms_cost = t_stamp.total_comms_cost_t2 - t_stamp.total_comms_cost_t1;
		timing_values[1][iteration] = t_stamp.diff_total_comms_cost;
		t_stamp.aggregate_comm += t_stamp.diff_total_comms_cost;

		if((iteration % file_write_per_unit_time) == 0)
		{
			t_stamp.write_t1 = MPI_Wtime();

			// HDF5 Start
			// HDF5 Start
			// HDF5 Start

			// Collect state variable data on writers.
			// TODO: Perhaps the ECs matrix is better stored as pointers in the grid?
			gather_EC_data(&grid, ec_buffer, ec);
			// TODO: Perhaps the SMCs matrix is better stored as pointers in the grid?
			gather_SMC_data(&grid, smc_buffer, smc);

			if(grid.rank_branch == 0)
			{
				// Write state variables to HDF5.
				write_EC_data_HDF5(&grid, ec_buffer, write_count, path);
				write_SMC_data_HDF5(&grid, smc_buffer, write_count, path);
			}

			// HDF5 End
			// HDF5 End
			// HDF5 End

			t_stamp.write_t2 = MPI_Wtime();
			t_stamp.diff_write = t_stamp.write_t2 - t_stamp.write_t1;
			timing_values[2][write_count] = t_stamp.diff_write;
			t_stamp.aggregate_write += t_stamp.diff_write;

			write_count++;
		}

		// checkpoint_timing_data(grid, check, tnow, t_stamp, iteration, file_offset_for_timing_data);
		initialize_t_stamp(&t_stamp);

		/// Increment the iteration as rksuite has finished solving between bounds tnow <= t <= tend.
		iteration++;
		tend += interval;
		rksuite.reset(tend);
	}

	// Release the EC and SMC buffers used for writing to HDF5.
	if(grid.rank_branch == 0)
	{
		free_EC_data_buffer(ec_buffer, 1);
		free_SMC_data_buffer(smc_buffer, 1);
	}
	else
	{
		free_EC_data_buffer(ec_buffer, 0);
		free_SMC_data_buffer(smc_buffer, 0);
	}

	// Prepare time profiling data.
	double tmp_array[write_count];

	for(int i = 0; i < write_count; i++)
	{
		tmp_array[i] = timing_values[2][i];
	}
	maximum(timing_values[0], iteration, &t_stamp.max_compute, &t_stamp.max_compute_index);
	maximum(timing_values[1], iteration, &t_stamp.max_comm, &t_stamp.max_comm_index);
	maximum(tmp_array, write_count, &t_stamp.max_write, &t_stamp.max_write_index);

	minimum(timing_values[0], iteration, &t_stamp.min_compute, &t_stamp.min_compute_index);
	minimum(timing_values[1], iteration, &t_stamp.min_comm, &t_stamp.min_comm_index);
	minimum(tmp_array, write_count, &t_stamp.min_write, &t_stamp.min_write_index);

	// Write time profiling data.
	// Time profiling data is lost in the event of a crash.
	checkpoint_coarse_time_profiling_data(grid, &t_stamp); //, my_IO_domain_info);
}
