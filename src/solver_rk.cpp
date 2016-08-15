#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>

#include "computelib.h"
#include "gather.h"
#include "writeHDF5.h"

#include "rksuite.h"

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
	compute(grid, smc, ec, cpl_cef, t, y, f);
}

// neq: redundant parameter.
// path: redundant parameter.
void rksuite_solver_CT(double tnow, double tfinal, double interval, double *y, double* yp,
		int neq, double TOL, double* thres, int file_write_per_unit_time, char* path)
{
	RKSUITE rksuite;

	// Solver method.
	int method = 2; // RK(4,5): Forth order evaluation and fifth order correction.
	double tend;
	int cflag = 0;
	int iteration = 0;
	int write_count = 0;
	int count = 0;
	char CTstr[] = "CT";

	tend = interval;

	rksuite.setup(neq, tnow, y, tend, TOL, thres, method, CTstr, false, 0.0, false);

	initialize_t_stamp(&t_stamp);

	// Exchange SMC and EC variables in the ghost cells.
	// Essential for restarts when data is loaded from a checkpoint.
	communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);

	// This is the way it is done in RKSUITE examples to ensure realistic values of y.
	// Do not need this for a restart on initialised run.
	// computeDerivatives(tnow, y, yp);
	MPI_Barrier(grid.universe);

	int totf, stpcst, stpsok;
	double waste, hnext;

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

	// Write JPLC for validation.
	if(grid.rank_branch == 0)
	{
		write_HDF5_JPLC(&grid, jplc_buffer, path);
		free(jplc_buffer);
	}

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

		double solver_start = MPI_Wtime();

		// The ct() function does not guarantee to advance all the way to the stop time. Keep stepping until it does.
		// The solver decides step magnitude depending on the stiffness of the problem.
		do
		{
			rksuite.ct(computeDerivatives, tnow, y, yp, cflag);

			// report_RK_suite_error(cflag, tnow, grid.universal_rank);
		}
		while(tnow < tend);

		/// rksuite.stat() routine gathers the statistics on the performance of the solver.
		/// Amongst other things it also informs about what the next step size should be.
		rksuite.stat(totf, stpcst, waste, stpsok, hnext);

		t_stamp.aggregate_compute += MPI_Wtime() - solver_start;

		/// Call for interprocessor communication.
		double comms_0_start = MPI_Wtime();
		communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);
		t_stamp.aggregate_comm += MPI_Wtime() - comms_0_start;

		if((iteration % file_write_per_unit_time) == 0)
		{
			// Collect state variable data on writers.

			double ec_gather = MPI_Wtime();
			gather_EC_data(&grid, ec_buffer, ec);
			t_stamp.aggregate_ec_gather += MPI_Wtime() - ec_gather;

			double smc_gather = MPI_Wtime();
			gather_SMC_data(&grid, smc_buffer, smc);
			t_stamp.aggregate_smc_gather += MPI_Wtime() - smc_gather;

			double ec_write = MPI_Wtime();
			if(grid.rank_branch == 0)
			{
				// Write state variables to HDF5.
				write_EC_data_HDF5(&grid, ec_buffer, write_count, path);
			}
			t_stamp.aggregate_ec_write += MPI_Wtime() - ec_write;

			double smc_write = MPI_Wtime();
			if(grid.rank_branch == 0)
			{
				// Write state variables to HDF5.
				write_SMC_data_HDF5(&grid, smc_buffer, write_count, path);
			}
			t_stamp.aggregate_smc_write += MPI_Wtime() - smc_write;

			write_count++;
		}

		// Read JPLC in if it is time to do so.
		if(tnow >= grid.stimulus_onset_time && !jplc_read_in)
		{
			read_lumanel_values(&grid, ec);
			jplc_read_in = true;
		}

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

	// Write time profiling data.
	// Time profiling data is lost in the event of a crash.
	dump_time_profiling(grid, &t_stamp);
}
