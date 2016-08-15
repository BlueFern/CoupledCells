#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>

#include "computelib.h"
#include "gather.h"
#include "writeHDF5.h"

using namespace boost::numeric::odeint;

extern conductance cpl_cef;
extern SMC_cell** smc;
extern EC_cell** ec;
extern double **sendbuf, **recvbuf;
extern grid_parms grid;
extern time_stamps t_stamp;

typedef std::vector< double > state_type;
typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;

void f(const state_type &y, state_type &f, const double t)
{
	compute(grid, smc, ec, cpl_cef, t, const_cast<double*> (&y[0]), &f[0]);
}

void odeint_solver(double tnow, double tfinal, double interval, double *yInitial, int neq, double relTol, double absTol,
		int file_write_per_unit_time, char* path)
{

	double tend = interval;
	int cflag = 0;
	int iteration = 0;
	int write_count = 0;
	int count = 0;
	initialize_t_stamp(&t_stamp);

	state_type y(yInitial, yInitial + neq); // Convert (copy) double* yinitial to std::vector y.

	// Exchange SMC and EC variables in the ghost cells.
	// Essential for restarts when data is loaded from a checkpoint.
	communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);

	MPI_Barrier(grid.universe);

	// Buffer for jplc values for the whole mesh.
	double *jplc_buffer = 0;

	// Allocate the jplc_buffer.
	if (grid.rank_branch == 0)
	{
		jplc_buffer = (double *)checked_malloc(grid.num_ranks_branch * grid.num_ec_axially * grid.num_ec_circumferentially * sizeof(double), SRC_LOC);
	}

	// Collect all jplc values in a single buffer on root node.
	gather_JPLC(&grid, jplc_buffer, ec);

	MPI_Barrier(MPI_COMM_WORLD);

	// Write jplc values to HDF5.
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
			ec[i][j].JPLC = grid.uniform_jplc; // agonist_profile((grid.stimulus_onset_time + 1), grid, i, j, ec[i][j].centeroid_point[1]);
		}
	}

	bool jplc_read_in = false;

	// Profiling.
	double palce_holder_for_timing_max_min[3][int(tfinal / interval)];

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

	while(tnow <= tfinal)
	{
		double solver_start = MPI_Wtime();

		// Read JPLC in if it is time to do so.
		if(tnow >= grid.stimulus_onset_time && !jplc_read_in)
		{
			read_lumanel_values(&grid, ec);
			jplc_read_in = true;
		}

		integrate_adaptive(make_controlled<error_stepper_type>( absTol , relTol), f , y , tnow , tend , interval);

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

		iteration++;
		tnow = tend;
		tend += interval;
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
