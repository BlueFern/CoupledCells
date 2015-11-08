#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <arkode/arkode.h>            /* prototypes for ARKode fcts., consts. */
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h>  /* def. of type 'realtype' */
#include <arkode/arkode_spgmr.h>
#include <arkode/arkode_bandpre.h>
#include "computelib.h"
#include "gather.h"
#include "writeHDF5.h"

#define EXPLICIT_ONLY 0

extern conductance cpl_cef;
extern SMC_cell** smc;
extern EC_cell** ec;
extern double **sendbuf, **recvbuf;
extern grid_parms grid;
extern time_stamps t_stamp;

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	N_VectorContent_Serial yContent = (N_VectorContent_Serial) y->content;
	N_VectorContent_Serial yDotContent = (N_VectorContent_Serial) ydot->content;
	compute(grid, smc, ec, cpl_cef, t, yContent->data, yDotContent->data);

	return 0;  // Success.
}

static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	N_VectorContent_Serial yContent = (N_VectorContent_Serial) y->content;
	N_VectorContent_Serial yDotContent = (N_VectorContent_Serial) ydot->content;

	compute_implicit(grid, smc, ec, cpl_cef, t, yContent->data, yDotContent->data);

	return 0;  // Success.
}

static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	N_VectorContent_Serial yContent = (N_VectorContent_Serial) y->content;
	N_VectorContent_Serial yDotContent = (N_VectorContent_Serial) ydot->content;

	compute_explicit(grid, smc, ec, cpl_cef, t, yContent->data, yDotContent->data);

	return 0;  // Success.
}

void ark_check_flag(int cflag, char *funcname, int rank, double tnow)
{
	switch (cflag)
	{
	case ARK_SUCCESS:
		break;

	case ARK_ROOT_RETURN:
		break;

	case ARK_TSTOP_RETURN:
		break;

	case ARK_MEM_NULL:
		fprintf(stderr, "[%d] in %s. arkode_mem is NULL at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);

	case ARK_MEM_FAIL:
		fprintf(stderr, "[%d] in %s. Memory allocation failed at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);

	case ARK_NO_MALLOC:
		fprintf(stderr, "[%d] in %s. arkode_mem was not allocated at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);

	case ARK_ILL_INPUT:
		fprintf(stderr, "[%d] in %s. Illegal input/s at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);

	case ARK_TOO_MUCH_WORK:
		fprintf(stderr, "[%d] in %s. Too much work at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		break;

	case ARK_TOO_MUCH_ACC:
		fprintf(stderr, "[%d] in %s. Required accuracy too high at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		break;

	case ARK_ERR_FAILURE:
		fprintf(stderr, "[%d] in %s. Error test failures at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		break;

	case ARK_CONV_FAILURE:
		fprintf(stderr, "[%d] in %s. Convergence test failures at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		break;

	case ARK_LINIT_FAIL:
		fprintf(stderr, "[%d] in %s. Linear solver’s initialization function failed at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);

	case ARK_LSETUP_FAIL:
		fprintf(stderr, "[%d] in %s. Linear solver’s setup routine failed at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);

	case ARK_LSOLVE_FAIL:
		fprintf(stderr, "[%d] in %s. Linear solver’s setup routine failed at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);

	case ARK_MASSINIT_FAIL:
		fprintf(stderr, "[%d] in %s. Mass matrix solver’s initialization function failed at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);

	case ARK_MASSSETUP_FAIL:
		fprintf(stderr, "[%d] in %s. Mass matrix solver’s setup routine failed at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);

	case ARK_MASSSOLVE_FAIL:
		fprintf(stderr, "[%d] in %s. Mass matrix solver’s setup routine failed at tnow = %f: %d\n", rank, funcname, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);

	default:
		fprintf(stderr, "[%d] Unexpected ARK Solver Error at tnow = %f: %d\n", rank, tnow, cflag);
		MPI_Abort(MPI_COMM_WORLD, 911);
	}
}

void arkode_solver(double tnow, double tfinal, double interval, double *yInitial, int neq, double relTOL, double absTOL,
		int file_write_per_unit_time, char* path)
{
	int iteration = 0;
	int write_count = 0;
	initialize_t_stamp(&t_stamp);

	/* general problem variables */
	int flag;                       /* reusable error-checking flag */
	N_Vector y = NULL;              /* empty vector for storing solution */
	void *arkode_mem = NULL;        /* empty ARKode memory structure */
	realtype t, tout;

	y = N_VMake_Serial(neq, yInitial);

	arkode_mem = ARKodeCreate();

#if EXPLICIT_ONLY
	flag = ARKodeInit(arkode_mem, f, NULL, tnow, y);
	ark_check_flag(flag, (char *)"ARKodeInit", grid.universal_rank, 0);

#else
	flag = ARKodeInit(arkode_mem, fe, fi, tnow, y);
	ark_check_flag(flag, (char *)"ARKodeInit", grid.universal_rank, 0);

	// Uses a scaled, preconditioned GMRES (Generalised Minimal Residual) solver without restarts.
	flag = ARKSpgmr(arkode_mem, PREC_RIGHT, 0);
	flag = ARKBandPrecInit(arkode_mem, neq, 1, 1);
#endif

	flag = ARKodeSStolerances(arkode_mem, relTOL, absTOL);
	ark_check_flag(flag, (char *)"ARKodeSVtolerances", grid.universal_rank, 0);

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
			ec[i][j].JPLC = grid.uniform_jplc;
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

	t = tnow;
	tout = tnow + interval;

	// ITERATION loop to go from INITIAL time to FINAL time.
	while (tfinal - t > 1.0e-15)
	{
		double solver_start = MPI_Wtime();

		// Read JPLC in if it is time to do so.
		if(t >= grid.stimulus_onset_time && !jplc_read_in)
		{
			read_init_ATP(&grid, ec);
			jplc_read_in = true;
		}

		flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);  // Call integrator.
		ark_check_flag(flag, (char *)"ARKode", grid.universal_rank, tnow);

		tout += interval;

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

		// checkpoint_timing_data(grid, check, tnow, t_stamp, iteration, file_offset_for_timing_data);
		initialize_t_stamp(&t_stamp);

		/// Increment the iteration as rksuite has finished solving between bounds tnow <= t <= tend.
		iteration++;
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
		tmp_array[i] = palce_holder_for_timing_max_min[2][i];
	}

	// Write time profiling data.
	// Time profiling data gets lost in the event of a crash.
	dump_time_profiling(grid, &t_stamp);
}
