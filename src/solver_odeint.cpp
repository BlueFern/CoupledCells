#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <boost/numeric/odeint.hpp>
#include <boost/range.hpp>
#include <list>

#include "computelib.h"
#include "gather.h"
#include "writeHDF5.h"

using namespace boost::numeric::odeint;
using namespace boost;

extern conductance cpl_cef;
extern SMC_cell** smc;
extern EC_cell** ec;
extern double **sendbuf, **recvbuf;
extern grid_parms grid;
extern time_stamps t_stamp;

typedef std::vector<double> state_type;
typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

bool jplc_read_in = false;

typedef iterator_range<std::vector<double>::iterator> time_range;

void f(const state_type &y, state_type &f, const double t)
{
	compute(grid, smc, ec, cpl_cef, t, const_cast<double*> (&y[0]), &f[0]);
	t_stamp.computeDerivatives_call_counter += 1;
}

struct observer
{
    int* iteration_;
    int* file_write_per_unit_time_;
    IO_domain_info *my_IO_domain_info_;
    double* palce_holder_for_timing_max_min_;
    // data_buffer *writer_buffer_;
    int* write_count_;
    int sizee_;
    checkpoint_handle *check_;
    ec_data_buffer *ec_buffer_;
	smc_data_buffer *smc_buffer_;
	char* path_;

    observer(double* palce_holder_for_timing_max_min, int* iteration, int* file_write_per_unit_time, IO_domain_info *my_IO_domain_info,
    		// data_buffer* writer_buffer,
			int* write_count, checkpoint_handle *check, ec_data_buffer *ec_buffer, smc_data_buffer *smc_buffer,
			char* path, int sizee) :
			palce_holder_for_timing_max_min_(palce_holder_for_timing_max_min) , iteration_(iteration),
			file_write_per_unit_time_(file_write_per_unit_time), my_IO_domain_info_ (my_IO_domain_info),
			// writer_buffer_ (writer_buffer),
			write_count_ (write_count), check_(check), ec_buffer_(ec_buffer),
			smc_buffer_ (smc_buffer), path_(path), sizee_ (sizee)
    { }

    void operator()( const state_type &y , const double tnow )
    {

    	// Read JPLC in if it is time to do so.
		if(tnow >= grid.stimulus_onset_time && !jplc_read_in)
		{
			read_init_ATP(&grid, ec);
			jplc_read_in = true;
		}

    	t_stamp.solver_t2 = MPI_Wtime();
		t_stamp.diff_solver = t_stamp.solver_t2 - t_stamp.solver_t1;

		palce_holder_for_timing_max_min_[0 * sizee_ + *iteration_] = t_stamp.diff_solver;
		t_stamp.aggregate_compute += t_stamp.diff_solver;

		/// Call for interprocessor communication.
		t_stamp.total_comms_cost_t1 = MPI_Wtime();

		communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);

		t_stamp.total_comms_cost_t2 = MPI_Wtime();
		t_stamp.diff_total_comms_cost = t_stamp.total_comms_cost_t2 - t_stamp.total_comms_cost_t1;
		palce_holder_for_timing_max_min_[1 * sizee_ + *iteration_] = t_stamp.diff_total_comms_cost;
		t_stamp.aggregate_comm += t_stamp.diff_total_comms_cost;

		if(((*iteration_) % (*file_write_per_unit_time_)) == 0)
		{

			t_stamp.write_t1 = MPI_Wtime();

#if 0
			// Geometry to be written.
			gather_smc_mesh_data_on_writers(&grid, my_IO_domain_info_, writer_buffer_, smc);
			gather_ec_mesh_data_on_writers(&grid, my_IO_domain_info_, writer_buffer_, ec);

			// State variables to be written as attributes.
			gather_smcData(&grid, my_IO_domain_info_, writer_buffer_, smc, *write_count_);
			gather_ecData(&grid, my_IO_domain_info_, writer_buffer_, ec, *write_count_);

			if (grid.rank == 0)
			{
				// Write out the meshes with attributes.
				initialise_time_wise_checkpoint(check_, grid, *write_count_, path_, my_IO_domain_info_);
				write_smc_and_ec_data(check_, &grid, tnow, smc, ec, *write_count_, my_IO_domain_info_, writer_buffer_);
				close_time_wise_checkpoints(check_);
			}
#endif
			// HDF5 Start
			// HDF5 Start
			// HDF5 Start

			// Collect state variable data on writers.
			// TODO: Perhaps the ECs matrix is better stored as pointers in the grid?
			gather_EC_data(&grid, ec_buffer_, ec);
			// TODO: Perhaps the SMCs matrix is better stored as pointers in the grid?
			gather_SMC_data(&grid, smc_buffer_, smc);

			if(grid.rank == 0)
			{
				// Write state variables to HDF5.
				write_EC_data_HDF5(&grid, ec_buffer_, *write_count_, path_);
				write_SMC_data_HDF5(&grid, smc_buffer_, *write_count_, path_);
			}

			// HDF5 End
			// HDF5 End
			// HDF5 End

			t_stamp.write_t2 = MPI_Wtime();
			t_stamp.diff_write = t_stamp.write_t2 - t_stamp.write_t1;
			palce_holder_for_timing_max_min_[2 * sizee_ + *write_count_] = t_stamp.diff_write;
			t_stamp.aggregate_write += t_stamp.diff_write;
			(*write_count_)++;
		}

		// checkpoint_timing_data(grid, check, tnow, t_stamp, iteration, file_offset_for_timing_data);
		initialize_t_stamp(&t_stamp);
		/// Increment the iteration as rksuite has finished solving between bounds tnow <= t <= tend.
		(*iteration_)++;
    }
};


void odeint_solver(double tnow, double tfinal, double interval, double *yInitial, int neq, double relTol, double absTol,
		int file_write_per_unit_time, checkpoint_handle *check, char* path, IO_domain_info *my_IO_domain_info)
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

	// This is the way it is done in RKSUITE examples to ensure realistic values of y.
	// Do not need this for a restart on initialised run.
	// computeDerivatives(tnow, y, yp);
	MPI_Barrier(grid.universe);

	// int file_offset_for_timing_data = determine_file_offset_for_timing_data(check, grid);
	int totf, stpcst, stpsok;
	double waste, hnext;
	int err;

#if 0
	// Aggregation of data for the writer.
	data_buffer* writer_buffer = (data_buffer*) checked_malloc(sizeof(data_buffer), SRC_LOC);
	writer_buffer->buffer_length = (int*) checked_malloc(12 * sizeof(int), SRC_LOC);
	writer_buffer->smc_stat_var_buffer_length = (int*) checked_malloc((grid.neq_smc) * sizeof(int), SRC_LOC);
	writer_buffer->ec_stat_var_buffer_length = (int*) checked_malloc((grid.neq_ec) * sizeof(int), SRC_LOC);
	writer_buffer->smc_cpl = (int*) checked_malloc((grid.num_coupling_species_smc) * sizeof(int), SRC_LOC);
	writer_buffer->ec_cpl = (int*) checked_malloc((grid.num_coupling_species_ec) * sizeof(int), SRC_LOC);

	// Dump MPI task mesh representation into vtk file to manifest task map.
	gather_tasks_mesh_point_data_on_writers(&grid, my_IO_domain_info, writer_buffer, smc, ec);
	if (grid.rank == 0)
	{
		// Validation, debugging.
		write_process_mesh(check, &grid, my_IO_domain_info, writer_buffer, path);
	}

	// Dump JPLC map on bifurcation into a vtk file.
	gather_ec_mesh_data_on_writers(&grid, my_IO_domain_info, writer_buffer, ec);
	gather_JPLC_map(&grid, my_IO_domain_info, writer_buffer, ec);

	if(grid.rank == 0)
	{
		// Initial concentration of JPLC in the EC cells.
		write_JPLC_map(check, &grid, my_IO_domain_info, writer_buffer, ec, path);
	}
#endif

	// Start HDF5 Output Prototyping.

	// Buffer for jplc values for the whole mesh.
	double *jplc_buffer = 0;

	// Allocate the jplc_buffer.
	if (grid.rank == 0)
	{
		jplc_buffer = (double *)checked_malloc(grid.tasks * grid.num_ec_axially * grid.num_ec_circumferentially * sizeof(double), SRC_LOC);
	}

	// Collect all jplc values in a single buffer on root node.
	gather_JPLC(&grid, jplc_buffer, ec);

	MPI_Barrier(MPI_COMM_WORLD);

	// Write jplc values to HDF5.
	if(grid.rank == 0)
	{
		write_HDF5_JPLC(&grid, jplc_buffer, path);
		free(jplc_buffer);
	}

	// End HDf5 Output Prototyping.

	// Reset JPLC to the uniform map.
	// The input file will have to be read later when the time is right.
	for (int i = 1; i <= grid.num_ec_circumferentially; i++)
	{
		for (int j = 1; j <= grid.num_ec_axially; j++)
		{
			ec[i][j].JPLC = grid.uniform_jplc; // agonist_profile((grid.stimulus_onset_time + 1), grid, i, j, ec[i][j].centeroid_point[1]);
		}
	}

	// Profiling.
	double palce_holder_for_timing_max_min[3][int(tfinal / interval)];

	// Data buffers for collecting and writing HDF5 output.
	ec_data_buffer *ec_buffer;
	smc_data_buffer *smc_buffer;

	// Allocate the EC and SMC buffers.
	if(grid.rank == 0)
	{
		ec_buffer = allocate_EC_data_buffer(grid.tasks, grid.num_ec_axially * grid.num_ec_circumferentially, 1);
		smc_buffer = allocate_SMC_data_buffer(grid.tasks, grid.num_smc_axially * grid.num_smc_circumferentially, 1);
	}
	else
	{
		ec_buffer = allocate_EC_data_buffer(grid.tasks, grid.num_ec_axially * grid.num_ec_circumferentially, 0);
		smc_buffer = allocate_SMC_data_buffer(grid.tasks, grid.num_smc_axially * grid.num_smc_circumferentially, 0);
	}


	int steps = (int) (((tfinal - tnow) / interval) + 1);

	std::vector<double> time_range;
	for (int i = 0; i < steps; i ++)
	{
		time_range.push_back((double) i * interval);
	}

	integrate_times(make_controlled< error_stepper_type >( absTol , relTol), f , y , time_range.begin() ,time_range.end(), interval,
			observer(&palce_holder_for_timing_max_min[0][0], &iteration, &file_write_per_unit_time, my_IO_domain_info,
			&write_count, check, ec_buffer, smc_buffer, path, int(tfinal / interval)));

	// Release the EC and SMC buffers used for writing to HDF5.
	if(grid.rank == 0)
	{
		free_EC_data_buffer(ec_buffer, 1);
		free_SMC_data_buffer(smc_buffer, 1);
	}
	else
	{
		free_EC_data_buffer(ec_buffer, 0);
		free_SMC_data_buffer(smc_buffer, 0);
	}

	//t_stamp.aggregate_compute = t_stamp.aggregate_compute / iteration;
	//t_stamp.aggregate_comm = t_stamp.aggregate_comm / iteration;
	//t_stamp.aggregate_write = t_stamp.aggregate_write / write_count;

	// Prepare time profiling data.
	double tmp_array[write_count];

	for(int i = 0; i < write_count; i++)
	{
		tmp_array[i] = palce_holder_for_timing_max_min[2][i];
	}
	//maximum(palce_holder_for_timing_max_min[0], iteration, &t_stamp.max_compute, &t_stamp.max_compute_index);
	//maximum(palce_holder_for_timing_max_min[1], iteration, &t_stamp.max_comm, &t_stamp.max_comm_index);
	//maximum(tmp_array, write_count, &t_stamp.max_write, &t_stamp.max_write_index);

	//minimum(palce_holder_for_timing_max_min[0], iteration, &t_stamp.min_compute, &t_stamp.min_compute_index);
	//minimum(palce_holder_for_timing_max_min[1], iteration, &t_stamp.min_comm, &t_stamp.min_comm_index);
	//minimum(tmp_array, write_count, &t_stamp.min_write, &t_stamp.min_write_index);

	// Write time profiling data.
	// Time profiling data gets lost in the event of a crash.
	checkpoint_coarse_time_profiling_data(grid, &t_stamp, my_IO_domain_info);
}

// **** Simpler version without observer. ****
#if 0
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/array.hpp>

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
	//printf("%f\n",(const_cast<double*> (&y[0]))[0]);
	// compute_with_time_profiling(&t_stamp, grid, smc, ec, cpl_cef, t, y, f);
	compute(grid, smc, ec, cpl_cef, t, const_cast<double*> (&y[0]), &f[0]);

	t_stamp.computeDerivatives_call_counter += 1;
}

void odeint_solver(double tnow, double tfinal, double interval, double *yInitial, int neq, double relTol, double absTol,
		int file_write_per_unit_time, checkpoint_handle *check, char* path, IO_domain_info *my_IO_domain_info)
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

	// This is the way it is done in RKSUITE examples to ensure realistic values of y.
	// Do not need this for a restart on initialised run.
	// computeDerivatives(tnow, y, yp);
	MPI_Barrier(grid.universe);

	// int file_offset_for_timing_data = determine_file_offset_for_timing_data(check, grid);
	int totf, stpcst, stpsok;
	double waste, hnext;
	int err;

	// Aggregation of data for the writer.
	data_buffer* writer_buffer = (data_buffer*) checked_malloc(sizeof(data_buffer), SRC_LOC);
	writer_buffer->buffer_length = (int*) checked_malloc(12 * sizeof(int), SRC_LOC);
	writer_buffer->smc_stat_var_buffer_length = (int*) checked_malloc((grid.neq_smc) * sizeof(int), SRC_LOC);
	writer_buffer->ec_stat_var_buffer_length = (int*) checked_malloc((grid.neq_ec) * sizeof(int), SRC_LOC);
	writer_buffer->smc_cpl = (int*) checked_malloc((grid.num_coupling_species_smc) * sizeof(int), SRC_LOC);
	writer_buffer->ec_cpl = (int*) checked_malloc((grid.num_coupling_species_ec) * sizeof(int), SRC_LOC);

	// Dump MPI task mesh representation into vtk file to manifest task map.
	gather_tasks_mesh_point_data_on_writers(&grid, my_IO_domain_info, writer_buffer, smc, ec);
	if (grid.rank == 0)
	{
		// Validation, debugging.
		write_process_mesh(check, &grid, my_IO_domain_info, writer_buffer, path);
	}

	// Dump JPLC map on bifurcation into a vtk file.
	gather_ec_mesh_data_on_writers(&grid, my_IO_domain_info, writer_buffer, ec);
	gather_JPLC_map(&grid, my_IO_domain_info, writer_buffer, ec);

	if(grid.rank == 0)
	{
		// Initial concentration of JPLC in the EC cells.
		write_JPLC_map(check, &grid, my_IO_domain_info, writer_buffer, ec, path);
	}

	// Start HDF5 Output Prototyping.
	// printf("* %d\t%d\t%d\t%d\t%d *\n", grid.universal_rank, grid.sub_universe_numtasks, grid.sub_universe_rank, grid.rank, grid.tasks);

	// Buffer for jplc values for the whole mesh.
	double *jplc_buffer = 0;

	// Allocate the jplc_buffer.
	if (grid.rank == 0)
	{
		jplc_buffer = (double *)checked_malloc(grid.tasks * grid.num_ec_axially * grid.num_ec_circumferentially * sizeof(double), SRC_LOC);
	}

	// Collect all jplc values in a single buffer on root node.
	gather_JPLC(&grid, jplc_buffer, ec);

	MPI_Barrier(MPI_COMM_WORLD);

	// Write jplc values to HDF5.
	if(grid.rank == 0)
	{
		write_HDF5_JPLC(&grid, jplc_buffer, path);
		free(jplc_buffer);
	}

	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("[%d] *** Pulling the plug in %s:%s\n", grid.universal_rank, __FILE__, __FUNCTION__);
	//MPI_Abort(MPI_COMM_WORLD, 911);

	// End HDf5 Output Prototyping.

	// printf("%s, grid.cart_comm: %p\n", __FUNCTION__, (void *)grid.cart_comm);

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
	if(grid.rank == 0)
	{
		ec_buffer = allocate_EC_data_buffer(grid.tasks, grid.num_ec_axially * grid.num_ec_circumferentially, 1);
		smc_buffer = allocate_SMC_data_buffer(grid.tasks, grid.num_smc_axially * grid.num_smc_circumferentially, 1);
	}
	else
	{
		ec_buffer = allocate_EC_data_buffer(grid.tasks, grid.num_ec_axially * grid.num_ec_circumferentially, 0);
		smc_buffer = allocate_SMC_data_buffer(grid.tasks, grid.num_smc_axially * grid.num_smc_circumferentially, 0);
	}

	while(tnow <= tfinal)
	{
		t_stamp.solver_t1 = MPI_Wtime();

		// Read JPLC in if it is time to do so.
		if(tnow >= grid.stimulus_onset_time && !jplc_read_in)
		{
			read_init_ATP(&grid, ec);
			jplc_read_in = true;
		}

		integrate_adaptive(make_controlled<error_stepper_type>( absTol , relTol), f , y , tnow , tend , interval);
		//integrate_const(stepper, f, y , tnow , tend , interval);

		t_stamp.solver_t2 = MPI_Wtime();
		t_stamp.diff_solver = t_stamp.solver_t2 - t_stamp.solver_t1;

		palce_holder_for_timing_max_min[0][iteration] = t_stamp.diff_solver;
		t_stamp.aggregate_compute += t_stamp.diff_solver;

		/// Call for interprocessor communication.
		t_stamp.total_comms_cost_t1 = MPI_Wtime();

		communication_async_send_recv(grid, sendbuf, recvbuf, smc, ec);

		t_stamp.total_comms_cost_t2 = MPI_Wtime();
		t_stamp.diff_total_comms_cost = t_stamp.total_comms_cost_t2 - t_stamp.total_comms_cost_t1;
		palce_holder_for_timing_max_min[1][iteration] = t_stamp.diff_total_comms_cost;
		t_stamp.aggregate_comm += t_stamp.diff_total_comms_cost;

		if((iteration % file_write_per_unit_time) == 0)
		{
			t_stamp.write_t1 = MPI_Wtime();

			// Geometry to be written.
			gather_smc_mesh_data_on_writers(&grid, my_IO_domain_info, writer_buffer, smc);
			gather_ec_mesh_data_on_writers(&grid, my_IO_domain_info, writer_buffer, ec);

			// State variables to be written as attributes.
			gather_smcData(&grid, my_IO_domain_info, writer_buffer, smc, write_count);
			gather_ecData(&grid, my_IO_domain_info, writer_buffer, ec, write_count);

			if (grid.rank == 0)
			{
				// Write out the meshes with attributes.
				initialise_time_wise_checkpoint(check, grid, write_count, path, my_IO_domain_info);
				write_smc_and_ec_data(check, &grid, tnow, smc, ec, write_count, my_IO_domain_info, writer_buffer);
				close_time_wise_checkpoints(check);
			}

			// HDF5 Start
			// HDF5 Start
			// HDF5 Start

			// Collect state variable data on writers.
			// TODO: Perhaps the ECs matrix is better stored as pointers in the grid?
			gather_EC_data(&grid, ec_buffer, ec);
			// TODO: Perhaps the SMCs matrix is better stored as pointers in the grid?
			gather_SMC_data(&grid, smc_buffer, smc);

			if(grid.rank == 0)
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
			palce_holder_for_timing_max_min[2][write_count] = t_stamp.diff_write;
			t_stamp.aggregate_write += t_stamp.diff_write;
			write_count++;
		}

		// checkpoint_timing_data(grid, check, tnow, t_stamp, iteration, file_offset_for_timing_data);
		initialize_t_stamp(&t_stamp);
		/// Increment the iteration as rksuite has finished solving between bounds tnow <= t <= tend.
		iteration++;
		tnow = tend;
		tend += interval;

	}

	// Release the EC and SMC buffers used for writing to HDF5.
	if(grid.rank == 0)
	{
		free_EC_data_buffer(ec_buffer, 1);
		free_SMC_data_buffer(smc_buffer, 1);
	}
	else
	{
		free_EC_data_buffer(ec_buffer, 0);
		free_SMC_data_buffer(smc_buffer, 0);
	}

	//t_stamp.aggregate_compute = t_stamp.aggregate_compute / iteration;
	//t_stamp.aggregate_comm = t_stamp.aggregate_comm / iteration;
	//t_stamp.aggregate_write = t_stamp.aggregate_write / write_count;

	// Prepare time profiling data.
	double tmp_array[write_count];

	for(int i = 0; i < write_count; i++)
	{
		tmp_array[i] = palce_holder_for_timing_max_min[2][i];
	}
	maximum(palce_holder_for_timing_max_min[0], iteration, &t_stamp.max_compute, &t_stamp.max_compute_index);
	maximum(palce_holder_for_timing_max_min[1], iteration, &t_stamp.max_comm, &t_stamp.max_comm_index);
	maximum(tmp_array, write_count, &t_stamp.max_write, &t_stamp.max_write_index);

	minimum(palce_holder_for_timing_max_min[0], iteration, &t_stamp.min_compute, &t_stamp.min_compute_index);
	minimum(palce_holder_for_timing_max_min[1], iteration, &t_stamp.min_comm, &t_stamp.min_comm_index);
	minimum(tmp_array, write_count, &t_stamp.min_write, &t_stamp.min_write_index);

	// Write time profiling data.
	// Time profiling data gets lost in the event of a crash.
	checkpoint_coarse_time_profiling_data(grid, &t_stamp, my_IO_domain_info);
}
#endif
