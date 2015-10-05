/*
 * This file declares the functions which, when given a group of processors, make either
 * a bifurcation or a straight segment, and recognise its remote and local neighbours.
 * These function are called in main.cpp.
 */

#include <assert.h>
#include <malloc.h>
#include <mpi.h>

#include "computelib.h"

void set_task_parameters(grid_parms *grid)
{
	// Sanity check: make sure the number of requested cores equals the number of cores required by the domain configuration parameters.

	// If the subdomain type is a straight segment.
	if(grid->domain_params[0][DOMAIN_TYPE] == STRSEG)
	{
		assert(grid->num_ranks == grid->domain_params[0][AX_QUADS] * grid->domain_params[0][CR_QUADS]); // Total number of tasks mapping to the straight segment are m x n.
	}
	// If the subdomain is a part of a bifurcation.
	else if(grid->domain_params[0][DOMAIN_TYPE] == BIF)
	{
		assert(grid->num_ranks == 3 * grid->domain_params[0][AX_QUADS] * grid->domain_params[0][CR_QUADS]); // Total number of tasks mapping to the bifurcation are 3 x m x n.
	}

	grid->m = grid->domain_params[0][AX_QUADS];
	grid->n = grid->domain_params[0][CR_QUADS];

	grid->domain_index = grid->domain_params[0][DOMAIN_NUM];
	grid->domain_type = grid->domain_params[0][DOMAIN_TYPE];

	// Each tasks now calculates the number of ECs per node.
	// Each tasks now calculates the number of SMCs per node.
	// Topological information of a functional block of coupled cells.
	// This is the minimum required to simulate a relevant coupled topology.

	grid->num_ec_axially = grid->domain_params[0][AX_ECS] * grid->num_ec_fundblk_axially;
	grid->num_smc_axially = grid->num_ec_axially * grid->num_smc_fundblk_axially;
	grid->num_smc_circumferentially = grid->domain_params[0][CR_SMCS] * grid->num_smc_fundblk_circumferentially;
	grid->num_ec_circumferentially = grid->num_smc_circumferentially * grid->num_ec_fundblk_circumferentially;

	grid->neq_ec_axially = grid->num_ec_axially * grid->neq_ec;
	grid->neq_smc_axially = grid->num_smc_axially * grid->neq_smc;

	grid->NEQ = grid->neq_smc * (grid->num_smc_axially * grid->num_smc_circumferentially) + grid->neq_ec * (grid->num_ec_axially * grid->num_ec_circumferentially);

	for (int i = 0; i < 4; i++) {
		grid->nbrs[local][i] = MPI_PROC_NULL;
		grid->nbrs[remote][i] = MPI_PROC_NULL;
	}
	for (int i = 0; i < 4; i++) {
		grid->flip_array[i] = 0;
	}
}

/**
 * Create a Cartesian grid for the components of a bifurcation, find the send and receive addresses for the edge tasks.
 */
// TODO: Pass a pointer, not a copy.
grid_parms make_bifucation_cart_grids(grid_parms grid)
{
	// Since there are 3 branches, there needs to be three values of a variable colour,
	// to identify the grouping of a rank to a particular sub-universe partitioned out of MPI_COMM_WORLD.
	int color = int(grid.universal_rank / (grid.m * grid.n));
	int key = 0;

	MPI_Comm split_comm;

	CHECK_MPI_ERROR(MPI_Comm_split(grid.universe, color, key, &split_comm));

	/// Parameters for cart create call.
	int ndims, nbrs[4], dims[2], periodic[2], reorder = 0, coords[2];
	ndims = 2;
	dims[0] = grid.m;
	dims[1] = grid.n;
	periodic[0] = 0;
	periodic[1] = 1;
	reorder = 0;

	CHECK_MPI_ERROR(MPI_Cart_create(split_comm, ndims, dims, periodic, reorder, &grid.cart_comm));
	CHECK_MPI_ERROR(MPI_Comm_rank(grid.cart_comm, &grid.rank_branch));
	CHECK_MPI_ERROR(MPI_Comm_size(grid.cart_comm, &grid.num_ranks_branch));

	// The inverse mapping, rank-to-coordinates translation.
	CHECK_MPI_ERROR(MPI_Cart_coords(grid.cart_comm, grid.rank_branch, ndims, grid.coords));

	CHECK_MPI_ERROR(MPI_Cart_shift(grid.cart_comm, 0, 1, &grid.nbrs[local][UP], &grid.nbrs[local][DOWN]));
	CHECK_MPI_ERROR(MPI_Cart_shift(grid.cart_comm, 1, 1, &grid.nbrs[local][LEFT], &grid.nbrs[local][RIGHT]));

	// Identifying remote neighbours.
	grid.offset_P = 0;
	grid.offset_L = (grid.m * grid.n) + ((grid.m - 1) * grid.n);
	grid.offset_R = 2 * (grid.m * grid.n) + ((grid.m - 1) * grid.n);

	if (color == 0)
	{
		grid.branch_tag = P;
	}
	else if (color == 1)
	{
		grid.branch_tag = L;
	}
	else if (color == 2)
	{
		grid.branch_tag = R;
	}

	// Label the ranks on the subdomain edges of a STRAIGHT SEGMENT as top or bottom boundary.
	for (int i = 0; i < (3 * grid.m * grid.n); i++)
	{
		if (grid.branch_tag == P)
		{
			if ((grid.rank_branch >= ((grid.m - 1) * grid.n)) && (grid.rank_branch <= (grid.m * grid.n - 1)))
			{
				grid.boundary_tag = 'B';
			}
			else
			{
				grid.boundary_tag = 'N';
			}
		}
		else if ((grid.branch_tag == L) || (grid.branch_tag == R))
		{
			if ((grid.rank_branch >= 0) && (grid.rank_branch <= (grid.n - 1)))
			{
				grid.boundary_tag = 'T';
			}
			else
			{
				grid.boundary_tag = 'N';
			}
		}
	}

	// Check whether number of processors in circumferential direction are EVEN or ODD.
	grid.scheme = grid.n % 2;

	// If number of processors in circumferential dimension is EVEN.
	if (grid.scheme == 0)
	{
		// For parent branch edge.
		if ((grid.universal_rank >= 0) && (grid.universal_rank < grid.n))
		{
			grid.boundary_tag = 'I';

			//Top edge which couples to left/right child branch.
			if ((grid.universal_rank - grid.offset_P) < (grid.n / 2))
			{
				grid.nbrs[remote][UP] = grid.offset_L + (grid.universal_rank - grid.offset_P);
			}
			else if ((grid.universal_rank - grid.offset_P) >= (grid.n / 2))
			{
				grid.nbrs[remote][UP] = grid.offset_R + (grid.universal_rank - grid.offset_P);

			}
		}
		// For left daughter branch edge.
		else if ((grid.universal_rank >= grid.offset_L) && (grid.universal_rank < (grid.offset_L + grid.n)))
		{
			grid.boundary_tag = 'I';

			if ((grid.universal_rank - grid.offset_L) < (grid.n / 2))
			{
				grid.nbrs[remote][DOWN] = grid.universal_rank - grid.offset_L;
				//grid.my_domain.internal_info.half_marker = 1;

			} else if ((grid.universal_rank - grid.offset_L) >= (grid.n / 2)) {
				grid.nbrs[remote][DOWN] = (grid.offset_R + (grid.n - 1)) - (grid.universal_rank - grid.offset_L);

				grid.flip_array[DOWN] = 1;
				//grid.my_domain.internal_info.half_marker = 2;
			}
		}
		// For Right daughter branch edge.
		else if ((grid.universal_rank >= grid.offset_R) && (grid.universal_rank < (grid.offset_R + grid.n)))
		{
			grid.boundary_tag = 'I';
			if ((grid.universal_rank - grid.offset_R) < (grid.n / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.offset_L + (grid.n - 1)) - (grid.universal_rank - grid.offset_R);
				grid.flip_array[DOWN] = 1;
				//grid.my_domain.internal_info.half_marker = 2;

			}
			else if ((grid.universal_rank - grid.offset_R) >= (grid.n / 2))
			{
				grid.nbrs[remote][DOWN] = grid.universal_rank - grid.offset_R;
				//grid.my_domain.internal_info.half_marker = 1;
			}
		}
	}

	// If number of processors in circumferential dimension are ODD.
	if (grid.scheme != 0)
	{
		// The parent edge.
		if ((grid.universal_rank >= 0) && (grid.universal_rank < grid.n))
		{
			grid.boundary_tag = 'I';
			if ((grid.universal_rank - grid.offset_P) < ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][UP] = grid.offset_L + (grid.universal_rank - grid.offset_P);
				grid.nbrs[remote][UP] += grid.offset_L + (grid.universal_rank - grid.offset_P);

			}
			else if ((grid.universal_rank - grid.offset_P) > ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][UP] = grid.offset_R + (grid.universal_rank - grid.offset_P);
				grid.nbrs[remote][UP] += grid.offset_R + (grid.universal_rank - grid.offset_P);

			}
			else if ((grid.universal_rank - grid.offset_P) == ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][UP] = grid.offset_L + (grid.universal_rank - grid.offset_P);
				grid.nbrs[remote][UP] += grid.offset_R + (grid.universal_rank - grid.offset_P);
			}
		}
		//The left daughter edge
		else if ((grid.universal_rank >= grid.offset_L) && (grid.universal_rank < grid.offset_L + grid.n))
		{
			grid.boundary_tag = 'I';
			if ((grid.universal_rank - grid.offset_L) < ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.universal_rank - grid.offset_L);
				grid.nbrs[remote][DOWN] += (grid.universal_rank - grid.offset_L);
				//grid.my_domain.internal_info.half_marker = 1;

			}
			else if ((grid.universal_rank - grid.offset_L) > ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.offset_R + (grid.n - 1)) - (grid.universal_rank - grid.offset_L);
				grid.nbrs[remote][DOWN] += (grid.offset_R + (grid.n - 1)) - (grid.universal_rank - grid.offset_L);
				grid.flip_array[DOWN] = 1;
				//grid.my_domain.internal_info.half_marker = 2;

			}
			else if ((grid.universal_rank - grid.offset_L) == ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.universal_rank - grid.offset_L);
				grid.nbrs[remote][DOWN] += (grid.offset_R + (grid.n - 1)) - (grid.universal_rank - grid.offset_L);
				grid.flip_array[DOWN] = 1;
				//grid.my_domain.internal_info.half_marker = 3;
			}
		}
		// The right daughter edge.
		else if ((grid.universal_rank >= grid.offset_R) && (grid.universal_rank < grid.offset_R + grid.n))
		{
			grid.boundary_tag = 'I';
			if ((grid.universal_rank - grid.offset_R) < ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.offset_L + (grid.n - 1)) - (grid.universal_rank - grid.offset_R);
				grid.nbrs[remote][DOWN] += (grid.offset_L + (grid.n - 1)) - (grid.universal_rank - grid.offset_R);
				grid.flip_array[DOWN] = 1;
				//grid.my_domain.internal_info.half_marker = 2;

			}
			else if ((grid.universal_rank - grid.offset_R) > ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = grid.universal_rank - grid.offset_R;
				grid.nbrs[remote][DOWN] += grid.universal_rank - grid.offset_R;
				//grid.my_domain.internal_info.half_marker = 1;

			}
			else if ((grid.universal_rank - grid.offset_R) == ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.offset_L + (grid.n - 1)) - (grid.universal_rank - grid.offset_R);
				grid.nbrs[remote][DOWN] += grid.universal_rank - grid.offset_R;
				grid.flip_array[DOWN] = 1;
				//grid.my_domain.internal_info.half_marker = 3;
			}
		}
	}

	// If I am a parent branch in my domain.
	if (grid.branch_tag == P)
	{
#if 0
		// If a parent domain exists for me.
		if (grid.my_domain.parent.domain_index >= 0)
		{
			// if I am a bottom row in my m x n cart grid.
			if ((grid.rank_branch >= ((grid.m - 1) * grid.n)) && (grid.rank_branch <= (grid.m * grid.n - 1))) {
				int stride = grid.rank_branch - ((grid.m - 1) * grid.n);
				grid.nbrs[remote][DOWN] = grid.my_domain.parent.domain_start + stride;
			}
		}
#endif
	}
	// If I am a left daughter branch in my domain.
	else if (grid.branch_tag == L)
	{
#if 0
		// If a child exists from me.
		if (grid.my_domain.left_child.domain_index >= 0)
		{
			// If I am top row in my m x n cart grid.
			if ((grid.rank_branch >= 0) && (grid.rank_branch <= (grid.n - 1)))
			{
				int stride = grid.rank_branch;
				grid.nbrs[remote][UP] = grid.my_domain.left_child.domain_start + stride;
			}
		}
#endif
	}
	// If I am a right daughter branch in my domain.
	else if (grid.branch_tag == R)
	{
#if 0
		// If a child exists from me.
		if (grid.my_domain.right_child.domain_index >= 0)
		{
			// If I am top row in my m x n cart grid.
			if ((grid.rank_branch >= 0) && (grid.rank_branch <= (grid.n - 1)))
			{
				int stride = grid.rank_branch;
				grid.nbrs[remote][UP] = grid.my_domain.right_child.domain_start + stride;
			}
		}
#endif
	}

	MPI_Comm_free(&split_comm);

	// Why do we need to return the grid, if it is passed as the argument?
	return grid;
}

/**
 * Create a Cartesian grid for a tube segment, find the send and receive addresses for the edge tasks.
 */
// TODO: Pass a pointer, not a copy.
grid_parms make_straight_segment_cart_grids(grid_parms grid)
{
	// Since there no branch, all processors have same colour.
	int color = 0;
	int key = 0;

	MPI_Comm split_comm;

	CHECK_MPI_ERROR(MPI_Comm_split(grid.universe, color, key, &split_comm));

	// Global variables that are to be read by each processor.
	int ndims, nbrs[4], dims[2], periodic[2], reorder = 0, coords[2];
	ndims = 2;
	dims[0] = grid.m;
	dims[1] = grid.n;
	periodic[0] = 0;
	periodic[1] = 1;
	reorder = 0;

	CHECK_MPI_ERROR(MPI_Cart_create(split_comm, ndims, dims, periodic, reorder, &grid.cart_comm));
	CHECK_MPI_ERROR(MPI_Comm_rank(grid.cart_comm, &grid.rank_branch));
	CHECK_MPI_ERROR(MPI_Comm_size(grid.cart_comm, &grid.num_ranks_branch));

	// The inverse mapping, rank-to-coordinates translation.
	CHECK_MPI_ERROR(MPI_Cart_coords(grid.cart_comm, grid.rank_branch, ndims, grid.coords));

	CHECK_MPI_ERROR(MPI_Cart_shift(grid.cart_comm, 0, 1, &grid.nbrs[local][UP], &grid.nbrs[local][DOWN]));
	CHECK_MPI_ERROR(MPI_Cart_shift(grid.cart_comm, 1, 1, &grid.nbrs[local][LEFT], &grid.nbrs[local][RIGHT]));

	grid.branch_tag = P;

	// Label the ranks on the subdomain edges of a STRAIGHT SEGMENT as top (T) or bottom boundary (B) or none (N).
	for (int i = 0; i < (grid.m * grid.n); i++)
	{
		if ((grid.rank_branch >= ((grid.m - 1) * grid.n)) && (grid.rank_branch <= (grid.m * grid.n - 1)))
		{
			grid.boundary_tag = 'B';
		}
		else if ((grid.rank_branch >= 0) && (grid.rank_branch <= (grid.n - 1)))
		{
			grid.boundary_tag = 'T';
		}
		else
		{
			grid.boundary_tag = 'N';
		}
	}

	// Find remote nearest neighbours on remote domains.

#if 0
	// If a parent domain exists for this subdomain.
	if (grid.my_domain.parent.domain_index >= 0)
	{
		// If we are in the bottom row in our m x n cart grid.
		if ((grid.rank_branch >= ((grid.m - 1) * grid.n)) && (grid.rank_branch <= (grid.m * grid.n - 1))) {
			int stride = grid.rank_branch - ((grid.m - 1) * grid.n);
			grid.nbrs[remote][DOWN] = grid.my_domain.parent.domain_start + stride;
		}
	}

	// If a child domain exists for this subdomain.
	if (grid.my_domain.left_child.domain_index >= 0) {
		// If we are in the top row in our m x n cart grid.
		if ((grid.rank_branch >= 0) && (grid.rank_branch <= (grid.n - 1))) {
			int stride = grid.rank_branch;
			grid.nbrs[remote][UP] = grid.my_domain.left_child.domain_start + stride;
		}
	}
#endif

	MPI_Comm_free(&split_comm);

	// Why do we need to return the grid, if it is passed as the argument?
	return grid;
}

