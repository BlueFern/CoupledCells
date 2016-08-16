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

	int scaling = grid->domain_params[0][AX_SCALE] * grid->domain_params[0][CR_SCALE];
	// Sanity check: make sure the number of requested cores equals the number of cores required by the domain configuration parameters.

	// If the subdomain type is a straight segment.
	if(grid->domain_params[0][DOMAIN_TYPE] == STRSEG)
	{
		assert(grid->num_ranks == (grid->domain_params[0][AX_QUADS] * grid->domain_params[0][CR_QUADS]) / scaling); // Total number of tasks mapping to the straight segment are m x n.
	}
	// If the subdomain is a part of a bifurcation.
	else if(grid->domain_params[0][DOMAIN_TYPE] == BIF)
	{
		assert(grid->num_ranks == 3 * (grid->domain_params[0][AX_QUADS] * grid->domain_params[0][CR_QUADS]) / scaling); // Total number of tasks mapping to the bifurcation are 3 x m x n.
	}

	grid->m = grid->domain_params[0][AX_QUADS] / grid->domain_params[0][AX_SCALE];
	grid->n = grid->domain_params[0][CR_QUADS] / grid->domain_params[0][CR_SCALE];

	grid->domain_index = grid->domain_params[0][DOMAIN_NUM];
	grid->domain_type = grid->domain_params[0][DOMAIN_TYPE];

	// Each tasks now calculates the number of ECs per node.
	// Each tasks now calculates the number of SMCs per node.
	// Topological information of a functional block of coupled cells.
	// This is the minimum required to simulate a relevant coupled topology.

	grid->num_ec_axially = grid->base_ec_axially * grid->num_ec_fundblk_axially * grid->domain_params[0][AX_SCALE];
	grid->num_smc_axially = grid->num_ec_axially * grid->num_smc_fundblk_axially;
	grid->num_smc_circumferentially = grid->base_smc_circumferentially * grid->num_smc_fundblk_circumferentially * grid->domain_params[0][CR_SCALE];
	grid->num_ec_circumferentially = grid->num_smc_circumferentially * grid->num_ec_fundblk_circumferentially;

	grid->neq_ec_axially = grid->num_ec_axially * grid->neq_ec;
	grid->neq_smc_axially = grid->num_smc_axially * grid->neq_smc;

	grid->NEQ = grid->neq_smc * (grid->num_smc_axially * grid->num_smc_circumferentially) + grid->neq_ec * (grid->num_ec_axially * grid->num_ec_circumferentially);
	printf("%d, %d, %d, %d, %d, %d, %d, %d, %d\n",grid->NEQ,grid->num_ec_circumferentially,grid->num_smc_circumferentially,grid->num_smc_axially,grid->num_ec_axially,
			grid->domain_type,grid->domain_index,grid->n,grid->m);

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
void make_bifucation_cart_grids(grid_parms *grid)
{
	// Since there are 3 branches, there needs to be three values of a variable colour,
	// to identify the grouping of a rank to a particular sub-universe partitioned out of MPI_COMM_WORLD.
	int color = int(grid->universal_rank / (grid->m * grid->n));
	int key = 0;

	MPI_Comm split_comm;

	CHECK_MPI_ERROR(MPI_Comm_split(grid->universe, color, key, &split_comm));

	/// Parameters for cart create call.
	int ndims, nbrs[4], dims[2], periodic[2], reorder = 0, coords[2];
	ndims = 2;
	dims[0] = grid->m;
	dims[1] = grid->n;
	periodic[0] = 0;
	periodic[1] = 1;
	reorder = 0;

	CHECK_MPI_ERROR(MPI_Cart_create(split_comm, ndims, dims, periodic, reorder, &grid->cart_comm));
	CHECK_MPI_ERROR(MPI_Comm_rank(grid->cart_comm, &grid->rank_branch));
	CHECK_MPI_ERROR(MPI_Comm_size(grid->cart_comm, &grid->num_ranks_branch));

	// The inverse mapping, rank-to-coordinates translation.
	CHECK_MPI_ERROR(MPI_Cart_coords(grid->cart_comm, grid->rank_branch, ndims, grid->coords));

	CHECK_MPI_ERROR(MPI_Cart_shift(grid->cart_comm, 0, 1, &grid->nbrs[local][UP], &grid->nbrs[local][DOWN]));
	CHECK_MPI_ERROR(MPI_Cart_shift(grid->cart_comm, 1, 1, &grid->nbrs[local][LEFT], &grid->nbrs[local][RIGHT]));

	if (color == 0)
	{
		grid->branch_tag = P;
	}
	else if (color == 1)
	{
		grid->branch_tag = L;
	}
	else if (color == 2)
	{
		grid->branch_tag = R;
	}

	// Label the ranks on the subdomain edges of a STRAIGHT SEGMENT as top or bottom boundary.
	if (grid->branch_tag == P)
	{
		// Last row.
		if ((grid->rank_branch >= ((grid->m - 1) * grid->n)) && (grid->rank_branch <= (grid->m * grid->n - 1)))
		{
			grid->boundary_tag = 'B';
		}
		// Other.
		else
		{
			grid->boundary_tag = 'N';
		}
	}
	else if ((grid->branch_tag == L) || (grid->branch_tag == R))
	{
		// Last row.
		if ((grid->rank_branch >= 0) && (grid->rank_branch <= (grid->n - 1)))
		{
			grid->boundary_tag = 'T';
		}
		// Other.
		else
		{
			grid->boundary_tag = 'N';
		}
	}

	// Identifying remote neighbours in universal ranks.
	grid->offset_P = 0;
	grid->offset_L = (grid->m * grid->n) + ((grid->m - 1) * grid->n);
	grid->offset_R = 2 * (grid->m * grid->n) + ((grid->m - 1) * grid->n);

	// For parent branch edge.
	if ((grid->universal_rank >= 0) && (grid->universal_rank < grid->n))
	{
		grid->boundary_tag = 'I';

		// Top edge which couples to left/right child branch.
		if ((grid->universal_rank - grid->offset_P) < (grid->n / 2))
		{
			grid->nbrs[remote][UP] = grid->offset_L + (grid->universal_rank - grid->offset_P);
		}
		else if ((grid->universal_rank - grid->offset_P) >= (grid->n / 2))
		{
			grid->nbrs[remote][UP] = grid->offset_R + (grid->universal_rank - grid->offset_P);
		}
	}
	// For left daughter branch edge.
	else if ((grid->universal_rank >= grid->offset_L) && (grid->universal_rank < (grid->offset_L + grid->n)))
	{
		grid->boundary_tag = 'I';

		if ((grid->universal_rank - grid->offset_L) < (grid->n / 2))
		{
			grid->nbrs[remote][DOWN] = grid->universal_rank - grid->offset_L;

		} else if ((grid->universal_rank - grid->offset_L) >= (grid->n / 2)) {
			grid->nbrs[remote][DOWN] = (grid->offset_R + (grid->n - 1)) - (grid->universal_rank - grid->offset_L);
			grid->flip_array[DOWN] = 1;
		}
	}
	// For Right daughter branch edge.
	else if ((grid->universal_rank >= grid->offset_R) && (grid->universal_rank < (grid->offset_R + grid->n)))
	{
		grid->boundary_tag = 'I';
		if ((grid->universal_rank - grid->offset_R) < (grid->n / 2))
		{
			grid->nbrs[remote][DOWN] = (grid->offset_L + (grid->n - 1)) - (grid->universal_rank - grid->offset_R);
			grid->flip_array[DOWN] = 1;

		}
		else if ((grid->universal_rank - grid->offset_R) >= (grid->n / 2))
		{
			grid->nbrs[remote][DOWN] = grid->universal_rank - grid->offset_R;
		}
	}

	MPI_Comm_free(&split_comm);
}

/**
 * Create a Cartesian grid for a tube segment, find the send and receive addresses for the edge tasks.
 */
void make_straight_cart_grid(grid_parms *grid)
{
	// Global variables that are to be read by each processor.
	int ndims, nbrs[4], dims[2], periodic[2], reorder = 0, coords[2];
	ndims = 2;
	dims[0] = grid->m;
	dims[1] = grid->n;
	periodic[0] = 0;
	periodic[1] = 1;
	reorder = 0;

	CHECK_MPI_ERROR(MPI_Cart_create(grid->universe, ndims, dims, periodic, reorder, &grid->cart_comm));
	CHECK_MPI_ERROR(MPI_Comm_rank(grid->cart_comm, &grid->rank_branch));
	CHECK_MPI_ERROR(MPI_Comm_size(grid->cart_comm, &grid->num_ranks_branch));

	// The inverse mapping, rank-to-coordinates translation.
	CHECK_MPI_ERROR(MPI_Cart_coords(grid->cart_comm, grid->rank_branch, ndims, grid->coords));

	CHECK_MPI_ERROR(MPI_Cart_shift(grid->cart_comm, 0, 1, &grid->nbrs[local][UP], &grid->nbrs[local][DOWN]));
	CHECK_MPI_ERROR(MPI_Cart_shift(grid->cart_comm, 1, 1, &grid->nbrs[local][LEFT], &grid->nbrs[local][RIGHT]));

	grid->branch_tag = P;

	// Label the ranks on the subdomain edges of a STRAIGHT SEGMENT as top (T) or bottom boundary (B) or none (N).
	for (int i = 0; i < (grid->m * grid->n); i++)
	{
		// Last row.
		if ((grid->rank_branch >= ((grid->m - 1) * grid->n)) && (grid->rank_branch <= (grid->m * grid->n - 1)))
		{
			grid->boundary_tag = 'B';
		}
		// First row.
		else if ((grid->rank_branch >= 0) && (grid->rank_branch <= (grid->n - 1)))
		{
			grid->boundary_tag = 'T';
		}
		// Other.
		else
		{
			grid->boundary_tag = 'N';
		}
	}
}

