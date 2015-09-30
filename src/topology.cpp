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
 * Calculate parameters of subdomains, including the parameters of parent and child subdomains.
 */
void configure_subdomains_topology(grid_parms *grid)
{

	int **subdomain_extents;

	// Just one subdomain in the case of a bifurcation or a straight segment.
	subdomain_extents = (int**) checked_malloc(grid->num_domains * sizeof(int*), SRC_LOC);

	for (int i = 0; i < grid->num_domains; i++) {
		subdomain_extents[i] = (int*) checked_malloc(3 * sizeof(int), SRC_LOC);
	}

	// Element 1: offset, Element 2: Start universal_rank, Element 3: End universal_rank.
	for (int i = 0; i < grid->num_domains; i++) {
		subdomain_extents[i][0] = 0;
		subdomain_extents[i][1] = 0;
		subdomain_extents[i][2] = 0;
	}

	for(int i = 0; i < grid->num_domains; i++)
	{
		// Setting up the offset bit.
		if(i == 0)
		{
			// Always zero offset on the first domain.
			subdomain_extents[i][0] = 0;
		}
		else
		{
			// WARNING: This code has never been verified to work because we always have only one domain.
			if((grid->domain_params[i - 1][DOMAIN_TYPE] == STRSEG))
			{
				// The offset for the current domain is the offset of the previous domain plus the size of the previous domain.
				subdomain_extents[i][0] = subdomain_extents[i - 1][0] + (grid->domain_params[i - 1][AX_QUADS] * grid->domain_params[i - 1][CR_QUADS]);
			}
			// If this is a bifurcation.
			else if(grid->domain_params[i - 1][DOMAIN_TYPE] == BIF)
			{
				// The offset for the current domain is the offset of the previous domain plus the size of the previous domain.
				subdomain_extents[i][0] = subdomain_extents[i - 1][0] + (3 * grid->domain_params[i - 1][AX_QUADS] * grid->domain_params[i - 1][CR_QUADS]);
			}
		}

		subdomain_extents[i][1] = subdomain_extents[i][0]; // Start universal_rank in MPI_COMM_WORLD.
		subdomain_extents[i][2] = subdomain_extents[i][0] + grid->num_ranks - 1; // End universal_rank in MPI_COMM_WORLD.
	}

	/* Now all processors have the information where each domain starts and ends. Using this information, each processor can identify which domain
	 * it belongs and can mark a colour (0 to num_domains - 1). This colour can now be used to split the MPI_COMM_WORLD into sub_domains.
	 *
	 * Identify the new reordered ranks in grid.sub_universe_ranks in these new communicators recorded in grid.sub_universe and update the size of this sub_domain in
	 * grid.sub_universe_numtasks.
	 *
	 * Since each processor has the information of its parent and child domains in domains[][] array, use this to update the my_tree structure.
	 * Update remote nearest neighbour locations accordingly.
	 */

	for(int i = 0; i < grid->num_domains; i++)
	{
		if((grid->universal_rank >= subdomain_extents[i][1]) && (grid->universal_rank <= subdomain_extents[i][2]))
		{
			//grid->my_domain_color = i;
			//grid->my_domain_key = 0;

			grid->m = grid->domain_params[i][AX_QUADS];
			grid->n = grid->domain_params[i][CR_QUADS];

			grid->my_domain.internal_info.domain_index = grid->domain_params[i][DOMAIN_NUM];
			grid->my_domain.internal_info.domain_type = grid->domain_params[i][DOMAIN_TYPE];
			grid->my_domain.internal_info.domain_start = subdomain_extents[i][1];
			grid->my_domain.internal_info.domain_end = subdomain_extents[i][2];
			grid->my_domain.internal_info.parent_branch_case_bifurcation = -1;

			grid->my_domain.parent.domain_index = grid->domain_params[i][PARENT_DOMAIN_NUM];
			grid->my_domain.parent.domain_type = -1;
			grid->my_domain.parent.domain_start = -1;
			grid->my_domain.parent.domain_end = -1;

			grid->my_domain.left_child.domain_index = grid->domain_params[i][LEFT_DOMAIN_NUM];
			grid->my_domain.left_child.domain_type = -1;
			grid->my_domain.left_child.domain_start = -1;
			grid->my_domain.left_child.domain_end = -1;

			grid->my_domain.right_child.domain_index = grid->domain_params[i][RIGHT_DOMAIN_NUM];
			grid->my_domain.right_child.domain_type = -1;
			grid->my_domain.right_child.domain_start = -1;
			grid->my_domain.right_child.domain_end = -1;
		}
	}

	for (int i = 0; i < grid->num_domains; i++) {
		free(subdomain_extents[i]);
	}
	free(subdomain_extents);

	// Since we always have only one domain, the sub-universe communicator will always be identical to the universe communicator.

	// Do the domain splitting to make subdomains.
	//CHECK_MPI_ERROR(MPI_Comm_split(grid->universe, grid->my_domain_color, grid->my_domain_key, &grid->sub_universe));

	// Reveal information of myself and size of grid.sub_universe.
	//CHECK_MPI_ERROR(MPI_Comm_rank(grid->sub_universe, &grid->sub_universe_rank));
	//CHECK_MPI_ERROR(MPI_Comm_size(grid->sub_universe, &grid->sub_universe_numtasks));
}

/**
 * Create a Cartesian grid for the components of a bifurcation, find the send and receive addresses for the edge tasks.
 */
grid_parms make_bifucation_cart_grids(grid_parms grid)
{
	// Since there are 3 branches, there needs to be three values of a variable colour,
	// to identify the grouping of a rank to a particular sub-universe partitioned out of MPI_COMM_WORLD.
	grid.color = int(grid.universal_rank / (grid.m * grid.n));
	grid.key = 0;

	MPI_Comm split_comm;

	CHECK_MPI_ERROR(MPI_Comm_split(grid.universe, grid.color, grid.key, &split_comm));

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

	if (grid.color == 0)
	{
		grid.branch_tag = P;
	}
	else if (grid.color == 1)
	{
		grid.branch_tag = L;
	}
	else if (grid.color == 2)
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
				grid.my_domain.internal_info.boundary_tag = 'B';
			}
			else
			{
				grid.my_domain.internal_info.boundary_tag = 'N';
			}
		}
		else if ((grid.branch_tag == L) || (grid.branch_tag == R))
		{
			if ((grid.rank_branch >= 0) && (grid.rank_branch <= (grid.n - 1)))
			{
				grid.my_domain.internal_info.boundary_tag = 'T';
			}
			else
			{
				grid.my_domain.internal_info.boundary_tag = 'N';
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
			grid.my_domain.internal_info.boundary_tag = 'I';

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
			grid.my_domain.internal_info.boundary_tag = 'I';

			if ((grid.universal_rank - grid.offset_L) < (grid.n / 2))
			{
				grid.nbrs[remote][DOWN] = grid.universal_rank - grid.offset_L;
				grid.my_domain.internal_info.half_marker = 1;

			} else if ((grid.universal_rank - grid.offset_L) >= (grid.n / 2)) {
				grid.nbrs[remote][DOWN] = (grid.offset_R + (grid.n - 1)) - (grid.universal_rank - grid.offset_L);

				grid.flip_array[DOWN] = 1;
				grid.my_domain.internal_info.half_marker = 2;
			}
		}
		// For Right daughter branch edge.
		else if ((grid.universal_rank >= grid.offset_R) && (grid.universal_rank < (grid.offset_R + grid.n)))
		{
			grid.my_domain.internal_info.boundary_tag = 'I';
			if ((grid.universal_rank - grid.offset_R) < (grid.n / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.offset_L + (grid.n - 1)) - (grid.universal_rank - grid.offset_R);
				grid.flip_array[DOWN] = 1;
				grid.my_domain.internal_info.half_marker = 2;

			}
			else if ((grid.universal_rank - grid.offset_R) >= (grid.n / 2))
			{
				grid.nbrs[remote][DOWN] = grid.universal_rank - grid.offset_R;
				grid.my_domain.internal_info.half_marker = 1;
			}
		}
	}

	// If number of processors in circumferential dimension are ODD.
	if (grid.scheme != 0)
	{
		// The parent edge.
		if ((grid.universal_rank >= 0) && (grid.universal_rank < grid.n))
		{
			grid.my_domain.internal_info.boundary_tag = 'I';
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
			grid.my_domain.internal_info.boundary_tag = 'I';
			if ((grid.universal_rank - grid.offset_L) < ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.universal_rank - grid.offset_L);
				grid.nbrs[remote][DOWN] += (grid.universal_rank - grid.offset_L);
				grid.my_domain.internal_info.half_marker = 1;

			}
			else if ((grid.universal_rank - grid.offset_L) > ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.offset_R + (grid.n - 1)) - (grid.universal_rank - grid.offset_L);
				grid.nbrs[remote][DOWN] += (grid.offset_R + (grid.n - 1)) - (grid.universal_rank - grid.offset_L);
				grid.flip_array[DOWN] = 1;
				grid.my_domain.internal_info.half_marker = 2;

			}
			else if ((grid.universal_rank - grid.offset_L) == ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.universal_rank - grid.offset_L);
				grid.nbrs[remote][DOWN] += (grid.offset_R + (grid.n - 1)) - (grid.universal_rank - grid.offset_L);
				grid.flip_array[DOWN] = 1;
				grid.my_domain.internal_info.half_marker = 3;
			}
		}
		// The right daughter edge.
		else if ((grid.universal_rank >= grid.offset_R) && (grid.universal_rank < grid.offset_R + grid.n))
		{
			grid.my_domain.internal_info.boundary_tag = 'I';
			if ((grid.universal_rank - grid.offset_R) < ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.offset_L + (grid.n - 1)) - (grid.universal_rank - grid.offset_R);
				grid.nbrs[remote][DOWN] += (grid.offset_L + (grid.n - 1)) - (grid.universal_rank - grid.offset_R);
				grid.flip_array[DOWN] = 1;
				grid.my_domain.internal_info.half_marker = 2;

			}
			else if ((grid.universal_rank - grid.offset_R) > ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = grid.universal_rank - grid.offset_R;
				grid.nbrs[remote][DOWN] += grid.universal_rank - grid.offset_R;
				grid.my_domain.internal_info.half_marker = 1;

			}
			else if ((grid.universal_rank - grid.offset_R) == ((grid.n - 1) / 2))
			{
				grid.nbrs[remote][DOWN] = (grid.offset_L + (grid.n - 1)) - (grid.universal_rank - grid.offset_R);
				grid.nbrs[remote][DOWN] += grid.universal_rank - grid.offset_R;
				grid.flip_array[DOWN] = 1;
				grid.my_domain.internal_info.half_marker = 3;
			}
		}
	}

	// If I am a parent branch in my domain.
	if (grid.branch_tag == P)
	{
		// If a parent domain exists for me.
		if (grid.my_domain.parent.domain_index >= 0)
		{
			// if I am a bottom row in my m x n cart grid.
			if ((grid.rank_branch >= ((grid.m - 1) * grid.n)) && (grid.rank_branch <= (grid.m * grid.n - 1))) {
				int stride = grid.rank_branch - ((grid.m - 1) * grid.n);
				grid.nbrs[remote][DOWN] = grid.my_domain.parent.domain_start + stride;
			}
		}
	}
	// If I am a left daughter branch in my domain.
	else if (grid.branch_tag == L)
	{
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
	}
	// If I am a right daughter branch in my domain.
	else if (grid.branch_tag == R)
	{
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
	}

	MPI_Comm_free(&split_comm);

	// Why do we need to return the grid, if it is passed as the argument?
	return grid;
}

/**
 * Create a Cartesian grid for a tube segment, find the send and receive addresses for the edge tasks.
 */
grid_parms make_straight_segment_cart_grids(grid_parms grid)
{
	// Since there no branch, all processors have same colour.
	grid.color = 0;
	grid.key = 0;

	MPI_Comm split_comm;

	CHECK_MPI_ERROR(MPI_Comm_split(grid.universe, grid.color, grid.key, &split_comm));

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

	// Label the ranks on the subdomain edges of a STRAIGHT SEGMENT as top (T) or bottom boundary (B) or none (N).
	for (int i = 0; i < (grid.m * grid.n); i++)
	{
		if ((grid.rank_branch >= ((grid.m - 1) * grid.n)) && (grid.rank_branch <= (grid.m * grid.n - 1)))
		{
			grid.my_domain.internal_info.boundary_tag = 'B';
		}
		else if ((grid.rank_branch >= 0) && (grid.rank_branch <= (grid.n - 1)))
		{
			grid.my_domain.internal_info.boundary_tag = 'T';
		}
		else
		{
			grid.my_domain.internal_info.boundary_tag = 'N';
		}
	}

	// Find remote nearest neighbours on remote domains.

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

	MPI_Comm_free(&split_comm);

	// Why do we need to return the grid, if it is passed as the argument?
	return grid;
}

