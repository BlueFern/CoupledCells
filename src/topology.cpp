/* This file declears the functions which, when given a group of processors, make either
 *  a bifurcation or a straight segment, and recognize its remote and local neighbours.
 *  These function are called in main.cpp.
 */

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"

grid_parms make_subdomains(grid_parms grid, int num_subdomains, int** domains){
    int **subdomain_extents; //Element 1: offset , Element 2: Start universal_rank, Element 3: End universal_rank;
    subdomain_extents = (int**) checked_malloc(num_subdomains * sizeof(int*),
	     "Subdomain information allocation");
    for (int i = 0; i < num_subdomains; i++)
	{
	subdomain_extents[i] = (int*) checked_malloc(3 * sizeof(int),
		"Subdomains array elements allocation");
	}
    for (int i = 0; i < num_subdomains; i++)
	{
	subdomain_extents[i][0] = 0;
	subdomain_extents[i][1] = 0;
	subdomain_extents[i][2] = 0;
	}

    int a;
    for (int i = 0; i < num_subdomains; i++)
	{
	if ((domains[i][1] == 0) || (domains[i][1] == 2)
		|| (domains[i][3] == 3))
	    {
	    a = domains[i][2] * domains[i][3];
	    }
	else if (domains[i][1] == 1)
	    {
	    a = 3 * domains[i][2] * domains[i][3];
	    }
	//setting up the offset bit
	if (i == 0)
	    {
	    subdomain_extents[i][0] = 0;
	    }
	else
	    {
	    if ((domains[i - 1][1] == 0) || (domains[i - 1][1] == 2)
		    || (domains[i - 1][3] == 3))
		{
		subdomain_extents[i][0] = subdomain_extents[i - 1][0]
			+ (domains[i - 1][2] * domains[i - 1][3]);
		}
	    else if (domains[i - 1][1] == 1)
		{
		subdomain_extents[i][0] = subdomain_extents[i - 1][0]
			+ (3 * domains[i - 1][2] * domains[i - 1][3]);
		}
	    }

	subdomain_extents[i][1] = subdomain_extents[i][0]; //start universal_rank in MPI_COMM_WORLD
	subdomain_extents[i][2] = subdomain_extents[i][0] + a - 1; //end universal_rank in MPI_COMM_WORLD
	}

    /* Now all processors have the information where each domain starts and ends. Using this information, each processor can identify which domain
     * it belongs and can mark a color (0 to num_domains-1). This color can now be used to split the MPI_COMM_WORLD into sub_domains.
     * Identify the new reordered ranks in grid.sub_universe_ranks in these new communicators recorded in grid.sub_universe and update the size of this sub_domain in
     * grid.sub_universe_numtasks.
     * Since each processor has the information of its parent and child(ren) domains in domain[][] array, use this to update the my_tree structure.
     * Update remote nearest neighbour locations accordingly.
     */

    for (int i = 0; i < num_subdomains; i++)
	{
	if ((grid.universal_rank >= subdomain_extents[i][1])
		&& (grid.universal_rank <= subdomain_extents[i][2]))
	    {
	    grid.my_domain_color = i;
	    grid.my_domain_key = 0;

	    grid.m = domains[i][2];
	    grid.n = domains[i][3];

	    grid.my_domain.internal_info.domain_index = domains[i][0];
	    grid.my_domain.internal_info.domain_type = domains[i][1];
	    grid.my_domain.internal_info.domain_start = subdomain_extents[i][1];
	    grid.my_domain.internal_info.domain_end = subdomain_extents[i][2];
	    grid.my_domain.internal_info.parent_branch_case_bifurcation = -1;

	    grid.my_domain.parent.domain_index = domains[i][4];
	    //If I have a parent domain
	    if (grid.my_domain.parent.domain_index >= 0)
		{
		grid.my_domain.parent.domain_type =
			domains[grid.my_domain.parent.domain_index][1];
		//Now decide which Branch in the parent domain (in case of a bifurcation) do I belong to as a child
		if (grid.my_domain.parent.domain_type == BIF)
		    {
		    if (grid.my_domain.internal_info.domain_index
			    == domains[grid.my_domain.parent.domain_index][5])
			{
			grid.my_domain.internal_info.parent_branch_case_bifurcation =
				L;
			}
		    else if (grid.my_domain.internal_info.domain_index
			    == domains[grid.my_domain.parent.domain_index][6])
			{
			grid.my_domain.internal_info.parent_branch_case_bifurcation =
				R;
			}
		    }
		else if (grid.my_domain.parent.domain_type == STRSEG)
		    {
		    grid.my_domain.internal_info.parent_branch_case_bifurcation =
			    -1;
		    }
		//If my parent is a bifurcation decide accordingly if I am a child from left branch or right. The ranks stored for the domains will be of MPI_COMM_WORLD
		if (grid.my_domain.internal_info.parent_branch_case_bifurcation
			>= 0)
		    {
		    if (grid.my_domain.internal_info.parent_branch_case_bifurcation
			    == L)
			{
			grid.my_domain.parent.domain_start =
				subdomain_extents[grid.my_domain.parent.domain_index][1]
					+ (grid.m * grid.n);
			grid.my_domain.parent.domain_end =
				grid.my_domain.parent.domain_start
					+ (grid.n - 1);
			}
		    else if (grid.my_domain.internal_info.parent_branch_case_bifurcation
			    == R)
			{
			grid.my_domain.parent.domain_start =
				subdomain_extents[grid.my_domain.parent.domain_index][1]
					+ 2 * (grid.m * grid.n);
			grid.my_domain.parent.domain_end =
				grid.my_domain.parent.domain_start
					+ (grid.n - 1);
			}
		    }
		else if ((grid.my_domain.internal_info.parent_branch_case_bifurcation
			== -1) && (grid.my_domain.parent.domain_type == STRSEG))
		    {
		    grid.my_domain.parent.domain_start =
			    subdomain_extents[grid.my_domain.parent.domain_index][1];
		    grid.my_domain.parent.domain_end =
			    grid.my_domain.parent.domain_start + (grid.n - 1);
		    }
		}
	    else
		{
		grid.my_domain.parent.domain_type = -1;
		grid.my_domain.parent.domain_start = -1;
		grid.my_domain.parent.domain_end = -1;
		}

	    grid.my_domain.left_child.domain_index = domains[i][5];
	    //In case I have a child domain
	    if (grid.my_domain.left_child.domain_index >= 0)
		{
		grid.my_domain.left_child.domain_type =
			domains[grid.my_domain.left_child.domain_index][1];
		//Irrespective of the domain type, the last row of my child's m by n grid is of my interest
		grid.my_domain.left_child.domain_start =
			subdomain_extents[grid.my_domain.left_child.domain_index][1]
				+ ((grid.m - 1) * grid.n);
		grid.my_domain.left_child.domain_end =
			grid.my_domain.left_child.domain_start + (grid.n - 1);
		}
	    else
		{
		grid.my_domain.left_child.domain_type = -1;
		grid.my_domain.left_child.domain_start = -1;
		grid.my_domain.left_child.domain_end = -1;
		}

	    grid.my_domain.right_child.domain_index = domains[i][6];
	    if (grid.my_domain.right_child.domain_index >= 0)
		{
		grid.my_domain.right_child.domain_type =
			domains[grid.my_domain.right_child.domain_index][1];
		//Irrespective of the domain type, the last row of my child's m by n grid is of my interest
		grid.my_domain.right_child.domain_start =
			subdomain_extents[grid.my_domain.right_child.domain_index][1]
				+ ((grid.m - 1) * grid.n);
		grid.my_domain.right_child.domain_end =
			grid.my_domain.right_child.domain_start + (grid.n - 1);
		}
	    else
		{
		grid.my_domain.right_child.domain_type = -1;
		grid.my_domain.right_child.domain_start = -1;
		grid.my_domain.right_child.domain_end = -1;
		}

	    }
	}
    check_flag(
	    MPI_Comm_split(grid.universe, grid.my_domain_color,
		    grid.my_domain_key, &grid.sub_universe),
	    "Comm-split failed at subdomain level.");

    //Reveal information of myself and size of grid.sub_universe
    check_flag(MPI_Comm_rank(grid.sub_universe, &grid.sub_universe_rank),
	    "error retrieving Subdomain_rank");
    check_flag(MPI_Comm_size(grid.sub_universe, &grid.sub_universe_numtasks),
	    "error retrieving Subdomain_size");
    return grid;
    } //end of make_subdomains()

grid_parms set_geometry_parameters(grid_parms grid, int e, int s){
	///Each tasks now calculates the number of ECs per node.
	//	if (grid.m != (grid.numtasks / grid.n))
	//		e = grid.m / grid.numtasks;

	///Each tasks now calculates the number of ECs per node.
		///topological information of a functional block of coupled cells. This is the minimum required to simulate a relevant coupled topology.
		grid.num_smc_fundblk_circumferentially = 1; grid.num_ec_fundblk_circumferentially =
				5; grid.num_smc_fundblk_axially = 13; grid.num_ec_fundblk_axially =
				1;

		grid.num_ghost_cells = 2;

		grid.num_fluxes_smc = 12;///number of SMC Ioinic currents to be evaluated for eval of LHS of the d/dt terms of the ODEs.
		grid.num_fluxes_ec = 12;///number of EC Ioinic currents to be evaluated for eval of LHS of the d/dt terms of the ODEs.

		grid.num_coupling_species_smc = 3;///number of SMC coupling species homogenic /heterogenic
		grid.num_coupling_species_ec = 3;///number of SMC coupling species homogenic /heterogenic

		grid.neq_smc = 5;			/// number of SMC ODEs for a single cell
		grid.neq_ec = 4;			/// number of EC ODEs for a single cell

		grid.num_ec_axially = e * 1;
		grid.num_smc_axially = e * 13;
		grid.num_ec_circumferentially = s * 5;
		grid.num_smc_circumferentially = s * 1;

		grid.neq_ec_axially = grid.num_ec_axially * grid.neq_ec;
		grid.neq_smc_axially = grid.num_smc_axially * grid.neq_smc;

		for (int i=0; i<4; i++){
				grid.nbrs[local][i]	= MPI_PROC_NULL;
				grid.nbrs[remote][i] = MPI_PROC_NULL;
			}
		for(int i=0; i<4; i++){
				grid.flip_array[i]	=	0;
			}
		return grid;
}// end of set_geometry_parameters()

grid_parms make_bifucation(grid_parms grid)
    {
    //Since there are 3 branches, there needs to be three values of a variable color, to identify association of a rank to a particular sub-universe partitioned out of MPI_COMM_WORLD.

    grid.color = int(grid.sub_universe_rank / (grid.m * grid.n));
    grid.key = 0;			//grid.color * ((grid.m * grid.n) - 1);

    check_flag(
	    MPI_Comm_split(grid.sub_universe, grid.color, grid.key,
		    &grid.split_comm), "Comm-split failed");

    ///Global variables that are to be read by each processor
    int ndims, nbrs[4], dims[2], periodic[2], reorder = 0, coords[2];
    ndims = 2;
    dims[0] = grid.m;
    dims[1] = grid.n;
    periodic[0] = 0;
    periodic[1] = 1;
    reorder = 0;

    check_flag(
	    MPI_Cart_create(grid.split_comm, ndims, dims, periodic, reorder,
		    &grid.cart_comm),"failed at cart create");
    check_flag(MPI_Comm_rank(grid.cart_comm, &grid.rank), "failed at comm rank");
    check_flag(MPI_Cart_coords(grid.cart_comm, grid.rank, ndims, grid.coords),
	     "failed at cart coords");

    check_flag(
	    MPI_Cart_shift(grid.cart_comm, 0, 1, &grid.nbrs[local][UP],
		    &grid.nbrs[local][DOWN]), "failed at cart shift up down");
    check_flag(
	    MPI_Cart_shift(grid.cart_comm, 1, 1, &grid.nbrs[local][LEFT],
		    &grid.nbrs[local][RIGHT]),"failed at cart left right");

    //Identifying remote neighbours

    grid.offset_P = 0;
    grid.offset_L = (grid.m * grid.n) + ((grid.m - 1) * grid.n);
    grid.offset_R = 2 * (grid.m * grid.n) + ((grid.m - 1) * grid.n);
    
  

    //check whether number of processors in circumferential direction are EVEN or ODD.
    grid.scheme = grid.n % 2;
    if (grid.color == 0){
	grid.branch_tag = P;
    }
    else if (grid.color == 1){
	grid.branch_tag	= L;
    }
    else if (grid.color == 2){
	grid.branch_tag = R;
    }
  
    //Label the ranks on the subdomain edges of at STRAIGHT SEGMENT as top or bottom boundary.
	for (int i = 0; i < (3 * grid.m * grid.n); i++) {
		if (grid.branch_tag == P) {
			if ((grid.rank >= ((grid.m - 1) * grid.n))
					&& (grid.rank <= (grid.m * grid.n - 1))) {
				grid.my_domain.internal_info.boundary_tag = 'B';
			} else {
				grid.my_domain.internal_info.boundary_tag = 'N';
			}
		} else if ((grid.branch_tag == L) || (grid.branch_tag == R)) {
			if ((grid.rank >= 0) && (grid.rank <= (grid.n - 1))) {
				grid.my_domain.internal_info.boundary_tag = 'T';
			} else {
				grid.my_domain.internal_info.boundary_tag = 'N';
			}
		}
	}

    //If number of processors in circumferentail dimension are EVEN
    if (grid.scheme == 0)
	{
	//For parent branch edge
	if ((grid.sub_universe_rank >= 0) && (grid.sub_universe_rank < grid.n))
	    {
		grid.my_domain.internal_info.boundary_tag = 'I';
	    if ((grid.sub_universe_rank - grid.offset_P) < (grid.n / 2))
		{
		grid.nbrs[remote][UP1] = grid.offset_L
			+ (grid.sub_universe_rank - grid.offset_P);
		grid.nbrs[remote][UP2] = grid.offset_L
			+ (grid.sub_universe_rank - grid.offset_P);
		}
	    else if ((grid.sub_universe_rank - grid.offset_P) >= (grid.n / 2))
		{
		grid.nbrs[remote][UP1] = grid.offset_R
			+ (grid.sub_universe_rank - grid.offset_P);
		grid.nbrs[remote][UP2] = grid.offset_R
			+ (grid.sub_universe_rank - grid.offset_P);
		}
	    //For Left daughter branch edge
	    }
	else if ((grid.sub_universe_rank >= grid.offset_L)
		&& (grid.sub_universe_rank < (grid.offset_L + grid.n)))
	    {
		grid.my_domain.internal_info.boundary_tag = 'I';
	    if ((grid.sub_universe_rank - grid.offset_L) < (grid.n / 2))
		{
		grid.nbrs[remote][DOWN1] = grid.sub_universe_rank
			- grid.offset_L;
		grid.nbrs[remote][DOWN2] = grid.sub_universe_rank
			- grid.offset_L;
		}
	    else if ((grid.sub_universe_rank - grid.offset_L) >= (grid.n / 2))
		{
		grid.nbrs[remote][DOWN1] = (grid.offset_R + (grid.n - 1))
			- (grid.sub_universe_rank - grid.offset_L);
		grid.nbrs[remote][DOWN2] = (grid.offset_R + (grid.n - 1))
			- (grid.sub_universe_rank - grid.offset_L);
		grid.flip_array[DOWN1] = 1;
		grid.flip_array[DOWN2] = 1;
		}
	    }
	//For Right daughter branch edge
	else if ((grid.sub_universe_rank >= grid.offset_R)
		&& (grid.sub_universe_rank < (grid.offset_R + grid.n)))
	    {
		grid.my_domain.internal_info.boundary_tag = 'I';
	    if ((grid.sub_universe_rank - grid.offset_R) < (grid.n / 2))
		{
		grid.nbrs[remote][DOWN1] = (grid.offset_L + (grid.n - 1))
			- (grid.sub_universe_rank - grid.offset_R);
		grid.nbrs[remote][DOWN2] = (grid.offset_L + (grid.n - 1))
			- (grid.sub_universe_rank - grid.offset_R);
		grid.flip_array[DOWN1] = 1;
		grid.flip_array[DOWN2] = 1;
		}
	    else if ((grid.sub_universe_rank - grid.offset_R) >= (grid.n / 2))
		{
		grid.nbrs[remote][DOWN1] = grid.sub_universe_rank
			- grid.offset_R;
		grid.nbrs[remote][DOWN2] = grid.sub_universe_rank
			- grid.offset_R;
		}
	    }
	}

    //In the case of n being ODD

    if (grid.scheme != 0)
	{
	//The parent artery edge
	if ((grid.sub_universe_rank >= 0) && (grid.sub_universe_rank < grid.n))
	    {
		grid.my_domain.internal_info.boundary_tag = 'I';
	    if ((grid.sub_universe_rank - grid.offset_P) < ((grid.n - 1) / 2))
		{
		grid.nbrs[remote][UP1] = grid.offset_L
			+ (grid.sub_universe_rank - grid.offset_P);
		grid.nbrs[remote][UP2] = grid.offset_L
			+ (grid.sub_universe_rank - grid.offset_P);
		}
	    else if ((grid.sub_universe_rank - grid.offset_P)
		    > ((grid.n - 1) / 2))
		{
		grid.nbrs[remote][UP1] = grid.offset_R
			+ (grid.sub_universe_rank - grid.offset_P);
		grid.nbrs[remote][UP2] = grid.offset_R
			+ (grid.sub_universe_rank - grid.offset_P);
		}
	    else if ((grid.sub_universe_rank - grid.offset_P)
		    == ((grid.n - 1) / 2))
		{
		grid.nbrs[remote][UP1] = grid.offset_L
			+ (grid.sub_universe_rank - grid.offset_P);
		grid.nbrs[remote][UP2] = grid.offset_R
			+ (grid.sub_universe_rank - grid.offset_P);
		}
	    }
	//The left daughter artery edge
	else if ((grid.sub_universe_rank >= grid.offset_L)
		&& (grid.sub_universe_rank < grid.offset_L + grid.n))
	    {
		grid.my_domain.internal_info.boundary_tag = 'I';
	    if ((grid.sub_universe_rank - grid.offset_L) < ((grid.n - 1) / 2))
		{
		grid.nbrs[remote][DOWN1] = (grid.sub_universe_rank
			- grid.offset_L);
		grid.nbrs[remote][DOWN2] = (grid.sub_universe_rank
			- grid.offset_L);
		}
	    else if ((grid.sub_universe_rank - grid.offset_L)
		    > ((grid.n - 1) / 2))
		{
		grid.nbrs[remote][DOWN1] = (grid.offset_R + (grid.n - 1))
			- (grid.sub_universe_rank - grid.offset_L);
		grid.nbrs[remote][DOWN2] = (grid.offset_R + (grid.n - 1))
			- (grid.sub_universe_rank - grid.offset_L);
		grid.flip_array[DOWN1] = 1;
		grid.flip_array[DOWN2] = 1;
		}
	    else if ((grid.sub_universe_rank - grid.offset_L)
		    == ((grid.n - 1) / 2))
		{
		grid.nbrs[remote][DOWN1] = (grid.sub_universe_rank
			- grid.offset_L);
		grid.nbrs[remote][DOWN2] = (grid.offset_R + (grid.n - 1))
			- (grid.sub_universe_rank - grid.offset_L);
		grid.flip_array[DOWN1] = 0;
		grid.flip_array[DOWN2] = 1;
		}
	    }
	//The right daughter artery edge
	else if ((grid.sub_universe_rank >= grid.offset_R)
		&& (grid.sub_universe_rank < grid.offset_R + grid.n))
	    {
		grid.my_domain.internal_info.boundary_tag = 'I';
	    if ((grid.sub_universe_rank - grid.offset_R) < ((grid.n - 1) / 2))
		{
		grid.nbrs[remote][DOWN1] = (grid.offset_L + (grid.n - 1))
			- (grid.sub_universe_rank - grid.offset_R);
		grid.nbrs[remote][DOWN2] = (grid.offset_L + (grid.n - 1))
			- (grid.sub_universe_rank - grid.offset_R);
		grid.flip_array[DOWN1] = 1;
		grid.flip_array[DOWN2] = 1;
		}
	    else if ((grid.sub_universe_rank - grid.offset_R)
		    > ((grid.n - 1) / 2))
		{
		grid.nbrs[remote][DOWN1] = grid.sub_universe_rank
			- grid.offset_R;
		grid.nbrs[remote][DOWN2] = grid.sub_universe_rank
			- grid.offset_R;
		}
	    else if ((grid.sub_universe_rank - grid.offset_R)
		    == ((grid.n - 1) / 2))
		{
		grid.nbrs[remote][DOWN1] = (grid.offset_L + (grid.n - 1))
			- (grid.sub_universe_rank - grid.offset_R);
		grid.nbrs[remote][DOWN2] = grid.sub_universe_rank
			- grid.offset_R;
		grid.flip_array[DOWN1] = 1;
		grid.flip_array[DOWN2] = 0;
		}
	    }
	}

    //If I am a parent branch in my domain
    if (grid.branch_tag == P)
	{
	//if a parent domain exists for me
	if (grid.my_domain.parent.domain_index >= 0)
	    {
	    //if I am a bottom row in my m x n cart grid
	    if ((grid.rank >= ((grid.m - 1) * grid.n))
		    && (grid.rank <= (grid.m * grid.n - 1)))
		{
		int stride = grid.rank - ((grid.m - 1) * grid.n);
		grid.nbrs[remote][DOWN1] = grid.my_domain.parent.domain_start
			+ stride;
		grid.nbrs[remote][DOWN2] = grid.my_domain.parent.domain_start
			+ stride;
		}
	    }
	}
    //If I am a left daughter branch in my domain
    else if (grid.branch_tag == L)
	{
	//if a child exists from me
	if (grid.my_domain.left_child.domain_index >= 0)
	    {
	    //if I am top row in my m x n cart grid
	    if ((grid.rank >= 0) && (grid.rank <= (grid.n - 1)))
		{
		int stride = grid.rank;
		grid.nbrs[remote][UP1] = grid.my_domain.left_child.domain_start
			+ stride;
		grid.nbrs[remote][UP2] = grid.my_domain.left_child.domain_start
			+ stride;
		}
	    }
	}
    //If I am a right daughter branch in my domain
    else if (grid.branch_tag == R)
	{
	//if a child exists from me
	if (grid.my_domain.right_child.domain_index >= 0)
	    {
	    //if I am top row in my m x n cart grid
	    if ((grid.rank >= 0) && (grid.rank <= (grid.n - 1)))
		{
		int stride = grid.rank;
		grid.nbrs[remote][UP1] = grid.my_domain.right_child.domain_start
			+ stride;
		grid.nbrs[remote][UP2] = grid.my_domain.right_child.domain_start
			+ stride;
		}
	    }
	}
    return grid;
    }			// end of make_bifurcation


grid_parms make_straight_segment(grid_parms grid)
    {
    //Since there no branch, all processors have same color.

    grid.color = 0;
    grid.key = 0;

    check_flag(
	    MPI_Comm_split(grid.sub_universe, grid.color, grid.key,
		    &grid.split_comm), "Comm-split failed");

    ///Global variables that are to be read by each processor
    int ndims, nbrs[4], dims[2], periodic[2], reorder = 0, coords[2];
    ndims = 2;
    dims[0] = grid.m;
    dims[1] = grid.n;
    periodic[0] = 0;
    periodic[1] = 1;
    reorder = 0;

    check_flag(
	    MPI_Cart_create(grid.split_comm, ndims, dims, periodic, reorder,
		    &grid.cart_comm), "failed at cart create");
    check_flag(MPI_Comm_rank(grid.cart_comm, &grid.rank),"failed at comm rank");
    check_flag(MPI_Cart_coords(grid.cart_comm, grid.rank, ndims, grid.coords),"failed at cart coords");

    check_flag(
	    MPI_Cart_shift(grid.cart_comm, 0, 1, &grid.nbrs[local][UP],
		    &grid.nbrs[local][DOWN]),"failed at cart shift up down");
    check_flag(
	    MPI_Cart_shift(grid.cart_comm, 1, 1, &grid.nbrs[local][LEFT],
		    &grid.nbrs[local][RIGHT]),"failed at cart left right");

    //Label the ranks on the subdomain edges of at STRAIGHT SEGMENT as top or bottom boundary.
	for (int i = 0; i < (grid.m * grid.n); i++) {
		if ((grid.rank >= ((grid.m - 1) * grid.n))
				&& (grid.rank <= (grid.m * grid.n - 1))) {
			grid.my_domain.internal_info.boundary_tag = 'B';
		} else if ((grid.rank >= 0) && (grid.rank <= (grid.n - 1))) {
			grid.my_domain.internal_info.boundary_tag = 'T';
		} else {
			grid.my_domain.internal_info.boundary_tag = 'N';
		}
	}
    
    
    //Find remote nearest neighbours on remote domains
	
	//if a parent domain exists for me
	if (grid.my_domain.parent.domain_index >= 0) {
		//if I am a bottom row in my m x n cart grid
		if ((grid.rank >= ((grid.m - 1) * grid.n))
				&& (grid.rank <= (grid.m * grid.n - 1))) {
			int stride = grid.rank - ((grid.m - 1) * grid.n);
			grid.nbrs[remote][DOWN1] = grid.my_domain.parent.domain_start
					+ stride;
			grid.nbrs[remote][DOWN2] = grid.my_domain.parent.domain_start
					+ stride;
		}
	}
	//if a child exists from me
	else if (grid.my_domain.left_child.domain_index >= 0) {
		//if I am top row in my m x n cart grid
		if ((grid.rank >= 0) && (grid.rank <= (grid.n - 1))) {
			int stride = grid.rank;
			grid.nbrs[remote][UP1] = grid.my_domain.left_child.domain_start
					+ stride;
			grid.nbrs[remote][UP2] = grid.my_domain.left_child.domain_start
					+ stride;
		}
	}
    return grid;
    }			//end of make_straight_segment()


