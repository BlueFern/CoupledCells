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

grid_parms make_subdomains(grid_parms grid, int num_subdomains, int** domains) {
	int **subdomain_extents; //Element 1: offset , Element 2: Start universal_rank, Element 3: End universal_rank;
	subdomain_extents = (int**) checked_malloc(num_subdomains * sizeof(int*), "Subdomain information allocation");
	for (int i = 0; i < num_subdomains; i++) {
		subdomain_extents[i] = (int*) checked_malloc(3 * sizeof(int), "Subdomains array elements allocation");
	}
	for (int i = 0; i < num_subdomains; i++) {
		subdomain_extents[i][0] = 0;
		subdomain_extents[i][1] = 0;
		subdomain_extents[i][2] = 0;
	}

	int a;
	for (int i = 0; i < num_subdomains; i++)	/// if the subdomain type is Straight segment
			{
		if ((domains[i][1] == 0) || (domains[i][1] == 2) || (domains[i][3] == 3)) {
			a = domains[i][2] * domains[i][3];		///total number of tasks mapping to the straight segment are m x n
		} else if (domains[i][1] == 1) {
			a = 3 * domains[i][2] * domains[i][3]; ///total number of tasks mapping to the bifurcation are 3 x m x n
		}
		//setting up the offset bit
		if (i == 0) {
			subdomain_extents[i][0] = 0;
		} else {
			if ((domains[i - 1][1] == 0) || (domains[i - 1][1] == 2) || (domains[i - 1][3] == 3)) {
				subdomain_extents[i][0] = subdomain_extents[i - 1][0] + (domains[i - 1][2] * domains[i - 1][3]);
			} else if (domains[i - 1][1] == 1) {
				subdomain_extents[i][0] = subdomain_extents[i - 1][0] + (3 * domains[i - 1][2] * domains[i - 1][3]);
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

	for (int i = 0; i < num_subdomains; i++) {
		if ((grid.universal_rank >= subdomain_extents[i][1]) && (grid.universal_rank <= subdomain_extents[i][2])) {
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
			if (grid.my_domain.parent.domain_index >= 0) {
				grid.my_domain.parent.domain_type = domains[grid.my_domain.parent.domain_index][1];
				grid.my_domain.parent.m = domains[grid.my_domain.parent.domain_index][2];
				grid.my_domain.parent.n = domains[grid.my_domain.parent.domain_index][3];
				//Now decide which Branch in the parent domain (in case of a bifurcation) do I belong to as a child
				if (grid.my_domain.parent.domain_type == BIF) {
					if (grid.my_domain.internal_info.domain_index == domains[grid.my_domain.parent.domain_index][5]) {
						grid.my_domain.internal_info.parent_branch_case_bifurcation = L;
					} else if (grid.my_domain.internal_info.domain_index == domains[grid.my_domain.parent.domain_index][6]) {
						grid.my_domain.internal_info.parent_branch_case_bifurcation = R;
					}
				} else if (grid.my_domain.parent.domain_type == STRSEG) {
					grid.my_domain.internal_info.parent_branch_case_bifurcation = -1;
				}
				//If my parent is a bifurcation decide accordingly if I am a child from left branch or right.
				//The ranks stored for the domains will be of MPI_COMM_WORLD.
				if (grid.my_domain.internal_info.parent_branch_case_bifurcation >= 0) {
					if (grid.my_domain.internal_info.parent_branch_case_bifurcation == L) {
						grid.my_domain.parent.domain_start = subdomain_extents[grid.my_domain.parent.domain_index][1]
								+ (grid.my_domain.parent.m * grid.my_domain.parent.n);
						grid.my_domain.parent.domain_end = grid.my_domain.parent.domain_start + (grid.my_domain.parent.n - 1);
					} else if (grid.my_domain.internal_info.parent_branch_case_bifurcation == R) {
						grid.my_domain.parent.domain_start = subdomain_extents[grid.my_domain.parent.domain_index][1]
								+ 2 * (grid.my_domain.parent.m * grid.my_domain.parent.n);
						grid.my_domain.parent.domain_end = grid.my_domain.parent.domain_start + (grid.my_domain.parent.n - 1);
					}
				} else if ((grid.my_domain.internal_info.parent_branch_case_bifurcation == -1) && (grid.my_domain.parent.domain_type == STRSEG)) {
					grid.my_domain.parent.domain_start = subdomain_extents[grid.my_domain.parent.domain_index][1];
					grid.my_domain.parent.domain_end = grid.my_domain.parent.domain_start + (grid.my_domain.parent.n - 1);
				}
			} else {
				grid.my_domain.parent.domain_type = -1;
				grid.my_domain.parent.domain_start = -1;
				grid.my_domain.parent.domain_end = -1;
			}

			grid.my_domain.left_child.domain_index = domains[i][5];
			//In case I have a child domain
			if (grid.my_domain.left_child.domain_index >= 0) {
				grid.my_domain.left_child.domain_type = domains[grid.my_domain.left_child.domain_index][1];
				grid.my_domain.left_child.m = domains[grid.my_domain.left_child.domain_index][2];
				grid.my_domain.left_child.n = domains[grid.my_domain.left_child.domain_index][3];
				//Irrespective of the domain type, the last row of my child's m by n grid is of my interest
				grid.my_domain.left_child.domain_start = subdomain_extents[grid.my_domain.left_child.domain_index][1]
						+ ((grid.my_domain.left_child.m - 1) * grid.my_domain.left_child.n);
				grid.my_domain.left_child.domain_end = grid.my_domain.left_child.domain_start + (grid.my_domain.left_child.n - 1);
			} else {
				grid.my_domain.left_child.domain_type = -1;
				grid.my_domain.left_child.domain_start = -1;
				grid.my_domain.left_child.domain_end = -1;
			}

			grid.my_domain.right_child.domain_index = domains[i][6];
			if (grid.my_domain.right_child.domain_index >= 0) {
				grid.my_domain.right_child.domain_type = domains[grid.my_domain.right_child.domain_index][1];
				grid.my_domain.right_child.m = domains[grid.my_domain.right_child.domain_index][2];
				grid.my_domain.right_child.n = domains[grid.my_domain.right_child.domain_index][3];

				//Irrespective of the domain type, the last row of my child's m by n grid is of my interest
				grid.my_domain.right_child.domain_start = subdomain_extents[grid.my_domain.right_child.domain_index][1]
						+ ((grid.my_domain.right_child.m - 1) * grid.my_domain.right_child.n);
				grid.my_domain.right_child.domain_end = grid.my_domain.right_child.domain_start + (grid.my_domain.right_child.n - 1);
			} else {
				grid.my_domain.right_child.domain_type = -1;
				grid.my_domain.right_child.domain_start = -1;
				grid.my_domain.right_child.domain_end = -1;
			}

		}
	}

	///Do the domain spliting to make subdomains.
	check_flag(MPI_Comm_split(grid.universe, grid.my_domain_color, grid.my_domain_key, &grid.sub_universe), "Comm-split failed at subdomain level.");

	//Reveal information of myself and size of grid.sub_universe
	check_flag(MPI_Comm_rank(grid.sub_universe, &grid.sub_universe_rank), "error retrieving Subdomain_rank");
	check_flag(MPI_Comm_size(grid.sub_universe, &grid.sub_universe_numtasks), "error retrieving Subdomain_size");

	return grid;
} //end of make_subdomains()

grid_parms set_geometry_parameters(grid_parms grid) {
	///Each tasks now calculates the number of ECs per node.

	///Each tasks now calculates the number of ECs per node.
	///topological information of a functional block of coupled cells. This is the minimum required to simulate a relevant coupled topology.
	grid.num_smc_fundblk_circumferentially = 1;
	grid.num_ec_fundblk_circumferentially = 5;
	grid.num_smc_fundblk_axially = 13;
	grid.num_ec_fundblk_axially = 1;

	grid.num_ghost_cells = 2;

	grid.num_fluxes_smc = 12;		///number of SMC Ioinic currents to be evaluated for eval of LHS of the d/dt terms of the ODEs.
	grid.num_fluxes_ec = 12;		///number of EC Ioinic currents to be evaluated for eval of LHS of the d/dt terms of the ODEs.

	grid.num_coupling_species_smc = 3;		///number of SMC coupling species homogenic /heterogenic
	grid.num_coupling_species_ec = 3;		///number of SMC coupling species homogenic /heterogenic

	grid.neq_smc = 5;			/// number of SMC ODEs for a single cell
	grid.neq_ec = 4;			/// number of EC ODEs for a single cell

	for (int i = 0; i < grid.num_domains; i++) {
		for (int j = 0; j < 9; j++) {
			grid.num_ec_axially = grid.domains[i][7] * 1;
			grid.num_smc_axially = grid.num_ec_axially * 13;
			grid.num_smc_circumferentially = grid.domains[i][8] * 1;
			grid.num_ec_circumferentially = grid.num_smc_circumferentially * 5;
		}
	}

	grid.neq_ec_axially = grid.num_ec_axially * grid.neq_ec;
	grid.neq_smc_axially = grid.num_smc_axially * grid.neq_smc;

	for (int i = 0; i < 4; i++) {
		grid.nbrs[local][i] = MPI_PROC_NULL;
		grid.nbrs[remote][i] = MPI_PROC_NULL;
	}
	for (int i = 0; i < 4; i++) {
		grid.flip_array[i] = 0;
	}

	grid.num_parameters = 2;
	return grid;
}			// end of set_geometry_parameters()

grid_parms make_bifucation(grid_parms grid) {
	//Since there are 3 branches, there needs to be three values of a variable color, to identify association of a rank to a particular sub-universe partitioned out of MPI_COMM_WORLD.

	grid.color = int(grid.sub_universe_rank / (grid.m * grid.n));
	grid.key = 0;			//grid.color * ((grid.m * grid.n) - 1);

	check_flag(MPI_Comm_split(grid.sub_universe, grid.color, grid.key, &grid.split_comm), "Comm-split failed");

	///Global variables that are to be read by each processor
	int ndims, nbrs[4], dims[2], periodic[2], reorder = 0, coords[2];
	ndims = 2;
	dims[0] = grid.m;
	dims[1] = grid.n;
	periodic[0] = 0;
	periodic[1] = 1;
	reorder = 0;

	check_flag(MPI_Cart_create(grid.split_comm, ndims, dims, periodic, reorder, &grid.cart_comm), "failed at cart create");
	check_flag(MPI_Comm_rank(grid.cart_comm, &grid.rank), "failed at cart comm rank");
	check_flag(MPI_Comm_size(grid.cart_comm, &grid.tasks), "failed at cart comm tasks");
	check_flag(MPI_Cart_coords(grid.cart_comm, grid.rank, ndims, grid.coords), "failed at cart coords");

	check_flag(MPI_Cart_shift(grid.cart_comm, 0, 1, &grid.nbrs[local][UP], &grid.nbrs[local][DOWN]), "failed at cart shift up down");
	check_flag(MPI_Cart_shift(grid.cart_comm, 1, 1, &grid.nbrs[local][LEFT], &grid.nbrs[local][RIGHT]), "failed at cart left right");

	//Identifying remote neighbours

	grid.offset_P = 0;
	grid.offset_L = (grid.m * grid.n) + ((grid.m - 1) * grid.n);
	grid.offset_R = 2 * (grid.m * grid.n) + ((grid.m - 1) * grid.n);

	//check whether number of processors in circumferential direction are EVEN or ODD.
	grid.scheme = grid.n % 2;
	if (grid.color == 0) {
		grid.branch_tag = P;
	} else if (grid.color == 1) {
		grid.branch_tag = L;
	} else if (grid.color == 2) {
		grid.branch_tag = R;
	}

	//Label the ranks on the subdomain edges of at STRAIGHT SEGMENT as top or bottom boundary.
	for (int i = 0; i < (3 * grid.m * grid.n); i++) {
		if (grid.branch_tag == P) {
			if ((grid.rank >= ((grid.m - 1) * grid.n)) && (grid.rank <= (grid.m * grid.n - 1))) {
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

	//If number of processors in circumferential dimension are EVEN
	if (grid.scheme == 0) {
		//For parent branch edge
		if ((grid.sub_universe_rank >= 0) && (grid.sub_universe_rank < grid.n)) {
			grid.my_domain.internal_info.boundary_tag = 'I';			//Top edge which couples to Left/Right child branch
			if ((grid.sub_universe_rank - grid.offset_P) < (grid.n / 2)) {
				grid.nbrs[remote][UP1] = grid.offset_L + (grid.sub_universe_rank - grid.offset_P);
				grid.nbrs[remote][UP2] = grid.offset_L + (grid.sub_universe_rank - grid.offset_P);
			} else if ((grid.sub_universe_rank - grid.offset_P) >= (grid.n / 2)) {
				grid.nbrs[remote][UP1] = grid.offset_R + (grid.sub_universe_rank - grid.offset_P);
				grid.nbrs[remote][UP2] = grid.offset_R + (grid.sub_universe_rank - grid.offset_P);
			}
			//For Left daughter branch edge
		} else if ((grid.sub_universe_rank >= grid.offset_L) && (grid.sub_universe_rank < (grid.offset_L + grid.n))) {
			grid.my_domain.internal_info.boundary_tag = 'I';
			if ((grid.sub_universe_rank - grid.offset_L) < (grid.n / 2)) {
				grid.nbrs[remote][DOWN1] = grid.sub_universe_rank - grid.offset_L;
				grid.nbrs[remote][DOWN2] = grid.sub_universe_rank - grid.offset_L;
				grid.my_domain.internal_info.half_marker = 1;
			} else if ((grid.sub_universe_rank - grid.offset_L) >= (grid.n / 2)) {
				grid.nbrs[remote][DOWN1] = (grid.offset_R + (grid.n - 1)) - (grid.sub_universe_rank - grid.offset_L);
				grid.nbrs[remote][DOWN2] = (grid.offset_R + (grid.n - 1)) - (grid.sub_universe_rank - grid.offset_L);
				grid.flip_array[DOWN1] = 1;
				grid.flip_array[DOWN2] = 1;
				grid.my_domain.internal_info.half_marker = 2;
			}
		}
		//For Right daughter branch edge
		else if ((grid.sub_universe_rank >= grid.offset_R) && (grid.sub_universe_rank < (grid.offset_R + grid.n))) {
			grid.my_domain.internal_info.boundary_tag = 'I';
			if ((grid.sub_universe_rank - grid.offset_R) < (grid.n / 2)) {
				grid.nbrs[remote][DOWN1] = (grid.offset_L + (grid.n - 1)) - (grid.sub_universe_rank - grid.offset_R);
				grid.nbrs[remote][DOWN2] = (grid.offset_L + (grid.n - 1)) - (grid.sub_universe_rank - grid.offset_R);
				grid.flip_array[DOWN1] = 1;
				grid.flip_array[DOWN2] = 1;
				grid.my_domain.internal_info.half_marker = 2;
			} else if ((grid.sub_universe_rank - grid.offset_R) >= (grid.n / 2)) {
				grid.nbrs[remote][DOWN1] = grid.sub_universe_rank - grid.offset_R;
				grid.nbrs[remote][DOWN2] = grid.sub_universe_rank - grid.offset_R;
				grid.my_domain.internal_info.half_marker = 1;
			}
		}
	}

	//If number of processors in circumferential dimension are ODD

	if (grid.scheme != 0) {
		//The parent artery edge
		if ((grid.sub_universe_rank >= 0) && (grid.sub_universe_rank < grid.n)) {
			grid.my_domain.internal_info.boundary_tag = 'I';
			if ((grid.sub_universe_rank - grid.offset_P) < ((grid.n - 1) / 2)) {
				grid.nbrs[remote][UP1] = grid.offset_L + (grid.sub_universe_rank - grid.offset_P);
				grid.nbrs[remote][UP2] = grid.offset_L + (grid.sub_universe_rank - grid.offset_P);
			} else if ((grid.sub_universe_rank - grid.offset_P) > ((grid.n - 1) / 2)) {
				grid.nbrs[remote][UP1] = grid.offset_R + (grid.sub_universe_rank - grid.offset_P);
				grid.nbrs[remote][UP2] = grid.offset_R + (grid.sub_universe_rank - grid.offset_P);
			} else if ((grid.sub_universe_rank - grid.offset_P) == ((grid.n - 1) / 2)) {
				grid.nbrs[remote][UP1] = grid.offset_L + (grid.sub_universe_rank - grid.offset_P);
				grid.nbrs[remote][UP2] = grid.offset_R + (grid.sub_universe_rank - grid.offset_P);
			}
		}
		//The left daughter artery edge
		else if ((grid.sub_universe_rank >= grid.offset_L) && (grid.sub_universe_rank < grid.offset_L + grid.n)) {
			grid.my_domain.internal_info.boundary_tag = 'I';
			if ((grid.sub_universe_rank - grid.offset_L) < ((grid.n - 1) / 2)) {
				grid.nbrs[remote][DOWN1] = (grid.sub_universe_rank - grid.offset_L);
				grid.nbrs[remote][DOWN2] = (grid.sub_universe_rank - grid.offset_L);
				grid.my_domain.internal_info.half_marker = 1;
			} else if ((grid.sub_universe_rank - grid.offset_L) > ((grid.n - 1) / 2)) {
				grid.nbrs[remote][DOWN1] = (grid.offset_R + (grid.n - 1)) - (grid.sub_universe_rank - grid.offset_L);
				grid.nbrs[remote][DOWN2] = (grid.offset_R + (grid.n - 1)) - (grid.sub_universe_rank - grid.offset_L);
				grid.flip_array[DOWN1] = 1;
				grid.flip_array[DOWN2] = 1;
				grid.my_domain.internal_info.half_marker = 2;
			} else if ((grid.sub_universe_rank - grid.offset_L) == ((grid.n - 1) / 2)) {
				grid.nbrs[remote][DOWN1] = (grid.sub_universe_rank - grid.offset_L);
				grid.nbrs[remote][DOWN2] = (grid.offset_R + (grid.n - 1)) - (grid.sub_universe_rank - grid.offset_L);
				grid.flip_array[DOWN1] = 0;
				grid.flip_array[DOWN2] = 1;
				grid.my_domain.internal_info.half_marker = 3;
			}
		}
		//The right daughter artery edge
		else if ((grid.sub_universe_rank >= grid.offset_R) && (grid.sub_universe_rank < grid.offset_R + grid.n)) {
			grid.my_domain.internal_info.boundary_tag = 'I';
			if ((grid.sub_universe_rank - grid.offset_R) < ((grid.n - 1) / 2)) {
				grid.nbrs[remote][DOWN1] = (grid.offset_L + (grid.n - 1)) - (grid.sub_universe_rank - grid.offset_R);
				grid.nbrs[remote][DOWN2] = (grid.offset_L + (grid.n - 1)) - (grid.sub_universe_rank - grid.offset_R);
				grid.flip_array[DOWN1] = 1;
				grid.flip_array[DOWN2] = 1;
				grid.my_domain.internal_info.half_marker = 2;
			} else if ((grid.sub_universe_rank - grid.offset_R) > ((grid.n - 1) / 2)) {
				grid.nbrs[remote][DOWN1] = grid.sub_universe_rank - grid.offset_R;
				grid.nbrs[remote][DOWN2] = grid.sub_universe_rank - grid.offset_R;
				grid.my_domain.internal_info.half_marker = 1;
			} else if ((grid.sub_universe_rank - grid.offset_R) == ((grid.n - 1) / 2)) {
				grid.nbrs[remote][DOWN1] = (grid.offset_L + (grid.n - 1)) - (grid.sub_universe_rank - grid.offset_R);
				grid.nbrs[remote][DOWN2] = grid.sub_universe_rank - grid.offset_R;
				grid.flip_array[DOWN1] = 1;
				grid.flip_array[DOWN2] = 0;
				grid.my_domain.internal_info.half_marker = 3;
			}
		}
	}

	//If I am a parent branch in my domain
	if (grid.branch_tag == P) {
		//if a parent domain exists for me
		if (grid.my_domain.parent.domain_index >= 0) {
			//if I am a bottom row in my m x n cart grid
			if ((grid.rank >= ((grid.m - 1) * grid.n)) && (grid.rank <= (grid.m * grid.n - 1))) {
				int stride = grid.rank - ((grid.m - 1) * grid.n);
				grid.nbrs[remote][DOWN1] = grid.my_domain.parent.domain_start + stride;
				grid.nbrs[remote][DOWN2] = grid.my_domain.parent.domain_start + stride;
			}
		}
	}
	//If I am a left daughter branch in my domain
	else if (grid.branch_tag == L) {
		//if a child exists from me
		if (grid.my_domain.left_child.domain_index >= 0) {
			//if I am top row in my m x n cart grid
			if ((grid.rank >= 0) && (grid.rank <= (grid.n - 1))) {
				int stride = grid.rank;
				grid.nbrs[remote][UP1] = grid.my_domain.left_child.domain_start + stride;
				grid.nbrs[remote][UP2] = grid.my_domain.left_child.domain_start + stride;
			}
		}
	}
	//If I am a right daughter branch in my domain
	else if (grid.branch_tag == R) {
		//if a child exists from me
		if (grid.my_domain.right_child.domain_index >= 0) {
			//if I am top row in my m x n cart grid
			if ((grid.rank >= 0) && (grid.rank <= (grid.n - 1))) {
				int stride = grid.rank;
				grid.nbrs[remote][UP1] = grid.my_domain.right_child.domain_start + stride;
				grid.nbrs[remote][UP2] = grid.my_domain.right_child.domain_start + stride;
			}
		}
	}
	return grid;
}			// end of make_bifurcation

grid_parms make_straight_segment(grid_parms grid) {
	//Since there no branch, all processors have same color.

	grid.color = 0;
	grid.key = 0;

	check_flag(MPI_Comm_split(grid.sub_universe, grid.color, grid.key, &grid.split_comm), "Comm-split failed");

	///Global variables that are to be read by each processor
	int ndims, nbrs[4], dims[2], periodic[2], reorder = 0, coords[2];
	ndims = 2;
	dims[0] = grid.m;
	dims[1] = grid.n;
	periodic[0] = 0;
	periodic[1] = 1;
	reorder = 0;

	check_flag(MPI_Cart_create(grid.split_comm, ndims, dims, periodic, reorder, &grid.cart_comm), "failed at cart create");
	check_flag(MPI_Comm_rank(grid.cart_comm, &grid.rank), "failed at comm rank");
	check_flag(MPI_Comm_size(grid.cart_comm, &grid.tasks), "failed at cart comm tasks");
	check_flag(MPI_Cart_coords(grid.cart_comm, grid.rank, ndims, grid.coords), "failed at cart coords");

	check_flag(MPI_Cart_shift(grid.cart_comm, 0, 1, &grid.nbrs[local][UP], &grid.nbrs[local][DOWN]), "failed at cart shift up down");
	check_flag(MPI_Cart_shift(grid.cart_comm, 1, 1, &grid.nbrs[local][LEFT], &grid.nbrs[local][RIGHT]), "failed at cart left right");

	//Label the ranks on the subdomain edges of at STRAIGHT SEGMENT as top (T) or bottom boundary (B) or none (N).
	for (int i = 0; i < (grid.m * grid.n); i++) {
		if ((grid.rank >= ((grid.m - 1) * grid.n)) && (grid.rank <= (grid.m * grid.n - 1))) {
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
		if ((grid.rank >= ((grid.m - 1) * grid.n)) && (grid.rank <= (grid.m * grid.n - 1))) {
			int stride = grid.rank - ((grid.m - 1) * grid.n);
			grid.nbrs[remote][DOWN1] = grid.my_domain.parent.domain_start + stride;
			grid.nbrs[remote][DOWN2] = grid.my_domain.parent.domain_start + stride;
		}
	}
	//if a child exists from me
	if (grid.my_domain.left_child.domain_index >= 0) {
		//if I am top row in my m x n cart grid
		if ((grid.rank >= 0) && (grid.rank <= (grid.n - 1))) {
			int stride = grid.rank;
			grid.nbrs[remote][UP1] = grid.my_domain.left_child.domain_start + stride;
			grid.nbrs[remote][UP2] = grid.my_domain.left_child.domain_start + stride;
		}
	}
	return grid;
}			//end of make_straight_segment()

/**************** Update global subdomain information ***************/
grid_parms update_global_subdomain_information(grid_parms grid, int num_subdomains, int** domains)
/********************************************************************/
{
	grid.global_domain_info.num_subdomains = num_subdomains;

	grid.global_domain_info.m = (int*) checked_malloc((grid.global_domain_info.num_subdomains * sizeof(int)), "aa");
	grid.global_domain_info.n = (int*) checked_malloc((grid.global_domain_info.num_subdomains * sizeof(int)), "bb");

	grid.global_domain_info.list_type_subdomains = (int*) checked_malloc((grid.global_domain_info.num_subdomains * sizeof(int)), "cc");
	grid.global_domain_info.list_num_ec_axially_per_domain = (int*) checked_malloc((grid.global_domain_info.num_subdomains * sizeof(int)), "dd");
	grid.global_domain_info.list_domain_z_coord_index = (double**) checked_malloc((grid.global_domain_info.num_subdomains * sizeof(double*)), "ee");
	for (int i = 0; i < grid.global_domain_info.num_subdomains; i++) {
		grid.global_domain_info.list_domain_z_coord_index[i] = (double*) checked_malloc((4 * sizeof(double)), "ff");
	}

	grid.my_domain.z_offset_start = 0.0;
	grid.my_domain.z_offset_end = 0.0;
	double theta = 3.1415 / 4;
	grid = z_coord_exchange(grid, theta);
	return (grid);
}
/********************************************************/
grid_parms z_coord_exchange(grid_parms grid, double theta)
/********************************************************/
{
	int tag = 101;
	MPI_Status status;
	int root;
	///If  there is no parent to me
	if (grid.my_domain.parent.domain_index == none) {
		root = 0;
		grid = my_z_offset(grid, theta);
		check_flag((MPI_Bcast(&grid.my_domain.z_offset_end, 1, MPI_DOUBLE, root, grid.cart_comm)), "broadcast Case no parent");

		if (grid.my_domain.internal_info.domain_type == STRSEG) {
			if (grid.my_domain.internal_info.boundary_tag == 'T') {
				check_flag((MPI_Send(&grid.my_domain.z_offset_end, 1, MPI_DOUBLE, grid.nbrs[remote][UP1], tag, grid.universe)),
						"z_info_no_parent_type_STRSEG");
			}
		} else if (grid.my_domain.internal_info.domain_type == BIF) {
			if (grid.branch_tag == P) {
				if (grid.my_domain.internal_info.boundary_tag == 'I') {
					check_flag((MPI_Send(&grid.my_domain.z_offset_end, 1, MPI_DOUBLE, grid.nbrs[remote][UP1], tag, grid.sub_universe)),
							"z_info_no_parent_type_BIF");
				}
			} else if (((grid.branch_tag == L) || (grid.branch_tag == R))) {
				if ((grid.my_domain.internal_info.boundary_tag == 'I') && (grid.my_domain.internal_info.half_marker == 1)) {
					check_flag((MPI_Recv(&grid.my_domain.z_offset_start, 1, MPI_DOUBLE, grid.nbrs[remote][DOWN1], tag, grid.sub_universe, &status)),
							"z_info_no_parent_type_BIF_daughter_branches");
				}
				if (grid.branch_tag == L) {
					root = (grid.m - 1) * grid.n;
				} else if (grid.branch_tag == R) {
					root = (grid.m * grid.n) - 1;
				}
				check_flag((MPI_Bcast(&grid.my_domain.z_offset_start, 1, MPI_DOUBLE, root, grid.cart_comm)), "broadcast Case no parent");
				grid = my_z_offset(grid, theta);
				if (grid.my_domain.internal_info.boundary_tag == 'T') {
					check_flag((MPI_Send(&grid.my_domain.z_offset_end, 1, MPI_DOUBLE, grid.nbrs[remote][UP1], tag, grid.universe)),
							"z_info_no_parent_type_BIF");
				}
			}
		}
	}			///If there is a parent domain do the following.
	else if (grid.my_domain.parent.domain_index > none) {
		///If my domain type is STRAIGHT SEGMENT do the following....
		if (grid.my_domain.internal_info.domain_type == STRSEG) {
			if (grid.my_domain.internal_info.boundary_tag == 'B') {
				check_flag((MPI_Recv(&grid.my_domain.z_offset_start, 1, MPI_DOUBLE, grid.nbrs[remote][DOWN1], tag, grid.universe, &status)),
						"z_info_parent_type_STRSEG");
			}
			root = (grid.m - 1) * grid.n;
			check_flag((MPI_Bcast(&grid.my_domain.z_offset_start, 1, MPI_DOUBLE, root, grid.cart_comm)), "broadcast Case with parent in STRSEG");
			grid = my_z_offset(grid, theta);
			if (grid.my_domain.internal_info.boundary_tag == 'T') {
				check_flag((MPI_Send(&grid.my_domain.z_offset_end, 1, MPI_DOUBLE, grid.nbrs[remote][UP1], tag, grid.universe)),
						"z_info_no_parent_type_BIF");
			}
		}
		/// Otherwise if my domain type is BIFURCATION do the following....
		else if (grid.my_domain.internal_info.domain_type == BIF) {
			if (grid.branch_tag == P) {
				if (grid.my_domain.internal_info.boundary_tag == 'B') {
					check_flag((MPI_Recv(&grid.my_domain.z_offset_start, 1, MPI_DOUBLE, grid.nbrs[remote][DOWN1], tag, grid.universe, &status)),
							"z_info_parent_type_BIF");
				}
				root = (grid.m - 1) * grid.n;

				check_flag((MPI_Bcast(&grid.my_domain.z_offset_start, 1, MPI_DOUBLE, root, grid.cart_comm)), "broadcast Case with parent in BIF");
				grid = my_z_offset(grid, theta);
				if (grid.my_domain.internal_info.boundary_tag == 'I') {
					check_flag((MPI_Send(&grid.my_domain.z_offset_end, 1, MPI_DOUBLE, grid.nbrs[remote][UP1], tag, grid.sub_universe)),
							"z_info_with_parent_type_BIF sending to daughter segments");
				}
			} else if ((grid.branch_tag == L) || (grid.branch_tag == R)) {
				if ((grid.my_domain.internal_info.boundary_tag == 'I') && (grid.my_domain.internal_info.half_marker == 1)) {
					check_flag((MPI_Recv(&grid.my_domain.z_offset_start, 1, MPI_DOUBLE, grid.nbrs[remote][DOWN1], tag, grid.sub_universe, &status)),
							"in type_BIF, z_info_recv from parent segment with parent domain");
				}
				if (grid.branch_tag == L) {
					root = (grid.m - 1) * grid.n;
				} else if (grid.branch_tag == R) {
					root = (grid.m * grid.n) - 1;
				}
				check_flag((MPI_Bcast(&grid.my_domain.z_offset_start, 1, MPI_DOUBLE, root, grid.cart_comm)), "broadcast Case with parent");
				grid = my_z_offset(grid, theta);
				if (grid.my_domain.internal_info.boundary_tag == 'T') {
					check_flag((MPI_Send(&grid.my_domain.z_offset_end, 1, MPI_DOUBLE, grid.nbrs[remote][UP1], tag, grid.universe)),
							"z_info_with_parent_type_BIF sending to daughter segments");
				}
			}
		}
	}

	double local_z_start, local_z_end;
	int rem1 = (grid.m - (int) (grid.rank / grid.n)) - 1, rem2 = grid.m - (int) (grid.rank / grid.n);

	grid.my_domain.local_z_end = grid.my_domain.z_offset_start
			+ rem1 * (grid.my_domain.z_offset_end - grid.my_domain.z_offset_start) / ((double) (grid.m));
	grid.my_domain.local_z_start = grid.my_domain.z_offset_start
			+ rem2 * (grid.my_domain.z_offset_end - grid.my_domain.z_offset_start) / ((double) (grid.m));
	return (grid);
}

/************************************************/
grid_parms my_z_offset(grid_parms grid, double theta)
/************************************************/
{
	if (grid.my_domain.internal_info.domain_type == STRSEG) {
		grid.my_domain.z_offset_end = grid.my_domain.z_offset_start + (double) (grid.m * grid.num_ec_axially) * 65e-6;
	} else if (grid.my_domain.internal_info.domain_type == BIF) {
		if (grid.branch_tag == P) {
			grid.my_domain.z_offset_end = grid.my_domain.z_offset_start + (double) (grid.m * grid.num_ec_axially) * 65e-6;
		} else if ((grid.branch_tag == L) || (grid.branch_tag == R)) {
			grid.my_domain.z_offset_end = grid.my_domain.z_offset_start + (double) (grid.m * grid.num_ec_axially) * 65e-6 * 0.707;
		}
	}
	return (grid);
}

/*************************************************/
IO_domain_info* make_io_domains(grid_parms* grid)
/*************************************************/
{
	IO_domain_info *my_IO_domain_info = (IO_domain_info*) checked_malloc((sizeof(IO_domain_info)),
			"allocation of IO_domain info structure failed in topology.cpp");
	if (grid->rank == 0) {
		my_IO_domain_info->my_IO_domain_color = WRITER_COLOR;
		my_IO_domain_info->my_IO_domain_key = WRITER_KEY;
	} else {
		my_IO_domain_info->my_IO_domain_color = COMPUTE_ONLY_COLOR;
		my_IO_domain_info->my_IO_domain_key = COMPUTE_ONLY_KEY;
	}

	MPI_Comm tmp_comm;

	/// Spliting MPI_COMM_WORLD into communicators aggregating Writers (i.e. all bridge processes located near IO nodes)
	/// and non-Writers
	check_flag(MPI_Comm_split(grid->universe, my_IO_domain_info->my_IO_domain_color, my_IO_domain_info->my_IO_domain_key, &tmp_comm),
			"error constructing IO_comm for writer nodes.");
	if (my_IO_domain_info->my_IO_domain_color == WRITER_COLOR) {
		check_flag(MPI_Comm_dup(tmp_comm, &my_IO_domain_info->writer_comm),
				"Comm_duplicate failed while processing the name change of the io-handlers communicator.");
	}
	check_flag(MPI_Comm_free(&tmp_comm), "Error freeing MPI_Comm tmp_comm that was created from MPI_COMM_WORLD.");

	///Reveal and record rank information of members of writer domains.
	for (int i = 0; i < grid->numtasks; i++) {
		my_IO_domain_info->writer_rank = -1;
		my_IO_domain_info->writer_tasks = -1;
	}

	if (my_IO_domain_info->my_IO_domain_color == WRITER_COLOR) {
		check_flag(MPI_Comm_rank(my_IO_domain_info->writer_comm, &my_IO_domain_info->writer_rank),
				"Error revealing the rank of writer-domain members.");
		check_flag(MPI_Comm_size(my_IO_domain_info->writer_comm, &my_IO_domain_info->writer_tasks),
				"Error revealing the size of writer-domain members.");
	}

	return (my_IO_domain_info);
}


