#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"
#include "rksuite.h"


using namespace std;
/*******************************************************************************************/
int couplingParms(int CASE,conductance* cpl_cef)
/*******************************************************************************************/
{
	if (CASE == 1) {
		cpl_cef->Vm_hm_smc = 1000.00;
		cpl_cef->Vm_hm_ec = 1000.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.05;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.00;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.00;
		cpl_cef->Ca_ht_ec = 0.00;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	} else if (CASE == 2) {
		cpl_cef->Vm_hm_smc = 1000.00;
		cpl_cef->Vm_hm_ec = 1000.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.05;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.00;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.05;
		cpl_cef->Ca_ht_ec = 0.05;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	} else if (CASE == 3) {
		cpl_cef->Vm_hm_smc = 1000.00;
		cpl_cef->Vm_hm_ec = 1000.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.05;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.05;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.05;
		cpl_cef->Ca_ht_ec = 0.05;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	} else if (CASE == 4) {
		cpl_cef->Vm_hm_smc = 1000.00;
		cpl_cef->Vm_hm_ec = 0.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.00;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.05;

		cpl_cef->Vm_ht_smc = 0.0;
		cpl_cef->Vm_ht_ec = 0.0;

		cpl_cef->Ca_ht_smc = 0.0;
		cpl_cef->Ca_ht_ec = 0.0;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	} else if (CASE == 5)// Simulating for experiments suggested by Dr. James Kazloski (IBM Watson Centre).
			// The homocellular Ca coupling between SMCs is changed to investigate the effects of
			//strreangth on coupling on the propagation speed of the spatial waves.
			{
		cpl_cef->Vm_hm_smc = 0.00;
		cpl_cef->Vm_hm_ec = 0.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.05;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.00;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.00;
		cpl_cef->Ca_ht_ec = 0.00;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	} else if (CASE == 6)// Simulating for experiments suggested by Dr. James Kazloski (IBM Watson Centre)
			{
		cpl_cef->Vm_hm_smc = 4000.00;
		cpl_cef->Vm_hm_ec = 1000.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.05;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.00;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.00;
		cpl_cef->Ca_ht_ec = 0.00;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	} else if (CASE == 7)// Simulating for experiments suggested by Dr. James Kazloski (IBM Watson Centre)
			{
		cpl_cef->Vm_hm_smc = 6000.00;
		cpl_cef->Vm_hm_ec = 1000.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.05;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.00;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.00;
		cpl_cef->Ca_ht_ec = 0.00;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	} else if (CASE == 8) {
		cpl_cef->Vm_hm_smc = 0.00;
		cpl_cef->Vm_hm_ec = 0.00;

		cpl_cef->Ca_hm_smc = 0.0;
		cpl_cef->Ca_hm_ec = 0.0;

		cpl_cef->IP3_hm_smc = 0.0;
		cpl_cef->IP3_hm_ec = 0.0;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.00;
		cpl_cef->Ca_ht_ec = 0.00;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	}

	return 0;
}

/*******************************************************************************************/
void Initialize_koeingsberger_smc(grid_parms grid, double y[], celltype1** smc)
/*******************************************************************************************/
{
	int k = 0, offset;

	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;
			y[k + ((j - 1) * grid.neq_smc) + smc_Ca] = 0.01;
			y[k + ((j - 1) * grid.neq_smc) + smc_SR] = 0.01;
			y[k + ((j - 1) * grid.neq_smc) + smc_Vm] = 0.01;
			y[k + ((j - 1) * grid.neq_smc) + smc_w] = 0.01;
			y[k + ((j - 1) * grid.neq_smc) + smc_IP3] = 0.01;
			/*if ((grid.branch_tag = P) || (grid.branch_tag = L)
					|| (grid.branch_tag = R)) {
				y[k + ((j - 1) * grid.neq_smc) + smc_Ca] =
						(double) (grid.universal_rank) + ((double)(i)*0.01) + ((double)(j)*0.0001) + 0.0000001;
				y[k + ((j - 1) * grid.neq_smc) + smc_SR] =
						(double) (grid.universal_rank) + ((double)(i)*0.01) + ((double)(j)*0.0001) + 0.0000002;
				y[k + ((j - 1) * grid.neq_smc) + smc_Vm] =
						(double) (grid.universal_rank) + ((double)(i)*0.01) + ((double)(j)*0.0001) + 0.0000003;
				y[k + ((j - 1) * grid.neq_smc) + smc_w] =
						(double) (grid.universal_rank) + ((double)(i)*0.01) + ((double)(j)*0.0001) + 0.0000004;
				y[k + ((j - 1) * grid.neq_smc) + smc_IP3] =
						(double) (grid.universal_rank) + ((double)(i)*0.01) + ((double)(j)*0.0001) + 0.0000005;
			//}*/
		}
	}

	for (int i = 0; i < (grid.num_smc_circumferentially + grid.num_ghost_cells);
			i++) {
		for (int j = 0; j < (grid.num_smc_axially + grid.num_ghost_cells);
				j++) {
			smc[i][j].p[smc_Ca] = 0.0;			//0.23464;
			smc[i][j].p[smc_SR] = 0.0;			//1.2816;
			smc[i][j].p[smc_Vm] = 0.0;			//-27.031;
			smc[i][j].p[smc_w] = 0.0;			//0.29697;
			smc[i][j].p[smc_IP3] = 0.0;			//0.099996;

			for (int k = 1; k <= grid.num_fluxes_smc; k++) {
				smc[i][j].A[k - 1] = 0.0;
			}
			for (int k = 1; k <= grid.num_coupling_species_smc; k++) {
				smc[i][j].B[k - 1] = 0.0;
				smc[i][j].C[k - 1] = 0.0;
			}
		}
	}
}
/*******************************************************************************************/
void Initialize_koeingsberger_ec(grid_parms grid, double y[], celltype2** ec)
/*******************************************************************************************/
{
int	k,offset = (grid.neq_smc * grid.num_smc_circumferentially
					* grid.num_smc_axially);

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			if (i > 1)
				k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
				k = offset + 0;
			y[k + ((j - 1) * grid.neq_ec) + ec_Ca] = 0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_SR] = 0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_Vm] = 0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_IP3] = 0.01;
			/*if ((grid.branch_tag = P) || (grid.branch_tag = L)
					|| (grid.branch_tag = R)) {
				y[k + ((j - 1) * grid.neq_ec) + ec_Ca] =
						(double) (grid.universal_rank) + ((double)(i)*0.01) + ((double)(j)*0.0001) + 0.0000006;
				y[k + ((j - 1) * grid.neq_ec) + ec_SR] =
						(double) (grid.universal_rank) + ((double)(i)*0.01) + ((double)(j)*0.0001) + 0.0000007;
				y[k + ((j - 1) * grid.neq_ec) + ec_Vm] =
						(double) (grid.universal_rank) + ((double)(i)*0.01) + ((double)(j)*0.0001) + 0.0000008;
				y[k + ((j - 1) * grid.neq_ec) + ec_IP3] =
						(double) (grid.universal_rank) + ((double)(i)*0.01) + ((double)(j)*0.0001) + 0.0000009;
			//}*/
		}
	}
	for (int i = 0; i < (grid.num_ec_circumferentially + grid.num_ghost_cells);
			i++) {
		for (int j = 0; j < (grid.num_ec_axially + grid.num_ghost_cells);
				j++) {
			ec[i][j].q[ec_Ca] = 0.0;//0.12988;
			ec[i][j].q[ec_SR] = 0.0; //0.27952;
			ec[i][j].q[ec_Vm] = 0.0;//-38.713;
			ec[i][j].q[ec_IP3] = 0.0;//0.099996;

			for (int k = 1; k <= grid.num_fluxes_ec; k++) {
				ec[i][j].A[k - 1] = 0.0;
			}
			for (int k = 1; k <= grid.num_coupling_species_ec; k++) {
				ec[i][j].B[k - 1] = 0.0;
				ec[i][j].C[k - 1] = 0.0;
			}
		}
	}
}
/*******************************************************************************************/
void map_solver_to_cells(grid_parms grid,double y[], celltype1** smc, celltype2** ec)
/*******************************************************************************************/
{
	///Mapping state vectors of each cell (all except the ghost cells) to the corresponding locations on the solver's solution array y[].
	int k = 0, offset;
/*
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;
			smc[i][j].p = &y[k + ((j - 1) * grid.neq_smc) + 0];
		}
	}
	offset = (grid.neq_smc * grid.num_smc_circumferentially
			* grid.num_smc_axially);

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			if (i > 1)
				k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
				k = offset + 0;
			ec[i][j].q = &y[k + ((j - 1) * grid.neq_ec) + 0];
		}
	}*/

	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;
			smc[i][j].p[smc_Ca] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca];
			smc[i][j].p[smc_SR] = y[k + ((j - 1) * grid.neq_smc) + smc_SR];
			smc[i][j].p[smc_Vm] = y[k + ((j - 1) * grid.neq_smc) + smc_Vm];
			smc[i][j].p[smc_w] = y[k + ((j - 1) * grid.neq_smc) + smc_w];
			smc[i][j].p[smc_IP3] = y[k + ((j - 1) * grid.neq_smc) + smc_IP3];
		}
	}
	offset = (grid.neq_smc * grid.num_smc_circumferentially
			* grid.num_smc_axially);

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			if (i > 1)
				k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
				k = offset + 0;
			ec[i][j].q[ec_Ca] = y[k + ((j - 1) * grid.neq_ec) + ec_Ca];
			ec[i][j].q[ec_SR] = y[k + ((j - 1) * grid.neq_ec) + ec_SR];
			ec[i][j].q[ec_Vm] = y[k + ((j - 1) * grid.neq_ec) + ec_Vm];
			ec[i][j].q[ec_IP3] = y[k + ((j - 1) * grid.neq_ec) + ec_IP3];
		}
	}
}
/*******************************************************************************************/
void map_GhostCells_to_cells(celltype1** smc, celltype2** ec, grid_parms grid)
/*******************************************************************************************/
{
	///Allocating arrays of appropriate lengths for holding values of ghost cell variable to which the relevant smc[i][j] and ec[i][j] members will point to.
		///The rest of the computational domain is to point to the solution vector y of the solver.
		double **ghost_cells_smc_circumferentially,
				**ghost_cells_ec_circumferentially,
				**ghost_cells_smc_axially,
				**ghost_cells_ec_axially;

		ghost_cells_smc_circumferentially = (double**) checked_malloc(
				sizeof(double*), stdout, "ghost cell in circ direction for smc");
		for (int i=0; i<2; i++){
			ghost_cells_smc_circumferentially[i] = (double*) checked_malloc(
				(grid.num_ghost_cells + grid.num_smc_circumferentially)
						* grid.neq_smc * sizeof(double), stdout,
						"ghost cell in circ direction for smc");
		}

		ghost_cells_ec_circumferentially = (double**) checked_malloc(
					sizeof(double*), stdout, "ghost cell in circ direction for ec");
		for (int i=0; i<2; i++){
				ghost_cells_ec_circumferentially[i] = (double*) checked_malloc(
					(grid.num_ghost_cells + grid.num_ec_circumferentially)
							* grid.neq_ec * sizeof(double), stdout,
							"ghost cell in circ direction for smc");
			}

		ghost_cells_smc_axially = (double**) checked_malloc(
					sizeof(double*), stdout, "ghost cell in axial direction for smc");
		for (int i=0; i<2; i++){
			ghost_cells_smc_axially[i] = (double*) checked_malloc(
				(grid.num_ghost_cells + grid.num_smc_axially)
						* grid.neq_smc * sizeof(double), stdout,
						"ghost cell in axial direction for smc");
		}

		ghost_cells_ec_axially = (double**) checked_malloc(
					sizeof(double*), stdout, "ghost cell in axial direction for ec");
		for (int i=0; i<2; i++){
				ghost_cells_ec_axially[i] = (double*) checked_malloc(
					(grid.num_ghost_cells + grid.num_ec_axially)
							* grid.neq_ec * sizeof(double), stdout,
							"ghost cell in axial direction for smc");
			}
/*
	for (int i = 0; i < (grid.num_ghost_cells + grid.num_smc_circumferentially);
			i++) {
		int j;
		j = 0;
			smc[i][j].p = &ghost_cells_smc_circumferentially[0][i * grid.neq_smc + 0];
		j = grid.num_ghost_cells + grid.num_smc_axially - 1;
			smc[i][j].p = &ghost_cells_smc_circumferentially[1][i * grid.neq_smc + 0];
	}
	for (int j = 0; j < (grid.num_ghost_cells + grid.num_smc_axially); j++) {
		int i;
		i = 0;
			smc[i][j].p = &ghost_cells_smc_axially[0][j * grid.neq_smc + 0];
		i = grid.num_ghost_cells + grid.num_smc_circumferentially - 1;
			smc[i][j].p = &ghost_cells_smc_axially[1][j * grid.neq_smc + 0];
	}

	for (int i = 0; i < (grid.num_ghost_cells + grid.num_ec_circumferentially);
			i++) {
		int j;
		j = 0;
			ec[i][j].q = &ghost_cells_ec_circumferentially[0][i * grid.neq_ec + 0];
		j = grid.num_ghost_cells + grid.num_ec_axially - 1;
			ec[i][j].q = &ghost_cells_ec_circumferentially[1][i * grid.neq_ec + 0];
	}
	for (int j = 0; j < (grid.num_ghost_cells + grid.num_ec_axially); j++) {
		int i;
		i = 0;
			ec[i][j].q = &ghost_cells_ec_axially[0][j * grid.neq_ec + 0];
		i = grid.num_ghost_cells + grid.num_ec_circumferentially - 1;
			ec[i][j].q = &ghost_cells_ec_axially[1][j * grid.neq_ec + 0];
	}
*/
/*	for (int i = 0; i < (grid.num_ghost_cells + grid.num_smc_circumferentially);
			i++) {
		int j;
		j = 0;
		smc[i][j].p[0]
				= ghost_cells_smc_circumferentially[0][i * grid.neq_smc + 0];
		smc[i][j].p[1]
				= ghost_cells_smc_circumferentially[0][i * grid.neq_smc + 1];
		smc[i][j].p[2]
				= ghost_cells_smc_circumferentially[0][i * grid.neq_smc + 2];
		smc[i][j].p[3]
				= ghost_cells_smc_circumferentially[0][i * grid.neq_smc + 3];
		smc[i][j].p[3]
				= ghost_cells_smc_circumferentially[0][i * grid.neq_smc + 4];
		j = grid.num_ghost_cells + grid.num_smc_axially - 1;
		smc[i][j].p[0]
				= ghost_cells_smc_circumferentially[1][i * grid.neq_smc + 0];
		smc[i][j].p[1]
				= ghost_cells_smc_circumferentially[1][i * grid.neq_smc + 1];
		smc[i][j].p[2]
				= ghost_cells_smc_circumferentially[1][i * grid.neq_smc + 2];
		smc[i][j].p[3]
				= ghost_cells_smc_circumferentially[1][i * grid.neq_smc + 3];
		smc[i][j].p[4]
				= ghost_cells_smc_circumferentially[1][i * grid.neq_smc + 4];
	}
	for (int j = 0; j < (grid.num_ghost_cells + grid.num_smc_axially); j++) {
		int i;
		i = 0;
		smc[i][j].p[0] = ghost_cells_smc_axially[0][j * grid.neq_smc + 0];
		smc[i][j].p[1] = ghost_cells_smc_axially[0][j * grid.neq_smc + 1];
		smc[i][j].p[2] = ghost_cells_smc_axially[0][j * grid.neq_smc + 2];
		smc[i][j].p[3] = ghost_cells_smc_axially[0][j * grid.neq_smc + 3];
		smc[i][j].p[4] = ghost_cells_smc_axially[0][j * grid.neq_smc + 4];
		i = grid.num_ghost_cells + grid.num_smc_circumferentially - 1;
		smc[i][j].p[0] = ghost_cells_smc_axially[1][j * grid.neq_smc + 0];
		smc[i][j].p[1] = ghost_cells_smc_axially[1][j * grid.neq_smc + 1];
		smc[i][j].p[2] = ghost_cells_smc_axially[1][j * grid.neq_smc + 2];
		smc[i][j].p[3] = ghost_cells_smc_axially[1][j * grid.neq_smc + 3];
		smc[i][j].p[4] = ghost_cells_smc_axially[1][j * grid.neq_smc + 4];

	}

	for (int i = 0; i < (grid.num_ghost_cells + grid.num_ec_circumferentially);
			i++) {
		int j;
		j = 0;
		ec[i][j].q[0]
				= ghost_cells_ec_circumferentially[0][i * grid.neq_ec + 0];
		ec[i][j].q[1]
				= ghost_cells_ec_circumferentially[0][i * grid.neq_ec + 1];
		ec[i][j].q[2]
				= ghost_cells_ec_circumferentially[0][i * grid.neq_ec + 2];
		ec[i][j].q[3]
				= ghost_cells_ec_circumferentially[0][i * grid.neq_ec + 3];
		j = grid.num_ghost_cells + grid.num_ec_axially - 1;
		ec[i][j].q[0]
				= ghost_cells_ec_circumferentially[0][i * grid.neq_ec + 0];
		ec[i][j].q[1]
				= ghost_cells_ec_circumferentially[0][i * grid.neq_ec + 1];
		ec[i][j].q[2]
				= ghost_cells_ec_circumferentially[0][i * grid.neq_ec + 2];
		ec[i][j].q[3]
				= ghost_cells_ec_circumferentially[0][i * grid.neq_ec + 3];
	}
	for (int j = 0; j < (grid.num_ghost_cells + grid.num_ec_axially); j++) {
		int i;
		i = 0;
		ec[i][j].q[0] = ghost_cells_ec_axially[0][j * grid.neq_ec + 0];
		ec[i][j].q[1] = ghost_cells_ec_axially[0][j * grid.neq_ec + 1];
		ec[i][j].q[2] = ghost_cells_ec_axially[0][j * grid.neq_ec + 2];
		ec[i][j].q[3] = ghost_cells_ec_axially[0][j * grid.neq_ec + 3];
		i = grid.num_ghost_cells + grid.num_ec_circumferentially - 1;
		ec[i][j].q[0] = ghost_cells_ec_axially[0][j * grid.neq_ec + 0];
		ec[i][j].q[1] = ghost_cells_ec_axially[0][j * grid.neq_ec + 1];
		ec[i][j].q[2] = ghost_cells_ec_axially[0][j * grid.neq_ec + 2];
		ec[i][j].q[3] = ghost_cells_ec_axially[0][j * grid.neq_ec + 3];
	}*/
}

/*******************************************************************************************/
void single_cell(double t, double y[], grid_parms grid,celltype1** smc, celltype2** ec)
/*******************************************************************************************/
{
	/* Constants for homogenically coupled SMCs*/
	const double Fi = 0.23;
	const double Kri = 1.00;
	const double GCai = 0.00129;
	const double vCa1 = 100.00;
	const double vCa2 = -24.00;
	const double RCai = 8.50;
	const double GNaCai = 0.00316;
	const double cNaCai = 0.5;
	const double vNaCai = -30.00;
	const double Bi = 2.025;
	const double cbi = 1.0;
	const double CICRi = 55.00;
	const double sci = 2.0;
	const double cci = 0.9;
	const double Di = 0.24;
	const double vdi = -100.00;
	const double Rdi = 250.00;
	const double Li = 0.025;
	const double gama = 1970.00;
	const double FNaK = 0.0432;
	const double GCli = 0.00134;
	const double vCl = -25.0;
	const double GKi = 0.0046;
	const double vKi = -94.00;
	const double lambda = 45.00;
	const double cwi = 0.0;
	const double beta = 0.13;
	const double vCa3 = -27.0;
	const double RKi = 12.0;
	const double ki = 0.1;
	/* Constants for homogenically coupled ECs*/
	const double Fj = 0.23;
	const double Krj = 1.00;
	const double Bj = 0.5;
	const double cbj = 1.0;
	const double CICRj = 5.0;
	const double scj = 2.0;
	const double ccj = 0.9;
	const double Dj = 0.24;
	const double Lj = 0.025;
	const double kj = 0.1;
	const double Gcatj = 0.66 * 1e-3;
	const double ECa = 50.00;
	const double m3cat = -0.18;			//-6.18;
	const double m4cat = 0.37;
	const double J0j = 0.029;
	const double Cmj = 25.8;
	const double Gtot = 6927;
	const double vKj = -80.0;
	const double a1j = 53.3;
	const double a2j = 53.3;
	const double bj = -80.8;
	const double c1j = -0.4; //-6.4;
	const double m3b = 1.32e-3;
	const double m4b = 0.30;
	const double m3s = -0.28;
	const double m4s = 0.389;
	const double GRj = 955;
	const double vrestj = -31.10;

	/*Intracellular calcium buffering*/
	const double k6 = 100.00;
	const double k7 = 300.00;
	const double BT = 120.00;

	int i, j, k;

//EVALUATING SINGLE CELL FLUXES :::::::For SMC::::::::
	for (i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (j = 1; j <= grid.num_smc_axially; j++) {

			//JIP3
			smc[i][j].A[J_IP3] = (Fi * P2(smc[i][j].p[smc_IP3]))
					/ (P2(Kri) + P2(smc[i][j].p[smc_IP3]));
			//JSRuptake
			smc[i][j].A[J_SERCA] = (Bi * P2(smc[i][j].p[smc_Ca]))
					/ (P2(smc[i][j].p[smc_Ca]) + P2(cbi));
			//Jcicr
			smc[i][j].A[J_CICR] =
					(CICRi
							* (P2(smc[i][j].p[smc_SR])
									* P4(smc[i][j].p[smc_Ca])))
							/ ((P2(sci) + P2(smc[i][j].p[smc_SR]))* ( P4(cci) + P4(smc[i][j].p[smc_Ca]) ) );
			//Jextrusion
			smc[i][j].A[J_Extrusion] = Di * smc[i][j].p[smc_Ca]
					* (1 + ((smc[i][j].p[smc_Vm] - vdi) / Rdi));

			//Jleak
			smc[i][j].A[J_Leak] = Li * smc[i][j].p[smc_SR];
			//Jvocc
			smc[i][j].A[J_VOCC] = GCai * (smc[i][j].p[smc_Vm] - vCa1)
					/ (1
							+ ((double) (exp(
									((-1) * (smc[i][j].p[smc_Vm] - vCa2))
											/ RCai))));
			//J Na/Ca
			smc[i][j].A[J_Na_Ca] = GNaCai * smc[i][j].p[smc_Ca]
					* (smc[i][j].p[smc_Vm] - vNaCai)
					/ (smc[i][j].p[smc_Ca] + cNaCai);
			//JNa/K
			smc[i][j].A[J_Na_K] = FNaK;
			//J Cl
			smc[i][j].A[J_Cl] = GCli * (smc[i][j].p[smc_Vm] - vCl);
			//JK
			smc[i][j].A[J_K] = GKi * smc[i][j].p[smc_w]
					* (smc[i][j].p[smc_Vm] - vKi);
			//Kactivation
			smc[i][j].A[K_activation] = P2(smc[i][j].p[smc_Ca] + cwi)
					/ (P2(smc[i][j].p[smc_Ca] + cwi)
							+ (beta
									* ((double) exp(
											(-1) * (smc[i][j].p[smc_Vm] - vCa3)
													/ RKi))));
			//Jdegradation
			smc[i][j].A[J_IP3_deg] = ki * smc[i][j].p[smc_IP3];

		} //end for j
	} //end for i

//EVALUATING SINGLE CELL FLUXES :::::::For EC::::::::

	for (i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (j = 1; j <= grid.num_ec_axially; j++) {
			//JIP3
			ec[i][j].A[J_IP3] = (Fj * P2(ec[i][j].q[ec_IP3]))
					/ (P2(Krj) + P2(ec[i][j].q[ec_IP3]));
			//JSRuptake
			ec[i][j].A[J_SERCA] = (Bj * P2(ec[i][j].q[ec_Ca]))
					/ (P2(ec[i][j].q[ec_Ca]) + P2(cbj));
			//Jcicr
			ec[i][j].A[J_CICR] =
					(CICRj * P2(ec[i][j].q[ec_SR]) * P4(ec[i][j].q[ec_Ca]))
							/ ((P2(ec[i][j].q[ec_SR]) + P2(scj))*(P4(ec[i][j].q[ec_Ca])+P4(ccj)));
			//Jextrusion
			ec[i][j].A[J_Extrusion] = Dj * ec[i][j].q[ec_Ca];
			//Jleak
			ec[i][j].A[J_Leak] = Lj * ec[i][j].q[ec_SR];
			//J_NonSelective Cation channels
			ec[i][j].A[J_NSC] = (Gcatj * (ECa - ec[i][j].q[ec_Vm]) * 0.5)
					* (1
							+ ((double) (tanh(
									(double) (((double) (log10(
											(double) (ec[i][j].q[ec_Ca])))
											- m3cat) / m4cat)))));
			//BK_channels
			ec[i][j].A[J_BK_Ca] =
					(0.4 / 2)
							* (1
									+ (double) (tanh(
											(double) ((((((double) (log10(
													(double) (ec[i][j].q[ec_Ca])))
													- c1j)
													* (ec[i][j].q[ec_Vm] - bj))
													- a1j)
													/ ((m3b
															* (P2(ec[i][j].q[ec_Vm]+(a2j*((double) (log10 ((double) (ec[i][j].q[ec_Ca])))-c1j))-bj)))+m4b))))));
			//SK_channels
			ec[i][j].A[J_SK_Ca] =
					(0.6 / 2)
							* (1
									+ (double) (tanh(
											(double) (((double) (log10(
													(double) (ec[i][j].q[ec_Ca])))
													- m3s) / m4s))));
			//Total K flux
			ec[i][j].A[J_Ktot] = Gtot * (ec[i][j].q[ec_Vm] - vKj)
					* (ec[i][j].A[J_BK_Ca] + ec[i][j].A[J_SK_Ca]);
			//Residual currents
			ec[i][j].A[J_Residual] = GRj * (ec[i][j].q[ec_Vm] - vrestj);
			//IP3 degradation
			ec[i][j].A[J_IP3_deg] = kj * ec[i][j].q[ec_IP3];
			//Grouping all other trivial Ca fluxes
			ec[i][j].A[J_trivial_Ca] = J0j;

		} //end for j
	} //end for i

}	//End of function single_cell

/*******************************************************************************************/
void coupling(double t,double y[], grid_parms grid,celltype1** smc, celltype2** ec,conductance cpl_cef)
/*******************************************************************************************/
{

	int i, j, k, l;

////******************** HOMOCELLULAR COUPLING *********************/
	for (i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (j = 1; j <= grid.num_smc_axially; j++) {
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
			smc[i][j].B[cpl_Ca] = -cpl_cef.Ca_hm_smc
					* ((smc[i][j].p[smc_Ca] - smc[i][up].p[smc_Ca])
							+ (smc[i][j].p[smc_Ca] - smc[i][down].p[smc_Ca])
							+ (smc[i][j].p[smc_Ca] - smc[left][j].p[smc_Ca])
							+ (smc[i][j].p[smc_Ca] - smc[right][j].p[smc_Ca]));
			smc[i][j].B[cpl_Vm] =
			-cpl_cef.Vm_hm_smc
			* ((smc[i][j].p[smc_Vm] - smc[i][up].p[smc_Vm])
					+ (smc[i][j].p[smc_Vm] - smc[i][down].p[smc_Vm])
					+ (smc[i][j].p[smc_Vm] - smc[left][j].p[smc_Vm])
					+ (smc[i][j].p[smc_Vm] - smc[right][j].p[smc_Vm]));
			smc[i][j].B[cpl_IP3] =
					-cpl_cef.IP3_hm_smc
							* ((smc[i][j].p[smc_IP3] - smc[i][up].p[smc_IP3])
									+ (smc[i][j].p[smc_IP3]
											- smc[i][down].p[smc_IP3])
									+ (smc[i][j].p[smc_IP3]
											- smc[left][j].p[smc_IP3])
									+ (smc[i][j].p[smc_IP3]
											- smc[right][j].p[smc_IP3]));
		}	//end j
	}	//end i

	for (i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (j = 1; j <= grid.num_ec_axially; j++) {
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
			ec[i][j].B[cpl_Ca] = -cpl_cef.Ca_hm_ec
					* ((ec[i][j].q[ec_Ca] - ec[i][up].q[ec_Ca])
							+ (ec[i][j].q[ec_Ca] - ec[i][down].q[ec_Ca])
							+ (ec[i][j].q[ec_Ca] - ec[left][j].q[ec_Ca])
							+ (ec[i][j].q[ec_Ca] - ec[right][j].q[ec_Ca]));
			ec[i][j].B[cpl_Vm] =
			-cpl_cef.Vm_hm_ec
			* ((ec[i][j].q[ec_Vm] - ec[i][up].q[ec_Vm])
					+ (ec[i][j].q[ec_Vm] - ec[i][down].q[ec_Vm])
					+ (ec[i][j].q[ec_Vm] - ec[left][j].q[ec_Vm])
					+ (ec[i][j].q[ec_Vm] - ec[right][j].q[ec_Vm]));
			ec[i][j].B[cpl_IP3] = -cpl_cef.IP3_hm_ec
					* ((ec[i][j].q[ec_IP3] - ec[i][up].q[ec_IP3])
							+ (ec[i][j].q[ec_IP3] - ec[i][down].q[ec_IP3])
							+ (ec[i][j].q[ec_IP3] - ec[left][j].q[ec_IP3])
							+ (ec[i][j].q[ec_IP3] - ec[right][j].q[ec_IP3]));

		}	//end j
	}	//end i

////******************** HETROCELLULAR COUPLING *********************/
	int offset_smc_circumferentially, offset_ec_axially;

	i = 0;
	j = 0;
	k = 0;
	l = 0;

	for (i = 1; i <= grid.num_smc_circumferentially; i++) {
		l = 1;
		for (j = 1; j <= grid.num_smc_axially; j++) {
			double dummy_smc[3] = { 0.0, 0.0, 0.0 };
			for (k = 1 + (i - 1) * 5; k <= i * 5; k++) {
				dummy_smc[cpl_Ca] = dummy_smc[cpl_Ca]
						+ (smc[i][j].p[smc_Ca] - ec[k][l].q[ec_Ca]);
				dummy_smc[cpl_Vm] = dummy_smc[cpl_Vm]
						+ (smc[i][j].p[smc_Vm] - ec[k][l].q[ec_Vm]);
				dummy_smc[cpl_IP3] = dummy_smc[cpl_IP3]
						+ (smc[i][j].p[smc_IP3] - ec[k][l].q[ec_IP3]);
			}
			if ((j % grid.num_smc_fundblk_axially) == 0) {
				l++;
			}
			smc[i][j].C[cpl_Ca] = -cpl_cef.Ca_ht_smc * dummy_smc[cpl_Ca];
			smc[i][j].C[cpl_Vm] = -cpl_cef.Vm_ht_smc * dummy_smc[cpl_Vm];
			smc[i][j].C[cpl_IP3] = -cpl_cef.IP3_ht_smc * dummy_smc[cpl_IP3];
		}
	}
	i = 0;
	j = 0;
	k = 0;
	l = 0;
	for (i = 1; i <= grid.num_ec_circumferentially; i++) {
		if ((i - 1) % 5 == 0)
			k++;
		for (j = 1; j <= grid.num_ec_axially; j++) {
			double dummy_ec[3] = { 0.0, 0.0, 0.0 };

			for (l = 1 + (j - 1) * 13; l <= j * 13; l++) {
				dummy_ec[cpl_Ca] = dummy_ec[cpl_Ca]
						+ (ec[i][j].q[ec_Ca] - smc[k][l].p[smc_Ca]);
				dummy_ec[cpl_Vm] = dummy_ec[cpl_Vm]
						+ (ec[i][j].q[ec_Vm] - smc[k][l].p[smc_Vm]);
				dummy_ec[cpl_IP3] = dummy_ec[cpl_IP3]
						+ (ec[i][j].q[ec_IP3] - smc[k][l].p[smc_IP3]);
			}
			ec[i][j].C[cpl_Ca] = -cpl_cef.Ca_ht_ec * dummy_ec[cpl_Ca];
			ec[i][j].C[cpl_Vm] = -cpl_cef.Vm_ht_ec * dummy_ec[cpl_Vm];
			ec[i][j].C[cpl_IP3] = -cpl_cef.IP3_ht_ec * dummy_ec[cpl_IP3];

		}
	}

}	//end of coupling()

int allocate_memory(grid_parms grid, celltype1** smc, celltype2** ec, double** sendbuf,double** recvbuf, checkpoint_handle* check)
{
	smc 	= (celltype1**) checked_malloc((grid.num_smc_circumferentially+grid.num_ghost_cells)* sizeof(celltype1*), stdout, "smc");
	for (int i=0; i<(grid.num_smc_circumferentially+grid.num_ghost_cells); i++){
		smc[i]	= (celltype1*) checked_malloc((grid.num_smc_axially+grid.num_ghost_cells)* sizeof(celltype1), stdout, "smc column dimension");

	}
	ec 	= (celltype2**) checked_malloc((grid.num_ec_circumferentially+grid.num_ghost_cells)* sizeof(celltype2*), stdout, "ec");
	for (int i=0; i<(grid.num_ec_circumferentially+grid.num_ghost_cells); i++){
		ec[i]	= (celltype2*) checked_malloc((grid.num_ec_axially+grid.num_ghost_cells)* sizeof(celltype2), stdout, "ec column dimension");
	}



///Memory allocation for state vector, the single cell evaluation placeholders (The RHS of the ODEs for each cell) and coupling fluxes is implemented in this section.
///In ghost cells, only the state vector array for each type of cells exists including all other cells.
///The memory is allocated for all the cells except the ghost cells, hence the ranges 1 to grid.num_ec_circumferentially(inclusive).
	///SMC domain
	for (int i = 0; i < (grid.num_smc_circumferentially+grid.num_ghost_cells); i++) {
		for (int j = 0; j < (grid.num_smc_axially+grid.num_ghost_cells); j++) {

			smc[i][j].p	=	(double*)checked_malloc(grid.neq_smc*sizeof(double),stdout,"allocation of array for state variables failed");
			smc[i][j].A = (double*) checked_malloc(
					grid.num_fluxes_smc * sizeof(double), stdout,
					"matrix A in smc");
			smc[i][j].B = (double*) checked_malloc(
					grid.num_coupling_species_smc * sizeof(double), stdout,
					"matrix B in smc");
			smc[i][j].C = (double*) checked_malloc(
					grid.num_coupling_species_smc * sizeof(double), stdout,
					"matrix C in smc");
		}
	}

	///EC domain
	for (int i = 0; i < (grid.num_ec_circumferentially+grid.num_ghost_cells); i++) {
		for (int j = 0; j < (grid.num_ec_axially+grid.num_ghost_cells); j++) {
			ec[i][j].q	=	(double*)checked_malloc(grid.neq_ec*sizeof(double),stdout,"allocation of array for state variables failed");

			ec[i][j].A = (double*) checked_malloc(
					grid.num_fluxes_ec * sizeof(double), stdout,
					"matrix A in ec");
			ec[i][j].B = (double*) checked_malloc(
					grid.num_coupling_species_ec * sizeof(double), stdout,
					"matrix B in ec");
			ec[i][j].C = (double*) checked_malloc(
					grid.num_coupling_species_ec * sizeof(double), stdout,
					"matrix C in ec");
		}
	}




	///Allocating memory space for coupling data to be sent and received by MPI communication routines.

	///sendbuf and recvbuf are 2D arrays having up,down,left and right directions as their first dimension.
	///Each dimension is broken down into two segments, e.g. up1,up2,down1 & down2,etc..
	///The length of the second dimension is equal to half the number of cells for which the information is to be sent and received.
	///Thus each communicating pair will exchange data twice to get the full lenght.

	sendbuf = (double**) checked_malloc(8 * sizeof(double*), stdout,
			"sendbuf dimension 1");
	recvbuf = (double**) checked_malloc(8 * sizeof(double*), stdout,
			"recvbuf dimension 1");

	///Each processor now allocates the memory for send and recv buffers those will hold the coupling information.
	///Since sendbuf must contain the information of number of SMCs and ECs being sent in the directions,
	///the first two elements contain the total count of SMCs located on the rank in the relevant dimension (circumferential or axial) and the count of SMCs for which
	///information is being sent, respectively.
	///The next two elements contain the same information for ECs.

	int extent_s, extent_e;	///Variables to calculate the length of the prospective buffer based on number of cell in either orientations (circumferential or axial).

	grid.added_info_in_send_buf = 4;///Number of elements containing additional information at the beginning of the send buffer.
	int seg_config_s, seg_config_e;	///Integers to decided whether the row or column being sent is overlapping or exactly divisible into two halfs.
	/// data to send to the neighbour in UP1 direction
	extent_s = (int) (ceil((double) (grid.num_smc_circumferentially) / 2));
	extent_e = (int) (ceil((double) (grid.num_ec_circumferentially) / 2));

	/// The seg_config variables are to recording the configuration of the split of the buffering in each direction.
	/// If the total number of cell in on a face (UP, DOWN, LEFT or RIGHT) are EVEN (i.e. seg_congif=0), the split will be non-overlapping
	/// (eg. if total SMCs are 26 in UP direction, UP1 buffer will send 13 and UP2 will send the other 13 to corresponding nbrs.
	/// If the total number of cells is ODD (i.e. seg_congif=1), then the split will be over lapping.
	/// (eg. if total SMCs are 13 (or any multiple of 13) UP1 will send elements from 0 - 6 and UP2 will send 6 - 12  to corresponding nbrs.
	/// These variables are used in send and recv buffers update before and after the MPI-communication routine is called.
	seg_config_s = grid.num_smc_circumferentially % 2;
	seg_config_e = grid.num_ec_circumferentially % 2;

	///Recording the number of elements in Send buffer in Up direction (UP1 or UP2) for use in update routine as count of elements.
	grid.num_elements_send_up = grid.added_info_in_send_buf
			+ (grid.num_coupling_species_smc * extent_s
					+ grid.num_coupling_species_ec * extent_e);

	sendbuf[UP1] = (double*) checked_malloc(
			grid.num_elements_send_up * sizeof(double), stdout,
			"sendbuf[UP1] dimension 2");

	sendbuf[UP1][0] = (double) (1); //Start of the 1st segment of SMC array in  in UP direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[UP1][1] = (double) (extent_s); //End of the 1st segment of SMC array in UP direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[UP1][2] = (double) (1); //Start of the 1st segment of EC array in UP direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[UP1][3] = (double) (extent_e); ///End of the 1st segment of EC array in UP direction (circumferential direction) to be sent to neighbouring processor

	sendbuf[UP2] = (double*) checked_malloc(
			grid.num_elements_send_up * sizeof(double), stdout,
			"sendbuf[UP2] dimension 2");

	if (seg_config_s != 0) {
		sendbuf[UP2][0] = (double) (extent_s); //Start of the 2nd segment of SMC array in UP direction (circumferential direction) to be sent to neighbouring processor
	} else if (seg_config_s == 0) {
		sendbuf[UP2][0] = (double) (extent_s + 1); //Start of the 2nd segment of SMC array in UP direction (circumferential direction) to be sent to neighbouring processor
	}
	sendbuf[UP2][1] = (double) (grid.num_smc_circumferentially); //End of the 2nd segment of SMC array in  in UP direction (circumferential direction) to be sent to neighbouring processor

	if (seg_config_e != 0) {
		sendbuf[UP2][2] = (double) (extent_e); //Start of the 2nd segment of EC array in  in UP direction (circumferential direction) to be sent to neighbouring processor
	} else if (seg_config_e == 0) {
		sendbuf[UP2][2] = (double) (extent_e + 1); //Start of the 2nd segment of EC array in  in UP direction (circumferential direction) to be sent to neighbouring processor
	}
	sendbuf[UP2][3] = (double) (grid.num_ec_circumferentially); //End of the 2nd segment of EC array in  in UP direction (circumferential direction) to be sent to neighbouring processor


	/// data to send to the neighbour in DOWN direction
	///Recording the number of elements in Send buffer in DOWN direction (DOWN1 or DOWN2) for use in update routine as count of elements.
	grid.num_elements_send_down = grid.added_info_in_send_buf
			+ (grid.num_coupling_species_smc * extent_s
					+ grid.num_coupling_species_ec * extent_e);

	sendbuf[DOWN1] = (double*) checked_malloc(
			grid.num_elements_send_down * sizeof(double), stdout,
			"sendbuf[DOWN1] dimension 2");

	sendbuf[DOWN1][0] = (double) (1); //Start of the 1st segment of SMC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[DOWN1][1] = (double) (extent_s); //End of the 1st segment of SMC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[DOWN1][2] = (double) (1); //Start of the 1st segment of EC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[DOWN1][3] = (double) (extent_e); ///End of the 1st segment of EC array in DOWN direction (circumferential direction) to be sent to neighbouring processor

	sendbuf[DOWN2] = (double*) checked_malloc(
			grid.num_elements_send_down * sizeof(double), stdout,
			"sendbuf[DOWN2] dimension 2");
	if (seg_config_s != 0) {
		sendbuf[DOWN2][0] = (double) (extent_s); //Start of the 2nd segment of SMC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
	} else if (seg_config_s == 0) {
		sendbuf[DOWN2][0] = (double) (extent_s + 1); //Start of the 2nd segment of SMC array in DOWN direction (circumferential direction) to be sent to neighbouring processor
	}

	sendbuf[DOWN2][1] = (double) (grid.num_smc_circumferentially); //End of the 2nd segment of SMC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor
	if (seg_config_e != 0) {
		sendbuf[DOWN2][2] = (double) (extent_e); //Start of the 2nd segment of EC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor
	} else if (seg_config_e == 0) {
		sendbuf[DOWN2][2] = (double) (extent_e + 1); //Start of the 2nd segment of EC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor
	}
	sendbuf[DOWN2][3] = (double) (grid.num_ec_circumferentially); //End of the 2nd segment of EC array in  in DOWN direction (circumferential direction) to be sent to neighbouring processor

	/// data to send to the neighbour in LEFT direction
	extent_s = (int) (ceil((double) (grid.num_smc_axially) / 2));
	extent_e = (int) (ceil((double) (grid.num_ec_axially) / 2));
	seg_config_s = grid.num_smc_axially % 2;
	seg_config_e = grid.num_ec_axially % 2;

	///Recording the number of elements in Send buffer in Left direction (LEFT1 or LEFT2) for use in update routine as count of elements.
	grid.num_elements_send_left = grid.added_info_in_send_buf
			+ (grid.num_coupling_species_smc * extent_s
					+ grid.num_coupling_species_ec * extent_e);
	sendbuf[LEFT1] = (double*) checked_malloc(
			grid.num_elements_send_left * sizeof(double), stdout,
			"sendbuf[LEFT1] dimension 2");
	sendbuf[LEFT1][0] = (double) (1); //Start of the 1st segment of SMC array in  in LEFT direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[LEFT1][1] = (double) (extent_s); //END of the 1st segment of SMC array in  in LEFT direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[LEFT1][2] = (double) (1); //Start of the 1st segment of EC array in  in LEFT direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[LEFT1][3] = (double) (extent_e); //END of the 1st segment of EC array in  in LEFT direction (circumferential direction) to be sent to neighbouring processor

	sendbuf[LEFT2] = (double*) checked_malloc(
			grid.num_elements_send_left * sizeof(double), stdout,
			"sendbuf[LEFT2] dimension 2");
	if (seg_config_s != 0) {
		sendbuf[LEFT2][0] = (double) (extent_s); //Start of the 2nd segment of SMC array in LEFT direction (circumferential direction) to be sent to neighbouring processor
	} else if (seg_config_s == 0) {
		sendbuf[LEFT2][0] = (double) (extent_s + 1); //Start of the 2nd segment of SMC array in LEFT direction (circumferential direction) to be sent to neighbouring processor
	}
	sendbuf[LEFT2][1] = (double) (grid.num_smc_axially); //END of the 2nd segment of SMC array in LEFT direction (circumferential direction) to be sent to neighbouring processor
	if (seg_config_e != 0) {
		sendbuf[LEFT2][2] = (double) (extent_e); //Start of the 2nd segment of EC array in LEFT direction (circumferential direction) to be sent to neighbouring processor
	} else if (seg_config_e == 0) {
		sendbuf[LEFT2][2] = (double) (extent_e + 1); //Start of the 2nd segment of EC array in LEFT direction (circumferential direction) to be sent to neighbouring processor
	}
	sendbuf[LEFT2][3] = (double) (grid.num_ec_axially); //END of the 2nd segment of EC array in LEFT direction (circumferential direction) to be sent to neighbouring processor

	/// data to send to the neighbour in RIGHT direction

	///Recording the number of elements in Send buffer in RIGHT direction (RIGHT1 or RIGHT2) for use in update routine as count of elements.
	grid.num_elements_send_right = grid.added_info_in_send_buf
			+ (grid.num_coupling_species_smc * extent_s
					+ grid.num_coupling_species_ec * extent_e);

	sendbuf[RIGHT1] = (double*) checked_malloc(
			grid.num_elements_send_right * sizeof(double), stdout,
			"sendbuf[RIGHT1] dimension 2");
	sendbuf[RIGHT1][0] = (double) (1); //Start of the 1st segment of SMC array in  in RIGHT direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[RIGHT1][1] = (double) (extent_s); //END of the 1st segment of SMC array in  in RIGHT direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[RIGHT1][2] = (double) (1); //Start of the 1st segment of EC array in  in RIGHT direction (circumferential direction) to be sent to neighbouring processor
	sendbuf[RIGHT1][3] = (double) (extent_e); //END of the 1st segment of EC array in  in RIGHT direction (circumferential direction) to be sent to neighbouring processor

	sendbuf[RIGHT2] = (double*) checked_malloc(
			grid.num_elements_send_right * sizeof(double), stdout,
			"sendbuf[RIGHT2] dimension 2");
	if (seg_config_s != 0) {
		sendbuf[RIGHT2][0] = (double) (extent_s); //Start of the 2nd segment of SMC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor
	} else if (seg_config_s == 0) {
		sendbuf[RIGHT2][0] = (double) (extent_s + 1); //Start of the 2nd segment of SMC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor
	}
	sendbuf[RIGHT2][1] = (double) (grid.num_smc_axially); //END of the 2nd segment of SMC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor
	if (seg_config_e != 0) {
		sendbuf[RIGHT2][2] = (double) (extent_e); //Start of the 2nd segment of EC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor
	} else if (seg_config_e == 0) {
		sendbuf[RIGHT2][2] = (double) (extent_e + 1); //Start of the 2nd segment of EC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor
	}
	sendbuf[RIGHT2][3] = (double) (grid.num_ec_axially); //END of the 2nd segment of EC array in RIGHT direction (circumferential direction) to be sent to neighbouring processor


	///Call communication to the number of elements to be recieved by neighbours and allocate memory of recvbuf for each direction accordingly.
		grid = communicate_num_recv_elements_to_nbrs(check->logptr,grid);
		///memory allocation

		fprintf(check->logptr,"recvbuf (u,d,l,r): (%d,%d,%d,%d)\n",
				grid.num_elements_recv_up,grid.num_elements_recv_down,grid.num_elements_recv_left,grid.num_elements_recv_right);


	/// data to receive from the neighbour in UP direction
	recvbuf[UP1] = (double*) checked_malloc(
			grid.num_elements_recv_up * sizeof(double), stdout,
			"recvbuf[UP1] dimension 2");
	recvbuf[UP2] = (double*) checked_malloc(
			grid.num_elements_recv_up * sizeof(double), stdout,
			"recvbuf[UP2] dimension 2");

	/// data to recv from the neighbour in DOWN direction
	recvbuf[DOWN1] = (double*) checked_malloc(
			grid.num_elements_recv_down * sizeof(double), stdout,
			"recvbuf[DOWN1] dimension 2");
	recvbuf[DOWN2] = (double*) checked_malloc(
			grid.num_elements_recv_down * sizeof(double), stdout,
			"recvbuf[DOWN2] dimension 2");

	/// data to receive from the neighbour in LEFT direction
	recvbuf[LEFT1] = (double*) checked_malloc(
			grid.num_elements_recv_left * sizeof(double), stdout,
			"recvbuf[LEFT1] dimension 2");
	recvbuf[LEFT2] = (double*) checked_malloc(
			grid.num_elements_recv_left * sizeof(double), stdout,
			"recvbuf[LEFT2] dimension 2");

	/// data to receive from the neighbour in RIGHT direction
	recvbuf[RIGHT1] = (double*) checked_malloc(
			grid.num_elements_recv_right * sizeof(double), stdout,
			"recvbuf[RIGHT1] dimension 2");
	recvbuf[RIGHT2] = (double*) checked_malloc(
			grid.num_elements_recv_right * sizeof(double), stdout,
			"recvbuf[RIGHT2] dimension 2");
return (success);
}
