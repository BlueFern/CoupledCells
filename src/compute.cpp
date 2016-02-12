#include <malloc.h>

#include "computelib.h"
#include "koenigsberger_model.h"
#include "tsoukias_model.h"

using namespace std;
time_stamps t_stamp;

/**
 * Wrapper around malloc to catch failed memory allocation. If allocation fails MPI_Abort is called.
 *
 * \param bytes Size of requested memory.
 * \param errmsg Message produced in the event of failed memory allocation.
 */
void* checked_malloc(size_t bytes, const char* errmsg)
{
	void *pval = malloc(bytes);
	if (pval == NULL)
	{
		fprintf(stderr, "MEMORY ALLOCATION ERROR: %s\n", errmsg);
		MPI_Abort(MPI_COMM_WORLD, 911);
	}
	return pval;
}

void set_coupling_parms(int CASE, conductance* cpl_cef)
{
	cpl_cef->ec_diffusion[0] = 0.5;
	cpl_cef->ec_diffusion[1] = 0.5;
	cpl_cef->ec_diffusion[2] = 0.5;
	cpl_cef->ec_diffusion[3] = 0.5;
	cpl_cef->ec_diffusion[4] = 1.0;
	cpl_cef->ec_diffusion[5] = 1.0;

	cpl_cef->smc_diffusion[0] = 1.0;
	cpl_cef->smc_diffusion[1] = 1.0;
	cpl_cef->smc_diffusion[2] = 0.5;
	cpl_cef->smc_diffusion[3] = 0.5;
	cpl_cef->smc_diffusion[4] = 0.5;
	cpl_cef->smc_diffusion[5] = 0.5;

	if(CASE == 1)
	{
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
	}
	else if(CASE == 2)
	{
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
	}
	else if(CASE == 3)
	{
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
	}
	else if(CASE == 4)
	{
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
	}
	else if(CASE == 5)
	{
		// Simulating for experiments suggested by Dr. James Kozloski (IBM Watson Centre).
		// The homocellular Ca coupling between SMCs is changed to investigate the effects of
		// strength on coupling on the propagation speed of the spatial waves.
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
	}
	else if(CASE == 6)
	{
		// Simulating for experiments suggested by Dr. James Kozloski (IBM Watson Centre).
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
	}
	else if(CASE == 7)
	{
		// Simulating for experiments suggested by Dr. James Kazloski (IBM Watson Centre).
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
	}
	else if(CASE == 8)
	{
		// What is this one?
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
}

// TODO: Move the Tsoukias code to the appropriate location.
// Mapping from state variable vector to cells.
int map_solver_output_to_cells(grid_parms grid, double* y, SMC_cell** smc, EC_cell** ec)
{
	int err = 0;
	switch (grid.smc_model)
	{
	case (TSK): {
		int k = 0, offset;
		for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
			for (int j = 1; j <= grid.num_smc_axially; j++) {
				if (i > 1)
					k = ((i - 1) * grid.neq_smc_axially);
				else if (i == 1)
					k = 0;
				smc[i][j].vars[smc_Vm] = y[k + ((j - 1) * grid.neq_smc) + smc_Vm];
				smc[i][j].vars[smc_d_L] = y[k + ((j - 1) * grid.neq_smc) + smc_d_L];
				smc[i][j].vars[smc_f_L] = y[k + ((j - 1) * grid.neq_smc) + smc_f_L];
				smc[i][j].vars[smc_p_f] = y[k + ((j - 1) * grid.neq_smc) + smc_p_f];
				smc[i][j].vars[smc_p_s] = y[k + ((j - 1) * grid.neq_smc) + smc_p_s];
				smc[i][j].vars[smc_q_1] = y[k + ((j - 1) * grid.neq_smc) + smc_q_1];
				smc[i][j].vars[smc_q_2] = y[k + ((j - 1) * grid.neq_smc) + smc_q_2];
				smc[i][j].vars[smc_p_K] = y[k + ((j - 1) * grid.neq_smc) + smc_p_K];
				smc[i][j].vars[smc_Ca_u] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca_u];
				smc[i][j].vars[smc_Ca_r] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca_r];
				smc[i][j].vars[smc_R_10] = y[k + ((j - 1) * grid.neq_smc) + smc_R_10];
				smc[i][j].vars[smc_R_11] = y[k + ((j - 1) * grid.neq_smc) + smc_R_11];
				smc[i][j].vars[smc_R_01] = y[k + ((j - 1) * grid.neq_smc) + smc_R_01];
				smc[i][j].vars[smc_h_IP3] = y[k + ((j - 1) * grid.neq_smc) + smc_h_IP3];
				smc[i][j].vars[smc_R_S_G] = y[k + ((j - 1) * grid.neq_smc) + smc_R_S_G];
				smc[i][j].vars[smc_R_S_P_G] = y[k + ((j - 1) * grid.neq_smc) + smc_R_S_P_G];
				smc[i][j].vars[smc_G] = y[k + ((j - 1) * grid.neq_smc) + smc_G];
				smc[i][j].vars[smc_IP3] = y[k + ((j - 1) * grid.neq_smc) + smc_IP3];
				smc[i][j].vars[smc_PIP2] = y[k + ((j - 1) * grid.neq_smc) + smc_PIP2];
				smc[i][j].vars[smc_V_cGMP] = y[k + ((j - 1) * grid.neq_smc) + smc_V_cGMP];
				smc[i][j].vars[smc_cGMP_i] = y[k + ((j - 1) * grid.neq_smc) + smc_cGMP_i];
				smc[i][j].vars[smc_Ca] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca];
				smc[i][j].vars[smc_Na_i] = y[k + ((j - 1) * grid.neq_smc) + smc_Na_i];
				smc[i][j].vars[smc_K_i] = y[k + ((j - 1) * grid.neq_smc) + smc_K_i];
				smc[i][j].vars[smc_Cl_i] = y[k + ((j - 1) * grid.neq_smc) + smc_Cl_i];
				smc[i][j].vars[smc_DAG] = y[k + ((j - 1) * grid.neq_smc) + smc_DAG];
			}
		}
		break;
	}
	case (KNBGR): {
		int k = 0, offset;
		for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
			for (int j = 1; j <= grid.num_smc_axially; j++) {
				if (i > 1)
					k = ((i - 1) * grid.neq_smc_axially);
				else if (i == 1)
					k = 0;
				smc[i][j].vars[smc_Ca] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca];
				smc[i][j].vars[smc_SR] = y[k + ((j - 1) * grid.neq_smc) + smc_SR];
				smc[i][j].vars[smc_Vm] = y[k + ((j - 1) * grid.neq_smc) + smc_Vm];
				smc[i][j].vars[smc_w] = y[k + ((j - 1) * grid.neq_smc) + smc_w];
				smc[i][j].vars[smc_IP3] =y[k + ((j - 1) * grid.neq_smc) + smc_IP3];
			}
		}
		break;
	}
	default: {
		err = 1;
		break;
	}
	}
	switch (grid.ec_model) {
	case (TSK): {
		int k, offset = (grid.neq_smc * grid.num_smc_circumferentially
				* grid.num_smc_axially);

		for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
			for (int j = 1; j <= grid.num_ec_axially; j++) {
				if (i > 1)
					k = offset + ((i - 1) * grid.neq_ec_axially);
				else if (i == 1)
					k = offset + 0;
				ec[i][j].vars[ec_Ca] = y[k + ((j - 1) * grid.neq_ec) + ec_Ca];
				ec[i][j].vars[ec_SR] = y[k + ((j - 1) * grid.neq_ec) + ec_SR];
				ec[i][j].vars[ec_Vm] = y[k + ((j - 1) * grid.neq_ec) + ec_Vm];
				ec[i][j].vars[ec_IP3] = y[k + ((j - 1) * grid.neq_ec) + ec_IP3];
			}
		}
		break;
	}
	case (KNBGR): {
		int k, offset = (grid.neq_smc * grid.num_smc_circumferentially
				* grid.num_smc_axially);

		for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
			for (int j = 1; j <= grid.num_ec_axially; j++) {
				if (i > 1)
					k = offset + ((i - 1) * grid.neq_ec_axially);
				else if (i == 1)
					k = offset + 0;
				ec[i][j].vars[ec_Ca] = y[k + ((j - 1) * grid.neq_ec) + ec_Ca];
				ec[i][j].vars[ec_SR] = y[k + ((j - 1) * grid.neq_ec) + ec_SR];
				ec[i][j].vars[ec_Vm] = y[k + ((j - 1) * grid.neq_ec) + ec_Vm];
				ec[i][j].vars[ec_IP3] = y[k + ((j - 1) * grid.neq_ec) + ec_IP3];
			}
		}
		break;
	}
	default: {
		err = 1;
		break;
	}
	}

	return (err);
}

// Assigns the corner cells the values of the adjacent (same row) cell.
void copy_smc_corner_coupling_species(grid_parms* grid, SMC_cell** smc, int species)
{
	smc[0][0].vars[species] = smc[1][1].vars[species];
	smc[grid->num_smc_circumferentially + 1][0].vars[species] = smc[grid->num_smc_circumferentially][0].vars[species];
	smc[0][grid->num_smc_axially + 1].vars[species] = smc[1][grid->num_smc_axially + 1].vars[species];
	smc[grid->num_smc_circumferentially + 1][grid->num_smc_axially + 1].vars[species] = smc[grid->num_smc_circumferentially][grid->num_smc_axially + 1].vars[species];
}

// Assigns the corner cells the values of the adjacent (same row) cell.
void copy_ec_corner_coupling_species(grid_parms* grid, EC_cell** ec, int species)
{
	ec[0][0].vars[species] = ec[1][1].vars[species];
	ec[grid->num_ec_circumferentially + 1][0].vars[species] = ec[grid->num_ec_circumferentially][0].vars[species];
	ec[0][grid->num_ec_axially + 1].vars[species] = ec[1][grid->num_ec_axially + 1].vars[species];
	ec[grid->num_ec_circumferentially + 1][grid->num_ec_axially + 1].vars[species] = ec[grid->num_ec_circumferentially][grid->num_ec_axially + 1].vars[species];
}


void coupling_implicit(double t, double y[], grid_parms grid, SMC_cell** smc, EC_cell** ec, conductance cpl_cef)
{
	int i, j, k, l;


	copy_smc_corner_coupling_species(&grid, smc, smc_Vm);
	copy_ec_corner_coupling_species(&grid, ec, ec_Vm);


////******************** HOMOCELLULAR COUPLING *********************/
	int offset = 0;
	for (i = 1; i <= grid.num_smc_circumferentially; i++) {

		for (j = 1; j <= grid.num_smc_axially; j++) {
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
			smc[i][j].homo_fluxes[cpl_Vm] = -cpl_cef.Vm_hm_smc * (
					(cpl_cef.smc_diffusion[0] * (smc[i][j].vars[smc_Vm] - smc[i][up].vars[smc_Vm])) + // Top.
					(cpl_cef.smc_diffusion[1] * (smc[i][j].vars[smc_Vm] - smc[i][down].vars[smc_Vm])) + // Bottom.

					(cpl_cef.smc_diffusion[2] * (smc[i][j].vars[smc_Vm] - smc[left][j - offset].vars[smc_Vm])) + // Top left.
					(cpl_cef.smc_diffusion[3] * (smc[i][j].vars[smc_Vm] - smc[left][j - offset + 1].vars[smc_Vm])) + // Bottom Left.
					(cpl_cef.smc_diffusion[4] * (smc[i][j].vars[smc_Vm] - smc[right][j - offset].vars[smc_Vm])) + // Top Right.
					(cpl_cef.smc_diffusion[5] * (smc[i][j].vars[smc_Vm] - smc[right][j - offset + 1].vars[smc_Vm]))); // Bottom right.
		offset = 1 - offset;


		}	//end j
	}	//end i

	for (i = 1; i <= grid.num_ec_circumferentially; i++) {
		// Every second row needs a small offset to see the correct cells above and below.
		// On daughter branches the EC brick pattern is inverted (for a correct tessellation)
		if (grid.branch_tag != P)
		{
			offset = 0;
		}
		else
		{
			offset = 1;
		}

		for (j = 1; j <= grid.num_ec_axially; j++) {
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
			ec[i][j].homo_fluxes[cpl_Vm] = -cpl_cef.Vm_hm_ec * (
					(cpl_cef.ec_diffusion[0] * (ec[i][j].vars[ec_Vm] - ec[i - offset][up].vars[ec_Vm])) + // Top left.
					(cpl_cef.ec_diffusion[1] * (ec[i][j].vars[ec_Vm] - ec[i - offset + 1][up].vars[ec_Vm])) + // Top right.
					(cpl_cef.ec_diffusion[2] * (ec[i][j].vars[ec_Vm] - ec[i - offset][down].vars[ec_Vm])) + // Bottom left.
					(cpl_cef.ec_diffusion[3] * (ec[i][j].vars[ec_Vm] - ec[i - offset + 1][down].vars[ec_Vm])) + // Bottom right.

					(cpl_cef.ec_diffusion[4] * (ec[i][j].vars[ec_Vm] - ec[left][j].vars[ec_Vm])) + // Left.
					(cpl_cef.ec_diffusion[5] * (ec[i][j].vars[ec_Vm] - ec[right][j].vars[ec_Vm]))); // Right.

			offset = 1 - offset; // Flip between 0 and 1.

		}	//end j
	}	//end i

////******************** HETROCELLULAR COUPLING *********************/
	int offset_smc_circumferentially, offset_ec_axially;

	i = 0; // x dim.
	j = 0; // y dim.
	k = 0;
	l = 0;

	for (i = 1; i <= grid.num_smc_circumferentially; i++)
	{
		l = 1;
		for (j = 1; j <= grid.num_smc_axially; j++)
		{
			double dummy_smc[3] = { 0.0, 0.0, 0.0 };
			for (k = 1 + (i - 1) * 5; k <= i * 5; k++)
			{
				dummy_smc[cpl_Vm] = dummy_smc[cpl_Vm] + (smc[i][j].vars[smc_Vm] - ec[k][l].vars[ec_Vm]);
			}
			if ((j % grid.num_smc_fundblk_axially) == 0) {
				l++;
			}
			smc[i][j].hetero_fluxes[cpl_Vm] = -cpl_cef.Vm_ht_smc * dummy_smc[cpl_Vm];
		}
	}

	i = 0;
	j = 0;
	k = 0;
	l = 0;

	for (i = 1; i <= grid.num_ec_circumferentially; i++)
	{
		if ((i - 1) % 5 == 0)
			k++;
		for (j = 1; j <= grid.num_ec_axially; j++)
		{
			double dummy_ec[3] = { 0.0, 0.0, 0.0 };
			for (l = 1 + (j - 1) * 13; l <= j * 13; l++)
			{
				dummy_ec[cpl_Vm] = dummy_ec[cpl_Vm] + (ec[i][j].vars[ec_Vm] - smc[k][l].vars[smc_Vm]);
			}
			ec[i][j].hetero_fluxes[cpl_Vm] = -cpl_cef.Vm_ht_ec * dummy_ec[cpl_Vm];
		}
	}
}

void coupling_explicit(double t, double y[], grid_parms grid, SMC_cell** smc, EC_cell** ec, conductance cpl_cef)
{
	int i, j, k, l;

	copy_smc_corner_coupling_species(&grid, smc, smc_Ca);
	copy_smc_corner_coupling_species(&grid, smc, smc_IP3);
	copy_ec_corner_coupling_species(&grid, ec, ec_Ca);
	copy_ec_corner_coupling_species(&grid, ec, ec_IP3);


////******************** HOMOCELLULAR COUPLING *********************/
	int offset = 0;
	for (i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (j = 1; j <= grid.num_smc_axially; j++) {
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
			smc[i][j].homo_fluxes[cpl_Ca] = -cpl_cef.Ca_hm_smc * (
					(cpl_cef.smc_diffusion[0] * (smc[i][j].vars[smc_Ca] - smc[i][up].vars[smc_Ca])) + // Top.
					(cpl_cef.smc_diffusion[1] * (smc[i][j].vars[smc_Ca] - smc[i][down].vars[smc_Ca])) + // Bottom.
					(cpl_cef.smc_diffusion[2] * (smc[i][j].vars[smc_Ca] - smc[left][j - offset].vars[smc_Ca])) + // Top left.
					(cpl_cef.smc_diffusion[3] * (smc[i][j].vars[smc_Ca] - smc[left][j - offset + 1].vars[smc_Ca])) + // Bottom left.
					(cpl_cef.smc_diffusion[4] * (smc[i][j].vars[smc_Ca] - smc[right][j - offset].vars[smc_Ca])) + // Top right.
					(cpl_cef.smc_diffusion[5] * (smc[i][j].vars[smc_Ca] - smc[right][j - offset + 1].vars[smc_Ca]))); // Bottom right.
			// Identical logic/process for IP3 instead of Ca.
			smc[i][j].homo_fluxes[cpl_IP3] = -cpl_cef.IP3_hm_smc * (
					(cpl_cef.smc_diffusion[0] * (smc[i][j].vars[smc_IP3] - smc[i][up].vars[smc_IP3])) +
					(cpl_cef.smc_diffusion[1] * (smc[i][j].vars[smc_IP3] - smc[i][down].vars[smc_IP3])) +
					(cpl_cef.smc_diffusion[2] * (smc[i][j].vars[smc_IP3] - smc[left][j - offset].vars[smc_IP3])) +
					(cpl_cef.smc_diffusion[3] * (smc[i][j].vars[smc_IP3] - smc[left][j - offset + 1].vars[smc_IP3])) +
					(cpl_cef.smc_diffusion[4] * (smc[i][j].vars[smc_IP3] - smc[right][j - offset].vars[smc_IP3])) +
					(cpl_cef.smc_diffusion[5] * (smc[i][j].vars[smc_IP3] - smc[right][j - offset + 1].vars[smc_IP3])));
		offset = 1 - offset; // Flip between 0 and 1.
		}	//end j
	}	//end i

	for (i = 1; i <= grid.num_ec_circumferentially; i++) {
		// Every second row needs a small offset to see the correct cells above and below.
		// On daughter branches the EC brick pattern is inverted (for a correct tessellation)
		if (grid.branch_tag != P)
		{
			offset = 0;
		}
		else
		{
			offset = 1;
		}

		for (j = 1; j <= grid.num_ec_axially; j++) {
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
			ec[i][j].homo_fluxes[cpl_Ca] = -cpl_cef.Ca_hm_ec * (
					(cpl_cef.ec_diffusion[0] * (ec[i][j].vars[ec_Ca] - ec[i - offset][up].vars[ec_Ca])) + // Top left.
					(cpl_cef.ec_diffusion[1] * (ec[i][j].vars[ec_Ca] - ec[i - offset + 1][up].vars[ec_Ca])) + // Top right.
					(cpl_cef.ec_diffusion[2] * (ec[i][j].vars[ec_Ca] - ec[i - offset][down].vars[ec_Ca])) + // Bottom left.
					(cpl_cef.ec_diffusion[3] * (ec[i][j].vars[ec_Ca] - ec[i - offset + 1][down].vars[ec_Ca])) + // Bottom right.
					(cpl_cef.ec_diffusion[4] * (ec[i][j].vars[ec_Ca] - ec[left][j].vars[ec_Ca])) +
					(cpl_cef.ec_diffusion[5] * (ec[i][j].vars[ec_Ca] - ec[right][j].vars[ec_Ca])));

			// Identical logic/process for IP3 instead of Ca.
			ec[i][j].homo_fluxes[cpl_IP3] = -cpl_cef.IP3_hm_ec * (
					(cpl_cef.ec_diffusion[0] * (ec[i][j].vars[ec_IP3] - ec[i - offset][up].vars[ec_IP3])) +
					(cpl_cef.ec_diffusion[1] * (ec[i][j].vars[ec_IP3] - ec[i - offset + 1][up].vars[ec_IP3])) +
					(cpl_cef.ec_diffusion[2] * (ec[i][j].vars[ec_IP3] - ec[i - offset][down].vars[ec_IP3])) +
					(cpl_cef.ec_diffusion[3] * (ec[i][j].vars[ec_IP3] - ec[i - offset + 1][down].vars[ec_IP3])) +
					(cpl_cef.ec_diffusion[4] * (ec[i][j].vars[ec_IP3] - ec[left][j].vars[ec_IP3])) +
					(cpl_cef.ec_diffusion[5] * (ec[i][j].vars[ec_IP3] - ec[right][j].vars[ec_IP3])));
			offset = 1 - offset; // Flip between 0 and 1.

		}	//end j
	}	//end i

////******************** HETROCELLULAR COUPLING *********************/
	int offset_smc_circumferentially, offset_ec_axially;

	i = 0; // x dim.
	j = 0; // y dim.
	k = 0;
	l = 0;

	for (i = 1; i <= grid.num_smc_circumferentially; i++)
	{
		l = 1;
		for (j = 1; j <= grid.num_smc_axially; j++)
		{
			double dummy_smc[3] = { 0.0, 0.0, 0.0 };
			for (k = 1 + (i - 1) * 5; k <= i * 5; k++)
			{
				dummy_smc[cpl_Ca] = dummy_smc[cpl_Ca] + (smc[i][j].vars[smc_Ca] - ec[k][l].vars[ec_Ca]);
				dummy_smc[cpl_IP3] = dummy_smc[cpl_IP3] + (smc[i][j].vars[smc_IP3] - ec[k][l].vars[ec_IP3]);
			}
			if ((j % grid.num_smc_fundblk_axially) == 0) {
				l++;
			}
			smc[i][j].hetero_fluxes[cpl_Ca] = -cpl_cef.Ca_ht_smc * dummy_smc[cpl_Ca];
			smc[i][j].hetero_fluxes[cpl_IP3] = -cpl_cef.IP3_ht_smc * dummy_smc[cpl_IP3];
		}
	}

	i = 0;
	j = 0;
	k = 0;
	l = 0;

	for (i = 1; i <= grid.num_ec_circumferentially; i++)
	{
		if ((i - 1) % 5 == 0)
			k++;
		for (j = 1; j <= grid.num_ec_axially; j++)
		{
			double dummy_ec[3] = { 0.0, 0.0, 0.0 };
			for (l = 1 + (j - 1) * 13; l <= j * 13; l++)
			{
				dummy_ec[cpl_Ca] = dummy_ec[cpl_Ca] + (ec[i][j].vars[ec_Ca] - smc[k][l].vars[smc_Ca]);
				dummy_ec[cpl_IP3] = dummy_ec[cpl_IP3] + (ec[i][j].vars[ec_IP3] - smc[k][l].vars[smc_IP3]);
			}
			ec[i][j].hetero_fluxes[cpl_Ca] = -cpl_cef.Ca_ht_ec * dummy_ec[cpl_Ca];
			ec[i][j].hetero_fluxes[cpl_IP3] = -cpl_cef.IP3_ht_ec * dummy_ec[cpl_IP3];
		}
	}
}

void coupling(double t, double y[], grid_parms grid, SMC_cell** smc, EC_cell** ec, conductance cpl_cef)
{
	coupling_implicit(t, y, grid, smc, ec, cpl_cef);
	coupling_explicit(t, y, grid, smc, ec, cpl_cef);
}


void compute(grid_parms grid, SMC_cell** smc, EC_cell** ec, conductance cpl_cef, double t, double* y, double* f)
{
	int err;

	map_solver_output_to_cells(grid, y, smc, ec);

#if CELL_MODEL == TSK

	tsoukias_smc(grid, smc);
	koenigsberger_ec(grid, ec);

	coupling(t, y, grid, smc, ec, cpl_cef);

	tsoukias_smc_derivatives(f, grid, smc);
	koenigsberger_ec_derivatives(t, f, grid, ec);

#elif CELL_MODEL == KNBGR

	koenigsberger_smc(grid, smc);
	koenigsberger_ec(grid, ec);

	coupling(t, y, grid, smc, ec, cpl_cef);

	koenigsberger_smc_derivatives(f, grid, smc);
	koenigsberger_ec_derivatives(t, f, grid, ec);

#endif
}

void compute_implicit(grid_parms grid, SMC_cell** smc, EC_cell** ec, conductance cpl_cef, double t, double* y, double* f)
{
	map_solver_output_to_cells(grid, y, smc, ec);

#if CELL_MODEL == TSK

	tsoukias_smc(grid, smc);
	koenigsberger_ec(grid, ec);

	coupling(t, y, grid, smc, ec, cpl_cef);

	tsoukias_smc_derivatives(f, grid, smc);
	koenigsberger_ec_derivatives(t, f, grid, ec);

#elif CELL_MODEL == KNBGR

	koenigsberger_smc_implicit(grid, smc);
	koenigsberger_ec_implicit(grid, ec);

	coupling_implicit(t, y, grid, smc, ec, cpl_cef);

	koenigsberger_smc_derivatives_implicit(f, grid, smc);
	koenigsberger_ec_derivatives_implicit(t, f, grid, ec);

#endif
}

void compute_explicit(grid_parms grid, SMC_cell** smc, EC_cell** ec, conductance cpl_cef, double t, double* y, double* f)
{
	map_solver_output_to_cells(grid, y, smc, ec);

#if CELL_MODEL == TSK

	tsoukias_smc(grid, smc);
	koenigsberger_ec(grid, ec);

	coupling(t, y, grid, smc, ec, cpl_cef);

	tsoukias_smc_derivatives(f, grid, smc);
	koenigsberger_ec_derivatives(t, f, grid, ec);

#elif CELL_MODEL == KNBGRs

	koenigsberger_smc_explicit(grid, smc);
	koenigsberger_ec_explicit(grid, ec);

	coupling_explicit(t, y, grid, smc, ec, cpl_cef);

	koenigsberger_smc_derivatives_explicit(f, grid, smc);
	koenigsberger_ec_derivatives_explicit(t, f, grid, ec);

#endif
}

#if 0
// TODO: These functions are to be moved to a util_func.cpp or something of that sort.
/************************************************************/
void minimum(double* table, int size, double *value, int *index) {
	///For evaluating minimum of an array.
	*value = table[0];
	for (int i = 0; i < size; i++) {
		if (*value > table[i]) {
			*value = table[i];
			*index = i;
		}
	}
}

/************************************************************/
void maximum(double* table, int size, double *value, int *index) {
	///For evaluating maximum of an array.
	*value = table[0];
	for (int i = 0; i < size; i++) {
		if (*value < table[i]) {
			*value = table[i];
			*index = i;
		}
	}
}

/************************************************************/
void average(double* table, int size, double *value) {
	///For evaluating average of an array.
	*value = 0;
	for (int i = 0; i < size; i++) {
		*value += table[i];
	}
	*value = *value / (double) (size);
}
#endif


