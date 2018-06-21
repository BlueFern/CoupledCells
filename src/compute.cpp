#include <malloc.h>
#include <stdlib.h> // malloc

#include "computelib.h"
#include "koenigsberger_model.h"

using namespace std;
time_stamps t_stamp;

// For electro-diffusion.
#define NUM_SPECIES 4


// Electro-diffusion coupling parameters. Most from Jacobsen et al. 2007.
double
	PCa = 5e-8, absT = 293.15, zCa = 2, r_cons = 8345, ec_seg = 4.16,
	sigma = 0.15, f_cons = 96487, PK = 5e-8, PNa = 4e-8, PCl = 3e-8,
	Cl_cyt = 59400, K_cyt = 140000, Na_cyt = 8400, zCl = -1, zNa = 1, zK = 1,
	A_smc = 200, A_ec = 500, A_cell = 350, V_smc = 1, smc_scaling = 25, ec_scaling = 35;


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
	// These values are honestly just trial and error guesses at what makes nice waves...
	// Really need to have a think about what they should be.
	// Logic is to prioritise waves moving in the direction of larger cell contact (smcs axially, ecs circumferentially)
	cpl_cef->ec_diffusion[0] = 0.66;
	cpl_cef->ec_diffusion[1] = 0.66;
	cpl_cef->ec_diffusion[2] = 1.33;
	cpl_cef->ec_diffusion[3] = 1.33;

	cpl_cef->smc_diffusion[0] = 1.33;
	cpl_cef->smc_diffusion[1] = 1.33;
	cpl_cef->smc_diffusion[2] = 0.66;
	cpl_cef->smc_diffusion[3] = 0.66;


	// Similar principal to above. From jacobsen et al. 2007. Stands for contact ratio between cells.
	cpl_cef->ec_rho[0] = 0.15;
	cpl_cef->ec_rho[1] = 0.15;
	cpl_cef->ec_rho[2] = 0.3;
	cpl_cef->ec_rho[3] = 0.3;

	cpl_cef->smc_rho[0] = 0.3;
	cpl_cef->smc_rho[1] = 0.3;
	cpl_cef->smc_rho[2] = 0.15;
	cpl_cef->smc_rho[3] = 0.15;

	// 1/20 (0.05) is estimated contact between EC and SMC based on our EC:SMC ratio per quad/unit area.
	cpl_cef->ec_smc_rho = 0.05;


	if(CASE == 1)
	{
		cpl_cef->Vm_hm_smc = 1000.00;
		cpl_cef->Vm_hm_ec = 1000.00;

		cpl_cef->Ca_hm_smc = 0.1;
		cpl_cef->Ca_hm_ec = 0.1;

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

		cpl_cef->Ca_hm_smc = 0.01;
		cpl_cef->Ca_hm_ec = 0.01;

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

// Mapping from state variable vector to cells.
int map_solver_output_to_cells(const grid_parms& grid, double* y, SMC_cell** smc, EC_cell** ec)
{
	int err = 0;

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
			smc[i][j].vars[smc_IP3] = y[k + ((j - 1) * grid.neq_smc) + smc_IP3];
		}
	}

	k, offset = (grid.neq_smc * grid.num_smc_circumferentially
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

	k, offset = (grid.neq_smc * grid.num_smc_circumferentially
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
			ec[i][j].vars[ec_Gprot] = y[k + ((j - 1) * grid.neq_ec) + ec_Gprot];
		}
	}


	return (err);
}


/*
if using um^2 and uM then we are returning femto amperes
*/
double single_electro_diffusion(double average, double gradient, double vm_gradient, double permeability, int valence, double cell_area)
{
		return permeability * cell_area * f_cons * (gradient + ((valence * f_cons * average * vm_gradient) / (r_cons * absT)));

}

double gap_junction_current(double cell_a_Ca, double cell_b_Ca, double cell_a_vm, double cell_b_vm, double cell_area)
{
	double permabilities[NUM_SPECIES] = {PCa, PK, PCl, PNa};
	int valences[NUM_SPECIES] = {zCa, zK, zCl, zNa};
	double gradients[NUM_SPECIES] = {cell_a_Ca - cell_b_Ca, 0, 0, 0};
	double averages[NUM_SPECIES] = {(cell_a_Ca + cell_b_Ca) / 2.0, K_cyt, Cl_cyt, Na_cyt};

	double vm_gradient = cell_a_vm - cell_b_vm;

	double current_coupling = 0;
	for (int i = 0; i < NUM_SPECIES; i++)
	{
		current_coupling += single_electro_diffusion(averages[i], gradients[i], vm_gradient, permabilities[i], valences[i], cell_area);
	}
	return current_coupling;

}


void coupling_implicit(double t, double y[], const grid_parms& grid, SMC_cell** smc, EC_cell** ec, const conductance& cpl_cef)
{
	int i, j, k, l;

	l = 0;
	for (i = 1; i <= grid.num_smc_circumferentially; i++) {
		l = 1;
		for (j = 1; j <= grid.num_smc_axially; j++) {
////******************** HOMOCELLULAR COUPLING SMCs *********************/
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;

			// 0.3 and 0.15 are the fractions of contact between cells. From Jacobsen et al. 2007.
			smc[i][j].homo_fluxes[cpl_Vm] =
					(cpl_cef.smc_rho[0] * gap_junction_current(smc[i][j].vars[smc_Ca], smc[i][up].vars[smc_Ca], smc[i][j].vars[smc_Vm], smc[i][up].vars[smc_Vm], A_smc)
					+ cpl_cef.smc_rho[1] * gap_junction_current(smc[i][j].vars[smc_Ca], smc[i][down].vars[smc_Ca], smc[i][j].vars[smc_Vm], smc[i][down].vars[smc_Vm], A_smc)
					+ cpl_cef.smc_rho[2] * gap_junction_current(smc[i][j].vars[smc_Ca], smc[left][j].vars[smc_Ca], smc[i][j].vars[smc_Vm], smc[left][j].vars[smc_Vm], A_smc)
					+ cpl_cef.smc_rho[3] * gap_junction_current(smc[i][j].vars[smc_Ca], smc[right][j].vars[smc_Ca], smc[i][j].vars[smc_Vm], smc[right][j].vars[smc_Vm], A_smc));

			// Divide by picofarads to get millivolts/s
			smc[i][j].homo_fluxes[cpl_Vm] /= -Cmj;

////******************** HETEROCELLULAR COUPLING SMCs *********************/
			double dummy_smc[3] = { 0.0, 0.0, 0.0 };
			for (k = 1 + (i - 1) * 5; k <= i * 5; k++)
			{
				dummy_smc[cpl_Vm] += gap_junction_current(smc[i][j].vars[smc_Ca], ec[k][l].vars[ec_Ca], smc[i][j].vars[smc_Vm], ec[k][l].vars[ec_Vm], A_cell);
			}
			if ((j % grid.num_smc_fundblk_axially) == 0) {
				l++;
			}

			// Divide by picofarads to get millivolts/s
			smc[i][j].hetero_fluxes[cpl_Vm] = cpl_cef.ec_smc_rho * dummy_smc[cpl_Vm] / -Cmj;
		}	//end j
	}	//end i

	k = 0;
	for (i = 1; i <= grid.num_ec_circumferentially; i++) {
		if ((i - 1) % 5 == 0)
			k++;
		for (j = 1; j <= grid.num_ec_axially; j++) {
////******************** HOMOCELLULAR COUPLING ECs *********************/
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;

			// 0.3 and 0.15 are the fractions of contact between cells. From Jacobsen et al. 2007.
			ec[i][j].homo_fluxes[cpl_Vm] =
					(cpl_cef.ec_rho[0] * gap_junction_current(ec[i][j].vars[ec_Ca], ec[i][up].vars[ec_Ca], ec[i][j].vars[ec_Vm], ec[i][up].vars[ec_Vm], A_ec)
					+ cpl_cef.ec_rho[1] * gap_junction_current(ec[i][j].vars[ec_Ca], ec[i][down].vars[ec_Ca], ec[i][j].vars[ec_Vm], ec[i][down].vars[ec_Vm], A_ec)
					+ cpl_cef.ec_rho[2] * gap_junction_current(ec[i][j].vars[ec_Ca], ec[left][j].vars[ec_Ca], ec[i][j].vars[ec_Vm], ec[left][j].vars[ec_Vm], A_ec)
					+ cpl_cef.ec_rho[3] * gap_junction_current(ec[i][j].vars[ec_Ca], ec[right][j].vars[ec_Ca], ec[i][j].vars[ec_Vm], ec[right][j].vars[ec_Vm], A_ec));

			// Divide by picofarads to get millivolts/s
			ec[i][j].homo_fluxes[cpl_Vm] /= -Cmj;

////******************** HETEROCELLULAR COUPLING ECs *********************/
			double dummy_ec[3] = { 0.0, 0.0, 0.0 };
			for (l = 1 + (j - 1) * 13; l <= j * 13; l++)
			{
				dummy_ec[cpl_Vm] += gap_junction_current(ec[i][j].vars[ec_Ca], smc[k][l].vars[smc_Ca], ec[i][j].vars[ec_Vm], smc[k][l].vars[smc_Vm], A_cell);
			}


			// Divide by picofarads to get millivolts/s
			ec[i][j].hetero_fluxes[cpl_Vm] = cpl_cef.ec_smc_rho * dummy_ec[cpl_Vm] / -Cmj;
		}	//end j
	}	//end i

}
void coupling_explicit(double t, double y[], const grid_parms& grid, SMC_cell** smc, EC_cell** ec, const conductance& cpl_cef)
{
	int i, j, k, l;
	double average;
	double grad;
	double vm_grad;
	for (i = 1; i <= grid.num_smc_circumferentially; i++) {
		l = 1;
		for (j = 1; j <= grid.num_smc_axially; j++) {
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
////******************** HOMOCELLULAR COUPLING SMCs *********************/


//			smc[i][j].homo_fluxes[cpl_Ca] = -cpl_cef.Ca_hm_smc
//											* (cpl_cef.smc_diffusion[0] * ((smc[i][j].vars[smc_Ca] - smc[i][up].vars[smc_Ca]))
//											   + (cpl_cef.smc_diffusion[1] * (smc[i][j].vars[smc_Ca] - smc[i][down].vars[smc_Ca]))
//											   + (cpl_cef.smc_diffusion[2] * (smc[i][j].vars[smc_Ca] - smc[left][j].vars[smc_Ca]))
//											   + (cpl_cef.smc_diffusion[3] * (smc[i][j].vars[smc_Ca] - smc[right][j].vars[smc_Ca])));

			smc[i][j].homo_fluxes[cpl_Ca] = 0.0;

			// coupling with above cell
			average = (smc[i][j].vars[smc_Ca] + smc[i][up].vars[smc_Ca]) / 2.0;
			grad = (smc[i][j].vars[smc_Ca] - smc[i][up].vars[smc_Ca]);
			vm_grad = (smc[i][j].vars[smc_Vm] - smc[i][up].vars[smc_Vm]);
			smc[i][j].homo_fluxes[cpl_Ca]  += cpl_cef.smc_diffusion[0] * single_electro_diffusion(average, grad, vm_grad, PCa, 2, A_smc);

			// coupling with cell below
			average = (smc[i][j].vars[smc_Ca] + smc[i][down].vars[smc_Ca]) / 2.0;
			grad = (smc[i][j].vars[smc_Ca] -  smc[i][down].vars[smc_Ca]);
			vm_grad = (smc[i][j].vars[smc_Vm] - smc[i][down].vars[smc_Vm]);
			smc[i][j].homo_fluxes[cpl_Ca]  += cpl_cef.smc_diffusion[1] * single_electro_diffusion(average, grad, vm_grad, PCa, 2, A_smc);

			// coupling to the left
			average = (smc[i][j].vars[smc_Ca] + smc[left][j].vars[smc_Ca]) / 2.0;
			grad = (smc[i][j].vars[smc_Ca] - smc[left][j].vars[smc_Ca]);
			vm_grad = (smc[i][j].vars[smc_Vm] - smc[left][j].vars[smc_Vm]);
			smc[i][j].homo_fluxes[cpl_Ca]  += cpl_cef.smc_diffusion[2] * single_electro_diffusion(average, grad, vm_grad, PCa, 2, A_smc);

			// coupling to the right
			average = (smc[i][j].vars[smc_Ca] + smc[right][j].vars[smc_Ca]) / 2.0;
			grad = (smc[i][j].vars[smc_Ca] - smc[right][j].vars[smc_Ca]);
			vm_grad = (smc[i][j].vars[smc_Vm] - smc[right][j].vars[smc_Vm]);
			smc[i][j].homo_fluxes[cpl_Ca]  += cpl_cef.smc_diffusion[3] * single_electro_diffusion(average, grad, vm_grad, PCa, 2, A_smc);


			// now convert fA to mM/s, then x 1000 to uM/s
			smc[i][j].homo_fluxes[cpl_Ca] = -1000 * (smc[i][j].homo_fluxes[cpl_Ca] / (f_cons * 2 * V_smc)) * smc_scaling;

			// Not electro diffused give IP3 has no charge.
			smc[i][j].homo_fluxes[cpl_IP3] = -cpl_cef.IP3_hm_smc
					* (cpl_cef.smc_diffusion[0] * ((smc[i][j].vars[smc_IP3] - smc[i][up].vars[smc_IP3]))
					+ (cpl_cef.smc_diffusion[1] * (smc[i][j].vars[smc_IP3] - smc[i][down].vars[smc_IP3]))
					+ (cpl_cef.smc_diffusion[2] * (smc[i][j].vars[smc_IP3] - smc[left][j].vars[smc_IP3]))
					+ (cpl_cef.smc_diffusion[3] * (smc[i][j].vars[smc_IP3] - smc[right][j].vars[smc_IP3])));

////******************** HETROCELLULAR COUPLING SMCs *********************/
			double dummy_smc[3] = { 0.0, 0.0, 0.0 };
			for (k = 1 + (i - 1) * 5; k <= i * 5; k++)
			{
				dummy_smc[cpl_Ca] += smc[i][j].vars[smc_Ca] - ec[k][l].vars[ec_Ca];
				dummy_smc[cpl_IP3] += smc[i][j].vars[smc_IP3] - ec[k][l].vars[ec_IP3];
			}
			if ((j % grid.num_smc_fundblk_axially) == 0) {
				l++;
			}
			smc[i][j].hetero_fluxes[cpl_Ca] = -cpl_cef.Ca_ht_smc * dummy_smc[cpl_Ca];
			smc[i][j].hetero_fluxes[cpl_IP3] = -cpl_cef.IP3_ht_smc * dummy_smc[cpl_IP3];
		}	//end j
	}	//end i

	k = 0;
	for (i = 1; i <= grid.num_ec_circumferentially; i++) {
		if ((i - 1) % 5 == 0)
			k++;
		for (j = 1; j <= grid.num_ec_axially; j++) {
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
////******************** HOMOCELLULAR COUPLING ECs *********************/


//			ec[i][j].homo_fluxes[cpl_Ca] = -cpl_cef.Ca_hm_ec
//										   * cpl_cef.ec_diffusion[0] * ((ec[i][j].vars[ec_Ca] - ec[i][up].vars[ec_Ca])
//																		+ cpl_cef.ec_diffusion[1] * (ec[i][j].vars[ec_Ca] - ec[i][down].vars[ec_Ca])
//																		+ cpl_cef.ec_diffusion[2] * (ec[i][j].vars[ec_Ca] - ec[left][j].vars[ec_Ca])
//																		+ cpl_cef.ec_diffusion[3] * (ec[i][j].vars[ec_Ca] - ec[right][j].vars[ec_Ca]));

			ec[i][j].homo_fluxes[cpl_Ca] = 0.0;

			// coupling with above cell
			average = (ec[i][j].vars[ec_Ca] + ec[i][up].vars[ec_Ca]) / 2.0;
			grad = (ec[i][j].vars[ec_Ca] - ec[i][up].vars[ec_Ca]);
			vm_grad = (ec[i][j].vars[ec_Vm] - ec[i][up].vars[ec_Vm]);
			ec[i][j].homo_fluxes[cpl_Ca]  += cpl_cef.ec_diffusion[0] * single_electro_diffusion(average, grad, vm_grad, PCa, 2, A_ec);

			// coupling with cell below
			average = (ec[i][j].vars[ec_Ca] + ec[i][down].vars[ec_Ca]) / 2.0;
			grad = (ec[i][j].vars[ec_Ca] -  ec[i][down].vars[ec_Ca]);
			vm_grad = (ec[i][j].vars[ec_Vm] - ec[i][down].vars[ec_Vm]);
			ec[i][j].homo_fluxes[cpl_Ca]  += cpl_cef.ec_diffusion[1] * single_electro_diffusion(average, grad, vm_grad, PCa, 2, A_ec);

			// coupling to the left
			average = (ec[i][j].vars[ec_Ca] + ec[left][j].vars[ec_Ca]) / 2.0;
			grad = (ec[i][j].vars[ec_Ca] - ec[left][j].vars[ec_Ca]);
			vm_grad = (ec[i][j].vars[ec_Vm] - ec[left][j].vars[ec_Vm]);
			ec[i][j].homo_fluxes[cpl_Ca]  += cpl_cef.ec_diffusion[2] * single_electro_diffusion(average, grad, vm_grad, PCa, 2, A_ec);

			// coupling to the right
			average = (ec[i][j].vars[ec_Ca] + ec[right][j].vars[ec_Ca]) / 2.0;
			grad = (ec[i][j].vars[ec_Ca] - ec[right][j].vars[ec_Ca]);
			vm_grad = (ec[i][j].vars[ec_Vm] - ec[right][j].vars[ec_Vm]);
			ec[i][j].homo_fluxes[cpl_Ca]  += cpl_cef.ec_diffusion[3] * single_electro_diffusion(average, grad, vm_grad, PCa, 2, A_ec);

			ec[i][j].homo_fluxes[cpl_Ca] = -1000 * (ec[i][j].homo_fluxes[cpl_Ca] / (f_cons * 2 * V_smc)) * ec_scaling;

			// Not electro diffused give IP3 has no charge.
			ec[i][j].homo_fluxes[cpl_IP3] = -cpl_cef.IP3_hm_ec
					* cpl_cef.ec_diffusion[0] * ((ec[i][j].vars[ec_IP3] - ec[i][up].vars[ec_IP3])
					+ cpl_cef.ec_diffusion[1] * (ec[i][j].vars[ec_IP3] - ec[i][down].vars[ec_IP3])
					+ cpl_cef.ec_diffusion[2] * (ec[i][j].vars[ec_IP3] - ec[left][j].vars[ec_IP3])
					+ cpl_cef.ec_diffusion[3] * (ec[i][j].vars[ec_IP3] - ec[right][j].vars[ec_IP3]));
////******************** HETEROCELLULAR COUPLING ECs *********************/
			double dummy_ec[3] = { 0.0, 0.0, 0.0 };
			for (l = 1 + (j - 1) * 13; l <= j * 13; l++)
			{
				dummy_ec[cpl_Ca] += ec[i][j].vars[ec_Ca] - smc[k][l].vars[smc_Ca];
				dummy_ec[cpl_IP3] += ec[i][j].vars[ec_IP3] - smc[k][l].vars[smc_IP3];
			}
			ec[i][j].hetero_fluxes[cpl_Ca] = -cpl_cef.Ca_ht_ec * dummy_ec[cpl_Ca];
			ec[i][j].hetero_fluxes[cpl_IP3] = -cpl_cef.IP3_ht_ec * dummy_ec[cpl_IP3];

		}	//end j
	}	//end i

}

void coupling(double t, double y[], const grid_parms& grid, SMC_cell** smc, EC_cell** ec, const conductance& cpl_cef)
{
	coupling_implicit(t, y, grid, smc, ec, cpl_cef);
	coupling_explicit(t, y, grid, smc, ec, cpl_cef);
}


void compute(const grid_parms& grid, SMC_cell** smc, EC_cell** ec, const conductance& cpl_cef, double t, double* y, double* f)
{
	int err;

	map_solver_output_to_cells(grid, y, smc, ec);

	koenigsberger_smc(grid, smc);
	koenigsberger_ec(grid, ec);

	coupling(t, y, grid, smc, ec, cpl_cef);

	koenigsberger_smc_derivatives(f, grid, smc);
	koenigsberger_ec_derivatives(t, f, grid, ec);

}

void compute_implicit(const grid_parms& grid, SMC_cell** smc, EC_cell** ec, const conductance& cpl_cef, double t, double* y, double* f)
{
	map_solver_output_to_cells(grid, y, smc, ec);

	koenigsberger_smc_implicit(grid, smc);
	koenigsberger_ec_implicit(grid, ec);

	coupling_implicit(t, y, grid, smc, ec, cpl_cef);

	koenigsberger_smc_derivatives_implicit(f, grid, smc);
	koenigsberger_ec_derivatives_implicit(t, f, grid, ec);
}

void compute_explicit(const grid_parms& grid, SMC_cell** smc, EC_cell** ec, const conductance& cpl_cef, double t, double* y, double* f)
{
	map_solver_output_to_cells(grid, y, smc, ec);

	koenigsberger_smc_explicit(grid, smc);
	koenigsberger_ec_explicit(grid, ec);

	coupling_explicit(t, y, grid, smc, ec, cpl_cef);

	koenigsberger_smc_derivatives_explicit(f, grid, smc);
	koenigsberger_ec_derivatives_explicit(t, f, grid, ec);
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
