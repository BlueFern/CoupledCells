#include <malloc.h>
#include <stdlib.h> // malloc

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
	// These values are honestly just trial and error guesses at what makes nice waves...
	// Really need to have a think about what they should be.
	cpl_cef->ec_diffusion[0] = 1;
	cpl_cef->ec_diffusion[1] = 1;
	cpl_cef->ec_diffusion[2] = 1;
	cpl_cef->ec_diffusion[3] = 1;

	cpl_cef->smc_diffusion[0] = 1;
	cpl_cef->smc_diffusion[1] = 1;
	cpl_cef->smc_diffusion[2] = 1;
	cpl_cef->smc_diffusion[3] = 1;


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

// TODO: Move the Tsoukias code to the appropriate location.
// Mapping from state variable vector to cells.
int map_solver_output_to_cells(const grid_parms& grid, 
                               double* y, SMC_cell** smc, EC_cell** ec)
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
				ec[i][j].vars[ec_Gprot] = y[k + ((j - 1) * grid.neq_ec) + ec_Gprot];
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

void coupling_implicit(double t, double y[], 
                       const grid_parms& grid, SMC_cell** smc, EC_cell** ec, 
                       const conductance& cpl_cef)
{
////******************** HOMOCELLULAR COUPLING *********************/
#pragma ivdep
#pragma omp parallel for
        for (int ij = 0; ij < grid.num_smc_circumferentially * grid.num_smc_axially; ij++) {
	        int i = ij / grid.num_smc_axially + 1;
                int j = ij % grid.num_smc_axially + 1;
		const double vSmc_Vm = smc[i][j].vars[smc_Vm];
		int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
                const double* __restrict__ vars = smc[i][j].vars;
                const double* __restrict__ upVars = smc[i][up].vars;
                const double* __restrict__ downVars = smc[i][down].vars;
                const double* __restrict__ leftVars = smc[left][j].vars;
                const double* __restrict__ rightVars = smc[right][j].vars;
                const double var = vars[smc_Vm];
                const double upVar = upVars[smc_Vm];
                const double downVar = downVars[smc_Vm];
                const double leftVar =  leftVars[smc_Vm];
                const double rightVar = rightVars[smc_Vm];
                const double* __restrict__ diffusion = cpl_cef.smc_diffusion;
                double* __restrict__ homo_fluxes = smc[i][j].homo_fluxes;

		homo_fluxes[cpl_Vm] = -cpl_cef.Vm_hm_smc
			  * (diffusion[0] * (var - upVar)
			  +  diffusion[1] * (var - downVar)
			  +  diffusion[2] * (var - leftVar)
			  +  diffusion[3] * (var - rightVar));

	}	//end ij

#pragma ivdep
#pragma omp parallel for
        for (int ij = 0; ij < grid.num_ec_circumferentially * grid.num_ec_axially; ij++) {
                int i = ij / grid.num_ec_axially + 1;
                int j = ij % grid.num_ec_axially + 1;
		int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
                const double* __restrict__ vars =  ec[i][j].vars;
                const double* __restrict__ upVars = ec[i][up].vars;
                const double* __restrict__ downVars = ec[i][down].vars;
                const double* __restrict__ leftVars = ec[left][j].vars;
                const double* __restrict__ rightVars = ec[right][j].vars;
                const double var = vars[ec_Vm];
                const double upVar = upVars[ec_Vm];
                const double downVar = downVars[ec_Vm];
                const double leftVar = leftVars[ec_Vm];
                const double rightVar = rightVars[ec_Vm];
                const double* __restrict__ diffusion = cpl_cef.ec_diffusion;
                double* __restrict__ homo_fluxes = ec[i][j].homo_fluxes;

		ec[i][j].homo_fluxes[cpl_Vm] = -cpl_cef.Vm_hm_ec
			  * (diffusion[0] * (var - upVar)
			  +  diffusion[1] * (var - downVar)
			  +  diffusion[2] * (var - leftVar)
			  +  diffusion[3] * (var - rightVar));

	}	//end ij

////******************** HETEROCELLULAR COUPLING *********************/
	int offset_smc_circumferentially, offset_ec_axially;

#pragma omp parallel for
        for (int ij = 0; ij < grid.num_smc_circumferentially * grid.num_smc_axially; ij++) {
	        int i = ij / grid.num_smc_axially + 1; // x dim
                int j = ij % grid.num_smc_axially + 1; // y dim
		int l = (j - 1) / grid.num_smc_fundblk_axially + 1;
                const double smcIJ = smc[i][j].vars[smc_Vm];
		double dummy_smc[3] = { 0.0, 0.0, 0.0 };
#pragma ivdep
		for (int k = 1 + (i - 1) * 5; k <= i * 5; k++) {
		      dummy_smc[cpl_Vm] = dummy_smc[cpl_Vm] + (smcIJ - ec[k][l].vars[ec_Vm]);
		}
		smc[i][j].hetero_fluxes[cpl_Vm] = -cpl_cef.Vm_ht_smc * dummy_smc[cpl_Vm];
	}

#pragma omp parallel for
        for (int ij = 0; ij < grid.num_ec_circumferentially * grid.num_ec_axially; ++ij) {
	        int i = ij / grid.num_ec_axially + 1;
		int j = ij % grid.num_ec_axially + 1;
                int k = (i - 1) / 5 + 1;
                const double ecIJ = ec[i][j].vars[ec_Vm];
		double dummy_ec[3] = { 0.0, 0.0, 0.0 };
#pragma ivdep
		for (int l = 1 + (j - 1) * 13; l <= j * 13; l++) {
		     dummy_ec[cpl_Vm] = dummy_ec[cpl_Vm] + (ecIJ - smc[k][l].vars[smc_Vm]);
		}
		ec[i][j].hetero_fluxes[cpl_Vm] = -cpl_cef.Vm_ht_ec * dummy_ec[cpl_Vm];
	}
}

void coupling_explicit(double t, double y[], 
                       const grid_parms& grid, 
                       SMC_cell** __restrict__ smc, 
                       EC_cell** __restrict__ ec, 
                       const conductance& cpl_cef)
{
////******************** HOMOCELLULAR COUPLING *********************/
#pragma ivdep
#pragma omp parallel for
	for (int ij = 0; ij < grid.num_smc_circumferentially * grid.num_smc_axially; ij++) {
	        int i = ij / grid.num_smc_axially + 1;
                int j = ij % grid.num_smc_axially + 1;
		int up = j - 1, down = j + 1, left = i - 1, right = i + 1;

                const double* __restrict__ vars = smc[i][j].vars;
                const double* __restrict__ upVars = smc[i][up].vars;
                const double* __restrict__ downVars = smc[i][down].vars;
                const double* __restrict__ leftVars = smc[left][j].vars;
                const double* __restrict__ rightVars = smc[right][j].vars;
                double* __restrict__ homo_fluxes = smc[i][j].homo_fluxes;

		homo_fluxes[cpl_Ca] = -cpl_cef.Ca_hm_smc
		                        * (cpl_cef.smc_diffusion[0] * ((vars[smc_Ca] - upVars[smc_Ca]))
					+ (cpl_cef.smc_diffusion[1] * (vars[smc_Ca] - downVars[smc_Ca]))
					+ (cpl_cef.smc_diffusion[2] * (vars[smc_Ca] - leftVars[smc_Ca]))
					+ (cpl_cef.smc_diffusion[3] * (vars[smc_Ca] - rightVars[smc_Ca])));

		homo_fluxes[cpl_IP3] = -cpl_cef.IP3_hm_smc
					* (cpl_cef.smc_diffusion[0] * ((vars[smc_IP3] - upVars[smc_IP3]))
					+ (cpl_cef.smc_diffusion[1] * (vars[smc_IP3] - downVars[smc_IP3]))
					+ (cpl_cef.smc_diffusion[2] * (vars[smc_IP3] - leftVars[smc_IP3]))
					+ (cpl_cef.smc_diffusion[3] * (vars[smc_IP3] - rightVars[smc_IP3])));
	}	//end ij

#pragma omp parallel for
	for (int ij = 0; ij < grid.num_ec_circumferentially * grid.num_ec_axially; ij++) {
                int i = ij / grid.num_ec_axially + 1;
                int j = ij % grid.num_ec_axially + 1;
                int up = j - 1, down = j + 1, left = i - 1, right = i + 1;

		const double* __restrict__ vars = ec[i][j].vars;
		const double* __restrict__ upVars = ec[i][up].vars;
                const double* __restrict__ downVars = ec[i][down].vars;
                const double* __restrict__ leftVars = ec[left][j].vars;
                const double* __restrict__ rightVars = ec[right][j].vars;
                double* __restrict__ homo_fluxes = ec[i][j].homo_fluxes;

		homo_fluxes[cpl_Ca] = -cpl_cef.Ca_hm_ec
					* cpl_cef.ec_diffusion[0] * ((vars[ec_Ca] - upVars[ec_Ca])
					+ cpl_cef.ec_diffusion[1] * (vars[ec_Ca] - downVars[ec_Ca])
					+ cpl_cef.ec_diffusion[2] * (vars[ec_Ca] - leftVars[ec_Ca])
					+ cpl_cef.ec_diffusion[3] * (vars[ec_Ca] - rightVars[ec_Ca]));

		homo_fluxes[cpl_IP3] = -cpl_cef.IP3_hm_ec
					* cpl_cef.ec_diffusion[0] * ((vars[ec_IP3] - upVars[ec_IP3])
					+ cpl_cef.ec_diffusion[1] * (vars[ec_IP3] - downVars[ec_IP3])
					+ cpl_cef.ec_diffusion[2] * (vars[ec_IP3] - leftVars[ec_IP3])
					+ cpl_cef.ec_diffusion[3] * (vars[ec_IP3] - rightVars[ec_IP3]));

	}	//end ij

////******************** HETEROCELLULAR COUPLING *********************/
	int offset_smc_circumferentially, offset_ec_axially;

#pragma omp parallel for
        for (int ij = 0; ij < grid.num_smc_circumferentially * grid.num_smc_axially; ij++) {
	        int i = ij / grid.num_smc_axially + 1; // x dim
	        int j = ij % grid.num_smc_axially + 1; // y dim
                int l = (j - 1) / grid.num_smc_fundblk_axially + 1;
	        double dummy_smc[3] = { 0.0, 0.0, 0.0 };
                const double* __restrict__ smcIJVars = smc[i][j].vars;
                double* __restrict__ hetero_fluxes = smc[i][j].hetero_fluxes;
		for (int k = 1 + (i - 1) * 5; k <= i * 5; k++) {
                        const double* __restrict__ ecKLVars = ec[k][l].vars;
		        dummy_smc[cpl_Ca] = dummy_smc[cpl_Ca] + (smcIJVars[smc_Ca] - ecKLVars[ec_Ca]);
			dummy_smc[cpl_IP3] = dummy_smc[cpl_IP3] + (smcIJVars[smc_IP3] - ecKLVars[ec_IP3]);
		}
		hetero_fluxes[cpl_Ca] = -cpl_cef.Ca_ht_smc * dummy_smc[cpl_Ca];
		hetero_fluxes[cpl_IP3] = -cpl_cef.IP3_ht_smc * dummy_smc[cpl_IP3];
	}

#pragma omp parallel for
        for (int ij = 0; ij < grid.num_ec_circumferentially * grid.num_ec_axially; ij++) {
	        int i = ij / grid.num_ec_axially + 1;
                int j = ij % grid.num_ec_axially + 1;
		int k = (i - 1) / 5 + 1;
                const double* __restrict__ ecIJVars = ec[i][j].vars;
                double* __restrict__ hetero_fluxes = ec[i][j].hetero_fluxes;
		double dummy_ec[3] = { 0.0, 0.0, 0.0 };
		for (int l = 1 + (j - 1) * 13; l <= j * 13; l++) {
		        const double* __restrict__ smcKLVars = smc[k][l].vars;
		        dummy_ec[cpl_Ca] = dummy_ec[cpl_Ca] + (ecIJVars[ec_Ca] - smcKLVars[smc_Ca]);
			dummy_ec[cpl_IP3] = dummy_ec[cpl_IP3] + (ecIJVars[ec_IP3] - smcKLVars[smc_IP3]);
		}
		hetero_fluxes[cpl_Ca] = -cpl_cef.Ca_ht_ec * dummy_ec[cpl_Ca];
	        hetero_fluxes[cpl_IP3] = -cpl_cef.IP3_ht_ec * dummy_ec[cpl_IP3];
	}
}

void coupling(double t, double y[], const grid_parms& grid, 
              SMC_cell** smc, EC_cell** ec, 
              const conductance& cpl_cef)
{
	coupling_implicit(t, y, grid, smc, ec, cpl_cef);
	coupling_explicit(t, y, grid, smc, ec, cpl_cef);
}


void compute(const grid_parms& grid, 
             SMC_cell** smc, 
             EC_cell** ec, 
             const conductance& cpl_cef, 
             double t, double* y, double* f)
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

#if PLOTTING && EXPLICIT_ONLY
	bufferPos = 0;
#endif

	koenigsberger_smc_derivatives(f, grid, smc);
	koenigsberger_ec_derivatives(t, f, grid, ec);

#if PLOTTING && EXPLICIT_ONLY

	if (grid.universal_rank == RANK)
	{
		for (int i = 0; i < OUTPUT_PLOTTING_SIZE; i++)
		{
			fprintf(var_file, "%f,", plotttingBuffer[i]);
		}
		fprintf(var_file, "%f\n", t);
	}
#endif


#endif
}

void compute_implicit(const grid_parms& grid, 
                      SMC_cell** smc, EC_cell** ec, 
                      const conductance& cpl_cef, 
                      double t, double* y, double* f)
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

void compute_explicit(const grid_parms& grid, 
                      SMC_cell** smc, EC_cell** ec, 
                      const conductance& cpl_cef, 
                      double t, double* y, double* f)
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
