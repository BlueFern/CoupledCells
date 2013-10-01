/*
 * koenigsberger_model.cpp
 *
 *  Created on: 12/06/2013
 *      Author: mohsinshaikh
 */
#include "computelib.h"
#include "koenigsberger_constants.h"

/*******************************************************************************************/
void Initialize_koeingsberger_smc(grid_parms grid, double* y, celltype1** smc)
/*******************************************************************************************/
{
	int k = 0, offset;

	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;
			/*y[k + ((j - 1) * grid.neq_smc) + smc_Ca] = 0.01;
			 y[k + ((j - 1) * grid.neq_smc) + smc_SR] = 0.01;
			 y[k + ((j - 1) * grid.neq_smc) + smc_Vm] = 0.01;
			 y[k + ((j - 1) * grid.neq_smc) + smc_w] = 0.01;
			 y[k + ((j - 1) * grid.neq_smc) + smc_IP3] = 0.01;*/

			y[k + ((j - 1) * grid.neq_smc) + smc_Ca] = (double)(grid.universal_rank)+(double) (i) * 1e-1
					+ (double) (j) * 1e-4 + (smc_Ca+1) * 1e-7;
			y[k + ((j - 1) * grid.neq_smc) + smc_SR] = (double)(grid.universal_rank)+(double) (i) * 1e-1
					+ (double) (j) * 1e-4 + (smc_SR+1) * 1e-7;
			y[k + ((j - 1) * grid.neq_smc) + smc_Vm] = (double)(grid.universal_rank)+(double) (i) * 1e-1
					+ (double) (j) * 1e-4 + (smc_Vm+1) * 1e-7;
			y[k + ((j - 1) * grid.neq_smc) + smc_w] = (double)(grid.universal_rank)+(double) (i) * 1e-1
					+ (double) (j) * 1e-4 + (smc_w+1) * 1e-7;
			y[k + ((j - 1) * grid.neq_smc) + smc_IP3] = (double)(grid.universal_rank)+(double) (i) * 1e-1
					+ (double) (j) * 1e-4 + (smc_IP3+1) * 1e-7;

		}
	}

	for (int i = 0; i < (grid.num_smc_circumferentially + grid.num_ghost_cells);
			i++) {
		for (int j = 0; j < (grid.num_smc_axially + grid.num_ghost_cells);
				j++) {
			smc[i][j].p[smc_Ca] = 0.0;
			smc[i][j].p[smc_SR] = 0.0;
			smc[i][j].p[smc_Vm] = 0.0;
			smc[i][j].p[smc_w] = 0.0;
			smc[i][j].p[smc_IP3] = 0.0;

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
void Initialize_koeingsberger_ec(grid_parms grid, double* y, celltype2** ec)
/*******************************************************************************************/
{
	int k, offset = (grid.neq_smc * grid.num_smc_circumferentially
			* grid.num_smc_axially);

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			if (i > 1)
				k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
				k = offset + 0;
			/*y[k + ((j - 1) * grid.neq_ec) + ec_Ca] = 0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_SR] = 0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_Vm] = 0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_IP3] = 0.01;*/

			y[k + ((j - 1) * grid.neq_ec) + ec_Ca] = -0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_SR] = -0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_Vm] = -0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_IP3] = -0.01;
		}
	}
	for (int i = 0; i < (grid.num_ec_circumferentially + grid.num_ghost_cells);
			i++) {
		for (int j = 0; j < (grid.num_ec_axially + grid.num_ghost_cells); j++) {
			ec[i][j].q[ec_Ca] = 0.0;
			ec[i][j].q[ec_SR] = 0.0;
			ec[i][j].q[ec_Vm] = 0.0;
			ec[i][j].q[ec_IP3] = 0.0;

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
void koenigsberger_smc(grid_parms grid, celltype1** smc)
/*******************************************************************************************/
{

	int i, j, k;
//#pragma omp for private(j)
//EVALUATING SINGLE CELL FLUXES :::::::For SMC::::::::
	for (i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (j = 1; j <= grid.num_smc_axially; j++) {

			//JIP3
			smc[i][j].A[J_IP3] = (Fi * P2(smc[i][j].p[smc_IP3]))
					/ (P2(Kri) + P2(smc[i][j].p[smc_IP3]));
			//JSRuptake
			smc[i][j].A[J_SERCA] = (Bi * P2(smc[i][j].p[smc_Ca])) / (P2(
					smc[i][j].p[smc_Ca]) + P2(cbi));
			//Jcicr
			smc[i][j].A[J_CICR] = (CICRi * (P2(smc[i][j].p[smc_SR]) * P4(
					smc[i][j].p[smc_Ca])))
					/ ((P2(sci) + P2(smc[i][j].p[smc_SR]))
					* (P4(cci)+ P4(
									smc[i][j].p[smc_Ca])));
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
					/ (P2(
							smc[i][j].p[smc_Ca] + cwi)
							+ (beta
									* ((double) exp(
											(-1) * (smc[i][j].p[smc_Vm] - vCa3)
													/ RKi))));
//Jdegradation
			smc[i][j].A[J_IP3_deg] = ki * smc[i][j].p[smc_IP3];

		} //end for j
	} //end for i
}
/***************************************************************************/
void koenigsberger_smc_derivatives(double* f, grid_parms grid,
		celltype1** smc) {
	/***************************************************************************/

	int k;
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;
			f[k + ((j - 1) * grid.neq_smc) + smc_Ca] = smc[i][j].A[J_IP3]
					- smc[i][j].A[J_SERCA] + smc[i][j].A[J_CICR]
					- smc[i][j].A[J_Extrusion] + smc[i][j].A[J_Leak]
					- smc[i][j].A[J_VOCC] + smc[i][j].A[J_Na_Ca]
					+ smc[i][j].B[cpl_Ca] + smc[i][j].C[cpl_Ca];

			f[k + ((j - 1) * grid.neq_smc) + smc_SR] = smc[i][j].A[J_SERCA]
					- smc[i][j].A[J_CICR] - smc[i][j].A[J_Leak];

			f[k + ((j - 1) * grid.neq_smc) + smc_Vm] = gama
					* (-smc[i][j].A[J_Na_K] - smc[i][j].A[J_Cl]
							- (2 * smc[i][j].A[J_VOCC]) - smc[i][j].A[J_Na_Ca]
							- smc[i][j].A[J_K]) + smc[i][j].B[cpl_Vm]
					+ smc[i][j].C[cpl_Vm];

			f[k + ((j - 1) * grid.neq_smc) + smc_w] = lambda
					* (smc[i][j].A[K_activation] - smc[i][j].p[smc_w]);

			f[k + ((j - 1) * grid.neq_smc) + smc_IP3] = -smc[i][j].A[J_IP3_deg]
					+ smc[i][j].B[cpl_IP3] + smc[i][j].C[cpl_IP3];
		}
	}
}

/************************************************************/
void koenigsberger_ec(grid_parms grid, celltype2** ec)
/*************************************************************                                                                                                                                                                         *
 * This is the multicell version for evaluating single cell ionic
 * currents in an EC
 *************************************************************/
{
	//EVALUATING SINGLE CELL FLUXES For EC
	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
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
}

/***************************************************************************/
void koenigsberger_ec_derivatives(double t, double* f, grid_parms grid,
		celltype2** ec) {
	/***************************************************************************/
	int k, offset = (grid.neq_smc * grid.num_smc_circumferentially
			* grid.num_smc_axially);
	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			if (i > 1)
				k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
				k = offset + 0;

			ec[i][j].JPLC = agonist_profile(t, grid, i, j, ec[i][j].z_coord);

			f[k + ((j - 1) * grid.neq_ec) + ec_Ca] = ec[i][j].A[J_IP3]
					- ec[i][j].A[J_SERCA] + ec[i][j].A[J_CICR]
					- ec[i][j].A[J_Extrusion] + ec[i][j].A[J_Leak]
					+ ec[i][j].A[J_NSC] + ec[i][j].A[J_trivial_Ca]
					+ ec[i][j].B[cpl_Ca] + ec[i][j].C[cpl_Ca];

			f[k + ((j - 1) * grid.neq_ec) + ec_SR] = ec[i][j].A[J_SERCA]
					- ec[i][j].A[J_CICR] - ec[i][j].A[J_Leak];

			f[k + ((j - 1) * grid.neq_ec) + ec_Vm] = ((-1 / Cmj)
					* (ec[i][j].A[J_Ktot] + ec[i][j].A[J_Residual]))
					+ ec[i][j].B[cpl_Vm] + ec[i][j].C[cpl_Vm];

			f[k + ((j - 1) * grid.neq_ec) + ec_IP3] = ec[i][j].JPLC
					- ec[i][j].A[J_IP3_deg] + ec[i][j].B[cpl_IP3]
					+ ec[i][j].C[cpl_IP3];

		}
	}
}

