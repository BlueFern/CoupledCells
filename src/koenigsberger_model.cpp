/*
 * koenigsberger_model.cpp
 *
 *  Created on: 12/06/2013
 *      Author: mohsinshaikh
 */

#include <math.h>
#include <cstdlib>
#include "computelib.h"
#include "koenigsberger_model.h"

double
	/* Constants for homogenically coupled SMCs. */
      Fi = 0.23,  Kri = 1.00,  GCai = 0.00129,  vCa1 = 100.00,
	  vCa2 = -24.00,  RCai = 8.50,  GNaCai = 0.00316,  cNaCai = 0.5,
	  vNaCai = -30.00,  Bi = 2.025,  cbi = 1.0,  CICRi = 55.00,
	  sci = 2.0,  cci = 0.9,  Di = 0.24,  vdi = -100.00,
	  Rdi = 250.00,  Li = 0.025, gama = 1970.00,
	  FNaK = 0.0432,  GCli = 0.00134,  vCl = -25.0,
	  GKi = 0.0046,  vKi = -94.00,  lambda = 45.00,
	  cwi = 0.0,  beta = 0.13,  vCa3 = -27.0,
	  RKi = 12.0,  ki = 0.1,
	/* Constants for homogenically coupled ECs. */
	  Fj = 0.23,  Krj = 1.00,  Bj = 0.5,
	  cbj = 1.0,  CICRj = 5.0,  scj = 2.0,
	  ccj = 0.9,  Dj = 0.24,  Lj = 0.025,
	  kj = 0.1,  Gcatj = 0.66 * 1e-3,  ECa = 50.00,
	  m3cat = -0.18,
	  m4cat = 0.37,  J0j = 0.029,  Cmj = 25.8,
	  Gtot = 6927,   vKj = -80.0,  a1j = 53.3,
	  a2j = 53.3,     bj = -80.8,  c1j = -0.4,
	  m3b = 1.32e-3, m4b = 0.30,   m3s = -0.28,
	  m4s = 0.389,   GRj = 955, vrestj = -31.10,
	/*Intracellular calcium buffering*/
	  // k6 = 100.00,  k7 = 300.00,  BT = 120.00;

	// Lemon at al. Constants
	kATP = 2, alpha_j = 2.781e-5, KCa = 0.4, R_PIP2_r = 10,
	PIP2_tot = 4e7, K_G_prot_act = 0.017, delta = 1.234e-3,
	G_prot_tot = 1e5, K_G_prot_deact = 0.15, N_a = 6.02252e23,
	V_ec = 1.17e3, unitcon_a = 1e21, PIP2_C = 5e7, kDeg = 1.25,
	cons_PIP2 = 3.989e7;


/// Initial values found by running a sufficiently long simulation and recording state values
/// after they have reached a steady state.
void initialize_koenigsberger_smc(const grid_parms& grid, double* y, SMC_cell** smc)
{
	int k = 0, offset;
	srand(grid.universal_rank);

	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;

			if (grid.random)
			{
				y[k + ((j - 1) * grid.neq_smc) + smc_Ca] = (float)rand() / (float)(RAND_MAX / 0.1) + 0.15;
				y[k + ((j - 1) * grid.neq_smc) + smc_SR] = (float)rand() / (float)(RAND_MAX / 0.2) + 1.3;
				y[k + ((j - 1) * grid.neq_smc) + smc_Vm] = (float)rand() / (float)(RAND_MAX / 5.0) - 68.0;
				y[k + ((j - 1) * grid.neq_smc) + smc_w] = (float)rand() / (float)(RAND_MAX / 0.005) + 0.01;
				y[k + ((j - 1) * grid.neq_smc) + smc_IP3] = (float)rand() / (float)(RAND_MAX / 0.2) + 0.6;
			}
			else
			{
				y[k + ((j - 1) * grid.neq_smc) + smc_Ca] = 0.15;
				y[k + ((j - 1) * grid.neq_smc) + smc_SR] = 1.3;
				y[k + ((j - 1) * grid.neq_smc) + smc_Vm] = - 68.0;
				y[k + ((j - 1) * grid.neq_smc) + smc_w] = 0.01;
				y[k + ((j - 1) * grid.neq_smc) + smc_IP3] = 0.6;
			}

		}
	}

	for (int i = 0; i < (grid.num_smc_circumferentially + grid.num_ghost_cells); i++) {
		for (int j = 0; j < (grid.num_smc_axially + grid.num_ghost_cells); j++) {
			smc[i][j].vars[smc_Ca] = 0.0;
			smc[i][j].vars[smc_SR] = 0.0;
			smc[i][j].vars[smc_Vm] = 0.0;
			smc[i][j].vars[smc_w] = 0.0;
			smc[i][j].vars[smc_IP3] = 0.0;

			// TODO: remove initialising fluxes for lemon model even when not used....
			for (int k = 1; k <= grid.num_fluxes_smc; k++) {
				smc[i][j].fluxes[k - 1] = 0.1;
			}
			for (int k = 1; k <= grid.num_coupling_species_smc; k++) {
				smc[i][j].homo_fluxes[k - 1] = 0.1;
				smc[i][j].hetero_fluxes[k - 1] = 0.1;
			}
		}
	}
}

/// Initial values found by running a sufficiently long simulation and recording state values
/// after they have reached a steady state.
void initialize_koenigsberger_ec(const grid_parms& grid, double* y, EC_cell** ec)
{
	int k, offset = (grid.neq_smc * grid.num_smc_circumferentially * grid.num_smc_axially);
	srand(grid.universal_rank);

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			if (i > 1)
				k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
				k = offset + 0;

			if (grid.random)
			{
				y[k + ((j - 1) * grid.neq_ec) + ec_Ca] = (float)rand() / (float)(RAND_MAX / 0.3) + 0.7;
				y[k + ((j - 1) * grid.neq_ec) + ec_SR] = (float)rand() / (float)(RAND_MAX / 0.2) + 0.5;
				y[k + ((j - 1) * grid.neq_ec) + ec_Vm] = (float)rand() / (float)(RAND_MAX / 5.0) - 68.0;
				y[k + ((j - 1) * grid.neq_ec) + ec_IP3] = (float)rand() / (float)(RAND_MAX / 0.1) + 0.9;
				y[k + ((j - 1) * grid.neq_ec) + ec_Gprot] = (float)rand() / (float)(RAND_MAX / 50.0) + 1450;
			}
			else
			{
				y[k + ((j - 1) * grid.neq_ec) + ec_Ca] = 0.825;
				y[k + ((j - 1) * grid.neq_ec) + ec_SR] = 0.63;
				y[k + ((j - 1) * grid.neq_ec) + ec_Vm] = -66.7;
				y[k + ((j - 1) * grid.neq_ec) + ec_IP3] = 1.057;
				y[k + ((j - 1) * grid.neq_ec) + ec_Gprot] = 1470.305;
			}

		}
	}
	for (int i = 0; i < (grid.num_ec_circumferentially + grid.num_ghost_cells); i++) {
		for (int j = 0; j < (grid.num_ec_axially + grid.num_ghost_cells); j++) {
			ec[i][j].vars[ec_Ca] = 0.0;
			ec[i][j].vars[ec_SR] = 0.0;
			ec[i][j].vars[ec_Vm] = 0.0;
			ec[i][j].vars[ec_IP3] = 0.0;
			ec[i][j].vars[ec_Gprot] = 0.0;

			for (int k = 1; k <= grid.num_fluxes_ec; k++) {
				ec[i][j].fluxes[k - 1] = 0.1;
			}
			for (int k = 1; k <= grid.num_coupling_species_ec; k++) {
				ec[i][j].homo_fluxes[k - 1] = 0.1;
				ec[i][j].hetero_fluxes[k - 1] = 0.1;
			}
		}
	}
}

void koenigsberger_smc(const grid_parms& grid, SMC_cell** smc)
{
	koenigsberger_smc_implicit(grid, smc);
	koenigsberger_smc_explicit(grid, smc);
}

void koenigsberger_smc_derivatives(double* f, const grid_parms& grid, SMC_cell** smc)
{
	koenigsberger_smc_derivatives_implicit(f, grid, smc);
	koenigsberger_smc_derivatives_explicit(f, grid, smc);
}

void koenigsberger_ec(const grid_parms& grid, EC_cell** ec, int atp_timestep)
{
	koenigsberger_ec_implicit(grid, ec);
	koenigsberger_ec_explicit(grid, ec, atp_timestep);
}

void koenigsberger_ec_derivatives(double t, double* f, const grid_parms& grid, EC_cell** ec)
{
	koenigsberger_ec_derivatives_implicit(t, f, grid, ec);
	koenigsberger_ec_derivatives_explicit(t, f, grid, ec);
}

void koenigsberger_smc_implicit(const grid_parms& grid, SMC_cell** smc)
{
	// Evaluate single cell fluxes.
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {

			//Jvocc
			smc[i][j].fluxes[J_VOCC] = GCai * (smc[i][j].vars[smc_Vm] - vCa1) / (1 + ((double) (exp(((-1) * (smc[i][j].vars[smc_Vm] - vCa2)) / RCai))));
			//J Na/Ca
			smc[i][j].fluxes[J_Na_Ca] = GNaCai * smc[i][j].vars[smc_Ca] * (smc[i][j].vars[smc_Vm] - vNaCai) / (smc[i][j].vars[smc_Ca] + cNaCai);
			//JNa/K
			smc[i][j].fluxes[J_Na_K] = FNaK;
			//J Cl
			smc[i][j].fluxes[J_Cl] = GCli * (smc[i][j].vars[smc_Vm] - vCl);
			//JK
			smc[i][j].fluxes[J_K] = GKi * smc[i][j].vars[smc_w] * (smc[i][j].vars[smc_Vm] - vKi);
		}
	}
}

void koenigsberger_smc_explicit(const grid_parms& grid, SMC_cell** smc)
{
	// Evaluate single cell fluxes.
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {

			//JIP3
			smc[i][j].fluxes[J_IP3] = (Fi * P2(smc[i][j].vars[smc_IP3])) / (P2(Kri) + P2(smc[i][j].vars[smc_IP3]));
			//JSRuptake
			smc[i][j].fluxes[J_SERCA] = (Bi * P2(smc[i][j].vars[smc_Ca])) / (P2(
					smc[i][j].vars[smc_Ca]) + P2(cbi));
			//Jcicr
			smc[i][j].fluxes[J_CICR] = (CICRi * (P2(smc[i][j].vars[smc_SR]) * P4(
					smc[i][j].vars[smc_Ca]))) / ((P2(sci) + P2(smc[i][j].vars[smc_SR]))
			* (P4(cci)+ P4(
							smc[i][j].vars[smc_Ca])));
			//Jextrusion
			smc[i][j].fluxes[J_Extrusion] = Di * smc[i][j].vars[smc_Ca] * (1 + ((smc[i][j].vars[smc_Vm] - vdi) / Rdi));

			//Jleak
			smc[i][j].fluxes[J_Leak] = Li * smc[i][j].vars[smc_SR];
			//Jvocc
			smc[i][j].fluxes[J_VOCC] = GCai * (smc[i][j].vars[smc_Vm] - vCa1) / (1 + ((double) (exp(((-1) * (smc[i][j].vars[smc_Vm] - vCa2)) / RCai))));
			//J Na/Ca
			smc[i][j].fluxes[J_Na_Ca] = GNaCai * smc[i][j].vars[smc_Ca] * (smc[i][j].vars[smc_Vm] - vNaCai) / (smc[i][j].vars[smc_Ca] + cNaCai);
			//Kactivation
			smc[i][j].fluxes[K_activation] = P2(smc[i][j].vars[smc_Ca] + cwi) / (P2(
					smc[i][j].vars[smc_Ca] + cwi) + (beta * ((double) exp((-1) * (smc[i][j].vars[smc_Vm] - vCa3) / RKi))));
			//Jdegradation
			smc[i][j].fluxes[J_IP3_deg] = ki * smc[i][j].vars[smc_IP3];
		}
	}
}

void koenigsberger_smc_derivatives_implicit(double* f, const grid_parms& grid, SMC_cell** smc)
{
	int k;
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;

			f[k + ((j - 1) * grid.neq_smc) + smc_Vm] = gama
					* (-smc[i][j].fluxes[J_Na_K] - smc[i][j].fluxes[J_Cl] - (2 * smc[i][j].fluxes[J_VOCC]) - smc[i][j].fluxes[J_Na_Ca] - smc[i][j].fluxes[J_K])
					+ smc[i][j].homo_fluxes[cpl_Vm] + smc[i][j].hetero_fluxes[cpl_Vm];

#if ! EXPLICIT_ONLY
			f[k + ((j - 1) * grid.neq_smc) + smc_Ca] = 0.0;
			f[k + ((j - 1) * grid.neq_smc) + smc_SR] = 0.0;
			f[k + ((j - 1) * grid.neq_smc) + smc_w] = 0.0;
			f[k + ((j - 1) * grid.neq_smc) + smc_IP3] = 0.0;
#endif
		}
	}
}

void koenigsberger_smc_derivatives_explicit(double* f, const grid_parms& grid, SMC_cell** smc)
{
	int k;
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;
			f[k + ((j - 1) * grid.neq_smc) + smc_Ca] = smc[i][j].fluxes[J_IP3] - smc[i][j].fluxes[J_SERCA] + smc[i][j].fluxes[J_CICR] - smc[i][j].fluxes[J_Extrusion]
					+ smc[i][j].fluxes[J_Leak] - smc[i][j].fluxes[J_VOCC] + smc[i][j].fluxes[J_Na_Ca] + smc[i][j].homo_fluxes[cpl_Ca] + smc[i][j].hetero_fluxes[cpl_Ca];

			f[k + ((j - 1) * grid.neq_smc) + smc_SR] = smc[i][j].fluxes[J_SERCA] - smc[i][j].fluxes[J_CICR] - smc[i][j].fluxes[J_Leak];

			f[k + ((j - 1) * grid.neq_smc) + smc_w] = lambda * (smc[i][j].fluxes[K_activation] - smc[i][j].vars[smc_w]);

			f[k + ((j - 1) * grid.neq_smc) + smc_IP3] = -smc[i][j].fluxes[J_IP3_deg] + smc[i][j].homo_fluxes[cpl_IP3] + smc[i][j].hetero_fluxes[cpl_IP3];

#if ! EXPLICIT_ONLY
			f[k + ((j - 1) * grid.neq_smc) + smc_Vm] = 0.0;
#endif

		}
	}
}

void koenigsberger_ec_implicit(const grid_parms& grid, EC_cell** ec)
{
	// Evaluate single cell fluxes.
	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			//Total K flux
			ec[i][j].fluxes[J_Ktot] = Gtot * (ec[i][j].vars[ec_Vm] - vKj) * (ec[i][j].fluxes[J_BK_Ca] + ec[i][j].fluxes[J_SK_Ca]);
			//Residual currents
			ec[i][j].fluxes[J_Residual] = GRj * (ec[i][j].vars[ec_Vm] - vrestj);
		}
	}
}

void koenigsberger_ec_explicit(const grid_parms& grid, EC_cell** ec, int atp_timestep)
{
	// Evaluate single cell fluxes.
	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			//JIP3
			ec[i][j].fluxes[J_IP3] = (Fj * P2(ec[i][j].vars[ec_IP3])) / (P2(Krj) + P2(ec[i][j].vars[ec_IP3]));
			//JSRuptake
			ec[i][j].fluxes[J_SERCA] = (Bj * P2(ec[i][j].vars[ec_Ca])) / (P2(ec[i][j].vars[ec_Ca]) + P2(cbj));
			//Jcicr
			ec[i][j].fluxes[J_CICR] = (CICRj * P2(ec[i][j].vars[ec_SR]) * P4(ec[i][j].vars[ec_Ca]))
					/ ((P2(ec[i][j].vars[ec_SR]) + P2(scj))*(P4(ec[i][j].vars[ec_Ca])+P4(ccj)));
			//Jextrusion
			ec[i][j].fluxes[J_Extrusion] = Dj * ec[i][j].vars[ec_Ca];
			//Jleak
			ec[i][j].fluxes[J_Leak] = Lj * ec[i][j].vars[ec_SR];

			//IP3 degradation
			ec[i][j].fluxes[J_IP3_deg] = kDeg * ec[i][j].vars[ec_IP3];

			//J_NonSelective Cation channels
			ec[i][j].fluxes[J_NSC] = (Gcatj * (ECa - ec[i][j].vars[ec_Vm]) * 0.5)
					* (1 + ((double) (tanh((double) (((double) (log10((double) (ec[i][j].vars[ec_Ca]))) - m3cat) / m4cat)))));
			//BK_channels
			ec[i][j].fluxes[J_BK_Ca] =
					(0.4 / 2)
							* (1
									+ (double) (tanh(
											(double) ((((((double) (log10((double) (ec[i][j].vars[ec_Ca]))) - c1j) * (ec[i][j].vars[ec_Vm] - bj)) - a1j)
													/ ((m3b * (P2(ec[i][j].vars[ec_Vm]+(a2j*((double) (log10 ((double) (ec[i][j].vars[ec_Ca]))) - c1j)) - bj))) + m4b))))));
			//SK_channels
			ec[i][j].fluxes[J_SK_Ca] = (0.6 / 2) * (1 + (double) (tanh((double) (((double) (log10((double) (ec[i][j].vars[ec_Ca]))) - m3s) / m4s))));
			//Grouping all other trivial Ca fluxes
			ec[i][j].fluxes[J_trivial_Ca] = J0j;
			// Ratio of bound to total P2Y
			ec[i][j].fluxes[L_P_P2Y] = ec[i][j].JPLC[atp_timestep] / (ec[i][j].JPLC[atp_timestep] + kATP); // TODO: Does this need to be calculated every time? is JPLC constant?

			// Rate of PIP2 hydrolysis.
			ec[i][j].fluxes[R_PIP2_H] = alpha_j * (ec[i][j].vars[ec_Ca] / (ec[i][j].vars[ec_Ca] + KCa)) * ec[i][j].vars[ec_Gprot];

			// Induced IP3 influx
			ec[i][j].fluxes[J_ind_I] = (ec[i][j].fluxes[R_PIP2_H] * cons_PIP2 * unitcon_a) / (N_a * V_ec);

		}
	}
}

void koenigsberger_ec_derivatives_implicit(double t, double* f, const grid_parms& grid, EC_cell** ec)
{
	int k, offset = (grid.neq_smc * grid.num_smc_circumferentially * grid.num_smc_axially);
	for(int i = 1; i <= grid.num_ec_circumferentially; i++)
	{
		for(int j = 1; j <= grid.num_ec_axially; j++)
		{
			if (i > 1)
				k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
				k = offset + 0;

			f[k + ((j - 1) * grid.neq_ec) + ec_Vm] =
					((-1 / Cmj) * (ec[i][j].fluxes[J_Ktot] + ec[i][j].fluxes[J_Residual])) + ec[i][j].homo_fluxes[cpl_Vm] + ec[i][j].hetero_fluxes[cpl_Vm];


#if !EXPLICIT_ONLY

			f[k + ((j - 1) * grid.neq_ec) + ec_Ca] = 0.0;
			f[k + ((j - 1) * grid.neq_ec) + ec_SR] = 0.0;
			f[k + ((j - 1) * grid.neq_ec) + ec_IP3] = 0.0;
			f[k + ((j - 1) * grid.neq_ec) + ec_Gprot] = 0.0;
#endif
		}
	}
}

void koenigsberger_ec_derivatives_explicit(double t, double* f, const grid_parms& grid, EC_cell** ec)
{
	int k, offset = (grid.neq_smc * grid.num_smc_circumferentially * grid.num_smc_axially);
	for(int i = 1; i <= grid.num_ec_circumferentially; i++)
	{
		for(int j = 1; j <= grid.num_ec_axially; j++)
		{
			if (i > 1)
				k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
				k = offset + 0;

			f[k + ((j - 1) * grid.neq_ec) + ec_Ca] =
					ec[i][j].fluxes[J_IP3] - ec[i][j].fluxes[J_SERCA] + ec[i][j].fluxes[J_CICR] - ec[i][j].fluxes[J_Extrusion]
					+ ec[i][j].fluxes[J_Leak] + ec[i][j].fluxes[J_NSC] + ec[i][j].fluxes[J_trivial_Ca] + ec[i][j].homo_fluxes[cpl_Ca] + ec[i][j].hetero_fluxes[cpl_Ca];

			f[k + ((j - 1) * grid.neq_ec) + ec_SR] =
					ec[i][j].fluxes[J_SERCA] - ec[i][j].fluxes[J_CICR] - ec[i][j].fluxes[J_Leak];

			f[k + ((j - 1) * grid.neq_ec) + ec_IP3] =
						ec[i][j].fluxes[J_ind_I] - ec[i][j].fluxes[J_IP3_deg] + ec[i][j].homo_fluxes[cpl_IP3] + ec[i][j].hetero_fluxes[cpl_IP3];


			f[k + ((j - 1) * grid.neq_ec) + ec_Gprot] =
						(K_G_prot_act * (delta + ec[i][j].fluxes[L_P_P2Y]) * (G_prot_tot - ec[i][j].vars[ec_Gprot]))
						- (K_G_prot_deact * ec[i][j].vars[ec_Gprot]);

#if ! EXPLICIT_ONLY
			f[k + ((j - 1) * grid.neq_ec) + ec_Vm] = 0.0;
#endif

		}
	}
}
