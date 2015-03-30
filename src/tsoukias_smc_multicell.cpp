/*
 * tsoukias_smc_multicell.cpp
 *
 *  Created on: 16/12/2012
 *      Author: mohsinshaikh
 */
#include "computelib.h"
#include "tsoukias_constants.h"
/*******************************************************************************************/
void Initialize_tsoukias_smc(grid_parms grid, double y[], SMC_cell** smc)
/*******************************************************************************************/
{

	int k = 0, offset;

		for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
			for (int j = 1; j <= grid.num_smc_axially; j++) {
				if (i > 1)
					k = ((i - 1) * grid.neq_smc_axially);
				else if (i == 1)
					k = 0;

				y[k + ((j - 1) * grid.neq_smc) + smc_Vm] 	= -59.4;		//mV
				y[k + ((j - 1) * grid.neq_smc) + smc_Ca]	= 68-6; 		//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_Ca_u]	= 0.66;			//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_Ca_r]	= 0.57;			//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_Na_i]	= 8.4;			//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_K_i]	= 140.0;		//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_Cl_i] 	= 59.4;			//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_IP3]	= 0.0;
				y[k + ((j - 1) * grid.neq_smc) + smc_DAG]	= 0.0;
				y[k + ((j - 1) * grid.neq_smc) + smc_cGMP_i]= 0.0;
				y[k + ((j - 1) * grid.neq_smc) + smc_V_cGMP]= 0.0;
				y[k + ((j - 1) * grid.neq_smc) + smc_d_L]   = 1 / (1 + exp(-y[k + ((j - 1) * grid.neq_smc) + smc_Vm] / 8.3));
				y[k + ((j - 1) * grid.neq_smc) + smc_f_L]   = 1 / (1 + exp((y[k + ((j - 1) * grid.neq_smc) + smc_Vm] + 42.0) / 9.1));

				smc[i][j].A[V_half_KCa] = -41.7 * log10(y[k + ((j-1) * grid.neq_smc) + smc_Ca]) - 128.2;
				y[k + ((j - 1) * grid.neq_smc) + smc_p_f]   = 1 / (1 + exp(-(y[k + ((j - 1) * grid.neq_smc) + smc_Vm] - smc[i][j].A[V_half_KCa]) / 18.25));


				y[k + ((j - 1) * grid.neq_smc) + smc_p_K] 	= 1 / (1 + exp(-(y[k + ((j - 1) * grid.neq_smc) + smc_Vm] + 11.0) / 15.0));
				y[k + ((j - 1) * grid.neq_smc) + smc_q_1] 	= 1 / (1 + exp((y[k + ((j - 1) * grid.neq_smc) + smc_Vm] + 40.0) / 14.0));
				y[k + ((j - 1) * grid.neq_smc) + smc_q_2] 	= 1 / (1 + exp((y[k + ((j - 1) * grid.neq_smc) + smc_Vm] + 40.0) / 14.0));
				y[k + ((j - 1) * grid.neq_smc) + smc_R_01]  = 0.995;
				y[k + ((j - 1) * grid.neq_smc) + smc_R_10]  = 0.0033;
				y[k + ((j - 1) * grid.neq_smc) + smc_R_11]  = 4.0e-6;
				y[k + ((j - 1) * grid.neq_smc) + smc_h_IP3] = K_inhIP3 / (y[k + ((j - 1) * grid.neq_smc) + smc_Ca] + K_inhIP3);

				y[k + ((j - 1) * grid.neq_smc) + smc_PIP2] = PIP2_T - (1 + (k_deg_G / r_r_G)) * gama_G * y[k + ((j - 1) * grid.neq_smc) + smc_IP3];

				smc[i][j].A[r_h_G] 	=   k_deg_G * gama_G * y[k + ((j - 1) * grid.neq_smc) + smc_IP3] / y[k + ((j - 1) * grid.neq_smc) + smc_PIP2];

				y[k + ((j - 1) * grid.neq_smc) + smc_R_S_G] = R_T_G * ksi_G;
				y[k + ((j - 1) * grid.neq_smc) + smc_R_S_P_G] = 0.0;
				y[k + (j * grid.neq_smc) + smc_G] = smc[i][j].A[r_h_G] * (K_c_G *  y[k + ((j - 1) * grid.neq_smc) + smc_Ca]) / (alpha_G * y[k + ((j - 1) * grid.neq_smc) + smc_Ca]);

				smc[i][j].A[delta_G] =   k_d_G * y[k + ((j - 1) * grid.neq_smc) + smc_G] / (k_a_G * (G_T_G - y[k + ((j - 1) * grid.neq_smc) + smc_G]));

				smc[i][j].NE = 0.0;	//mM
				smc[i][j].NO = 0.0;	//nM
				smc[i][j].I_stim = 0.0;	//pA
			}
		}

		for (int i = 0; i < (grid.num_smc_circumferentially + grid.num_ghost_cells);
				i++) {
			for (int j = 0; j < (grid.num_smc_axially + grid.num_ghost_cells);
					j++) {
				smc[i][j].p[smc_Vm]				= 0.0;
				smc[i][j].p[smc_d_L] 			= 0.0;
				smc[i][j].p[smc_f_L] 		 	= 0.0;
				smc[i][j].p[smc_p_f] 			= 0.0;
				smc[i][j].p[smc_p_s] 			= 0.0;
				smc[i][j].p[smc_q_1] 			= 0.0;
				smc[i][j].p[smc_q_2] 			= 0.0;
				smc[i][j].p[smc_p_K] 			= 0.0;
				smc[i][j].p[smc_Ca_u] 			= 0.0;
				smc[i][j].p[smc_Ca_r] 			= 0.0;
				smc[i][j].p[smc_R_10] 			= 0.0;
				smc[i][j].p[smc_R_11] 			= 0.0;
				smc[i][j].p[smc_R_01] 			= 0.0;
				smc[i][j].p[smc_h_IP3] 			= 0.0;
				smc[i][j].p[smc_R_S_G] 			= 0.0;
				smc[i][j].p[smc_R_S_P_G]		= 0.0;
				smc[i][j].p[smc_G] 				= 0.0;
				smc[i][j].p[smc_IP3] 			= 0.0;
				smc[i][j].p[smc_PIP2] 			= 0.0;
				smc[i][j].p[smc_V_cGMP]			= 0.0;
				smc[i][j].p[smc_cGMP_i]			= 0.0;
				smc[i][j].p[smc_Ca]				= 0.0;
				smc[i][j].p[smc_Na_i] 			= 0.0;
				smc[i][j].p[smc_K_i] 			= 0.0;
				smc[i][j].p[smc_Cl_i]			= 0.0;
				smc[i][j].p[smc_DAG]			= 0.0;

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
/************************************************************************************************/
void tsoukias_smc(grid_parms grid, SMC_cell** smc)
/*																				   		      	*
 * This is the multicell version of the function tosukias_smc() used for single   		  	  	*
 * EC-SMC Unit simulation.														  		  		*
 * int i and j are the indicies of the SMC being simulated.						  		  		*
 * Int NO_Path and cGMP_Path are the two switches to enable or disable the		  		  		*
 * corresponding pathways.														  		  		*
 ************************************************************************************************/
{

// --------------------------------------------------------------------------------------------------------
//   SYSTEM OF EQNS OF A SINGLE VASCULAR SMOOTH MUSCLE CELL   E L E C T R O P H Y S I O L O G Y   M O D E L
// --------------------------------------------------------------------------------------------------------
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {

			smc[i][j].p[smc_DAG] = smc[i][j].p[smc_IP3];

//SMC single cell ionic currents

			smc[i][j].A[E_K] = (R_const * T / (z_K * F)) * log(K_o / smc[i][j].p[smc_K_i]);
			smc[i][j].A[E_Ca] = (R_const * T / (z_Ca * F)) * log(Ca_o / smc[i][j].p[smc_Ca]);
			smc[i][j].A[E_Na] = (R_const * T / (z_Na * F)) * log(Na_o / smc[i][j].p[smc_Na_i]);
			smc[i][j].A[E_Cl] = (R_const * T / (z_Cl * F)) * log(Cl_o / smc[i][j].p[smc_Cl_i]);

/**** Membrane Electrophysiology ****/

			//L-type voltage operated Ca channels

			smc[i][j].A[dbar_L] = 1 / (1 + exp(-smc[i][j].p[smc_Vm] / 8.3));
			smc[i][j].A[tau_dL] = 1.15 + 2.5 * exp(-P2((smc[i][j].p[smc_Vm] + 40.0) / 30.0));
			smc[i][j].A[fbar_L] = 1 / (1 + exp((smc[i][j].p[smc_Vm] + 42.0) / 9.1));
			smc[i][j].A[tau_fL] = 45.0 + 65.0 * exp(-P2((smc[i][j].p[smc_Vm] + 35.0) / 25.0));
			smc[i][j].A[I_VOCC] =
			    1e6 * Am * P_VOCC * smc[i][j].p[smc_d_L] *
			    smc[i][j].p[smc_f_L] * smc[i][j].p[smc_Vm] *
			    (P2(z_Ca * F) / (R_const * T)) *
			    ((Ca_o -
			      (smc[i][j].p[smc_Ca] *
			       exp((smc[i][j].p[smc_Vm] * z_Ca * F) /
				   (R_const * T)))) / (1 - exp((smc[i][j].p[smc_Vm] * z_Ca * F) / (R_const * T))));

			//Large conductance Ca activated potassium channels

			smc[i][j].A[R_NO] = grid.NO_path * smc[i][j].NO / (smc[i][j].NO + 200.00);
			smc[i][j].A[R_cGMP] = grid.cGMP_path * P2(smc[i][j].p[smc_cGMP_i]) / (P2(smc[i][j].p[smc_cGMP_i]) + P2(1.5));
			smc[i][j].A[V_half_KCa] = -41.7 * log10(smc[i][j].p[smc_Ca]) - 128.2;

			smc[i][j].A[pbar_o] = 1 / (1 + exp(-(smc[i][j].p[smc_Vm] - smc[i][j].A[V_half_KCa]) / 18.25));
			smc[i][j].A[P_KCa] = 0.17 * smc[i][j].p[smc_p_f] + 0.83 * smc[i][j].p[smc_p_s];
			smc[i][j].A[i_KCa] =
			    1e6 * P_BKCa * smc[i][j].p[smc_Vm] * (P2(F) /
								   (R_const * T)) *
			    ((K_o - (smc[i][j].p[smc_K_i] * exp(smc[i][j].p[smc_Vm] * F / (R_const * T)))) / (1 - exp(smc[i]
														 [j].p[smc_Vm]
														 * F / (R_const * T))));
			smc[i][j].A[I_BKCa] = Am * N_BKCa * smc[i][j].A[P_KCa] * smc[i][j].A[i_KCa];

			//Voltage dependent K channels

			smc[i][j].A[qbar] = 1.0 / (1 + exp((smc[i][j].p[smc_Vm] + 40.0) / 14.0));
			smc[i][j].A[tau_pK] = 61.5 * exp(-0.027 * smc[i][j].p[smc_Vm]);
			smc[i][j].A[pbar_K] = 1.0 / (1 + exp(-(smc[i][j].p[smc_Vm] + 11.0) / 15.0));
			smc[i][j].A[I_Kv] =
			    g_Kv * smc[i][j].p[smc_p_K] * (0.45 *
							   smc[i][j].p[smc_q_1] +
							   0.55 * smc[i][j].p[smc_q_2]) * (smc[i][j].p[smc_Vm] - smc[i][j].A[E_K]);

			//Unspecified K Leak channels
			smc[i][j].A[I_Kleak] = g_K_leak * (smc[i][j].p[smc_Vm] - smc[i][j].A[E_K]);

			//Non-selective cation channels
			smc[i][j].A[Po_NSC] = 0.4344 + ((1 - 0.4344) / (1 + exp(-(smc[i][j].p[smc_Vm] - 47.12)
										/ 24.24)));
			smc[i][j].A[INa_NSC] =
			    1e6 * Am *
			    ((smc[i][j].p[smc_DAG] /
			      (smc[i][j].p[smc_DAG] + K_NSC)) +
			     d_NSCmin) * smc[i][j].A[Po_NSC] * P_NaNSC *
			    smc[i][j].p[smc_Vm] * (P2(F) / (R_const * T)) *
			    ((Na_o - (smc[i][j].p[smc_Na_i] * exp(smc[i][j].p[smc_Vm] * F / (R_const * T)))) / (1 - exp(smc[i]
														   [j].p[smc_Vm]
														   * F / (R_const * T))));
			smc[i][j].A[IK_NSC] =
			    1e6 * Am *
			    ((smc[i][j].p[smc_DAG] /
			      (smc[i][j].p[smc_DAG] + K_NSC)) +
			     d_NSCmin) * smc[i][j].A[Po_NSC] * P_KNSC *
			    smc[i][j].p[smc_Vm] * (P2(F) / (R_const * T)) *
			    ((K_o - (smc[i][j].p[smc_K_i] * exp(smc[i][j].p[smc_Vm] * F / (R_const * T)))) / (1 - exp(smc[i]
														 [j].p[smc_Vm]
														 * F / (R_const * T))));
			smc[i][j].A[ICa_NSC] =
			    1e6 * Am * d_NSCmin * smc[i][j].A[Po_NSC] *
			    P_CaNSC * smc[i][j].p[smc_Vm] * (P2(z_Ca * F) /
							      (R_const * T)) *
			    ((Ca_o -
			      (smc[i][j].p[smc_Ca] *
			       exp(smc[i][j].p[smc_Vm] * z_Ca * F /
				   (R_const * T)))) / (1 - exp(smc[i][j].p[smc_Vm] * z_Ca * F / (R_const * T))));
			smc[i][j].A[I_NSC] = smc[i][j].A[INa_NSC] + smc[i][j].A[IK_NSC] + smc[i][j].A[ICa_NSC];

			//Store operated non selective cation channels
			smc[i][j].A[P_SOC] = 1.0 / (1 + (smc[i][j].p[smc_Ca_u] / K_SOC));
			smc[i][j].A[ICa_SOC] = g_SOCCa * smc[i][j].A[P_SOC] * (smc[i][j].p[smc_Vm] - smc[i][j].A[E_Ca]);
			smc[i][j].A[INa_SOC] = g_SOCNa * smc[i][j].A[P_SOC] * (smc[i][j].p[smc_Vm] - smc[i][j].A[E_Na]);
			smc[i][j].A[I_SOC] = smc[i][j].A[ICa_SOC] + smc[i][j].A[INa_SOC];

			//Calcium activated Chloride channels
			smc[i][j].A[alpha_Cl] =
			    pow(smc[i][j].p[smc_cGMP_i],
				n_ClcGMP) / (pow(smc[i][j].p[smc_cGMP_i], n_ClcGMP) + pow(K_ClcGMP, n_ClcGMP));
			smc[i][j].A[K_ClCacGMP] = (1 - 0.9 * smc[i][j].A[alpha_Cl]) * 400.00;
			smc[i][j].A[P_Cl] =
			    (0.0132 * pow(smc[i][j].p[smc_Ca] * 1e6, n_ClCa) /
			     (pow(smc[i][j].p[smc_Ca] * 1e6, n_ClCa) +
			      pow(K_ClCa,
				  n_ClCa))) +
			    (smc[i][j].A[alpha_Cl] *
			     pow(smc[i][j].p[smc_Ca] * 1e6,
				 n_ClCa) / (pow(smc[i][j].p[smc_Ca] * 1e6, n_ClCa) + pow(smc[i][j].A[K_ClCacGMP], n_ClCa)));
			smc[i][j].A[I_ClCa] = C_m * g_ClCa * smc[i][j].A[P_Cl] * (smc[i][j].p[smc_Vm] - smc[i][j].A[E_Cl]);

			//Plasma membrane Ca pump
			smc[i][j].A[I_PMCA] = I_bar_PMCA * (smc[i][j].p[smc_Ca]) / ((smc[i][j].p[smc_Ca]) + K_mPMCA);

			//Plasma membrane Na-Ca exchanger
			smc[i][j].A[phi_F] = exp((F * gama_NCX * smc[i][j].p[smc_Vm]) / (R_const * T));
			smc[i][j].A[phi_R] = exp((F * (gama_NCX - 1) * smc[i][j].p[smc_Vm]) / (R_const * T));
			smc[i][j].A[R_NCX_cGMP] = 1 + 0.55 * (smc[i][j].p[smc_cGMP_i] / (smc[i][j].p[smc_cGMP_i] + 45.0));
			smc[i][j].A[I_NCX] =
			    g_NCX * smc[i][j].A[R_NCX_cGMP] *
			    (((P3(smc[i][j].p[smc_Na_i]) * Ca_o *
			       smc[i][j].A[phi_F]) - (P3(Na_o) * smc[i][j].p[smc_Ca] * smc[i][j].A[phi_R]))
			     / (1 + d_NCX * (P3(Na_o) * smc[i][j].p[smc_Ca] + P3(smc[i][j].p[smc_Na_i]) * Ca_o)));
			//Sodium potassium pump
			smc[i][j].A[Q_10] = 1.87;
			smc[i][j].A[Q] = pow(smc[i][j].A[Q_10], ((T - 309.15) / 10.0));
			smc[i][j].A[I_NaK] =
			    C_m * I_bar_NaK * smc[i][j].A[Q] *
			    (pow(K_o, n_HKo) /
			     (pow(K_o, n_HKo) +
			      pow(K_dKo, n_HKo))) * (pow(smc[i][j].p[smc_Na_i],
							 n_HNai) /
						     (pow
						      (smc[i][j].p[smc_Na_i],
						       n_HNai) + pow(K_dNai,
								     n_HNai))) *
			    ((smc[i][j].p[smc_Vm] + 150.0) / (smc[i][j].p[smc_Vm] + 200.0));

			//Sodium potassium chloride cotransport
			smc[i][j].A[R_NaKCl_cGMP] = 1.0 + 3.5 * (smc[i][j].p[smc_cGMP_i] / (smc[i][j].p[smc_cGMP_i] + 6.4));
			smc[i][j].A[I_Cl_NaKCl] =
			    -(smc[i][j].A[R_NaKCl_cGMP] * z_Cl * Am * L_NaKCl *
			      R_const * F * T) * log((Na_o / smc[i][j].p[smc_Na_i]) *
					       (K_o / smc[i][j].p[smc_K_i]) * P2(Cl_o / smc[i][j].p[smc_Cl_i]));
			smc[i][j].A[I_Na_NaKCl] = -0.5 * smc[i][j].A[I_Cl_NaKCl];
			smc[i][j].A[I_K_NaKCl] = -0.5 * smc[i][j].A[I_Cl_NaKCl];

	/**** Sarcoplasmic Reticulum ****/
			smc[i][j].A[I_SERCA] = I_bar_SERCA * (smc[i][j].p[smc_Ca] * 1e3 / (smc[i][j].p[smc_Ca] * 1e3 + K_mUp));
			smc[i][j].A[I_tr] = (smc[i][j].p[smc_Ca_u] - smc[i][j].p[smc_Ca_r]) * z_Ca * F * vol_u / tau_tr;
			smc[i][j].A[I_rel] =
			    (P2(smc[i][j].p[smc_R_10]) +
			     R_leak) * (smc[i][j].p[smc_Ca_r] - smc[i][j].p[smc_Ca]) * z_Ca * F * vol_r / tau_rel;

			// Ryanodine receptor
			smc[i][j].A[R_00] = 1.0 - smc[i][j].p[smc_R_01] - smc[i][j].p[smc_R_10] - smc[i][j].p[smc_R_11];

			// IP3 Receptor
			smc[i][j].A[I_IP3] =
			    I_bar_IP3 * z_Ca * vol_Ca * F *
			    P3(smc[i][j].p[smc_h_IP3] *
			       (smc[i][j].p[smc_IP3] /
				(smc[i][j].p[smc_IP3] +
				 K_IP3)) * (smc[i][j].p[smc_Ca] /
					    (smc[i][j].p[smc_Ca] + K_actIP3))) * (smc[i][j].p[smc_Ca_u] - smc[i][j].p[smc_Ca]);

	/**** Alpha Adrenoreceptor Activation and IP3 formation ****/
			smc[i][j].A[r_h_G] =
			    alpha_G * (smc[i][j].p[smc_Ca] / (smc[i][j].p[smc_Ca] + K_c_G)) * smc[i][j].p[smc_G];
			smc[i][j].A[rho_r_G] = smc[i][j].NE * smc[i][j].p[smc_R_S_G] / (ksi_G * R_T_G * (K_1_G + smc[i][j].NE));
			smc[i][j].A[delta_G] = 0.0;	//deltaG_0;//0.0;

	/**** sGC activation and cGMP formation ****/
			smc[i][j].A[A0_sGC] = (((k_1_sGC + k2_sGC) * kD_sGC) + (k_1_sGC * k_2_sGC)) / (k1_sGC * k3_sGC);
			smc[i][j].A[A1_sGC] = ((k1_sGC + k3_sGC) * kD_sGC + (k2_sGC + k_2_sGC) * k1_sGC) / (k1_sGC * k3_sGC);
			smc[i][j].A[V_bar_cGMP] =
			    V_cGMPmax * (B5_sGC * smc[i][j].NO * 1e-6 +
					 P2(smc[i][j].NO * 1e-6)) /
			    (smc[i][j].A[A0_sGC] + smc[i][j].A[A1_sGC] * smc[i][j].NO * 1e-6 + P2(smc[i][j].NO * 1e-6));
			if ((smc[i][j].A[V_bar_cGMP] - smc[i][j].p[smc_V_cGMP]) >= 0) {
				smc[i][j].A[tau_sGC] = tau_a_sGC;
			} else {
				smc[i][j].A[tau_sGC] = tau_d_sGC;
			}

	/**** Ionic Balances ****/
			smc[i][j].A[I_Catotm] =
			    smc[i][j].A[ICa_SOC] + smc[i][j].A[I_VOCC] -
			    2 * smc[i][j].A[I_NCX] + smc[i][j].A[I_PMCA] + smc[i][j].A[ICa_NSC];
			smc[i][j].A[I_Natotm] =
			    smc[i][j].A[I_Na_NaKCl] + smc[i][j].A[INa_SOC] +
			    3 * smc[i][j].A[I_NaK] + 3 * smc[i][j].A[I_NCX] + smc[i][j].A[INa_NSC];
			smc[i][j].A[I_Ktotm] =
			    smc[i][j].A[I_K_NaKCl] + smc[i][j].A[I_Kv] +
			    smc[i][j].A[I_BKCa] + smc[i][j].A[IK_NSC] + smc[i][j].A[I_Kleak] - 2 * smc[i][j].A[I_NaK];
			smc[i][j].A[I_Cltotm] = smc[i][j].A[I_Cl_NaKCl] + smc[i][j].A[I_ClCa];

		}		//end j

	}			//end i

}
/**************************************************************************************/
/**/ void tsoukias_smc_derivatives(double* f, grid_parms grid,SMC_cell** smc) /**/
/**************************************************************************************/
{
	int k;
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;

			f[k + ((j - 1) * grid.neq_smc) + smc_Vm] = (-1e3
					* (smc[i][j].A[I_VOCC] + smc[i][j].A[I_Kv]
							+ smc[i][j].A[I_BKCa] + smc[i][j].A[I_Kleak]
							+ smc[i][j].A[I_NSC] + smc[i][j].A[I_SOC]
							+ smc[i][j].A[I_ClCa] + smc[i][j].A[I_PMCA]
							+ smc[i][j].A[I_NaK] + smc[i][j].A[I_NCX]) / C_m)
					+ smc[i][j].B[cpl_Vm] + smc[i][j].C[cpl_Vm];

			f[k + ((j - 1) * grid.neq_smc) + smc_d_L] = 1e3
					* (smc[i][j].A[dbar_L] - smc[i][j].p[smc_d_L])
					/ smc[i][j].A[tau_dL];
			f[k + ((j - 1) * grid.neq_smc) + smc_f_L] = 1e3
					* (smc[i][j].A[fbar_L] - smc[i][j].p[smc_f_L])
					/ smc[i][j].A[tau_fL];

			f[k + ((j - 1) * grid.neq_smc) + smc_p_f] = 1e3
					* (smc[i][j].A[pbar_o] - smc[i][j].p[smc_p_f]) / tau_pf;
			f[k + ((j - 1) * grid.neq_smc) + smc_p_s] = 1e3
					* (smc[i][j].A[pbar_o] - smc[i][j].p[smc_p_s]) / tau_ps;

			f[k + ((j - 1) * grid.neq_smc) + smc_q_1] = 1e3
					* (smc[i][j].A[qbar] - smc[i][j].p[smc_q_1]) / tau_q1;
			f[k + ((j - 1) * grid.neq_smc) + smc_q_2] = 1e3
					* (smc[i][j].A[qbar] - smc[i][j].p[smc_q_2]) / tau_q2;
			f[k + ((j - 1) * grid.neq_smc) + smc_p_K] = 1e3
					* (smc[i][j].A[pbar_K] - smc[i][j].p[smc_p_K])
					/ smc[i][j].A[tau_pK];

			f[k + ((j - 1) * grid.neq_smc) + smc_Ca_u] = 1e3
					* (smc[i][j].A[I_SERCA] - smc[i][j].A[I_tr]
							- smc[i][j].A[I_IP3]) / (z_Ca * F * vol_u);
			f[k + ((j - 1) * grid.neq_smc) + smc_Ca_r] = 1e3
					* ((smc[i][j].A[I_tr] - smc[i][j].A[I_rel])
							/ (z_Ca * F * vol_r))
					/ ((1
							+ (CSQN * K_CSQN
									/ P2(K_CSQN + smc[i][j].p[smc_Ca_r]))));

			f[k + ((j - 1) * grid.neq_smc) + smc_R_10] =
					1e3
							* ((Kr1 * P2(smc[i][j].p[smc_Ca])
									* smc[i][j].A[R_00])
									- ((K_r1 + Kr2 * smc[i][j].p[smc_Ca])
											* smc[i][j].p[smc_R_10])
									+ (K_r2 * smc[i][j].p[smc_R_11]));
			f[k + ((j - 1) * grid.neq_smc) + smc_R_11] =
					1e3
							* ((Kr2 * smc[i][j].p[smc_Ca]
									* smc[i][j].p[smc_R_10])
									- ((K_r1 + K_r2) * smc[i][j].p[smc_R_11])
									+ (Kr1 * P2(smc[i][j].p[smc_Ca])
											* smc[i][j].p[smc_R_01]));
			f[k + ((j - 1) * grid.neq_smc) + smc_R_01] =
					1e3	* ((Kr2 * smc[i][j].p[smc_Ca] * smc[i][j].A[R_00])
									+ (K_r1 * smc[i][j].p[smc_R_11])
									- ((K_r2 + Kr1 * P2(smc[i][j].p[smc_Ca]))* smc[i][j].p[smc_R_01]));
			f[k + ((j - 1) * grid.neq_smc) + smc_h_IP3] = 1e3 * k_on_IP3
					* (K_inhIP3
							- (smc[i][j].p[smc_Ca] + K_inhIP3)
									* smc[i][j].p[smc_h_IP3]);

			f[k + ((j - 1) * grid.neq_smc) + smc_R_S_G] = 1e3
					* (k_r_G * ksi_G * R_T_G
							- (smc[i][j].p[smc_R_S_G]
									* (k_r_G
											+ (k_p_G * smc[i][j].NE
													/ (K_1_G + smc[i][j].NE))))
							- (k_r_G * smc[i][j].p[smc_R_S_P_G]));
			f[k + ((j - 1) * grid.neq_smc) + smc_R_S_P_G] = 1e3 * smc[i][j].NE
					* ((k_p_G * smc[i][j].p[smc_R_S_G] / (K_1_G + smc[i][j].NE))
							- (k_e_G * smc[i][j].p[smc_R_S_P_G]
									/ (K_2_G + smc[i][j].NE)));

			f[k + ((j - 1) * grid.neq_smc) + smc_G] = 1e3
					* (k_a_G
							* ((smc[i][j].A[delta_G] + smc[i][j].A[rho_r_G])
									* (G_T_G - smc[i][j].p[smc_G]))
							- (k_d_G * smc[i][j].p[smc_G]));
			f[k + ((j - 1) * grid.neq_smc) + smc_IP3] = 1e3
					* (((smc[i][j].A[r_h_G] / gama_G) * smc[i][j].p[smc_PIP2])
							- (k_deg_G * smc[i][j].p[smc_IP3]))
					+ smc[i][j].B[cpl_IP3] + smc[i][j].C[cpl_IP3];
			f[k + ((j - 1) * grid.neq_smc) + smc_PIP2] = 1e3
					* ((r_r_G * PIP2_T)
							- (smc[i][j].A[r_h_G] + r_r_G)
									* smc[i][j].p[smc_PIP2]
							- (r_r_G * gama_G * smc[i][j].p[smc_IP3]));

			f[k + ((j - 1) * grid.neq_smc) + smc_V_cGMP] = 1e3
					* (smc[i][j].A[V_bar_cGMP] - smc[i][j].p[smc_V_cGMP])
					/ smc[i][j].A[tau_sGC];
			f[k + ((j - 1) * grid.neq_smc) + smc_cGMP_i] = 1e3
					* (smc[i][j].p[smc_V_cGMP]
							- k_pde_cGMP * P2(smc[i][j].p[smc_cGMP_i])
									/ (smc[i][j].p[smc_cGMP_i] + K_m_pde));

			f[k + ((j - 1) * grid.neq_smc) + smc_Ca] =
					1e3
							* ((-(smc[i][j].A[I_Catotm] + smc[i][j].A[I_SERCA]
									- smc[i][j].A[I_rel] - smc[i][j].A[I_IP3])
									/ (z_Ca * F * vol_Ca))
									/ (1
											+ ((S_CM * K_d)
													/ P2(K_d + smc[i][j].p[smc_Ca]))+
											((B_F * K_dB) / P2(K_dB + smc[i][j].p[smc_Ca])))) + smc[i][j].B[cpl_Ca] + smc[i][j].C[cpl_Ca];

			f[k + ((j - 1) * grid.neq_smc) + smc_Na_i] = 1e3
					* (-smc[i][j].A[I_Natotm] / (F * vol_i));
			f[k + ((j - 1) * grid.neq_smc) + smc_K_i] = 1e3
					* (-smc[i][j].A[I_Ktotm] / (F * vol_i));
			f[k + ((j - 1) * grid.neq_smc) + smc_Cl_i] = 1e3
					* (-smc[i][j].A[I_Cltotm] / (z_Cl * F * vol_i));
			f[k + ((j - 1) * grid.neq_smc) + smc_DAG] = f[k
					+ ((j - 1) * grid.neq_smc) + smc_IP3];
		}
	}
}
