#ifndef KOENIGSBERGER_CONSTANTS_H_
#define KOENIGSBERGER_CONSTANTS_H_

/*
 * koenigsberger_constants.h
 *
 * FIXED MODEL PARAMETERS FOR THE KOENIGSBERGER MODEL
 *
 */

//SMC state variables
#define     smc_Ca              0
#define     smc_SR              1
#define     smc_Vm              2
#define		smc_w				3
#define     smc_IP3             4
//EC state variables
#define     ec_Ca               0
#define     ec_SR               1
#define     ec_Vm               2
#define     ec_IP3              3
// Lemon et al. specific
#define 	ec_Gprot 			4

//Ionic currents common to both celltypes
#define     J_IP3               0
#define     J_SERCA             1
#define     J_CICR              2
#define     J_Extrusion         3
#define     J_Leak              4
#define     J_IP3_deg           5
//SMC specific currents
#define     J_VOCC              6
#define     J_Na_Ca             7
#define     J_Na_K              8
#define     J_Cl	            9
#define     J_K                 10
#define     K_activation        11
//EC specific currents
#define     J_NSC               6
#define     J_BK_Ca             7
#define     J_SK_Ca             8
#define     J_Ktot              9
#define     J_Residual          10
#define     J_trivial_Ca        11

// Lemon et al. specific
#define		J_ind_I 			12
#define 	L_P_P2Y 			13
#define 	R_PIP2_H			14
// Bennett et al. specific
#define		J_ind_I 			12
#define 	B_P_P2Y 			13
#define 	J_Gprot 			14


#define	    cpl_Ca		        0
#define	    cpl_Vm		        1
#define	    cpl_IP3		        2

void initialize_koenigsberger_smc(const grid_parms&, double*, SMC_cell**);
void initialize_koenigsberger_ec(const grid_parms&, double*, EC_cell**);

void koenigsberger_smc(const grid_parms&, SMC_cell**);
void koenigsberger_smc_derivatives(double*, const grid_parms&, SMC_cell**);
void koenigsberger_ec(const grid_parms&, EC_cell**);
void koenigsberger_ec_derivatives(double, double*, const grid_parms&, EC_cell**);

void koenigsberger_smc_implicit(const grid_parms&, SMC_cell**);
void koenigsberger_smc_derivatives_implicit(double*, const grid_parms&, SMC_cell**);
void koenigsberger_ec_implicit(const grid_parms&, EC_cell**);
void koenigsberger_ec_derivatives_implicit(double, double*, const grid_parms&, EC_cell**);

void koenigsberger_smc_explicit(const grid_parms&, SMC_cell**);
void koenigsberger_smc_derivatives_explicit(double*, const grid_parms&, SMC_cell**);
void koenigsberger_ec_explicit(const grid_parms&, EC_cell**);
void koenigsberger_ec_derivatives_explicit(double, double*, const grid_parms&, EC_cell**);

#endif /* KOENIGSBERGER_CONSTANTS_H_ */
