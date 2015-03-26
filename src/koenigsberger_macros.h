// SMC state variables.
#define     smc_Ca              0
#define     smc_SR              1
#define     smc_Vm              2
#define		smc_w				3
#define     smc_IP3             4
// EC state variables.
#define     ec_Ca               0
#define     ec_SR               1
#define     ec_Vm               2
#define     ec_IP3              3

// Ionic currents common to both cell types.
#define     J_IP3               0
#define     J_SERCA             1
#define     J_CICR              2
#define     J_Extrusion         3
#define     J_Leak              4
#define     J_IP3_deg           5

// SMC specific currents.
#define     J_VOCC              6
#define     J_Na_Ca             7
#define     J_Na_K              8
#define     J_Cl	            9
#define     J_K                 10
#define     K_activation        11

// EC specific currents.
#define     J_NSC               6
#define     J_BK_Ca             7
#define     J_SK_Ca             8
#define     J_Ktot              9
#define     J_Residual          10
#define     J_trivial_Ca        11

// Coupling constants.
#define	    cpl_Ca		         0
#define	    cpl_Vm		         1
#define	    cpl_IP3		         2
