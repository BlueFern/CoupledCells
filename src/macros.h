
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define SRC_LOC __FILE__ ":" TOSTRING(__LINE__)

#define CHECK_MPI_ERROR(fn) \
	{ \
	int errcode = (fn); \
	if(errcode != MPI_SUCCESS) \
	{ \
		fprintf(stderr, "MPI ERROR: %d; %s.\n", errcode, SRC_LOC); \
		MPI_Abort(MPI_COMM_WORLD, 911); \
	} \
	}

/****** marcos for identifying models ******/
#define 		KNBGR			0
#define 		TSK				1

// WARNING, THESE MACROS WERE/ARE INCONSISTENT WITH KOENIGSBERGER_MACROS.H
// TODO: Remove duplicate macros and put the rest in tsoukias.
/*******Koenigsberger model*******/
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
#define     J_Cl	             9
#define     J_K                 10
#define     K_activation        11
//EC specific currents
#define     J_NSC               6
#define     J_BK_Ca             7
#define     J_SK_Ca             8
#define     J_Ktot              9
#define     J_Residual          10
#define     J_trivial_Ca        11



#define	    cpl_Ca		         0
#define	    cpl_Vm		         1
#define	    cpl_IP3		         2

/*******Tsoukias model*******/
//SMC state variables in Tsoukias model
#define		smc_q_1		    	smc_w		/// Kv fast inact q1 gate probability
#define		smc_q_2		    	5			/// Kv slow inact q2 gate probability
#define 	smc_d_L 	    	6			/// L-type Ca d gate probability
#define		smc_f_L         	7			/// L-type Ca f gate probability
#define		smc_p_f     		8			/// BKCa fast p gate probability
#define		smc_p_s         	9			/// BKCa slow p gate probability
#define		smc_p_K		    	10			/// Kv activating p gate probability
#define		smc_h_IP3	    	11			/// IIP3 recept h gate probability
#define		smc_Ca_u			smc_SR		/// SR uptake Ca concentration
#define		smc_Ca_r	    	12			/// SR release Ca concentration
#define		smc_R_10	    	13			/// RyR w/ act Ca bnd state
#define		smc_R_11	    	14			/// RyR w/act-inact Ca bnd state
#define		smc_R_01	    	15			/// RyR w/ inact Ca bnd state
#define		smc_V_cGMP  		16			/// Rate of cGMP formation
#define		smc_cGMP_i			17			/// Cytosolic cGMP concentration
#define		smc_Na_i			18			/// Cytosolic Na+ concentration
#define		smc_K_i				19			/// Cytosolic K+ concentration
#define		smc_Cl_i			20			/// Cytosolic Cl- concentration
#define 	smc_DAG				21			/// Cytosolic DAG concentration
#define 	smc_PIP2        	22			/// Membrane bound PIP2 concentration
#define 	smc_R_S_G       	23
#define 	smc_R_S_P_G     	24
#define 	smc_G           	25

#define		E_K					0	// K+ Nernst potential
#define 	E_Na				1	// Na+ Nernst potential
#define 	E_Ca				2	// Ca2+ Nernst potential
#define 	E_Cl				3	// Cl- Nernst potential
#define		dbar_L	   			4		// L-type Ca SS d gate prob
#define 	fbar_L	    		5		// L-type Ca SS f gate prob
#define 	tau_dL				6	// L-type Ca d gate time const
#define 	tau_fL				7	// L-type Ca f gate time const
#define 	I_VOCC				8	// L-type Ca total channel current
#define		P_KCa				9	// BKCa total fast/slow gate prob
#define		i_KCa				10	// BKCa single channel current
#define		R_NO				11	// BKCa NO V50KCa shift fact
#define 	R_cGMP	   			12		// BKCa cGMP V50KCa shift fact
#define		V_half_KCa			13	// BKCa 1/2 max act membr potential
#define		pbar_o	    		14		// BKCa SS fast/slow gate prob
#define		I_BKCa				15	// BKCa total channel current
#define		pbar_K	    		16		// Kv SS act p gate probability
#define		qbar				17	// Kv SS inact q1/q2 gate probability
#define		tau_pK				18	// Kv act p gate time constant
#define		I_Kv				19	// Kv total channel current
#define		I_Kleak				20	// Kleak total channel current
#define		Po_NSC	    		21		// NSC act gate probability
#define		INa_NSC				22	// NSC total Na+ current
#define		IK_NSC				23	// NSC total K+ current
#define		ICa_NSC				24	// NSC total Ca2+ current
#define		I_NSC				25	// NSC total channel current
#define		P_SOC	    		26		// SOC act gate probability
#define 	INa_SOC				27	// SOC total Na+ current
#define 	ICa_SOC				28	// SOC total Ca2+ current
#define 	I_SOC				29	// SOC total channel current
#define		alpha_Cl			30	// ClCa cGMP-dep factor
#define 	K_ClCacGMP			31	// ClCa cGMP 1/2 max act concern
#define 	P_Cl				32	// ClCa act gate probability
#define 	I_ClCa				33	// ClCa total channel current
#define		I_PMCA				34	// PMCA total pump current
#define		R_NCX_cGMP			35	// NCX cGMP-dep act factor
#define 	phi_F				36	// NCX V-dep forward flux partition
#define 	phi_R				37	// NCX V-dep reverse flux partition
#define 	I_NCX				38	// NCX total exchanger current
#define		I_NaK				39	// NaK total pump current
#define		R_NaKCl_cGMP		40	// NaKCl cGMP-dep act factor
#define 	I_Cl_NaKCl			41	// NaKCl total Cl- current
#define 	I_Na_NaKCl			42	// NaKCl total Na+ current
#define 	I_K_NaKCl			43	// NaKCl total K+ current
#define		I_IP3				44	// IP3 receptor total current
#define		I_SERCA				45	// SERCA total current
#define 	I_tr				46	// Diff current btwn SR_u and SR_r
#define 	I_rel				47	// RyR total current
#define		R_00				48	// RyR state with no Ca2+ bound
#define		V_bar_cGMP			49	// cGMP max rate of formation.
#define 	A0_sGC      		50
#define		A1_sGC     			51
#define		tau_sGC				52	// cGMP sGC-dep time constant
#define		I_Catotm			53	// Total Ca2+ membrane current
#define 	I_Natotm			54	// Total Na+ membrane current
#define 	I_Ktotm				55	// Total K+ membrane current
#define 	I_Cltotm			56	// Total Cl- membrane current
#define 	Q           		57
#define 	Q_10        		58
#define 	r_h_G       		59
#define 	rho_r_G     		60
#define 	delta_G     		61


