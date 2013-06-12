// FIXED MODEL PARAMETERS FOR TSOUKIAS MODEL

/**** Standard model parameter values ****/
double z_K = 1.0;
double z_Na = 1.0;
double z_Ca = 2.0;
double z_Cl = -1.0;
//Avogadro's number
double N_Avogadro = 6.022e23;
//Faraday's constant (C/mol)
double F = 96487.0;
//Universal Gas constant (mJ/mol*K)
double R_const = 8.341e3;
//Temperature (K)
double T = 293.0;
//Cell membrane capacitance (pF)
double C_m = 25.0, Am = C_m / 1e6;	//cm^2
//Extracellular concentrations of Ca; Na; K and Cl (mM)
double Ca_o = 2.0, Na_o = 140.0, Cl_o = 129.0, K_o = 5.0;
//Cell volumes (pL)
//cytosolic
double vol_i = 1.0;
//available to free Ca
double vol_Ca = 0.7;
/**** L-type voltage operated Ca2+ channels ****/
//Whole cell permeability (cm/s)
double P_VOCC = 1.88e-5;
/**** Large conductance Ca2+ activated potassium channels ****/
//Single channel permeability (cm^3/s)
double P_BKCa = 3.9e-13;
//Channel density
double N_BKCa = 6.6e6;
//Fast and slow activation time constants (ms)
double tau_pf = 0.84, tau_ps = 35.9;
//Maximum NO-induced V_1/2KCa shift (mV)
double dV_half_KCaNO = 46.3;
//Maximum cGMP-induced V_1/2KCa shift (mV)
double dV_half_KCacGMP = 76.0;
/**** Voltage dependent K channels ****/
//Maximum whole cell conductance (nS)
double g_Kv = 1.35;
//Fast and slow inactivation time constants (ms)
double tau_q1 = 371.00, tau_q2 = 2884.00;
/**** ATP sensetive K channels	****/
//Whole cell conductance (nS)
double g_K_leak = 0.067;
/**** Nonselective cation channels ****/
//DAG concentration for half maximum activation (nM)
double K_NSC = 3000e-6;		//corrected and compared with jsim model
//and changed to mM
//Whole cell Na; K ; Ca permeabilities (cm/s)
double P_NaNSC = 5.11e-7, P_KNSC = 1.06 * P_NaNSC, P_CaNSC = 4.54 * P_NaNSC;
//Constituent activation
double d_NSCmin = 0.0244;
/**** Store operated nonselective cation channels ****/
//Whole cell conductances to Ca and Na (nM)
double g_SOCCa = 0.0083, g_SOCNa = 0.0575;
//Ca concentration in SR for half activation (nM)
double K_SOC = 100e-6;		//changed to mM
/**** Calcium activated chloride channels ****/
//Maximum conductance (nS/pF)
double g_ClCa = 0.23;
//cGMP independent component
double R_ClcGMPmin = 0.0132;
//Hill coefficient
double n_ClCa = 2;
//EC50 for Ca dependent activation (nM)
double K_ClCa = 365.00;
//Hill coefficient;
double n_ClcGMP = 3.3;
//EC50 for cGMP dependent activation (uM)
double K_ClcGMP = 6.4;
/**** Plasma membrane Ca ATPase ****/
//Maximum current (pA)
double I_bar_PMCA = 5.37;
//Michaelis constant (nM)
double K_mPMCA = 170e-6;		//changed to mM
/**** Plasma membrane Na-Ca exchange ****/
//Scaling factor (pS)
double g_NCX = 4.87e-4;		//0.0487;
double gama_NCX = 0.45, d_NCX = 3e-4;
/**** Sodium potasium pump ****/
//Maximum current density (pA/pF)
double I_bar_NaK = 2.3083;
//Binding constants and Hill coefficients (mM) for K_out and Na_intracellular
double K_dKo = 1.6,
		K_dNai = 22.0,
		n_HKo = 1.1,
		n_HNai = 1.7;
/**** Sodium potassium chloride cotransport ****/
//Cotransport coefficient (nmol^2/(J*s*cm^2))
double L_NaKCl = 1.79e-11;		//changed to (mmol^2/(J*s*cm^2))
/**** Sarcoplasmic reticulum ****/
/** (parameters as in Yang(2003) ) **/
//Time constant of the internal diffusion (ms)
double tau_tr = 1000.00;
//Time constant of the diffusion from SR release compartment (ms)
double tau_rel = 0.0333;
//Activation rate constant of RyR (1/ (ms * mM^2))
double Kr1 = 2500.00;
//Inactivation rate constant of RyR (1/ (ms * mM))
double Kr2 = 1.05;
//Unbinding rate constant from activation (1/ms)
double K_r1 = 0.0076;
//Unbinding rate constant from inactivation (1/ms)
double K_r2 = 0.084;
/** (Rest of the parameters as in Tsoukias(2008) ) **/
//Michaelis constant of SR Ca pump (uM)
double K_mUp = 1.0;
//Oscillation and non oscillation parameter
double k_leak = 4.0;
//Maximum SR uptake current (pA)
double I_bar_SERCA = 20.0;//(3.34)*(k_leak+1);               //ranging between 6.68 - 20 pA
//Leakage parameter
double R_leak = 5.35e-5;//1.07e-5*(k_leak);         //ranging between 1.07e-5-5.35e-5
//Binding affinity of calsequestrin (mM)
double K_CSQN = 0.8;
//Concentration of calsequestrin in the release compartment (mM)
double CSQN = 15.0;
//Rate constant of Ca release by IP3R (1/ms)
double I_bar_IP3 = 2880e-6;
//Dissociation constants for Ca activation and inhibitory sites (nM)
double K_actIP3 = 170e-6;		//changed to uM
double K_inhIP3 = 100e-6;		//changed to uM
//Dissociation constant for IP3 binding to IP3R (nM)
double K_IP3 = 120e-6;		//changed to uM
//Rate of Ca binding to the inhibitory site (1/(mM*ms))
double k_on_IP3 = 1.4;
//Volume of uptake compartment  (pl)
double vol_u = 0.07;
//Volume of release compartment (pl)
double vol_r = 0.007;
/**** alpha 1 adrenoreceptor activation and IP3 formation ****/
/** (parameters as in Bennett et al (2005) ) **/
//Total number of receptors
double R_T_G = 2.0e4;
//Unphosphorylated receptor dissociation constanst (mM)
double K_1_G = 0.01;
//Phosphorylated receptor dissociation constant (mM)
double K_2_G = 0.2;
//Receptor recycling rate (1/ms)
double k_r_G = 1.75e-7;
/*kpG declared below */
//Receptor endocytosis rate (1/ms)
double k_e_G = 6e-6;
//Fraction of mobile receptors
double ksi_G = 0.85;
//Total number of G protein molecules
double G_T_G = 1e5;
//IP3 degradation rate (1/ms)
double k_deg_G = 1.25e-3;
//G protein activation rate (1/ms)
double k_a_G = 0.17e-3;
//G protein deactivation rate (1/ms)
double k_d_G = 1.5e-3;
//Total PIP2 molecules
double PIP2_T = 5e7;
//PIP2 replenishment rate (1/ms)
double r_r_G = 0.015e-3;
//Dissociation constant for Ca binding to PLC (mM)
double K_c_G = 0.4e-3;
//Effective signal gain parameter (1/ms)
double alpha_G = 2.781e-8;
//parameter
double gama_G = N_Avogadro * vol_i * 1e-15;
/** (Rest of the parameters as in Tsoukias(2008) ) **/
//Receptor phosphorylation rate (1/s)
double k_p_G = 0.0;	//set it to 0.1(1/ms) while replicating the results in Fig3 of (Tsoukias 2008)
double deltaG_0 = 0;
/**** sGC activation and cGMP formation ****/
//sGC activation and inactivation time constants (s)
double tau_a_sGC = 0.23e3;		//changed to ms
double tau_d_sGC = 10.0e3;		//changed to ms
//Maximum cGMP formation rate (mM/ms)
double V_cGMPmax = 1.26e-7;
//Michaelis Menten Constant (mM)
double K_m_pde = 1e-3, k1_sGC = 2e3,	//Units (1/mM*ms)
		k2_sGC = 0.64e-5,		//Units (1/ms)
		k3_sGC = 4.2,		//Units (1/mM*ms)
		kD_tau_sGC = 0.1e-3,	//Units (1/ms)
		k_1_sGC = 15e-3,		//Units (1/ms)
		k_2_sGC = 0.1e-6,		//Units (1/ms)
		kD_sGC = 0.4e-3,		//Units (1/ms)
		k_4_sGC = 0.1e-3,		//Units (1/ms)
		B5_sGC = k2_sGC / k3_sGC, k_pde_cGMP = 0.0695e-3;
/**** Ca buffering in the cytosol ****/
//Concentration of calmodulin and other buffers (mM) and their dissociation constants (nM)
double S_CM = 0.1,			//mM
		B_F = 0.1,			//mM
		K_d = 260.0e-6,		//nM
		K_dB = 530.0e-6;		//nM
