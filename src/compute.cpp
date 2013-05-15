//#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"


using namespace std;
time_stamps		t_stamp;


/*******************************************************************************************/
int couplingParms(int CASE,conductance* cpl_cef)
/*******************************************************************************************/
{
	if (CASE == 1) {
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
	} else if (CASE == 2) {
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
	} else if (CASE == 3) {
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
	} else if (CASE == 4) {
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
	} else if (CASE == 5)// Simulating for experiments suggested by Dr. James Kazloski (IBM Watson Centre).
			// The homocellular Ca coupling between SMCs is changed to investigate the effects of
			//strreangth on coupling on the propagation speed of the spatial waves.
			{
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
	} else if (CASE == 6)// Simulating for experiments suggested by Dr. James Kazloski (IBM Watson Centre)
			{
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
	} else if (CASE == 7)// Simulating for experiments suggested by Dr. James Kazloski (IBM Watson Centre)
			{
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
	} else if (CASE == 8) {
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

	return 0;
}

/*******************************************************************************************/
void Initialize_koeingsberger_smc(grid_parms grid, double y[], celltype1** smc)
/*******************************************************************************************/
{
	int k = 0, offset;

	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;
			y[k + ((j - 1) * grid.neq_smc) + smc_Ca] = 0.01;
			y[k + ((j - 1) * grid.neq_smc) + smc_SR] = 0.01;
			y[k + ((j - 1) * grid.neq_smc) + smc_Vm] = 0.01;
			y[k + ((j - 1) * grid.neq_smc) + smc_w] = 0.01;
			y[k + ((j - 1) * grid.neq_smc) + smc_IP3] = 0.01;
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
void Initialize_koeingsberger_ec(grid_parms grid, double y[], celltype2** ec)
/*******************************************************************************************/
{
int	k,offset = (grid.neq_smc * grid.num_smc_circumferentially
					* grid.num_smc_axially);

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			if (i > 1)
				k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
				k = offset + 0;
			y[k + ((j - 1) * grid.neq_ec) + ec_Ca] = 0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_SR] = 0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_Vm] = 0.01;
			y[k + ((j - 1) * grid.neq_ec) + ec_IP3] = 0.01;
		}
	}
	for (int i = 0; i < (grid.num_ec_circumferentially + grid.num_ghost_cells);
			i++) {
		for (int j = 0; j < (grid.num_ec_axially + grid.num_ghost_cells);
				j++) {
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
void map_solver_to_cells(grid_parms grid,double y[], celltype1** smc, celltype2** ec)
/*******************************************************************************************/
{
	///Mapping state vectors of each cell (all except the ghost cells) to the corresponding locations on the solver's solution array y[].
	int k = 0, offset;
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;
			smc[i][j].p[smc_Ca] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca];
			smc[i][j].p[smc_SR] = y[k + ((j - 1) * grid.neq_smc) + smc_SR];
			smc[i][j].p[smc_Vm] = y[k + ((j - 1) * grid.neq_smc) + smc_Vm];
			smc[i][j].p[smc_w] = y[k + ((j - 1) * grid.neq_smc) + smc_w];
			smc[i][j].p[smc_IP3] = y[k + ((j - 1) * grid.neq_smc) + smc_IP3];
		}
	}
	offset = (grid.neq_smc * grid.num_smc_circumferentially
			* grid.num_smc_axially);

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			if (i > 1)
				k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
				k = offset + 0;
			ec[i][j].q[ec_Ca] = y[k + ((j - 1) * grid.neq_ec) + ec_Ca];
			ec[i][j].q[ec_SR] = y[k + ((j - 1) * grid.neq_ec) + ec_SR];
			ec[i][j].q[ec_Vm] = y[k + ((j - 1) * grid.neq_ec) + ec_Vm];
			ec[i][j].q[ec_IP3] = y[k + ((j - 1) * grid.neq_ec) + ec_IP3];
		}
	}
}
/*******************************************************************************************/
void map_GhostCells_to_cells(celltype1** smc, celltype2** ec, grid_parms grid)
/*******************************************************************************************/
{
	///Allocating arrays of appropriate lengths for holding values of ghost cell variable to which the relevant smc[i][j] and ec[i][j] members will point to.
		///The rest of the computational domain is to point to the solution vector y of the solver.
		double **ghost_cells_smc_circumferentially,
				**ghost_cells_ec_circumferentially,
				**ghost_cells_smc_axially,
				**ghost_cells_ec_axially;

		ghost_cells_smc_circumferentially = (double**) checked_malloc(
				sizeof(double*), "ghost cell in circ direction for smc");
		for (int i=0; i<2; i++){
			ghost_cells_smc_circumferentially[i] = (double*) checked_malloc(
				(grid.num_ghost_cells + grid.num_smc_circumferentially)
						* grid.neq_smc * sizeof(double),
						"ghost cell in circ direction for smc");
		}

		ghost_cells_ec_circumferentially = (double**) checked_malloc(
					sizeof(double*), "ghost cell in circ direction for ec");
		for (int i=0; i<2; i++){
				ghost_cells_ec_circumferentially[i] = (double*) checked_malloc(
					(grid.num_ghost_cells + grid.num_ec_circumferentially)
							* grid.neq_ec * sizeof(double),
							"ghost cell in circ direction for smc");
			}

		ghost_cells_smc_axially = (double**) checked_malloc(
					sizeof(double*), "ghost cell in axial direction for smc");
		for (int i=0; i<2; i++){
			ghost_cells_smc_axially[i] = (double*) checked_malloc(
				(grid.num_ghost_cells + grid.num_smc_axially)
						* grid.neq_smc * sizeof(double),
						"ghost cell in axial direction for smc");
		}

		ghost_cells_ec_axially = (double**) checked_malloc(
					sizeof(double*),  "ghost cell in axial direction for ec");
		for (int i=0; i<2; i++){
				ghost_cells_ec_axially[i] = (double*) checked_malloc(
					(grid.num_ghost_cells + grid.num_ec_axially)
							* grid.neq_ec * sizeof(double),
							"ghost cell in axial direction for smc");
			}
}

/*******************************************************************************************/
void single_cell(double t, double y[], grid_parms grid,celltype1** smc, celltype2** ec)
/*******************************************************************************************/
{
	/* Constants for homogenically coupled SMCs*/
	const double Fi = 0.23;
	const double Kri = 1.00;
	const double GCai = 0.00129;
	const double vCa1 = 100.00;
	const double vCa2 = -24.00;
	const double RCai = 8.50;
	const double GNaCai = 0.00316;
	const double cNaCai = 0.5;
	const double vNaCai = -30.00;
	const double Bi = 2.025;
	const double cbi = 1.0;
	const double CICRi = 55.00;
	const double sci = 2.0;
	const double cci = 0.9;
	const double Di = 0.24;
	const double vdi = -100.00;
	const double Rdi = 250.00;
	const double Li = 0.025;
	const double gama = 1970.00;
	const double FNaK = 0.0432;
	const double GCli = 0.00134;
	const double vCl = -25.0;
	const double GKi = 0.0046;
	const double vKi = -94.00;
	const double lambda = 45.00;
	const double cwi = 0.0;
	const double beta = 0.13;
	const double vCa3 = -27.0;
	const double RKi = 12.0;
	const double ki = 0.1;
	/* Constants for homogenically coupled ECs*/
	const double Fj = 0.23;
	const double Krj = 1.00;
	const double Bj = 0.5;
	const double cbj = 1.0;
	const double CICRj = 5.0;
	const double scj = 2.0;
	const double ccj = 0.9;
	const double Dj = 0.24;
	const double Lj = 0.025;
	const double kj = 0.1;
	const double Gcatj = 0.66 * 1e-3;
	const double ECa = 50.00;
	const double m3cat = -0.18;			//-6.18;
	const double m4cat = 0.37;
	const double J0j = 0.029;
	const double Cmj = 25.8;
	const double Gtot = 6927;
	const double vKj = -80.0;
	const double a1j = 53.3;
	const double a2j = 53.3;
	const double bj = -80.8;
	const double c1j = -0.4; //-6.4;
	const double m3b = 1.32e-3;
	const double m4b = 0.30;
	const double m3s = -0.28;
	const double m4s = 0.389;
	const double GRj = 955;
	const double vrestj = -31.10;

	/*Intracellular calcium buffering*/
	const double k6 = 100.00;
	const double k7 = 300.00;
	const double BT = 120.00;

	int i, j, k;
//#pragma omp for private(j)
//EVALUATING SINGLE CELL FLUXES :::::::For SMC::::::::
	for (i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (j = 1; j <= grid.num_smc_axially; j++) {

			//JIP3
			smc[i][j].A[J_IP3] = (Fi * P2(smc[i][j].p[smc_IP3]))
					/ (P2(Kri) + P2(smc[i][j].p[smc_IP3]));
			//JSRuptake
			smc[i][j].A[J_SERCA] = (Bi * P2(smc[i][j].p[smc_Ca]))
					/ (P2(smc[i][j].p[smc_Ca]) + P2(cbi));
			//Jcicr
			smc[i][j].A[J_CICR] =
					(CICRi
							* (P2(smc[i][j].p[smc_SR])
									* P4(smc[i][j].p[smc_Ca])))
							/ ((P2(sci) + P2(smc[i][j].p[smc_SR]))* ( P4(cci) + P4(smc[i][j].p[smc_Ca]) ) );
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
					/ (P2(smc[i][j].p[smc_Ca] + cwi)
							+ (beta
									* ((double) exp(
											(-1) * (smc[i][j].p[smc_Vm] - vCa3)
													/ RKi))));
			//Jdegradation
			smc[i][j].A[J_IP3_deg] = ki * smc[i][j].p[smc_IP3];

		} //end for j
	} //end for i

//#pragma omp for private(j)
//EVALUATING SINGLE CELL FLUXES For EC
	for (i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (j = 1; j <= grid.num_ec_axially; j++) {
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

}	//End of function single_cell

/*******************************************************************************************/
void coupling(double t,double y[], grid_parms grid,celltype1** smc, celltype2** ec,conductance cpl_cef)
/*******************************************************************************************/
{

	int i, j, k, l;

////******************** HOMOCELLULAR COUPLING *********************/
	for (i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (j = 1; j <= grid.num_smc_axially; j++) {
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
			smc[i][j].B[cpl_Ca] = -cpl_cef.Ca_hm_smc
					* ((smc[i][j].p[smc_Ca] - smc[i][up].p[smc_Ca])
							+ (smc[i][j].p[smc_Ca] - smc[i][down].p[smc_Ca])
							+ (smc[i][j].p[smc_Ca] - smc[left][j].p[smc_Ca])
							+ (smc[i][j].p[smc_Ca] - smc[right][j].p[smc_Ca]));
			smc[i][j].B[cpl_Vm] =
			-cpl_cef.Vm_hm_smc
			* ((smc[i][j].p[smc_Vm] - smc[i][up].p[smc_Vm])
					+ (smc[i][j].p[smc_Vm] - smc[i][down].p[smc_Vm])
					+ (smc[i][j].p[smc_Vm] - smc[left][j].p[smc_Vm])
					+ (smc[i][j].p[smc_Vm] - smc[right][j].p[smc_Vm]));
			smc[i][j].B[cpl_IP3] =
					-cpl_cef.IP3_hm_smc
							* ((smc[i][j].p[smc_IP3] - smc[i][up].p[smc_IP3])
									+ (smc[i][j].p[smc_IP3]
											- smc[i][down].p[smc_IP3])
									+ (smc[i][j].p[smc_IP3]
											- smc[left][j].p[smc_IP3])
									+ (smc[i][j].p[smc_IP3]
											- smc[right][j].p[smc_IP3]));
		}	//end j
	}	//end i

	for (i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (j = 1; j <= grid.num_ec_axially; j++) {
			int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
			ec[i][j].B[cpl_Ca] = -cpl_cef.Ca_hm_ec
					* ((ec[i][j].q[ec_Ca] - ec[i][up].q[ec_Ca])
							+ (ec[i][j].q[ec_Ca] - ec[i][down].q[ec_Ca])
							+ (ec[i][j].q[ec_Ca] - ec[left][j].q[ec_Ca])
							+ (ec[i][j].q[ec_Ca] - ec[right][j].q[ec_Ca]));
			ec[i][j].B[cpl_Vm] =
			-cpl_cef.Vm_hm_ec
			* ((ec[i][j].q[ec_Vm] - ec[i][up].q[ec_Vm])
					+ (ec[i][j].q[ec_Vm] - ec[i][down].q[ec_Vm])
					+ (ec[i][j].q[ec_Vm] - ec[left][j].q[ec_Vm])
					+ (ec[i][j].q[ec_Vm] - ec[right][j].q[ec_Vm]));
			ec[i][j].B[cpl_IP3] = -cpl_cef.IP3_hm_ec
					* ((ec[i][j].q[ec_IP3] - ec[i][up].q[ec_IP3])
							+ (ec[i][j].q[ec_IP3] - ec[i][down].q[ec_IP3])
							+ (ec[i][j].q[ec_IP3] - ec[left][j].q[ec_IP3])
							+ (ec[i][j].q[ec_IP3] - ec[right][j].q[ec_IP3]));

		}	//end j
	}	//end i

////******************** HETROCELLULAR COUPLING *********************/
	int offset_smc_circumferentially, offset_ec_axially;

	i = 0;
	j = 0;
	k = 0;
	l = 0;

	for (i = 1; i <= grid.num_smc_circumferentially; i++) {
		l = 1;
		for (j = 1; j <= grid.num_smc_axially; j++) {
			double dummy_smc[3] = { 0.0, 0.0, 0.0 };
			for (k = 1 + (i - 1) * 5; k <= i * 5; k++) {
				dummy_smc[cpl_Ca] = dummy_smc[cpl_Ca]
						+ (smc[i][j].p[smc_Ca] - ec[k][l].q[ec_Ca]);
				dummy_smc[cpl_Vm] = dummy_smc[cpl_Vm]
						+ (smc[i][j].p[smc_Vm] - ec[k][l].q[ec_Vm]);
				dummy_smc[cpl_IP3] = dummy_smc[cpl_IP3]
						+ (smc[i][j].p[smc_IP3] - ec[k][l].q[ec_IP3]);
			}
			if ((j % grid.num_smc_fundblk_axially) == 0) {
				l++;
			}
			smc[i][j].C[cpl_Ca] = -cpl_cef.Ca_ht_smc * dummy_smc[cpl_Ca];
			smc[i][j].C[cpl_Vm] = -cpl_cef.Vm_ht_smc * dummy_smc[cpl_Vm];
			smc[i][j].C[cpl_IP3] = -cpl_cef.IP3_ht_smc * dummy_smc[cpl_IP3];
		}
	}
	i = 0;
	j = 0;
	k = 0;
	l = 0;
	for (i = 1; i <= grid.num_ec_circumferentially; i++) {
		if ((i - 1) % 5 == 0)
			k++;
		for (j = 1; j <= grid.num_ec_axially; j++) {
			double dummy_ec[3] = { 0.0, 0.0, 0.0 };

			for (l = 1 + (j - 1) * 13; l <= j * 13; l++) {
				dummy_ec[cpl_Ca] = dummy_ec[cpl_Ca]
						+ (ec[i][j].q[ec_Ca] - smc[k][l].p[smc_Ca]);
				dummy_ec[cpl_Vm] = dummy_ec[cpl_Vm]
						+ (ec[i][j].q[ec_Vm] - smc[k][l].p[smc_Vm]);
				dummy_ec[cpl_IP3] = dummy_ec[cpl_IP3]
						+ (ec[i][j].q[ec_IP3] - smc[k][l].p[smc_IP3]);
			}
			ec[i][j].C[cpl_Ca] = -cpl_cef.Ca_ht_ec * dummy_ec[cpl_Ca];
			ec[i][j].C[cpl_Vm] = -cpl_cef.Vm_ht_ec * dummy_ec[cpl_Vm];
			ec[i][j].C[cpl_IP3] = -cpl_cef.IP3_ht_ec * dummy_ec[cpl_IP3];

		}
	}

}	//end of coupling()

/*******************************************************************************************/
double agonist_profile(double t, grid_parms grid, int i, int j){

	double JPLC;
	if (t > grid.stimulus_onset_time) {

				JPLC =grid.min_jplc+ (grid.max_jplc/
						(1 + exp(-grid.gradient * ( ((j-1)+grid.num_ec_axially*floor(grid.rank/grid.n)) -(grid.m*grid.num_ec_axially / 2) )) ) );

			} else if (t <= grid.stimulus_onset_time){
				JPLC = grid.uniform_jplc;
			}
	return JPLC;
}


/*******************************************************************************************/
//double z_coordinate(grid_parms gird, double domains[][], int num_subdomains, int i, int j)
/*******************************************************************************************/



/*******************************************************************************************/
void initialize_t_stamp(time_stamps t_stamp){
	t_stamp.diff_async_comm_calls	=	0.0;
	t_stamp.diff_async_comm_calls_wait=	0.0;
	t_stamp.diff_barrier_in_solver_before_comm=0.0;
	t_stamp.diff_map_function =0.0;
	t_stamp.diff_single_cell_fluxes=0.0;
	t_stamp.diff_coupling_fluxes=0.0;
}
