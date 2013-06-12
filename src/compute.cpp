//#include <omp.h>
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

/*************************************************************************/
int map_solver_to_cells(grid_parms grid, double* y, int smc_model,
		celltype1** smc, int ec_model, celltype2** ec) {
/*************************************************************************/
	int err = 1;
	switch (smc_model) {
		case (TSK): {
			int k = 0, offset;
			for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
				for (int j = 1; j <= grid.num_smc_axially; j++) {
					if (i > 1)
					k = ((i - 1) * grid.neq_smc_axially);
					else if (i == 1)
					k = 0;
					smc[i][j].p[smc_Vm] = y[k + ((j - 1) * grid.neq_smc) + smc_Vm];
					smc[i][j].p[smc_d_L] = y[k + ((j - 1) * grid.neq_smc) + smc_d_L];
					smc[i][j].p[smc_f_L] = y[k + ((j - 1) * grid.neq_smc) + smc_f_L];
					smc[i][j].p[smc_p_f] = y[k + ((j - 1) * grid.neq_smc) + smc_p_f];
					smc[i][j].p[smc_p_s] = y[k + ((j - 1) * grid.neq_smc) + smc_p_s];
					smc[i][j].p[smc_q_1] = y[k + ((j - 1) * grid.neq_smc) + smc_q_1];
					smc[i][j].p[smc_q_2] = y[k + ((j - 1) * grid.neq_smc) + smc_q_2];
					smc[i][j].p[smc_p_K] = y[k + ((j - 1) * grid.neq_smc) + smc_p_K];
					smc[i][j].p[smc_Ca_u] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca_u];
					smc[i][j].p[smc_Ca_r] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca_r];
					smc[i][j].p[smc_R_10] = y[k + ((j - 1) * grid.neq_smc) + smc_R_10];
					smc[i][j].p[smc_R_11] = y[k + ((j - 1) * grid.neq_smc) + smc_R_11];
					smc[i][j].p[smc_R_01] = y[k + ((j - 1) * grid.neq_smc) + smc_R_01];
					smc[i][j].p[smc_h_IP3] = y[k + ((j - 1) * grid.neq_smc) + smc_h_IP3];
					smc[i][j].p[smc_R_S_G] = y[k + ((j - 1) * grid.neq_smc) + smc_R_S_G];
					smc[i][j].p[smc_R_S_P_G] = y[k + ((j - 1) * grid.neq_smc) + smc_R_S_P_G];
					smc[i][j].p[smc_G] = y[k + ((j - 1) * grid.neq_smc) + smc_G];
					smc[i][j].p[smc_IP3] = y[k + ((j - 1) * grid.neq_smc) + smc_IP3];
					smc[i][j].p[smc_PIP2] = y[k + ((j - 1) * grid.neq_smc) + smc_PIP2];
					smc[i][j].p[smc_V_cGMP] = y[k + ((j - 1) * grid.neq_smc) + smc_V_cGMP];
					smc[i][j].p[smc_cGMP_i] = y[k + ((j - 1) * grid.neq_smc) + smc_cGMP_i];
					smc[i][j].p[smc_Ca] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca];
					smc[i][j].p[smc_Na_i] = y[k + ((j - 1) * grid.neq_smc) + smc_Na_i];
					smc[i][j].p[smc_K_i] = y[k + ((j - 1) * grid.neq_smc) + smc_K_i];
					smc[i][j].p[smc_Cl_i] = y[k + ((j - 1) * grid.neq_smc) + smc_Cl_i];
					smc[i][j].p[smc_DAG] = y[k + ((j - 1) * grid.neq_smc) + smc_DAG];
				}
			}
			break;
		}case (KNBGR): {
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
			break;
		}default: {
			err =1;
			break;
		}
	}
	switch (ec_model) {
		case (TSK): {
			int k, offset = (grid.neq_smc * grid.num_smc_circumferentially
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
		}case (KNBGR): {
			int k, offset = (grid.neq_smc * grid.num_smc_circumferentially
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
		}default: {
			err =1;
			break;
		}
	}
return (err);
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
	return (JPLC);
}
/*******************************************************************************************/
void initialize_t_stamp(time_stamps t_stamp){
/*******************************************************************************************/
	t_stamp.diff_async_comm_calls	=	0.0;
	t_stamp.diff_async_comm_calls_wait=	0.0;
	t_stamp.diff_barrier_in_solver_before_comm=0.0;
	t_stamp.diff_map_function =0.0;
	t_stamp.diff_single_cell_fluxes=0.0;
	t_stamp.diff_coupling_fluxes=0.0;
}

/*****************************************************************************/
int compute(time_stamps t_stamp, grid_parms grid, int smc_model,
		celltype1** smc, int ec_model, celltype2** ec, conductance cpl_cef, double t, double* y,
		double* f, int NO_Path, int cGMP_Path) {
/*****************************************************************************/

	int err;
	t_stamp.computeDerivatives_call_counter =
			t_stamp.computeDerivatives_call_counter + 1;

	t_stamp.map_function_t1 = MPI_Wtime();
	map_solver_to_cells(grid, y, smc_model, smc, ec_model, ec);
	t_stamp.map_function_t2 = MPI_Wtime();
	t_stamp.diff_map_function = t_stamp.diff_map_function
			+ (t_stamp.map_function_t2 - t_stamp.map_function_t1);

	t_stamp.single_cell_fluxes_t1 = MPI_Wtime();
	switch (smc_model) {
	case (TSK): {
		tsoukias_smc(grid, smc, NO_Path, cGMP_Path);
		break;
	}
	case (KNBGR): {
		koenigsberger_smc(grid, smc);
		break;
	}
	default: {
		err = 1;
		break;
	}
	}
	switch (ec_model) {
	case (TSK): {
		koenigsberger_ec(grid, ec);
		break;
	}
	case (KNBGR): {
		koenigsberger_ec(grid, ec);
		break;
	}
	default: {
		err = 1;
		break;
	}
	}
	t_stamp.single_cell_fluxes_t2 = MPI_Wtime();
	t_stamp.diff_single_cell_fluxes = t_stamp.diff_single_cell_fluxes
			+ (t_stamp.single_cell_fluxes_t2 - t_stamp.single_cell_fluxes_t1);

	t_stamp.coupling_fluxes_t1 = MPI_Wtime();
	coupling(t, y, grid, smc, ec, cpl_cef);
	t_stamp.coupling_fluxes_t2 = MPI_Wtime();
	t_stamp.diff_coupling_fluxes = t_stamp.diff_coupling_fluxes
			+ (t_stamp.coupling_fluxes_t2 - t_stamp.coupling_fluxes_t1);

	tsoukias_smc_derivatives(f, grid, smc);
	koenigsberger_ec_derivatives(t, f, grid, ec);

	return (err);
}
