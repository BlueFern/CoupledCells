//#include <omp.h>
#include "computelib.h"
using namespace std;
time_stamps t_stamp;

/**
 * Wrapper around malloc to catch failed memory allocation. If allocation fails
 * MPI_Abort is called.
 *
 * \param bytes is the size of requested memory.
 * \param errmsg is the message produced in the event of failed memory allocation.
 *
 * \todo Perhaps the error messages should be printed to stderr.
 */
void* checked_malloc(size_t bytes, const char* errmsg) {
	void *pval = malloc(bytes);

	if (pval == NULL) {
		fprintf(stdout, "%s", errmsg);
		MPI_Abort(MPI_COMM_WORLD, 100);
	}

	return pval;
}

/*******************************************************************************************/
int couplingParms(int CASE, conductance* cpl_cef)
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
	} else if (CASE == 5) // Simulating for experiments suggested by Dr. James Kazloski (IBM Watson Centre).
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
int map_solver_to_cells(grid_parms grid, double* y, celltype1** smc,
		celltype2** ec) {
	/*************************************************************************/
	int err = 0;
	switch (grid.smc_model) {
	case (TSK): {
		int k = 0, offset;
		for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
			for (int j = 1; j <= grid.num_smc_axially; j++) {
				if (i > 1)
					k = ((i - 1) * grid.neq_smc_axially);
				else if (i == 1)
					k = 0;
				smc[i][j].p[smc_Vm] = y[k + ((j - 1) * grid.neq_smc) + smc_Vm];
				smc[i][j].p[smc_d_L] =
						y[k + ((j - 1) * grid.neq_smc) + smc_d_L];
				smc[i][j].p[smc_f_L] =
						y[k + ((j - 1) * grid.neq_smc) + smc_f_L];
				smc[i][j].p[smc_p_f] =
						y[k + ((j - 1) * grid.neq_smc) + smc_p_f];
				smc[i][j].p[smc_p_s] =
						y[k + ((j - 1) * grid.neq_smc) + smc_p_s];
				smc[i][j].p[smc_q_1] =
						y[k + ((j - 1) * grid.neq_smc) + smc_q_1];
				smc[i][j].p[smc_q_2] =
						y[k + ((j - 1) * grid.neq_smc) + smc_q_2];
				smc[i][j].p[smc_p_K] =
						y[k + ((j - 1) * grid.neq_smc) + smc_p_K];
				smc[i][j].p[smc_Ca_u] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_Ca_u];
				smc[i][j].p[smc_Ca_r] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_Ca_r];
				smc[i][j].p[smc_R_10] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_R_10];
				smc[i][j].p[smc_R_11] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_R_11];
				smc[i][j].p[smc_R_01] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_R_01];
				smc[i][j].p[smc_h_IP3] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_h_IP3];
				smc[i][j].p[smc_R_S_G] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_R_S_G];
				smc[i][j].p[smc_R_S_P_G] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_R_S_P_G];
				smc[i][j].p[smc_G] = y[k + ((j - 1) * grid.neq_smc) + smc_G];
				smc[i][j].p[smc_IP3] =
						y[k + ((j - 1) * grid.neq_smc) + smc_IP3];
				smc[i][j].p[smc_PIP2] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_PIP2];
				smc[i][j].p[smc_V_cGMP] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_V_cGMP];
				smc[i][j].p[smc_cGMP_i] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_cGMP_i];
				smc[i][j].p[smc_Ca] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca];
				smc[i][j].p[smc_Na_i] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_Na_i];
				smc[i][j].p[smc_K_i] =
						y[k + ((j - 1) * grid.neq_smc) + smc_K_i];
				smc[i][j].p[smc_Cl_i] = y[k + ((j - 1) * grid.neq_smc)
						+ smc_Cl_i];
				smc[i][j].p[smc_DAG] =
						y[k + ((j - 1) * grid.neq_smc) + smc_DAG];
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
				smc[i][j].p[smc_Ca] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca];
				smc[i][j].p[smc_SR] = y[k + ((j - 1) * grid.neq_smc) + smc_SR];
				smc[i][j].p[smc_Vm] = y[k + ((j - 1) * grid.neq_smc) + smc_Vm];
				smc[i][j].p[smc_w] = y[k + ((j - 1) * grid.neq_smc) + smc_w];
				smc[i][j].p[smc_IP3] =y[k + ((j - 1) * grid.neq_smc) + smc_IP3];

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
				ec[i][j].q[ec_Ca] = y[k + ((j - 1) * grid.neq_ec) + ec_Ca];
				ec[i][j].q[ec_SR] = y[k + ((j - 1) * grid.neq_ec) + ec_SR];
				ec[i][j].q[ec_Vm] = y[k + ((j - 1) * grid.neq_ec) + ec_Vm];
				ec[i][j].q[ec_IP3] = y[k + ((j - 1) * grid.neq_ec) + ec_IP3];
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
				ec[i][j].q[ec_Ca] = y[k + ((j - 1) * grid.neq_ec) + ec_Ca];
				ec[i][j].q[ec_SR] = y[k + ((j - 1) * grid.neq_ec) + ec_SR];
				ec[i][j].q[ec_Vm] = y[k + ((j - 1) * grid.neq_ec) + ec_Vm];
				ec[i][j].q[ec_IP3] = y[k + ((j - 1) * grid.neq_ec) + ec_IP3];
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

/*******************************************************************************************/
void map_GhostCells_to_cells(celltype1** smc, celltype2** ec, grid_parms grid)
/*******************************************************************************************/
{
	///Allocating arrays of appropriate lengths for holding values of ghost cell variable to which the relevant smc[i][j] and ec[i][j] members will point to.
	///The rest of the computational domain is to point to the solution vector y of the solver.
	double **ghost_cells_smc_circumferentially,
			**ghost_cells_ec_circumferentially, **ghost_cells_smc_axially,
			**ghost_cells_ec_axially;

	ghost_cells_smc_circumferentially = (double**) checked_malloc(
			sizeof(double*), "ghost cell in circ direction for smc");
	for (int i = 0; i < 2; i++) {
		ghost_cells_smc_circumferentially[i] = (double*) checked_malloc(
				(grid.num_ghost_cells + grid.num_smc_circumferentially)
						* grid.neq_smc * sizeof(double),
				"ghost cell in circ direction for smc");
	}

	ghost_cells_ec_circumferentially = (double**) checked_malloc(
			sizeof(double*), "ghost cell in circ direction for ec");
	for (int i = 0; i < 2; i++) {
		ghost_cells_ec_circumferentially[i] = (double*) checked_malloc(
				(grid.num_ghost_cells + grid.num_ec_circumferentially)
						* grid.neq_ec * sizeof(double),
				"ghost cell in circ direction for smc");
	}

	ghost_cells_smc_axially = (double**) checked_malloc(sizeof(double*),
			"ghost cell in axial direction for smc");
	for (int i = 0; i < 2; i++) {
		ghost_cells_smc_axially[i] = (double*) checked_malloc(
				(grid.num_ghost_cells + grid.num_smc_axially) * grid.neq_smc
						* sizeof(double),
				"ghost cell in axial direction for smc");
	}

	ghost_cells_ec_axially = (double**) checked_malloc(sizeof(double*),
			"ghost cell in axial direction for ec");
	for (int i = 0; i < 2; i++) {
		ghost_cells_ec_axially[i] = (double*) checked_malloc(
				(grid.num_ghost_cells + grid.num_ec_axially) * grid.neq_ec
						* sizeof(double),
				"ghost cell in axial direction for smc");
	}
}

/*******************************************************************************************/
void coupling(double t, double y[], grid_parms grid, celltype1** smc,
		celltype2** ec, conductance cpl_cef)
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
			smc[i][j].B[cpl_Vm] = -cpl_cef.Vm_hm_smc
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
			ec[i][j].B[cpl_Vm] = -cpl_cef.Vm_hm_ec
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


/**
 * Provide a JPLC value for a given position along the axial dimension of the
 * current domain. The value is provided by a sigmoid function.
 *
 * \todo Provide a formula for this equation.
 *
 * \param t
 * \param grid
 * \param i
 * \param j
 * \param axial_coordinate
 * \return
 */
double agonist_profile(double t, grid_parms grid, int i, int j, double axial_coordinate)
{
	double JPLC;
	if (t > grid.stimulus_onset_time) {

		/*	JPLC =grid.min_jplc+ (grid.max_jplc/
		 (1 + exp(-grid.gradient * ( ((j-1)+grid.num_ec_axially*floor(grid.rank/grid.n)) -(grid.m*grid.num_ec_axially / 2) )) ) );
		 */
		JPLC = grid.min_jplc
				+ (grid.max_jplc / (1.0 + exp(-grid.gradient * axial_coordinate)));
	} else if (t <= grid.stimulus_onset_time) {
		JPLC = grid.uniform_jplc;
	}
	return JPLC;
}

/**********************************************************************/
celltype2** ith_ec_z_coordinate(grid_parms grid, celltype2** ec)
/**********************************************************************/
{
	/*	double array[2 * grid.num_ec_axially];

	 for (int i = 0; i <= 2 * grid.num_ec_axially; i++) {
	 array[i] = grid.my_domain.local_z_end
	 + (i)
	 * ((grid.my_domain.local_z_start
	 - grid.my_domain.local_z_end)
	 / (2 * grid.num_ec_axially));
	 }

	 for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
	 int indx = 2 * grid.num_ec_axially - 1;
	 for (int j = 1; j <= grid.num_ec_axially; j++) {
	 ec[i][j].z_coord = array[indx];
	 indx -= 2;
	 }
	 }
	 return (ec);*/
	double adapted_ec_length = grid.hx_ec;//-(grid.my_domain.local_z_end	- grid.my_domain.local_z_start) / grid.num_ec_axially;
	double array[grid.num_ec_axially];
	for (int i = 0; i < grid.num_ec_axially; i++) {
		array[i] = grid.my_domain.local_z_end
				+ ((double) (i) * adapted_ec_length) + (adapted_ec_length / 2);
	}

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			ec[i][j].z_coord = array[grid.num_ec_axially - i];
		}
	}

	return (ec);
}
/**********************************************************************/
void initialize_t_stamp(time_stamps* t_stamp) {
	t_stamp->diff_async_comm_calls = 0.0;
	t_stamp->diff_async_comm_calls_wait = 0.0;
	t_stamp->diff_barrier_in_solver_before_comm = 0.0;
	t_stamp->diff_map_function = 0.0;
	t_stamp->diff_single_cell_fluxes = 0.0;
	t_stamp->diff_coupling_fluxes = 0.0;
}
/********************************************************************************/
int recognize_end_of_file_index(checkpoint_handle* check, grid_parms grid) {
	/*******************************************************************************/
	MPI_Offset disp;
	MPI_Status status;
	int index;

	disp = grid.universal_rank * sizeof(int);
	check_flag(
			MPI_File_read_at(check->line_number, disp, &index, 1, MPI_INT,
					&status),
			"error reading the line number in recognize_end_of_file_index.");
	return (index);
}
/********************************************************************************/
double reinitialize_time(checkpoint_handle* check, int line_index,
		grid_parms grid) {
	/************************************************************************/
	MPI_Offset disp;
	MPI_Status status;

	int elements = 1;
	double time;

	disp = ((line_index - 1) * grid.tasks * elements * sizeof(double))
			+ (grid.rank * elements * sizeof(double));
	check_flag(
			MPI_File_read_at(check->Time, disp, &time, elements, MPI_DOUBLE,
					&status), "error read in the reinit data in reinit_smc.");
	return (time);
}
/********************************************************************************/
double* reinitialize_koenigsberger_smc(checkpoint_handle* check, int line_index,
		grid_parms grid, double* y, celltype1** smc) {
	/***********************************************************************************/
	MPI_Offset disp;
	MPI_Status status;

	int elements = grid.num_smc_circumferentially * grid.num_smc_axially;

	double** buffer;
	buffer = (double**) checked_malloc(grid.neq_smc * sizeof(double*),
			"error allocating memory for read buffer in reinit_smc.");
	for (int i = 0; i < grid.neq_smc; i++) {
		buffer[i] = (double*) checked_malloc(elements * sizeof(double),
				"error allocating memory for read buffer in reinit_smc.");
	}

	disp = ((line_index - 1) * grid.tasks * elements * sizeof(double))
			+ (grid.rank * elements * sizeof(double));
	check_flag(
			MPI_File_read_at(check->ci, disp, buffer[smc_Ca], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit Ca data in reinit_smc.");
	check_flag(
			MPI_File_read_at(check->si, disp, buffer[smc_SR], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit SR data in reinit_smc.");
	check_flag(
			MPI_File_read_at(check->vi, disp, buffer[smc_Vm], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit Vm data in reinit_smc.");
	check_flag(
			MPI_File_read_at(check->wi, disp, buffer[smc_w], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit K_Ca open channel probability data in reinit_smc.");
	check_flag(
			MPI_File_read_at(check->Ii, disp, buffer[smc_IP3], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit IP3 data in reinit_smc.");

	int k = 0, offset;
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;
			y[k + ((j - 1) * grid.neq_smc) + smc_Ca] = buffer[smc_Ca][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_SR] = buffer[smc_SR][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_Vm] = buffer[smc_Vm][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_w] = buffer[smc_w][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_IP3] = buffer[smc_IP3][(i - 1)
					* grid.num_smc_axially + (j - 1)];
		}
	}
	return (y);
}
/********************************************************************************/
double* reinitialize_tsoukias_smc(checkpoint_handle* check, int line_index,
		grid_parms grid, double* y, celltype1** smc) {
	/***********************************************************************************/
	MPI_Offset disp;
	MPI_Status status;

	int elements = grid.num_smc_circumferentially * grid.num_smc_axially;

	double** buffer;
	buffer = (double**) checked_malloc(grid.neq_smc * sizeof(double*),
			"error allocating memory for read buffer in reinit_smc.");
	for (int i = 0; i < grid.neq_smc; i++) {
		buffer[i] = (double*) checked_malloc(elements * sizeof(double),
				"error allocating memory for read buffer in reinit_smc.");
	}

	disp = ((line_index - 1) * grid.tasks * elements * sizeof(double))
			+ (grid.rank * elements * sizeof(double));
	check_flag(
			MPI_File_read_at(check->ci, disp, buffer[smc_Ca], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit Ca data in reinit_smc.");
	check_flag(
			MPI_File_read_at(check->si, disp, buffer[smc_SR], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit SR data in reinit_smc.");
	check_flag(
			MPI_File_read_at(check->vi, disp, buffer[smc_Vm], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit Vm data in reinit_smc.");
	check_flag(
			MPI_File_read_at(check->wi, disp, buffer[smc_w], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit K_Ca open channel probability data in reinit_smc.");
	check_flag(
			MPI_File_read_at(check->Ii, disp, buffer[smc_IP3], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit IP3 data in reinit_smc.");

	int k = 0, offset;
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;

			y[k + ((j - 1) * grid.neq_smc) + smc_Vm] = buffer[smc_Vm][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_Ca] = buffer[smc_Ca][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_Ca_u] =
					buffer[smc_Ca_u][(i - 1) * grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_Ca_r] =
					buffer[smc_Ca_r][(i - 1) * grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_Na_i] =
					buffer[smc_Na_i][(i - 1) * grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_K_i] = buffer[smc_K_i][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_Cl_i] =
					buffer[smc_Cl_i][(i - 1) * grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_IP3] = buffer[smc_IP3][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_DAG] = buffer[smc_DAG][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_cGMP_i] = buffer[smc_cGMP_i][(i
					- 1) * grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_V_cGMP] = buffer[smc_V_cGMP][(i
					- 1) * grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_d_L] = buffer[smc_d_L][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_f_L] = buffer[smc_f_L][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_p_f] = buffer[smc_p_f][(i - 1)
					* grid.num_smc_axially + (j - 1)];

			y[k + ((j - 1) * grid.neq_smc) + smc_p_K] = buffer[smc_p_K][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_q_1] = buffer[smc_q_1][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_q_2] = buffer[smc_q_2][(i - 1)
					* grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_R_01] =
					buffer[smc_R_01][(i - 1) * grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_R_10] =
					buffer[smc_R_10][(i - 1) * grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_R_11] =
					buffer[smc_R_11][(i - 1) * grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_h_IP3] = buffer[smc_h_IP3][(i
					- 1) * grid.num_smc_axially + (j - 1)];

			y[k + ((j - 1) * grid.neq_smc) + smc_PIP2] =
					buffer[smc_PIP2][(i - 1) * grid.num_smc_axially + (j - 1)];

			y[k + ((j - 1) * grid.neq_smc) + smc_R_S_G] = buffer[smc_R_S_G][(i
					- 1) * grid.num_smc_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_smc) + smc_R_S_P_G] =
					buffer[smc_R_S_P_G][(i - 1) * grid.num_smc_axially + (j - 1)];
			y[k + (j * grid.neq_smc) + smc_G] = buffer[smc_G][(i - 1)
					* grid.num_smc_axially + (j - 1)];
		}
	}
	return (y);
}
/************************************************************************************/
double* reinitialize_koenigsberger_ec(checkpoint_handle* check, int line_index,
		grid_parms grid, double* y, celltype2** ec) {
	/************************************************************************************/
	MPI_Offset disp;
	MPI_Status status;

	int elements = grid.num_ec_circumferentially * grid.num_ec_axially;

	double** buffer;
	buffer = (double**) checked_malloc(grid.neq_ec * sizeof(double*),
			"error allocating memory for read buffer in reinit_ec.");
	for (int i = 0; i < grid.neq_ec; i++) {
		buffer[i] = (double*) checked_malloc(elements * sizeof(double),
				"error allocating memory for read buffer in reinit_ec.");
	}

	disp = ((line_index - 1) * grid.tasks * elements * sizeof(double))
			+ (grid.rank * elements * sizeof(double));

	check_flag(
			MPI_File_read_at(check->cj, disp, buffer[ec_Ca], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit Ca data in reinit_ec.");
	check_flag(
			MPI_File_read_at(check->sj, disp, buffer[ec_SR], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit SR data in reinit_ec.");
	check_flag(
			MPI_File_read_at(check->vj, disp, buffer[ec_Vm], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit Vm data in reinit_ec.");
	check_flag(
			MPI_File_read_at(check->Ij, disp, buffer[ec_IP3], elements,
					MPI_DOUBLE, &status),
			"error read in the reinit IP3 data in reinit_ec.");

	int k = 0, offset = (grid.neq_smc * grid.num_smc_circumferentially
			* grid.num_smc_axially);

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			if (i > 1)
				k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
				k = offset + 0;
			y[k + ((j - 1) * grid.neq_ec) + ec_Ca] = buffer[ec_Ca][(i - 1)
					* grid.num_ec_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_ec) + ec_SR] = buffer[ec_SR][(i - 1)
					* grid.num_ec_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_ec) + ec_Vm] = buffer[ec_Vm][(i - 1)
					* grid.num_ec_axially + (j - 1)];
			y[k + ((j - 1) * grid.neq_ec) + ec_IP3] = buffer[ec_IP3][(i - 1)
					* grid.num_ec_axially + (j - 1)];
		}
	}
	return (y);
}

/**
 *
 * \param t_stamp
 * \param grid
 * \param smc
 * \param ec
 * \param cpl_cef
 * \param t
 * \param y
 * \param f
 * \return
 */
int compute_with_time_profiling(time_stamps* t_stamp, grid_parms grid,
		celltype1** smc, celltype2** ec, conductance cpl_cef, double t,
		double* y, double* f) {
	int err;

	t_stamp->map_function_t1 = MPI_Wtime();
	map_solver_to_cells(grid, y, smc, ec);
	t_stamp->map_function_t2 = MPI_Wtime();
	t_stamp->diff_map_function = t_stamp->diff_map_function
			+ (t_stamp->map_function_t2 - t_stamp->map_function_t1);

	t_stamp->single_cell_fluxes_t1 = MPI_Wtime();
	switch (grid.smc_model) {
	case (TSK): {
		tsoukias_smc(grid, smc);
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
	switch (grid.ec_model) {
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
	t_stamp->single_cell_fluxes_t2 = MPI_Wtime();
	t_stamp->diff_single_cell_fluxes = t_stamp->diff_single_cell_fluxes
			+ (t_stamp->single_cell_fluxes_t2 - t_stamp->single_cell_fluxes_t1);

	t_stamp->coupling_fluxes_t1 = MPI_Wtime();
	coupling(t, y, grid, smc, ec, cpl_cef);
	t_stamp->coupling_fluxes_t2 = MPI_Wtime();
	t_stamp->diff_coupling_fluxes = t_stamp->diff_coupling_fluxes
			+ (t_stamp->coupling_fluxes_t2 - t_stamp->coupling_fluxes_t1);

	//tsoukias_smc_derivatives(f, grid, smc);
	//koenigsberger_ec_derivatives(t, f, grid, ec);
	switch (grid.smc_model) {
	case (TSK): {
		tsoukias_smc_derivatives(f, grid, smc);
		break;
	}
	case (KNBGR): {
		koenigsberger_smc_derivatives(f, grid, smc);
		break;
	}
	default: {
		err = 1;
		break;
	}
	}
	switch (grid.ec_model) {
	case (TSK): {
		koenigsberger_ec_derivatives(t, f, grid, ec);
		break;
	}
	case (KNBGR): {
		koenigsberger_ec_derivatives(t, f, grid, ec);
		break;
	}
	default: {
		err = 1;
		break;
	}
	}
	return (err);
}

/*****************************************************************************/
int compute(grid_parms grid, celltype1** smc, celltype2** ec,
		conductance cpl_cef, double t, double* y, double* f) {
	/*****************************************************************************/

	int err;

	map_solver_to_cells(grid, y, smc, ec);

	switch (grid.smc_model) {
	case (TSK): {
		tsoukias_smc(grid, smc);
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
	switch (grid.ec_model) {
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
	coupling(t, y, grid, smc, ec, cpl_cef);

	switch (grid.smc_model) {
	case (TSK): {
		tsoukias_smc_derivatives(f, grid, smc);
		break;
	}
	case (KNBGR): {
		koenigsberger_smc_derivatives(f, grid, smc);
		break;
	}
	default: {
		err = 1;
		break;
	}
	}
	switch (grid.ec_model) {
	case (TSK): {
		koenigsberger_ec_derivatives(t, f, grid, ec);
		break;
	}
	case (KNBGR): {
		koenigsberger_ec_derivatives(t, f, grid, ec);
		break;
	}
	default: {
		err = 1;
		break;
	}
	}
	return (err);
}

/*****************************************************************************/
void Total_cells_in_computational_domain(grid_parms grid) {
/// Gathers local information on each CPU about the number of ECs and SMCs and sends to the Root
/// which evaluates the total number of ECs and SMCs constituting the global computational domain.
	/****************************************************************************/
	int sendcount = 2, recvcount = 2;
	int root = 0;
	int sendarray[2], recvarray[sendcount * grid.numtasks];
	int sumEC = 0, sumSMC = 0;
	sendarray[0] = grid.num_ec_axially * grid.num_ec_circumferentially;
	sendarray[1] = grid.num_smc_axially * grid.num_smc_circumferentially;

	check_flag(
			MPI_Gather(sendarray, sendcount, MPI_INT, recvarray, recvcount,
					MPI_INT, root, grid.universe), "error MPI_Gather");
	if (grid.universal_rank == 0) {
		for (int i = 0; i < grid.numtasks * recvcount; i += 2) {
			sumEC += recvarray[i];
		}
		for (int i = 1; i < grid.numtasks * recvcount; i += 2) {
			sumSMC += recvarray[i];
		}
		printf("Total number of ECs = %d\n", sumEC);
		printf("Total number of SMCs = %d\n", sumSMC);
		printf("Total number of ECs+SMCs = %d\n", sumEC + sumSMC);
	}
}

/****************************************************************************/
void process_time_profiling_data(grid_parms grid, double** time_profiler,
		int count) {
///This function is suppose to evaluate min, max, average of the stored time profiling data
// for every processor.
	double processed_data[11][3]; ///This will hold min max average of each of the 11 arrays.
	int index, root = 0, sendcount = 3, recvcount = grid.numtasks * sendcount;
	for (int n = 0; n < 11; n++) {
		minimum(time_profiler[n], count, &processed_data[n][0], &index);
		maximum(time_profiler[n], count, &processed_data[n][1], &index);
		average(time_profiler[n], count, &processed_data[n][2]);
	}

	/*printf("[%d] %2.15lf\t%2.15lf\t%2.15lf\n", grid.universal_rank,
	 processed_data[3][0], processed_data[3][1],
	 processed_data[3][2]);*/

	double *sendarray = (double*) malloc(sendcount * sizeof(double));
	double *recvarray = (double*) malloc(recvcount * sizeof(double));
	sendarray[0] = processed_data[3][0];
	sendarray[1] = processed_data[3][1];
	sendarray[2] = processed_data[3][2];
	check_flag(
			MPI_Gather(sendarray, sendcount, MPI_DOUBLE, recvarray, recvcount,
					MPI_DOUBLE, root, grid.universe), "error MPI_Gather");

	if (grid.universal_rank == 0) {
		FILE* fw;
		fw = fopen("time_profile_data.txt", "w+");
		for (int i = 0; i < grid.numtasks; i++) {
			fprintf(fw, "%d %lf %lf %lf\n", i, recvarray[i * sendcount],
					recvarray[i * sendcount + 1], recvarray[i * sendcount + 2]);
		}
		fclose(fw);
	}
	MPI_Barrier(grid.universe);
}
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
