#ifndef CVODE
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"

extern "C" {
#include "rksuite.h"
}

extern conductance cpl_cef;
extern celltype1** smc;
extern celltype2** ec;
extern double **sendbuf, **recvbuf;
extern grid_parms grid;


///***************************************************************************************/
///************ComputeDerivates(tnow, state_variables[], first_derivatives[])*************/
///***************************************************************************************/
void computeDerivatives(double t, double y[], double f[]) {


	const double gama = 1970.00;
	const double lambda = 45.00;
	const double Cmj = 25.8;
	int k = 0, talk = 0;
	map_solver_to_cells(grid,y,smc,ec);

	single_cell( t, y,  grid, smc,  ec);
	coupling(t,y, grid, smc, ec,cpl_cef);
	///In previous version of code (agonist_variation_main.cpp) the indexing in this section was
	///handelled differently but the end result is same.
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
			k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
			k = 0;
			f[k + ((j - 1) * grid.neq_smc) +smc_Ca]= smc[i][j].A[J_IP3] - smc[i][j].A[J_SERCA]
			+ smc[i][j].A[J_CICR] - smc[i][j].A[J_Extrusion]
			+ smc[i][j].A[J_Leak] - smc[i][j].A[J_VOCC]
			+ smc[i][j].A[J_Na_Ca] + smc[i][j].B[cpl_Ca]
			+ smc[i][j].C[cpl_Ca];

			f[k + ((j - 1) * grid.neq_smc) +smc_SR]= smc[i][j].A[J_SERCA] - smc[i][j].A[J_CICR]
			- smc[i][j].A[J_Leak];

			f[k + ((j - 1) * grid.neq_smc) +smc_Vm]= gama
			* (-smc[i][j].A[J_Na_K] - smc[i][j].A[J_Cl]
					- (2 * smc[i][j].A[J_VOCC]) - smc[i][j].A[J_Na_Ca]
					- smc[i][j].A[J_K]) + smc[i][j].B[cpl_Vm]
			+ smc[i][j].C[cpl_Vm];

			f[k + ((j - 1) * grid.neq_smc) +smc_w] = lambda
			* (smc[i][j].A[K_activation] - smc[i][j].p[smc_w]);

			f[k + ((j - 1) * grid.neq_smc) +smc_IP3] = -smc[i][j].A[J_IP3_deg]
			+ smc[i][j].B[cpl_IP3] + smc[i][j].C[cpl_IP3];
		}
	}

	int offset = (grid.neq_smc * grid.num_smc_circumferentially
			* grid.num_smc_axially);

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			if (i > 1)
			k = offset + ((i - 1) * grid.neq_ec_axially);
			else if (i == 1)
			k = offset + 0;

			ec[i][j].JPLC = agonist_profile(t,grid,i,j,ec[i][j].z_coord);

			f[k + ((j - 1) * grid.neq_ec) +ec_Ca] = ec[i][j].A[J_IP3]
			- ec[i][j].A[J_SERCA] + ec[i][j].A[J_CICR]
			- ec[i][j].A[J_Extrusion] + ec[i][j].A[J_Leak]
			+ ec[i][j].A[J_NSC] + ec[i][j].A[J_trivial_Ca]
			+ ec[i][j].B[cpl_Ca] + ec[i][j].C[cpl_Ca];

			f[k + ((j - 1) * grid.neq_ec) +ec_SR] = ec[i][j].A[J_SERCA]
			- ec[i][j].A[J_CICR] - ec[i][j].A[J_Leak];

			f[k + ((j - 1) * grid.neq_ec) +ec_Vm] = ((-1 / Cmj)
					* (ec[i][j].A[J_Ktot] + ec[i][j].A[J_Residual]))
			+ ec[i][j].B[cpl_Vm] + ec[i][j].C[cpl_Vm];

			f[k + ((j - 1) * grid.neq_ec) +ec_IP3] =	 ec[i][j].JPLC
			- ec[i][j].A[J_IP3_deg] + ec[i][j].B[cpl_IP3]
			+ ec[i][j].C[cpl_IP3];

		}
	}

}    //end of computeDerivatives()
void rksuite_solver_CT(double tnow, double tfinal, double interval, double *y, double* yp,
		int total, double TOL, double* thres, int file_write_per_unit_time,
		checkpoint_handle *check) {

	RKSUITE rksuite;
	//Solver method
	int method = 2;		//RK(4,5)
	double tend;
	int cflag = 0;
	int itteration = 0;
	int write_count=0;
	int write_once=0;

	tend= interval;
	rksuite.setup(total, tnow, y, tend, TOL, thres, method, "CT", false,
			0.0, false);

	communication_async_send_recv(grid,sendbuf,recvbuf,smc,ec);
	computeDerivatives(tnow, y, yp);
	MPI_Barrier (grid.universe);
	while (tnow <= tfinal) {
		// the ct() function does not guarantee to advance all the
		// way to the stop time.  Keep stepping until it does.
		do {
			rksuite.ct(computeDerivatives, tnow, y, yp, cflag);
			if (cflag >= 5) {
				fprintf(stdout,
						"[%d] \t RKSUITE failed with error flag %d at t=%lf\n\n",
						grid.rank, cflag, tnow);
				MPI_Abort(MPI_COMM_WORLD, 300);
			}
		} while (tnow < tend);

		///Increament the itteration as rksuite has finished solving between bounds tnow<= t <= tend.
		itteration++;
		/// Call for interprocessor communication
		MPI_Barrier(grid.universe);
		communication_async_send_recv(grid,sendbuf,recvbuf,smc,ec);


		/*if (itteration == 5) {
			dump_JPLC(grid, ec, check, "Local agonist before t=100s\n");
		}*/

		if ((write_once <= 1) && (tnow >= grid.stimulus_onset_time)) {
			write_once++;
			if (grid.rank % grid.n == 0) {
				dump_JPLC(grid, ec, check, "Local agonist after t=100s");
			}
		}
		if ((itteration % file_write_per_unit_time) == 0) {
			checkpoint(check, grid, tnow, smc, ec,write_count);
		write_count++;
		}		//end itteration
		//MPI_Barrier(grid.universe);
		tend += interval;
		rksuite.reset(tend);
	}			//end while()

}

void rksuite_solver_UT(double tnow, double tfinal, double interval, double *y, double* yp,
		int total, double TOL, double* thres, int file_write_per_unit_time,
		checkpoint_handle *check) {
	printf("[%d]: I have called RKSUITE\n",grid.universal_rank);

	RKSUITE rksuite;
	//Solver method
	int method = 2;		//RK(4,5)
	//Error Flag
	int uflag = 0;
	int itteration = 0;
	double tend;
	int write_count=0;
	int write_once=0;
	double* ymax = (double*) checked_malloc(grid.NEQ * sizeof(double),
			"Solver array ymax for RKSUITE");
	
	communication_async_send_recv(grid,sendbuf,recvbuf,smc,ec);
	computeDerivatives(tnow, y, yp);
	rksuite.setup(grid.NEQ, tnow, y, tfinal, TOL, thres, method, "UT", false,
			0.00, false);


	///Iterative  calls to the solver start here.
	for (tend = interval; tend < tfinal; tend += interval) {

		/// RKSUITE UT Call
		/// prototype (f,twant,tgot,ygot,ypgot,ymax,work,uflag)
		rksuite.ut(computeDerivatives, tend, tnow, y, yp, ymax, uflag);
		if (uflag >= 5) {
			fprintf(stdout,
					"[%d] \t RKSUITE failed with error flag %d at t=%lf\n\n",
					grid.rank, uflag, tnow);
			MPI_Abort(MPI_COMM_WORLD, 300);
		}

		///Increament the itteration as rksuite has finished solving between bounds tnow<= t <= tend.
		itteration++;
		/// Call for interprocessor communication
		MPI_Barrier(grid.universe);
		communication_async_send_recv(grid,sendbuf,recvbuf,smc,ec);


		/*if (itteration == 5) {
			dump_JPLC(grid, ec, check, "Local agonist before t=100s");
		}

		if ((itteration == int(grid.stimulus_onset_time/interval)) && (write_once<=1)) {
			write_once++;
			dump_JPLC(grid, ec, check, "Local agonist after t=100s");
		}*/

		if ((itteration % file_write_per_unit_time) == 0) {
					checkpoint(check, grid, tnow, smc, ec,write_count);
				write_count++;
				}		//end itteration
		//MPI_Barrier(grid.universe);
	}		//end of for loop on TEND

}
#endif
