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

			if (t > 100.00) {

				ec[i][j].JPLC =grid.min_jplc+ (grid.max_jplc/
						(1 + exp(-grid.gradient * ( ((j-1)+grid.num_ec_axially*floor(grid.rank/grid.n)) -(grid.m*grid.num_ec_axially / 2) )) ) );

			} else if (t <= 100.00)
			ec[i][j].JPLC = grid.uniform_jplc;

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
	tend= interval;
	rksuite.setup(total, tnow, y, tend, TOL, thres, method, "CT", false,
			0.0, false);

	communication_async_send_recv(check->logptr,grid,sendbuf,recvbuf,smc,ec);
	computeDerivatives(tnow, y, yp);
	MPI_Barrier (MPI_COMM_WORLD);
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
		communication_async_send_recv(check->logptr,grid,sendbuf,recvbuf,smc,ec);
		MPI_Barrier(MPI_COMM_WORLD);

		if (itteration == 5) {
			dump_JPLC(grid, ec, check, "Local agonist before t=100s\n");
		}

		if (itteration == 1e5) {
			dump_JPLC(grid, ec, check, "Local agonist after t=100s\n");
		}

		if ((itteration % file_write_per_unit_time) == 0) {
			checkpoint(check, grid, tnow, smc, ec,write_count);
		write_count++;
		}		//end itteration
		tend += interval;
		rksuite.reset(tend);
	}			//end while()

}
/*
void rksuite_solver_UT(double tnow, double tfinal, double interval, double *y, double* yp,
		int total, double TOL, double* thres, int file_write_per_unit_time,
		checkpoint_handle *check) {

	RKSUITE rksuite;
	//Solver method
	int method = 2;		//RK(4,5)
	//Error Flag
	int uflag = 0;
	int itteration = 0;
	double tend;

	double* ymax = (double*) checked_malloc(grid.NEQ * sizeof(double), stdout,
			"Solver array ymax for RKSUITE");
	communication_async_send_recv(check->logptr);

	rksuite.setup(grid.NEQ, tnow, y, tfinal, TOL, thres, method, "UT", false,
			0.00, false);

	double t1 = MPI_Wtime();
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
		MPI_Barrier (MPI_COMM_WORLD);
		communication_async_send_recv(check->logptr);

		if (itteration == 5) {
			dump_JPLC(grid, ec, check, "Local agonist before t=100s");
		}

		if (itteration == 1e5) {
			dump_JPLC(grid, ec, check, "Local agonist after t=100s");
		}

		if ((itteration % file_write_per_unit_time) == 0) {
			checkpoint_with_ghost_cells(check, grid, tnow, smc, ec);
		}		//end itteration
	}		//end of for loop on TEND

}*/
