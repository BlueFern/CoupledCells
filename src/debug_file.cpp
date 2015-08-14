#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"


void print_domains(FILE* logptr, grid_parms grid, SMC_cell** smc,
		EC_cell** ec) {

	fprintf(logptr, "****** SMC Domain ***********\n");
	for (int i = 0; i < (grid.num_smc_circumferentially + grid.num_ghost_cells);
			i++) {
		fprintf(logptr, " ------ i = %d -------------\n", i);
		for (int j = 0; j < (grid.num_smc_axially + grid.num_ghost_cells);
				j++) {
			fprintf(logptr,
					"[%d,%d]\t%2.5lf\t%2.5lf\t%2.5lf\t%2.5lf\t%2.5lf\n", i, j,
					smc[i][j].p[smc_Ca], smc[i][j].p[smc_SR],
					smc[i][j].p[smc_Vm], smc[i][j].p[smc_w],
					smc[i][j].p[smc_IP3]);
		}
	}
	fprintf(logptr, "\n****** EC Domain ***********\n");
	for (int i = 0; i < (grid.num_ec_circumferentially + grid.num_ghost_cells);
			i++) {
		fprintf(logptr, " ------ i = %d -------------\n", i);
		for (int j = 0; j < (grid.num_ec_axially + grid.num_ghost_cells); j++) {
			fprintf(logptr, "[%d,%d]\t%2.5lf\t%2.5lf\t%2.5lf\t%2.5lf\n", i, j,
					ec[i][j].q[ec_Ca], ec[i][j].q[ec_SR], ec[i][j].q[ec_Vm],
					ec[i][j].q[ec_IP3]);
		}
	}

/*
	fprintf(logptr,"****** SMC Domain ***********\n");
	for (int i = 0; i < (grid.num_smc_circumferentially + grid.num_ghost_cells);
			i++) {
		fprintf(logptr," ------ i = %d -------------\n",i);
		for (int j = 0; j < (grid.num_smc_axially + grid.num_ghost_cells);
				j++) {
			fprintf(logptr, "[%d,%d]\t %2.5lf\t%2.5lf\t%2.5lf\t%2.5lf\t%2.5lf\t\t \t%2.5lf\t%2.5lf\t%2.5lf\n", i,
					j, smc[i][j].p[smc_Ca], smc[i][j].p[smc_SR],
					smc[i][j].p[smc_Vm], smc[i][j].p[smc_w],
					smc[i][j].p[smc_IP3],smc[i][j].B[cpl_Ca],smc[i][j].B[cpl_Vm],smc[i][j].B[cpl_IP3]);
		}
	}


	fprintf(logptr,"****** EC Domain ***********\n");
	for (int i = 0; i < (grid.num_ec_circumferentially + grid.num_ghost_cells);
			i++) {
		fprintf(logptr," ------ i = %d -------------\n",i);
		for (int j = 0; j < (grid.num_ec_axially + grid.num_ghost_cells);
				j++) {
			fprintf(logptr, "[%d,%d]\t %2.5lf\t%2.5lf\t%2.5lf\t%2.5lf\t\t \t%2.5lf\t%2.5lf\t%2.5lf\n", i,
					j, ec[i][j].q[ec_Ca], ec[i][j].q[ec_SR],
					ec[i][j].q[ec_Vm],ec[i][j].q[ec_IP3],ec[i][j].B[cpl_Ca],ec[i][j].B[cpl_Vm],ec[i][j].B[cpl_IP3]);
		}
	}*/
}

void print_send_buffer(FILE* logptr,grid_parms grid, double** sendbuf)
{
	fprintf(logptr, "***Up direction***\n");
	for (int i = 0; i < grid.num_elements_send_up; i++) {
		fprintf(logptr, "[%d]\t %lf\n", i, sendbuf[UP][i]);
	}
	fprintf(logptr, "***Down direction***\n");
	for (int i = 0; i < grid.num_elements_send_down; i++) {
		fprintf(logptr, "[%d]\t %lf\n", i, sendbuf[DOWN][i]);
	}
	fprintf(logptr, "***Left direction***\n");
	for (int i = 0; i < grid.num_elements_send_left; i++) {
		fprintf(logptr, "[%d]\t %lf\n", i, sendbuf[LEFT][i]);
	}
	fprintf(logptr, "***Right direction***\n");
	for (int i = 0; i < grid.num_elements_send_right; i++) {
		fprintf(logptr, "[%d]\t %lf\n", i, sendbuf[RIGHT][i]);
	}
}


void print_recv_buffer(FILE* logptr,grid_parms grid, double** recvbuf)
{
fprintf(logptr, "***Up direction***\n");
	for (int i = 0; i < grid.num_elements_recv_up; i++) {
		fprintf(logptr, "[%d]\t %2.8lf\n", i, recvbuf[UP][i]);
	}
	fprintf(logptr, "***Down direction***\n");
	for (int i = 0; i < grid.num_elements_recv_down; i++) {
		fprintf(logptr, "[%d]\t %2.8lf\n", i, recvbuf[DOWN][i]);
	}
	fprintf(logptr, "***Left direction***\n");
	for (int i = 0; i < grid.num_elements_recv_left; i++) {
		fprintf(logptr, "[%d]\t %2.8lf\n", i, recvbuf[LEFT][i]);
	}
	fprintf(logptr, "***Right direction***\n");
	for (int i = 0; i < grid.num_elements_recv_right; i++) {
		fprintf(logptr, "[%d]\t %2.8lf\n", i, recvbuf[RIGHT][i]);
	}
}

void print_compare(FILE* logptr, double t, double y[], grid_parms grid, SMC_cell** smc, EC_cell** ec)
{
	if(grid.rank==0){
		fprintf(logptr,"*** t = %lf***\n--------------SMC DOMAIN -------------\n",t);
			int kk, off;
			for (int i = 1; i <= grid.num_smc_circumferentially; i++) {

				for (int j = 1; j <= grid.num_smc_axially; j++) {
					if (i > 1)
					kk = ((i - 1) * grid.neq_smc_axially);
					else if (i == 1)
					kk = 0;
					/*fprintf(logptr,
							"[%d,%d]\t %2.3lf\t %2.3lf\t %2.3lf\t\t [%d] %2.3lf\t [%d] %2.3lf\t [%d] %2.3lf\n",
							i, j, smc[i][j].p[smc_Ca],smc[i][j].p[smc_Vm],smc[i][j].p[smc_IP3],
							kk + ((j - 1) * grid.neq_smc) + smc_Ca,y[kk + ((j - 1) * grid.neq_smc) + smc_Ca],
							kk + ((j - 1) * grid.neq_smc) + smc_Vm,y[kk + ((j - 1) * grid.neq_smc) + smc_Vm],
							kk + ((j - 1) * grid.neq_smc) + smc_IP3,y[kk + ((j - 1) * grid.neq_smc) + smc_IP3]);*/

					fprintf(logptr,
												"SMC : [%d,%d]\t %2.3lf\t %2.3lf\t %2.3lf\t %2.3lf\t\n\n",
												i, j, smc[i][j].p[smc_Vm],
												smc[i][j].B[cpl_Ca],smc[i][j].B[cpl_Vm],smc[i][j].B[cpl_IP3]);
					/*fprintf(logptr,
										"[%d,%d]\t %2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\t%2.3lf\n",
										i, j, smc[i][j].A[J_IP3], smc[i][j].A[J_SERCA],
										smc[i][j].A[J_CICR],smc[i][j].A[J_Extrusion],smc[i][j].A[J_Leak],smc[i][j].A[J_IP3_deg],
										smc[i][j].A[J_VOCC],smc[i][j].A[J_Na_Ca],smc[i][j].A[J_Na_K],smc[i][j].A[J_Cl],smc[i][j].A[J_K],smc[i][j].A[K_activation]);*/
				}
			}
			off = (grid.neq_smc * grid.num_smc_circumferentially
					* grid.num_smc_axially);
		//	fprintf(logptr,"--------------EC DOMAIN -------------\n");
			for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
				for (int j = 1; j <= grid.num_ec_axially; j++) {
					if (i > 1)
					kk = off + ((i - 1) * grid.neq_ec_axially);
					else if (i == 1)
					kk = off + 0;
			/*		fprintf(logptr,
							"[%d,%d]\t %2.3lf\t %2.3lf\t %2.3lf\t\t [%d] %2.3lf\t [%d] %2.3lf\t [%d] %2.3lf\n",
							i, j, ec[i][j].q[ec_Ca],ec[i][j].q[ec_Vm],ec[i][j].q[ec_IP3],
							kk + ((j - 1) * grid.neq_ec) + ec_Ca,y[kk + ((j - 1) * grid.neq_ec) + ec_Ca],
							kk + ((j - 1) * grid.neq_ec) + ec_Vm,y[kk + ((j - 1) * grid.neq_ec) + ec_Vm],
							kk + ((j - 1) * grid.neq_ec) + ec_IP3,y[kk + ((j - 1) * grid.neq_ec) + ec_IP3]);*/
				fprintf(logptr,
						"EC : [%d,%d]\t %2.3lf\t %2.3lf\t %2.3lf\t %2.3lf\t\n\n",
						i, j, ec[i][j].q[ec_Vm], ec[i][j].B[cpl_Ca],
						ec[i][j].B[cpl_Vm], ec[i][j].B[cpl_IP3]);
				}
			}
	}
}


