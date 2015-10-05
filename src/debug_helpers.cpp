#include <mpi.h>
#include <malloc.h>

#include "computelib.h"
#include "koenigsberger_model.h"

/**
 * Dump info/debug output to a log file.
 */
void dump_rank_info(conductance cpl_cef, grid_parms grid)
{
	MPI_Status status;
	MPI_Offset displacement = 0;
	char *buffer = (char*)checked_malloc(2 * 1024 * sizeof(char), SRC_LOC);
	int root = 0;
	char filename[50];
	int logfile_displacements = 0;
	char *logfile_write_buffer = NULL;
	int length =
			sprintf(buffer,
					"BRANCH_TAG	= %d\n[Universal_Rank, Cart_Rank= (%d,%d)] \tcoords= %d,%d\t nbrs: local (u,d,l,r)=(%d %d %d %d)\t "
					"remote: (up,down)=(%d %d)\nflip_array: (%d,%d,%d,%d)\n"
					"Boundary_tag = %c\n(T = Top\t B= Bottom\t I=Interior edges of the bifurcation segments, parent or children\t N=Interior of the subdomain)\n"
					"COUPLING COEFFICIENTS\n"
					"Vm_hm_smc=%2.5lf\nVm_hm_ec=%2.5lf\nCa_hm_smc=%2.5lf\nCa_hm_ec=%2.5lf\nIP3_hm_smc=%2.5lf\n"
					"IP3_hm_ec=%2.5lf\nVm_ht_smc=%2.5lf\nVm_ht_ec=%2.5lf\nCa_ht_smc=%2.5lf\nCa_ht_ec=%2.5lf\n"
					"IP3_ht_smc=%2.5lf\nIP3_ht_ec=%2.5lf\n\n"
					"Total Tasks=%d\n"
					"Number of grid points in axial direction =%d\n"
					"Number of grid points in circumferential direction =%d\n"
					"Number of ECs per node (axially) =%d\n"
					"Number of SMCs per node (circumferentially) =%d\n"
					"Total ECs on this node =%d\n"
					"Total SMCs on this node =%d\n"
					"Total number of cells on this node =%d\n"
					"\nTotal ECs in the full computational domain =%d\n"
					"Total SMCs in the full computational domain =%d\n"
					"Total number of cells in the full computational domain =%d\n"
					"Total number of equations in the full computational domain =%d\n"
					"------------------------------------------------------------------\n\n",
					grid.branch_tag, grid.universal_rank, grid.rank_branch, grid.coords[0], grid.coords[1], grid.nbrs[local][UP], grid.nbrs[local][DOWN],
					grid.nbrs[local][LEFT], grid.nbrs[local][RIGHT], grid.nbrs[remote][UP], grid.nbrs[remote][DOWN],
					grid.flip_array[0], grid.flip_array[1], grid.flip_array[2], grid.flip_array[3],
					grid.boundary_tag, cpl_cef.Vm_hm_smc, cpl_cef.Vm_hm_ec, cpl_cef.Ca_hm_smc, cpl_cef.Ca_hm_ec,
					cpl_cef.IP3_hm_smc, cpl_cef.IP3_hm_ec, cpl_cef.Vm_ht_smc, cpl_cef.Vm_ht_ec, cpl_cef.Ca_ht_smc, cpl_cef.Ca_ht_ec,
					cpl_cef.IP3_ht_smc, cpl_cef.IP3_ht_ec, grid.num_ranks, grid.m,
					grid.n, grid.num_ec_axially, grid.num_smc_circumferentially, grid.num_ec_axially * grid.num_ec_circumferentially,
					grid.num_smc_axially * grid.num_smc_circumferentially,
					(grid.num_ec_axially * grid.num_ec_circumferentially) + (grid.num_smc_axially * grid.num_smc_circumferentially),
					(grid.num_ec_circumferentially * grid.num_ec_axially * grid.num_ranks),
					(grid.num_smc_circumferentially * grid.num_smc_axially * grid.num_ranks),
					((grid.num_ec_axially * grid.num_ec_circumferentially) + (grid.num_smc_axially * grid.num_smc_circumferentially)) * grid.num_ranks,
					grid.NEQ * grid.num_ranks);

	int *recv_count = (int*) checked_malloc(grid.num_ranks_branch * sizeof(int), SRC_LOC);
	int *disp = (int*) checked_malloc(grid.num_ranks_branch * sizeof(int), SRC_LOC);

	// Gathering and summing the length of all the CHARs contained in every send_buffer containing coordinates from each MPI process.
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid.cart_comm));

	for (int i = 0; i < grid.num_ranks_branch; i++)
	{
		disp[i] = logfile_displacements;
		logfile_displacements += recv_count[i];
	}

	if (grid.rank_branch == 0)
	{
		logfile_write_buffer = (char*) checked_malloc(logfile_displacements * sizeof(char), SRC_LOC);
	}
	CHECK_MPI_ERROR(MPI_Gatherv(buffer, length, MPI_CHAR, logfile_write_buffer, recv_count, disp, MPI_CHAR, root, grid.cart_comm));

	if (grid.rank_branch == 0)
	{

		sprintf(filename, "Logfile_%d_%d.txt", grid.domain_index, grid.branch_tag);

		printf("[%d] Writing %s\n", grid.universal_rank, filename);
		MPI_File rank_info_file;
		CHECK_MPI_ERROR(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &rank_info_file));
		CHECK_MPI_ERROR(MPI_File_write_at(rank_info_file, displacement, logfile_write_buffer, logfile_displacements, MPI_CHAR, &status));
		MPI_File_close(&rank_info_file);
	}

	if (grid.rank_branch == 0)
	{
		free(logfile_write_buffer);
	}

	free(buffer);
	free(recv_count);
	free(disp);
}

void print_domains(FILE* logptr, grid_parms grid, SMC_cell** smc, EC_cell** ec)
{
	fprintf(logptr, "****** SMC Domain ***********\n");
	for (int i = 0; i < (grid.num_smc_circumferentially + grid.num_ghost_cells);
			i++) {
		fprintf(logptr, " ------ i = %d -------------\n", i);
		for (int j = 0; j < (grid.num_smc_axially + grid.num_ghost_cells);
				j++) {
			fprintf(logptr,
					"[%d,%d]\t%2.5lf\t%2.5lf\t%2.5lf\t%2.5lf\t%2.5lf\n", i, j,
					smc[i][j].vars[smc_Ca], smc[i][j].vars[smc_SR],
					smc[i][j].vars[smc_Vm], smc[i][j].vars[smc_w],
					smc[i][j].vars[smc_IP3]);
		}
	}
	fprintf(logptr, "\n****** EC Domain ***********\n");
	for (int i = 0; i < (grid.num_ec_circumferentially + grid.num_ghost_cells);
			i++) {
		fprintf(logptr, " ------ i = %d -------------\n", i);
		for (int j = 0; j < (grid.num_ec_axially + grid.num_ghost_cells); j++) {
			fprintf(logptr, "[%d,%d]\t%2.5lf\t%2.5lf\t%2.5lf\t%2.5lf\n", i, j,
					ec[i][j].vars[ec_Ca], ec[i][j].vars[ec_SR], ec[i][j].vars[ec_Vm],
					ec[i][j].vars[ec_IP3]);
		}
	}
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
	if(grid.rank_branch == 0)
	{
		fprintf(logptr,"*** t = %lf***\n--------------SMC DOMAIN -------------\n",t);
		int kk, off;
		for (int i = 1; i <= grid.num_smc_circumferentially; i++) {

			for (int j = 1; j <= grid.num_smc_axially; j++) {
				if (i > 1)
				kk = ((i - 1) * grid.neq_smc_axially);
				else if (i == 1)
				kk = 0;

				fprintf(logptr,
						"SMC : [%d,%d]\t %2.3lf\t %2.3lf\t %2.3lf\t %2.3lf\t\n\n",
						i, j, smc[i][j].vars[smc_Vm],
						smc[i][j].homo_fluxes[cpl_Ca],smc[i][j].homo_fluxes[cpl_Vm],smc[i][j].homo_fluxes[cpl_IP3]);
			}
		}
		off = (grid.neq_smc * grid.num_smc_circumferentially
				* grid.num_smc_axially);

		for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
			for (int j = 1; j <= grid.num_ec_axially; j++) {
				if (i > 1)
				kk = off + ((i - 1) * grid.neq_ec_axially);
				else if (i == 1)
				kk = off + 0;
				fprintf(logptr,
						"EC : [%d,%d]\t %2.3lf\t %2.3lf\t %2.3lf\t %2.3lf\t\n\n",
						i, j, ec[i][j].vars[ec_Vm], ec[i][j].homo_fluxes[cpl_Ca],
						ec[i][j].homo_fluxes[cpl_Vm], ec[i][j].homo_fluxes[cpl_IP3]);
			}
		}
	}
}


