#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"
#define CHECK(fn) {int errcode; errcode = (fn); if (errcode != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD,1); }

#ifdef PARALLEL_IO
checkpoint_handle* initialise_checkpoint(int rank){
	checkpoint_handle *check = (checkpoint_handle*)malloc(sizeof(checkpoint_handle));
	char filename[30];
	int err;

	err=sprintf(filename,"time.txt");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Time));

	err=sprintf(filename,"Log_file.txt");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->logptr));

	err=sprintf(filename,"smc_Ca.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->ci));


	err=sprintf(filename,"ec_Ca.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cj));

	err=sprintf(filename,"smc_V.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->vi));

	err=sprintf(filename,"ec_V.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->vj));

	err=sprintf(filename,"smc_IP3.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Ii));

	err=sprintf(filename,"ec_IP3.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Ij));

	err=sprintf(filename,"smc_cpc.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpCi));

	err=sprintf(filename,"ec_cpc.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpCj));

	err=sprintf(filename,"smc_cpV.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpVi));

	err=sprintf(filename,"ec_cpV.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpVj));

	err=sprintf(filename,"smc_cpIP3.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpIi));

	err=sprintf(filename,"ec_cpIP3.vtk");
	CHECK(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpIj));

	return check;
}

void dump_smc(grid_parms grid, celltype1 **smc, checkpoint_handle *check, int write_count){
	MPI_Status	status;
	MPI_Offset	disp;
	int 	write_element_count,time_offset_in_file;
	double *buffer;


	write_element_count = grid.num_smc_axially;
	time_offset_in_file = write_count*write_element_count*grid.numtasks*sizeof(double);
	buffer = (double*) malloc(write_element_count* sizeof(double));

	int i=0;
	disp = time_offset_in_file + (grid.rank * write_element_count*sizeof(double));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = smc[i][j].p[smc_Ca];
	}
	CHECK(MPI_File_write_at_all(check->ci, disp, buffer, write_element_count, MPI_DOUBLE, &status));


	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = smc[i][j].p[smc_SR];
	}
	CHECK(MPI_File_write_at_all(check->si, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = smc[i][j].p[smc_Vm];
	}
	CHECK(MPI_File_write_at_all(check->vi, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = smc[i][j].p[smc_w];
	}
	CHECK(MPI_File_write_at_all(check->wi, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = smc[i][j].p[smc_IP3];
	}
	CHECK(MPI_File_write_at_all(check->Ii, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = smc[i][j].B[cpl_Ca];
	}
	CHECK(MPI_File_write_at_all(check->cpCi, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = smc[i][j].B[cpl_Vm];
	}
	CHECK(MPI_File_write_at_all(check->cpVi, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = smc[i][j].B[cpl_IP3];
	}
	CHECK(MPI_File_write_at_all(check->cpIi, disp, buffer, write_element_count, MPI_DOUBLE, &status));


	free(buffer);
}

void dump_ec(grid_parms grid, celltype2 **ec, checkpoint_handle *check, int write_count){
	MPI_Status	status;
	MPI_Offset	disp;
	int 	write_element_count,time_offset_in_file;
	double *buffer;

	write_element_count = grid.num_ec_axially;
	time_offset_in_file = write_count*write_element_count*grid.numtasks*sizeof(double);
	buffer = (double*) malloc(write_element_count* sizeof(double));

	int i=0;
	disp = time_offset_in_file + (grid.rank*write_element_count*sizeof(double));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = ec[i][j].q[smc_Ca];
	}
	CHECK(MPI_File_write_at_all(check->cj, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	disp = time_offset_in_file + (grid.rank*write_element_count*sizeof(double));
	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = ec[i][j].q[smc_SR];
	}
	CHECK(MPI_File_write_at_all(check->sj, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = ec[i][j].q[smc_Vm];
	}
	CHECK(MPI_File_write_at_all(check->vj, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = ec[i][j].q[smc_IP3];
	}
	CHECK(MPI_File_write_at_all(check->Ij, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = ec[i][j].B[cpl_Ca];
	}
	CHECK(MPI_File_write_at_all(check->cpCj, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = ec[i][j].B[cpl_Vm];
	}
	CHECK(MPI_File_write_at_all(check->cpVj, disp, buffer, write_element_count, MPI_DOUBLE, &status));

	for(int j=1; j<=grid.num_smc_axially;j++) {
		buffer[j] = ec[i][j].B[cpl_IP3];
	}
	CHECK(MPI_File_write_at_all(check->cpIj, disp, buffer, write_element_count, MPI_DOUBLE, &status));

}

void checkpoint(checkpoint_handle* check, grid_parms grid, double tnow, celltype1** smc, celltype2** ec, int write_count)
{
	MPI_Status	status;
	MPI_Offset	disp;
	disp = (write_count*1*grid.numtasks*sizeof(double))+ (grid.rank*sizeof(double));
	CHECK(MPI_File_write_at_all(check->Time, disp, &tnow, 1, MPI_DOUBLE, &status));

	dump_smc(grid, smc, check, write_count);
	dump_ec(grid, ec, check, write_count);

}


void final_checkpoint(checkpoint_handle *check, grid_parms grid,double t1, double t2){
	MPI_Status	status;
	MPI_Offset	disp;
	double diff = t2-t1;

	disp = grid.rank*sizeof(double);
	CHECK(MPI_File_write_at_all(check->logptr, disp, &diff, 1, MPI_DOUBLE, &status));

	MPI_File_close(&check->Time);
	MPI_File_close(&check->logptr);
	MPI_File_close(&check->ci);
	MPI_File_close(&check->cj);
	MPI_File_close(&check->si);
	MPI_File_close(&check->sj);
	MPI_File_close(&check->vi);
	MPI_File_close(&check->vj);
	MPI_File_close(&check->wi);
	MPI_File_close(&check->Ii);
	MPI_File_close(&check->Ij);

	MPI_File_close(&check->cpCi);
	MPI_File_close(&check->cpCj);
	MPI_File_close(&check->cpVi);
	MPI_File_close(&check->cpVj);
	MPI_File_close(&check->cpIi);
	MPI_File_close(&check->cpIj);
}




#else
checkpoint_handle* initialise_checkpoint(int rank){
	checkpoint_handle *check = (checkpoint_handle*)malloc(sizeof(checkpoint_handle));
	char filename[30];
	int err;
	err = sprintf(filename, "logfile%d.txt", rank);
	check->logptr=fopen(filename,"w+");

	err = sprintf(filename, "time%d.txt", rank);
	check->Time= fopen(filename,"w+");

	err = sprintf(filename, "smc_c%d.txt", rank);
	check->ci= fopen(filename,"w+");

	err = sprintf(filename, "smc_s%d.txt", rank);
	check->si= fopen(filename,"w+");

	err = sprintf(filename, "smc_v%d.txt", rank);
	check->vi= fopen(filename,"w+");

	err = sprintf(filename, "smc_w%d.txt", rank);
	check->wi= fopen(filename,"w+");

	err = sprintf(filename, "smc_I%d.txt", rank);
	check->Ii= fopen(filename,"w+");

	err = sprintf(filename, "ec_c%d.txt", rank);
	check->cj= fopen(filename,"w+");

	err = sprintf(filename, "ec_s%d.txt", rank);
	check->sj= fopen(filename,"w+");

	err = sprintf(filename, "ec_v%d.txt", rank);
	check->vj= fopen(filename,"w+");

	err = sprintf(filename, "ec_I%d.txt", rank);
	check->Ij= fopen(filename,"w+");

	err = sprintf(filename, "smc_cpC%d.txt", rank);
	check->cpCi= fopen(filename,"w+");

	err = sprintf(filename, "ec_cpC%d.txt", rank);
	check->cpCj= fopen(filename,"w+");

	err = sprintf(filename, "smc_cpV%d.txt", rank);
	check->cpVi= fopen(filename,"w+");

	err = sprintf(filename, "ec_cpV%d.txt", rank);
	check->cpVj= fopen(filename,"w+");

	err = sprintf(filename, "smc_cpI%d.txt", rank);
	check->cpIi= fopen(filename,"w+");

	err = sprintf(filename, "ec_cpI%d.txt", rank);
	check->cpIj= fopen(filename,"w+");

	return check;
}

void dump_smc(grid_parms grid, celltype1 **smc, checkpoint_handle *check, int write_count){

	for (int i =1; i <= grid.num_smc_circumferentially; i++) {
		for (int j =1; j <= grid.num_smc_axially; j++) {
			fprintf(check->ci , "%2.6lf\t",smc[i][j].p[smc_Ca] );
			fprintf(check->si , "%2.6lf\t",smc[i][j].p[smc_SR] );
			fprintf(check->vi , "%2.6lf\t",smc[i][j].p[smc_Vm] );
			fprintf(check->wi , "%2.6lf\t",smc[i][j].p[smc_w] );
			fprintf(check->Ii , "%2.6lf\t",smc[i][j].p[smc_IP3] );
			fprintf(check->cpCi , "%2.6lf\t",smc[i][j].B[cpl_Ca] );
			fprintf(check->cpVi , "%2.6lf\t",smc[i][j].B[cpl_Vm] );
			fprintf(check->cpIi , "%2.6lf\t",smc[i][j].B[cpl_IP3] );
		}		//end j
	}		//end i
	fprintf(check->ci ,"\n");
	fprintf(check->si ,"\n");
	fprintf(check->vi ,"\n");
	fprintf(check->wi ,"\n");
	fprintf(check->Ii ,"\n");
	fprintf(check->cpCi ,"\n");
	fprintf(check->cpVi ,"\n");
	fprintf(check->cpIi ,"\n");

}
void dump_JPLC(grid_parms grid, celltype2 **ec, checkpoint_handle *check, const char *message){
	fprintf(check->logptr, message);
	for (int p =1; p <=grid.num_ec_axially; p++) {
		fprintf(check->logptr, "[%d]\t%lf\n", p, ec[1][p].JPLC);
	}
}

void dump_ec(grid_parms grid, celltype2 **ec, checkpoint_handle *check, int write_count){

	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
			for (int j = 1; j <= grid.num_ec_axially; j++) {
				fprintf(check->cj , "%2.6lf\t",ec[i][j].q[ec_Ca] );
				fprintf(check->sj , "%2.6lf\t",ec[i][j].q[ec_SR] );
				fprintf(check->vj , "%2.6lf\t",ec[i][j].q[ec_Vm] );
				fprintf(check->Ij , "%2.6lf\t",ec[i][j].q[ec_IP3] );
				fprintf(check->cpCj , "%2.6lf\t",ec[i][j].B[cpl_Ca] );
				fprintf(check->cpVj , "%2.6lf\t",ec[i][j].B[cpl_Vm] );
				fprintf(check->cpIj , "%2.6lf\t",ec[i][j].B[cpl_IP3] );
			}		//end j
		}		//end i
		fprintf(check->cj ,"\n");
		fprintf(check->sj ,"\n");
		fprintf(check->vj ,"\n");
		fprintf(check->Ij ,"\n");
		fprintf(check->cpCj ,"\n");
		fprintf(check->cpVj ,"\n");
		fprintf(check->cpIj ,"\n");


}

void checkpoint(checkpoint_handle* check, grid_parms grid, double tnow, celltype1** smc, celltype2** ec, int write_count)
{
	fprintf(check->Time, "%2.6lf\n", tnow);
	dump_smc(grid, smc, check,write_count);
	dump_ec(grid, ec, check,write_count);

}


void final_checkpoint(checkpoint_handle *check, grid_parms, double t1, double t2){

	fprintf(check->logptr, "Elapsed time =%2.5lf\n\n", (t2 - t1));

	fclose(check->Time);
	fclose(check->logptr);
	fclose(check->ci);
	fclose(check->cj);
	fclose(check->si);
	fclose(check->sj);
	fclose(check->vi);
	fclose(check->vj);
	fclose(check->wi);
	fclose(check->Ii);
	fclose(check->Ij);

	fclose(check->cpCi);
	fclose(check->cpCj);
	fclose(check->cpVi);
	fclose(check->cpVj);
	fclose(check->cpIi);
	fclose(check->cpIj);

}

void dump_rank_info(checkpoint_handle *check, conductance cpl_cef,
		grid_parms grid) {

	fprintf(check->logptr,
			"BRANCH_TAG	=	%d\n(Universal_Rank, Cart_Rank= (%d,%d) \tcoords= %d,%d\t nbrs: local (u,d,l,r)=(%d %d %d %d) \t remote: (up1,up2,down1,down2)=(%d %d %d %d)\n\n flip_array: (%d,%d,%d,%d)\n\n",
			grid.branch_tag, grid.universal_rank, grid.rank, grid.coords[0],
			grid.coords[1], grid.nbrs[local][UP], grid.nbrs[local][DOWN],
			grid.nbrs[local][LEFT], grid.nbrs[local][RIGHT],
			grid.nbrs[remote][UP1], grid.nbrs[remote][UP2],
			grid.nbrs[remote][DOWN1], grid.nbrs[remote][DOWN2],
			grid.flip_array[0], grid.flip_array[1], grid.flip_array[2],
			grid.flip_array[3]);
	fprintf(check->logptr, "Boundary_tag = %c\n(T = Top\t B= Bottom\t N=Interior of the subdomain)\n",grid.my_domain.internal_info.boundary_tag);
	///Write my local information in my rank's logfile
	fprintf(check->logptr,
			"COUPLING COEFFICIENTS\nVm_hm_smc=%2.5lf\nVm_hm_ec=%2.5lf\nCa_hm_smc=%2.5lf\nCa_hm_ec=%2.5lf\nIP3_hm_smc=%2.5lf\nIP3_hm_ec=%2.5lf\nVm_ht_smc=%2.5lf\nVm_ht_ec=%2.5lf\nCa_ht_smc=%2.5lf\nCa_ht_ec=%2.5lf\nIP3_ht_smc=%2.5lf\nIP3_ht_ec=%2.5lf\n\n",
			cpl_cef.Vm_hm_smc, cpl_cef.Vm_hm_ec, cpl_cef.Ca_hm_smc,
			cpl_cef.Ca_hm_ec, cpl_cef.IP3_hm_smc, cpl_cef.IP3_hm_ec,
			cpl_cef.Vm_ht_smc, cpl_cef.Vm_ht_ec, cpl_cef.Ca_ht_smc,
			cpl_cef.Ca_ht_ec, cpl_cef.IP3_ht_smc, cpl_cef.IP3_ht_ec);

	fprintf(check->logptr,
			"Spatial Gradient info:\nUniform JPLC\t=%2.5lf\nMinimum JPLC\t=%2.5lf\nMaximum JPLC\t=%2.5lf\nGradient\t=%2.5lf\n",
			grid.uniform_jplc, grid.min_jplc, grid.max_jplc, grid.gradient);

	fprintf(check->logptr, "Total Tasks=%d\n", grid.numtasks);
	fprintf(check->logptr, "Number of grid points in axial direction =%d\n ",
			grid.m);
	fprintf(check->logptr,
			"Number of grid points in circumferential direction =%d\n ",
			grid.n);
	fprintf(check->logptr, "Number of ECs per node (axially) =%d\n ",
			grid.num_ec_axially);
	fprintf(check->logptr, "Number of SMCs per node (circumferentially) =%d\n ",
			grid.num_smc_circumferentially);
	fprintf(check->logptr, "Total ECs on this node =%d\n ",
			(grid.num_ec_axially * grid.num_ec_circumferentially));
	fprintf(check->logptr, "Total SMCs on this node =%d\n ",
			(grid.num_smc_axially * grid.num_smc_circumferentially));
	fprintf(check->logptr, "Total number of cells on this node =%d\n",
			(grid.num_ec_axially * grid.num_ec_circumferentially)
					+ (grid.num_smc_axially * grid.num_smc_circumferentially));
	fprintf(check->logptr,
			"Total number of cells in the full computational domain =%d\n ",
			((grid.num_ec_axially * grid.num_ec_circumferentially)
					+ (grid.num_smc_axially * grid.num_smc_circumferentially))
					* grid.numtasks);
	fprintf(check->logptr,
			"Total number of equations in the full computational domain =%d\n ",
			grid.NEQ * grid.numtasks);
}


void dump_smc_with_ghost_cells(grid_parms grid, celltype1 **smc, checkpoint_handle *check, int write_count){

	for (int i =0; i <(grid.num_smc_circumferentially+grid.num_ghost_cells); i++) {
		for (int j =0; j <(grid.num_smc_axially+grid.num_ghost_cells); j++) {
			fprintf(check->ci , "%2.6lf\t",smc[i][j].p[smc_Ca] );
			fprintf(check->si , "%2.6lf\t",smc[i][j].p[smc_SR] );
			fprintf(check->vi , "%2.6lf\t",smc[i][j].p[smc_Vm] );
			fprintf(check->wi , "%2.6lf\t",smc[i][j].p[smc_w] );
			fprintf(check->Ii , "%2.6lf\t",smc[i][j].p[smc_IP3] );

		}		//end j
	}		//end i
	fprintf(check->ci ,"\n");
	fprintf(check->si ,"\n");
	fprintf(check->vi ,"\n");
	fprintf(check->wi ,"\n");
	fprintf(check->Ii ,"\n");

}

void dump_ec_with_ghost_cells(grid_parms grid, celltype2 **ec, checkpoint_handle *check, int write_count){

	for (int i =0; i <(grid.num_ec_circumferentially+grid.num_ghost_cells); i++) {
			for (int j =0; j <(grid.num_ec_axially+grid.num_ghost_cells); j++) {
				fprintf(check->cj , "%2.6lf\t",ec[i][j].q[ec_Ca] );
				fprintf(check->sj , "%2.6lf\t",ec[i][j].q[ec_SR] );
				fprintf(check->vj , "%2.6lf\t",ec[i][j].q[ec_Vm] );
				fprintf(check->Ij , "%2.6lf\t",ec[i][j].q[ec_IP3] );
			}		//end j
		}		//end i
		fprintf(check->cj ,"\n");
		fprintf(check->sj ,"\n");
		fprintf(check->vj ,"\n");
		fprintf(check->Ij ,"\n");
}

void checkpoint_with_ghost_cells(checkpoint_handle* check, grid_parms grid, double tnow, celltype1** smc, celltype2** ec, int write_count)
{
	fprintf(check->Time, "%2.6lf\n", tnow);
	dump_smc_with_ghost_cells(grid, smc, check,write_count);
	dump_ec_with_ghost_cells(grid, ec, check,write_count);

}
#endif






