#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"
#define CHECK(fn) {int errcode; errcode = (fn); if (errcode != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD,1); }
//#ifdef PARALLEL_IO
checkpoint_handle* initialise_checkpoint(grid_parms grid){

	checkpoint_handle *check = (checkpoint_handle*)malloc(sizeof(checkpoint_handle));
	char filename[50],suffix[10];
	int err;
	

	//Prepare the suffix which indicates my subdomain information and if I am a bifurcation, then also tells about which branch do I belong to
	int subdomain, branch;

	if (grid.my_domain.internal_info.domain_type == STRSEG){
		subdomain	=	grid.my_domain.internal_info.domain_index;
		err=sprintf(suffix,"%d.txt",subdomain);
	}else if (grid.my_domain.internal_info.domain_type == BIF){
		subdomain	=	grid.my_domain.internal_info.domain_index;
		branch		=	grid.branch_tag;
		err=sprintf(suffix,"%d_%d.txt",subdomain,branch);
	}

	err=sprintf(filename,"time_%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Time));

	err=sprintf(filename,"Log_file%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->logptr));

	err=sprintf(filename,"smc_Ca%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->ci));

	err=sprintf(filename,"ec_Ca%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cj));

	err=sprintf(filename,"smc_SERCA%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->si));

	err=sprintf(filename,"ec_SERCA%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->sj));

	err=sprintf(filename,"smc_V%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->vi));

	err=sprintf(filename,"ec_V%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->vj));

	err=sprintf(filename,"smc_KCa%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->wi));

	err=sprintf(filename,"smc_IP3%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Ii));

	err=sprintf(filename,"ec_IP3%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Ij));

	err=sprintf(filename,"smc_cpc%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpCi));

	err=sprintf(filename,"ec_cpc%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpCj));

	err=sprintf(filename,"smc_cpV%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpVi));

	err=sprintf(filename,"ec_cpV%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpVj));

	err=sprintf(filename,"smc_cpIP3%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpIi));

	err=sprintf(filename,"ec_cpIP3%s",suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpIj));
	
	err=sprintf(filename,"Elasped_time%s",suffix);
        CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->elapsed_time));
	return (check);
}

void dump_smc(grid_parms grid, celltype1 **smc, checkpoint_handle *check,
		int write_count) {
	MPI_Status status[8];
	MPI_Offset disp;
	int write_element_count, time_offset_in_file;

	write_element_count = grid.num_smc_axially * grid.num_smc_circumferentially;
	time_offset_in_file = write_count * write_element_count * grid.tasks
			* sizeof(double);

	double b1[grid.num_smc_circumferentially * grid.num_smc_axially],
			b2[grid.num_smc_circumferentially * grid.num_smc_axially],
			b3[grid.num_smc_circumferentially * grid.num_smc_axially],
			b4[grid.num_smc_circumferentially * grid.num_smc_axially],
			b5[grid.num_smc_circumferentially * grid.num_smc_axially],
			b6[grid.num_smc_circumferentially * grid.num_smc_axially],
			b7[grid.num_smc_circumferentially * grid.num_smc_axially],
			b8[grid.num_smc_circumferentially * grid.num_smc_axially];
	int k;

	k = 0;
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			b1[k] = smc[i][j].p[smc_Ca];
			b2[k] = smc[i][j].p[smc_SR];
			b3[k] = smc[i][j].p[smc_Vm];
			b4[k] = smc[i][j].p[smc_w];
			b5[k] = smc[i][j].p[smc_IP3];
			b6[k] = smc[i][j].B[cpl_Ca];
			b7[k] = smc[i][j].B[cpl_Vm];
			b8[k] = smc[i][j].B[cpl_IP3];
			k++;
		}
	}
	disp = time_offset_in_file
			+ (grid.rank * write_element_count * sizeof(double));
	CHECK(
			MPI_File_write_at(check->ci, disp, &b1, write_element_count, MPI_DOUBLE, &status[0]));
	CHECK(
			MPI_File_write_at(check->si, disp, &b2, write_element_count, MPI_DOUBLE, &status[1]));
	CHECK(
			MPI_File_write_at(check->vi, disp, &b3, write_element_count, MPI_DOUBLE, &status[2]));
	CHECK(
			MPI_File_write_at(check->wi, disp, &b4, write_element_count, MPI_DOUBLE, &status[3]));
	CHECK(
			MPI_File_write_at(check->Ii, disp, &b5, write_element_count, MPI_DOUBLE, &status[4]));
	CHECK(
			MPI_File_write_at(check->cpCi, disp, &b6, write_element_count, MPI_DOUBLE, &status[5]));
	CHECK(
			MPI_File_write_at(check->cpVi, disp, &b7, write_element_count, MPI_DOUBLE, &status[6]));
	CHECK(
			MPI_File_write_at(check->cpIi, disp, &b8, write_element_count, MPI_DOUBLE, &status[7]));
}

void dump_ec(grid_parms grid, celltype2 **ec, checkpoint_handle *check, int write_count){
	MPI_Status status[8];
	MPI_Offset disp;
	int write_element_count, time_offset_in_file;

	write_element_count = grid.num_ec_axially * grid.num_ec_circumferentially;
	time_offset_in_file = write_count * write_element_count * grid.tasks
			* sizeof(double);

	double b1[grid.num_ec_circumferentially * grid.num_ec_axially],
			b2[grid.num_ec_circumferentially * grid.num_ec_axially],
			b3[grid.num_ec_circumferentially * grid.num_ec_axially],
			b4[grid.num_ec_circumferentially * grid.num_ec_axially],
			b5[grid.num_ec_circumferentially * grid.num_ec_axially],
			b6[grid.num_ec_circumferentially * grid.num_ec_axially],
			b7[grid.num_ec_circumferentially * grid.num_ec_axially];
	int k;

	k = 0;
	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			b1[k] = ec[i][j].q[ec_Ca];
			b2[k] = ec[i][j].q[ec_SR];
			b3[k] = ec[i][j].q[ec_Vm];
			b4[k] = ec[i][j].q[ec_IP3];
			b5[k] = ec[i][j].B[cpl_Ca];
			b6[k] = ec[i][j].B[cpl_Vm];
			b7[k] = ec[i][j].B[cpl_IP3];
			k++;
		}
	}
	disp = time_offset_in_file
			+ (grid.rank * write_element_count * sizeof(double));
	CHECK(
			MPI_File_write_at(check->cj, disp, &b1, write_element_count, MPI_DOUBLE, &status[0]));
	CHECK(
			MPI_File_write_at(check->sj, disp, &b2, write_element_count, MPI_DOUBLE, &status[1]));
	CHECK(
			MPI_File_write_at(check->vj, disp, &b3, write_element_count, MPI_DOUBLE, &status[2]));
	CHECK(
			MPI_File_write_at(check->Ij, disp, &b4, write_element_count, MPI_DOUBLE, &status[3]));
	CHECK(
			MPI_File_write_at(check->cpCj, disp, &b5, write_element_count, MPI_DOUBLE, &status[4]));
	CHECK(
			MPI_File_write_at(check->cpVj, disp, &b6, write_element_count, MPI_DOUBLE, &status[5]));
	CHECK(
			MPI_File_write_at(check->cpIj, disp, &b7, write_element_count, MPI_DOUBLE, &status[6]));
}


/*******************/
/*Asnyc MPI-IO test*/
/*******************/
void dump_smc_async(grid_parms grid, celltype1 **smc, checkpoint_handle *check,
		int write_count) {
	MPI_Status status;
	MPI_Request request[8];
	MPI_Offset disp;
	int write_element_count, time_offset_in_file;

	write_element_count = grid.num_smc_axially * grid.num_smc_circumferentially;
	time_offset_in_file = write_count * write_element_count * grid.tasks
			* sizeof(double);

	double b1[grid.num_smc_circumferentially * grid.num_smc_axially],
			b2[grid.num_smc_circumferentially * grid.num_smc_axially],
			b3[grid.num_smc_circumferentially * grid.num_smc_axially],
			b4[grid.num_smc_circumferentially * grid.num_smc_axially],
			b5[grid.num_smc_circumferentially * grid.num_smc_axially],
			b6[grid.num_smc_circumferentially * grid.num_smc_axially],
			b7[grid.num_smc_circumferentially * grid.num_smc_axially],
			b8[grid.num_smc_circumferentially * grid.num_smc_axially];
	int k;

	k = 0;
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			b1[k] = smc[i][j].p[smc_Ca];
			b2[k] = smc[i][j].p[smc_SR];
			b3[k] = smc[i][j].p[smc_Vm];
			b4[k] = smc[i][j].p[smc_w];
			b5[k] = smc[i][j].p[smc_IP3];
			b6[k] = smc[i][j].B[cpl_Ca];
			b7[k] = smc[i][j].B[cpl_Vm];
			b8[k] = smc[i][j].B[cpl_IP3];
			k++;
		}
	}
	disp = time_offset_in_file
			+ (grid.rank * write_element_count * sizeof(double));
	CHECK(
			MPI_File_iwrite_at(check->ci, disp, &b1, write_element_count, MPI_DOUBLE, &request[0]));
	CHECK(
			MPI_File_iwrite_at(check->si, disp, &b2, write_element_count, MPI_DOUBLE, &request[1]));
	CHECK(
			MPI_File_iwrite_at(check->vi, disp, &b3, write_element_count, MPI_DOUBLE, &request[2]));
	CHECK(
			MPI_File_iwrite_at(check->wi, disp, &b4, write_element_count, MPI_DOUBLE, &request[3]));
	CHECK(
			MPI_File_iwrite_at(check->Ii, disp, &b5, write_element_count, MPI_DOUBLE, &request[4]));
	CHECK(
			MPI_File_iwrite_at(check->cpCi, disp, &b6, write_element_count, MPI_DOUBLE, &request[5]));
	CHECK(
			MPI_File_iwrite_at(check->cpVi, disp, &b7, write_element_count, MPI_DOUBLE, &request[6]));
	CHECK(
			MPI_File_iwrite_at(check->cpIi, disp, &b8, write_element_count, MPI_DOUBLE, &request[7]));
}

void dump_ec_async(grid_parms grid, celltype2 **ec, checkpoint_handle *check, int write_count){
	MPI_Status status;
	MPI_Request request[8];
	MPI_Offset disp;
	int write_element_count, time_offset_in_file;

	write_element_count = grid.num_ec_axially * grid.num_ec_circumferentially;
	time_offset_in_file = write_count * write_element_count * grid.tasks
			* sizeof(double);

	double b1[grid.num_ec_circumferentially * grid.num_ec_axially],
			b2[grid.num_ec_circumferentially * grid.num_ec_axially],
			b3[grid.num_ec_circumferentially * grid.num_ec_axially],
			b4[grid.num_ec_circumferentially * grid.num_ec_axially],
			b5[grid.num_ec_circumferentially * grid.num_ec_axially],
			b6[grid.num_ec_circumferentially * grid.num_ec_axially],
			b7[grid.num_ec_circumferentially * grid.num_ec_axially];
	int k;

	k = 0;
	for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid.num_ec_axially; j++) {
			b1[k] = ec[i][j].q[ec_Ca];
			b2[k] = ec[i][j].q[ec_SR];
			b3[k] = ec[i][j].q[ec_Vm];
			b4[k] = ec[i][j].q[ec_IP3];
			b5[k] = ec[i][j].B[cpl_Ca];
			b6[k] = ec[i][j].B[cpl_Vm];
			b7[k] = ec[i][j].B[cpl_IP3];
			k++;
		}
	}
	disp = time_offset_in_file
			+ (grid.rank * write_element_count * sizeof(double));

	CHECK(
			MPI_File_iwrite_at(check->cj, disp, &b1, write_element_count, MPI_DOUBLE, &request[0]));
	CHECK(
			MPI_File_iwrite_at(check->sj, disp, &b2, write_element_count, MPI_DOUBLE, &request[1]));
	CHECK(
			MPI_File_iwrite_at(check->vj, disp, &b3, write_element_count, MPI_DOUBLE, &request[2]));
	CHECK(
			MPI_File_iwrite_at(check->Ij, disp, &b4, write_element_count, MPI_DOUBLE, &request[3]));
	CHECK(
			MPI_File_iwrite_at(check->cpCj, disp, &b5, write_element_count, MPI_DOUBLE, &request[4]));
	CHECK(
			MPI_File_iwrite_at(check->cpVj, disp, &b6, write_element_count, MPI_DOUBLE, &request[5]));
	CHECK(
			MPI_File_iwrite_at(check->cpIj, disp, &b7, write_element_count, MPI_DOUBLE, &request[6]));
}

/*******************/


void checkpoint(checkpoint_handle* check, grid_parms grid, double tnow, celltype1** smc, celltype2** ec, int write_count)
{
	MPI_Status	status;
	MPI_Offset	disp;
	disp = (write_count*1*grid.numtasks*sizeof(double))+ (grid.rank*sizeof(double));
	CHECK(MPI_File_write_at(check->Time, disp, &tnow, 1, MPI_DOUBLE, &status));

	dump_smc(grid, smc, check, write_count);
	dump_ec(grid, ec, check, write_count);

}

void dump_rank_info(checkpoint_handle *check, conductance cpl_cef,
		grid_parms grid) {
	MPI_Status status;
	MPI_Offset disp;
	int bytes = 2 * 1024; 			// 2kB space
	char* buffer;
	buffer = (char*) checked_malloc(bytes,
			"allocation for logfile segment space\n");

sprintf(buffer,
			"BRANCH_TAG	=	%d\n(Universal_Rank, Cart_Rank= (%d,%d) \tcoords= %d,%d\t nbrs: local (u,d,l,r)=(%d %d %d %d) \t "
					"remote: (up1,up2,down1,down2)=(%d %d %d %d)\n\n flip_array: (%d,%d,%d,%d)\n\n"
					"Boundary_tag = %c\n(T = Top\t B= Bottom\t N=Interior of the subdomain)\n"
					"COUPLING COEFFICIENTS\n"
					"Vm_hm_smc=%2.5lf\nVm_hm_ec=%2.5lf\nCa_hm_smc=%2.5lf\nCa_hm_ec=%2.5lf\nIP3_hm_smc=%2.5lf\n"
					"IP3_hm_ec=%2.5lf\nVm_ht_smc=%2.5lf\nVm_ht_ec=%2.5lf\nCa_ht_smc=%2.5lf\nCa_ht_ec=%2.5lf\n"
					"IP3_ht_smc=%2.5lf\nIP3_ht_ec=%2.5lf\n\n"
					"Spatial Gradient info:\nUniform JPLC\t=%2.5lf\nMinimum JPLC\t=%2.5lf\nMaximum JPLC\t=%2.5lf\nGradient\t=%2.5lf\n"
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
					"Total number of equations in the full computational domain =%d\n "
					"------------------------------------------------------------------",

			grid.branch_tag, grid.universal_rank, grid.rank, grid.coords[0],
			grid.coords[1], grid.nbrs[local][UP], grid.nbrs[local][DOWN],
			grid.nbrs[local][LEFT], grid.nbrs[local][RIGHT],
			grid.nbrs[remote][UP1], grid.nbrs[remote][UP2],
			grid.nbrs[remote][DOWN1], grid.nbrs[remote][DOWN2],
			grid.flip_array[0], grid.flip_array[1], grid.flip_array[2],
			grid.flip_array[3], grid.my_domain.internal_info.boundary_tag,
			cpl_cef.Vm_hm_smc, cpl_cef.Vm_hm_ec, cpl_cef.Ca_hm_smc,
			cpl_cef.Ca_hm_ec, cpl_cef.IP3_hm_smc, cpl_cef.IP3_hm_ec,
			cpl_cef.Vm_ht_smc, cpl_cef.Vm_ht_ec, cpl_cef.Ca_ht_smc,
			cpl_cef.Ca_ht_ec, cpl_cef.IP3_ht_smc, cpl_cef.IP3_ht_ec,
			grid.uniform_jplc, grid.min_jplc, grid.max_jplc, grid.gradient,
			grid.numtasks,grid.m,grid.n,
			grid.num_ec_axially, grid.num_smc_circumferentially,
			grid.num_ec_axially * grid.num_ec_circumferentially,
			grid.num_smc_axially * grid.num_smc_circumferentially,
			(grid.num_ec_axially * grid.num_ec_circumferentially)
					+ (grid.num_smc_axially * grid.num_smc_circumferentially),
			(grid.num_ec_circumferentially * grid.num_ec_axially * grid.numtasks),
			(grid.num_smc_circumferentially * grid.num_smc_axially
					* grid.numtasks),
			((grid.num_ec_axially * grid.num_ec_circumferentially)
					+ (grid.num_smc_axially * grid.num_smc_circumferentially))
					* grid.numtasks, grid.NEQ * grid.numtasks);        






	disp	=	grid.rank * bytes;

	CHECK(MPI_File_write_at(check->logptr, disp, buffer, bytes, MPI_CHAR, &status));

}

void final_checkpoint(grid_parms grid, checkpoint_handle *check,double t1, double t2){
	MPI_Status	status;
	MPI_Offset	disp;
	double diff = t2-t1;

	disp = grid.rank*sizeof(double);
	CHECK(MPI_File_write_at_all(check->elapsed_time, disp, &diff, 1, MPI_DOUBLE, &status));

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
	MPI_File_close(&check->elapsed_time);
	
}


/*

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
*/





