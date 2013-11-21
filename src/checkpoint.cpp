#include "computelib.h"
#define CHECK(fn) {int errcode; errcode = (fn);if (errcode != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD,1); }
//#ifdef PARALLEL_IO
checkpoint_handle* initialise_checkpoint(grid_parms grid) {

	checkpoint_handle *check = (checkpoint_handle*) malloc(sizeof(checkpoint_handle));
	char filename[50], suffix[10];
	int err;

	open_common_checkpoint(check, grid);

	return (check);
}
/*********************************************************************/
void open_common_checkpoint(checkpoint_handle* check, grid_parms grid) {
	/*********************************************************************/

	int err;
	char filename[50];
	err = sprintf(filename, "Log_file%s.txt", grid.suffix);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->logptr));

	err = sprintf(filename, "Elasped_time%s.txt", grid.suffix);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->elapsed_time));

	err = sprintf(filename, "JPLC%s.txt", grid.suffix);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->jplc));

	err = sprintf(filename, "coords%s.txt", grid.suffix);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR,MPI_INFO_NULL, &check->coords));
}

/*********************************************************************/
void open_koenigsberger_smc_checkpoint(checkpoint_handle* check, grid_parms grid, int write_count, char* path) {
	/*********************************************************************/
	int err;
	char filename[50];
	err = sprintf(filename, "%s/smc_Ca%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->ci));
	err = sprintf(filename, "%s/smc_SERCA%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->si));

	err = sprintf(filename, "%s/smc_V%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->vi));

	err = sprintf(filename, "%s/smc_KCa%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->wi));

	err = sprintf(filename, "%s/smc_IP3%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Ii));
}

/*************************************************************************/
void open_koenigsberger_ec_checkpoint(checkpoint_handle* check, grid_parms grid, int write_count, char* path) {
	/*************************************************************************/
	int err;
	char filename[50];
	err = sprintf(filename, "%s/ec_Ca%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cj));

	err = sprintf(filename, "%s/ec_SERCA%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->sj));
	err = sprintf(filename, "%s/ec_V%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->vj));
	err = sprintf(filename, "%s/ec_IP3%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Ij));
}

/*************************************************************************/
void open_coupling_data_checkpoint(checkpoint_handle* check, grid_parms grid, int write_count, char* path) {
	/*************************************************************************/
	int err;
	char filename[50];
	err = sprintf(filename, "%s/smc_cpc%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpCi));

	err = sprintf(filename, "%s/ec_cpc%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpCj));

	err = sprintf(filename, "%s/smc_cpV%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpVi));

	err = sprintf(filename, "%s/ec_cpV%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpVj));

	err = sprintf(filename, "%s/smc_cpIP3%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpIi));

	err = sprintf(filename, "%s/ec_cpIP3%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpIj));

}
/***************************************************************************/
void dump_smc(grid_parms grid, celltype1 **smc, checkpoint_handle *check, int line_number, int write_count) {
	/***************************************************************************/
	MPI_Status status[8];
	MPI_Request req[8];

	MPI_Status status_tmp[3];
	MPI_Request req_tmp[3];

	MPI_Offset disp;
	int write_element_count, time_offset_in_file, file_offset;

	write_element_count = grid.num_smc_axially * grid.num_smc_circumferentially;
	time_offset_in_file = write_count * write_element_count * grid.tasks * sizeof(double);
	file_offset = (line_number * grid.tasks * write_element_count * sizeof(double));

	double b1[grid.num_smc_circumferentially * grid.num_smc_axially], b2[grid.num_smc_circumferentially * grid.num_smc_axially],
			b3[grid.num_smc_circumferentially * grid.num_smc_axially], b4[grid.num_smc_circumferentially * grid.num_smc_axially],
			b5[grid.num_smc_circumferentially * grid.num_smc_axially], b6[grid.num_smc_circumferentially * grid.num_smc_axially],
			b7[grid.num_smc_circumferentially * grid.num_smc_axially], b8[grid.num_smc_circumferentially * grid.num_smc_axially];
	int k;

	k = 0;
	// An amendment in the orientation of data written, so that it is aligned with the vtk geomerty
	for (int j = 1; j <= grid.num_smc_axially; j++) {
		for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
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
	disp = /*file_offset + time_offset_in_file + */(grid.rank * write_element_count * sizeof(double));
	CHECK( MPI_File_write_at(check->ci, disp, &b1, write_element_count, MPI_DOUBLE, &status[0]));
	CHECK( MPI_File_write_at(check->vi, disp, &b3, write_element_count, MPI_DOUBLE, &status[2]));
	CHECK( MPI_File_write_at(check->Ii, disp, &b5, write_element_count, MPI_DOUBLE, &status[4]));
	CHECK( MPI_File_write_at(check->si, disp, &b2, write_element_count, MPI_DOUBLE, &status[1]));

	CHECK( MPI_File_write_at(check->wi, disp, &b4, write_element_count, MPI_DOUBLE, &status[3]));

	CHECK( MPI_File_write_at(check->cpCi, disp, &b6, write_element_count, MPI_DOUBLE, &status[5]));
	CHECK( MPI_File_write_at(check->cpVi, disp, &b7, write_element_count, MPI_DOUBLE, &status[6]));
	CHECK( MPI_File_write_at(check->cpIi, disp, &b8, write_element_count, MPI_DOUBLE, &status[7]));
}
/***********************************************************************/
void dump_ec(grid_parms grid, celltype2 **ec, checkpoint_handle *check, int line_number, int write_count) {
	/***********************************************************************/
	MPI_Status status[8];

	MPI_Status status_tmp[3];
	MPI_Request req_tmp[3];

	MPI_Offset disp;
	int write_element_count, time_offset_in_file, file_offset;

	write_element_count = grid.num_ec_axially * grid.num_ec_circumferentially;
	time_offset_in_file = write_count * write_element_count * grid.tasks * sizeof(double);
	file_offset = (line_number * grid.tasks * write_element_count * sizeof(double));

	double b1[grid.num_ec_circumferentially * grid.num_ec_axially], b2[grid.num_ec_circumferentially * grid.num_ec_axially],
			b3[grid.num_ec_circumferentially * grid.num_ec_axially], b4[grid.num_ec_circumferentially * grid.num_ec_axially],
			b5[grid.num_ec_circumferentially * grid.num_ec_axially], b6[grid.num_ec_circumferentially * grid.num_ec_axially],
			b7[grid.num_ec_circumferentially * grid.num_ec_axially];
	int k;

	k = 0;
// An amendment in the orientation of data written, so that it is aligned with the vtk geomerty
	/*for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
	 for (int j = 1; j <= grid.num_ec_axially; j++) {*/
	for (int j = 1; j <= grid.num_ec_axially; j++) {
		for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
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
	disp = /*file_offset + time_offset_in_file + */(grid.rank * write_element_count * sizeof(double));
	CHECK( MPI_File_write_at(check->cj, disp, &b1, write_element_count, MPI_DOUBLE, &status[0]));

	CHECK( MPI_File_write_at(check->vj, disp, &b3, write_element_count, MPI_DOUBLE, &status[2]));
	CHECK( MPI_File_write_at(check->Ij, disp, &b4, write_element_count, MPI_DOUBLE, &status[3]));
	CHECK( MPI_File_write_at(check->sj, disp, &b2, write_element_count, MPI_DOUBLE, &status[1]));
	CHECK( MPI_File_write_at(check->cpCj, disp, &b5, write_element_count, MPI_DOUBLE, &status[4]));
	CHECK( MPI_File_write_at(check->cpVj, disp, &b6, write_element_count, MPI_DOUBLE, &status[5]));
	CHECK( MPI_File_write_at(check->cpIj, disp, &b7, write_element_count, MPI_DOUBLE, &status[6]));
}

/*******************/
/*Asnyc MPI-IO test*/
/*******************/

void dump_data(checkpoint_handle* check, grid_parms grid, int line_number, double tnow, celltype1** smc, celltype2** ec, int write_count) {
	MPI_Status status;
	MPI_Offset disp;

	int write_element_count = 1, time_offset_in_file = (write_count * 1 * grid.tasks * sizeof(double)), file_offset = (line_number * grid.tasks * 1
			* sizeof(double));

	disp = file_offset + time_offset_in_file + (grid.rank * sizeof(double));
	CHECK(MPI_File_write_at(check->Time, disp, &tnow, 1, MPI_DOUBLE, &status));

	dump_smc(grid, smc, check, line_number, write_count);
	dump_ec(grid, ec, check, line_number, write_count);

}

void update_line_number(checkpoint_handle* check, grid_parms grid, int line_number) {
	MPI_Status status;
	MPI_Offset disp;

	disp = grid.universal_rank * sizeof(int);
	CHECK( MPI_File_write_at(check->line_number, disp, &line_number, 1, MPI_INT, &status));
}

void dump_rank_info(checkpoint_handle *check, conductance cpl_cef, grid_parms grid) {
	MPI_Status status;
	MPI_Offset disp;
	int bytes = 2 * 1024; 			// 2kB space
	char* buffer;
	buffer = (char*) checked_malloc(bytes, "allocation for logfile segment space\n");

	sprintf(buffer,
			"BRANCH_TAG	=	%d\n[Universal_Rank, Cart_Rank= (%d,%d)] \tcoords= %d,%d\t nbrs: local (u,d,l,r)=(%d %d %d %d) \t "
					"remote: (up1,up2,down1,down2)=(%d %d %d %d)\n\n flip_array: (%d,%d,%d,%d)\n\n"
					"Boundary_tag = %c\n(T = Top\t B= Bottom\t I=Interior edges of the bifurcation segmensts, parent or children\t N=Interior of the subdomain)\n"
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
					"z_coordinates:       start = %lf     end = %lf\n local_z_start = %lf  local_z_end = %lf\n"

					"------------------------------------------------------------------",

			grid.branch_tag, grid.universal_rank, grid.rank, grid.coords[0], grid.coords[1], grid.nbrs[local][UP], grid.nbrs[local][DOWN],
			grid.nbrs[local][LEFT], grid.nbrs[local][RIGHT], grid.nbrs[remote][UP1], grid.nbrs[remote][UP2], grid.nbrs[remote][DOWN1],
			grid.nbrs[remote][DOWN2], grid.flip_array[0], grid.flip_array[1], grid.flip_array[2], grid.flip_array[3],
			grid.my_domain.internal_info.boundary_tag, cpl_cef.Vm_hm_smc, cpl_cef.Vm_hm_ec, cpl_cef.Ca_hm_smc, cpl_cef.Ca_hm_ec, cpl_cef.IP3_hm_smc,
			cpl_cef.IP3_hm_ec, cpl_cef.Vm_ht_smc, cpl_cef.Vm_ht_ec, cpl_cef.Ca_ht_smc, cpl_cef.Ca_ht_ec, cpl_cef.IP3_ht_smc, cpl_cef.IP3_ht_ec,
			grid.uniform_jplc, grid.min_jplc, grid.max_jplc, grid.gradient, grid.numtasks, grid.m, grid.n, grid.num_ec_axially,
			grid.num_smc_circumferentially, grid.num_ec_axially * grid.num_ec_circumferentially,
			grid.num_smc_axially * grid.num_smc_circumferentially,
			(grid.num_ec_axially * grid.num_ec_circumferentially) + (grid.num_smc_axially * grid.num_smc_circumferentially),
			(grid.num_ec_circumferentially * grid.num_ec_axially * grid.numtasks),
			(grid.num_smc_circumferentially * grid.num_smc_axially * grid.numtasks),
			((grid.num_ec_axially * grid.num_ec_circumferentially) + (grid.num_smc_axially * grid.num_smc_circumferentially)) * grid.numtasks,
			grid.NEQ * grid.numtasks, grid.my_domain.z_offset_start, grid.my_domain.z_offset_end, grid.my_domain.local_z_start,
			grid.my_domain.local_z_end);

	disp = grid.rank * bytes;

	CHECK( MPI_File_write_at(check->logptr, disp, buffer, bytes, MPI_CHAR, &status));

}

void dump_JPLC(grid_parms grid, celltype2 **ec, checkpoint_handle *check, const char *message) {

	/*	MPI_Status	status;
	 MPI_Offset	disp;
	 int write_element_count	=	grid.num_ec_circumferentially * grid.num_ec_axially;
	 double buffer[write_element_count];

	 int k = 0;
	 for (int i = 1; i <= grid.num_ec_circumferentially; i++) {
	 for (int j = 1; j <= grid.num_ec_axially; j++) {
	 buffer[k]	=	ec[i][j].JPLC;
	 k++;
	 }
	 }

	 int k = 0;
	 int i=1;
	 for (int j = 1; j <= grid.num_ec_axially; j++) {
	 buffer[k]	=	ec[i][j].JPLC;
	 k++;
	 }
	 int offset = grid.rank/grid.n;
	 //disp = (grid.rank * write_element_count * sizeof(double));
	 disp = (offset * write_element_count * sizeof(double));
	 CHECK(MPI_File_write_at(check->jplc, disp, &buffer, write_element_count, MPI_DOUBLE, &status));*/

	MPI_Status status;
	MPI_Offset disp;
	int write_element_count = grid.num_ec_axially;
	double buffer[write_element_count];

	int k = 0;
	int i = 1;
	/*for (int j = 1; j <= grid.num_ec_axially; j++) {
	 buffer[k] = ec[i][j].JPLC;
	 k++;
	 }
	 int offset = grid.rank / grid.n;
	 disp = (offset * write_element_count * sizeof(double));*/
	for (int j = grid.num_ec_axially; j >= 1; j--) {
		buffer[k] = ec[i][j].JPLC;
		k++;
	}
	disp = grid.num_ec_axially * (grid.m - ((grid.rank + 1) / grid.n)) * sizeof(double);
	CHECK( MPI_File_write_at(check->jplc, disp, &buffer, write_element_count, MPI_DOUBLE, &status));
}

void dump_coords(grid_parms grid, celltype2** ec, checkpoint_handle* check, const char* message) {

	MPI_Status status;
	MPI_Offset disp;
	int write_element_count = grid.num_ec_axially;
	double buffer[write_element_count];

	int k = 0;
	int i = 1;
	for (int j = grid.num_ec_axially; j >= 1; j--) {
		buffer[k] = ec[i][j].z_coord;
		k++;
	}

	int offset;

//offset = grid.rank / grid.n;
//disp = (offset * write_element_count * sizeof(double));

	disp = grid.num_ec_axially * (grid.m - ((grid.rank + 1) / grid.n)) * sizeof(double);
	CHECK( MPI_File_write_at(check->coords, disp, &buffer, write_element_count, MPI_DOUBLE, &status));
}

void checkpoint_timing_data(grid_parms grid, checkpoint_handle* check, double tnow, time_stamps t_stamp, int itteration, int append_point) {

	MPI_Status status;
	MPI_Offset disp_itteration, disp_write;
	int n = 16;
	double buffer[n];

	buffer[0] = tnow;
	buffer[1] = t_stamp.diff_async_comm_calls;
	buffer[2] = t_stamp.diff_async_comm_calls_wait;
	buffer[3] = t_stamp.diff_barrier_in_solver_before_comm;
	buffer[4] = t_stamp.diff_map_function;
	buffer[5] = t_stamp.diff_single_cell_fluxes;
	buffer[6] = t_stamp.diff_coupling_fluxes;
	buffer[7] = t_stamp.diff_solver;
	buffer[8] = t_stamp.diff_write;
	buffer[9] = (double) (t_stamp.computeDerivatives_call_counter);
	buffer[10] = (double) (itteration);

	buffer[11] = t_stamp.diff_remote_async_comm_calls;
	buffer[12] = t_stamp.diff_remote_async_comm_calls_wait;
	buffer[13] = t_stamp.diff_update_sendbuf;
	buffer[14] = t_stamp.diff_update_recvbuf;
	buffer[15] = t_stamp.diff_total_comms_cost;
	int write_element_count, time_offset_in_file, file_offset;

	file_offset = (append_point * grid.tasks * 1 * sizeof(double));

	write_element_count = 1;
	time_offset_in_file = itteration * write_element_count * grid.tasks * sizeof(double);

	disp_write = time_offset_in_file + file_offset + (grid.rank * write_element_count * sizeof(double));
	disp_itteration = (grid.rank * write_element_count * sizeof(int));

	CHECK( MPI_File_write_at_all(check->time_profiling, disp_write, &buffer[0], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->async_calls, disp_write, &buffer[1], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->async_wait, disp_write, &buffer[2], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->barrier_before_comm, disp_write, &buffer[3], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->map_function, disp_write, &buffer[4], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->single_cell_fluxes, disp_write, &buffer[5], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->coupling_fluxes, disp_write, &buffer[6], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->solver, disp_write, &buffer[7], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->writer_func, disp_write, &buffer[8], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->derivative_calls, disp_write, &buffer[9], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->itter_count, disp_write, &buffer[10], 1, MPI_DOUBLE, &status));
/// Write Comms time profiling data...
	CHECK( MPI_File_write_at_all(check->remote_async_calls, disp_write, &buffer[11], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->remote_async_wait, disp_write, &buffer[12], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->send_buf_update, disp_write, &buffer[13], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->recv_buf_update, disp_write, &buffer[14], 1, MPI_DOUBLE, &status));
	CHECK( MPI_File_write_at_all(check->total_comms_cost, disp_write, &buffer[15], 1, MPI_DOUBLE, &status));

}

void Record_timing_data_in_arrays(grid_parms grid, double tnow, time_stamps t_stamp, int itteration, double** time_profiler) {
	time_profiler[0][itteration] = tnow;
	time_profiler[1][itteration] = t_stamp.diff_async_comm_calls;
	time_profiler[2][itteration] = t_stamp.diff_async_comm_calls_wait;
	time_profiler[3][itteration] = t_stamp.diff_barrier_in_solver_before_comm;
	time_profiler[4][itteration] = t_stamp.diff_map_function;
	time_profiler[5][itteration] = t_stamp.diff_single_cell_fluxes;
	time_profiler[6][itteration] = t_stamp.diff_coupling_fluxes;
	time_profiler[7][itteration] = t_stamp.diff_solver;
	time_profiler[8][itteration] = t_stamp.diff_write;
	time_profiler[9][itteration] = (double) (t_stamp.computeDerivatives_call_counter);
	time_profiler[10][itteration] = (double) (itteration);
}

void final_checkpoint(checkpoint_handle *check, grid_parms grid) {

	MPI_Barrier(grid.universe);
	close_common_checkpoints(check);
//	close_time_wise_checkpoints(check);
//	close_time_profiling_checkpoints(check);
}

void close_common_checkpoints(checkpoint_handle* check) {
	MPI_File_close(&check->logptr);
	MPI_File_close(&check->elapsed_time);
	MPI_File_close(&check->jplc);
	MPI_File_close(&check->coords);
}
void close_time_wise_checkpoints(checkpoint_handle* check) {
	MPI_File_close(&check->Time);

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
void close_time_profiling_checkpoints(checkpoint_handle* check) {

	MPI_File_close(&check->time_profiling);
	MPI_File_close(&check->async_calls);
	MPI_File_close(&check->async_wait);
	MPI_File_close(&check->barrier_before_comm);
	MPI_File_close(&check->map_function);
	MPI_File_close(&check->single_cell_fluxes);
	MPI_File_close(&check->coupling_fluxes);
	MPI_File_close(&check->solver);
	MPI_File_close(&check->writer_func);
	MPI_File_close(&check->derivative_calls);
	MPI_File_close(&check->itter_count);
	MPI_File_close(&check->coords);
	MPI_File_close(&check->line_number);

	MPI_File_close(&check->remote_async_calls);
	MPI_File_close(&check->remote_async_wait);
	MPI_File_close(&check->send_buf_update);
	MPI_File_close(&check->recv_buf_update);
	MPI_File_close(&check->total_comms_cost);
}
/******************************************************************************/
int checkpoint(checkpoint_handle* check, grid_parms grid, double* tnow, double* y, celltype1** smc, celltype2** ec) {
	/******************************************************************************/
/// After when the MPI_IO files have been opened, check whether their current instance is first or did they previously existed.
/// This is checked by retrieving the file size of the file recording line number of the the timefile.
/// If the file is empty, then it is assumed to be the first instance of simulation (starting from t=0)
/// otherwise the linenumber indicates where the last complete result was written in the file and is used as
/// an offset to read initial values for the following simulation as well as a displacement for writing/appending new
/// data into the file.
	MPI_Offset line_number, file_offset;
	int line_index, err;
	check_flag(MPI_File_get_size(check->line_number, &line_number), "error reading the file size");
	if (line_number > 0) {
		line_index = recognize_end_of_file_index(check, grid);
		*tnow = reinitialize_time(check, line_index, grid);
		switch (grid.smc_model) {
		case (TSK): {
			y = reinitialize_tsoukias_smc(check, line_index, grid, y, smc);
			break;
		}
		case (KNBGR): {
			y = reinitialize_koenigsberger_smc(check, line_index, grid, y, smc);
			break;
		}
		default: {
			err = 1;
			break;
		}
		}
		switch (grid.ec_model) {
		case (TSK): {
			reinitialize_koenigsberger_ec(check, line_index, grid, y, ec);
			break;
		}
		case (KNBGR): {
			reinitialize_koenigsberger_ec(check, line_index, grid, y, ec);
			break;
		}
		default: {
			err = 1;
			break;
		}
		}

		y = reinitialize_koenigsberger_ec(check, line_index, grid, y, ec);
	} else if ((line_number == 0) || (line_number == NULL)) {
		switch (grid.smc_model) {
		case (TSK): {
			Initialize_tsoukias_smc(grid, y, smc);
			break;
		}
		case (KNBGR): {
			Initialize_koeingsberger_smc(grid, y, smc);
			break;
		}
		default: {
			err = 1;
			break;
		}
		}
		switch (grid.ec_model) {
		case (TSK): {
			Initialize_koeingsberger_ec(grid, y, ec);
			break;
		}
		case (KNBGR): {
			Initialize_koeingsberger_ec(grid, y, ec);
			break;
		}
		default: {
			err = 1;
			break;
		}
		}
	} else if (line_number < 0) {
		check_flag(MPI_Abort(grid.universe, 911), "error in locating the offset of for continued initial conditions for next time set.");
	}

	return (line_index);
}

/*********************************************************************/
int read_domain_info(int rank, char* filename, grid_parms* grid) {
	/*********************************************************************/
/// Read data from the domain_info.txt to retrieve  the information related to how the domain is set up.
/// All tasks open the same file to read.
/// Every task has the same displacement so each start to read data from the same position
/// Each task decide the read buffer size to be allocated by looking at the file size of the file opened
/// for read operation.
/// The data is read and sorted into delimiters and numbers and stored into corresponding place holder,
/// domains[][] in the structure grid_parms grid.
	int err;
	MPI_File data;
	MPI_Offset file_size, disp;
	MPI_Status status;
	int chunk, count;
	char* buffer;
	int *p;
	CHECK( MPI_File_open(grid->universe, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &data));

	err = MPI_File_get_size(data, &file_size);
	chunk = (int) (file_size / sizeof(char));
	buffer = (char*) malloc(chunk * sizeof(char));
	p = (int*) malloc(chunk * sizeof(int));
	if (p == NULL) {
		printf("[%d]: error allocating p\n", grid->universal_rank);
		MPI_Abort(grid->universe, 10001);
	}
	disp = 0;
	err = MPI_File_read_at(data, disp, buffer, chunk, MPI_CHAR, &status);
	err = MPI_Get_count(&status, MPI_CHAR, &count);
	char *token;
	char *delimiter = ";,\n";
	token = strtok(buffer, delimiter);
	int index = 0;
	while (token != NULL) {
		p[index] = atoi(token);
		token = strtok(NULL, delimiter);
		index++;
	}
	grid->num_domains = p[0];

	grid->domains = (int**) checked_malloc(grid->num_domains * sizeof(int*), "Subdomain information allocation");

/// second coordinate has following information:
/// Element 0: 	Key_val or serial number of the subdomain
/// Element 1:	Subdomain Type (2 possibilities and their values)
/// 				1. Straight Segment (STRSEG)					(0)
/// 				2. Bifurcation	(BIF)							(1)
/// Element 2:	Axial extent of processor of current key_val
/// Element 3: 	circumferential extent of processors of current key_val
/// Element 4:	Parent subdomain key_val of current Key_val.
/// Element 5: 	Left Child subdomain key_val of the current Key_val.
/// Element 6:   Right Child subdomain key_val of the current Key_val.
/// Element 7:   Required number of Endothelial cells axially per node.
/// Element 8:   Required number of Smooth muscle cells circumferentially per node.
///
/// In the case of elements 6 & 7, if subdomain type of current key_val is a straight segment, left Child is positive or zero, and right Child is negative.
/// If subdomain type of current key_val is a bifurcation, then both right and left child subdomains are non-negative.
////
	for (int i = 0; i < grid->num_domains; i++) {
		grid->domains[i] = (int*) checked_malloc(9 * sizeof(int), "Subdomains array elements allocation");
	}
	for (int i = 0; i < grid->num_domains; i++) {
		for (int j = 0; j < 9; j++) {
			grid->domains[i][j] = p[1 + (i * 9) + j];
		}
	}
	MPI_File_close(&data);
	return (0);
}

void update_elapsed_time(checkpoint_handle* check, grid_parms grid, time_keeper* elps_t) {

	MPI_Offset disp;
	MPI_Status status;
	char filename[50];

	disp = grid.rank * sizeof(double);

	elps_t->t_new = MPI_Wtime();
	elps_t->elapsed_time = elps_t->t_new - elps_t->t_old;

	if (elps_t > 0) {
		double time_from_file;
		///if elapsed time is non-zero then open the elapsed time file and read the previously written data
		///and updated it with new time; which is new elapsed time + previously written elapsed time.
		/*	int err = sprintf(filename, "Elasped_time%s", grid.suffix);
		 CHECK(
		 MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->elapsed_time));
		 */

		check_flag(MPI_File_read_at(check->elapsed_time, disp, &time_from_file, 1, MPI_DOUBLE, &status), "error read elapsed time from file.");

		time_from_file = time_from_file + elps_t->elapsed_time;

		check_flag(MPI_File_write_at(check->elapsed_time, disp, &time_from_file, 1, MPI_DOUBLE, &status), "error read elapsed time from file.");
	}
	elps_t->t_old = elps_t->t_new;
}

void naming_convention(grid_parms* grid) {
//Prepare the suffix which indicates my subdomain information and if I am a bifurcation, then also tells about which branch do I belong to
	int subdomain, branch, err;

	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		subdomain = grid->my_domain.internal_info.domain_index;
		err = sprintf(grid->suffix, "%d", subdomain);
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		subdomain = grid->my_domain.internal_info.domain_index;
		branch = grid->branch_tag;
		err = sprintf(grid->suffix, "%d_%d", subdomain, branch);
	}
}

int determine_file_offset_for_timing_data(checkpoint_handle* check, grid_parms grid) {
///This function implements checkpointing for time profiling data.
///If the simulation somehow crashes or gets discontinued, the file offset will be determined
///by reading in the last iteration number written in the file.
///This is to avoid overwriting on the pre-existing data and returned.

	MPI_Offset disp;
	MPI_Status status;
	char filename[50];
//int err = sprintf(filename, "itter_count%s", grid.suffix);
	int file_offset = 0;
	disp = grid.rank * sizeof(int);
	/*CHECK(
	 MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->itter_count));
	 */
	CHECK( MPI_File_read_at(check->itter_count, disp, &file_offset, 1, MPI_INT, &status));

	return (file_offset);

}

void jplc_plot_data(grid_parms grid, checkpoint_handle* check) {

	MPI_Offset disp;
	MPI_Status status;
	double z, jplc;
	FILE *fh;
	char filename[50];
	int err = sprintf(filename, "agonist%s.txt", grid.suffix);

	fh = fopen(filename, "w+");
	for (int i = 0; i < (grid.m * grid.num_ec_axially); i++) {
		disp = i * sizeof(double);
		int count = grid.num_ec_axially * grid.m;
		CHECK( MPI_File_read_at(check->coords, disp, &z, 1, MPI_DOUBLE, &status));
		CHECK( MPI_File_read_at(check->jplc, disp, &jplc, 1,MPI_DOUBLE, &status));
		fprintf(fh, "%2.8lf\t%2.8lf\n", z, jplc);
	}
	fclose(fh);

}

checkpoint_handle* initialise_time_wise_checkpoint(checkpoint_handle* check, grid_parms grid, int write_count, char* path) {
	int err;
	char filename[50];

	err = sprintf(filename, "%s/time_%s_t_%d.txt", path, grid.suffix, write_count);
	CHECK( MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Time));

	open_koenigsberger_smc_checkpoint(check, grid, write_count, path);
	open_koenigsberger_ec_checkpoint(check, grid, write_count, path);
	open_coupling_data_checkpoint(check, grid, write_count, path);
	return (check);
}

