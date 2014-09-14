#include "computelib.h"
#define CHECK(fn) {int errcode; errcode = (fn);if (errcode != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD,1); }
//#ifdef PARALLEL_IO
checkpoint_handle* initialise_checkpoint(grid_parms grid) {

	checkpoint_handle *check = (checkpoint_handle*) malloc(sizeof(checkpoint_handle));

	return (check);
}
/*********************************************************************/
void open_common_checkpoint(checkpoint_handle* check, grid_parms grid) {
	/*********************************************************************/

	int err;
	char filename[50];
	err = sprintf(filename, "Log_file%s.txt", grid.suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->logptr));

	err = sprintf(filename, "Elasped_time%s.txt", grid.suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->elapsed_time));

	err = sprintf(filename, "JPLC%s.txt", grid.suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->jplc));

	err = sprintf(filename, "coords%s.txt", grid.suffix);
	CHECK(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->coords));
}

/*******************************************************************************************************************************************************/
void open_koenigsberger_smc_checkpoint(checkpoint_handle* check, grid_parms grid, int write_count, char* path, IO_domain_info* my_IO_domain_info)
/***************************************************************************************************************************************************/
{
	int err;
	char filename[50];
	err = sprintf(filename, "%s/smc_Data_t_%d.vtk", path, write_count);
	CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->smc_data_file));
	/*err = sprintf(filename, "%s/smc_Ca_t_%d.vtk", path, write_count);
	 CHECK( MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->ci));
	 err = sprintf(filename, "%s/smc_SERCA_t_%d.vtk", path, write_count);
	 CHECK( MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->si));
	 err = sprintf(filename, "%s/smc_V_t_%d.vtk", path, write_count);
	 CHECK( MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->vi));
	 err = sprintf(filename, "%s/smc_KCa_t_%d.vtk", path, write_count);
	 CHECK( MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->wi));
	 err = sprintf(filename, "%s/smc_IP3_t_%d.vtk", path, write_count);
	 CHECK( MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Ii));*/
}

/*************************************************************************/
void open_koenigsberger_ec_checkpoint(checkpoint_handle* check, grid_parms grid, int write_count, char* path, IO_domain_info* my_IO_domain_info) {
	/*************************************************************************/
	int err;
	char filename[50];
	err = sprintf(filename, "%s/ec_Data_t_%d.vtk", path, write_count);
	CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->ec_data_file));
	/*err = sprintf(filename, "%s/ec_Ca_t_%d.vtk", path, write_count);
	 CHECK( MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cj));
	 err = sprintf(filename, "%s/ec_SERCA_t_%d.vtk", path, write_count);
	 CHECK( MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->sj));
	 err = sprintf(filename, "%s/ec_V_t_%d.vtk", path, write_count);
	 CHECK( MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->vj));
	 err = sprintf(filename, "%s/ec_IP3_t_%d.vtk", path, write_count);
	 CHECK( MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Ij));*/
}

/*************************************************************************/
void open_coupling_data_checkpoint(checkpoint_handle* check, grid_parms grid, int write_count, char* path, IO_domain_info* my_IO_domain_info) {
	/*************************************************************************/
	int err;
	char filename[50];
	err = sprintf(filename, "%s/smc_cpc_t_%d.vtk", path, write_count);
	CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpCi));

	err = sprintf(filename, "%s/ec_cpc_t_%d.vtk", path, write_count);
	CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpCj));

	err = sprintf(filename, "%s/smc_cpV_t_%d.vtk", path, write_count);
	CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpVi));

	err = sprintf(filename, "%s/ec_cpV_t_%d.vtk", path, write_count);
	CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpVj));

	err = sprintf(filename, "%s/smc_cpIP3_t_%d.vtk", path, write_count);
	CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpIi));

	err = sprintf(filename, "%s/ec_cpIP3_t_%d.vtk", path, write_count);
	CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpIj));

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
	CHECK(MPI_File_write_at(check->ci, disp, &b1, write_element_count, MPI_DOUBLE, &status[0]));
	CHECK(MPI_File_write_at(check->vi, disp, &b3, write_element_count, MPI_DOUBLE, &status[2]));
	CHECK(MPI_File_write_at(check->Ii, disp, &b5, write_element_count, MPI_DOUBLE, &status[4]));
	CHECK(MPI_File_write_at(check->si, disp, &b2, write_element_count, MPI_DOUBLE, &status[1]));

	CHECK(MPI_File_write_at(check->wi, disp, &b4, write_element_count, MPI_DOUBLE, &status[3]));

	CHECK(MPI_File_write_at(check->cpCi, disp, &b6, write_element_count, MPI_DOUBLE, &status[5]));
	CHECK(MPI_File_write_at(check->cpVi, disp, &b7, write_element_count, MPI_DOUBLE, &status[6]));
	CHECK(MPI_File_write_at(check->cpIi, disp, &b8, write_element_count, MPI_DOUBLE, &status[7]));
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
	CHECK(MPI_File_write_at(check->cj, disp, &b1, write_element_count, MPI_DOUBLE, &status[0]));

	CHECK(MPI_File_write_at(check->vj, disp, &b3, write_element_count, MPI_DOUBLE, &status[2]));
	CHECK(MPI_File_write_at(check->Ij, disp, &b4, write_element_count, MPI_DOUBLE, &status[3]));
	CHECK(MPI_File_write_at(check->sj, disp, &b2, write_element_count, MPI_DOUBLE, &status[1]));
	CHECK(MPI_File_write_at(check->cpCj, disp, &b5, write_element_count, MPI_DOUBLE, &status[4]));
	CHECK(MPI_File_write_at(check->cpVj, disp, &b6, write_element_count, MPI_DOUBLE, &status[5]));
	CHECK(MPI_File_write_at(check->cpIj, disp, &b7, write_element_count, MPI_DOUBLE, &status[6]));
}

void dump_data(checkpoint_handle* check, grid_parms* grid, int line_number, double tnow, celltype1** smc, celltype2** ec, int write_count,
		IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer) {

	dump_smc_data(check, grid, my_IO_domain_info, writer_buffer, smc, write_count);
	dump_ec_data(check, grid, my_IO_domain_info, writer_buffer, ec, write_count);
}
/*********************************************************************************************************************************/
void dump_process_data(checkpoint_handle* check, grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, char* path)
/*********************************************************************************************************************************/
{
	MPI_Status status;
	MPI_Offset disp;
	char filename[50];
	int err = sprintf(filename, "%s/task_mesh.vtk", path);
	CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->task_mesh));

	int write_element_count = 1, header_offset[3] = { 0, 0, 0 }, point_offset = 0, cell_offset = 0, celltype_offset = 0, count = 0;
	char* header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0.\n");
	int branches;
	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		branches = 1;
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		branches = 3;
	}

	/*************** Writing VTK header **************/
	sprintf(header, "# vtk DataFile Version 2.0\n"
			"Task mesh show how MPI processes are connected\n"
			"ASCII\n"
			"DATASET UNSTRUCTURED_GRID\n"
			"POINTS %d double\n", grid->info[ProcessMesh][TOTAL_CELLS] * 4 * branches);

	header_offset[0] = strlen(header);
	count = header_offset[0];
	disp = 0;
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->task_mesh, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}

	/*************** Writing Point data **************/
	int buffer_lengths[4] = { 0, 0, 0, 0 };
	check_flag(MPI_Allgather(&writer_buffer->buffer_length[ProcessMesh], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths");

	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}

	disp = (header_offset[0] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[ProcessMesh];
	check_flag(MPI_File_write_at(check->task_mesh, disp, writer_buffer->process_mesh_points, count, MPI_CHAR, &status),
			"error writing the coordinates in task_mesh.\n");

	for (int i = 0; i < 4; i++) {
		point_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->process_mesh_points);
	/*************** Writing cell data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0.\n");
	header_offset[1] = sprintf(header, "CELLS %d %d\n", branches * grid->info[ProcessMesh][TOTAL_CELLS],
			5 * 3 * grid->info[ProcessMesh][TOTAL_CELLS]);

	count = header_offset[1];
	disp = (header_offset[0] + point_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->task_mesh, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.\n");
	}

	check_flag(MPI_Allgather(&writer_buffer->buffer_length[ProcessCell], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + disp + header_offset[1]) * sizeof(char);
	count = writer_buffer->buffer_length[ProcessCell];
	check_flag(MPI_File_write_at(check->task_mesh, disp, writer_buffer->process_mesh_cells, count, MPI_CHAR, &status),
			"error writing the cells in task_mesh.\n");
	for (int i = 0; i < 4; i++) {
		cell_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->process_mesh_cells);
	/*************** Writing cell type data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0.\n");
	header_offset[2] = sprintf(header, "CELL_TYPES %d\n", branches * grid->info[ProcessMesh][TOTAL_CELLS]);

	count = header_offset[2];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->task_mesh, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.\n");
	}
	check_flag(
			MPI_Allgather(&writer_buffer->buffer_length[ProcessCellType], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths\n");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[ProcessCellType];
	check_flag(MPI_File_write_at(check->task_mesh, disp, writer_buffer->process_mesh_type, count, MPI_CHAR, &status),
			"error writing the cell_type in task_mesh.\n");
	for (int i = 0; i < 4; i++) {
		celltype_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->process_mesh_type);

	MPI_File_close(&check->task_mesh);
}
/**************************************************************************************/
void dump_smc_data(checkpoint_handle* check, grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, celltype1** smc,
		int write_count) {
	MPI_Status status;
	MPI_Offset disp;

	int write_element_count = 1, *header_offset, point_offset = 0, cell_offset = 0, celltype_offset = 0, *smcDataOffset;
	header_offset = (int*) checked_malloc((3 + grid->neq_smc + grid->num_coupling_species_smc) * sizeof(int),
			"allocation failed for header_offset array in dump_smc_data.");
	smcDataOffset = (int*) checked_malloc((grid->neq_smc + grid->num_coupling_species_smc) * sizeof(int),
			"allocation failed for smcDataOffset array in dump_smc_data.");
	for (int i = 0; i < (3 + grid->neq_smc + grid->num_coupling_species_smc); i++) {
		header_offset[i] = 0;
	}
	for (int i = 0; i < (grid->neq_smc + grid->num_coupling_species_smc); i++) {
		smcDataOffset[i] = 0;
	}

	int count = 0;
	char* header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	int branches;
	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		branches = 1;
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		branches = 3;
	}

	/*************** Writing VTK header **************/
	header_offset[0] = sprintf(header, "# vtk DataFile Version 2.0\n"
			"Time file at t = %d seconds\n"
			"ASCII\n"
			"DATASET UNSTRUCTURED_GRID\n"
			"POINTS %d double\n", write_count, grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * 4 * branches);
	count = header_offset[0];
	disp = 0;
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status),
				"error writing into time file by writer_rank 0.\n");
	}
	/*************** Writing Point data **************/
	int buffer_lengths[4] = { 0, 0, 0, 0 };
	check_flag(MPI_Allgather(&writer_buffer->buffer_length[smcMesh], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths.\n");

	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}

	disp = (header_offset[0] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[smcMesh];
	check_flag(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->smc_mesh_points, count, MPI_CHAR, &status),
			"error writing the smc coordinates in file.\n");

	for (int i = 0; i < 4; i++) {
		point_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->smc_mesh_points);

	/*************** Writing cell data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[1] = sprintf(header, "CELLS %d %d\n", (grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches),
			5 * (grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches));

	count = header_offset[1];
	disp = (header_offset[0] + point_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}

	check_flag(MPI_Allgather(&writer_buffer->buffer_length[smcCell], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + disp + header_offset[1]) * sizeof(char);
	count = writer_buffer->buffer_length[smcCell];
	check_flag(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->smc_mesh_cells, count, MPI_CHAR, &status),
			"error writing the smc cells in file.");
	for (int i = 0; i < 4; i++) {
		cell_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->smc_mesh_cells);

	/*************** Writing cell type data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[2] = sprintf(header, "CELL_TYPES %d\n", (grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches));

	count = header_offset[2];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(MPI_Allgather(&writer_buffer->buffer_length[smcCellType], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths");
	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[smcCellType];
	check_flag(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->smc_mesh_type, count, MPI_CHAR, &status),
			"error writing the cell_type in file.");
	for (int i = 0; i < 4; i++) {
		celltype_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->smc_mesh_type);

	/***************************************************************/
	/********		Writing Field 1 : SMC Ca data 			********/
	/***************************************************************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[3] = sprintf(header, "CELL_DATA %d\nFIELD SMC_Data %d\n"
			"SMC_Ca %d %d float\n", grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches,
			grid->neq_smc + grid->num_coupling_species_smc, 1, grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches);

	count = header_offset[3];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(
			MPI_Allgather(&writer_buffer->smc_stat_var_buffer_length[smc_Ca], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT,
					my_IO_domain_info->writer_comm), "error in all gather called for smc Data buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + disp)
			* sizeof(char);
	count = writer_buffer->smc_stat_var_buffer_length[smc_Ca];
	check_flag(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->ci, count, MPI_CHAR, &status),
			"error writing the smc calcium in file.\n");
	for (int i = 0; i < 4; i++) {
		smcDataOffset[0] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ci);

	/*************** Writing Field 2 : SMC SR data ***************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[4] = sprintf(header, "SMC_SR %d %d float\n", 1, grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches);

	count = header_offset[4];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}

	check_flag(
			MPI_Allgather(&writer_buffer->smc_stat_var_buffer_length[smc_SR], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT,
					my_IO_domain_info->writer_comm), "error in all gather called for smc Data buffer lengths\n");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + disp) * sizeof(char);
	count = writer_buffer->smc_stat_var_buffer_length[smc_SR];
	check_flag(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->si, count, MPI_CHAR, &status), "error writing the smc SR Ca in file.\n");
	for (int i = 0; i < 4; i++) {
		smcDataOffset[1] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->si);

	/*************** Writing Field 3 : SMC Vm data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[5] = sprintf(header, "SMC_Vm %d %d float\n", 1, grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches);

	count = header_offset[5];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(
			MPI_Allgather(&writer_buffer->smc_stat_var_buffer_length[smc_Vm], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT,
					my_IO_domain_info->writer_comm), "error in all gather called for smc Data buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + disp) * sizeof(char);
	count = writer_buffer->smc_stat_var_buffer_length[smc_Vm];
	check_flag(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->vi, count, MPI_CHAR, &status), "error writing the smc Vm in file.\n");
	for (int i = 0; i < 4; i++) {
		smcDataOffset[2] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->vi);
	/*************** Writing Field 4 : SMC w data ***************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[6] = sprintf(header, "SMC_w %d %d float\n", 1, grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches);

	count = header_offset[6];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(
			MPI_Allgather(&writer_buffer->smc_stat_var_buffer_length[smc_w], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT,
					my_IO_domain_info->writer_comm), "error in all gather called for smc Data buffer lengths.\n");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + disp) * sizeof(char);
	count = writer_buffer->smc_stat_var_buffer_length[smc_w];
	check_flag(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->wi, count, MPI_CHAR, &status), "error writing the smc KCa in file.\n");
	for (int i = 0; i < 4; i++) {
		smcDataOffset[3] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->wi);
	/*************** Writing Field 5 : SMC I data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0.\n");
	header_offset[7] = sprintf(header, "SMC_IP3 %d %d float\n", 1,
			grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches);

	count = header_offset[7];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3])
			* sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(
			MPI_Allgather(&writer_buffer->smc_stat_var_buffer_length[smc_IP3], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT,
					my_IO_domain_info->writer_comm), "error in all gather called for smc Data buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + disp) * sizeof(char);
	count = writer_buffer->smc_stat_var_buffer_length[smc_IP3];
	check_flag(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->Ii, count, MPI_CHAR, &status), "error writing the smc IP3 in file.\n");
	for (int i = 0; i < 4; i++) {
		smcDataOffset[4] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->Ii);
	/*************** Writing Field 6 : SMC Ca coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[8] = sprintf(header, "SMC_Ca_coupling %d %d float\n", 1,
			grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches);

	count = header_offset[8];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(MPI_Allgather(&writer_buffer->smc_cpl[cpl_Ca], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for smc Ca cpl Data buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4] + header_offset[8] + disp) * sizeof(char);
	count = writer_buffer->smc_cpl[cpl_Ca];
	check_flag(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->cpCi, count, MPI_CHAR, &status),
			"error writing the smc Ca coupling in file.\n");
	for (int i = 0; i < 4; i++) {
		smcDataOffset[5] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpCi);
	/*************** Writing Field 7 : SMC Vm coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[9] = sprintf(header, "SMC_Vm_coupling %d %d float\n", 1,
			grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches);

	count = header_offset[9];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4] + header_offset[8] + smcDataOffset[5]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(MPI_Allgather(&writer_buffer->smc_cpl[cpl_Vm], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for smc Vm cpl Data buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4] + header_offset[8] + smcDataOffset[5] + header_offset[9] + disp) * sizeof(char);
	count = writer_buffer->smc_cpl[cpl_Vm];
	check_flag(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->cpVi, count, MPI_CHAR, &status),
			"error writing the smc Vm coupling in file.\n");
	for (int i = 0; i < 4; i++) {
		smcDataOffset[6] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpVi);

	/*************** Writing Field 8 : SMC IP3 coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[10] = sprintf(header, "SMC_IP3_coupling %d %d float\n", 1,
			grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * branches);

	count = header_offset[10];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4] + header_offset[8] + smcDataOffset[5] + header_offset[9] + smcDataOffset[6]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(MPI_Allgather(&writer_buffer->smc_cpl[cpl_IP3], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for smc IP3 cpl Data buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4] + header_offset[8] + smcDataOffset[5] + header_offset[9] + smcDataOffset[6] + header_offset[10]
			+ disp) * sizeof(char);
	count = writer_buffer->smc_cpl[cpl_IP3];
	check_flag(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->cpIi, count, MPI_CHAR, &status),
			"error writing the smc IP3 coupling in file.\n");
	for (int i = 0; i < 4; i++) {
		smcDataOffset[5] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpIi);
}

/**************************************************************************************/
void dump_agonists_map(checkpoint_handle* check, grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, celltype2** ec,
		char* path) {
	MPI_Status status;
	MPI_Offset disp;
	char filename[50];
	int err = sprintf(filename, "%s/Agonist_map.vtk", path);
	CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->ec_agonist_file));

	int write_element_count = 1, *header_offset, point_offset = 0, cell_offset = 0, celltype_offset = 0, *ecDataOffset;
	header_offset = (int*) checked_malloc((3 + grid->num_parameters) * sizeof(int),
			"allocation failed for header_offset array in dump_agonists_map.");
	ecDataOffset = (int*) checked_malloc((grid->num_parameters) * sizeof(int), "allocation failed for ecDataOffset array in dump_agonists_map.");
	for (int i = 0; i < (3 + grid->num_parameters); i++) {
		header_offset[i] = 0;
	}
	for (int i = 0; i < (grid->num_parameters); i++) {
		ecDataOffset[i] = 0;
	}

	int count = 0;
	char* header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	int branches;
	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		branches = 1;
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		branches = 3;
	}

	/*************** Writing VTK header **************/
	header_offset[0] = sprintf(header, "# vtk DataFile Version 2.0\n"
			"Agonist Map file\n"
			"ASCII\n"
			"DATASET UNSTRUCTURED_GRID\n"
			"POINTS %d double\n", grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * 4 * branches);
	count = header_offset[0];
	disp = 0;
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_agonist_file, disp, header, count, MPI_CHAR, &status),
				"error writing into time file by writer_rank 0.\n");
	}
	/*************** Writing Point data **************/
	int buffer_lengths[4] = { 0, 0, 0, 0 };
	check_flag(MPI_Allgather(&writer_buffer->buffer_length[ecMesh], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths\n");

	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}

	disp = (header_offset[0] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[ecMesh];
	check_flag(MPI_File_write_at(check->ec_agonist_file, disp, writer_buffer->ec_mesh_points, count, MPI_CHAR, &status),
			"error writing the ec coordinates in file.\n");

	for (int i = 0; i < 4; i++) {
		point_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ec_mesh_points);

	/*************** Writing cell data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[1] = sprintf(header, "CELLS %d %d\n", (grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches),
			5 * (grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches));

	count = header_offset[1];
	disp = (header_offset[0] + point_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_agonist_file, disp, header, count, MPI_CHAR, &status),
				"error writing into time file by writer_rank 0.\n");
	}

	check_flag(MPI_Allgather(&writer_buffer->buffer_length[ecCell], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths\n");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + disp + header_offset[1]) * sizeof(char);
	count = writer_buffer->buffer_length[ecCell];
	check_flag(MPI_File_write_at(check->ec_agonist_file, disp, writer_buffer->ec_mesh_cells, count, MPI_CHAR, &status),
			"error writing the ec cells in ec_agonist_file.\n");
	for (int i = 0; i < 4; i++) {
		cell_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ec_mesh_cells);

	/*************** Writing cell type data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[2] = sprintf(header, "CELL_TYPES %d\n", (grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches));

	count = header_offset[2];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_agonist_file, disp, header, count, MPI_CHAR, &status),
				"error writing into time file by writer_rank 0.\n");
	}
	check_flag(MPI_Allgather(&writer_buffer->buffer_length[ecCellType], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths.\n");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[ecCellType];
	check_flag(MPI_File_write_at(check->ec_agonist_file, disp, writer_buffer->ec_mesh_type, count, MPI_CHAR, &status),
			"error writing the ec cell type in ec_agonist_file.\n");
	for (int i = 0; i < 4; i++) {
		celltype_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ec_mesh_type);

	/***************************************************************/
	/********		Writing Field 1 : JPLC Data 			********/
	/***************************************************************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[3] = sprintf(header, "CELL_DATA %d\nFIELD ec_Data %d\n"
			"JPLC %d %d float\n", grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches,
	/*grid->num_parameters*/1, 1, grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches);

	count = header_offset[3];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_agonist_file, disp, header, count, MPI_CHAR, &status),
				"error writing into agonist file by writer_rank 0.\n");
	}
	check_flag(MPI_Allgather(&writer_buffer->jplc_buffer_length, 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for ec Data buffer lengths.\n");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + disp)
			* sizeof(char);
	count = writer_buffer->jplc_buffer_length;
	check_flag(MPI_File_write_at(check->ec_agonist_file, disp, writer_buffer->jplc, count, MPI_CHAR, &status),
			"error writing the ec JPLC in ec_agonist_file.\n");
	for (int i = 0; i < 4; i++) {
		ecDataOffset[0] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->jplc);

	MPI_File_close(&check->ec_agonist_file);
}
/*************************************************************************************/
/**************************************************************************************/
void dump_ec_data(checkpoint_handle* check, grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, celltype2** ec,
		int write_count) {
	MPI_Status status;
	MPI_Offset disp;

	int write_element_count = 1, *header_offset, point_offset = 0, cell_offset = 0, celltype_offset = 0, *ecDataOffset;
	header_offset = (int*) checked_malloc((3 + grid->neq_ec + grid->num_coupling_species_ec) * sizeof(int),
			"allocation failed for header_offset array in dump_ec_data.");
	ecDataOffset = (int*) checked_malloc((grid->neq_ec + grid->num_coupling_species_ec) * sizeof(int),
			"allocation failed for ecDataOffset array in dump_ec_data.\n");
	for (int i = 0; i < (3 + grid->neq_ec + grid->num_coupling_species_ec); i++) {
		header_offset[i] = 0;
	}
	for (int i = 0; i < (grid->neq_ec + grid->num_coupling_species_ec); i++) {
		ecDataOffset[i] = 0;
	}

	int count = 0;
	char* header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	int branches;
	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		branches = 1;
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		branches = 3;
	}

	/*************** Writing VTK header **************/
	header_offset[0] = sprintf(header, "# vtk DataFile Version 2.0\n"
			"Time file at t = %d seconds\n"
			"ASCII\n"
			"DATASET UNSTRUCTURED_GRID\n"
			"POINTS %d double\n", write_count, grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * 4 * branches);
	count = header_offset[0];
	disp = 0;
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	/*************** Writing Point data **************/
	int buffer_lengths[4] = { 0, 0, 0, 0 };
	check_flag(MPI_Allgather(&writer_buffer->buffer_length[ecMesh], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths");

	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}

	disp = (header_offset[0] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[ecMesh];
	check_flag(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->ec_mesh_points, count, MPI_CHAR, &status),
			"error writing the ec coordinates in file.\n");

	for (int i = 0; i < 4; i++) {
		point_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ec_mesh_points);

	/*************** Writing cell data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[1] = sprintf(header, "CELLS %d %d\n", (grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches),
			5 * (grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches));

	count = header_offset[1];
	disp = (header_offset[0] + point_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}

	check_flag(MPI_Allgather(&writer_buffer->buffer_length[ecCell], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths.\n");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + disp + header_offset[1]) * sizeof(char);
	count = writer_buffer->buffer_length[ecCell];
	check_flag(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->ec_mesh_cells, count, MPI_CHAR, &status),
			"error writing the ec cells in file.\n");
	for (int i = 0; i < 4; i++) {
		cell_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ec_mesh_cells);

	/*************** Writing cell type data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[2] = sprintf(header, "CELL_TYPES %d\n", (grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches));

	count = header_offset[2];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(MPI_Allgather(&writer_buffer->buffer_length[ecCellType], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[ecCellType];
	check_flag(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->ec_mesh_type, count, MPI_CHAR, &status),
			"error writing the cell type in file.\n");
	for (int i = 0; i < 4; i++) {
		celltype_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ec_mesh_type);

	/***************************************************************/
	/********		Writing Field 1 : ec Ca data 			********/
	/***************************************************************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[3] = sprintf(header, "CELL_DATA %d\nFIELD ec_Data %d\n"
			"ec_Ca %d %d float\n", grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches,
			grid->neq_ec + grid->num_coupling_species_ec, 1, grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches);

	count = header_offset[3];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(
			MPI_Allgather(&writer_buffer->ec_stat_var_buffer_length[ec_Ca], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT,
					my_IO_domain_info->writer_comm), "error in all gather called for ec Data buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + disp)
			* sizeof(char);
	count = writer_buffer->ec_stat_var_buffer_length[ec_Ca];
	check_flag(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->cj, count, MPI_CHAR, &status), "error writing the ec Ca in file.\n");
	for (int i = 0; i < 4; i++) {
		ecDataOffset[0] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cj);

	/*************** Writing Field 2 : ec SR data ***************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[4] = sprintf(header, "ec_SR %d %d float\n", 1, grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches);

	count = header_offset[4];
	disp =
			(header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
					+ ecDataOffset[0]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}

	check_flag(
			MPI_Allgather(&writer_buffer->ec_stat_var_buffer_length[ec_SR], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT,
					my_IO_domain_info->writer_comm), "error in all gather called for ec Data buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + disp) * sizeof(char);
	count = writer_buffer->ec_stat_var_buffer_length[ec_SR];
	check_flag(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->sj, count, MPI_CHAR, &status), "error writing the ec SR Ca in file.\n");
	for (int i = 0; i < 4; i++) {
		ecDataOffset[1] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->sj);

	/*************** Writing Field 3 : ec Vm data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[5] = sprintf(header, "ec_Vm %d %d float\n", 1, grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches);

	count = header_offset[5];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(
			MPI_Allgather(&writer_buffer->ec_stat_var_buffer_length[ec_Vm], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT,
					my_IO_domain_info->writer_comm), "error in all gather called for ec Data buffer lengths.\n");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + disp) * sizeof(char);
	count = writer_buffer->ec_stat_var_buffer_length[ec_Vm];
	check_flag(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->vj, count, MPI_CHAR, &status), "error writing the ec Vm  in file.\n");
	for (int i = 0; i < 4; i++) {
		ecDataOffset[2] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->vj);
	/*************** Writing Field 4 : ec I data ***************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[6] = sprintf(header, "ec_IP3 %d %d float\n", 1, grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches);

	count = header_offset[6];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(
			MPI_Allgather(&writer_buffer->ec_stat_var_buffer_length[ec_IP3], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT,
					my_IO_domain_info->writer_comm), "error in all gather called for ec Data buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + disp) * sizeof(char);
	count = writer_buffer->ec_stat_var_buffer_length[ec_IP3];
	check_flag(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->Ij, count, MPI_CHAR, &status), "error writing the ec IP3 in file.\n");
	for (int i = 0; i < 4; i++) {
		ecDataOffset[3] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->Ij);
	/*************** Writing Field 5 : ec Ca coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[7] = sprintf(header, "ec_Ca_coupling %d %d float\n", 1,
			grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches);

	count = header_offset[7];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(MPI_Allgather(&writer_buffer->ec_cpl[cpl_Ca], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for ec Data buffer lengths");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3] + header_offset[7] + disp)
			* sizeof(char);
	count = writer_buffer->ec_cpl[cpl_Ca];
	check_flag(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->cpCj, count, MPI_CHAR, &status),
			"error writing the ec Ca coupling in file.\n");
	for (int i = 0; i < 4; i++) {
		ecDataOffset[4] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpCj);
	/*************** Writing Field 6 : ec Vm coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0.\n");
	header_offset[8] = sprintf(header, "ec_Vm_coupling %d %d float\n", 1,
			grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches);

	count = header_offset[8];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3] + header_offset[7]
			+ ecDataOffset[4]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status),
				"error writing into time file by writer_rank 0.\n");
	}
	check_flag(MPI_Allgather(&writer_buffer->ec_cpl[cpl_Vm], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for ec Vm cpl Data buffer lengths.\n");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3] + header_offset[7]
			+ ecDataOffset[4] + header_offset[8] + disp) * sizeof(char);
	count = writer_buffer->ec_cpl[cpl_Vm];
	check_flag(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->cpVj, count, MPI_CHAR, &status),
			"error writing the ec Vm coupling in file.\n");
	for (int i = 0; i < 4; i++) {
		ecDataOffset[5] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpVj);
	/*************** Writing Field 7 : ec IP3 coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[9] = sprintf(header, "ec_IP3_coupling %d %d float\n", 1,
			grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * branches);

	count = header_offset[9];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3] + header_offset[7]
			+ ecDataOffset[4] + header_offset[8] + ecDataOffset[5]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}
	check_flag(MPI_Allgather(&writer_buffer->ec_cpl[cpl_IP3], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for ec Vm cpl Data buffer lengths.\n");

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3] + header_offset[7]
			+ ecDataOffset[4] + header_offset[8] + ecDataOffset[5] + header_offset[9] + disp) * sizeof(char);
	count = writer_buffer->ec_cpl[cpl_IP3];
	check_flag(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->cpIj, count, MPI_CHAR, &status),
			"error writing the ec IP3 coupling in file.\n");
	for (int i = 0; i < 4; i++) {
		ecDataOffset[6] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpIj);
}

/*************************************************************************************/
void update_line_number(checkpoint_handle* check, grid_parms grid, int line_number) {
	MPI_Status status;
	MPI_Offset disp;

	disp = grid.universal_rank * sizeof(int);
	CHECK(MPI_File_write_at(check->line_number, disp, &line_number, 1, MPI_INT, &status));
}
/****************************************************************************************/
void dump_rank_info(checkpoint_handle* check, conductance cpl_cef, grid_parms grid, IO_domain_info* my_IO_domain_info) {
	MPI_Status status;
	MPI_Offset displacement = 0;
	char* buffer = (char*) checked_malloc(2 * 1024 * sizeof(char), "allocation for logfile segment space\n");
	int root = 0;
	char filename[50];
	int length =
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
					grid.my_domain.internal_info.boundary_tag, cpl_cef.Vm_hm_smc, cpl_cef.Vm_hm_ec, cpl_cef.Ca_hm_smc, cpl_cef.Ca_hm_ec,
					cpl_cef.IP3_hm_smc, cpl_cef.IP3_hm_ec, cpl_cef.Vm_ht_smc, cpl_cef.Vm_ht_ec, cpl_cef.Ca_ht_smc, cpl_cef.Ca_ht_ec,
					cpl_cef.IP3_ht_smc, cpl_cef.IP3_ht_ec, grid.uniform_jplc, grid.min_jplc, grid.max_jplc, grid.gradient, grid.numtasks, grid.m,
					grid.n, grid.num_ec_axially, grid.num_smc_circumferentially, grid.num_ec_axially * grid.num_ec_circumferentially,
					grid.num_smc_axially * grid.num_smc_circumferentially,
					(grid.num_ec_axially * grid.num_ec_circumferentially) + (grid.num_smc_axially * grid.num_smc_circumferentially),
					(grid.num_ec_circumferentially * grid.num_ec_axially * grid.numtasks),
					(grid.num_smc_circumferentially * grid.num_smc_axially * grid.numtasks),
					((grid.num_ec_axially * grid.num_ec_circumferentially) + (grid.num_smc_axially * grid.num_smc_circumferentially)) * grid.numtasks,
					grid.NEQ * grid.numtasks, grid.my_domain.z_offset_start, grid.my_domain.z_offset_end, grid.my_domain.local_z_start,
					grid.my_domain.local_z_end);

	int *recv_count = (int*) checked_malloc(grid.tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	int *disp = (int*) checked_malloc(grid.tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	/// Gathering and summing the length of all the CHARs contained in every send_buffer containing coordinates from each MPI process.
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid.cart_comm),
			"error in MPI_Gather gathering Logfile buffer length by each process.");
	grid.logfile_displacements = 0;
	for (int i = 0; i < grid.tasks; i++) {
		disp[i] = grid.logfile_displacements;
		grid.logfile_displacements += recv_count[i];
	}

	if (grid.rank == 0) {
		grid.logfile_write_buffer = (char*) checked_malloc(grid.logfile_displacements * sizeof(char),
				"allocation error for writer_buffer for Logfile.");
	}
	check_flag(MPI_Gatherv(buffer, length, MPI_CHAR, grid.logfile_write_buffer, recv_count, disp, MPI_CHAR, root, grid.cart_comm),
			"Error gathering Logfile data.");

	if (grid.universal_rank == 0) {
		int dir_status = mkdir("Logfiles", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);

		printf("dir status return =  %d\n", dir_status);
	}
	if (grid.rank == 0) {
		sprintf(filename, "Logfiles/Logfile_%s.txt", grid.suffix);
		MPI_Barrier(grid.universe);
		check_flag(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->logptr), "error opening Logfile");
		check_flag(MPI_File_write_at(check->logptr, displacement, grid.logfile_write_buffer, grid.logfile_displacements, MPI_CHAR, &status),
				"Error writing data into Logfile.");
		MPI_File_close(&check->logptr);
		free(grid.logfile_write_buffer);
	}

	free(recv_count);
	free(disp);

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
	CHECK(MPI_File_write_at(check->jplc, disp, &buffer, write_element_count, MPI_DOUBLE, &status));
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
	CHECK(MPI_File_write_at(check->coords, disp, &buffer, write_element_count, MPI_DOUBLE, &status));
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

	CHECK(MPI_File_write_at_all(check->time_profiling, disp_write, &buffer[0], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->async_calls, disp_write, &buffer[1], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->async_wait, disp_write, &buffer[2], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->barrier_before_comm, disp_write, &buffer[3], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->map_function, disp_write, &buffer[4], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->single_cell_fluxes, disp_write, &buffer[5], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->coupling_fluxes, disp_write, &buffer[6], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->solver, disp_write, &buffer[7], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->writer_func, disp_write, &buffer[8], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->derivative_calls, disp_write, &buffer[9], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->itter_count, disp_write, &buffer[10], 1, MPI_DOUBLE, &status));
/// Write Comms time profiling data...
	CHECK(MPI_File_write_at_all(check->remote_async_calls, disp_write, &buffer[11], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->remote_async_wait, disp_write, &buffer[12], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->send_buf_update, disp_write, &buffer[13], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->recv_buf_update, disp_write, &buffer[14], 1, MPI_DOUBLE, &status));
	CHECK(MPI_File_write_at_all(check->total_comms_cost, disp_write, &buffer[15], 1, MPI_DOUBLE, &status));

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
	MPI_File_close(&check->smc_data_file);
	MPI_File_close(&check->ec_data_file);

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
	CHECK(MPI_File_open(grid->universe, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &data));

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
/// Element 4:	Parent branch index of the parent subdomain.
/// Element 5:  Subdomain index of the parent subdomain in the Parent branch.
/// Element 6:  Global subdomain index of the parent subdomain.

/// Element 7:	Left child branch index of the parent subdomain.
/// Element 8:  Subdomain index of the left child subdomain in the Left child branch.
/// Element 9:  Global subdomain index of the left child subdomain.

/// Element 10:	Right child branch index of the parent subdomain.
/// Element 11: Subdomain index of the Right child subdomain in the Right child branch.
/// Element 12: Global subdomain index of the Right child subdomain.

/// Element 13: Requested number of ECs per mesh element in the subdomain
/// Element 14: Requested number of SMCs per mesh element in the subdomain

/// Element 15: My own branch index
/// Element 16: My own local subdomain index of the branch I belong to.

/// Element 17 - 19: Number of points to be read from each type of grid, i.e. Processor Mesh, SMC Mesh and EC Mesh respectively, for my subdomain.
/// Element 20 - 22: Number of cells to be read from each type of grid, i.e. Processor Mesh, SMC Mesh and EC Mesh respectively, for my subdomain. The last entry i.e. index 22 denotes bot number of EC quad and EC Centroids to be read.

/// In the case of elements (7,8,9) & (10,11,12), if subdomain type of current key_val is a straight segment, left Child is positive or zero, and right Child is negative.
/// If subdomain type of current key_val is a bifurcation, then both right and left child subdomains are non-negative.
////
	for (int i = 0; i < grid->num_domains; i++) {
		grid->domains[i] = (int*) checked_malloc(23 * sizeof(int), "Subdomains array elements allocation");
	}
	for (int i = 0; i < grid->num_domains; i++) {
		for (int j = 0; j < 17; j++) {
			grid->domains[i][j] = p[1 + (i * 23) + j];
		}
	}
	MPI_File_close(&data);
	return (0);
}

void update_elapsed_time(checkpoint_handle* check, grid_parms grid, time_keeper* elps_t, IO_domain_info* my_IO_domain_info) {

	MPI_Status status;
	char filename[50];
	int root = 0;
	char *buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * sizeof(char), "Allocation error in buffer for elapsed time data.\n");
	char* write_buffer;
	elps_t->t_new = MPI_Wtime();
	elps_t->elapsed_time = elps_t->t_new - elps_t->t_old;
	elps_t->t_old = elps_t->t_new;

	int length = sprintf(buffer, "%2.12lf\n", elps_t->elapsed_time);

	int *recv_count = (int*) checked_malloc(grid.tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	int *disp = (int*) checked_malloc(grid.tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	/// Gathering and summing the length of all the CHARs contained in every send_buffer containing coordinates from each MPI process.
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid.cart_comm),
			"error in MPI_Gather gathering Logfile buffer length by each process.");
	int total_buffer_length = 0;
	for (int i = 0; i < grid.tasks; i++) {
		disp[i] = total_buffer_length;
		total_buffer_length += recv_count[i];
	}
	if (grid.rank == 0) {
		write_buffer = (char*) checked_malloc(total_buffer_length * sizeof(char), "allocation error for writer_buffer for Elapsed_time file.");
	}
	check_flag(MPI_Gatherv(buffer, length, MPI_CHAR, write_buffer, recv_count, disp, MPI_CHAR, root, grid.cart_comm),
			"Error gathering Logfile data.");
	if (grid.rank == 0) {
		sprintf(filename, "Elapsed_time_%s.txt", grid.suffix);
		check_flag(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->elapsed_time),
				"error opening Logfile");
		check_flag(MPI_File_write_at(check->elapsed_time, 0, write_buffer, total_buffer_length, MPI_CHAR, &status),
				"Error writing data into Elapsed Time file.");
		MPI_File_close(&check->elapsed_time);
		free(write_buffer);
	}

	free(recv_count);

	free(disp);

}

void naming_convention(grid_parms* grid) {
//Prepare the suffix which indicates my subdomain information and if I am a bifurcation, then also tells about which branch do I belong to
	int subdomain, branch, err;

	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		err = sprintf(grid->suffix, "%d_%d", grid->my_domain.internal_info.domain_branch, grid->my_domain.internal_info.domain_local_subdomain);

	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		subdomain = grid->my_domain.internal_info.domain_index;
		branch = grid->branch_tag;
		err = sprintf(grid->suffix, "%d_%d_%d", grid->my_domain.internal_info.domain_branch, grid->my_domain.internal_info.domain_local_subdomain,
				grid->branch_tag);
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
	CHECK(MPI_File_read_at(check->itter_count, disp, &file_offset, 1, MPI_INT, &status));

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
		CHECK(MPI_File_read_at(check->coords, disp, &z, 1, MPI_DOUBLE, &status));
		CHECK(MPI_File_read_at(check->jplc, disp, &jplc, 1, MPI_DOUBLE, &status));
		fprintf(fh, "%2.8lf\t%2.8lf\n", z, jplc);
	}
	fclose(fh);

}

/*********************************************************************************************************/
checkpoint_handle* initialise_time_wise_checkpoint(checkpoint_handle* check, grid_parms grid, int write_count, char* path,
		IO_domain_info* my_IO_domain_info)
		/*********************************************************************************************************/
		{
	int err;
	char filename[50];
	open_koenigsberger_smc_checkpoint(check, grid, write_count, path, my_IO_domain_info);
	open_koenigsberger_ec_checkpoint(check, grid, write_count, path, my_IO_domain_info);
	//open_coupling_data_checkpoint(check, grid, write_count, path, my_IO_domain_info);
	return (check);
}

//Read coordinates from relevant geometry files for each ec and smc in the computational domain.
/************************************************************************************************/
int retrieve_topology_info(char* filename, grid_parms* grid, celltype1 **smc, celltype2 **ec)
/************************************************************************************************/
{
	/** This is a function to read the points data from the .txt files that are the output from the
	 * Surface mesh generator code. These points are effectively the vertices of the quadrilaterals
	 * representing either a processor/ MPI-Task, an SMC or an EC. The four vertices are in a set of
	 * four consecutive points (each with three components x,y & z) in the file.
	 */

	/// info is a 2d array in the structure grid which will have the information on the amount of data to be read by the reading processor and to decide which data belongs to which processor/MPI-task of the MPI_COMM_WORLD
	grid->info = (int**) checked_malloc(4 * sizeof(int*), "memory allocation for info failed"); /// 4 levels of information i.e. MPI-Tasks(ProcessMesh), SMC grid (smcMesh), EC grid (ecMesh), EC centroids(ecCentroids).
	for (int i = 0; i < 4; i++) {
		grid->info[i] = (int*) checked_malloc(2 * sizeof(int), "memory allocation for info failed d2"); /// Each information level has the information for Points and Cells (total, in axial direction, in circumferential direction) of that respective category.
	}

	grid->info[ProcessMesh][0] = grid->domains[grid->my_domain.internal_info.domain_index][17]; /// Record total number of points for processor/MPI-task grid to be read
	grid->info[ProcessMesh][1] = grid->domains[grid->my_domain.internal_info.domain_index][20]; /// Record total number of cells for processor/MPI-task grid to be read
	grid->info[smcMesh][0] = grid->domains[grid->my_domain.internal_info.domain_index][18]; /// Record total number of points for smcMesh grid to be read
	grid->info[smcMesh][1] = grid->domains[grid->my_domain.internal_info.domain_index][21]; /// Record total number of cells for smcMesh grid to be read
	grid->info[ecMesh][0] = grid->domains[grid->my_domain.internal_info.domain_index][19]; /// Record total number of points for ecMesh grid to be read
	grid->info[ecMesh][1] = grid->domains[grid->my_domain.internal_info.domain_index][22]; /// Record total number of cells for ecMesh grid to be read
	grid->info[ecCentroids][0] = grid->domains[grid->my_domain.internal_info.domain_index][22]; /// Record total number of points for ecMesh grid to be read
	grid->info[ecCentroids][1] = grid->domains[grid->my_domain.internal_info.domain_index][22]; /// Record total number of cells for ecMesh grid to be read

	int *disp, *send_count, recv_count, root = 0;
	double *send_points, *recv_points;
	int* indx = (int*) malloc(2 * sizeof(int));
	int num_tuple_components, num_tuples;
	int tuple_offset;

	/// Setting up file reading environment for reading in the ProcessMesh vertices from the process mesh.
	num_tuple_components = 3;
	num_tuples = 4;
	tuple_offset = num_tuple_components * num_tuples;

	/// coordinates is the array to store four vertices of the quad representing me,a processors, in the subdomain I am a member of.
	grid->coordinates = (double**) checked_malloc(num_tuples * sizeof(double*), "allocation error in D1 of grid coordinates.");
	for (int i = 0; i < num_tuples; i++) {
		grid->coordinates[i] = (double*) checked_malloc(num_tuple_components * sizeof(double), "allocation error in D2 of grid coordinates.");
	}

	send_count = (int*) malloc(grid->info[ProcessMesh][TOTAL_CELLS] * sizeof(int));
	disp = (int*) malloc(grid->info[ProcessMesh][TOTAL_CELLS] * sizeof(int));

	send_points = (double*) malloc(tuple_offset * grid->info[ProcessMesh][TOTAL_CELLS] * sizeof(double));

	/// Only Rank 0 of each subdomain is going to read the coordinates from the files. These coordinates will then be sent
	/// the other member processors of that subdomain, again to minimize the expense of reading operation on filesystem.
	if (grid->rank == 0) {
		vtk_info* process_mesh = (vtk_info*) malloc(sizeof(vtk_info));
		process_mesh->points = (double**) malloc(grid->info[ProcessMesh][TOTAL_POINTS] * sizeof(double*));
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_POINTS]; i++) {
			process_mesh->points[i] = (double*) malloc(3 * sizeof(double));
		}
		process_mesh->cells = (int**) malloc(grid->info[ProcessMesh][TOTAL_CELLS] * sizeof(int*));
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS]; i++) {
			process_mesh->cells[i] = (int*) malloc(5 * sizeof(int));
		}

		indx = read_coordinates(grid, grid->info, process_mesh, ProcessMesh, grid->info[ProcessMesh][TOTAL_POINTS],
				grid->info[ProcessMesh][TOTAL_CELLS]);
		for (int i = 0; i < indx[1]; i++) {
			send_points[(i * tuple_offset) + 0] = process_mesh->points[process_mesh->cells[i][1]][0];
			send_points[(i * tuple_offset) + 1] = process_mesh->points[process_mesh->cells[i][1]][1];
			send_points[(i * tuple_offset) + 2] = process_mesh->points[process_mesh->cells[i][1]][2];

			send_points[(i * tuple_offset) + 3] = process_mesh->points[process_mesh->cells[i][2]][0];
			send_points[(i * tuple_offset) + 4] = process_mesh->points[process_mesh->cells[i][2]][1];
			send_points[(i * tuple_offset) + 5] = process_mesh->points[process_mesh->cells[i][2]][2];

			send_points[(i * tuple_offset) + 6] = process_mesh->points[process_mesh->cells[i][3]][0];
			send_points[(i * tuple_offset) + 7] = process_mesh->points[process_mesh->cells[i][3]][1];
			send_points[(i * tuple_offset) + 8] = process_mesh->points[process_mesh->cells[i][3]][2];

			send_points[(i * tuple_offset) + 9] = process_mesh->points[process_mesh->cells[i][4]][0];
			send_points[(i * tuple_offset) + 10] = process_mesh->points[process_mesh->cells[i][4]][1];
			send_points[(i * tuple_offset) + 11] = process_mesh->points[process_mesh->cells[i][4]][2];
		}

		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_POINTS]; i++) {
			free(process_mesh->points[i]);
		}
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS]; i++) {
			free(process_mesh->cells[i]);
		}
		free(process_mesh->points);
		free(process_mesh->cells);
		free(process_mesh);
	}
	recv_points = (double*) malloc(tuple_offset * sizeof(double));

	for (int i = 0; i < grid->tasks; i++) {
		send_count[i] = tuple_offset;
		disp[i] = tuple_offset * i;
	}
	recv_count = tuple_offset;
	check_flag(MPI_Scatterv(send_points, send_count, disp, MPI_DOUBLE, recv_points, recv_count, MPI_DOUBLE, root, grid->cart_comm),
			"error in scatter of coordinates");
	for (int i = 0; i < num_tuples; i++) {
		grid->coordinates[i][0] = recv_points[(i * num_tuple_components) + 0];
		grid->coordinates[i][1] = recv_points[(i * num_tuple_components) + 1];
		grid->coordinates[i][2] = recv_points[(i * num_tuple_components) + 2];
	}
	free(send_points);
	free(recv_points);

//Reading in SMC mesh and communicating to each branch communicator memeber.
	send_points = (double*) malloc(tuple_offset * grid->info[smcMesh][TOTAL_CELLS] * grid->info[ProcessMesh][TOTAL_CELLS] * sizeof(double));

	if (grid->rank == 0) {

		vtk_info* smc_mesh = (vtk_info*) malloc(sizeof(vtk_info));
		smc_mesh->points = (double**) malloc(grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_POINTS] * sizeof(double*));
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_POINTS]; i++) {
			smc_mesh->points[i] = (double*) malloc(3 * sizeof(double));
		}
		smc_mesh->cells = (int**) malloc(grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS] * sizeof(int*));
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS]; i++) {
			smc_mesh->cells[i] = (int*) malloc(5 * sizeof(int));
		}

		indx = read_coordinates(grid, grid->info, smc_mesh, smcMesh, (grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_POINTS]),
				(grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS]));

		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS]; i++) {
			for (int j = 0; j < grid->info[smcMesh][TOTAL_CELLS]; j++) {
				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 0] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][1]][0];
				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 1] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][1]][1];
				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 2] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][1]][2];

				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 3] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][2]][0];
				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 4] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][2]][1];
				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 5] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][2]][2];

				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 6] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][3]][0];
				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 7] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][3]][1];
				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 8] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][3]][2];

				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 9] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][4]][0];
				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 10] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][4]][1];
				send_points[(i * grid->info[smcMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 11] = smc_mesh->points[smc_mesh->cells[(i
						* grid->info[smcMesh][TOTAL_CELLS]) + j][4]][2];
			}
		}
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_POINTS]; i++) {
			free(smc_mesh->points[i]);
		}
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[smcMesh][TOTAL_CELLS]; i++) {
			free(smc_mesh->cells[i]);
		}
		free(smc_mesh->points);
		free(smc_mesh->cells);
		free(smc_mesh);
	}

	recv_points = (double*) malloc(tuple_offset * grid->info[smcMesh][TOTAL_CELLS] * sizeof(double));
	for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS]; i++) {
		send_count[i] = tuple_offset * grid->info[smcMesh][TOTAL_CELLS];
		disp[i] = tuple_offset * grid->info[smcMesh][TOTAL_CELLS] * i;
	}
	recv_count = tuple_offset * grid->info[smcMesh][TOTAL_CELLS];

	check_flag(MPI_Scatterv(send_points, send_count, disp, MPI_DOUBLE, recv_points, recv_count, MPI_DOUBLE, root, grid->cart_comm),
			"error scattering SMC mesh coordinates.ß");

	int count = 0;
	for (int n = 1; n <= grid->num_smc_axially; n++) {
		for (int m = 1; m <= grid->num_smc_circumferentially; m++) {
			for (int l = 0; l < 4; l++) {
				smc[m][n].x_coordinate[l] = recv_points[count * tuple_offset + l * 3 + 0];
				smc[m][n].y_coordinate[l] = recv_points[count * tuple_offset + l * 3 + 1];
				smc[m][n].z_coordinate[l] = recv_points[count * tuple_offset + l * 3 + 2];
			}
			count++;
		}
	}
	free(send_points);
	free(recv_points);

//Reading in EC mesh and communicating to each branch communicator memeber.
	send_points = (double*) malloc(tuple_offset * grid->info[ecMesh][TOTAL_CELLS] * grid->info[ProcessMesh][TOTAL_CELLS] * sizeof(double));

	if (grid->rank == 0) {
		vtk_info* ec_mesh = (vtk_info*) malloc(sizeof(vtk_info));
		ec_mesh->points = (double**) malloc(grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_POINTS] * sizeof(double*));
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_POINTS]; i++) {
			ec_mesh->points[i] = (double*) malloc(3 * sizeof(double));
		}
		ec_mesh->cells = (int**) malloc(grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS] * sizeof(int*));
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS]; i++) {
			ec_mesh->cells[i] = (int*) malloc(5 * sizeof(int));
		}

		indx = read_coordinates(grid,grid->info, ec_mesh, ecMesh, (grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_POINTS]),
				(grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS]));

		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS]; i++) {
			for (int j = 0; j < grid->info[ecMesh][TOTAL_CELLS]; j++) {
				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 0] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][1]][0];
				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 1] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][1]][1];
				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 2] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][1]][2];

				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 3] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][2]][0];
				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 4] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][2]][1];
				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 5] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][2]][2];

				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 6] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][3]][0];
				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 7] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][3]][1];
				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 8] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][3]][2];

				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 9] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][4]][0];
				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 10] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][4]][1];
				send_points[(i * grid->info[ecMesh][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 11] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[ecMesh][TOTAL_CELLS]) + j][4]][2];
			}
		}
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_POINTS]; i++) {
			free(ec_mesh->points[i]);
		}
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecMesh][TOTAL_CELLS]; i++) {
			free(ec_mesh->cells[i]);
		}
		free(ec_mesh->points);
		free(ec_mesh->cells);
		free(ec_mesh);
	}
	recv_points = (double*) malloc(tuple_offset * grid->info[ecMesh][TOTAL_CELLS] * sizeof(double));
	for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS]; i++) {
		send_count[i] = tuple_offset * grid->info[ecMesh][TOTAL_CELLS];
		disp[i] = tuple_offset * grid->info[ecMesh][TOTAL_CELLS] * i;
	}
	recv_count = tuple_offset * grid->info[ecMesh][TOTAL_CELLS];

	check_flag(MPI_Scatterv(send_points, send_count, disp, MPI_DOUBLE, recv_points, recv_count, MPI_DOUBLE, root, grid->cart_comm),
			"error scattering EC mesh coordinates.");
	count = 0;
	for (int n = 1; n <= grid->num_ec_axially; n++) {
		for (int m = 1; m <= grid->num_ec_circumferentially; m++) {
			for (int l = 0; l < 4; l++) {
				ec[m][n].x_coordinate[l] = recv_points[count * tuple_offset + l * 3 + 0];
				ec[m][n].y_coordinate[l] = recv_points[count * tuple_offset + l * 3 + 1];
				ec[m][n].z_coordinate[l] = recv_points[count * tuple_offset + l * 3 + 2];
			}
			count++;
		}
	}
	free(send_points);
	free(recv_points);

//Reading in EC centroids and communicating to each branch communicator memeber.
	num_tuples = 1;
	tuple_offset = num_tuple_components * num_tuples;

	send_points = (double*) malloc(tuple_offset * grid->info[ecCentroids][TOTAL_CELLS] * grid->info[ProcessMesh][TOTAL_CELLS] * sizeof(double));

	if (grid->rank == 0) {
		vtk_info* ec_centroids = (vtk_info*) malloc(sizeof(vtk_info));
		ec_centroids->points = (double**) malloc(grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecCentroids][TOTAL_POINTS] * sizeof(double*));
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecCentroids][TOTAL_POINTS]; i++) {
			ec_centroids->points[i] = (double*) malloc(3 * sizeof(double));
		}
		ec_centroids->cells = (int**) malloc(grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecCentroids][TOTAL_CELLS] * sizeof(int*));
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecCentroids][TOTAL_CELLS]; i++) {
			ec_centroids->cells[i] = (int*) malloc(2 * sizeof(int));
		}

		indx = read_coordinates(grid, grid->info, ec_centroids, ecCentroids,
				(grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecCentroids][TOTAL_POINTS]),
				(grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecCentroids][TOTAL_CELLS]));

		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS]; i++) {
			for (int j = 0; j < grid->info[ecCentroids][TOTAL_CELLS]; j++) {
				send_points[(i * grid->info[ecCentroids][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 0] =
						ec_centroids->points[ec_centroids->cells[(i * grid->info[ecCentroids][TOTAL_CELLS]) + j][1]][0];
				send_points[(i * grid->info[ecCentroids][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 1] =
						ec_centroids->points[ec_centroids->cells[(i * grid->info[ecCentroids][TOTAL_CELLS]) + j][1]][1];
				send_points[(i * grid->info[ecCentroids][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 2] =
						ec_centroids->points[ec_centroids->cells[(i * grid->info[ecCentroids][TOTAL_CELLS]) + j][1]][2];
			}
		}
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecCentroids][TOTAL_POINTS]; i++) {
			free(ec_centroids->points[i]);
		}
		for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS] * grid->info[ecCentroids][TOTAL_CELLS]; i++) {
			free(ec_centroids->cells[i]);
		}
		free(ec_centroids->points);
		free(ec_centroids->cells);
		free(ec_centroids);
	}
	recv_points = (double*) malloc(tuple_offset * grid->info[ecCentroids][TOTAL_CELLS] * sizeof(double));
	for (int i = 0; i < grid->info[ProcessMesh][TOTAL_CELLS]; i++) {
		send_count[i] = tuple_offset * grid->info[ecCentroids][TOTAL_CELLS];
		disp[i] = tuple_offset * grid->info[ecCentroids][TOTAL_CELLS] * i;
	}

	recv_count = tuple_offset * grid->info[ecCentroids][TOTAL_CELLS];
	check_flag(MPI_Scatterv(send_points, send_count, disp, MPI_DOUBLE, recv_points, recv_count, MPI_DOUBLE, root, grid->cart_comm),
			"error scattering EC centroids coordinates.");

	count = 0;
	for (int n = 1; n <= grid->num_ec_axially; n++) {
		for (int m = 1; m <= grid->num_ec_circumferentially; m++) {
			ec[m][n].centeroid_point[0] = recv_points[count * tuple_offset + 0];
			ec[m][n].centeroid_point[1] = recv_points[count * tuple_offset + 1];
			ec[m][n].centeroid_point[2] = recv_points[count * tuple_offset + 2];
			count++;
		}
	}

	free(send_points);
	free(recv_points);

	MPI_Barrier(grid->universe);
	return (0);
}

/*******************************************************************************************************/
int* read_coordinates(grid_parms* grid, int** info, vtk_info* mesh, int mesh_type, int num_points, int num_cells)
/*******************************************************************************************************/
{

	/**
	 * This function reads coordinates from the files in accordance with the information passed by the calling function.
	 * vtk_info = Structure into which the points and cells are to be read in.
	 * Mesh type = ProcessMesh or smcMesh or ecMesh or ecCentroids
	 * num_points = Total number of points in the source file
	 * num_cells  = Total number of vtk cell entries in the source files
	 */
	FILE *fr;
	char filename_points[256], filename_cells[256];
	int* indx = (int*) malloc(2 * sizeof(int));
	indx[0] = 0;
	indx[1] = 0;
	sprintf(grid->input_file_path, "files");
	if (mesh_type == 0) {
		if (grid->my_domain.internal_info.domain_type == STRSEG) {
			sprintf(filename_points, "%s/branch_%d_subdomain_%d_points.txt", grid->input_file_path, grid->my_domain.internal_info.domain_branch,
					grid->my_domain.internal_info.domain_local_subdomain);
		} else if (grid->my_domain.internal_info.domain_type == BIF) {
			if (grid->branch_tag == P) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%d_points.txt", grid->input_file_path, grid->my_domain.internal_info.domain_branch,
						grid->my_domain.internal_info.domain_local_subdomain);
			} else if (grid->branch_tag == L) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%d_points.txt", grid->input_file_path, grid->my_domain.left_child.domain_branch,
						grid->my_domain.left_child.domain_local_subdomain - 1);
			} else if (grid->branch_tag == R) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%d_points.txt", grid->input_file_path, grid->my_domain.right_child.domain_branch,
						grid->my_domain.right_child.domain_local_subdomain - 1);
			}
		}
	} else if (mesh_type == 1) {
		if (grid->my_domain.internal_info.domain_type == STRSEG) {
			sprintf(filename_points, "%s/branch_%d_subdomain_%d_smc_mesh_points.txt", grid->input_file_path,
					grid->my_domain.internal_info.domain_branch, grid->my_domain.internal_info.domain_local_subdomain);
		} else if (grid->my_domain.internal_info.domain_type == BIF) {
			if (grid->branch_tag == P) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%d_smc_mesh_points.txt", grid->input_file_path,
						grid->my_domain.internal_info.domain_branch, grid->my_domain.internal_info.domain_local_subdomain);
			} else if (grid->branch_tag == L) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%d_smc_mesh_points.txt", grid->input_file_path,
						grid->my_domain.left_child.domain_branch, grid->my_domain.left_child.domain_local_subdomain - 1);
			} else if (grid->branch_tag == R) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%d_smc_mesh_points.txt", grid->input_file_path,
						grid->my_domain.right_child.domain_branch, grid->my_domain.right_child.domain_local_subdomain - 1);
			}
		}
	} else if (mesh_type == 2) {
		if (grid->my_domain.internal_info.domain_type == STRSEG) {
			sprintf(filename_points, "%s/branch_%d_subdomain_%d_ec_mesh_points.txt", grid->input_file_path,
					grid->my_domain.internal_info.domain_branch, grid->my_domain.internal_info.domain_local_subdomain);
		} else if (grid->my_domain.internal_info.domain_type == BIF) {
			if (grid->branch_tag == P) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%d_ec_mesh_points.txt", grid->input_file_path,
						grid->my_domain.internal_info.domain_branch, grid->my_domain.internal_info.domain_local_subdomain);
			} else if (grid->branch_tag == L) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%d_ec_mesh_points.txt", grid->input_file_path,
						grid->my_domain.left_child.domain_branch, grid->my_domain.left_child.domain_local_subdomain - 1);
			} else if (grid->branch_tag == R) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%d_ec_mesh_points.txt", grid->input_file_path,
						grid->my_domain.right_child.domain_branch, grid->my_domain.right_child.domain_local_subdomain - 1);
			}
		}
	} else if (mesh_type == 3) {
		if (grid->my_domain.internal_info.domain_type == STRSEG) {
			sprintf(filename_points, "%s/branch_%d_subdomain_%d_ec_centeroid_points.txt", grid->input_file_path,
					grid->my_domain.internal_info.domain_branch, grid->my_domain.internal_info.domain_local_subdomain);
		} else if (grid->my_domain.internal_info.domain_type == BIF) {
			if (grid->branch_tag == P) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%d_ec_centeroid_points.txt", grid->input_file_path,
						grid->my_domain.internal_info.domain_branch, grid->my_domain.internal_info.domain_local_subdomain);
			} else if (grid->branch_tag == L) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%_ec_centeroid_points.txt", grid->input_file_path,
						grid->my_domain.left_child.domain_branch, grid->my_domain.left_child.domain_local_subdomain - 1);
			} else if (grid->branch_tag == R) {
				sprintf(filename_points, "%s/branch_%d_subdomain_%d_ec_centeroid_points.txt", grid->input_file_path,
						grid->my_domain.right_child.domain_branch, grid->my_domain.right_child.domain_local_subdomain - 1);
			}
		}
	}

	fr = fopen(filename_points, "r+");
	for (int i = 0; i < num_points; i++) {
		fscanf(fr, "%lf", &mesh->points[i][0]);
		fscanf(fr, "%lf", &mesh->points[i][1]);
		fscanf(fr, "%lf", &mesh->points[i][2]);
		indx[0]++;
	}
	fclose(fr);

	fr = fopen(filename_cells, "r+");
	if (mesh_type < 3) {
		for (int i = 0; i < num_cells; i++) {
			fscanf(fr, "%d", &mesh->cells[i][0]);
			fscanf(fr, "%d", &mesh->cells[i][1]);
			fscanf(fr, "%d", &mesh->cells[i][2]);
			fscanf(fr, "%d", &mesh->cells[i][3]);
			fscanf(fr, "%d", &mesh->cells[i][4]);
			indx[1]++;
		}
	} else if (mesh_type == 3) {
		for (int i = 0; i < num_cells; i++) {
			fscanf(fr, "%d", &mesh->cells[i][0]);
			fscanf(fr, "%d", &mesh->cells[i][1]);
			indx[1]++;
		}
	}
	fclose(fr);

	return (indx);
}
/************************************************/
void gather_tasks_mesh_point_data_on_writers(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, celltype1** smc,
		celltype2** ec) {

	int branch;
	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		branch = P - 1;
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		branch = grid->branch_tag - 1;
	}

	/************* Point type data *************/
	int num_tuple_components = 3, num_tuples = 4;
	int root = 0;
	int *recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	int *disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	char *send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");
	int length = 0;
	for (int i = 0; i < num_tuples; i++) {
		length += sprintf(send_buffer + length, "%2.9lf\t%2.9lf\t%2.9lf\n", grid->coordinates[i][0], grid->coordinates[i][1],
				grid->coordinates[i][2]);
	}

	/// Gathering and summing the length of all the CHARs contained in every send_buffer containing coordinates from each MPI process.
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[ProcessMesh] = 0;
	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ProcessMesh];
		writer_buffer->buffer_length[ProcessMesh] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->process_mesh_points = (char*) checked_malloc(writer_buffer->buffer_length[ProcessMesh] * sizeof(char),
				"allocation error for writer_buffer member process mesh.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->process_mesh_points, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for task coordinates.");
	free(recv_count);
	free(send_buffer);
	free(disp);
	/************* Cell data *************/
	num_tuple_components = 5;
	num_tuples = 1;
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");

	length = sprintf(send_buffer, "%d %d %d %d %d\n", 4, (branch * grid->tasks * 4) + 4 * grid->rank + 0,
			(branch * grid->tasks * 4) + 4 * grid->rank + 1, (branch * grid->tasks * 4) + 4 * grid->rank + 2,
			(branch * grid->tasks * 4) + 4 * grid->rank + 3);

	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[ProcessCell] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ProcessCell];
		writer_buffer->buffer_length[ProcessCell] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->process_mesh_cells = (char*) checked_malloc(writer_buffer->buffer_length[ProcessCell] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->process_mesh_cells, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for task cells.");
	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* Cell type data *************/
	num_tuple_components = 1;
	num_tuples = 1;
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");

	length = sprintf(send_buffer, "%d\n", 9);

	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[ProcessCellType] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ProcessCellType];
		writer_buffer->buffer_length[ProcessCellType] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->process_mesh_type = (char*) checked_malloc(writer_buffer->buffer_length[ProcessCellType] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->process_mesh_type, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for task cells.");
	free(recv_count);
	free(send_buffer);
	free(disp);
}

/*****************************************************************************************/
void gather_smc_mesh_data_on_writers(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, celltype1** smc) {
	int branch;
	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		branch = P - 1;
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		branch = grid->branch_tag - 1;
	}

	/************* Point data *************/
	int num_tuple_components = 3, num_tuples = 4 * grid->num_smc_axially * grid->num_smc_circumferentially;
	int root = 0;
	int *recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	int *disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	char *send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");
	int length = 0;

	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			for (int k = 0; k < 4; k++) {
				length += sprintf(send_buffer + length, "%2.9lf\t%2.9lf\t%2.9lf\n", smc[i][j].x_coordinate[k], smc[i][j].y_coordinate[k],
						smc[i][j].z_coordinate[k]);
			}
		}
	}

/// Gathering and summing the length of all the CHARs contained in every send_buffer containing coordinates from each MPI process.
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[smcMesh] = 0;
	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[smcMesh];
		writer_buffer->buffer_length[smcMesh] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->smc_mesh_points = (char*) checked_malloc(writer_buffer->buffer_length[smcMesh] * sizeof(char),
				"allocation error for writer_buffer member process mesh.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->smc_mesh_points, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for task coordinates.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* Cell data *************/
	int my_branch_smc_offset = 0;
	num_tuple_components = 5;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");

	int* global_smc_point_offset_info = (int*) checked_malloc(grid->numtasks * sizeof(int),
			"error in global info array allocation in function gather_smc_mesh_data_on_writers");
	int my_smc_point_offset_info = 4 * grid->num_smc_axially * grid->num_smc_circumferentially;
	check_flag(MPI_Allgather(&my_smc_point_offset_info, 1, MPI_INT, global_smc_point_offset_info, 1, MPI_INT, grid->universe),
			"error in mpi_allgatherv in function gather_smc_mesh_data_on_writers");

	length = 0;
	int offset = 0;
	for (int p = grid->universal_rank; p > 0; p--) {
		offset += global_smc_point_offset_info[p];
	}
	for (int i = 0; i < grid->num_smc_axially; i++) {
		for (int j = 0; j < grid->num_smc_circumferentially; j++) {
			length += sprintf(send_buffer + length, "%d %d %d %d %d\n", 4, offset + i * 4 * grid->num_smc_circumferentially + j * 4 + 0,
					offset + i * 4 * grid->num_smc_circumferentially + j * 4 + 1, offset + i * 4 * grid->num_smc_circumferentially + j * 4 + 2,
					offset + i * 4 * grid->num_smc_circumferentially + j * 4 + 3);

		}
	}
	free(global_smc_point_offset_info);
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[smcCell] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[smcCell];
		writer_buffer->buffer_length[smcCell] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->smc_mesh_cells = (char*) checked_malloc(writer_buffer->buffer_length[smcCell] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->smc_mesh_cells, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for task cells.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* Celltype data *************/
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");
	length = 0;
	for (int i = 0; i < grid->num_smc_axially; i++) {
		for (int j = 0; j < grid->num_smc_circumferentially; j++) {
			length += sprintf(send_buffer + length, "%d\n", 9);
		}
	}
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[smcCellType] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[smcCellType];
		writer_buffer->buffer_length[smcCellType] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->smc_mesh_type = (char*) checked_malloc(writer_buffer->buffer_length[smcCellType] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->smc_mesh_type, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for task cells.");
	free(recv_count);
	free(send_buffer);
	free(disp);
}
/*****************************************************************************************/
void gather_ec_mesh_data_on_writers(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, celltype2** ec) {
	int branch;
	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		branch = P - 1;
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		branch = grid->branch_tag - 1;
	}

	/************* Point data *************/
	int num_tuple_components = 3, num_tuples = 4 * grid->num_ec_axially * grid->num_ec_circumferentially;
	int root = 0;
	int *recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	int *disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	char *send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");
	int length = 0;

	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			for (int k = 0; k < 4; k++) {
				length += sprintf(send_buffer + length, "%2.9lf\t%2.9lf\t%2.9lf\n", ec[i][j].x_coordinate[k], ec[i][j].y_coordinate[k],
						ec[i][j].z_coordinate[k]);
			}
		}
	}

/// Gathering and summing the length of all the CHARs contained in every send_buffer containing coordinates from each MPI process.
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[ecMesh] = 0;
	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ecMesh];
		writer_buffer->buffer_length[ecMesh] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->ec_mesh_points = (char*) checked_malloc(writer_buffer->buffer_length[ecMesh] * sizeof(char),
				"allocation error for writer_buffer member process mesh.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->ec_mesh_points, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for task coordinates.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* Cell data *************/
	int my_branch_ec_offset = 0;
	num_tuple_components = 5;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");

	int* global_ec_point_offset_info = (int*) checked_malloc(grid->numtasks * sizeof(int),
			"error in global info array allocation in function gather_ec_mesh_data_on_writers");
	int my_ec_point_offset_info = 4 * grid->num_ec_axially * grid->num_ec_circumferentially;
	check_flag(MPI_Allgather(&my_ec_point_offset_info, 1, MPI_INT, global_ec_point_offset_info, 1, MPI_INT, grid->universe),
			"error in mpi_allgatherv in function gather_ec_mesh_data_on_writers");

	length = 0;
	int offset = 0;
	for (int p = grid->universal_rank; p > 0; p--) {
		offset += global_ec_point_offset_info[p];
	}
	for (int i = 0; i < grid->num_ec_axially; i++) {
		for (int j = 0; j < grid->num_ec_circumferentially; j++) {
			length += sprintf(send_buffer + length, "%d %d %d %d %d\n", 4, offset + i * 4 * grid->num_ec_circumferentially + j * 4 + 0,
					offset + i * 4 * grid->num_ec_circumferentially + j * 4 + 1, offset + i * 4 * grid->num_ec_circumferentially + j * 4 + 2,
					offset + i * 4 * grid->num_ec_circumferentially + j * 4 + 3);

		}
	}
	free(global_ec_point_offset_info);
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[ecCell] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ecCell];
		writer_buffer->buffer_length[ecCell] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->ec_mesh_cells = (char*) checked_malloc(writer_buffer->buffer_length[ecCell] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->ec_mesh_cells, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for task cells.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* Celltype data *************/
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");
	length = 0;
	for (int i = 0; i < grid->num_ec_axially; i++) {
		for (int j = 0; j < grid->num_ec_circumferentially; j++) {
			length += sprintf(send_buffer + length, "%d\n", 9);
		}
	}
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[ecCellType] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ecCellType];
		writer_buffer->buffer_length[ecCellType] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->ec_mesh_type = (char*) checked_malloc(writer_buffer->buffer_length[ecCellType] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->ec_mesh_type, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for task cells.");
	free(recv_count);
	free(send_buffer);
	free(disp);
}

/********************************************************************************************************/
void gather_smcData(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, celltype1** smc, int write_count) {
	int branch;
	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		branch = P - 1;
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		branch = grid->branch_tag - 1;
	}
	/************* smcCa field data *************/
	int num_tuple_components = 1, num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;
	int root = 0;
	int *recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	int *disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	int my_branch_smc_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	char* send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_smcData.");

	int length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].p[smc_Ca]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->smc_stat_var_buffer_length[smc_Ca] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_stat_var_buffer_length[smc_Ca];
		writer_buffer->smc_stat_var_buffer_length[smc_Ca] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->ci = (char*) checked_malloc(writer_buffer->smc_stat_var_buffer_length[smc_Ca] * sizeof(char),
				"allocation error for writer_buffer member ci.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->ci, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for smc Ca.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* smcSR field data *************/
	num_tuple_components = 1, num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	my_branch_smc_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_smcData.");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].p[smc_SR]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->smc_stat_var_buffer_length[smc_SR] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_stat_var_buffer_length[smc_SR];
		writer_buffer->smc_stat_var_buffer_length[smc_SR] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->si = (char*) checked_malloc(writer_buffer->smc_stat_var_buffer_length[smc_SR] * sizeof(char),
				"allocation error for writer_buffer member si.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->si, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for smc SR.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* smcVm field data *************/
	num_tuple_components = 1, num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	my_branch_smc_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_smcData.");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].p[smc_Vm]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->smc_stat_var_buffer_length[smc_Vm] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_stat_var_buffer_length[smc_Vm];
		writer_buffer->smc_stat_var_buffer_length[smc_Vm] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->vi = (char*) checked_malloc(writer_buffer->smc_stat_var_buffer_length[smc_Vm] * sizeof(char),
				"allocation error for writer_buffer member vi.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->vi, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for smc Vm.");

	free(recv_count);
	free(send_buffer);
	free(disp);
	/************* smc_w field data *************/
	num_tuple_components = 1, num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	my_branch_smc_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_smcData.");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].p[smc_w]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->smc_stat_var_buffer_length[smc_w] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_stat_var_buffer_length[smc_w];
		writer_buffer->smc_stat_var_buffer_length[smc_w] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->wi = (char*) checked_malloc(writer_buffer->smc_stat_var_buffer_length[smc_w] * sizeof(char),
				"allocation error for writer_buffer member wi.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->wi, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for smc w.");

	free(recv_count);
	free(send_buffer);
	free(disp);
	/************* smc_IP3 field data *************/
	num_tuple_components = 1, num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	my_branch_smc_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_smcData.");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].p[smc_IP3]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->smc_stat_var_buffer_length[smc_IP3] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_stat_var_buffer_length[smc_IP3];
		writer_buffer->smc_stat_var_buffer_length[smc_IP3] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->Ii = (char*) checked_malloc(writer_buffer->smc_stat_var_buffer_length[smc_IP3] * sizeof(char),
				"allocation error for writer_buffer member Ii.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->Ii, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for smc IP3.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* smc_cpl_Ca field data *************/
	my_branch_smc_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_smcData.");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].B[cpl_Ca]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->smc_cpl[cpl_Ca] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_cpl[cpl_Ca];
		writer_buffer->smc_cpl[cpl_Ca] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpCi = (char*) checked_malloc(writer_buffer->smc_cpl[cpl_Ca] * sizeof(char),
				"allocation error for writer_buffer member smc_cplCa.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpCi, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for smc_cplCa.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* smc_cpl_Vm field data *************/
	my_branch_smc_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_smcData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].B[cpl_Vm]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->smc_cpl[cpl_Vm] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_cpl[cpl_Vm];
		writer_buffer->smc_cpl[cpl_Vm] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpVi = (char*) checked_malloc(writer_buffer->smc_cpl[cpl_Vm] * sizeof(char),
				"allocation error for writer_buffer member smc_cplVm.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpVi, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for smc_cplVm.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* smc_cpl_IP3 field data *************/
	my_branch_smc_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_smcData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].B[cpl_IP3]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->smc_cpl[cpl_IP3] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_cpl[cpl_IP3];
		writer_buffer->smc_cpl[cpl_IP3] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpIi = (char*) checked_malloc(writer_buffer->smc_cpl[cpl_IP3] * sizeof(char),
				"allocation error for writer_buffer member smc_cplIi.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpIi, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for smc_cplIi.");

	free(recv_count);
	free(send_buffer);
	free(disp);
}
/********************************************************************************************************/
void gather_ecData(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, celltype2** ec, int write_count) {
	int branch;
	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		branch = P - 1;
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		branch = grid->branch_tag - 1;
	}

	int num_tuple_components = 1, num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;
	int root = 0;
	int *recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	int *disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	/************* ecCa field data *************/
	int my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	char* send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_ecData.");

	int length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].q[ec_Ca]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->ec_stat_var_buffer_length[ec_Ca] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_stat_var_buffer_length[ec_Ca];
		writer_buffer->ec_stat_var_buffer_length[ec_Ca] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cj = (char*) checked_malloc(writer_buffer->ec_stat_var_buffer_length[ec_Ca] * sizeof(char),
				"allocation error for writer_buffer member cj.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cj, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for ec Ca.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* ecSR field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].q[ec_SR]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->ec_stat_var_buffer_length[ec_SR] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_stat_var_buffer_length[ec_SR];
		writer_buffer->ec_stat_var_buffer_length[ec_SR] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->sj = (char*) checked_malloc(writer_buffer->ec_stat_var_buffer_length[ec_SR] * sizeof(char),
				"allocation error for writer_buffer member sj.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->sj, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for ec SR.");

	free(recv_count);
	free(send_buffer);
	free(disp);
	/************* ecVm field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].q[ec_Vm]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->ec_stat_var_buffer_length[ec_Vm] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_stat_var_buffer_length[ec_Vm];
		writer_buffer->ec_stat_var_buffer_length[ec_Vm] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->vj = (char*) checked_malloc(writer_buffer->ec_stat_var_buffer_length[ec_Vm] * sizeof(char),
				"allocation error for writer_buffer member cj.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->vj, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for ec Vm.");

	free(recv_count);
	free(send_buffer);
	free(disp);
	/************* ec_IP3 field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].q[ec_IP3]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->ec_stat_var_buffer_length[ec_IP3] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_stat_var_buffer_length[ec_IP3];
		writer_buffer->ec_stat_var_buffer_length[ec_IP3] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->Ij = (char*) checked_malloc(writer_buffer->ec_stat_var_buffer_length[ec_Vm] * sizeof(char),
				"allocation error for writer_buffer member Ij.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->Ij, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for ec Ii.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* ec_cpl_Ca field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].B[cpl_Ca]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->ec_cpl[cpl_Ca] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_cpl[cpl_Ca];
		writer_buffer->ec_cpl[cpl_Ca] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpCj = (char*) checked_malloc(writer_buffer->ec_cpl[cpl_Ca] * sizeof(char),
				"allocation error for writer_buffer member ec_cplCa.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpCj, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for ec_cplCa.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* ec_cpl_Vm field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].B[cpl_Vm]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->ec_cpl[cpl_Vm] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_cpl[cpl_Vm];
		writer_buffer->ec_cpl[cpl_Vm] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpVj = (char*) checked_malloc(writer_buffer->ec_cpl[cpl_Vm] * sizeof(char),
				"allocation error for writer_buffer member ec_cplVm.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpVj, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for ec_cplVm.");

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* ec_cpl_IP3 field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].B[cpl_IP3]);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->ec_cpl[cpl_IP3] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_cpl[cpl_IP3];
		writer_buffer->ec_cpl[cpl_IP3] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpIj = (char*) checked_malloc(writer_buffer->ec_cpl[cpl_IP3] * sizeof(char),
				"allocation error for writer_buffer member ec_cplIi.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpIj, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for ec_cplIi.");

	free(recv_count);
	free(send_buffer);
	free(disp);
}
void gather_JPLC_map(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, celltype2** ec) {
	int branch;
	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		branch = P - 1;
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		branch = grid->branch_tag - 1;
	}

	int num_tuple_components = 1, num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;
	int root = 0;
	int *recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	int *disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	/************* JPLC field data *************/
	int my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	char* send_buffer = (char*) checked_malloc(
	NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), "error allocating send_buffer in gather_ecData.");

	int length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].JPLC);
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");

	writer_buffer->jplc_buffer_length = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->jplc_buffer_length;
		writer_buffer->jplc_buffer_length += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->jplc = (char*) checked_malloc(writer_buffer->jplc_buffer_length * sizeof(char),
				"allocation error for writer_buffer member jplc.");
	}
	check_flag(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->jplc, recv_count, disp, MPI_CHAR, root, grid->cart_comm),
			"Error gathering data for jplc.");

	free(recv_count);
	free(send_buffer);
	free(disp);
}

void checkpoint_coarse_time_profiling_data(grid_parms grid, time_stamps* t_stamp, IO_domain_info* my_IO_domain_info) {

	push_coarse_timing_data_to_file("aggregated_compute_time", grid, t_stamp->aggregate_compute, my_IO_domain_info);
	push_coarse_timing_data_to_file("aggregated_comm_time", grid, t_stamp->aggregate_comm, my_IO_domain_info);
	push_coarse_timing_data_to_file("aggregated_write_time", grid, t_stamp->aggregate_write, my_IO_domain_info);

	push_task_wise_min_max_of_time_profile("min_max_of_aggregate_compute", grid, t_stamp->aggregate_compute, my_IO_domain_info);
	push_task_wise_min_max_of_time_profile("min_max_of_aggregate_comm", grid, t_stamp->aggregate_comm, my_IO_domain_info);
	push_task_wise_min_max_of_time_profile("min_max_of_aggregate_write", grid, t_stamp->aggregate_write, my_IO_domain_info);
}
void push_coarse_timing_data_to_file(char* file_prefix, grid_parms grid, double field, IO_domain_info* my_IO_domain_info) {
	MPI_Status status;
	MPI_Offset displacement = 0;
	MPI_File fw;
	char* buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * sizeof(char), "allocation for logfile segment space\n");
	char* write_buffer;
	int root = 0;
	char filename[50];

	int length = sprintf(buffer, "%2.12lf\n", field);

	int *recv_count = (int*) checked_malloc(grid.tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	int *disp = (int*) checked_malloc(grid.tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	/// Gathering and summing the length of all the CHARs contained in every send_buffer containing coordinates from each MPI process.
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid.cart_comm),
			"error in MPI_Gather gathering Logfile buffer length by each process.");
	int total_buffer_length = 0;
	for (int i = 0; i < grid.tasks; i++) {
		disp[i] = total_buffer_length;
		total_buffer_length += recv_count[i];
	}

	if (grid.rank == 0) {
		write_buffer = (char*) checked_malloc(total_buffer_length * sizeof(char),
				"allocation error for writer_buffer for time_profiling data in function push_coarse_timing_data_to_file().\n");
	}
	check_flag(MPI_Gatherv(buffer, length, MPI_CHAR, write_buffer, recv_count, disp, MPI_CHAR, root, grid.cart_comm),
			"Error gathering time_profiling data.");

	if (grid.rank == 0) {
		sprintf(filename, "%s/%s_%s.txt", grid.time_profiling_dir, file_prefix, grid.suffix);
		check_flag(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fw), "error opening Time_profiling data");
		check_flag(MPI_File_write_at(fw, 0, write_buffer, total_buffer_length, MPI_CHAR, &status), "Error writing data into Time profiling.");
		MPI_File_close(&fw);
		free(write_buffer);
	}
	free(recv_count);
	free(buffer);
	free(disp);
}

void push_task_wise_min_max_of_time_profile(char* file_prefix, grid_parms grid, double field, IO_domain_info* my_IO_domain_info) {

	MPI_Status status;
	MPI_Offset displacement = 0;
	MPI_File fw;
	char* buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * sizeof(char), "allocation for logfile segment space\n");
	char* write_buffer;
	int root = 0;
	char filename[50];
	int length = 0;
	double max, min;
	int max_ind, min_ind;
	double array[grid.tasks];
	int *recv_count = (int*) checked_malloc(grid.tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	int *disp = (int*) checked_malloc(grid.tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	check_flag(MPI_Gather(&field, 1, MPI_DOUBLE, array, 1, MPI_DOUBLE, root, grid.cart_comm),
			"error in MPI_Gather gathering aggregate_time_profiling data in function push_task_wise_min_max_of_time_profile().\n");
	if (grid.rank == 0) {
		maximum(array, grid.tasks, &max, &max_ind);
		minimum(array, grid.tasks, &min, &min_ind);
		write_buffer = (char*) checked_malloc(1024 * sizeof(char), "error allocating write_buffer for pushing max min values for time profiling.\n");
		sprintf(filename, "%s/%s_%s.txt", grid.time_profiling_dir, file_prefix, grid.suffix);
		check_flag(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fw), "error opening Time_profiling data");

		length = sprintf(write_buffer, "Maximum = %lf\t by Rank = %d\nMinimum = %lf\t by Rank = %d\n", max, max_ind, min, min_ind);
		check_flag(MPI_File_write_at(fw, 0, write_buffer, length, MPI_CHAR, &status), "Error writing data into Time profiling.");
		MPI_File_close(&fw);
		free(write_buffer);
	}
	free(recv_count);
	free(buffer);
	free(disp);
}

