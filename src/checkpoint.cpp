#include <assert.h>
#include "computelib.h"

#define NUM_CONFIG_ELEMENTS 9

checkpoint_handle* initialise_checkpoint(grid_parms grid)
{
	checkpoint_handle *check = (checkpoint_handle*) malloc(sizeof(checkpoint_handle));
	return (check);
}

#if 0
void open_common_checkpoint(checkpoint_handle* check, grid_parms grid) {

	int err;
	char filename[50];
	err = sprintf(filename, "Log_file%s.txt", grid.suffix);
	CHECK_MPI_ERROR(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->logptr));

	err = sprintf(filename, "Elasped_time%s.txt", grid.suffix);
	CHECK_MPI_ERROR(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->elapsed_time));

	err = sprintf(filename, "JPLC%s.txt", grid.suffix);
	CHECK_MPI_ERROR(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->jplc));

	err = sprintf(filename, "coords%s.txt", grid.suffix);
	CHECK_MPI_ERROR(MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR,MPI_INFO_NULL, &check->coords));
}
#endif

void open_koenigsberger_smc_checkpoint(checkpoint_handle* check, grid_parms grid, int write_count, char* path, IO_domain_info* my_IO_domain_info)
{
	int err;
	char filename[50];
	err = sprintf(filename, "%s/smc_Data_t_%d.vtk", path, write_count);
	CHECK_MPI_ERROR(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->smc_data_file));
	/*err = sprintf(filename, "%s/smc_Ca_t_%d.vtk", path, write_count);
	 CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->ci));
	 err = sprintf(filename, "%s/smc_SERCA_t_%d.vtk", path, write_count);
	 CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->si));
	 err = sprintf(filename, "%s/smc_V_t_%d.vtk", path, write_count);
	 CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->vi));
	 err = sprintf(filename, "%s/smc_KCa_t_%d.vtk", path, write_count);
	 CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->wi));
	 err = sprintf(filename, "%s/smc_IP3_t_%d.vtk", path, write_count);
	 CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Ii));*/
}

void open_koenigsberger_ec_checkpoint(checkpoint_handle* check, grid_parms grid, int write_count, char* path, IO_domain_info* my_IO_domain_info) {
	int err;
	char filename[50];
	err = sprintf(filename, "%s/ec_Data_t_%d.vtk", path, write_count);
	CHECK_MPI_ERROR(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->ec_data_file));
	/*err = sprintf(filename, "%s/ec_Ca_t_%d.vtk", path, write_count);
	 CHECK_MPI_ERROR(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cj));
	 err = sprintf(filename, "%s/ec_SERCA_t_%d.vtk", path, write_count);
	 CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->sj));
	 err = sprintf(filename, "%s/ec_V_t_%d.vtk", path, write_count);
	 CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->vj));
	 err = sprintf(filename, "%s/ec_IP3_t_%d.vtk", path, write_count);
	 CHECK(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->Ij));*/
}

#if 0
void open_coupling_data_checkpoint(checkpoint_handle* check, grid_parms grid, int write_count, char* path, IO_domain_info* my_IO_domain_info)
{
	int err;
	char filename[50];
	err = sprintf(filename, "%s/smc_cpc_t_%d.vtk", path, write_count);
	CHECK_MPI_ERROR(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpCi));

	err = sprintf(filename, "%s/ec_cpc_t_%d.vtk", path, write_count);
	CHECK_MPI_ERROR(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpCj));

	err = sprintf(filename, "%s/smc_cpV_t_%d.vtk", path, write_count);
	CHECK_MPI_ERROR(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpVi));

	err = sprintf(filename, "%s/ec_cpV_t_%d.vtk", path, write_count);
	CHECK_MPI_ERROR(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpVj));

	err = sprintf(filename, "%s/smc_cpIP3_t_%d.vtk", path, write_count);
	CHECK_MPI_ERROR(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpIi));

	err = sprintf(filename, "%s/ec_cpIP3_t_%d.vtk", path, write_count);
	CHECK_MPI_ERROR(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->cpIj));
}
#endif

#if 0
void dump_smc(grid_parms grid, SMC_cell **smc, checkpoint_handle *check, int line_number, int write_count)
{
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
	CHECK_MPI_ERROR(MPI_File_write_at(check->ci, disp, &b1, write_element_count, MPI_DOUBLE, &status[0]));
	CHECK_MPI_ERROR(MPI_File_write_at(check->vi, disp, &b3, write_element_count, MPI_DOUBLE, &status[2]));
	CHECK_MPI_ERROR(MPI_File_write_at(check->Ii, disp, &b5, write_element_count, MPI_DOUBLE, &status[4]));
	CHECK_MPI_ERROR(MPI_File_write_at(check->si, disp, &b2, write_element_count, MPI_DOUBLE, &status[1]));

	CHECK_MPI_ERROR(MPI_File_write_at(check->wi, disp, &b4, write_element_count, MPI_DOUBLE, &status[3]));

	CHECK_MPI_ERROR(MPI_File_write_at(check->cpCi, disp, &b6, write_element_count, MPI_DOUBLE, &status[5]));
	CHECK_MPI_ERROR(MPI_File_write_at(check->cpVi, disp, &b7, write_element_count, MPI_DOUBLE, &status[6]));
	CHECK_MPI_ERROR(MPI_File_write_at(check->cpIi, disp, &b8, write_element_count, MPI_DOUBLE, &status[7]));
}

void dump_ec(grid_parms grid, EC_cell **ec, checkpoint_handle *check, int line_number, int write_count)
{
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
	CHECK_MPI_ERROR(MPI_File_write_at(check->cj, disp, &b1, write_element_count, MPI_DOUBLE, &status[0]));

	CHECK_MPI_ERROR(MPI_File_write_at(check->vj, disp, &b3, write_element_count, MPI_DOUBLE, &status[2]));
	CHECK_MPI_ERROR(MPI_File_write_at(check->Ij, disp, &b4, write_element_count, MPI_DOUBLE, &status[3]));
	CHECK_MPI_ERROR(MPI_File_write_at(check->sj, disp, &b2, write_element_count, MPI_DOUBLE, &status[1]));
	CHECK_MPI_ERROR(MPI_File_write_at(check->cpCj, disp, &b5, write_element_count, MPI_DOUBLE, &status[4]));
	CHECK_MPI_ERROR(MPI_File_write_at(check->cpVj, disp, &b6, write_element_count, MPI_DOUBLE, &status[5]));
	CHECK_MPI_ERROR(MPI_File_write_at(check->cpIj, disp, &b7, write_element_count, MPI_DOUBLE, &status[6]));
}
#endif

void write_smc_and_ec_data(checkpoint_handle* check, grid_parms* grid, double tnow, SMC_cell** smc, EC_cell** ec, int write_count,
		IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer) {

	dump_smc_data(check, grid, my_IO_domain_info, writer_buffer, smc, write_count);
	dump_ec_data(check, grid, my_IO_domain_info, writer_buffer, ec, write_count);
}
/*********************************************************************************************************************************/
void write_process_mesh(checkpoint_handle* check, grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, char* path)
/*********************************************************************************************************************************/
{
	MPI_Status status;
	MPI_Offset disp;
	char filename[50];

	int err = sprintf(filename, "%s/task_mesh.vtk", path);
	printf("[%d] ======>>>>>> Entering %s:%s to write %s \n", grid->universal_rank, __FILE__, __FUNCTION__, filename);

	CHECK_MPI_ERROR(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->task_mesh));

	int write_element_count = 1, header_offset[3] = { 0, 0, 0 }, point_offset = 0, cell_offset = 0, celltype_offset = 0, count = 0;
	char* header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
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
			"POINTS %d double\n", grid->info[PROCESS_MESH][TOTAL_CELLS] * 4 * branches);

	header_offset[0] = strlen(header);
	count = header_offset[0];
	disp = 0;
	if (my_IO_domain_info->writer_rank == 0)
	{
		CHECK_MPI_ERROR(MPI_File_write_at(check->task_mesh, disp, header, count, MPI_CHAR, &status));
	}

	/*************** Writing Point data **************/
	int buffer_lengths[4] = { 0, 0, 0, 0 };
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[PROCESS_MESH], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}

	disp = (header_offset[0] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[PROCESS_MESH];
	CHECK_MPI_ERROR(MPI_File_write_at(check->task_mesh, disp, writer_buffer->process_mesh_points, count, MPI_CHAR, &status));

	for (int i = 0; i < 4; i++) {
		point_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->process_mesh_points);
	/*************** Writing cell data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[1] = sprintf(header, "CELLS %d %d\n",
			branches * grid->info[PROCESS_MESH][TOTAL_CELLS],
			5 * branches * grid->info[PROCESS_MESH][TOTAL_CELLS]);

	count = header_offset[1];
	disp = (header_offset[0] + point_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->task_mesh, disp, header, count, MPI_CHAR, &status));
	}

	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[ProcessCell], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + disp + header_offset[1]) * sizeof(char);
	count = writer_buffer->buffer_length[ProcessCell];
	CHECK_MPI_ERROR(MPI_File_write_at(check->task_mesh, disp, writer_buffer->process_mesh_cells, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		cell_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->process_mesh_cells);
	/*************** Writing cell type data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[2] = sprintf(header, "CELL_TYPES %d\n", branches * grid->info[PROCESS_MESH][TOTAL_CELLS]);

	count = header_offset[2];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0)
	{
		CHECK_MPI_ERROR(MPI_File_write_at(check->task_mesh, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[ProcessCellType], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[ProcessCellType];
	CHECK_MPI_ERROR(MPI_File_write_at(check->task_mesh, disp, writer_buffer->process_mesh_type, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++)
	{
		celltype_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->process_mesh_type);

	MPI_File_close(&check->task_mesh);
}

/**************************************************************************************/
void dump_smc_data(checkpoint_handle* check, grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, SMC_cell** smc,
		int write_count) {
	MPI_Status status;
	MPI_Offset disp;

	printf("[%d] ------>>>>>> Entering %s:%s to write %d\n", grid->universal_rank, __FILE__, __FUNCTION__, check->smc_data_file);

	int write_element_count = 1, *header_offset, point_offset = 0, cell_offset = 0, celltype_offset = 0, *smcDataOffset;
	header_offset = (int*) checked_malloc((3 + grid->neq_smc + grid->num_coupling_species_smc) * sizeof(int), SRC_LOC);
	smcDataOffset = (int*) checked_malloc((grid->neq_smc + grid->num_coupling_species_smc) * sizeof(int), SRC_LOC);

	for (int i = 0; i < (3 + grid->neq_smc + grid->num_coupling_species_smc); i++)
	{
		header_offset[i] = 0;
	}
	for (int i = 0; i < (grid->neq_smc + grid->num_coupling_species_smc); i++)
	{
		smcDataOffset[i] = 0;
	}

	int count = 0;
	char* header = (char*) checked_malloc(1024 * sizeof(char), SRC_LOC);
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
			"POINTS %d double\n", write_count, grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * 4 * branches);
	count = header_offset[0];
	disp = 0;
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status));
	}
	/*************** Writing Point data **************/
	int buffer_lengths[4] = { 0, 0, 0, 0 };
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[SMC_MESH], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}

	disp = (header_offset[0] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[SMC_MESH];
	CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->smc_mesh_points, count, MPI_CHAR, &status));

	for (int i = 0; i < 4; i++) {
		point_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->smc_mesh_points);

	/*************** Writing cell data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[1] = sprintf(header, "CELLS %d %d\n", (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches),
			5 * (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches));

	count = header_offset[1];
	disp = (header_offset[0] + point_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status));
	}

	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[smcCell], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + disp + header_offset[1]) * sizeof(char);
	count = writer_buffer->buffer_length[smcCell];
	CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->smc_mesh_cells, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		cell_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->smc_mesh_cells);

	/*************** Writing cell type data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[2] = sprintf(header, "CELL_TYPES %d\n", (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches));

	count = header_offset[2];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[smcCellType], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));
	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[smcCellType];
	CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->smc_mesh_type, count, MPI_CHAR, &status));
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
			"SMC_Ca %d %d float\n", grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches,
			grid->neq_smc + grid->num_coupling_species_smc, 1, grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[3];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->smc_stat_var_buffer_length[smc_Ca], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + disp)
			* sizeof(char);
	count = writer_buffer->smc_stat_var_buffer_length[smc_Ca];
	CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->ci, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		smcDataOffset[0] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ci);

	/*************** Writing Field 2 : SMC SR data ***************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[4] = sprintf(header, "SMC_SR %d %d float\n", 1, grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[4];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status));
	}

	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->smc_stat_var_buffer_length[smc_SR], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + disp) * sizeof(char);
	count = writer_buffer->smc_stat_var_buffer_length[smc_SR];
	CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->si, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		smcDataOffset[1] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->si);

	/*************** Writing Field 3 : SMC Vm data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[5] = sprintf(header, "SMC_Vm %d %d float\n", 1, grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[5];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->smc_stat_var_buffer_length[smc_Vm], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + disp) * sizeof(char);
	count = writer_buffer->smc_stat_var_buffer_length[smc_Vm];
	CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->vi, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		smcDataOffset[2] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->vi);
	/*************** Writing Field 4 : SMC w data ***************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[6] = sprintf(header, "SMC_w %d %d float\n", 1, grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[6];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->smc_stat_var_buffer_length[smc_w], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + disp) * sizeof(char);
	count = writer_buffer->smc_stat_var_buffer_length[smc_w];
	CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->wi, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		smcDataOffset[3] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->wi);
	/*************** Writing Field 5 : SMC I data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[7] = sprintf(header, "SMC_IP3 %d %d float\n", 1,
			grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[7];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3])
			* sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->smc_stat_var_buffer_length[smc_IP3], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + disp) * sizeof(char);
	count = writer_buffer->smc_stat_var_buffer_length[smc_IP3];
	CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->Ii, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		smcDataOffset[4] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->Ii);
	/*************** Writing Field 6 : SMC Ca coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[8] = sprintf(header, "SMC_Ca_coupling %d %d float\n", 1,
			grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[8];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->smc_cpl[cpl_Ca], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4] + header_offset[8] + disp) * sizeof(char);
	count = writer_buffer->smc_cpl[cpl_Ca];
	CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->cpCi, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		smcDataOffset[5] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpCi);
	/*************** Writing Field 7 : SMC Vm coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[9] = sprintf(header, "SMC_Vm_coupling %d %d float\n", 1,
			grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[9];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4] + header_offset[8] + smcDataOffset[5]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->smc_cpl[cpl_Vm], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4] + header_offset[8] + smcDataOffset[5] + header_offset[9] + disp) * sizeof(char);
	count = writer_buffer->smc_cpl[cpl_Vm];
	CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->cpVi, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		smcDataOffset[6] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpVi);

	/*************** Writing Field 8 : SMC IP3 coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[10] = sprintf(header, "SMC_IP3_coupling %d %d float\n", 1,
			grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[10];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4] + header_offset[8] + smcDataOffset[5] + header_offset[9] + smcDataOffset[6]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->smc_cpl[cpl_IP3], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
			+ smcDataOffset[0] + header_offset[4] + smcDataOffset[1] + header_offset[5] + smcDataOffset[2] + header_offset[6] + smcDataOffset[3]
			+ header_offset[7] + smcDataOffset[4] + header_offset[8] + smcDataOffset[5] + header_offset[9] + smcDataOffset[6] + header_offset[10]
			+ disp) * sizeof(char);
	count = writer_buffer->smc_cpl[cpl_IP3];
	CHECK_MPI_ERROR(MPI_File_write_at(check->smc_data_file, disp, writer_buffer->cpIi, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		smcDataOffset[5] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpIi);

	printf("[%d] ------>>>>>> Leaving %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);
}

void write_JPLC_map(checkpoint_handle* check, grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, EC_cell** ec, char* path)
{
	MPI_Status status;
	MPI_Offset disp;
	char filename[50];
	int err = sprintf(filename, "%s/Agonist_map.vtk", path);
	printf("[%d] ======>>>>>> Entering %s:%s to write %s\n", grid->universal_rank, __FILE__, __FUNCTION__, filename);

	CHECK_MPI_ERROR(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->ec_agonist_file));

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
			"POINTS %d double\n", grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * 4 * branches);
	count = header_offset[0];
	disp = 0;
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_agonist_file, disp, header, count, MPI_CHAR, &status));
	}

	/*************** Writing Point data **************/
	int buffer_lengths[4] = { 0, 0, 0, 0 };
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[EC_MESH], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}

	disp = (header_offset[0] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[EC_MESH];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_agonist_file, disp, writer_buffer->ec_mesh_points, count, MPI_CHAR, &status));

	for (int i = 0; i < 4; i++) {
		point_offset += buffer_lengths[i];
	}

	free(header);
	free(writer_buffer->ec_mesh_points);

	/*************** Writing cell data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[1] = sprintf(header, "CELLS %d %d\n", (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches),
			5 * (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches));

	count = header_offset[1];
	disp = (header_offset[0] + point_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_agonist_file, disp, header, count, MPI_CHAR, &status));
	}

	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[ecCell], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + disp + header_offset[1]) * sizeof(char);
	count = writer_buffer->buffer_length[ecCell];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_agonist_file, disp, writer_buffer->ec_mesh_cells, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		cell_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ec_mesh_cells);

	/*************** Writing cell type data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "Allocation memory for writing header failed at MPI_COMM_WORLD Rank 0.");
	header_offset[2] = sprintf(header, "CELL_TYPES %d\n", (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches));

	count = header_offset[2];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_agonist_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[ecCellType], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[ecCellType];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_agonist_file, disp, writer_buffer->ec_mesh_type, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		celltype_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ec_mesh_type);


	/******** Writing Field 1 : JPLC Data ********/
	header = (char*) checked_malloc(1024 * sizeof(char), "Allocation memory for writing header failed at MPI_COMM_WORLD Rank 0.");
	header_offset[3] = sprintf(header,
			"CELL_DATA %d\nFIELD ec_Data %d\n"
			"JPLC %d %d float\n", grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches,
	/*grid->num_parameters*/1, 1, grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[3];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0)
	{
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_agonist_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->jplc_buffer_length, 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + disp)
			* sizeof(char);
	count = writer_buffer->jplc_buffer_length;
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_agonist_file, disp, writer_buffer->jplc, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		ecDataOffset[0] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->jplc);

	MPI_File_close(&check->ec_agonist_file);

	// printf("[%d] Leaving %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);
}

void dump_ec_data(checkpoint_handle* check, grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, EC_cell** ec,
		int write_count) {
	MPI_Status status;
	MPI_Offset disp;

	printf("[%d] ++++++>>>>>> Entering %s:%s to write %d\n", grid->universal_rank, __FILE__, __FUNCTION__, check->ec_data_file);

	int write_element_count = 1, *header_offset, point_offset = 0, cell_offset = 0, celltype_offset = 0, *ecDataOffset;
	header_offset = (int*) checked_malloc((3 + grid->neq_ec + grid->num_coupling_species_ec) * sizeof(int),
			"allocation failed for header_offset array in dump_ec_data.");
	ecDataOffset = (int*) checked_malloc((grid->neq_ec + grid->num_coupling_species_ec) * sizeof(int),
			"allocation failed for ecDataOffset array in dump_ec_data.");
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
			"POINTS %d double\n", write_count, grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * 4 * branches);
	count = header_offset[0];
	disp = 0;
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status));
	}
	/*************** Writing Point data **************/
	int buffer_lengths[4] = { 0, 0, 0, 0 };
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[EC_MESH], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}

	disp = (header_offset[0] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[EC_MESH];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->ec_mesh_points, count, MPI_CHAR, &status));

	for (int i = 0; i < 4; i++) {
		point_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ec_mesh_points);

	/*************** Writing cell data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[1] = sprintf(header, "CELLS %d %d\n", (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches),
			5 * (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches));

	count = header_offset[1];
	disp = (header_offset[0] + point_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status));
	}

	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[ecCell], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + disp + header_offset[1]) * sizeof(char);
	count = writer_buffer->buffer_length[ecCell];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->ec_mesh_cells, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		cell_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ec_mesh_cells);

	/*************** Writing cell type data **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[2] = sprintf(header, "CELL_TYPES %d\n", (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches));

	count = header_offset[2];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->buffer_length[ecCellType], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + disp) * sizeof(char);
	count = writer_buffer->buffer_length[ecCellType];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->ec_mesh_type, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		celltype_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->ec_mesh_type);

	/***************************************************************/
	/********		Writing Field 1 : ec Ca data 			********/
	/***************************************************************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[3] = sprintf(header, "CELL_DATA %d\nFIELD EC_Data %d\n"
			"EC_Ca %d %d float\n", grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches,
			grid->neq_ec + grid->num_coupling_species_ec, 1, grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[3];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->ec_stat_var_buffer_length[ec_Ca], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + disp)
			* sizeof(char);
	count = writer_buffer->ec_stat_var_buffer_length[ec_Ca];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->cj, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		ecDataOffset[0] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cj);

	/*************** Writing Field 2 : ec SR data ***************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[4] = sprintf(header, "EC_SR %d %d float\n", 1, grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[4];
	disp =
			(header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3]
					+ ecDataOffset[0]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status));
	}

	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->ec_stat_var_buffer_length[ec_SR], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + disp) * sizeof(char);
	count = writer_buffer->ec_stat_var_buffer_length[ec_SR];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->sj, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		ecDataOffset[1] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->sj);

	/*************** Writing Field 3 : ec Vm data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[5] = sprintf(header, "EC_Vm %d %d float\n", 1, grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[5];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->ec_stat_var_buffer_length[ec_Vm], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + disp) * sizeof(char);
	count = writer_buffer->ec_stat_var_buffer_length[ec_Vm];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->vj, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		ecDataOffset[2] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->vj);
	/*************** Writing Field 4 : ec I data ***************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[6] = sprintf(header, "EC_IP3 %d %d float\n", 1, grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[6];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->ec_stat_var_buffer_length[ec_IP3], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + disp) * sizeof(char);
	count = writer_buffer->ec_stat_var_buffer_length[ec_IP3];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->Ij, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		ecDataOffset[3] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->Ij);
	/*************** Writing Field 5 : ec Ca coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[7] = sprintf(header, "EC_Ca_coupling %d %d float\n", 1,
			grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[7];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->ec_cpl[cpl_Ca], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3] + header_offset[7] + disp)
			* sizeof(char);
	count = writer_buffer->ec_cpl[cpl_Ca];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->cpCj, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		ecDataOffset[4] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpCj);
	/*************** Writing Field 6 : ec Vm coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[8] = sprintf(header, "EC_Vm_coupling %d %d float\n", 1,
			grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[8];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3] + header_offset[7]
			+ ecDataOffset[4]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->ec_cpl[cpl_Vm], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3] + header_offset[7]
			+ ecDataOffset[4] + header_offset[8] + disp) * sizeof(char);
	count = writer_buffer->ec_cpl[cpl_Vm];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->cpVj, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		ecDataOffset[5] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpVj);
	/*************** Writing Field 7 : ec IP3 coupling data ***************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0");
	header_offset[9] = sprintf(header, "EC_IP3_coupling %d %d float\n", 1,
			grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * branches);

	count = header_offset[9];
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3] + header_offset[7]
			+ ecDataOffset[4] + header_offset[8] + ecDataOffset[5]) * sizeof(char);
	if (my_IO_domain_info->writer_rank == 0) {
		CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, header, count, MPI_CHAR, &status));
	}
	CHECK_MPI_ERROR(MPI_Allgather(&writer_buffer->ec_cpl[cpl_IP3], 1, MPI_INT, &buffer_lengths[1], 1, MPI_INT, my_IO_domain_info->writer_comm));

	disp = 0;
	for (int j = my_IO_domain_info->writer_rank; j > 0; j--) {
		disp = disp + buffer_lengths[j];
	}
	disp = (header_offset[0] + point_offset + header_offset[1] + cell_offset + header_offset[2] + celltype_offset + header_offset[3] + ecDataOffset[0]
			+ header_offset[4] + ecDataOffset[1] + header_offset[5] + ecDataOffset[2] + header_offset[6] + ecDataOffset[3] + header_offset[7]
			+ ecDataOffset[4] + header_offset[8] + ecDataOffset[5] + header_offset[9] + disp) * sizeof(char);
	count = writer_buffer->ec_cpl[cpl_IP3];
	CHECK_MPI_ERROR(MPI_File_write_at(check->ec_data_file, disp, writer_buffer->cpIj, count, MPI_CHAR, &status));
	for (int i = 0; i < 4; i++) {
		ecDataOffset[6] += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->cpIj);

	printf("[%d] ++++++>>>>>> Leaving %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);
}

/**
 * Dump info/debug output to a log file.
 */
void dump_rank_info(checkpoint_handle* check, conductance cpl_cef, grid_parms grid, IO_domain_info* my_IO_domain_info)
{
	MPI_Status status;
	MPI_Offset displacement = 0;
	char* buffer = (char*) checked_malloc(2 * 1024 * sizeof(char), "allocation for logfile segment space\n");
	int root = 0;
	char filename[50];
	int length =
			sprintf(buffer,
					"BRANCH_TAG	= %d\n[Universal_Rank, Cart_Rank= (%d,%d)] \tcoords= %d,%d\t nbrs: local (u,d,l,r)=(%d %d %d %d)\t "
							"remote: (up,down)=(%d %d)\n\n flip_array: (%d,%d,%d,%d)\n\n"
							"Boundary_tag = %c\n(T = Top\t B= Bottom\t I=Interior edges of the bifurcation segments, parent or children\t N=Interior of the subdomain)\n"
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
					grid.nbrs[local][LEFT], grid.nbrs[local][RIGHT], grid.nbrs[remote][UP], grid.nbrs[remote][DOWN],
					grid.flip_array[0], grid.flip_array[1], grid.flip_array[2], grid.flip_array[3],
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
			"Allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers.");
	int *disp = (int*) checked_malloc(grid.tasks * sizeof(int), "Allocation failed for disp array in gather_tasks_mesh_point_data_on_writers.");

	// Gathering and summing the length of all the CHARs contained in every send_buffer containing coordinates from each MPI process.
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid.cart_comm));

	grid.logfile_displacements = 0;
	for (int i = 0; i < grid.tasks; i++)
	{
		disp[i] = grid.logfile_displacements;
		grid.logfile_displacements += recv_count[i];
	}

	if (grid.rank == 0)
	{
		grid.logfile_write_buffer = (char*) checked_malloc(grid.logfile_displacements * sizeof(char),
				"Allocation error for writer_buffer for log file.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(buffer, length, MPI_CHAR, grid.logfile_write_buffer, recv_count, disp, MPI_CHAR, root, grid.cart_comm));

	if (grid.rank == 0)
	{
		sprintf(filename, "Logfile_%s.txt", grid.suffix);
		CHECK_MPI_ERROR(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->logptr));
		CHECK_MPI_ERROR(MPI_File_write_at(check->logptr, displacement, grid.logfile_write_buffer, grid.logfile_displacements, MPI_CHAR, &status));
		MPI_File_close(&check->logptr);
		free(grid.logfile_write_buffer);
	}

	free(recv_count);
	free(disp);
}

#if 0
void dump_coords(grid_parms grid, EC_cell** ec, checkpoint_handle* check, const char* message) {

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
	CHECK_MPI_ERROR(MPI_File_write_at(check->coords, disp, &buffer, write_element_count, MPI_DOUBLE, &status));
}
#endif

#if 0
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

	CHECK_MPI_ERROR(MPI_File_write_at_all(check->time_profiling, disp_write, &buffer[0], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->async_calls, disp_write, &buffer[1], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->async_wait, disp_write, &buffer[2], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->barrier_before_comm, disp_write, &buffer[3], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->map_function, disp_write, &buffer[4], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->single_cell_fluxes, disp_write, &buffer[5], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->coupling_fluxes, disp_write, &buffer[6], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->solver, disp_write, &buffer[7], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->writer_func, disp_write, &buffer[8], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->derivative_calls, disp_write, &buffer[9], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->itter_count, disp_write, &buffer[10], 1, MPI_DOUBLE, &status));
	/// Write Comms time profiling data...
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->remote_async_calls, disp_write, &buffer[11], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->remote_async_wait, disp_write, &buffer[12], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->send_buf_update, disp_write, &buffer[13], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->recv_buf_update, disp_write, &buffer[14], 1, MPI_DOUBLE, &status));
	CHECK_MPI_ERROR(MPI_File_write_at_all(check->total_comms_cost, disp_write, &buffer[15], 1, MPI_DOUBLE, &status));
}
#endif

#if 0
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
#endif

void final_checkpoint(checkpoint_handle *check, grid_parms grid)
{
	MPI_Barrier(grid.universe);
	close_common_checkpoints(check);
//	close_time_wise_checkpoints(check);
//	close_time_profiling_checkpoints(check);
}

void close_common_checkpoints(checkpoint_handle* check)
{
	MPI_File_close(&check->logptr);
	MPI_File_close(&check->elapsed_time);
	MPI_File_close(&check->jplc);
	MPI_File_close(&check->coords);
}

void close_time_wise_checkpoints(checkpoint_handle* check)
{
	MPI_File_close(&check->smc_data_file);
	MPI_File_close(&check->ec_data_file);
}

#if 0
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
#endif

#if 0
int checkpoint(checkpoint_handle* check, grid_parms grid, double* tnow, double* y, SMC_cell** smc, EC_cell** ec) {
/// After when the MPI_IO files have been opened, check whether their current instance is first or did they previously existed.
/// This is checked by retrieving the file size of the file recording line number of the the timefile.
/// If the file is empty, then it is assumed to be the first instance of simulation (starting from t=0)
/// otherwise the line number indicates where the last complete result was written in the file and is used as
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
	}
	else if ((line_number == 0) || (line_number == NULL))
	{
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
#endif

/// Read data from the domain_info.txt to retrieve the information related to how the domain is set up.
/// All tasks open the same file to read.
/// Every task has the same displacement so each start to read data from the same position.
/// Each task decides the read buffer size to be allocated by looking at the file size of the file opened for read operation.
/// The data is read and sorted into delimited strings; numbers and stored into corresponding place holders in the array
/// domains[][] in the structure grid_parms grid.
///
///	For a bifurcation as well as a straight segment there's only one domain.
///
/// The arrays in the domains[][] store the following information:
///
/// Element 0: 	Key_val or serial number of the subdomain
/// Element 1:	Subdomain Type (2 possibilities and their values)
/// 				1. Straight Segment (STRSEG): (0)
/// 				2. Bifurcation	(BIF): (1)
/// Element 2: Number of quads/tasks in the axial extent for the current key_val.
/// Element 3: Number of quads/tasks in the circumferential for the current key_val.
/// Element 4: Parent subdomain key_val of current key_val.
/// Element 5: Left child subdomain key_val of the current key_val.
/// Element 6: Right child subdomain key_val of the current key_val.
/// Element 7: Required number of endothelial cells axially per quad/task.
/// Element 8: Required number of smooth muscle cells circumferentially per quad/task.
///
/// In the case of elements 5 and 6, if subdomain type of current key_val is a straight segment,
/// left child is positive or zero, and right child is negative. ???
/// If subdomain type of current key_val is a bifurcation, then both right and left child subdomains are non-negative. ???

int read_config_file(int rank, char* filename, grid_parms* grid) {
	int err;
	MPI_File input_file;
	MPI_Offset file_size;
	MPI_Status status;
	const char *delimiter = ";,\n";
	char *buffer;
	char *token;
	int *values;

	CHECK_MPI_ERROR(MPI_File_open(grid->universe, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &input_file));

	err = MPI_File_get_size(input_file, &file_size);

	buffer = (char*)checked_malloc((int)file_size * sizeof(char) + 1, SRC_LOC); // Extra char is for the null character.

	err = MPI_File_read_all(input_file, buffer, file_size, MPI_CHAR, &status);
	buffer[file_size] = '\0'; // Null-terminate the damn string!

	// Parse the first value indicating the number of domains.
	token = strtok(buffer, delimiter);

	if(token != NULL)
	{
		grid->num_domains = atoi(token);
	}
	else
	{
		printf("[%d] Unable to read the number of domains value from file %s.\n", grid->universal_rank, filename);
		MPI_Abort(MPI_COMM_WORLD, 911);
	}

	// Allocate the array for the rest of the values in the config file.
	values = (int*)checked_malloc(NUM_CONFIG_ELEMENTS * grid->num_domains * sizeof(int), SRC_LOC);

	// Parse the rest of the values in the config file.
	int index = 0;
	token = strtok(NULL, delimiter);
	while(token != NULL)
	{
		values[index++] = atoi(token);
		token = strtok(NULL, delimiter);
	}

	// Error checking.
	if(index != (NUM_CONFIG_ELEMENTS * grid->num_domains))
	{
		printf("[%d] Insufficient number of values in the config file %s.\n", grid->universal_rank, filename);
		MPI_Abort(MPI_COMM_WORLD, 911);
	}

	// TODO: Define a function for allocating 2D arrays.
	// TODO: Check the allocated memory is released when appropriate.
	// TODO: Define a function for releasing 2D arrays.

	// Allocate first dimension array.
	grid->domains = (int**) checked_malloc(grid->num_domains * sizeof(int*), SRC_LOC);

	// Allocate second dimension arrays.
	for (int i = 0; i < grid->num_domains; i++) {
		grid->domains[i] = (int*) checked_malloc(NUM_CONFIG_ELEMENTS * sizeof(int), SRC_LOC);
	}

	// Copy the data into the domains array from the values array.
	for (int i = 0; i < grid->num_domains; i++) {
		for (int j = 0; j < NUM_CONFIG_ELEMENTS; j++) {
			grid->domains[i][j] = values[(i * NUM_CONFIG_ELEMENTS) + j];
		}
	}

	MPI_File_close(&input_file);

	free(buffer);
	free(values);

	return 0;
}

// Every cylinder root writes elapsed time to a file.
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
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid.cart_comm));
	int total_buffer_length = 0;
	for (int i = 0; i < grid.tasks; i++) {
		disp[i] = total_buffer_length;
		total_buffer_length += recv_count[i];
	}
	if (grid.rank == 0) {
		write_buffer = (char*) checked_malloc(total_buffer_length * sizeof(char), "allocation error for writer_buffer for Elapsed_time file.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(buffer, length, MPI_CHAR, write_buffer, recv_count, disp, MPI_CHAR, root, grid.cart_comm));
	if (grid.rank == 0) {
		sprintf(filename, "Elapsed_time_%s.txt", grid.suffix);
		CHECK_MPI_ERROR(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->elapsed_time));
		CHECK_MPI_ERROR(MPI_File_write_at(check->elapsed_time, 0, write_buffer, total_buffer_length, MPI_CHAR, &status));
		MPI_File_close(&check->elapsed_time);
		free(write_buffer);
	}

	free(recv_count);
	free(disp);
}

// Prepare the suffix which indicates our subdomain information, bifurcation or tube segment suffix, and the containing branch info.
void set_file_naming_strings(grid_parms* grid)
{
	int subdomain, branch, err;

	if(grid->my_domain.internal_info.domain_type == STRSEG)
	{
		subdomain = grid->my_domain.internal_info.domain_index;
		err = sprintf(grid->suffix, "%d", subdomain);
	}
	else if(grid->my_domain.internal_info.domain_type == BIF)
	{
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
	/*CHECK_MPI_ERROR(
	 MPI_File_open(grid.cart_comm, filename, MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &check->itter_count));
	 */
	CHECK_MPI_ERROR(MPI_File_read_at(check->itter_count, disp, &file_offset, 1, MPI_INT, &status));

	return (file_offset);

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

/// Read coordinates from relevant geometry files for each ec and smc in the computational domain.
int read_topology_info(char* filename, grid_parms* grid, SMC_cell **smc, EC_cell **ec)
{
	int buffer[24];
	int msg_tag_1 = 1, msg_tag_2 = 2, msg_tag_3 = 3, msg_tag_4 = 4;
	MPI_Status status;

	if(grid->universal_rank == 0)
	{
		FILE *fr, *fw;
		fr = fopen(filename, "r+");
		if(fr == NULL)
		{
			printf("[%d] Unable to open %s for reading.\n", grid->universal_rank, filename);
			MPI_Abort(grid->universe, 101);
		}
		int p = 0;
		while(fgetc(fr) != EOF)
		{
			int err = fscanf(fr, "%d", &buffer[p]);
			if (err > 0)
				p++;
		}
		fclose(fr);
	}

	CHECK_MPI_ERROR(MPI_Bcast(buffer, 24, MPI_INT, 0, grid->universe));

	/// These kinds of meshes are to be read. For each mesh type the points and cells are to be read
	/// into a structure vtk_info. Using the retrieved data the root will send coordinates to all processes,
	/// for their specific boundaries, and their constituent SMCs and ECs.

	/// The following ints are used as macros below:

	int branch;
	if (grid->my_domain.internal_info.domain_type == STRSEG)
	{
		branch = P;
	}
	else if (grid->my_domain.internal_info.domain_type == BIF)
	{
		branch = grid->branch_tag;
	}

	grid->info = (int**)checked_malloc(4 * sizeof(int*), "Memory allocation for info failed.");
	for (int i = 0; i < 4; i++)
	{
		grid->info[i] = (int*)checked_malloc(6 * sizeof(int), "Memory allocation for info failed.");
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 6; j++) {
			grid->info[i][j] = buffer[i * 6 + j];

		}
	}

	int *disp, *send_count, recv_count, root = 0;
	double *send_points, *recv_points;

	int num_tuple_components = 3, num_tuples = 4;
	int tuple_offset = num_tuple_components * num_tuples;

	grid->coordinates = (double**)checked_malloc(num_tuples * sizeof(double*), SRC_LOC);

	for(int i = 0; i < num_tuples; i++)
	{
		grid->coordinates[i] = (double*)checked_malloc(num_tuple_components * sizeof(double), SRC_LOC);
	}

	send_count = (int*)checked_malloc(grid->info[PROCESS_MESH][TOTAL_CELLS] * sizeof(int), SRC_LOC);

	disp = (int*)checked_malloc(grid->info[PROCESS_MESH][TOTAL_CELLS] * sizeof(int), SRC_LOC);

	send_points = (double*)checked_malloc(tuple_offset * grid->info[PROCESS_MESH][TOTAL_CELLS] * sizeof(double), SRC_LOC);

	if (grid->rank == 0)
	{
		vtk_info* process_mesh = (vtk_info*)checked_malloc(sizeof(vtk_info), SRC_LOC);
		process_mesh->points = (double**)checked_malloc(grid->info[PROCESS_MESH][TOTAL_POINTS] * sizeof(double*), SRC_LOC);
		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_POINTS]; i++)
		{
			process_mesh->points[i] = (double*)checked_malloc(3 * sizeof(double), SRC_LOC);
		}
		process_mesh->cells = (int**)checked_malloc(grid->info[PROCESS_MESH][TOTAL_CELLS] * sizeof(int*), SRC_LOC);
		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS]; i++)
		{
			process_mesh->cells[i] = (int*)checked_malloc(5 * sizeof(int), SRC_LOC);
		}

		int read_counts[2] = {0, 0};
		read_coordinates(grid->info, process_mesh, branch, PROCESS_MESH,
				grid->info[PROCESS_MESH][TOTAL_POINTS],
				grid->info[PROCESS_MESH][TOTAL_CELLS],
				read_counts);

		assert(read_counts[0] == grid->info[PROCESS_MESH][TOTAL_POINTS] &&
				read_counts[1] == grid->info[PROCESS_MESH][TOTAL_CELLS]);

		// Do we really want to iterate over the cells here?
		for(int i = 0; i < read_counts[1]; i++)
		{
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

		for(int i = 0; i < grid->info[PROCESS_MESH][TOTAL_POINTS]; i++)
		{
			free(process_mesh->points[i]);
		}
		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS]; i++)
		{
			free(process_mesh->cells[i]);
		}
		free(process_mesh->points);
		free(process_mesh->cells);
		free(process_mesh);
	}

	recv_points = (double*)checked_malloc(tuple_offset * sizeof(double), SRC_LOC);

	for(int i = 0; i < grid->tasks; i++)
	{
		send_count[i] = tuple_offset;
		disp[i] = tuple_offset * i;
	}

	recv_count = tuple_offset;

	CHECK_MPI_ERROR(MPI_Scatterv(send_points, send_count, disp, MPI_DOUBLE, recv_points, recv_count, MPI_DOUBLE, root, grid->cart_comm));

	for(int i = 0; i < num_tuples; i++)
	{
		grid->coordinates[i][0] = recv_points[(i * num_tuple_components) + 0];
		grid->coordinates[i][1] = recv_points[(i * num_tuple_components) + 1];
		grid->coordinates[i][2] = recv_points[(i * num_tuple_components) + 2];
	}

	free(send_points);
	free(recv_points);

	// Reading in SMC mesh and communicating to each branch communicator member.
	send_points = (double*)checked_malloc(tuple_offset * grid->info[SMC_MESH][TOTAL_CELLS] * grid->info[PROCESS_MESH][TOTAL_CELLS] * sizeof(double), SRC_LOC);

	if(grid->rank == 0)
	{

		int num_SMC_points = grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_POINTS];
		int num_SMC_cells = grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[SMC_MESH][TOTAL_CELLS];

		vtk_info* smc_mesh = (vtk_info*)checked_malloc(sizeof(vtk_info), SRC_LOC);

		smc_mesh->points = (double**)checked_malloc(num_SMC_points * sizeof(double*), SRC_LOC);
		for(int i = 0; i < num_SMC_points; i++)
		{
			smc_mesh->points[i] = (double*)checked_malloc(3 * sizeof(double), SRC_LOC);
		}

		smc_mesh->cells = (int**)checked_malloc(num_SMC_cells * sizeof(int*), SRC_LOC);
		for(int i = 0; i < num_SMC_cells; i++)
		{
			smc_mesh->cells[i] = (int*)checked_malloc(5 * sizeof(int), SRC_LOC);
		}

		int read_counts[2] = {0, 0};
		read_coordinates(grid->info, smc_mesh, branch, SMC_MESH,
				(num_SMC_points),
				(num_SMC_cells),
				read_counts);

		assert(read_counts[0] == (num_SMC_points) && read_counts[1] == (num_SMC_cells));

		for(int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS]; i++)
		{
			for(int j = 0; j < grid->info[SMC_MESH][TOTAL_CELLS]; j++)
			{
				// Declare two local variables to reduce the number of the initial additions and multiplications
				// to 1/12th because these variables can be used in the next 12 lines. Could reduce it further.
				int sp_off = (i * grid->info[SMC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset);
				int cell_pos = (i * grid->info[SMC_MESH][TOTAL_CELLS]) + j;

				// Sanity check.
				assert(cell_pos < num_SMC_cells);

				int pt_id_1 = smc_mesh->cells[cell_pos][1];
				int pt_id_2 = smc_mesh->cells[cell_pos][2];
				int pt_id_3 = smc_mesh->cells[cell_pos][3];
				int pt_id_4 = smc_mesh->cells[cell_pos][4];

				// Sanity check.
				assert(pt_id_1 < num_SMC_points);
				assert(pt_id_2 < num_SMC_points);
				assert(pt_id_3 < num_SMC_points);
				assert(pt_id_4 < num_SMC_points);

				send_points[sp_off + 0] = smc_mesh->points[pt_id_1][0];
				send_points[sp_off + 1] = smc_mesh->points[pt_id_1][1];
				send_points[sp_off + 2] = smc_mesh->points[pt_id_1][2];

				send_points[sp_off + 3] = smc_mesh->points[pt_id_2][0];
				send_points[sp_off + 4] = smc_mesh->points[pt_id_2][1];
				send_points[sp_off + 5] = smc_mesh->points[pt_id_2][2];

				send_points[sp_off + 6] = smc_mesh->points[pt_id_3][0];
				send_points[sp_off + 7] = smc_mesh->points[pt_id_3][1];
				send_points[sp_off + 8] = smc_mesh->points[pt_id_3][2];

				send_points[sp_off + 9] = smc_mesh->points[pt_id_4][0];
				send_points[sp_off + 10] = smc_mesh->points[pt_id_4][1];
				send_points[sp_off + 11] = smc_mesh->points[pt_id_4][2];
			}
		}

		for(int i = 0; i < num_SMC_points; i++)
		{
			free(smc_mesh->points[i]);
		}

		for(int i = 0; i < num_SMC_cells; i++)
		{
			free(smc_mesh->cells[i]);
		}
		free(smc_mesh->points);
		free(smc_mesh->cells);
		free(smc_mesh);
	}

	recv_points = (double*)checked_malloc(tuple_offset * grid->info[SMC_MESH][TOTAL_CELLS] * sizeof(double), SRC_LOC);
	for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS]; i++)
	{
		send_count[i] = tuple_offset * grid->info[SMC_MESH][TOTAL_CELLS];
		disp[i] = tuple_offset * grid->info[SMC_MESH][TOTAL_CELLS] * i;
	}
	recv_count = tuple_offset * grid->info[SMC_MESH][TOTAL_CELLS];

	CHECK_MPI_ERROR(MPI_Scatterv(send_points, send_count, disp, MPI_DOUBLE, recv_points, recv_count, MPI_DOUBLE, root, grid->cart_comm));

	int count = 0;
	for(int n = 1; n <= grid->num_smc_axially; n++)
	{
		for(int m = 1; m <= grid->num_smc_circumferentially; m++)
		{
			for(int l = 0; l < 4; l++)
			{
				smc[m][n].x_coordinate[l] = recv_points[count * tuple_offset + l * 3 + 0];
				smc[m][n].y_coordinate[l] = recv_points[count * tuple_offset + l * 3 + 1];
				smc[m][n].z_coordinate[l] = recv_points[count * tuple_offset + l * 3 + 2];
			}
			count++;
		}
	}
	free(send_points);
	free(recv_points);

	// Reading in EC mesh and communicating to each branch communicator member.
	send_points = (double*)checked_malloc(tuple_offset * grid->info[EC_MESH][TOTAL_CELLS] * grid->info[PROCESS_MESH][TOTAL_CELLS] * sizeof(double), SRC_LOC);

	if (grid->rank == 0) {
		vtk_info* ec_mesh = (vtk_info*)checked_malloc(sizeof(vtk_info), SRC_LOC);

		ec_mesh->points = (double**)checked_malloc(grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_POINTS] * sizeof(double*), SRC_LOC);
		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_POINTS]; i++)
		{
			ec_mesh->points[i] = (double*)checked_malloc(3 * sizeof(double), SRC_LOC);
		}

		ec_mesh->cells = (int**)checked_malloc(grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS] * sizeof(int*), SRC_LOC);
		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS]; i++)
		{
			ec_mesh->cells[i] = (int*)checked_malloc(5 * sizeof(int), SRC_LOC);
		}

		int read_counts[2] = {0, 0};
		read_coordinates(grid->info, ec_mesh, branch, EC_MESH,
				(grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_POINTS]),
				(grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS]),
				read_counts);

		assert(read_counts[0] == (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_POINTS]) &&
				read_counts[1] == (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS]));

		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS]; i++) {
			for (int j = 0; j < grid->info[EC_MESH][TOTAL_CELLS]; j++) {
				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 0] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][1]][0];
				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 1] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][1]][1];
				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 2] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][1]][2];

				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 3] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][2]][0];
				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 4] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][2]][1];
				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 5] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][2]][2];

				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 6] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][3]][0];
				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 7] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][3]][1];
				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 8] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][3]][2];

				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 9] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][4]][0];
				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 10] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][4]][1];
				send_points[(i * grid->info[EC_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 11] = ec_mesh->points[ec_mesh->cells[(i
						* grid->info[EC_MESH][TOTAL_CELLS]) + j][4]][2];
			}
		}
		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_POINTS]; i++) {
			free(ec_mesh->points[i]);
		}
		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_MESH][TOTAL_CELLS]; i++) {
			free(ec_mesh->cells[i]);
		}
		free(ec_mesh->points);
		free(ec_mesh->cells);
		free(ec_mesh);
	}

	recv_points = (double*)checked_malloc(tuple_offset * grid->info[EC_MESH][TOTAL_CELLS] * sizeof(double), SRC_LOC);
	for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS]; i++)
	{
		send_count[i] = tuple_offset * grid->info[EC_MESH][TOTAL_CELLS];
		disp[i] = tuple_offset * grid->info[EC_MESH][TOTAL_CELLS] * i;
	}
	recv_count = tuple_offset * grid->info[EC_MESH][TOTAL_CELLS];

	CHECK_MPI_ERROR(MPI_Scatterv(send_points, send_count, disp, MPI_DOUBLE, recv_points, recv_count, MPI_DOUBLE, root, grid->cart_comm));

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

	// Reading in EC centroids and communicating to each branch communicator member.
	num_tuples = 1;
	tuple_offset = num_tuple_components * num_tuples;

	send_points = (double*)checked_malloc(tuple_offset * grid->info[EC_CENT_MESH][TOTAL_CELLS] * grid->info[PROCESS_MESH][TOTAL_CELLS] * sizeof(double), SRC_LOC);

	if (grid->rank == 0) {
		vtk_info* ec_centroids = (vtk_info*)checked_malloc(sizeof(vtk_info), SRC_LOC);
		ec_centroids->points = (double**)checked_malloc(grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_CENT_MESH][TOTAL_POINTS] * sizeof(double*), SRC_LOC);
		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_CENT_MESH][TOTAL_POINTS]; i++) {
			ec_centroids->points[i] = (double*)checked_malloc(3 * sizeof(double), SRC_LOC);
		}
		ec_centroids->cells = (int**)checked_malloc(grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_CENT_MESH][TOTAL_CELLS] * sizeof(int*), SRC_LOC);
		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_CENT_MESH][TOTAL_CELLS]; i++) {
			ec_centroids->cells[i] = (int*)checked_malloc(2 * sizeof(int), SRC_LOC);
		}

		int read_counts[2] = {0, 0};
		read_coordinates(grid->info, ec_centroids, branch, EC_CENT_MESH,
				(grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_CENT_MESH][TOTAL_POINTS]),
				(grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_CENT_MESH][TOTAL_CELLS]),
				read_counts);

		assert(read_counts[0] == (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_CENT_MESH][TOTAL_POINTS]) &&
				read_counts[1] == (grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_CENT_MESH][TOTAL_CELLS]));

		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS]; i++) {
			for (int j = 0; j < grid->info[EC_CENT_MESH][TOTAL_CELLS]; j++) {
				send_points[(i * grid->info[EC_CENT_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 0] =
						ec_centroids->points[ec_centroids->cells[(i * grid->info[EC_CENT_MESH][TOTAL_CELLS]) + j][1]][0];
				send_points[(i * grid->info[EC_CENT_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 1] =
						ec_centroids->points[ec_centroids->cells[(i * grid->info[EC_CENT_MESH][TOTAL_CELLS]) + j][1]][1];
				send_points[(i * grid->info[EC_CENT_MESH][TOTAL_CELLS] * tuple_offset) + (j * tuple_offset) + 2] =
						ec_centroids->points[ec_centroids->cells[(i * grid->info[EC_CENT_MESH][TOTAL_CELLS]) + j][1]][2];
			}
		}
		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_CENT_MESH][TOTAL_POINTS]; i++) {
			free(ec_centroids->points[i]);
		}
		for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS] * grid->info[EC_CENT_MESH][TOTAL_CELLS]; i++) {
			free(ec_centroids->cells[i]);
		}
		free(ec_centroids->points);
		free(ec_centroids->cells);
		free(ec_centroids);
	}

	recv_points = (double*)checked_malloc(tuple_offset * grid->info[EC_CENT_MESH][TOTAL_CELLS] * sizeof(double), SRC_LOC);
	for (int i = 0; i < grid->info[PROCESS_MESH][TOTAL_CELLS]; i++) {
		send_count[i] = tuple_offset * grid->info[EC_CENT_MESH][TOTAL_CELLS];
		disp[i] = tuple_offset * grid->info[EC_CENT_MESH][TOTAL_CELLS] * i;
	}

	recv_count = tuple_offset * grid->info[EC_CENT_MESH][TOTAL_CELLS];
	CHECK_MPI_ERROR(MPI_Scatterv(send_points, send_count, disp, MPI_DOUBLE, recv_points, recv_count, MPI_DOUBLE, root, grid->cart_comm));

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

	free(disp);
	free(send_count);

	MPI_Barrier(grid->universe);

	return (0);
}

void read_init_ATP(grid_parms *grid, EC_cell **ECs)
{
	int branch;
	if (grid->my_domain.internal_info.domain_type == STRSEG)
	{
		branch = P;
	}
	else if (grid->my_domain.internal_info.domain_type == BIF)
	{
		branch = grid->branch_tag;
	}

	int jplc_per_task_count = grid->num_ec_circumferentially * grid->num_ec_axially;

	int jplc_in_size = jplc_per_task_count * grid->tasks;

	// This can be allocated only on the IO nodes.
	double *send_jplc = (double *)checked_malloc(jplc_in_size * sizeof(double), SRC_LOC);

	// Only the IO nodes read the input files.
	if (grid->rank == 0)
	{
		char jplc_file_name[64];
		switch(branch)
		{
			case P:
				sprintf(jplc_file_name, "files/parent_atp.txt");
				break;
			case L:
				sprintf(jplc_file_name, "files/left_daughter_atp.txt");
				break;
			case R:
				sprintf(jplc_file_name, "files/right_daughter_atp.txt");
				break;
			default:
				; // Do something sensible here otherwise all hell breaks loose...
		}

		FILE *fr = fopen(jplc_file_name, "r+");
		printf("Reading ATP from %s, FILE is %s...\n", jplc_file_name, fr == NULL ? "NULL" : "OK");

		int count_in = 0;
		double jplc_val = 0;
		while(fscanf(fr, "%lf", &jplc_val) == 1)
		{
			send_jplc[count_in] = jplc_val;
			count_in++;
		}
		if(feof(fr))
		{
			// Check the number of values read.
			assert(jplc_in_size == count_in);
		}
		fclose(fr);
	}

	int *send_jplc_counts = (int *)checked_malloc(grid->tasks * sizeof(int), SRC_LOC);
	int *send_jplc_offsets = (int *)checked_malloc(grid->tasks * sizeof(int), SRC_LOC);

	for(int task = 0; task < grid->tasks; task++)
	{
		send_jplc_counts[task] = jplc_per_task_count;
		send_jplc_offsets[task] = task * jplc_per_task_count;
	}

	int recv_jplc_count = jplc_per_task_count;
	double *recv_jplc = (double *)checked_malloc(recv_jplc_count * sizeof(double), SRC_LOC);

	int root = 0;

	// printf("%s, grid->cart_comm: %p\n", __FUNCTION__, (void *)grid->cart_comm);

	// Scatter JPLC values to the nodes in this Cartesian grid.
	CHECK_MPI_ERROR(MPI_Scatterv(send_jplc, send_jplc_counts, send_jplc_offsets, MPI_DOUBLE, recv_jplc, recv_jplc_count, MPI_DOUBLE, root, grid->cart_comm));

	// Assign received JPLC values to the cells.
	for(int m = 1; m <= grid->num_ec_circumferentially; m++)
	{
		for(int n = 1; n <= grid->num_ec_axially; n++)
		{
			// Fortran array referencing!
			ECs[m][n].JPLC = recv_jplc[(n - 1) * grid->num_ec_circumferentially + m - 1];
		}
	}

	free(send_jplc);
	free(send_jplc_counts);
	free(send_jplc_offsets);
	free(recv_jplc);
}

void read_coordinates(int** info, vtk_info* mesh, int branch, int mesh_type, int points, int cells, int *read_counts)
{
	FILE *fr;
	char filename_points[64], filename_cells[64];
	if (mesh_type == 0) {
		if (branch == P) {
			sprintf(filename_points, "files/parent_points.txt");
			sprintf(filename_cells, "files/parent_cells.txt");
		} else if (branch == L) {
			sprintf(filename_points, "files/left_daughter_points.txt");
			sprintf(filename_cells, "files/left_daughter_cells.txt");
		} else if (branch == R) {
			sprintf(filename_points, "files/right_daughter_points.txt");
			sprintf(filename_cells, "files/right_daughter_cells.txt");
		}
	} else if (mesh_type == 1) {
		if (branch == P) {
			sprintf(filename_points, "files/parent_smc_mesh_points.txt");
			sprintf(filename_cells, "files/parent_smc_mesh_cells.txt");
		} else if (branch == L) {
			sprintf(filename_points, "files/left_daughter_smc_mesh_points.txt");
			sprintf(filename_cells, "files/left_daughter_smc_mesh_cells.txt");
		} else if (branch == R) {
			sprintf(filename_points, "files/right_daughter_smc_mesh_points.txt");
			sprintf(filename_cells, "files/right_daughter_smc_mesh_cells.txt");
		}
	} else if (mesh_type == 2) {
		if (branch == P) {
			sprintf(filename_points, "files/parent_ec_mesh_points.txt");
			sprintf(filename_cells, "files/parent_ec_mesh_cells.txt");
		} else if (branch == L) {
			sprintf(filename_points, "files/left_daughter_ec_mesh_points.txt");
			sprintf(filename_cells, "files/left_daughter_ec_mesh_cells.txt");
		} else if (branch == R) {
			sprintf(filename_points, "files/right_daughter_ec_mesh_points.txt");
			sprintf(filename_cells, "files/right_daughter_ec_mesh_cells.txt");
		}
	} else if (mesh_type == 3) {
		if (branch == P) {
			sprintf(filename_points, "files/parent_ec_centeroid_points.txt");
			sprintf(filename_cells, "files/parent_ec_centeroid_cells.txt");
		} else if (branch == L) {
			sprintf(filename_points, "files/left_daughter_ec_centeroid_points.txt");
			sprintf(filename_cells, "files/left_daughter_ec_centeroid_cells.txt");
		} else if (branch == R) {
			sprintf(filename_points, "files/right_daughter_ec_centeroid_points.txt");
			sprintf(filename_cells, "files/right_daughter_ec_centeroid_cells.txt");
		}
	}

	fr = fopen(filename_points, "r+");
	printf("Reading points from %s, FILE is %s\n", filename_points, fr == NULL ? "NULL" : "OK");
	for (int i = 0; i < points; i++) {
		fscanf(fr, "%lf", &mesh->points[i][0]);
		fscanf(fr, "%lf", &mesh->points[i][1]);
		fscanf(fr, "%lf", &mesh->points[i][2]);
		read_counts[0]++;
	}
	fclose(fr);

	fr = fopen(filename_cells, "r+");
	printf("Reading cells from %s, FILE is %s\n", filename_cells, fr == NULL ? "NULL" : "OK");
	if (mesh_type < 3) {
		for (int i = 0; i < cells; i++) {
			fscanf(fr, "%d", &mesh->cells[i][0]);
			fscanf(fr, "%d", &mesh->cells[i][1]);
			fscanf(fr, "%d", &mesh->cells[i][2]);
			fscanf(fr, "%d", &mesh->cells[i][3]);
			fscanf(fr, "%d", &mesh->cells[i][4]);
			read_counts[1]++;
		}
	} else if (mesh_type == 3) {
		for (int i = 0; i < cells; i++) {
			fscanf(fr, "%d", &mesh->cells[i][0]);
			fscanf(fr, "%d", &mesh->cells[i][1]);
			read_counts[1]++;
		}
	}
	fclose(fr);
}

void gather_tasks_mesh_point_data_on_writers(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, SMC_cell** smc,
		EC_cell** ec)
{

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

	char *send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");
	int length = 0;
	for (int i = 0; i < num_tuples; i++) {
		length += sprintf(send_buffer + length, "%2.9lf\t%2.9lf\t%2.9lf\n", grid->coordinates[i][0], grid->coordinates[i][1],
				grid->coordinates[i][2]);
	}

	/// Gathering and summing the length of all the CHARs contained in every send_buffer containing coordinates from each MPI process.
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));
	writer_buffer->buffer_length[PROCESS_MESH] = 0;
	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[PROCESS_MESH];
		writer_buffer->buffer_length[PROCESS_MESH] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->process_mesh_points = (char*) checked_malloc(writer_buffer->buffer_length[PROCESS_MESH] * sizeof(char),
				"allocation error for writer_buffer member process mesh.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->process_mesh_points, recv_count, disp, MPI_CHAR, root, grid->cart_comm));
	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* Cell data *************/
	num_tuple_components = 5;
	num_tuples = 1;
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");

	length = sprintf(send_buffer, "%d %d %d %d %d\n", 4, (branch * grid->tasks * 4) + 4 * grid->rank + 0,
			(branch * grid->tasks * 4) + 4 * grid->rank + 1, (branch * grid->tasks * 4) + 4 * grid->rank + 2,
			(branch * grid->tasks * 4) + 4 * grid->rank + 3);

	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));
	writer_buffer->buffer_length[ProcessCell] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ProcessCell];
		writer_buffer->buffer_length[ProcessCell] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->process_mesh_cells = (char*) checked_malloc(writer_buffer->buffer_length[ProcessCell] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->process_mesh_cells, recv_count, disp, MPI_CHAR, root, grid->cart_comm));
	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* Cell type data *************/
	num_tuple_components = 1;
	num_tuples = 1;
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");

	length = sprintf(send_buffer, "%d\n", 9);

	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));
	writer_buffer->buffer_length[ProcessCellType] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ProcessCellType];
		writer_buffer->buffer_length[ProcessCellType] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->process_mesh_type = (char*) checked_malloc(writer_buffer->buffer_length[ProcessCellType] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->process_mesh_type, recv_count, disp, MPI_CHAR, root, grid->cart_comm));
	free(recv_count);
	free(send_buffer);
	free(disp);
}

void gather_smc_mesh_data_on_writers(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, SMC_cell** smc) {
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

	char *send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
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
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));
	writer_buffer->buffer_length[SMC_MESH] = 0;
	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[SMC_MESH];
		writer_buffer->buffer_length[SMC_MESH] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->smc_mesh_points = (char*) checked_malloc(writer_buffer->buffer_length[SMC_MESH] * sizeof(char),
				"allocation error for writer_buffer member process mesh.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->smc_mesh_points, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* Cell data *************/
	int my_branch_smc_offset = 0;
	num_tuple_components = 5;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");

	int* global_smc_point_offset_info = (int*) checked_malloc(grid->numtasks * sizeof(int),
			"error in global info array allocation in function gather_smc_mesh_data_on_writers");
	int my_smc_point_offset_info = 4 * grid->num_smc_axially * grid->num_smc_circumferentially;
	CHECK_MPI_ERROR(MPI_Allgather(&my_smc_point_offset_info, 1, MPI_INT, global_smc_point_offset_info, 1, MPI_INT, grid->universe));

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

	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));
	writer_buffer->buffer_length[smcCell] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[smcCell];
		writer_buffer->buffer_length[smcCell] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->smc_mesh_cells = (char*) checked_malloc(writer_buffer->buffer_length[smcCell] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->smc_mesh_cells, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* Celltype data *************/
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
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

	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));
	writer_buffer->buffer_length[smcCellType] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[smcCellType];
		writer_buffer->buffer_length[smcCellType] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->smc_mesh_type = (char*) checked_malloc(writer_buffer->buffer_length[smcCellType] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->smc_mesh_type, recv_count, disp, MPI_CHAR, root, grid->cart_comm));
	free(recv_count);
	free(send_buffer);
	free(disp);
}

void gather_ec_mesh_data_on_writers(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, EC_cell** ec) {
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

	char *send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
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
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));
	writer_buffer->buffer_length[EC_MESH] = 0;
	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[EC_MESH];
		writer_buffer->buffer_length[EC_MESH] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->ec_mesh_points = (char*) checked_malloc(writer_buffer->buffer_length[EC_MESH] * sizeof(char),
				"allocation error for writer_buffer member process mesh.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->ec_mesh_points, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* Cell data *************/
	int my_branch_ec_offset = 0;
	num_tuple_components = 5;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");

	int* global_ec_point_offset_info = (int*) checked_malloc(grid->numtasks * sizeof(int),
			"error in global info array allocation in function gather_ec_mesh_data_on_writers");
	int my_ec_point_offset_info = 4 * grid->num_ec_axially * grid->num_ec_circumferentially;
	CHECK_MPI_ERROR(MPI_Allgather(&my_ec_point_offset_info, 1, MPI_INT, global_ec_point_offset_info, 1, MPI_INT, grid->universe));

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

	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));
	writer_buffer->buffer_length[ecCell] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ecCell];
		writer_buffer->buffer_length[ecCell] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->ec_mesh_cells = (char*) checked_malloc(writer_buffer->buffer_length[ecCell] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->ec_mesh_cells, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* Celltype data *************/
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
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

	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));
	writer_buffer->buffer_length[ecCellType] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ecCellType];
		writer_buffer->buffer_length[ecCellType] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->ec_mesh_type = (char*) checked_malloc(writer_buffer->buffer_length[ecCellType] * sizeof(char),
				"allocation error for writer_buffer member process cells.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->ec_mesh_type, recv_count, disp, MPI_CHAR, root, grid->cart_comm));
	free(recv_count);
	free(send_buffer);
	free(disp);
}

void gather_smcData(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, SMC_cell** smc, int write_count)
{
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

	char* send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_smcData.");

	int length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].p[smc_Ca]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->smc_stat_var_buffer_length[smc_Ca] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_stat_var_buffer_length[smc_Ca];
		writer_buffer->smc_stat_var_buffer_length[smc_Ca] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->ci = (char*) checked_malloc(writer_buffer->smc_stat_var_buffer_length[smc_Ca] * sizeof(char),
				"allocation error for writer_buffer member ci.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->ci, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

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

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_smcData.");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].p[smc_SR]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->smc_stat_var_buffer_length[smc_SR] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_stat_var_buffer_length[smc_SR];
		writer_buffer->smc_stat_var_buffer_length[smc_SR] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->si = (char*) checked_malloc(writer_buffer->smc_stat_var_buffer_length[smc_SR] * sizeof(char),
				"allocation error for writer_buffer member si.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->si, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

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

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_smcData.");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].p[smc_Vm]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->smc_stat_var_buffer_length[smc_Vm] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_stat_var_buffer_length[smc_Vm];
		writer_buffer->smc_stat_var_buffer_length[smc_Vm] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->vi = (char*) checked_malloc(writer_buffer->smc_stat_var_buffer_length[smc_Vm] * sizeof(char),
				"allocation error for writer_buffer member vi.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->vi, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

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

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_smcData.");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].p[smc_w]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->smc_stat_var_buffer_length[smc_w] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_stat_var_buffer_length[smc_w];
		writer_buffer->smc_stat_var_buffer_length[smc_w] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->wi = (char*) checked_malloc(writer_buffer->smc_stat_var_buffer_length[smc_w] * sizeof(char),
				"allocation error for writer_buffer member wi.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->wi, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

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

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_smcData.");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].p[smc_IP3]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->smc_stat_var_buffer_length[smc_IP3] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_stat_var_buffer_length[smc_IP3];
		writer_buffer->smc_stat_var_buffer_length[smc_IP3] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->Ii = (char*) checked_malloc(writer_buffer->smc_stat_var_buffer_length[smc_IP3] * sizeof(char),
				"allocation error for writer_buffer member Ii.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->Ii, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

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

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_smcData.");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].B[cpl_Ca]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->smc_cpl[cpl_Ca] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_cpl[cpl_Ca];
		writer_buffer->smc_cpl[cpl_Ca] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpCi = (char*) checked_malloc(writer_buffer->smc_cpl[cpl_Ca] * sizeof(char),
				"allocation error for writer_buffer member smc_cplCa.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpCi, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* smc_cpl_Vm field data *************/
	my_branch_smc_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_smcData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].B[cpl_Vm]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->smc_cpl[cpl_Vm] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_cpl[cpl_Vm];
		writer_buffer->smc_cpl[cpl_Vm] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpVi = (char*) checked_malloc(writer_buffer->smc_cpl[cpl_Vm] * sizeof(char),
				"allocation error for writer_buffer member smc_cplVm.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpVi, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* smc_cpl_IP3 field data *************/
	my_branch_smc_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_smc_axially * grid->num_smc_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_smcData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid->num_smc_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", smc[i][j].B[cpl_IP3]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->smc_cpl[cpl_IP3] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->smc_cpl[cpl_IP3];
		writer_buffer->smc_cpl[cpl_IP3] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpIi = (char*) checked_malloc(writer_buffer->smc_cpl[cpl_IP3] * sizeof(char),
				"allocation error for writer_buffer member smc_cplIi.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpIi, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);
}

void gather_ecData(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, EC_cell** ec, int write_count)
{
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

	char* send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_ecData.");

	int length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].q[ec_Ca]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->ec_stat_var_buffer_length[ec_Ca] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_stat_var_buffer_length[ec_Ca];
		writer_buffer->ec_stat_var_buffer_length[ec_Ca] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cj = (char*) checked_malloc(writer_buffer->ec_stat_var_buffer_length[ec_Ca] * sizeof(char),
				"allocation error for writer_buffer member cj.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cj, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* ecSR field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].q[ec_SR]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->ec_stat_var_buffer_length[ec_SR] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_stat_var_buffer_length[ec_SR];
		writer_buffer->ec_stat_var_buffer_length[ec_SR] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->sj = (char*) checked_malloc(writer_buffer->ec_stat_var_buffer_length[ec_SR] * sizeof(char),
				"allocation error for writer_buffer member sj.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->sj, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* ecVm field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].q[ec_Vm]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->ec_stat_var_buffer_length[ec_Vm] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_stat_var_buffer_length[ec_Vm];
		writer_buffer->ec_stat_var_buffer_length[ec_Vm] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->vj = (char*) checked_malloc(writer_buffer->ec_stat_var_buffer_length[ec_Vm] * sizeof(char),
				"allocation error for writer_buffer member cj.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->vj, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* ec_IP3 field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].q[ec_IP3]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->ec_stat_var_buffer_length[ec_IP3] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_stat_var_buffer_length[ec_IP3];
		writer_buffer->ec_stat_var_buffer_length[ec_IP3] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->Ij = (char*) checked_malloc(writer_buffer->ec_stat_var_buffer_length[ec_Vm] * sizeof(char),
				"allocation error for writer_buffer member Ij.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->Ij, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* ec_cpl_Ca field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].B[cpl_Ca]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->ec_cpl[cpl_Ca] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_cpl[cpl_Ca];
		writer_buffer->ec_cpl[cpl_Ca] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpCj = (char*) checked_malloc(writer_buffer->ec_cpl[cpl_Ca] * sizeof(char),
				"allocation error for writer_buffer member ec_cplCa.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpCj, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* ec_cpl_Vm field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].B[cpl_Vm]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->ec_cpl[cpl_Vm] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_cpl[cpl_Vm];
		writer_buffer->ec_cpl[cpl_Vm] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpVj = (char*) checked_malloc(writer_buffer->ec_cpl[cpl_Vm] * sizeof(char),
				"allocation error for writer_buffer member ec_cplVm.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpVj, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);

	/************* ec_cpl_IP3 field data *************/
	my_branch_ec_offset = 0;
	num_tuple_components = 1;
	num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;

	send_buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char),
			"error allocating send_buffer in gather_ecData.");
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	length = 0;
	for (int i = 1; i <= grid->num_ec_circumferentially; i++) {
		for (int j = 1; j <= grid->num_ec_axially; j++) {
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].B[cpl_IP3]);
		}
	}
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->ec_cpl[cpl_IP3] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->ec_cpl[cpl_IP3];
		writer_buffer->ec_cpl[cpl_IP3] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->cpIj = (char*) checked_malloc(writer_buffer->ec_cpl[cpl_IP3] * sizeof(char),
				"allocation error for writer_buffer member ec_cplIi.");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->cpIj, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);
}

void gather_JPLC_map(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, EC_cell** ec)
{
	printf("[%d] ===>>> Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	int num_tuple_components = 1;
	int num_tuples = grid->num_ec_axially * grid->num_ec_circumferentially;
	int root = 0;

	int *recv_count = (int*) checked_malloc(grid->tasks * sizeof(int), SRC_LOC);
	int *disp = (int*) checked_malloc(grid->tasks * sizeof(int), SRC_LOC);

	/************* JPLC field data *************/
	int my_branch_ec_offset = 0;

	char* send_buffer = (char*)checked_malloc(NUM_DBL_TO_CHAR_BYTES * num_tuple_components * num_tuples * sizeof(char), SRC_LOC);

	int length = 0;
	// Put things in the buffer, keep track of the length.
	for(int i = 1; i <= grid->num_ec_circumferentially; i++)
	{
		for(int j = 1; j <= grid->num_ec_axially; j++)
		{
			length += sprintf(send_buffer + length, "%2.12lf\n", ec[i][j].JPLC);
		}
	}

	// Collect the length of buffers on all ranks.
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm));

	writer_buffer->jplc_buffer_length = 0;

	for (int i = 0; i < grid->tasks; i++)
	{
		disp[i] = writer_buffer->jplc_buffer_length;
		writer_buffer->jplc_buffer_length += recv_count[i];
	}

	if (grid->rank == 0)
	{
		writer_buffer->jplc = (char*)checked_malloc(writer_buffer->jplc_buffer_length * sizeof(char), SRC_LOC);
	}
	CHECK_MPI_ERROR(MPI_Gatherv(send_buffer, length, MPI_CHAR, writer_buffer->jplc, recv_count, disp, MPI_CHAR, root, grid->cart_comm));

	free(recv_count);
	free(send_buffer);
	free(disp);
}

void checkpoint_coarse_time_profiling_data(grid_parms grid, time_stamps* t_stamp, IO_domain_info* my_IO_domain_info) {

	write_timing((char *)"aggregated_compute_time", grid, t_stamp->aggregate_compute, my_IO_domain_info);
	write_timing((char *)"aggregated_comm_time", grid, t_stamp->aggregate_comm, my_IO_domain_info);
	write_timing((char *)"aggregated_write_time", grid, t_stamp->aggregate_write, my_IO_domain_info);

	write_min_max_timing((char *)"min_max_of_aggregate_compute" , grid, t_stamp->aggregate_compute, my_IO_domain_info);
	write_min_max_timing((char *)"min_max_of_aggregate_comm" , grid, t_stamp->aggregate_comm, my_IO_domain_info);
	write_min_max_timing((char *)"min_max_of_aggregate_write" , grid, t_stamp->aggregate_write, my_IO_domain_info);
}

void write_timing(char* file_prefix, grid_parms grid, double field, IO_domain_info* my_IO_domain_info) {
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
	CHECK_MPI_ERROR(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid.cart_comm));
	int total_buffer_length = 0;
	for (int i = 0; i < grid.tasks; i++) {
		disp[i] = total_buffer_length;
		total_buffer_length += recv_count[i];
	}

	if (grid.rank == 0) {
		write_buffer = (char*) checked_malloc(total_buffer_length * sizeof(char),
				"allocation error for writer_buffer for time_profiling data in function push_coarse_timing_data_to_file().\n");
	}
	CHECK_MPI_ERROR(MPI_Gatherv(buffer, length, MPI_CHAR, write_buffer, recv_count, disp, MPI_CHAR, root, grid.cart_comm));

	if (grid.rank == 0) {
		sprintf(filename, "%s/%s_%s.txt", grid.time_profiling_dir, file_prefix, grid.suffix);
		CHECK_MPI_ERROR(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fw));
		CHECK_MPI_ERROR(MPI_File_write_at(fw, 0, write_buffer, total_buffer_length, MPI_CHAR, &status));
		MPI_File_close(&fw);
		free(write_buffer);
	}
	free(recv_count);
	free(buffer);
	free(disp);
}

void write_min_max_timing(char *file_prefix, grid_parms grid, double field, IO_domain_info* my_IO_domain_info)
{
	MPI_Status status;
	MPI_Offset displacement = 0;
	MPI_File fw;
	char* buffer = (char*) checked_malloc(NUM_DBL_TO_CHAR_BYTES * sizeof(char), SRC_LOC);
	char* write_buffer;
	int root = 0;
	char filename[50];
	int length = 0;
	double max, min;
	int max_ind, min_ind;
	double array[grid.tasks];

	int *recv_count = (int*) checked_malloc(grid.tasks * sizeof(int), SRC_LOC);
	int *disp = (int*) checked_malloc(grid.tasks * sizeof(int), SRC_LOC);

	CHECK_MPI_ERROR(MPI_Gather(&field, 1, MPI_DOUBLE, array, 1, MPI_DOUBLE, root, grid.cart_comm));

	if (grid.rank == 0)
	{
		maximum(array, grid.tasks, &max, &max_ind);
		minimum(array, grid.tasks, &min, &min_ind);
		write_buffer = (char*) checked_malloc(1024 * sizeof(char), SRC_LOC);
		sprintf(filename, "%s/%s_%s.txt", grid.time_profiling_dir, file_prefix, grid.suffix);

		CHECK_MPI_ERROR(MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fw));

		length = sprintf(write_buffer, "Maximum = %lf\t by Rank = %d\nMinimum = %lf\t by Rank = %d\n", max, max_ind, min, min_ind);

		CHECK_MPI_ERROR(MPI_File_write_at(fw, 0, write_buffer, length, MPI_CHAR, &status));

		MPI_File_close(&fw);
		free(write_buffer);
	}

	free(recv_count);
	free(buffer);
	free(disp);
}

