#include "computelib.h"

/************************************************/
void gather_tasks_mesh_point_data_on_writers_ver2(grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer, celltype1** smc,
		celltype2** ec) {

	int branch;
	if (grid->my_domain.internal_info.domain_type == STRSEG) {
		branch = P - 1;
	} else if (grid->my_domain.internal_info.domain_type == BIF) {
		branch = grid->branch_tag - 1;
	}

	/************* Point data *************/
	int num_tuple_components = 3, num_tuples = 4;
	int root = 0;
	int *recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	int *disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	double *send_buffer_double = (double*) checked_malloc(num_tuple_components * num_tuples * sizeof(double),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");
	int length = 0;
	///Number of vtk_points per task of each subdomain. This information will be sent to each writer.
	check_flag(MPI_Reduce(&num_tuples, &my_IO_domain_info->num_vtk_points_mesh, 1, MPI_INT, MPI_SUM, root, grid->cart_comm),
			"Error in MPI_Reduce in function gather_task_mesh_point_data.");

	for (int i = 0; i < num_tuples; i++) {
		for (int j = 0; j < num_tuple_components; j++) {
			send_buffer_double[i * num_tuple_components + j] = grid->coordinates[i][j];
			length++;
		}
	}

	/// Gathering and summing the length of all the CHARs contained in every send_buffer containing coordinates from each MPI process.
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[ProcessMesh] = 0;
	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ProcessMesh];
		writer_buffer->buffer_length[ProcessMesh] += recv_count[i];
	}

	if (grid->rank == 0) {
		writer_buffer->process_mesh_points_double = (double*) checked_malloc(writer_buffer->buffer_length[ProcessMesh] * sizeof(double),
				"allocation error for writer_buffer member process mesh.");
	}
	check_flag(
			MPI_Gatherv(send_buffer_double, length, MPI_DOUBLE, writer_buffer->process_mesh_points_double, recv_count, disp, MPI_DOUBLE, root,
					grid->cart_comm), "Error gathering data for task coordinates.");
	free(recv_count);
	free(send_buffer_double);
	free(disp);
	/************* Cell data *************/

	num_tuple_components = 5;
	num_tuples = 1;

	check_flag(MPI_Reduce(&num_tuples, &my_IO_domain_info->num_vtk_cells_mesh, 1, MPI_INT, MPI_SUM, root, grid->cart_comm),
			"Error in MPI_Reduce in function gather_task_mesh_point_data while gathering num_vtk_cells_mesh.\n");

	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	int *send_buffer_int = (int*) checked_malloc(num_tuple_components * num_tuples * sizeof(int),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");
	length = 0;
	for (int i = 0; i < num_tuples; i++) {
		send_buffer_int[i * num_tuple_components + 0] = 4;
		length++;
		for (int j = 1; j < num_tuple_components; j++) {
			send_buffer_int[i * num_tuple_components + j] = i * (num_tuple_components - 1) + j - 1;
			length++;
		}
	}

	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[ProcessCell] = 0;
	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ProcessCell];
		writer_buffer->buffer_length[ProcessCell] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->process_mesh_cells_int = (int*) checked_malloc(writer_buffer->buffer_length[ProcessCell] * sizeof(int),
				"allocation error for writer_buffer member process cells.");
	}
	check_flag(MPI_Gatherv(send_buffer_int, length, MPI_INT, writer_buffer->process_mesh_cells_int, recv_count, disp, MPI_INT, root, grid->cart_comm),
			"Error gathering data for task cells.");
	free(recv_count);
	free(send_buffer_int);
	free(disp);

	/************* Cell type data *************/
	num_tuple_components = 1;
	num_tuples = 1;
	recv_count = (int*) checked_malloc(grid->tasks * sizeof(int),
			"allocation failed for recv_count array in gather_tasks_mesh_point_data_on_writers");
	disp = (int*) checked_malloc(grid->tasks * sizeof(int), "allocation failed for disp array in gather_tasks_mesh_point_data_on_writers");

	send_buffer_int = (int*) checked_malloc(num_tuple_components * num_tuples * sizeof(int),
			"error allocating send_buffer in gather_tasks_mesh_point_data_on_writers.");

	length = 0;
	for (int i = 0; i < num_tuples; i++) {
		for (int j = 0; j < num_tuple_components; j++) {
			send_buffer_int[i * num_tuple_components + j] = 9;
			length++;
		}
	}
	check_flag(MPI_Gather(&length, 1, MPI_INT, recv_count, 1, MPI_INT, root, grid->cart_comm), "error in MPI_Gather.");
	writer_buffer->buffer_length[ProcessCellType] = 0;

	for (int i = 0; i < grid->tasks; i++) {
		disp[i] = writer_buffer->buffer_length[ProcessCellType];
		writer_buffer->buffer_length[ProcessCellType] += recv_count[i];
	}
	if (grid->rank == 0) {
		writer_buffer->process_mesh_type_int = (int*) checked_malloc(writer_buffer->buffer_length[ProcessCellType] * sizeof(int),
				"allocation error for writer_buffer member process cells.");
	}
	check_flag(MPI_Gatherv(send_buffer_int, length, MPI_INT, writer_buffer->process_mesh_type_int, recv_count, disp, MPI_INT, root, grid->cart_comm),
			"Error gathering data for task cells.");
	free(recv_count);
	free(send_buffer_int);
	free(disp);
}

/*****************************************************************************************/
void dump_process_data_ver2(checkpoint_handle* check, grid_parms* grid, IO_domain_info* my_IO_domain_info, data_buffer* writer_buffer,
		static_info_of_geometry* static_info, char* path)
		/*********************************************************************************************************************************/
		{
	MPI_Status status;
	int disp;
	char filename[50];
	int err;
	gather_static_info_on_writers(grid, my_IO_domain_info, static_info);
	err = sprintf(filename, "%s/task_mesh.vtk", path);
	check_flag(MPI_File_open(my_IO_domain_info->writer_comm, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &check->task_mesh),
			"Error opening file for writing task map.\n");

	int write_element_count = 1, header_offset[3] = { 0, 0, 0 }, point_offset = 0, cell_offset = 0, celltype_offset = 0, count = 0;
	char* header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at MPI_COMM_WORLD Rank 0.\n");

	int root = 0;
/// Get the sum of the number of points to be written by each participating writer. Total_points is significant/used
/// only by Rank 0 of the writer communicator.
	int Total_points = sum_array(my_IO_domain_info->writer_tasks, static_info->num_mesh_points), Total_cells = sum_array(
			my_IO_domain_info->writer_tasks, static_info->num_mesh_cells);

	/*************** Writing VTK header **************/
	sprintf(header, "# vtk DataFile Version 2.0\n"
			"Task mesh show how MPI processes are connected\n"
			"BINARY\n"
			"DATASET UNSTRUCTURED_GRID\n"
			"POINTS %d double\n", Total_points);

	header_offset[0] = strlen(header);
	count = header_offset[0];
	disp = 0;
	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->task_mesh, disp, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.");
	}

	disp = header_offset[0] * sizeof(char);
	check_flag(MPI_File_set_view(check->task_mesh, disp, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL), "Error adjusting file pointers.\n");
	/*************** Writing Point data **************/
	int buffer_lengths[my_IO_domain_info->writer_tasks]; 	/// Array holding the total number of points that each writer is going to write.
	check_flag(MPI_Allgather(&writer_buffer->buffer_length[ProcessMesh], 1, MPI_INT, buffer_lengths, 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths");
	disp = 0;
	for (int j = 0; j < my_IO_domain_info->writer_rank; j++) {
		disp += (static_info->num_mesh_points[j] * 3);
	}
	disp = disp * sizeof(double);
	count = static_info->num_mesh_points[my_IO_domain_info->writer_rank] * 3;
	check_flag(MPI_File_write_at(check->task_mesh, disp, writer_buffer->process_mesh_points_double, count, MPI_DOUBLE, &status),
			"error writing the coordinates in task_mesh.\n");

	for (int i = 0; i < my_IO_domain_info->writer_tasks; i++) {
		point_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->process_mesh_points_double);
	/*************** Writing VTK Cell connectivity**************/

	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at writers.\n");
	header_offset[1] = sprintf(header, "CELLS %d %d\n", Total_cells, 5 * Total_cells);

	count = header_offset[1];

	disp = header_offset[0] * sizeof(char) + point_offset * sizeof(double);
	check_flag(MPI_File_set_view(check->task_mesh, disp, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL), "Error adjusting file pointers.\n");

	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->task_mesh, 0, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.\n");
	}

	check_flag(MPI_Allgather(&writer_buffer->buffer_length[ProcessCell], 1, MPI_INT, buffer_lengths, 1, MPI_INT, my_IO_domain_info->writer_comm),
			"error in all gather called for buffer lengths");

	disp = header_offset[0] * sizeof(char) + point_offset * sizeof(double) + header_offset[1] * sizeof(char);
	check_flag(MPI_File_set_view(check->task_mesh, disp, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL), "Error adjusting file pointers.\n");

	disp = 0;
	for (int j = 0; j < my_IO_domain_info->writer_rank; j++) {
		disp += (static_info->num_mesh_cells[j] * 5);
	}

	for (int i = 0; i < static_info->num_mesh_cells[my_IO_domain_info->writer_rank] * 5; i += 4) {
		for (int j = 1; j < 5; j++) {
			if (my_IO_domain_info->writer_rank != 0) {
				writer_buffer->process_mesh_cells_int[i + j] += (disp * 4);
			}
		}
	}
	disp = disp * sizeof(int);
	count = static_info->num_mesh_cells[my_IO_domain_info->writer_rank] * 5;
	check_flag(MPI_File_write_at(check->task_mesh, disp, writer_buffer->process_mesh_cells_int, count, MPI_INT, &status),
			"error writing the cells in task_mesh.\n");
	for (int i = 0; i < my_IO_domain_info->writer_tasks; i++) {
		cell_offset += buffer_lengths[i];
	}


	if (my_IO_domain_info->writer_rank == 0) {
			double buf[Total_cells * 5];
			FILE *p, *wrt;
			p = fopen(filename, "rb");
			fseek(p, 0, SEEK_SET);
			wrt = fopen("tmp.txt", "w+");
			fseek(p, header_offset[0] * sizeof(char) + point_offset * sizeof(double) + header_offset[1] * sizeof(char), SEEK_SET);
			fread(buf, sizeof(int), Total_cells * 5, p);
			for (int i = 0; i < Total_cells * 5; i += 5)
				fprintf(wrt, "%d %d %d %d %d\n", buf[i], buf[i + 1], buf[i + 2], buf[i + 3], buf[i + 4]);
			fclose(wrt);
			fclose(p);
		}
		char f_name[20];
		sprintf(f_name, "tmp%d.txt", my_IO_domain_info->writer_rank);
		FILE* wrt = fopen(f_name, "w+");
		for (int i = 0; i < static_info->num_mesh_cells[my_IO_domain_info->writer_rank] * 5; i += 5)
			fprintf(wrt, "%d %d %d %d %d\n", writer_buffer->process_mesh_cells_int[i], writer_buffer->process_mesh_cells_int[i + 1]
					, writer_buffer->process_mesh_cells_int[i + 2], writer_buffer->process_mesh_cells_int[i + 3], writer_buffer->process_mesh_cells_int[i + 4]);
		fclose(wrt);

	free(header);
	free(writer_buffer->process_mesh_cells_int);
	/*************** Writing VTK Cell type **************/
	header = (char*) checked_malloc(1024 * sizeof(char), "allocation memory for writing header failed at writers.\n");
	header_offset[2] = sprintf(header, "CELL_TYPES %d\n", Total_cells);

	count = header_offset[2];
	disp = header_offset[0] * sizeof(char) + point_offset * sizeof(double) + header_offset[1] * sizeof(char) + cell_offset * sizeof(int);
	check_flag(MPI_File_set_view(check->task_mesh, disp, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL), "Error adjusting file pointers.\n");

	if (my_IO_domain_info->writer_rank == 0) {
		check_flag(MPI_File_write_at(check->task_mesh, 0, header, count, MPI_CHAR, &status), "error writing into time file by writer_rank 0.\n");
	}
	disp = header_offset[0] * sizeof(char) + point_offset * sizeof(double) + header_offset[1] * sizeof(char) + cell_offset * sizeof(int)
			+ header_offset[2] * sizeof(char);
	check_flag(MPI_File_set_view(check->task_mesh, disp, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL), "Error adjusting file pointers.\n");

	disp = 0;
	for (int j = 0; j < my_IO_domain_info->writer_rank; j++) {
		disp += static_info->num_mesh_cells[j];
	}
	disp = disp * sizeof(int);
	count = static_info->num_mesh_cells[my_IO_domain_info->writer_rank];

	check_flag(MPI_File_write_at(check->task_mesh, disp, writer_buffer->process_mesh_type_int, count, MPI_INT, &status),
			"error writing the cells in task_mesh.\n");
	for (int i = 0; i < 4; i++) {
		celltype_offset += buffer_lengths[i];
	}
	free(header);
	free(writer_buffer->process_mesh_type_int);

	MPI_File_close(&check->task_mesh);
}
/**************************************************************************************/

void gather_static_info_on_writers(grid_parms* grid, IO_domain_info* my_IO_domain_info, static_info_of_geometry* static_info) {
	int root = 0;
	int recv_count[my_IO_domain_info->writer_tasks];
	for (int i = 0; i < my_IO_domain_info->writer_tasks; i++) {
		recv_count[i] = 1;
	}
	int disp[my_IO_domain_info->writer_tasks];
	for (int i = 0; i < my_IO_domain_info->writer_tasks; i++) {
		disp[i] = i;
	}

	static_info->num_mesh_points = (int*) checked_malloc(my_IO_domain_info->writer_tasks * sizeof(int),
			"Allaocation for array num_mesh_points failed in struct static_info in dump_process_data.\n");
	static_info->num_mesh_cells = (int*) checked_malloc(my_IO_domain_info->writer_tasks * sizeof(int),
			"Allaocation for array num_mesh_cells failed in struct static_info in dump_process_data.\n");
	static_info->num_smc_points = (int*) checked_malloc(my_IO_domain_info->writer_tasks * sizeof(int),
			"Allaocation for array num_smc_points failed in struct static_info in dump_process_data.\n");
	static_info->num_smc_cells = (int*) checked_malloc(my_IO_domain_info->writer_tasks * sizeof(int),
			"Allaocation for array num_smc_cells failed in struct static_info in dump_process_data.\n");
	static_info->num_ec_points = (int*) checked_malloc(my_IO_domain_info->writer_tasks * sizeof(int),
			"Allaocation for array num_ec_points failed in struct static_info in dump_process_data.\n");
	static_info->num_ec_cells = (int*) checked_malloc(my_IO_domain_info->writer_tasks * sizeof(int),
			"Allaocation for array num_ec_cells failed in struct static_info in dump_process_data.\n");

	check_flag(
			MPI_Allgatherv(&my_IO_domain_info->num_vtk_points_mesh, 1, MPI_INT, static_info->num_mesh_points, recv_count, disp, MPI_INT,
					my_IO_domain_info->writer_comm), "Error in All_gather_v on all writers for num_mesh_points.\n");
	check_flag(
			MPI_Allgatherv(&my_IO_domain_info->num_vtk_cells_mesh, 1, MPI_INT, static_info->num_mesh_cells, recv_count, disp, MPI_INT,
					my_IO_domain_info->writer_comm), "Error in All_gather_v on all writers for num_mesh_cells.\n");
	for (int i = 0; i < my_IO_domain_info->writer_tasks; i++) {
		printf("[%d,%d,%d] : points = %d   ---   cells = %d\n", grid->universal_rank, my_IO_domain_info->writer_rank, i,
				static_info->num_mesh_points[i], static_info->num_mesh_cells[i]);
	}
}
