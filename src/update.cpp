#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"

/*extern	   conductance		cpl_cef;
 extern	   celltype1** 		smc;
 extern	   celltype2** 		ec;
 extern 	   double       	**sendbuf,**recvbuf;
 extern	   grid_parms		grid;
 */
using namespace std;
extern time_stamps t_stamp;

///************************************/
///************ check_flag*************/
///************************************/
void check_flag(int err, const char* errmsg) {
	if (err != MPI_SUCCESS) {
		fprintf(stdout, "%s", errmsg);
		MPI_Abort(MPI_COMM_WORLD, 200);
	}
}

void determine_source_destination(grid_parms grid, int source[], int dest[]) {
	if (grid.nbrs[local][UP] >= 0) {
		dest[UP] = grid.nbrs[local][UP];
		source[UP] = grid.nbrs[local][UP];
	} else if (grid.nbrs[local][UP] < 0) {
		dest[UP] = grid.rank; //MPI_PROC_NULL;
		source[UP] = grid.rank; //MPI_PROC_NULL;
	}

	if (grid.nbrs[local][DOWN] >= 0) {
		dest[DOWN] = grid.nbrs[local][DOWN];
		source[DOWN] = grid.nbrs[local][DOWN];
	} else if (grid.nbrs[local][DOWN] < 0) {
		dest[DOWN] = grid.rank; //MPI_PROC_NULL;
		source[DOWN] = grid.rank; //MPI_PROC_NULL;
	}

	if (grid.nbrs[local][LEFT] >= 0) {
		dest[LEFT] = grid.nbrs[local][LEFT];
		source[LEFT] = grid.nbrs[local][LEFT];
	} else if (grid.nbrs[local][LEFT] < 0) {
		dest[LEFT] = grid.rank; //MPI_PROC_NULL;
		source[LEFT] = grid.rank; //MPI_PROC_NULL;
	}

	if (grid.nbrs[local][RIGHT] >= 0) {
		dest[RIGHT] = grid.nbrs[local][RIGHT];
		source[RIGHT] = grid.nbrs[local][RIGHT];
	} else if (grid.nbrs[local][RIGHT] < 0) {
		dest[RIGHT] = grid.rank; //MPI_PROC_NULL;
		source[RIGHT] = grid.rank; // MPI_PROC_NULL;
	}
}
/*******************************************************************************************/
grid_parms communicate_num_recv_elements_to_nbrs(grid_parms grid)
/*******************************************************************************************/
{
	int err;
	MPI_Request reqs[8];
	MPI_Status stats[8];

	int source[4], dest[4];
	int tag;
	int num[4];
	determine_source_destination(grid, source, dest);

	tag = 0;
	int tmp[4];
	check_flag(
			MPI_Irecv(&grid.num_elements_recv_up, 1, MPI_INT, source[UP], tag,
					grid.cart_comm, &reqs[4 + UP]), "Irecv");
	check_flag(
			MPI_Isend(&grid.num_elements_send_up, 1, MPI_INT, dest[UP], tag,
					grid.cart_comm, &reqs[UP]), "Isend");
	check_flag(
			MPI_Irecv(&grid.num_elements_recv_down, 1, MPI_INT, source[DOWN],
					tag, grid.cart_comm, &reqs[4 + DOWN]), "Irecv");
	check_flag(
			MPI_Isend(&grid.num_elements_send_down, 1, MPI_INT, dest[DOWN], tag,
					grid.cart_comm, &reqs[DOWN]), "Isend");
	check_flag(
			MPI_Irecv(&grid.num_elements_recv_left, 1, MPI_INT, source[LEFT],
					tag, grid.cart_comm, &reqs[4 + LEFT]), "Irecv");
	check_flag(
			MPI_Isend(&grid.num_elements_send_left, 1, MPI_INT, dest[LEFT], tag,
					grid.cart_comm, &reqs[LEFT]), "Isend");
	check_flag(
			MPI_Irecv(&grid.num_elements_recv_right, 1, MPI_INT, source[RIGHT],
					tag, grid.cart_comm, &reqs[4 + RIGHT]), "Irecv");
	check_flag(
			MPI_Isend(&grid.num_elements_send_right, 1, MPI_INT, dest[RIGHT],
					tag, grid.cart_comm, &reqs[RIGHT]), "Isend");

	MPI_Waitall(8, reqs, stats);

	return (grid);

}

// Carry out asynchronous communication send and receive of edge cell data.
void communication_async_send_recv(grid_parms grid, double** sendbuf,
		double** recvbuf, SMC_cell** smc, EC_cell** ec)
		/*******************************************************************************************/
		{
	/// For recording an error in communication. err can either be MPI_SUCCESS or !MPI_SUCCESS
	int err;
	/// 8 MPI variables of types Request and Status are declared. 4 of each are for sending information to 4 neighbours and the other 4 are to retrive the receive operation status.
	MPI_Request reqs[16];
	MPI_Status stats[16];

	/// Two array, source and dest, hold the ranks for the communicating tasks.
	/// dest has ranks to which my task will send a message to.
	/// source contains tasks from which I expect a message.
	int source[4], dest[4];
	/// Message tag 1 and 2 representing segment 1 and 2 send or received by processor.
	int tag_1 = 1, tag_2 = 2;

	// Get nearest neighbours indices.
	determine_source_destination(grid, source, dest);
	t_stamp.update_sendbuf_t1 = MPI_Wtime();

	// Prepare the buffer for exchanging edge cell data with ghost cells.
	communication_update_sendbuf(grid, sendbuf, smc, ec);

	t_stamp.update_sendbuf_t2 = MPI_Wtime();
	t_stamp.diff_update_sendbuf = t_stamp.update_sendbuf_t2
			- t_stamp.update_sendbuf_t1;
	t_stamp.async_comm_calls_t1 = MPI_Wtime();


	/// Communication block
	check_flag(
			MPI_Irecv(&recvbuf[UP1][0], grid.num_elements_recv_up, MPI_DOUBLE,
					source[UP], tag_1, grid.cart_comm, &reqs[8 + UP1]),
			"Receive message operation for UP1");
	check_flag(
			MPI_Irecv(&recvbuf[UP2][0], grid.num_elements_recv_up, MPI_DOUBLE,
					source[UP], tag_2, grid.cart_comm, &reqs[8 + UP2]),
			"Receive message operation for UP2");

	check_flag(
			MPI_Irecv(&recvbuf[DOWN1][0], grid.num_elements_recv_down,
					MPI_DOUBLE, source[DOWN], tag_1, grid.cart_comm,
					&reqs[8 + DOWN1]), "Receive message operation for DOWN1");
	check_flag(
			MPI_Irecv(&recvbuf[DOWN2][0], grid.num_elements_recv_down,
					MPI_DOUBLE, source[DOWN], tag_2, grid.cart_comm,
					&reqs[8 + DOWN2]), "Receive message operation for DOWN2");

	check_flag(
			MPI_Irecv(&recvbuf[LEFT1][0], grid.num_elements_recv_left,
					MPI_DOUBLE, source[LEFT], tag_1, grid.cart_comm,
					&reqs[8 + LEFT1]), "Receive message operation for LEFT1");
	check_flag(
			MPI_Irecv(&recvbuf[LEFT2][0], grid.num_elements_recv_left,
					MPI_DOUBLE, source[LEFT], tag_2, grid.cart_comm,
					&reqs[8 + LEFT2]), "Receive message operation for LEFT2");

	check_flag(
			MPI_Irecv(&recvbuf[RIGHT1][0], grid.num_elements_recv_right,
					MPI_DOUBLE, source[RIGHT], tag_1, grid.cart_comm,
					&reqs[8 + RIGHT1]), "Receive message operation for RIGHT1");
	check_flag(
			MPI_Irecv(&recvbuf[RIGHT2][0], grid.num_elements_recv_right,
					MPI_DOUBLE, source[RIGHT], tag_2, grid.cart_comm,
					&reqs[8 + RIGHT2]), "Receive message operation for RIGHT2");
	check_flag(
			MPI_Isend(sendbuf[UP1], grid.num_elements_send_up, MPI_DOUBLE,
					dest[UP], tag_1, grid.cart_comm, &reqs[UP1]),
			"Send message operation for UP1");
	check_flag(
			MPI_Isend(sendbuf[UP2], grid.num_elements_send_up, MPI_DOUBLE,
					dest[UP], tag_2, grid.cart_comm, &reqs[UP2]),
			"Send message operation for UP2");

	check_flag(
			MPI_Isend(sendbuf[DOWN1], grid.num_elements_send_down, MPI_DOUBLE,
					dest[DOWN], tag_1, grid.cart_comm, &reqs[DOWN1]),
			"Send message operation for DOWN1");
	check_flag(
			MPI_Isend(sendbuf[DOWN2], grid.num_elements_send_down, MPI_DOUBLE,
					dest[DOWN], tag_2, grid.cart_comm, &reqs[DOWN2]),
			"Send message operation for DOWN2");

	check_flag(
			MPI_Isend(sendbuf[LEFT1], grid.num_elements_send_left, MPI_DOUBLE,
					dest[LEFT], tag_1, grid.cart_comm, &reqs[LEFT1]),
			"Send message operation for LEFT1");
	check_flag(
			MPI_Isend(sendbuf[LEFT2], grid.num_elements_send_left, MPI_DOUBLE,
					dest[LEFT], tag_2, grid.cart_comm, &reqs[LEFT2]),
			"Send message operation for LEFT2");
	check_flag(
			MPI_Isend(sendbuf[RIGHT1], grid.num_elements_send_right, MPI_DOUBLE,
					dest[RIGHT], tag_1, grid.cart_comm, &reqs[RIGHT1]),
			"Send message operation for RIGHT1");
	check_flag(
			MPI_Isend(sendbuf[RIGHT2], grid.num_elements_send_right, MPI_DOUBLE,
					dest[RIGHT], tag_2, grid.cart_comm, &reqs[RIGHT2]),
			"Send message operation for RIGHT2");
	t_stamp.async_comm_calls_t2 = MPI_Wtime();

	t_stamp.async_comm_calls_wait_t1 = MPI_Wtime();
	MPI_Waitall(16, reqs, stats);
	t_stamp.async_comm_calls_wait_t2 = MPI_Wtime();

	t_stamp.diff_async_comm_calls = t_stamp.async_comm_calls_t2
			- t_stamp.async_comm_calls_t1;
	t_stamp.diff_async_comm_calls_wait = t_stamp.async_comm_calls_wait_t2
			- t_stamp.async_comm_calls_wait_t1;

	int remote_source[4], remote_dest[4];
	MPI_Request remote_req[8];
	MPI_Status remote_status[8];
	int tag_remote_1 = 3, tag_remote_2 = 4;

	if (grid.my_domain.internal_info.boundary_tag == 'I') {
		t_stamp.remote_async_comm_calls_t1 = MPI_Wtime();
		for (int i = 0; i < 4; i++) {
			remote_source[i] = grid.nbrs[remote][i];
			remote_dest[i] = grid.nbrs[remote][i];
		}

		check_flag(
				MPI_Irecv(recvbuf[UP1], grid.num_elements_recv_up, MPI_DOUBLE,
						remote_source[UP1], tag_remote_1, grid.sub_universe,
						&remote_req[4 + UP1]),
				"Receive message operation for remote UP1");
		check_flag(
				MPI_Irecv(recvbuf[UP2], grid.num_elements_recv_up, MPI_DOUBLE,
						remote_source[UP2], tag_remote_2, grid.sub_universe,
						&remote_req[4 + UP2]),
				"Receive message operation for remote UP2");
		check_flag(
				MPI_Irecv(recvbuf[DOWN1], grid.num_elements_recv_down,
						MPI_DOUBLE, remote_source[DOWN1], tag_remote_1,
						grid.sub_universe, &remote_req[4 + DOWN1]),
				"Receive message operation for remote DOWN1");
		check_flag(
				MPI_Irecv(recvbuf[DOWN2], grid.num_elements_recv_down,
						MPI_DOUBLE, remote_source[DOWN2], tag_remote_2,
						grid.sub_universe, &remote_req[4 + DOWN2]),
				"Receive message operation for remote DOWN2");

		check_flag(
				MPI_Isend(sendbuf[UP1], grid.num_elements_send_up, MPI_DOUBLE,
						remote_dest[UP1], tag_remote_1, grid.sub_universe,
						&remote_req[UP1]),
				"Send message operation for remote UP1");
		check_flag(
				MPI_Isend(sendbuf[UP2], grid.num_elements_send_up, MPI_DOUBLE,
						remote_dest[UP2], tag_remote_2, grid.sub_universe,
						&remote_req[UP2]),
				"Send message operation for remote UP2");
		check_flag(
				MPI_Isend(sendbuf[DOWN1], grid.num_elements_send_down,
						MPI_DOUBLE, remote_dest[DOWN1], tag_remote_1,
						grid.sub_universe, &remote_req[DOWN1]),
				"Send message operation for remote DOWN1");
		check_flag(
				MPI_Isend(sendbuf[DOWN2], grid.num_elements_send_down,
						MPI_DOUBLE, remote_dest[DOWN2], tag_remote_2,
						grid.sub_universe, &remote_req[DOWN2]),
				"Send message operation for remote DOWN2");
		t_stamp.remote_async_comm_calls_wait_t1 = MPI_Wtime();
		MPI_Waitall(8, remote_req, remote_status);
		t_stamp.remote_async_comm_calls_wait_t2 = MPI_Wtime();
		t_stamp.diff_remote_async_comm_calls_wait =
				t_stamp.remote_async_comm_calls_wait_t2
						- t_stamp.remote_async_comm_calls_wait_t1;

		t_stamp.remote_async_comm_calls_t2 = MPI_Wtime();
		t_stamp.diff_remote_async_comm_calls =
				t_stamp.remote_async_comm_calls_t2
						- t_stamp.remote_async_comm_calls_t1;

	} else if ((grid.my_domain.internal_info.boundary_tag == 'T')
			|| (grid.my_domain.internal_info.boundary_tag == 'B')) {
		t_stamp.remote_async_comm_calls_t1 = MPI_Wtime();
		for (int i = 0; i < 4; i++) {
			remote_source[i] = grid.nbrs[remote][i];
			remote_dest[i] = grid.nbrs[remote][i];
		}

		check_flag(
				MPI_Irecv(recvbuf[UP1], grid.num_elements_recv_up, MPI_DOUBLE,
						remote_source[UP1], tag_remote_1, grid.universe,
						&remote_req[4 + UP1]),
				"Receive message operation for remote UP1");
		check_flag(
				MPI_Irecv(recvbuf[UP2], grid.num_elements_recv_up, MPI_DOUBLE,
						remote_source[UP2], tag_remote_2, grid.universe,
						&remote_req[4 + UP2]),
				"Receive message operation for remote UP2");
		check_flag(
				MPI_Irecv(recvbuf[DOWN1], grid.num_elements_recv_down,
						MPI_DOUBLE, remote_source[DOWN1], tag_remote_1,
						grid.universe, &remote_req[4 + DOWN1]),
				"Receive message operation for remote DOWN1");
		check_flag(
				MPI_Irecv(recvbuf[DOWN2], grid.num_elements_recv_down,
						MPI_DOUBLE, remote_source[DOWN2], tag_remote_2,
						grid.universe, &remote_req[4 + DOWN2]),
				"Receive message operation for remote DOWN2");

		check_flag(
				MPI_Isend(sendbuf[UP1], grid.num_elements_send_up, MPI_DOUBLE,
						remote_dest[UP1], tag_remote_1, grid.universe,
						&remote_req[UP1]),
				"Send message operation for remote UP1");
		check_flag(
				MPI_Isend(sendbuf[UP2], grid.num_elements_send_up, MPI_DOUBLE,
						remote_dest[UP2], tag_remote_2, grid.universe,
						&remote_req[UP2]),
				"Send message operation for remote UP2");
		check_flag(
				MPI_Isend(sendbuf[DOWN1], grid.num_elements_send_down,
						MPI_DOUBLE, remote_dest[DOWN1], tag_remote_1,
						grid.universe, &remote_req[DOWN1]),
				"Send message operation for remote DOWN1");
		check_flag(
				MPI_Isend(sendbuf[DOWN2], grid.num_elements_send_down,
						MPI_DOUBLE, remote_dest[DOWN2], tag_remote_2,
						grid.universe, &remote_req[DOWN2]),
				"Send message operation for remote DOWN2");
		t_stamp.remote_async_comm_calls_wait_t1 = MPI_Wtime();
		MPI_Waitall(8, remote_req, remote_status);
		t_stamp.remote_async_comm_calls_wait_t2 = MPI_Wtime();
		t_stamp.diff_remote_async_comm_calls_wait =
				t_stamp.remote_async_comm_calls_wait_t2
						- t_stamp.remote_async_comm_calls_wait_t1;

		t_stamp.remote_async_comm_calls_t2 = MPI_Wtime();
			t_stamp.diff_remote_async_comm_calls = t_stamp.remote_async_comm_calls_t2
					- t_stamp.remote_async_comm_calls_t1;
	}

	t_stamp.update_recvbuf_t1 = MPI_Wtime();

	// Unpack received data int ghost cells.
	communication_update_recvbuf(grid, recvbuf, smc, ec);

	t_stamp.update_recvbuf_t2 = MPI_Wtime();
	t_stamp.diff_update_recvbuf = t_stamp.update_recvbuf_t2
			- t_stamp.update_recvbuf_t1;

}	//end of update_async()

// Prepare the sendbuf pub putting ECs's and SMCs's edge data to be sent to the neighbours's ghost cells.
void communication_update_sendbuf(grid_parms grid, double** sendbuf,
		SMC_cell** smc, EC_cell** ec)
/*******************************************************************************************/
		{
	int k, buf_offset;

///UP direction	///
	///UP1
	buf_offset = grid.added_info_in_send_buf;
	//Updating SEND BUFFERwith SMC information  to be sent in UP1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int i = (int) sendbuf[UP1][0]; i <= (int) sendbuf[UP1][1]; i++) {
		int j = 1;
		sendbuf[UP1][buf_offset + k + 0] = smc[i][j].p[smc_Ca];
		sendbuf[UP1][buf_offset + k + 1] = smc[i][j].p[smc_Vm];
		sendbuf[UP1][buf_offset + k + 2] = smc[i][j].p[smc_IP3];
		k += grid.num_coupling_species_smc;
	}

	//Setting up the offset to transfer the EC info into SEND BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (sendbuf[UP1][1] - sendbuf[UP1][0] + 1);
	//Updating SEND BUFFER with EC information to be sent in UP1 direction with corresponding downstream neighbour cells located locally.
	k = 0;

	for (int i = (int) sendbuf[UP1][2]; i <= (int) sendbuf[UP1][3]; i++) {
		int j = 1;
		sendbuf[UP1][buf_offset + k + 0] = ec[i][j].q[ec_Ca];
		sendbuf[UP1][buf_offset + k + 1] = ec[i][j].q[ec_Vm];
		sendbuf[UP1][buf_offset + k + 2] = ec[i][j].q[ec_IP3];
		k += grid.num_coupling_species_ec;

	}

	///UP2
	buf_offset = grid.added_info_in_send_buf;
	//Updating SEND BUFFER to be sent in UP2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int i = (int) sendbuf[UP2][0]; i <= (int) sendbuf[UP2][1]; i++) {
		int j = 1;
		sendbuf[UP2][buf_offset + k + 0] = smc[i][j].p[smc_Ca];
		sendbuf[UP2][buf_offset + k + 1] = smc[i][j].p[smc_Vm];
		sendbuf[UP2][buf_offset + k + 2] = smc[i][j].p[smc_IP3];
		k += grid.num_coupling_species_smc;
	}

	//Setting up the offset to transfer the EC info into SEND BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (sendbuf[UP2][1] - sendbuf[UP2][0] + 1);
	//Updating SEND BUFFER to be sent in UP2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int i = (int) sendbuf[UP2][2]; i <= (int) sendbuf[UP2][3]; i++) {
		int j = 1;
		sendbuf[UP2][buf_offset + k + 0] = ec[i][j].q[ec_Ca];
		sendbuf[UP2][buf_offset + k + 1] = ec[i][j].q[ec_Vm];
		sendbuf[UP2][buf_offset + k + 2] = ec[i][j].q[ec_IP3];
		k += grid.num_coupling_species_ec;
	}

///DOWN direction	///
	///DOWN1
	buf_offset = grid.added_info_in_send_buf;
	//Updating SEND BUFFERwith SMC information  to be sent in DOWN1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int i = (int) sendbuf[DOWN1][0]; i <= (int) sendbuf[DOWN1][1]; i++) {
		int j = grid.num_smc_axially;
		sendbuf[DOWN1][buf_offset + k + 0] = smc[i][j].p[smc_Ca];
		sendbuf[DOWN1][buf_offset + k + 1] = smc[i][j].p[smc_Vm];
		sendbuf[DOWN1][buf_offset + k + 2] = smc[i][j].p[smc_IP3];
		k += grid.num_coupling_species_smc;
	}
	//Setting up the offset to transfer the EC info into SEND BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (sendbuf[DOWN1][1] - sendbuf[DOWN1][0] + 1);
	//Updating SEND BUFFER with EC information to be sent in DOWN1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int i = (int) sendbuf[DOWN1][2]; i <= (int) sendbuf[DOWN1][3]; i++) {
		int j = grid.num_ec_axially;
		sendbuf[DOWN1][buf_offset + k + 0] = ec[i][j].q[ec_Ca];
		sendbuf[DOWN1][buf_offset + k + 1] = ec[i][j].q[ec_Vm];
		sendbuf[DOWN1][buf_offset + k + 2] = ec[i][j].q[ec_IP3];
		k += grid.num_coupling_species_ec;
	}

	///DOWN2
	buf_offset = grid.added_info_in_send_buf;
	//Updating SEND BUFFER to be sent in DOWN2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int i = (int) sendbuf[DOWN2][0]; i <= (int) sendbuf[DOWN2][1]; i++) {
		int j = grid.num_smc_axially;
		sendbuf[DOWN2][buf_offset + k + 0] = smc[i][j].p[smc_Ca];
		sendbuf[DOWN2][buf_offset + k + 1] = smc[i][j].p[smc_Vm];
		sendbuf[DOWN2][buf_offset + k + 2] = smc[i][j].p[smc_IP3];
		k += grid.num_coupling_species_smc;
	}

	//Setting up the offset to transfer the EC info into SEND BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (sendbuf[DOWN2][1] - sendbuf[DOWN2][0] + 1);
	//Updating SEND BUFFER to be sent in DOWN2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int i = (int) sendbuf[DOWN2][2]; i <= (int) sendbuf[DOWN2][3]; i++) {
		int j = grid.num_ec_axially;
		sendbuf[DOWN2][buf_offset + k + 0] = ec[i][j].q[ec_Ca];
		sendbuf[DOWN2][buf_offset + k + 1] = ec[i][j].q[ec_Vm];
		sendbuf[DOWN2][buf_offset + k + 2] = ec[i][j].q[ec_IP3];
		k += grid.num_coupling_species_ec;
	}

///LEFT direction	///
	///LEFT1
	buf_offset = grid.added_info_in_send_buf;
	//Updating SEND BUFFER with SMC information  to be sent in LEFT1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) sendbuf[LEFT1][0]; j <= (int) sendbuf[LEFT1][1]; j++) {
		int i = 1;
		sendbuf[LEFT1][buf_offset + k + 0] = smc[i][j].p[smc_Ca];
		sendbuf[LEFT1][buf_offset + k + 1] = smc[i][j].p[smc_Vm];
		sendbuf[LEFT1][buf_offset + k + 2] = smc[i][j].p[smc_IP3];
		k += grid.num_coupling_species_smc;
	}
	//Setting up the offset to transfer the EC info into SEND BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (sendbuf[LEFT1][1] - sendbuf[LEFT1][0] + 1);
	//Updating SEND BUFFER with EC information to be sent in DOWN1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) sendbuf[LEFT1][2]; j <= (int) sendbuf[LEFT1][3]; j++) {
		int i = 1;
		sendbuf[LEFT1][buf_offset + k + 0] = ec[i][j].q[ec_Ca];
		sendbuf[LEFT1][buf_offset + k + 1] = ec[i][j].q[ec_Vm];
		sendbuf[LEFT1][buf_offset + k + 2] = ec[i][j].q[ec_IP3];
		k += grid.num_coupling_species_ec;
	}

	///LEFT2
	buf_offset = grid.added_info_in_send_buf;
	//Updating SEND BUFFER to be sent in LEFT2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) sendbuf[LEFT2][0]; j <= (int) sendbuf[LEFT2][1]; j++) {
		int i = 1;
		sendbuf[LEFT2][buf_offset + k + 0] = smc[i][j].p[smc_Ca];
		sendbuf[LEFT2][buf_offset + k + 1] = smc[i][j].p[smc_Vm];
		sendbuf[LEFT2][buf_offset + k + 2] = smc[i][j].p[smc_IP3];
		k += grid.num_coupling_species_smc;
	}
	//Setting up the offset to transfer the EC info into SEND BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (sendbuf[LEFT2][1] - sendbuf[LEFT2][0] + 1);
	//Updating SEND BUFFER to be sent in DOWN2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) sendbuf[LEFT2][2]; j <= (int) sendbuf[LEFT2][3]; j++) {
		int i = 1;
		sendbuf[LEFT2][buf_offset + k + 0] = ec[i][j].q[ec_Ca];
		sendbuf[LEFT2][buf_offset + k + 1] = ec[i][j].q[ec_Vm];
		sendbuf[LEFT2][buf_offset + k + 2] = ec[i][j].q[ec_IP3];
		k += grid.num_coupling_species_ec;
	}

///RIGHT direction	///
	///RIGHT1
	buf_offset = grid.added_info_in_send_buf;
	//Updating SEND BUFFER with SMC information  to be sent in LEFT1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) sendbuf[RIGHT1][0]; j <= (int) sendbuf[RIGHT1][1]; j++) {
		int i = grid.num_smc_circumferentially;
		sendbuf[RIGHT1][buf_offset + k + 0] = smc[i][j].p[smc_Ca];
		sendbuf[RIGHT1][buf_offset + k + 1] = smc[i][j].p[smc_Vm];
		sendbuf[RIGHT1][buf_offset + k + 2] = smc[i][j].p[smc_IP3];
		k += grid.num_coupling_species_smc;
	}
	//Setting up the offset to transfer the EC info into SEND BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (sendbuf[RIGHT1][1] - sendbuf[RIGHT1][0] + 1);
	//Updating SEND BUFFER with EC information to be sent in RIGHT1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) sendbuf[RIGHT1][2]; j <= (int) sendbuf[RIGHT1][3]; j++) {
		int i = grid.num_ec_circumferentially;
		sendbuf[RIGHT1][buf_offset + k + 0] = ec[i][j].q[ec_Ca];
		sendbuf[RIGHT1][buf_offset + k + 1] = ec[i][j].q[ec_Vm];
		sendbuf[RIGHT1][buf_offset + k + 2] = ec[i][j].q[ec_IP3];
		k += grid.num_coupling_species_ec;
	}

	///RIGHT2
	buf_offset = grid.added_info_in_send_buf;
	//Updating SEND BUFFER to be sent in RIGHT2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) sendbuf[RIGHT2][0]; j <= (int) sendbuf[RIGHT2][1]; j++) {
		int i = grid.num_smc_circumferentially;
		sendbuf[RIGHT2][buf_offset + k + 0] = smc[i][j].p[smc_Ca];
		sendbuf[RIGHT2][buf_offset + k + 1] = smc[i][j].p[smc_Vm];
		sendbuf[RIGHT2][buf_offset + k + 2] = smc[i][j].p[smc_IP3];
		k += grid.num_coupling_species_smc;
	}

	//Setting up the offset to transfer the EC info into SEND BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (sendbuf[RIGHT2][1] - sendbuf[RIGHT2][0] + 1);
	//Updating SEND BUFFER to be sent in RIGHT2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) sendbuf[RIGHT2][2]; j <= (int) sendbuf[RIGHT2][3]; j++) {
		int i = grid.num_ec_circumferentially;

		sendbuf[RIGHT2][buf_offset + k + 0] = ec[i][j].q[ec_Ca];
		sendbuf[RIGHT2][buf_offset + k + 1] = ec[i][j].q[ec_Vm];
		sendbuf[RIGHT2][buf_offset + k + 2] = ec[i][j].q[ec_IP3];
		k += grid.num_coupling_species_ec;
	}

}	// end of communication_update_sendbuf()

// Unpacking received data into ghost cells.
void communication_update_recvbuf(grid_parms grid, double** recvbuf,
		SMC_cell** smc, EC_cell** ec)
		/*******************************************************************************************/
		{
	int k, buf_offset;
///UP direction	///
///UP1
	buf_offset = grid.added_info_in_send_buf;
//Updating RECEIVE BUFFER with SMC information  to be sent in UP1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int i = (int) recvbuf[UP1][0]; i <= (int) recvbuf[UP1][1]; i++) {
		int j = 0;
		smc[i][j].p[smc_Ca] = recvbuf[UP1][buf_offset + k + 0];
		smc[i][j].p[smc_Vm] = recvbuf[UP1][buf_offset + k + 1];
		smc[i][j].p[smc_IP3] = recvbuf[UP1][buf_offset + k + 2];
		k += grid.num_coupling_species_smc;
	}
//Setting up the offset to transfer the EC info into RECEIVE BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (recvbuf[UP1][1] - recvbuf[UP1][0] + 1);
//Updating RECEIVE BUFFER with EC information to be sent in UP1 direction with corresponding downstream neighbour cells located locally.
	k = 0;

	for (int i = (int) recvbuf[UP1][2]; i <= (int) recvbuf[UP1][3]; i++) {
		int j = 0;
		ec[i][j].q[ec_Ca] = recvbuf[UP1][buf_offset + k + 0];
		ec[i][j].q[ec_Vm] = recvbuf[UP1][buf_offset + k + 1];
		ec[i][j].q[ec_IP3] = recvbuf[UP1][buf_offset + k + 2];
		k += grid.num_coupling_species_ec;
	}

///UP2
	buf_offset = grid.added_info_in_send_buf;
//Updating RECEIVE BUFFER to be sent in UP2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int i = (int) recvbuf[UP2][0]; i <= (int) recvbuf[UP2][1]; i++) {
		int j = 0;
		smc[i][j].p[smc_Ca] = recvbuf[UP2][buf_offset + k + 0];
		smc[i][j].p[smc_Vm] = recvbuf[UP2][buf_offset + k + 1];
		smc[i][j].p[smc_IP3] = recvbuf[UP2][buf_offset + k + 2];
		k += grid.num_coupling_species_smc;
	}

//Setting up the offset to transfer the EC info into RECEIVE BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (recvbuf[UP2][1] - recvbuf[UP2][0] + 1);
//Updating RECEIVE BUFFER to be sent in UP2 direction with corresponding downstream neighbour cells located locally.
	k = 0;

	for (int i = (int) recvbuf[UP2][2]; i <= (int) recvbuf[UP2][3]; i++) {
		int j = 0;
		ec[i][j].q[ec_Ca] = recvbuf[UP2][buf_offset + k + 0];
		ec[i][j].q[ec_Vm] = recvbuf[UP2][buf_offset + k + 1];
		ec[i][j].q[ec_IP3] = recvbuf[UP2][buf_offset + k + 2];
		k += grid.num_coupling_species_ec;
	}

///DOWN direction	///
///DOWN1
	buf_offset = grid.added_info_in_send_buf;
//Updating RECEIVE BUFFERwith SMC information  to be sent in DOWN1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	if (grid.flip_array[DOWN1] == 0) {
		for (int i = (int) recvbuf[DOWN1][0]; i <= (int) recvbuf[DOWN1][1];
				i++) {
			int j = grid.num_smc_axially + 1;
			smc[i][j].p[smc_Ca] = recvbuf[DOWN1][buf_offset + k + 0];
			smc[i][j].p[smc_Vm] = recvbuf[DOWN1][buf_offset + k + 1];
			smc[i][j].p[smc_IP3] = recvbuf[DOWN1][buf_offset + k + 2];
			k += grid.num_coupling_species_smc;
		}
	} else if (grid.flip_array[DOWN1] == 1) {
		int start = (int) recvbuf[DOWN2][0], end = (int) recvbuf[DOWN2][1];
		for (int i = end; i >= start; i--) {
			int j = grid.num_smc_axially + 1;
			smc[i][j].p[smc_Ca] = recvbuf[DOWN1][buf_offset + k + 0];
			smc[i][j].p[smc_Vm] = recvbuf[DOWN1][buf_offset + k + 1];
			smc[i][j].p[smc_IP3] = recvbuf[DOWN1][buf_offset + k + 2];
			k += grid.num_coupling_species_smc;
		}
	}

//Setting up the offset to transfer the EC info into RECEIVE BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (recvbuf[DOWN1][1] - recvbuf[DOWN1][0] + 1);
//Updating RECEIVE BUFFER with EC information to be sent in DOWN1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	if (grid.flip_array[DOWN1] == 0) {
		for (int i = (int) recvbuf[DOWN1][2]; i <= (int) recvbuf[DOWN1][3];
				i++) {
			int j = grid.num_ec_axially + 1;
			ec[i][j].q[ec_Ca] = recvbuf[DOWN1][buf_offset + k + 0];
			ec[i][j].q[ec_Vm] = recvbuf[DOWN1][buf_offset + k + 1];
			ec[i][j].q[ec_IP3] = recvbuf[DOWN1][buf_offset + k + 2];
			k += grid.num_coupling_species_ec;
		}
	} else if (grid.flip_array[DOWN1] == 1) {
		int start = (int) recvbuf[DOWN2][2], end = (int) recvbuf[DOWN2][3];
		for (int i = end; i >= start; i--) {
			int j = grid.num_ec_axially + 1;
			ec[i][j].q[ec_Ca] = recvbuf[DOWN1][buf_offset + k + 0];
			ec[i][j].q[ec_Vm] = recvbuf[DOWN1][buf_offset + k + 1];
			ec[i][j].q[ec_IP3] = recvbuf[DOWN1][buf_offset + k + 2];
			k += grid.num_coupling_species_ec;
		}
	}

///DOWN2
	buf_offset = grid.added_info_in_send_buf;
//Updating RECEIVE BUFFER to be sent in DOWN2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	if (grid.flip_array[DOWN2] == 0) {
		for (int i = (int) recvbuf[DOWN2][0]; i <= (int) recvbuf[DOWN2][1];
				i++) {
			int j = grid.num_smc_axially + 1;
			smc[i][j].p[smc_Ca] = recvbuf[DOWN2][buf_offset + k + 0];
			smc[i][j].p[smc_Vm] = recvbuf[DOWN2][buf_offset + k + 1];
			smc[i][j].p[smc_IP3] = recvbuf[DOWN2][buf_offset + k + 2];
			k += grid.num_coupling_species_smc;
		}
	} else if (grid.flip_array[DOWN2] == 1) {
		int start = (int) recvbuf[DOWN1][0], end = (int) recvbuf[DOWN1][1];
		for (int i = end; i >= start; i--) {
			int j = grid.num_smc_axially + 1;
			smc[i][j].p[smc_Ca] = recvbuf[DOWN2][buf_offset + k + 0];
			smc[i][j].p[smc_Vm] = recvbuf[DOWN2][buf_offset + k + 1];
			smc[i][j].p[smc_IP3] = recvbuf[DOWN2][buf_offset + k + 2];
			k += grid.num_coupling_species_smc;
		}
	}

//Setting up the offset to transfer the EC info into RECEIVE BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (recvbuf[DOWN2][1] - recvbuf[DOWN2][0] + 1);
//Updating RECEIVE BUFFER to be sent in DOWN2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	if (grid.flip_array[DOWN2] == 0) {
		for (int i = (int) recvbuf[DOWN2][2]; i <= (int) recvbuf[DOWN2][3];
				i++) {
			int j = grid.num_ec_axially + 1;
			ec[i][j].q[ec_Ca] = recvbuf[DOWN2][buf_offset + k + 0];
			ec[i][j].q[ec_Vm] = recvbuf[DOWN2][buf_offset + k + 1];
			ec[i][j].q[ec_IP3] = recvbuf[DOWN2][buf_offset + k + 2];
			k += grid.num_coupling_species_ec;
		}
	} else if (grid.flip_array[DOWN2] == 1) {
		int start = (int) recvbuf[DOWN1][2], end = (int) recvbuf[DOWN1][3];
		for (int i = end; i >= start; i--) {
			int j = grid.num_ec_axially + 1;
			ec[i][j].q[ec_Ca] = recvbuf[DOWN2][buf_offset + k + 0];
			ec[i][j].q[ec_Vm] = recvbuf[DOWN2][buf_offset + k + 1];
			ec[i][j].q[ec_IP3] = recvbuf[DOWN2][buf_offset + k + 2];
			k += grid.num_coupling_species_ec;
		}
	}

///LEFT direction	///
///LEFT1
	buf_offset = grid.added_info_in_send_buf;
//Updating RECEIVE BUFFER with SMC information  to be sent in LEFT1 direction with corresponding downstream neighbour cells located locally.
	k = 0;

	for (int j = (int) recvbuf[LEFT1][0]; j <= (int) recvbuf[LEFT1][1]; j++) {
		int i = 0;
		smc[i][j].p[smc_Ca] = recvbuf[LEFT1][buf_offset + k + 0];
		smc[i][j].p[smc_Vm] = recvbuf[LEFT1][buf_offset + k + 1];
		smc[i][j].p[smc_IP3] = recvbuf[LEFT1][buf_offset + k + 2];
		k += grid.num_coupling_species_smc;
	}
//Setting up the offset to transfer the EC info into RECEIVE BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (recvbuf[LEFT1][1] - recvbuf[LEFT1][0] + 1);
//Updating RECEIVE BUFFER with EC information to be sent in DOWN1 direction with corresponding downstream neighbour cells located locally.
	k = 0;

	for (int j = (int) recvbuf[LEFT1][2]; j <= (int) recvbuf[LEFT1][3]; j++) {
		int i = 0;
		ec[i][j].q[ec_Ca] = recvbuf[LEFT1][buf_offset + k + 0];
		ec[i][j].q[ec_Vm] = recvbuf[LEFT1][buf_offset + k + 1];
		ec[i][j].q[ec_IP3] = recvbuf[LEFT1][buf_offset + k + 2];
		k += grid.num_coupling_species_ec;
	}

///LEFT2
	buf_offset = grid.added_info_in_send_buf;
//Updating RECEIVE BUFFER to be sent in LEFT2 direction with corresponding downstream neighbour cells located locally.
	k = 0;

	for (int j = (int) recvbuf[LEFT2][0]; j <= (int) recvbuf[LEFT2][1]; j++) {
		int i = 0;
		smc[i][j].p[smc_Ca] = recvbuf[LEFT2][buf_offset + k + 0];
		smc[i][j].p[smc_Vm] = recvbuf[LEFT2][buf_offset + k + 1];
		smc[i][j].p[smc_IP3] = recvbuf[LEFT2][buf_offset + k + 2];
		k += grid.num_coupling_species_smc;
	}

//Setting up the offset to transfer the EC info into RECEIVE BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (recvbuf[LEFT2][1] - recvbuf[LEFT2][0] + 1);
//Updating RECEIVE BUFFER to be sent in LEFT2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) recvbuf[LEFT2][2]; j <= (int) recvbuf[LEFT2][3]; j++) {
		int i = 0;
		ec[i][j].q[ec_Ca] = recvbuf[LEFT2][buf_offset + k + 0];
		ec[i][j].q[ec_Vm] = recvbuf[LEFT2][buf_offset + k + 1];
		ec[i][j].q[ec_IP3] = recvbuf[LEFT2][buf_offset + k + 2];
		k += grid.num_coupling_species_ec;
	}

///RIGHT direction	///
///RIGHT1
	buf_offset = grid.added_info_in_send_buf;
//Updating RECEIVE BUFFER with SMC information  to be sent in LEFT1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) recvbuf[RIGHT1][0]; j <= (int) recvbuf[RIGHT1][1]; j++) {
		int i = grid.num_smc_circumferentially + 1;
		smc[i][j].p[smc_Ca] = recvbuf[RIGHT1][buf_offset + k + 0];
		smc[i][j].p[smc_Vm] = recvbuf[RIGHT1][buf_offset + k + 1];
		smc[i][j].p[smc_IP3] = recvbuf[RIGHT1][buf_offset + k + 2];
		k += grid.num_coupling_species_smc;
	}
//Setting up the offset to transfer the EC info into RECEIVE BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (recvbuf[RIGHT1][1] - recvbuf[RIGHT1][0] + 1);
//Updating RECEIVE BUFFER with EC information to be sent in RIGHT1 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) recvbuf[RIGHT1][2]; j <= (int) recvbuf[RIGHT1][3]; j++) {
		int i = grid.num_ec_circumferentially + 1;
		ec[i][j].q[ec_Ca] = recvbuf[RIGHT1][buf_offset + k + 0];
		ec[i][j].q[ec_Vm] = recvbuf[RIGHT1][buf_offset + k + 1];
		ec[i][j].q[ec_IP3] = recvbuf[RIGHT1][buf_offset + k + 2];
		k += grid.num_coupling_species_ec;
	}

///RIGHT2
	buf_offset = grid.added_info_in_send_buf;
//Updating RECEIVE BUFFER to be sent in RIGHT2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) recvbuf[RIGHT2][0]; j <= (int) recvbuf[RIGHT2][1]; j++) {
		int i = grid.num_smc_circumferentially + 1;
		smc[i][j].p[smc_Ca] = recvbuf[RIGHT2][buf_offset + k + 0];
		smc[i][j].p[smc_Vm] = recvbuf[RIGHT2][buf_offset + k + 1];
		smc[i][j].p[smc_IP3] = recvbuf[RIGHT2][buf_offset + k + 2];
		k += grid.num_coupling_species_smc;
	}

//Setting up the offset to transfer the EC info into RECEIVE BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc
					* (recvbuf[RIGHT2][1] - recvbuf[RIGHT2][0] + 1);
//Updating RECEIVE BUFFER to be sent in RIGHT2 direction with corresponding downstream neighbour cells located locally.
	k = 0;
	for (int j = (int) recvbuf[RIGHT2][2]; j <= (int) recvbuf[RIGHT2][3]; j++) {
		int i = grid.num_ec_circumferentially + 1;
		ec[i][j].q[ec_Ca] = recvbuf[RIGHT2][buf_offset + k + 0];
		ec[i][j].q[ec_Vm] = recvbuf[RIGHT2][buf_offset + k + 1];
		ec[i][j].q[ec_IP3] = recvbuf[RIGHT2][buf_offset + k + 2];
		k += grid.num_coupling_species_ec;
	}
}	// end of communication_update_recvbuf()


