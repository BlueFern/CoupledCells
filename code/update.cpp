#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "computelib.h"

extern	   conductance		cpl_cef;
extern	   celltype1** 		smc;
extern	   celltype2** 		ec;
extern 	   double       	**sendbuf,**recvbuf;
extern	   grid_parms		grid;

using namespace std;

///************************************/
///************ check_flag*************/
///************************************/
void check_flag(int err, FILE* errfile, const char* errmsg){
	if (err !=MPI_SUCCESS) {
		fprintf(errfile,"MPI-Comm Error:%s\n",errmsg);
		MPI_Abort(MPI_COMM_WORLD, 200);
	}
}

void determin_source_destination(int source[], int dest[]){
	if (grid.nbrs[local][UP] >= 0) {
				dest[UP] = grid.nbrs[local][UP];
				source[UP] = grid.nbrs[local][UP];
			} else if (grid.nbrslocal][UP] < 0) {
				dest[UP] = grid.rank;//MPI_PROC_NULL;
				source[UP] = grid.rank;//MPI_PROC_NULL;
			}

			if (grid.nbrs[DOWN] >= 0) {
				dest[DOWN] = grid.nbrs[DOWN];
				source[DOWN] = grid.nbrs[DOWN];
			} else if (grid.nbrs[DOWN] < 0) {
				dest[DOWN] = grid.rank;//MPI_PROC_NULL;
				source[DOWN] = grid.rank;//MPI_PROC_NULL;
			}

			if (grid.nbrs[LEFT] >= 0) {
				dest[LEFT] = grid.nbrs[LEFT];
				source[LEFT] = grid.nbrs[LEFT];
			} else if (grid.nbrs[LEFT] < 0) {
				dest[LEFT] = grid.rank;//MPI_PROC_NULL;
				source[LEFT] = grid.rank;//MPI_PROC_NULL;
			}

			if (grid.nbrs[RIGHT] >= 0) {
				dest[RIGHT] = grid.nbrs[RIGHT];
				source[RIGHT] = grid.nbrs[RIGHT];
			} else if (grid.nbrs[RIGHT] < 0) {
				dest[RIGHT] = grid.rank;//MPI_PROC_NULL;
				source[RIGHT] =grid.rank;// MPI_PROC_NULL;
			}
}
/*******************************************************************************************/
void communicate_num_recv_elements_to_nbrs(FILE* logptr)
/*******************************************************************************************/
{
	int err;
	MPI_Request reqs[8], req[4];
	MPI_Status stats[8];

	int source[4], dest[4];
	int tag;
	int num[4];
	determin_source_destination(source,dest);

	tag = 0;
	check_flag(MPI_Irecv(&grid.num_elements_recv_up, 1, MPI_INT, source[UP], tag,
			grid.comm, &reqs[4 + UP]), logptr,"Irecv" );
	check_flag(MPI_Isend(&grid.num_elements_send_up, 1, MPI_INT, dest[UP], tag,
			grid.comm, &reqs[UP]), logptr,"Isend" );
	check_flag(MPI_Irecv(&grid.num_elements_recv_down, 1, MPI_INT, source[DOWN], tag,
			grid.comm, &reqs[4 + DOWN]), logptr,"Irecv" );
	check_flag(MPI_Isend(&grid.num_elements_send_down, 1, MPI_INT, dest[DOWN], tag,
			grid.comm, &reqs[DOWN]), logptr,"Isend" );
	check_flag(MPI_Irecv(&grid.num_elements_recv_left, 1, MPI_INT, source[LEFT], tag,
			grid.comm, &reqs[4 + LEFT]), logptr,"Irecv" );
	check_flag(MPI_Isend(&grid.num_elements_send_left, 1, MPI_INT, dest[LEFT], tag,
			grid.comm, &reqs[LEFT]), logptr,"Isend" );
	check_flag(MPI_Irecv(&grid.num_elements_recv_right, 1, MPI_INT, source[RIGHT], tag,
			grid.comm, &reqs[4 + RIGHT]), logptr,"Irecv" );
	check_flag(MPI_Isend(&grid.num_elements_send_right, 1, MPI_INT, dest[RIGHT], tag,
			grid.comm, &reqs[RIGHT]), logptr,"Isend" );

	MPI_Waitall(8, reqs, stats);

}
/*******************************************************************************************/
void communication_async_send_recv(FILE* logptr)
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
	int source[4],dest[4];
	/// Message tag 1 and 2 representing segment 1 and 2 send or received by processor.
	int tag_1=1,tag_2=2;
	determin_source_destination(source,dest);

communication_update_sendbuf(logptr);
	/// Communication block
check_flag(
		MPI_Irecv(&recvbuf[UP1][0], grid.num_elements_recv_up, MPI_DOUBLE,
				source[UP], tag_1, grid.comm, &reqs[8 + UP1]), logptr,
		"Receive message operation for UP1");
check_flag(
		MPI_Irecv(&recvbuf[UP2][0], grid.num_elements_recv_up, MPI_DOUBLE,
				source[UP], tag_2, grid.comm, &reqs[8 + UP2]), logptr,
		"Receive message operation for UP2");

check_flag(
		MPI_Irecv(&recvbuf[DOWN1][0], grid.num_elements_recv_down, MPI_DOUBLE,
				source[DOWN], tag_1, grid.comm, &reqs[8 + DOWN1]), logptr,
		"Receive message operation for DOWN1");
check_flag(
		MPI_Irecv(&recvbuf[DOWN2][0], grid.num_elements_recv_down, MPI_DOUBLE,
				source[DOWN], tag_2, grid.comm, &reqs[8 + DOWN2]), logptr,
		"Receive message operation for DOWN2");

check_flag(
		MPI_Irecv(&recvbuf[LEFT1][0], grid.num_elements_recv_left, MPI_DOUBLE,
				source[LEFT], tag_1, grid.comm, &reqs[8 + LEFT1]), logptr,
		"Receive message operation for LEFT1");
check_flag(
		MPI_Irecv(&recvbuf[LEFT2][0], grid.num_elements_recv_left, MPI_DOUBLE,
				source[LEFT], tag_2, grid.comm, &reqs[8 + LEFT2]), logptr,
		"Receive message operation for LEFT2");

check_flag(
		MPI_Irecv(&recvbuf[RIGHT1][0], grid.num_elements_recv_right, MPI_DOUBLE,
				source[RIGHT], tag_1, grid.comm, &reqs[8 + RIGHT1]), logptr,
		"Receive message operation for RIGHT1");
check_flag(
		MPI_Irecv(&recvbuf[RIGHT2][0], grid.num_elements_recv_right, MPI_DOUBLE,
				source[RIGHT], tag_2, grid.comm, &reqs[8 + RIGHT2]), logptr,
		"Receive message operation for RIGHT2");
check_flag(
		MPI_Isend(sendbuf[UP1], grid.num_elements_send_up, MPI_DOUBLE,
				dest[UP], tag_1, grid.comm, &reqs[UP1]), logptr,
		"Send message operation for UP1");
check_flag(
		MPI_Isend(sendbuf[UP2], grid.num_elements_send_up, MPI_DOUBLE,
				dest[UP], tag_2, grid.comm, &reqs[UP2]), logptr,
		"Send message operation for UP2");

check_flag(
		MPI_Isend(sendbuf[DOWN1], grid.num_elements_send_down, MPI_DOUBLE,
				dest[DOWN], tag_1, grid.comm, &reqs[DOWN1]), logptr,
		"Send message operation for DOWN1");
check_flag(
		MPI_Isend(sendbuf[DOWN2], grid.num_elements_send_down, MPI_DOUBLE,
				dest[DOWN], tag_2, grid.comm, &reqs[DOWN2]), logptr,
		"Send message operation for DOWN2");

check_flag(
		MPI_Isend(sendbuf[LEFT1], grid.num_elements_send_left, MPI_DOUBLE,
				dest[LEFT], tag_1, grid.comm, &reqs[LEFT1]), logptr,
		"Send message operation for LEFT1");
check_flag(
		MPI_Isend(sendbuf[LEFT2], grid.num_elements_send_left, MPI_DOUBLE,
				dest[LEFT], tag_2, grid.comm, &reqs[LEFT2]), logptr,
		"Send message operation for LEFT2");

check_flag(
		MPI_Isend(sendbuf[RIGHT1], grid.num_elements_send_right, MPI_DOUBLE,
				dest[RIGHT], tag_1, grid.comm, &reqs[RIGHT1]), logptr,
		"Send message operation for RIGHT1");
check_flag(
		MPI_Isend(sendbuf[RIGHT2], grid.num_elements_send_right, MPI_DOUBLE,
				dest[RIGHT], tag_2, grid.comm, &reqs[RIGHT2]), logptr,
		"Send message operation for RIGHT2");

	MPI_Waitall(16, reqs, stats);

int remote_source[4],remote_dest[4];
MPI_Request		remote_req[8];
MPI_Status		remote_status[8];

if ( (grid.branch_tag==P)||(grid.branch_tag==L)||(grid.branch_tag==R) ){
	for (int i=0; i<4 i++){
		remote_source[i]	=	nbrs[remote][i];
		remote_dest[i]		=	nbrs[remote][i];
	}

	check_flag(
			MPI_Irecv(&recvbuf[UP1][0], grid.num_elements_recv_right, MPI_DOUBLE,
					remote_source[UP1], tag_remote, grid.universe, &reqs[4 + UP1]), logptr,
			"Receive message operation for remote UP1");
	check_flag(
				MPI_Irecv(&recvbuf[UP2][0], grid.num_elements_recv_right, MPI_DOUBLE,
						remote_source[UP2], tag_remote, grid.universe, &reqs[4 + UP2]), logptr,
				"Receive message operation for remote UP2");
	check_flag(
				MPI_Irecv(&recvbuf[DOWN1][0], grid.num_elements_recv_right, MPI_DOUBLE,
						remote_source[DOWN1], tag_remote, grid.universe, &reqs[4 + DOWN1]), logptr,
				"Receive message operation for remote DOWN1");
	check_flag(
				MPI_Irecv(&recvbuf[DOWN2][0], grid.num_elements_recv_right, MPI_DOUBLE,
						remote_source[DOWN2], tag_remote, grid.universe, &reqs[4 + DOWN2]), logptr,
				"Receive message operation for remote DOWN2");


	check_flag(
			MPI_Isend(sendbuf[UP1], grid.num_elements_send_right, MPI_DOUBLE,
					remote_dest[UP1], tag_remote, grid.universe, &reqs[UP1]), logptr,
			"Send message operation for remote UP1");
	check_flag(
				MPI_Isend(sendbuf[UP2], grid.num_elements_send_right, MPI_DOUBLE,
						remote_dest[UP2], tag_remote, grid.universe, &reqs[UP2]), logptr,
				"Send message operation for remote UP2");
	check_flag(
				MPI_Isend(sendbuf[DOWN1], grid.num_elements_send_right, MPI_DOUBLE,
						remote_dest[DOWN1], tag_remote, grid.universe, &reqs[DOWN1]), logptr,
				"Send message operation for remote DOWN1");
	check_flag(
				MPI_Isend(sendbuf[DOWN2], grid.num_elements_send_right, MPI_DOUBLE,
						remote_dest[DOWN2], tag_remote, grid.universe, &reqs[DOWN2]), logptr,
				"Send message operation for remote DOWN2");

}

	communication_update_recvbuf(logptr);
}//end of update_async()


/*******************************************************************************************/
void communication_update_sendbuf(FILE* logptr)
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


/*******************************************************************************************/
void communication_update_recvbuf(FILE* logptr)
/*******************************************************************************************/
{
	int k, buf_offset;
///UP direction	///
///UP1
	buf_offset = grid.added_info_in_send_buf;
//Updating RECEIVE BUFFERwith SMC information  to be sent in UP1 direction with corresponding downstream neighbour cells located locally.
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
		int start = (int) recvbuf[DOWN2][0], end = (int) recvbuf[DOWN2][1]
		for (int i = end; i <= start; i--) {
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
		for (int i = end ; i <= start;	i--) {
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
if (grid.flip_array[DOWN2]==0){
	for (int i = (int) recvbuf[DOWN2][0]; i <= (int) recvbuf[DOWN2][1]; i++) {
		int j = grid.num_smc_axially + 1;
		smc[i][j].p[smc_Ca] = recvbuf[DOWN2][buf_offset + k + 0];
		smc[i][j].p[smc_Vm] = recvbuf[DOWN2][buf_offset + k + 1];
		smc[i][j].p[smc_IP3] = recvbuf[DOWN2][buf_offset + k + 2];
		k += grid.num_coupling_species_smc;
	}
}
else if (grid.flip_array[DOWN2]==1){
	int start = (int) recvbuf[DOWN1][0] , end = (int) recvbuf[DOWN1][1];
	for (int i = end ; i <= start; i--) {
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
if (grid.flip_array[DOWN2] == 0){
	for (int i = (int) recvbuf[DOWN2][2]; i <= (int) recvbuf[DOWN2][3]; i++) {
		int j = grid.num_ec_axially + 1;
		ec[i][j].q[ec_Ca] = recvbuf[DOWN2][buf_offset + k + 0];
		ec[i][j].q[ec_Vm] = recvbuf[DOWN2][buf_offset + k + 1];
		ec[i][j].q[ec_IP3] = recvbuf[DOWN2][buf_offset + k + 2];
		k += grid.num_coupling_species_ec;
	}
}
else if (grid.flip_array[DOWN2]==1){
	int start =  (int) recvbuf[DOWN1][2], end = (int) recvbuf[DOWN1][3];
	for (int i = end; i <= start; i--) {
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
		int 	i = grid.num_smc_circumferentially + 1;
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
		int 	i = grid.num_ec_circumferentially + 1;
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
		int 	i = grid.num_smc_circumferentially + 1;
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
		int 	i = grid.num_ec_circumferentially + 1;
		ec[i][j].q[ec_Ca] = recvbuf[RIGHT2][buf_offset + k + 0];
		ec[i][j].q[ec_Vm] = recvbuf[RIGHT2][buf_offset + k + 1];
		ec[i][j].q[ec_IP3] = recvbuf[RIGHT2][buf_offset + k + 2];
		k += grid.num_coupling_species_ec;
	}
}	// end of communication_update_recvbuf()


void print_domains(FILE* logptr){

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
*/

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
	}
}

void print_send_buffer(FILE* logptr){

	fprintf(logptr, "***Up direction***\n");
	for (int i = 0; i < grid.num_elements_send_up; i++) {
		fprintf(logptr, "[%d]\t %2.3lf \t %2.3lf\n", i, sendbuf[UP1][i],
				sendbuf[UP2][i]);
	}
	fprintf(logptr, "***Down direction***\n");
	for (int i = 0; i < grid.num_elements_send_down; i++) {
		fprintf(logptr, "[%d]\t %2.3lf \t %2.3lf\n", i, sendbuf[DOWN1][i],
				sendbuf[DOWN2][i]);
	}
	fprintf(logptr, "***Left direction***\n");
	for (int i = 0; i < grid.num_elements_send_left; i++) {
		fprintf(logptr, "[%d]\t %2.3lf \t %2.3lf\n", i, sendbuf[LEFT1][i],
				sendbuf[LEFT2][i]);
	}
	fprintf(logptr, "***Right direction***\n");
	for (int i = 0; i < grid.num_elements_send_right; i++) {
		fprintf(logptr, "[%d]\t %2.3lf \t %2.3lf\n", i, sendbuf[RIGHT1][i],
				sendbuf[RIGHT2][i]);
	}
}


void print_recv_buffer(FILE* logptr){

fprintf(logptr, "***Up direction***\n");
	for (int i = 0; i < grid.num_elements_send_up; i++) {
		fprintf(logptr, "[%d]\t %2.3lf \t %2.3lf\n", i, recvbuf[UP1][i],
				recvbuf[UP2][i]);
	}
	fprintf(logptr, "***Down direction***\n");
	for (int i = 0; i < grid.num_elements_send_down; i++) {
		fprintf(logptr, "[%d]\t %2.3lf \t %2.3lf\n", i, recvbuf[DOWN1][i],
				recvbuf[DOWN2][i]);
	}
	fprintf(logptr, "***Left direction***\n");
	for (int i = 0; i < grid.num_elements_send_left; i++) {
		fprintf(logptr, "[%d]\t %2.3lf \t %2.3lf\n", i, recvbuf[LEFT1][i],
				recvbuf[LEFT2][i]);
	}
	fprintf(logptr, "***Right direction***\n");
	for (int i = 0; i < grid.num_elements_send_right; i++) {
		fprintf(logptr, "[%d]\t %2.3lf \t %2.3lf\n", i, recvbuf[RIGHT1][i],
				recvbuf[RIGHT2][i]);
	}
}

void print_compare(FILE* logptr, double t, double y[],celltype1** smc, celltype2** ec){
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








