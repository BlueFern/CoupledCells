#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define		UP		0
#define		DOWN	1
#define 	LEFT	2
#define		RIGHT	3

void check_flag(int err, FILE* errfile, const char* errmsg, int rank){
	if (err !=MPI_SUCCESS) {
		fprintf(errfile,"[%d]\t MPI-Comm Error:%s\n",rank,errmsg);
		MPI_Abort(MPI_COMM_WORLD, 200);
	}
}
int main(int argc, char *argv[])
{

	/* Starts MPI processes ... */
	      MPI_Init(&argc, &argv);
	      int rank,size, new_rank,new_size, cart_rank,cart_size;
	  	int
	  		m	=	6,
	  		n	=	10;

	      int nbrs[4], dims[2], periods[2], reorder=0, coords[2];
	      int outbuf, inbuf[4]={MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL},
	      	source,dest;

	      MPI_Request	reqs[8];
	      MPI_Status	status[8];

	      MPI_Comm	newcomm,comm;
	      MPI_Comm_rank(MPI_COMM_WORLD, &rank);  /* get current process id */
	      MPI_Comm_size(MPI_COMM_WORLD, &size);    /* get number of processes */

	      char filename[30];
	      sprintf(filename, "log%d", rank);
	      FILE *log;
	      log = fopen(filename, "w+");

	      int color = int(rank/(m*n));
	      int key	=	color*((m*n)-1);

	     check_flag(MPI_Comm_split(MPI_COMM_WORLD,color,key,&newcomm),stdout, "Comm-split failed",rank);
	     MPI_Comm_rank(newcomm, &new_rank);  /* get current process id */
	     MPI_Comm_size(newcomm, &size);    /* get number of processes */
	     //     if((color==0)||(color==2))
		//if (color==1)
	     	//     fprintf(log,"[%d] color=%d, key=%d, comm_size = %d\t new_rank=%d\n",rank,color,key,size,new_rank);
	     //else if (color==1)
	    //	 	 fprintf(stdout,"[%d] color=%d, key=%d, comm_size = %d\t new_rank=%d\n",rank,color,key,size,new_rank);

	     	dims[0]		=m;
	     	dims[1]		=n;
	     	periods[0]	=0;
	     	periods[1]	=1;
	     	reorder		=0;

	     check_flag(MPI_Cart_create(newcomm, 2, dims, periods, reorder,
					&comm), stdout, "failed at cart create",rank);
	     check_flag(MPI_Comm_rank(comm, &cart_rank), stdout,
			"failed at comm rank",rank);
	     check_flag(MPI_Cart_coords(comm, cart_rank, 2, coords), stdout,
			"failed at cart coords",rank);


	     MPI_Barrier(MPI_COMM_WORLD);



	     check_flag(MPI_Cart_shift(comm, 0, 1, &nbrs[UP], &nbrs[DOWN]),stdout,"failed at cart shift up down",rank);
	     check_flag(MPI_Cart_shift(comm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]),stdout,"failed at cart left right",rank);

	     int tag	=	1;

	     for (int i=0; i<4; i++) {
	          dest = nbrs[i];
	          source = nbrs[i];
	          MPI_Isend(&outbuf, 1, MPI_INT, dest, tag, comm, &reqs[i]);
	          MPI_Irecv(&inbuf[i], 1, MPI_INT, source, tag, comm, &reqs[i+4]);
	          }
	       MPI_Waitall(8, reqs, status);
	       fprintf(log,"\n[%d -> %d -> %d]  coords (%d , %d)\t (u,d,l,r) (%d, %d, %d, %d)\n",rank,new_rank,cart_rank,coords[0],coords[1],nbrs[UP],nbrs[DOWN],nbrs[LEFT],nbrs[RIGHT]);

	     fclose(log);

	     MPI_Finalize();

}
