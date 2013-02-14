#include "mpi.h"
#include <stdio.h>

#define true	1
#define false  0
void check_flag(int err, FILE* errfile, const char* errmsg, int rank){
	if (err !=MPI_SUCCESS) {
		fprintf(errfile,"[%d]\t MPI-Comm Error:%s\n",rank,errmsg);
		MPI_Abort(MPI_COMM_WORLD, 200);
	}
}
int main( int argc, char *argv[] )
{
	int ndims = 3;
	int errs = 0;
    int size, dims[3], periods[3], remain[3];
    int result;
    int coords[2];
    int rank,new_rank;
    MPI_Comm comm, newcomm;
    MPI_Init( &argc, &argv );

    /* First, create a 1-dim cartesian communicator */
    periods[0] = 0; periods[1] = 0; periods[2] = 0;
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    dims[0] = 3;dims[1] = 2; dims[2] = 4;

    check_flag(MPI_Cart_create( MPI_COMM_WORLD, ndims, dims, periods, 0, &comm ),stdout,"Cart_create failed",rank);

    /* Now, extract a communicator with no dimensions */
    remain[0] = false;remain[1] = true; remain[2] = true;
    MPI_Cart_sub( comm, remain, &newcomm );

    /* This should be congruent to MPI_COMM_SELF
    MPI_Comm_compare( MPI_COMM_SELF, newcomm, &result );
    if (result != MPI_CONGRUENT) {
        errs++;
        printf( "cart sub to size 0 did not give self\n" );fflush(stdout);
    }*/
    check_flag(MPI_Comm_rank(newcomm,&new_rank),stdout,"Comm_new_rank failed",rank);

    check_flag(MPI_Cart_coords(newcomm, new_rank , 2, coords),stdout,"Cart_sub failed",rank);


    printf("[%d,%d] :coords(%d,%d)\n",rank,new_rank,coords[0],coords[1]);
    /* Free the new communicator */
    MPI_Comm_free( &newcomm );
    MPI_Comm_free( &comm );

    MPI_Finalize();
    return 0;
}
