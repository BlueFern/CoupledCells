#include <mpi.h>
int main(int argc, char *argv[])
{
      int mcol, irow, jcol, p;
      MPI_Comm row_comm, col_comm, comm2D;
      int Iam, row_id, col_id;
      int row_group, row_key, map[6];

/* Starts MPI processes ... */
      MPI_Init(&argc, &argv);                          /* starts MPI */
      MPI_Comm_rank(MPI_COMM_WORLD, &Iam);  /* get current process id */
      MPI_Comm_size(MPI_COMM_WORLD, &p);    /* get number of processes */
     
      map[0]=2; map[1]=1; map[2]=2; map[3]=1; map[4]=0; map[5]=1;
      mcol=2;   /* nrow = 3 */
      if(Iam == 0) {
        printf("\n");
        printf("Example of MPI_Comm_split Usage\n");
        printf("Split 3x2 grid into 2 different communicators\n");
        printf("which correspond to 3 rows and 2 columns.");
        printf("\n");
        printf("     Iam     irow     jcol   row-id   col-id\n");
      }
/* virtual topology with nrow rows and mcol columns */
      irow = Iam/mcol;        /* row number */
      jcol = Iam%mcol;        /* column number */
      comm2D = MPI_COMM_WORLD;
      MPI_Comm_split(comm2D, irow, jcol, &row_comm);
      MPI_Comm_split(comm2D, jcol, irow, &col_comm);

      MPI_Comm_rank(row_comm, &row_id);
      MPI_Comm_rank(col_comm, &col_id);
      MPI_Barrier(MPI_COMM_WORLD);

      printf("%8d %8d %8d %8d %8d\n",Iam,irow,jcol,row_id,col_id);
      MPI_Barrier(MPI_COMM_WORLD);

      if(Iam == 0) {
        printf("\n");
        printf("Next, create more general communicator\n");
        printf("which consists of two groups :\n");
        printf("Rows 1 and 2 belongs to group 1 and row 3 is group 2\n");
        printf("\n");
      }
/* MPI_Comm_split is more general than MPI_Cart_sub
   simple example of 6 processes divided into 2 groups;
   1st 4 belongs to group 1 and remaining two to group 2 */
      row_group = Iam/4;     /* this expression by no means general */
      row_key = Iam - row_group*4;    /* group1:0,1,2,3; group2:0,1 */
      MPI_Comm_split(comm2D, row_group, row_key, &row_comm);
      MPI_Comm_rank(row_comm, &row_id);
      printf("%8d %8d\n",Iam,row_id);
      MPI_Barrier(MPI_COMM_WORLD);

      if(Iam == 0) {
        printf("\n");
        printf("If two processes have same key, the ranks\n");
        printf("of these two processes in the new\n");
        printf("communicator will be ordered according'\n");
        printf("to their order in the old communicator\n");
        printf(" key = map[Iam]; map = (2,1,2,1,0,1)\n");
        printf("\n");
      }
/* MPI_Comm_split is more general than MPI_Cart_sub
   simple example of 6 processes dirowided into 2 groups;
   1st 4 belongs to group 1 and remaining two to group 2 */
      row_group = Iam/4;     /* this expression by no means general */
      row_key = map[Iam];
      MPI_Comm_split(comm2D, row_group, row_key, &row_comm);
      MPI_Comm_rank(row_comm, &row_id);
      MPI_Barrier(MPI_COMM_WORLD);
      printf("%8d %8d %8d\n",Iam,row_id,row_group);

      MPI_Finalize();                  /* let MPI finish up ...  */
}
