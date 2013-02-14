#include "mpi.h"
int main(int argc, char *argv[])
{
      int Iam, p;
      int Neven, Nodd, members[6], even_rank, odd_rank;
      MPI_Group group_world, even_group, odd_group;
/* Starts MPI processes ... */
      MPI_Init(&argc, &argv);                          /* starts MPI */
      MPI_Comm_rank(MPI_COMM_WORLD, &Iam);  /* get current process id */
      MPI_Comm_size(MPI_COMM_WORLD, &p);    /* get number of processes */
      Neven = (p + 1)/2;     /* All processes of MPI_COMM_WORLD are divided */
      Nodd = p - Neven;      /* into 2 groups, odd- and even-numbered groups */
      members[0] = 2;
      members[1] = 0;
      members[2] = 4;
      MPI_Comm_group(MPI_COMM_WORLD, &group_world);
      MPI_Group_incl(group_world, Neven, members, &even_group);
      MPI_Group_excl(group_world, Neven, members,  &odd_group);

      MPI_Barrier(MPI_COMM_WORLD);
      if(Iam == 0) {
        printf("MPI_Group_incl/excl Usage Example\n");
        printf("\n");
        printf("Number of processes is %d\n", p);
        printf("Number of odd processes is %d\n", Nodd);
        printf("Number of even processes is %d\n", Neven);
        printf("\n");
        printf("     Iam     even      odd\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);

      MPI_Group_rank(even_group, &even_rank);
      MPI_Group_rank( odd_group,  &odd_rank);
      printf("%8d %8d %8d\n",Iam, even_rank, odd_rank);

      MPI_Finalize();                  /* let MPI finish up ...  */
}     
