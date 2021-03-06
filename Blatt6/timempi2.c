#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <limits.h>

typedef struct {
  char hostname[HOST_NAME_MAX];
  struct timeval time;
} timempi_res_t;

int main(int argc, char** argv) {
  timempi_res_t package;
  struct timeval time;
  time_t curtime;
  char buffer[120];
  MPI_Init(NULL, NULL);
  
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // Get the rank of the process
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  long int usec = 0;
 
  if(rank == 0){
    for(int i = 1; i < size;i++){
      MPI_Recv(&package,sizeof(package),MPI_BYTE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      curtime=package.time.tv_sec;
      strftime(buffer,120,"%F %T",localtime(&curtime));
      printf("%s: %s.%ld\n",
             package.hostname, buffer,package.time.tv_usec);
    }
  } 
  else{
    int name_len;
    MPI_Get_processor_name(package.hostname, &name_len);
    gettimeofday(&time,NULL);
    package.time = time;
    MPI_Send(&package,sizeof(package),MPI_BYTE,0,0,MPI_COMM_WORLD);
  }
  MPI_Reduce(&package.time.tv_usec, &usec,1,MPI_LONG,MPI_MIN,0,MPI_COMM_WORLD);
  if(rank == 0){
    printf("%ld\n",usec);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Rang %d endet jetzt!\n", rank);
  // Finalize the MPI environment.
  MPI_Finalize();
}
