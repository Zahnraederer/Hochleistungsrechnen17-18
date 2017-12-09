#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

int*
init (int k,int rest,int rank)
{
  int* buf = malloc(sizeof(int) * (k+1));

  //printf("Init in Rank %d mit k %d und rest %d\n",rank,k,rest);
  srand(time(NULL));
  
  //All Ranks smaller then rest get an extra element to get N elements
  if(rest > rank){
    for (int i = 0; i < k+1; i++){
      //DO not modify % 25
      buf[i] = abs((rand()*rank+rand()) % 25);
    }
  }
  else{
    for (int i = 0; i < k+1; i++){
      if(i < k)
        // Do not modify % 25
        buf[i] = abs((rand()*rank+rand()) % 25);
      else
        //-1 prevents an element from beeing printed
        buf[i] = -1;
    }
  }

  return buf;
}

int*
circle (int* buf,int* newBuf,int size,int rank,int k,int rest)
{
  //the last rank must send to the first
  if(rank == size-1)
    MPI_Send(buf,k+rest,MPI_INT,0,0,MPI_COMM_WORLD);
  //all other ranks send to the next
  else
    MPI_Send(buf,k+rest,MPI_INT,rank+1,0,MPI_COMM_WORLD);
  
  //the first rank recevies from the last
  if(rank != 0)
    MPI_Recv(newBuf,k+rest,MPI_INT,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  //all other ranks receive from the previous rank
  else
    MPI_Recv(newBuf,k+rest,MPI_INT,size-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
  //Ensure that exactly one iteration was done
  MPI_Barrier(MPI_COMM_WORLD);

  return newBuf;
}

int
main (int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  
  char arg[256];
  int N;
  int rank;
  int size;
  int* buf;
  //Each rank has a placeholderspace for the new Elements
  int* newBuf;
  //The last rank must hold the first element to trigger the break
  int firstElement;
  //All other ranks must receive an abort signal
  int abort = 0;

  if (argc < 2)
  {
    printf("Arguments error!\n");
    return EXIT_FAILURE;
  }

  sscanf(argv[1], "%s", arg);

  //size
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if(size == 1){
    printf("Nothing to be done for a single Thread!\n");
    return EXIT_SUCCESS;
  }
  //myrank
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  // Array length
  N = atoi(arg);
  //Abgerundete Elementzahl
  int k = N/size;
  int rest = N%size;

  buf = init(k,rest,rank);
  newBuf = malloc(sizeof(int) * (k+1));
  //Rest will be 1 to be consistent with a previous implementation
  rest = 1;

  if(rank == 0)
    MPI_Send(buf,1,MPI_INT,size-1,0,MPI_COMM_WORLD);
  if(rank == size-1)
    MPI_Recv(&firstElement,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
 
  //First the first rank will print all buffers to the console
  if(rank == 0){
    printf("\nBEFORE\n");
    for (int i = 0; i < k+rest; i++){
      if(buf[i] != -1)
        printf ("rank %d: %d\n", rank, buf[i]);
    }
    for (int j = 1; j < size; j++){
      MPI_Recv(newBuf,k+rest,MPI_INT,j,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      for (int i = 0; i < k+rest; i++){
        if(newBuf[i] != -1)
          printf ("rank %d: %d\n", j, newBuf[i]);
      }
    }
  }
  else
    MPI_Send(buf,k+rest,MPI_INT,0,0,MPI_COMM_WORLD);

  //While the last rank does't hold the first element we iterate
  while(1){
    if(rank == size-1){
      if(buf[0] == firstElement){
        abort = 1;
        MPI_Bcast(&abort,1,MPI_INT,size-1,MPI_COMM_WORLD);
        break;
      }
      else{
        abort = 0;
        MPI_Bcast(&abort,1,MPI_INT,size-1,MPI_COMM_WORLD);
      }
    }
    else
      MPI_Bcast(&abort,1,MPI_INT,size-1,MPI_COMM_WORLD);
    if(abort == 1)
      break;
    buf = circle(buf,newBuf,size,rank,k,rest);
  }
  
  //Again the first rank will print all buffers to the console
  if(rank == 0){
    printf("\nAFTER\n");
    for (int i = 0; i < k+rest; i++)
    {
      if(buf[i] != -1)
        printf ("rank %d: %d\n", rank, buf[i]);
    }
    
    for (int j = 1; j < size; j++){
      MPI_Recv(newBuf,k+rest,MPI_INT,j,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      for (int i = 0; i < k+rest; i++)
      {
        if(newBuf[i] != -1)
          printf ("rank %d: %d\n", j, newBuf[i]);
      }
    }
  }
  else{
    MPI_Send(buf,k+rest,MPI_INT,0,1,MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
