#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

//Erstellt einen Speicherblock für double der Dimension N*k
double** createMatrix(uint64_t N, uint64_t k){
  printf("In createMatrix mit N:%ld k:%ld\n",N,k);
  uint64_t i,j;
  
  double* M = malloc( (N) * (k) * sizeof(double));
  if(M == NULL) exit(EXIT_FAILURE);
  
  double** Matrix = malloc((k) * sizeof(double*));
  if(Matrix == NULL) exit(EXIT_FAILURE);

  //printf("Vor der Adressschleife\n");
  for (j = 0; j < k; j++){
    //Matrix[j] = M[j*N];
    Matrix[j] = M + j*k;
  }
  
  //printf("Nach der Adressschleife vor der Fuellschleif\n");
  
  /* initialize matrix/matrices with zeros */

  for (i = 0; i < k; i++){
    for (j = 0; j < N; j++){
      Matrix[i][j] = 0.0;
    }
  }
  //printf("Nach der Fuellschleif\n");
  return Matrix;
}

double** init (uint64_t N,uint64_t k,uint64_t rest,int rank,uint64_t* newK){
  
  printf("In init mit N:%ld k:%ld rest:%ld rank:%d\n",N,k,rest,rank);
  double** Matrix;

  //printf("Init in Rank %d mit k %d und rest %d\n",rank,k,rest);
  //srand(time(NULL));
  
  //All Ranks smaller then rest get an extra line to get N elements
  if(rest > (uint64_t) rank)
    *newK = k+1;
    //Matrix = createMatrix(N,k+1);
  else
    *newK = k;
  
  Matrix = createMatrix(N,*newK);
  
  return Matrix;
}

int main (int argc, char** argv){
  MPI_Init(&argc,&argv);
  
  char arg[256];
  uint64_t N;
  uint64_t k;
  int rank;
  int size;
  double** Matrix;

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
  uint64_t segments = N/size;
  uint64_t rest = N%size;

  Matrix = init(N,segments,rest,rank,&k);
  printf("Matrix[0][0] = %f Matrix[k-1][N] = %f Mit k = %ld\n",Matrix[0][0],Matrix[k-1][N+1],k);
  
  free(Matrix[0]);
  free(Matrix);
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}
