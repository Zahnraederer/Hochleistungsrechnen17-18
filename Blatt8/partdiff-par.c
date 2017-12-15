#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

#include "partdiff-par.h"
//#include "partdiff-seq.c"

//#ifndef PI
//#define PI 			3.141592653589793
//#endif
//#define TWO_PI_SQUARE 		(2 * PI * PI)


//Erstellt einen Speicherblock für double der Dimension N*k
double** createMatrix(uint64_t N, uint64_t k){
  //printf("In createMatrix mit N:%ld k:%ld\n",N,k);
  uint64_t i,j;
  
  double* M = malloc( (N) * (k) * sizeof(double));
  if(M == NULL) exit(EXIT_FAILURE);
  
  double** Matrix = malloc((k) * sizeof(double*));
  if(Matrix == NULL) exit(EXIT_FAILURE);

  //printf("Vor der Adressschleife\n");
  for (j = 0; j < k; j++){
    //Matrix[j] = M[j*N];
    Matrix[j] = M + j*N;
  }
  
  //printf("Nach der Adressschleife vor der Fuellschleif\n");
  
  //int start = rank*(N*k);
  
  /* initialize matrix/matrices with zeros */

  for (i = 0; i < k*N; i++){
    //for (j = 0; j < N; j++){
      //Matrix[i][j] = 0.0;
    //}
    M[i] = i * 1.0;
  }
  //printf("Nach der Fuellschleif\n");
  return Matrix;
}

double** init (uint64_t N,uint64_t k,uint64_t rest,int rank,uint64_t* newK,int
               size){
  
  //printf("In init mit N:%ld k:%ld rest:%ld rank:%d\n",N,k,rest,rank);
  double** Matrix;

  //printf("Init in Rank %d mit k %d und rest %d\n",rank,k,rest);
  //srand(time(NULL));
  
  //All Ranks smaller then rest get an extra line to get N elements
  if(rest > (uint64_t) rank)
    *newK = k+1+2;
    //Matrix = createMatrix(N,k+1);
  else
    *newK = k+2;
  
  if(rank == 0 || rank == size-1)
    *newK -=1 ;
  
  Matrix = createMatrix(N,*newK);
  
  return Matrix;
}


int main (int argc, char** argv){
  MPI_Init(&argc,&argv);
  
  struct options options;
  //struct calculation_arguments arguments;
  //struct calculation_results results;
  
  //char arg[256];
  //char argIter[256];
  uint64_t N;
  uint64_t k;
  uint64_t iter;
  uint64_t i,j,l;
  int rank;
  int size;
  double** Matrix_In;
  double** Matrix_Out;
  double** tmp;
  double star;
  //double PI = 3.141592653589793;
  //double TWO_PI_SQUARE = (2 * PI * PI);
  MPI_Request req;

  //int* foo = allocateMemory(5);
  
  //size
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if(size == 1){
    //printf("Nothing to be done for a single Thread!\n");
    seq_Call(argc,argv);
    MPI_Finalize();
    return EXIT_SUCCESS;
  }
  
  
  //myrank
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //Haben wir weniger als 7 Parameter haben...
  if (argc < 7){
    //Muss der erste nachfragen
    if(rank == 0){
      AskParams(&options, argc, argv);
    }
  }
  //Anderenfalls kann AskParams benutzt werden zum Parsen in jedem Prozess
  else{
    AskParams(&options, argc, argv);
  }

  //Bei Gauß-Seidel nur sequentiell
  if(options.method == 1){
    if(rank == 0)
      seq_Call(argc,argv);

    MPI_Finalize();
    return EXIT_SUCCESS;
  }
  
  /*
  struct options{
  uint64_t number;         // Number of threads                              /
  uint64_t method;         // Gauss Seidel or Jacobi method of iteration     /
  uint64_t interlines;     // matrix size = interlines*8+9                   /
  uint64_t inf_func;       // inference function                             /
  uint64_t termination;    // termination condition                          /
  uint64_t term_iteration; // terminate if iteration number reached          /
  double   term_precision; // terminate if precision reached                 /
  };
  */
  
  
  // LineLength as usual N * 8 + 9 - 1
  N =  (options.interlines+1)*8;//atoi(arg);
  iter = options.term_iteration;//atoi(argIter);
  
  //Abgerundete Elementzahl
  uint64_t segments = N/size;
  uint64_t rest = N%size;

  Matrix_In = init(N,segments,rest,rank,&k,size);
  Matrix_Out = init(N,segments,rest,rank,&k,size);
  //Er kann auf elemente hinter der letzten Spalte zugreifen und die sind 0.0
  //printf("Matrix_In[1][1] = %f Matrix_In[k-2][N-2] = %f Mit k = %ld und N=%ld
  //\n",Matrix_In[1][1],Matrix_In[k-2][N-2],k,N);

  double h = 1.0/(double)N;
  double pih = PI * h;
  double fpisin = 0.25 * TWO_PI_SQUARE * h * h;
  //printf("Iter ist:%lu h:%f  pih:%f fpisin:%f \n",iter,h,pih,fpisin);
    
  for(i = 0; i <= iter; i++){
    
    //over all rows
    for (j = 1; j <= k-2; j++){
      double fpisin_i = 0.0;
      fpisin_i = fpisin * sin(pih * (double)j);
      
      //over all columns
      for (l = 1; l <= N-2; l++){
        star = 0.25 * (Matrix_In[j-1][l] + Matrix_In[j][l-1] + Matrix_In[j][l+1]
                       + Matrix_In[j+1][l]);
        star += fpisin_i * sin(pih * (double)l);
        //if(!(iter % 25))
          //printf("star: %f sin(x): %f\n",star,sin(pih * (double)j));
        Matrix_Out[j][l] = star;
      }
      
      if(j == 1){
        //Sende erste Zeile an vorgaenger wenn wir nicht 0 sind
        if(rank != 0){
          MPI_Isend(Matrix_Out[j],N,MPI_DOUBLE, rank-1,
                    0,MPI_COMM_WORLD, &req);
        //MPI_Send(buf,k+rest,MPI_INT,rank+1,0,MPI_COMM_WORLD);
          printf("Rank %d hat Nach oben gesendet!\n",rank);
        }
      }
      
      
      if(j == k-3){
        //empfange von nachfolger
        if(rank != size-1){
          MPI_Recv(Matrix_In[k-1],N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
        
          printf("Rank %d hat von unten empfangen!\n",rank);
        }
      }
      
      
      if(j == k-2){
        //Sende an Nachfolger
        if(rank != size-1){
          MPI_Isend(Matrix_Out[k-2],N,MPI_DOUBLE, rank+1,
                    0,MPI_COMM_WORLD, &req);
        
        
          printf("Rank %d hat nach unten gesendet!\n",rank);
        }
        if(rank != 0){
        //Empfange von Vorgaenger
          MPI_Recv(Matrix_In[0],N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
        
          printf("Rank %d hat von oben empfangen!\n",rank);
        }
      }
    }
    tmp = Matrix_In;
    Matrix_In = Matrix_Out;
    Matrix_Out = tmp;
  }
  
  printf("Matrix_Out[1][1] = %f Matrix_Out[k-2][N-2] = %f Mit k = %ld und N=%ld\n",Matrix_Out[1][1],Matrix_Out[k-2][N-2],k,N);
  
  
  free(Matrix_Out[0]);
  free(Matrix_Out);
  free(Matrix_In[0]);
  free(Matrix_In);
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}
