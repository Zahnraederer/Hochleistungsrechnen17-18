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


//------------------------------------------------------------------------------

void DisplayMatrix (/*struct calculation_arguments* arguments, struct calculation_results* results, */struct options* options, int rank, int size, int from, int to,double** Matrix_In)
{
  int const elements = 8 * options->interlines + 9;
  
  int x, y;
  double** Matrix = Matrix_In;
  MPI_Status status;
  
  // first line belongs to rank 0
  if (rank == 0){
    from--;
  }
  
  // last line belongs to rank size - 1
  if (rank + 1 == size){
    to++;
  }
  
  printf("Display Matrix rank = %d from = %d to = %d\n",rank,from,to);
  if (rank == 0)
    printf("Matrix:\n");
  
  for (y = 0; y < 9; y++)
  {
    int line = y * (options->interlines + 1);
    
    if (rank == 0)
    {
      // check whether this line belongs to rank 0
      if (line < from || line > to)
      {
        // use the tag to receive the lines in the correct order
        // the line is stored in Matrix[0], because we do not need it anymore
        MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
      }
    }
    else
    {
      if (line >= from && line <= to)
      {
        // if the line belongs to this process, send it to rank 0
        // (line - from + 1) is used to calculate the correct local address
        printf("Sending %f to %f  from rank %d\n",Matrix[line - from + 1][0],Matrix[line - from + 1][options->interlines],rank);

        MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
      }
    }
    
    if (rank == 0)
    {
      for (x = 0; x < 9; x++)
      {
       // int col = x * (options->interlines + 1);
        
        if (line >= from && line <= to)
        {
          // this line belongs to rank 0
         // printf("%7.4f", Matrix[line][col]);
        }
        else
        {
          // this line belongs to another rank and was received above
          //printf("%7.4f", Matrix[0][col]);
        }
      }
      
      printf("\n");
    }
  }
  
  fflush(stdout);
}

//------------------------------------------------------------------------------



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



void calculate_Gaus(uint64_t N,uint64_t iter,double** Matrix_In,uint64_t rows, struct options const* options){
  int initialized = 0;
  int residuum_iter = 0;
  int residuum_total_iter = 0;
  int residuum_rest = 0;
  double residuum;
  double maxresiduum;
  double globalmaxresiduum;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
  //options->term_precision
    
  uint64_t termination = options->termination;
  double star;
  
  MPI_Request SendFirstRowReq;
  MPI_Request SendLastRowReq;
  
  double h = 1.0/(double)N;
  double pih = 0.0;
  double fpisin = 0.0;
  
  uint64_t i,j,k;
  
  if (options->inf_func == FUNC_FPISIN)
  {
    pih = PI * h;
    fpisin = 0.25 * TWO_PI_SQUARE * h * h;
  }
  
  
  if(termination == TERM_ITER){
    for(k = 0; k < iter; k++){
      
      //Schon in der ersten Iteration müssen wir empfangen
      if(rank != 0){
        MPI_Recv(Matrix_In[0],N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      }
      if((initialized == 1) && (rank != 0)){
        //printf("Will jetzt auf den Send warten!\n");
        MPI_Wait(&SendFirstRowReq, MPI_STATUS_IGNORE);
        //printf("Habe korrekt auf die erste Zeilenversendung gewartet!\n");
      }
      
      //Ueber die Zeilen
      for(i = 1; i < rows-2; i++){
        double fpisin_i = 0.0;
        if (options->inf_func == FUNC_FPISIN)
          fpisin_i = fpisin * sin(pih * (double)i);
        //Ueber die Spalten
        for(j = 1; j < N-1; j++){
          star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);
          if (options->inf_func == FUNC_FPISIN)
            star += fpisin_i * sin(pih * (double)j);
        }
        if((i == 1) && (rank != 0)){
          //An den Vorgänger schicken
          //printf("Habe an den vorgaenger abgeschickt!\n");
          MPI_Isend(Matrix_In[1],N,MPI_DOUBLE, rank-1,
                    0,MPI_COMM_WORLD, &SendFirstRowReq);
        }
      }
      //Vor dem berechnen der letzten Reihe wenn nötig empfangen
      if( (initialized == 1) && (rank != size-1) ){
        MPI_Recv(Matrix_In[rows-1],N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Wait(&SendLastRowReq, MPI_STATUS_IGNORE);
      }
      double fpisin_i = 0.0;
      //Die letzte Reihe berechnen
      for(j = 1; j < N-1; j++){
        i = rows-2;
        star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);
        if (options->inf_func == FUNC_FPISIN)
        star += fpisin_i * sin(pih * (double)j);
      }
      if(rank != size-1)
        MPI_Isend(Matrix_In[rows-2],N,MPI_DOUBLE, rank+1,
                  0,MPI_COMM_WORLD, &SendLastRowReq);
      //Nach der ersten Iteration müssen wir immer die letzte reihe empfangen
      initialized = 1;
    }
  }
  
  else if(termination == TERM_PREC){
    //Solange rechnen bis die gleiche Iterationszahl mit richtigem res da ist
    while( residuum_rest != 1 || residuum_iter >= residuum_total_iter){
      
      //Schon in der ersten Iteration müssen wir empfangen
      if(rank != 0)
        MPI_Recv(Matrix_In[0],N,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
      
      if((initialized == 1) && (rank != 0))
        MPI_Wait(&SendFirstRowReq, MPI_STATUS_IGNORE);
      
      //Ueber die Zeilen
      for(i = 1; i < rows-2; i++){
        double fpisin_i = 0.0;
        if (options->inf_func == FUNC_FPISIN)
          fpisin_i = fpisin * sin(pih * (double)i);
        //Ueber die Spalten
        for(j = 1; j < N-1; j++){
          star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);
          if (options->inf_func == FUNC_FPISIN)
            star += fpisin_i * sin(pih * (double)j);
          residuum = Matrix_In[i][j] - star;
          residuum = (residuum < 0) ? -residuum : residuum;
          maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
          }
        if((i == 1) && (rank != 0))
          //An den Vorgänger schicken
          MPI_Isend(Matrix_In[1],N,MPI_DOUBLE, rank-1,
                    0,MPI_COMM_WORLD, &SendFirstRowReq);
      }
      
      //Vor dem berechnen der letzten Reihe wenn nötig empfangen
      if( (initialized == 1) && (rank != size-1) ){
        MPI_Recv(Matrix_In[rows-1],N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Wait(&SendLastRowReq, MPI_STATUS_IGNORE);
      }
      double fpisin_i = 0.0;
      //Die letzte Reihe berechnen
      for(j = 1; j < N-1; j++){
        i = rows-2;
        star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);
        if (options->inf_func == FUNC_FPISIN)
          star += fpisin_i * sin(pih * (double)j);
        residuum = Matrix_In[i][j] - star;
        residuum = (residuum < 0) ? -residuum : residuum;
        maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
      }
      if(rank != size-1)
        MPI_Isend(Matrix_In[rows-2],N,MPI_DOUBLE, rank+1,
                  0,MPI_COMM_WORLD, &SendLastRowReq);
      //Nach der ersten Iteration müssen wir immer die letzte reihe empfangen
      initialized = 1;
      
      residuum_iter++;
      if(!residuum_rest){
        residuum_total_iter++;
        //Wenn alle threads schon laufen
        if( (!residuum_rest) && residuum_iter >= size - rank){
          MPI_Allreduce(&maxresiduum, &globalmaxresiduum, 1, MPI_DOUBLE, MPI_MAX,
                        MPI_COMM_WORLD);
          if((residuum_iter%10000) == 0)
            printf("Warte auf Reduce in iteration: %d bin rank: %d mit maxes: %f\n",residuum_iter,rank,globalmaxresiduum);
          if(globalmaxresiduum < options->term_precision){
            MPI_Allreduce(&residuum_total_iter, &residuum_total_iter, 1, MPI_INT, MPI_MAX,
                          MPI_COMM_WORLD);
            residuum_rest = 1;
          }
        }
      }
    }   
  }
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

  //arguments.N = N;
  //arguments.num_matrices = (options.method == METH_JACOBI) ? 2 : 1;
  //arguments.h = 1.0 / N;
  //arguments.Matrix = &Matrix_In;
  
  //results.m = 0;
  //results.stat_iteration = 0;
  //results.stat_precision = 0;
  
  
  //Bei Gauß-Seidel nur sequentiell
  if(options.method == 1){
    if(rank == 0)
      printf("Gaus-Seidel benutzen\n");

    calculate_Gaus(N,iter,Matrix_In,k,&options);
    MPI_Barrier(MPI_COMM_WORLD);
    DisplayMatrix(/*&arguments,&results,*/&options,rank,size,1,options.interlines,Matrix_In);
    
    MPI_Finalize();
    return EXIT_SUCCESS;
  }
  
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
