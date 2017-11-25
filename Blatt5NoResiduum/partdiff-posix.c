/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and   **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <pthread.h>
#include <assert.h>

#include "partdiff-seq.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

struct thread_struct
{
  int start_index;  //start index of the loop
  int end_index;    //End index of the loop
  int M;
  double const h;
  struct options const* options;  //the options
  double** Matrix_Out;
  double** Matrix_In;
  double pih;
  double fpisin;
  double* residuum;
  int term_iteration;
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */

/* ************************************************************************ */
/* Threadfunction                                                           */
/* ************************************************************************ */

void* threadFunction(void* thread_struct){
  //printf("Thread gestartet\n");
  struct thread_struct* threadstruct = (struct thread_struct*) thread_struct;
  
  int i,k;
  int N = threadstruct->end_index;
  int M = threadstruct->M;
  struct options const* options = threadstruct->options;
  
  double** Matrix_Out = threadstruct->Matrix_Out;
  double** Matrix_In  = threadstruct->Matrix_In;
  
  
  int j;                                      /* local variables for loops */
  //int m1, m2;                                 /* used as indices for old and new matrices */
  double star;                                /* four times center value minus 4 neigh.b values */
  //double residuum = *(threadstruct->residuum);
  int term_iteration = threadstruct->term_iteration;
  //double residuum;                            /* residuum of current iteration */
  //double maxresiduum = *(threadstruct->residuum);                         /* maximum residuum value of a slave in iteration */

  //int const N = arguments->N;
  //double const h = threadstruct->h;

  double pih = threadstruct->pih;
  double fpisin = threadstruct->fpisin;
  
  //i-1 doesn't exist for i == 0 so we ensure i > 0
  if(threadstruct->start_index == 0)
    k = 1;
  else
    k = threadstruct->start_index;
  //printf("Laufe von %d bis %d\n",threadstruct->start_index,N);
  for (i = k; i < N; i++)
  {
    double fpisin_i = 0.0;
    
    if (options->inf_func == FUNC_FPISIN)
    {
      fpisin_i = fpisin * sin(pih * (double)i);
    }
    
    /* over all columns */
    for (j = 1; j < M; j++)
    {
      star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);
      
      if (options->inf_func == FUNC_FPISIN)
      {
        star += fpisin_i * sin(pih * (double)j);
      }
      
      if (options->termination == TERM_PREC || term_iteration == 1)
      {
        //residuum = Matrix_In[i][j] - star;
        //residuum = (residuum < 0) ? -residuum : residuum;
        //printf("Residuum in thread: %f\n",residuum);
        //maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
      }
      //printf("Maxresiduum in Thread: %f",maxresiduum);
      //*(threadstruct->residuum) = maxresiduum;
      Matrix_Out[i][j] = star;
    }
  }
  return NULL;
}

/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t g, i, j;                                /*  local variables for loops   */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 0; i <= N; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				Matrix[g][0][i] = 1.0 - (h * i);
				Matrix[g][N][i] = h * i;
			}

			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
  //printf("Bin in Calculate gewesen!\n");
	unsigned int i; //j;                                   /* local variables for loops */
	int m1, m2;                                 /* used as indices for old and new matrices */
	//double star;                                /* four times center value minus 4 neigh.b values */
	//double residuum;                            /* residuum of current iteration */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */
  
	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	struct thread_struct tstruct[options->number];
	pthread_t threads[options->number];
  int segments[(options->number)+1];
  //double* residuum[options->number];
  //printf("before residuum setting\n");
  //*residuum[0] = 5.0;
  //printf("after residuum setting\n");
  //int errorcode;
		
	//tstruct = malloc(sizeof(tstruct)*options->number);
	//threads = malloc(sizeof(pthread_t)*options->number);
	
	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}
  
  //printf("%d ist N\n",N);

  for(i = 0; i < options->number; i++){
    segments[i] = i*(N/options->number);
    //printf("%d Segment Nr: %d\n",segments[i],i);
  }
  segments[(options->number)] = N;
  //printf("%d Segment Nr: %d\n",segments[(int)options->number],(int)options->number);
  
  //exit(1);
  
	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxresiduum = 0;

		for(i = 0; i < options->number; i++){
      //printf("Bin im for mit i = %d\n",i);
		  tstruct[i].start_index = segments[i];
		  tstruct[i].end_index = segments[i+1];
		  tstruct[i].M = N;
		  //tstruct[i].h = h;
		  tstruct[i].options = options;
		  tstruct[i].Matrix_Out = Matrix_Out;
		  tstruct[i].Matrix_In = Matrix_In ;
		  tstruct[i].pih = pih;
		  tstruct[i].fpisin = fpisin;
      //*residuum[i] = maxresiduum;
      //printf("got here!\n");
      //tstruct[i].residuum = residuum[i];
      //printf("*residuum = %f  residuum[i] = %p\n",*residuum[i],(void*) residuum[i]);
      tstruct[i].term_iteration = term_iteration;
      
		  assert(!pthread_create(&threads[i],NULL,threadFunction,&tstruct[i]));
		}
		
		for(i = 0; i < options->number; i++){
      //printf("Will joinen\n");
		  pthread_join(threads[i], NULL);
		}
		/* over all rows 
		for (i = 1; i < N; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			//over all columns
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}
    */
    //for(i = 0; i < options->number; i++){
      //maxresiduum = (*residuum[i] < maxresiduum) ? maxresiduum : *residuum[i];
      //printf("Residuum = %f Maxresiduum = %f\n",*residuum[i],maxresiduum);
    //}
		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxresiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	//free(tstruct);
	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion DisplayMatrix:                               **/
/**                                                                        **/
/** Die Funktion DisplayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf ("\n");
	}

	fflush (stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

  //printf("AskParams\n");
	AskParams(&options, argc, argv);
  
  //printf("initVariables\n");
	initVariables(&arguments, &results, &options);
  
  //printf("allocateMatrices\n");
	allocateMatrices(&arguments);
  
  //printf("initMatricies\n");
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
  
  //printf("calculate\n");
	calculate(&arguments, &results, &options);
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	DisplayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);

	return 0;
}
