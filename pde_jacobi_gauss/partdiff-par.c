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
#include <mpi.h>

#include "partdiff-par.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	// New for MPI: local matrix isn't quadratic
	uint64_t  R;              /* number of spaces between rows (rows=R+1)     */
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

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->R = arguments->N;
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
	uint64_t const R = arguments->R;

	// New for MPI: local matrix isn't quadratic
	arguments->M = allocateMemory(arguments->num_matrices * (R + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((R + 1) * sizeof(double*));

		// Set pointers to row starts
		for (j = 0; j <= R; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (R + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options, int rank, int size, int from)
{
	uint64_t g, i, j;                                /*  local variables for loops   */

	uint64_t const N = arguments->N;
	uint64_t const R = arguments->R;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= R; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	// When looping over all rows:
	// row i in local memory is row (i + offset) in global matrix
	int offset = from - 1;

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			// Iterate over rows only
			for (i = 0; i <= R; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * (i + offset));
				Matrix[g][i][N] = h * (i + offset);
			}
			// Iterate over columns only
			for (i = 0; i <= N; i++)
			{
				// Do not handle process borders as matrix borders
				if (rank == 0)
					Matrix[g][0][i] = 1.0 - (h * i);
				if (rank == size - 1)
					Matrix[g][R][i] = h * i;
			}

			if (rank == size - 1)
				Matrix[g][R][0] = 0.0;
			if (rank == 0)
				Matrix[g][0][N] = 0.0;
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation - MPI Gauss-Seidel version                */
/* ************************************************************************ */
static
void
calculate_gs (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, int rank, int size, int from)
{
	int i, j;                                   /* local variables for loops */
	int m1;                                     /* used as indices for old and new matrices */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	// rare case of more processes than matrix rows
	// for gauss-seidel we can exit, as no multi-communications are used
	if (rank >= size)
		return;

	int const N = arguments->N;
	int const R = arguments->R;
	int offset = from - 1;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	double iter_maxres = 0;
	int prec_finish = 0;
	int test_finish = 0;

	/* initialize m1 */
	m1 = 0;

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m1];

		maxresiduum = 0;

		// Communication is a little spread, some before calc some in between and
		// some after.
		// From each row first and last item aren't changed and thus not exchanged.
		if (size > 1 && rank > 0)
		{
			MPI_Recv(&Matrix_In[0][1], N-1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (options->termination == TERM_PREC || term_iteration == 1)
			{
				// Recieve maxresiduum for precision and termination
				// from process with rank-1, use recieved max as base for own search
				MPI_Recv(&iter_maxres, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				maxresiduum = maxresiduum > iter_maxres ? maxresiduum : iter_maxres;
				MPI_Recv(&prec_finish, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}

		/* over all rows */
		for (i = 1; i < R; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)(i + offset));
			}

			/* over all columns */
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

			// process with rank-1 needs our first row ASAP, send it when calculated
			if (i == 1 && size > 1 && rank > 0)
			{
				MPI_Send(&Matrix_In[1][1], N-1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
			}
		}
		// Communication after calc
		if (size > 1 && rank < size - 1)
		{
			MPI_Send(&Matrix_In[R-1][1], N-1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
			if (options->termination == TERM_PREC || term_iteration == 1)
			{
				// Send maxresiduum for precision and termination
				MPI_Send(&maxresiduum, 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
				MPI_Send(&prec_finish, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
			}

			MPI_Recv(&Matrix_In[R][1], N-1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			// in case of size 1 do like unparallel version
			if (size == 1)
			{
				if (maxresiduum < options->term_precision)
					term_iteration = 0;
			}
			// Recieve signal from last process if precision was reached.
			// First message expected when pipeline full, so check iteration.
			// To avoid conflict with other messages use different tag 1.
			// Don't stop immediately, do another iteration to start transmission
			// chain of prec_finish to following process, achieved by having
			// "else if (prec_finish)" below.
			else if (rank == 0 && !prec_finish && results->stat_iteration >= (unsigned)size)
			{
				MPI_Recv(&prec_finish, 1, MPI_INT, size-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			// Last process has accumulated maxresiduum and can determine if wanted
			// precision is reached. Signal this to process 0, although he has gone
			// beyond this iteration. Process 0 will start message chain to signal
			// one process after the other to stop after all finished same iteration.
			// Last process must not stop immediate but wait until message arrives
			// from process chain.
			// But do not signal process 0 again as it will have stopped iterations.
			else if (rank == size-1 && !test_finish)
			{
				test_finish = maxresiduum < options->term_precision;
				MPI_Send(&test_finish, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
			}
			else if (prec_finish)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	// Exchange precision for output
	// process 0 only waits a long time for this after nothing else is to do
	if (size > 1)
	{
		if (rank == size-1)
		{
			MPI_Send(&maxresiduum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		else if (rank == 0)
		{
			MPI_Recv(&iter_maxres, 1, MPI_DOUBLE, size-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			maxresiduum = maxresiduum > iter_maxres ? maxresiduum : iter_maxres;

			results->stat_precision = maxresiduum;
		}
	}

	results->m = m1;
}

/* ************************************************************************ */
/* calculate: solves the equation - MPI Jacobi version                      */
/* ************************************************************************ */
static
void
calculate_j (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, int rank, int size, int from)
{
	int i, j;                                   /* local variables for loops */
	int m1, m2;                                 /* used as indices for old and new matrices */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	int const R = arguments->R;
	int offset = from - 1;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	m1 = 0;
	m2 = 1;

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxresiduum = 0;

		// communication alternating pattern
		// from each row first and last item aren't changed and thus not exchanged
		// rank < size may happen in rare case of less rows than processes where we reduced size variable manually
		if (size > 1 && rank < size)
		{
			if (rank % 2)
			{
				if (rank > 0)
					MPI_Recv(&Matrix_In[0][1], N-1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (rank < size - 1)
					MPI_Send(&Matrix_In[R-1][1], N-1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
				if (rank > 0)
					MPI_Send(&Matrix_In[1][1], N-1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
				if (rank < size - 1)
					MPI_Recv(&Matrix_In[R][1], N-1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			else
			{
				if (rank < size - 1)
					MPI_Send(&Matrix_In[R-1][1], N-1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
				if (rank > 0)
					MPI_Recv(&Matrix_In[0][1], N-1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (rank < size - 1)
					MPI_Recv(&Matrix_In[R][1], N-1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (rank > 0)
					MPI_Send(&Matrix_In[1][1], N-1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
			}
		}

		/* over all rows */
		for (i = 1; i < R; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)(i + offset));
			}

			/* over all columns */
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

		if (options->termination == TERM_PREC || term_iteration == 1)
		{
			// Exchange maxresiduum for precision and termination
			double globalmaxres;
			MPI_Allreduce(&maxresiduum, &globalmaxres, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			maxresiduum = globalmaxres;
		}

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

//include special mpi display version
#include "displaymatrix-mpi.c"

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	// Initialize MPI and get rank and size of running MPI processes
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	AskParams(&options, argc, argv, rank);

	initVariables(&arguments, &results, &options);

	// additional check for rare case of less rows than processes
	// adjust size to keep all checks simple
	// idle processes must of cource not do send/recv
	// but must not terminate as they have to participate in
	// broadcasts, reduce etc.
	if ((unsigned)size > arguments.N-1)
		size = arguments.N-1;

	// Domain decomposition
	int num_rows = (arguments.N-1) / size;
	int from;
	int to;
	if ((unsigned)rank < (arguments.N-1) % size)
	{
		num_rows++;
		from = rank * num_rows + 1;
		to = (rank+1) * num_rows;
	}
	else
	{
		from = rank * num_rows + 1 + (arguments.N-1) % size;
		to = (rank+1) * num_rows + (arguments.N-1) % size;
	}

	// Each process stores and manages only a part
	arguments.R = num_rows + 1;

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options, rank, size, from);

	gettimeofday(&start_time, NULL);

	if (options.method == METH_JACOBI)
	{
		calculate_j(&arguments, &results, &options, rank, size, from);
	}
	else
	{
		calculate_gs(&arguments, &results, &options, rank, size, from);
	}

	gettimeofday(&comp_time, NULL);

	if (rank == 0)
		displayStatistics(&arguments, &results, &options);
	if (rank < size)
		DisplayMatrix_mpi(&arguments, &results, &options, rank, size, from, to);

	freeMatrices(&arguments);

	// Cleanup MPI and exit
	MPI_Finalize();
	return 0;
}
