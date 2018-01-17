# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff-par
OBJS = partdiff-par.o askparams.o

# Targets ...
all: $(TGTS)

partdiff-par: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-par.o: partdiff-par.c Makefile
	$(CC) -c $(CFLAGS) $*.c

askparams.o: askparams.c Makefile

scorep:
	scorep $(CC) -c $(CFLAGS) askparams.c
	scorep $(CC) -c $(CFLAGS) partdiff-par.c
	scorep $(CC) $(LFLAGS) -o partdiff-par $(OBJS) $(LIBS)

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

pde_jacobi_f1_iter.log:
	command -v sbatch >/dev/null || mpiexec ./partdiff-par 1 2 0 1 2 82 >pde_jacobi_f1_iter.log

pde_jacobi_f2_iter.log:
	command -v sbatch >/dev/null || mpiexec ./partdiff-par 1 2 0 2 2 85 >pde_jacobi_f2_iter.log

pde_gaussseidel_f1_iter.log:
	command -v sbatch >/dev/null || mpiexec ./partdiff-par 1 1 0 1 2 47 >pde_gaussseidel_f1_iter.log

pde_gaussseidel_f2_iter.log:
	command -v sbatch >/dev/null || mpiexec ./partdiff-par 1 1 0 2 2 48 >pde_gaussseidel_f2_iter.log

pde_jacobi_f1_prec.log:
	command -v sbatch >/dev/null || mpiexec ./partdiff-par 1 2 0 1 1 1e-6 >pde_jacobi_f1_prec.log

pde_jacobi_f2_prec.log:
	command -v sbatch >/dev/null || mpiexec ./partdiff-par 1 2 0 2 1 1e-6 >pde_jacobi_f2_prec.log

# reference output is based on resulting number of iterations
# but this depends on processes used, due to pipeline principle
# so fix it to 24 like in slurm file
# f1 results to 83 iterations
pde_gaussseidel_f1_prec.log:
	command -v sbatch >/dev/null || mpiexec -n 24 ./partdiff-par 1 1 0 1 1 1e-6 >pde_gaussseidel_f1_prec.log

# f2 results to 84 iterations
pde_gaussseidel_f2_prec.log:
	command -v sbatch >/dev/null || mpiexec -n 24 ./partdiff-par 1 1 0 2 1 1e-6 >pde_gaussseidel_f2_prec.log

test: partdiff-par pde_jacobi_f1_iter.log pde_jacobi_f2_iter.log pde_gaussseidel_f1_iter.log pde_gaussseidel_f2_iter.log pde_jacobi_f1_prec.log pde_jacobi_f2_prec.log pde_gaussseidel_f1_prec.log pde_gaussseidel_f2_prec.log
	cat pde_jacobi_f1_iter.log | sed -e /^$$/d | tail -n 13 > out.log
	tail -n 13 referenz/Jacobi.f1 > ref.log
	diff out.log ref.log && echo "SUCCESS (1/8)"
	cat pde_jacobi_f2_iter.log | sed -e /^$$/d | tail -n 13 > out.log
	tail -n 13 referenz/Jacobi.f2 > ref.log
	diff out.log ref.log && echo "SUCCESS (2/8)"
	cat pde_gaussseidel_f1_iter.log | sed -e /^$$/d | tail -n 13 > out.log
	tail -n 13 referenz/GaussSeidel.f1 > ref.log
	diff out.log ref.log && echo "SUCCESS (3/8)"
	cat pde_gaussseidel_f2_iter.log | sed -e /^$$/d | tail -n 13 > out.log
	tail -n 13 referenz/GaussSeidel.f2 > ref.log
	diff out.log ref.log && echo "SUCCESS (4/8)"
	cat pde_jacobi_f1_prec.log | sed -e /^$$/d | tail -n 13 > out.log
	tail -n 13 referenz/Jacobi_prec.f1 > ref.log
	diff out.log ref.log && echo "SUCCESS (5/8)"
	cat pde_jacobi_f2_prec.log | sed -e /^$$/d | tail -n 13 > out.log
	tail -n 13 referenz/Jacobi_prec.f2 > ref.log
	diff out.log ref.log && echo "SUCCESS (6/8)"
	cat pde_gaussseidel_f1_prec.log | sed -e /^$$/d | tail -n 13 > out.log
	tail -n 13 referenz/GaussSeidel_prec.f1 > ref.log
	diff out.log ref.log && echo "SUCCESS (7/8)"
	cat pde_gaussseidel_f2_prec.log | sed -e /^$$/d | tail -n 13 > out.log
	tail -n 13 referenz/GaussSeidel_prec.f2 > ref.log
	diff out.log ref.log && echo "SUCCESS (8/8)"
	rm out.log ref.log

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
	$(RM) *.log
