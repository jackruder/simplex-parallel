#ifndef MATRIX_H
#define MATRIX_H
#include "matrix.h"
#endif /* ifndef MATRIX_H */

#ifndef PTHREAD_H
#define PTHREAD_H
#include "pthread.h"
#endif

typedef struct SplxTbl {
  int m; // number of constraints
  int n; // number of variables

  int *basic; // 0 or 1 array, 1 indicates var is basic

  Matrix *A;  // Constraints
  Matrix *b;  // RHS
  Matrix *ct; // objective function coefficient vector transpose
  Matrix *x;  // current BFS

  Matrix *B;
  Matrix *N;

} SplxTbl;

// thread argument structure
typedef struct dArg {

  SplxTbl *sx;
  int j;                  // index of nonbasic variable
  pthread_mutex_t *mutex; // in charge of below

  Matrix **v;
  double *min;
  int *best;

} dArg;

// callback funcion for findDirections_t
void *findDirectionThread(void *v);

/* compute simplex directions, returns index of entering variable, and returns
 * the simplex direction in newly allocated v
 */
int findDirections(SplxTbl *sx, Matrix **v);
int findDirections_t(SplxTbl *sx, Matrix **v, int n_threads);

/* updates the array of basic variables according to the current BFS,
 * then recalculates B and N, and then runs an LU decomposition on B.
 */
int updateBasic(SplxTbl *sx, int threaded);

// initializer, copy values and allocate memory. Pointers in arguments are
// reused
SplxTbl *initializeTbl(Matrix *A, Matrix *b, Matrix *c, Matrix *x);
// frees all memory
void cleanTbl(SplxTbl *tbl);

/*
 * Calculates the maximum step size, returning the maximum number. If negative,
 * the LP is unbounded
 */
double stepSize(Matrix *d, Matrix *x);

// updates the soution using the new direction and step size
void updateSolution(SplxTbl *sx, Matrix *d, double lambda);

// run the simplex method on the paramaters, given an initial BFS
int simplex(Matrix *A, Matrix *b, Matrix *ct, Matrix *x, int threaded);

// prints the table information
void printTbl(SplxTbl *tbl);
