#include "simplex.h"

#ifndef ALGEBRA_H
#define ALGEBRA_H
#include "algebra.h"
#endif

#ifndef MATRIXOPS_H
#define MATRIXOPS_H
#include "matrixops.h"
#endif /* ifndef MATRIXOPS_H */

#ifndef FLOAT_H
#define FLOAT_H
#include "float.h"
#endif /* ifndef FLOAT_H */

#ifndef STDLIB_H
#define STDLIB_H
#include "stdlib.h"
#endif /* ifndef STDLIB_H */

#ifndef STDIO_H
#define STDIO_H
#include "stdio.h"
#endif /* ifndef STDIO_H */

SplxTbl *initializeTbl(Matrix *A, Matrix *b, Matrix *ct, Matrix *x) {
  SplxTbl *tbl = (SplxTbl *)malloc(sizeof(SplxTbl));
  tbl->A = A;
  tbl->b = b;
  tbl->ct = ct;
  tbl->x = x;

  tbl->m = A->m;
  tbl->n = A->n;

  tbl->B = matrix_empty(tbl->m, tbl->m);
  tbl->N = matrix_empty(tbl->m, tbl->n - tbl->m);

  tbl->basic = (int *)malloc(sizeof(int) * tbl->n);
  return tbl;
}

void cleanTbl(SplxTbl *tbl) {
  delete_matrix(tbl->A);
  delete_matrix(tbl->b);
  delete_matrix(tbl->ct);
  delete_matrix(tbl->B);
  delete_matrix(tbl->N);
  free(tbl->x);
  free(tbl->basic);
  free(tbl);
}

int updateBasic(SplxTbl *sx, int threaded) {
  int *b_indices = (int *)malloc(sizeof(int) * sx->m);
  int *n_indices = (int *)malloc(sizeof(int) * (sx->n - sx->m));
  int b_index = 0;
  int n_index = 0;
  for (int j = 0; j < sx->n; j++) {
    if (*get_entry_address(sx->x, j, 0)) {
      sx->basic[j] = 1;
      b_indices[b_index] = j;
      b_index++;
    } else {
      if (n_index == sx->n - sx->m) {
        printf("degenerate");
        return -1; // degenerate
      }

      sx->basic[j] = 0;
      n_indices[n_index] = j;
      n_index++;
    }
  }

  matrix_from_columns(sx->A, sx->B, b_indices, sx->m);
  matrix_from_columns(sx->A, sx->N, n_indices, sx->n - sx->m);
  lu(sx->B);
  free(n_indices);
  free(b_indices);
  return 0;
}

int findDirections_t(SplxTbl *sx, Matrix **v, int n_threads) {
  pthread_mutex_t minmutex = PTHREAD_MUTEX_INITIALIZER;
  double min = DBL_MAX;

  int best = -1;
  pthread_t *threads = NULL;

  threads = (pthread_t *)malloc(sizeof(pthread_t) * n_threads);
  int n_nb = sx->n - sx->m;
  int nbpt = n_nb / n_threads; // number of nonbasic vars to calculate
  int extra = n_nb % n_threads;
  int nonbasic[n_nb];
  dArg *thread_args[n_threads];
  int index = 0;
  for (int j = 0; j < sx->n; j++) { // find the nonbasic indicies
    if (!sx->basic[j]) {
      nonbasic[index] = j;
      index++;
    }
  }

  index = 0;
  for (int x = 0; x < n_threads; x++) {
    int n_entries_calc;
    if (x < extra) {
      n_entries_calc = nbpt + 1;
    } else {
      n_entries_calc = nbpt;
    }

    int stop = index + n_entries_calc;
    thread_args[x] = (dArg *)malloc(sizeof(dArg) * n_entries_calc);

    int t = 0;
    while (index < stop) {
      dArg da;
      da.sx = sx;
      da.v = v;
      da.j = nonbasic[index];
      da.mutex = &minmutex;
      da.best = &best;
      da.min = &min;

      thread_args[x][t] = da;

      index++;
      t++;
    }
    pthread_create(&threads[x], NULL, findDirectionThread, thread_args[x]);
  }
  for (int x = 0; x < n_threads; x++) {
    pthread_join(threads[x], NULL);
    free(thread_args[x]);
  }
  free(threads);
  return best;
}

void *findDirectionThread(void *v) {
  dArg *da = (dArg *)(v);
  SplxTbl *sx = da->sx;
  int j = da->j;

  Matrix *a_k = get_column(sx->A, j);
  matrix_scalar(a_k, -1);
  solveLU(sx->B, a_k); // compute d^k_B = -B^-1 * A^k

  Matrix *d_k = matrix_empty(sx->n, 1); // initialize d^k_B
  int index = 0;
  for (int k = 0; k < sx->n; k++) {
    if (sx->basic[k]) {
      *get_entry_address(d_k, k, 0) =
          *get_entry_address(a_k, index, 0); // fill in basic var slots
      index++;
    } else {
      if (j == k) {
        *get_entry_address(d_k, k, 0) = 1; // move in this direction
      } else {
        *get_entry_address(d_k, k, 0) = 0; // don't move
      }
    }
  }
  free(a_k);
  double val = dot(sx->ct, d_k) / norm(d_k);
  if (val > 0) {                   // this is an improving direction
    pthread_mutex_lock(da->mutex); // enter critical section
    if (val < *da->min) {          // improving and best so fa
      *da->min = val;
      if (!(da->v)) {
        free(*da->v);
      }
      *da->v = d_k;
      *da->best = j;
    }
    pthread_mutex_unlock(da->mutex); // exit critical section
  } else {
    free(d_k);
  }
  return NULL;
}
int findDirections(SplxTbl *sx, Matrix **v) {
  double min = DBL_MAX;
  int best = -1;
  for (int j = 0; j < sx->n; j++) {
    if (!sx->basic[j]) {
      Matrix *a_k = get_column(sx->A, j);
      matrix_scalar(a_k, -1);
      solveLU(sx->B, a_k); // compute d^k_B = -B^-1 * A^k

      Matrix *d_k = matrix_empty(sx->n, 1); // initialize d^k_B
      int index = 0;
      for (int k = 0; k < sx->n; k++) {
        if (sx->basic[k]) {
          *get_entry_address(d_k, k, 0) =
              *get_entry_address(a_k, index, 0); // fill in basic var slots
          index++;
        } else {
          if (j == k) {
            *get_entry_address(d_k, k, 0) = 1; // move in this direction
          } else {
            *get_entry_address(d_k, k, 0) = 0; // don't move
          }
        }
      }
      free(a_k);
      double val = dot(sx->ct, d_k) / norm(d_k);

      if (val > 0) { // this is an improving direction
        if (val < min) {
          min = val;
          if (!(*v)) {
            free(*v);
          }
          *v = d_k;
          best = j;
        } else {
          free(d_k);
        }
      } else {
        free(d_k);
      } // otherwise continue
    }
  }
  return best; // if -1, we never found an improving direction, hence the
               // current BFS is optimal
}

double stepSize(Matrix *d, Matrix *x) {
  double min = DBL_MAX;
  int bounded = 0;
  for (int i = 0; i < d->m; i++) {
    double d_j = *get_entry_address(d, i, 0);
    if (d_j < 0) { // only negative components matter
      double val = -1 * *get_entry_address(x, i, 0) / d_j; // choose min ratio
      if (val < min) {
        min = val;
        bounded = 1;
      }
    }
  }
  if (bounded) {
    return min;
  } else {
    return -1;
  }
}

void updateSolution(SplxTbl *sx, Matrix *d, double lambda) {
  matrix_scalar(d, lambda);
  matrix_add(sx->x, d, sx->x);
}

int simplex(Matrix *A, Matrix *b, Matrix *ct, Matrix *x, int threaded) {
  SplxTbl *sx = initializeTbl(A, b, ct, x);

  updateBasic(sx, 0);
  Matrix **v = (Matrix **)malloc(sizeof(Matrix *));

  while (1) {
    int entering;
    if (threaded) {
      entering = findDirections_t(sx, v, threaded);
    } else {
      entering = findDirections(sx, v);
    }
    if (entering < 0) { // optimal
      printTbl(sx);
      if (!(*v)) {
        free(*v);
        free(v);
        cleanTbl(sx);
      }
      return 0;
    }
    printf("\nsimplex direction, choosing %d\n", entering);
    print(*v);

    double lambda = stepSize(*v, sx->x);
    printf("step size of  %f", lambda);
    if (lambda < 0) {
      printf("unbounded\n"); // if returned -1 then unbounded, exit
      free(*v);
      free(v);
      cleanTbl(sx);
    }
    updateSolution(sx, *v, lambda);
    if (updateBasic(sx, 0)) { // update basis
      printf("degenerate\n"); // if returned -1 then degeneracy, exit
      free(*v);
      free(v);
      cleanTbl(sx);

      return -1;
    }
    if (!(*v)) {
      free(*v);
    }
  }
}

void printTbl(SplxTbl *tbl) {

  printf("\n\nSOLUTION FOUND:\nCoefficient Matrix:\n");
  print(tbl->A);

  printf("\nRHS:\n");
  print(tbl->b);

  printf("\nObjective Function coefficients:\n");
  print(tbl->ct);

  printf("\nOptimal Solution:\n");
  print(tbl->x);
}
