#include "matrixops.h"

#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif /* ifndef STDLIB_H */

#ifndef PTHREAD_H
#define PTHREAD_H
#include <pthread.h>
#endif /* ifndef PTHREAD_H */

#ifndef LIMITS_H
#define LIMITS_H
#include <limits.h>
#endif

#ifndef MATH_H
#define MATH_H
#include <math.h>
#endif /* ifndef MATH_H */

int matrix_mult(Matrix *a, Matrix *b, Matrix *c) {
  int i, j;
  for (j = 0; j < a->n; j++) {   // for each row
    for (i = 0; i < b->m; i++) { // for each entry in row
      calc_matrix_mult_entry(    // calculate and store the result in c
          a, b, c, i, j);
    }
  }
  return 0;
}
int matrix_mult_t(Matrix *a, Matrix *b, Matrix *c, int n_threads) {
  int m, n, ents, ept, extra, x, stop, index = 0;
  pthread_t *threads = NULL;

  threads =
      (pthread_t *)malloc(sizeof(pthread_t) * n_threads); // create thread array
  m = c->m;                         // number of rows in product
  n = c->n;                         // number of cols in product
  ents = m * n;                     // number of entries to calculate
  ept = ents / n_threads;           // number of entries per thread
  extra = ents % n_threads;         // number of outlying threads
  tMultArg *thread_args[n_threads]; // array of args

  int i = 0;
  int j = 0;
  int e = 0;

  for (x = 0; x < n_threads; x++) { // for each thread
    int n_entries_calc;

    if (x < extra) {
      n_entries_calc = ept + 1; // add on extra entries to thread //TODO
    } else {
      n_entries_calc = ept; // otherwise just the number
    }

    thread_args[x] = (tMultArg *)malloc(
        sizeof(tMultArg) * (n_entries_calc)); // create memory for arguments
    stop = e + n_entries_calc;                // stop here
    index = 0;
    while (e < stop) { // we only went to 1131
      tMultArg ta;     // create argument struct
      ta.a = a;
      ta.b = b;
      ta.c = c;
      ta.i = i;
      ta.j = j;
      ta.n_calc = n_entries_calc;
      thread_args[x][index] = ta; // copy struct to memory

      i++;
      if (i == m) { // if we have reached last row
        i = 0;      // wrap
        j++;
      }
      e = i + n * j; // calculate e
      index++;
    }
    pthread_create(&threads[x], NULL, calc_matrix_mult_entry_thread,
                   thread_args[x]);
  }

  for (int x = 0; x < n_threads; x++) {
    pthread_join(threads[x], NULL);
    free(thread_args[x]);
  }
  free(threads);
  return 0;
}

int calc_matrix_mult_entry(Matrix *a, Matrix *b, Matrix *c, int i, int j) {
  int s, k;
  double ent_a, ent_b;

  s = 0;
  for (k = 0; k < a->n; k++) {
    ent_a = *get_entry_address(a, i, k); // a_ik
    ent_b = *get_entry_address(b, k, j); // b_kj
    s += ent_a * ent_b;
  }
  *get_entry_address(c, i, j) = s;
  return 0;
}
void *calc_matrix_mult_entry_thread(void *v) {
  tMultArg *t = (tMultArg *)(v);
  int n_calc = (t[0]).n_calc;
  int i, j, s, x, k;
  double ent_a, ent_b;
  Matrix *a, *b, *c;

  for (x = 0; x < n_calc; x++) {
    a = t[x].a;
    b = t[x].b;
    c = t[x].c;
    i = t[x].i;
    j = t[x].j;
    s = 0;
    for (k = 0; k < a->n; k++) {
      ent_a = *get_entry_address(a, i, k); // a_ik
      ent_b = *get_entry_address(b, k, j); // b_kj
      s += ent_a * ent_b;
    }
    *get_entry_address(c, i, j) = s;
  }
  return NULL;
}

int matrix_add(Matrix *a, Matrix *b, Matrix *c) {
  int i, j;

  for (j = 0; j < a->n; j++) {   // for each column
    for (i = 0; i < b->m; i++) { // for each entry in column
      calc_matrix_add_entry(     // calculate and store the result in c
          a, b, c, i, j);
    }
  }
  return 0;
}
int matrix_add_t(Matrix *a, Matrix *b, Matrix *c, int n_threads) {
  int m, n, ents, ept, extra, x, stop, index = 0;
  pthread_t *threads = NULL;

  threads =
      (pthread_t *)malloc(sizeof(pthread_t) * n_threads); // create thread array
  m = c->m;                        // number of rows in product
  n = c->n;                        // number of cols in product
  ents = m * n;                    // number of entries to calculate
  ept = ents / n_threads;          // number of entries per thread
  extra = ents % n_threads;        // number of outlying threads
  tAddArg *thread_args[n_threads]; // array of args

  for (x = 0; x < n_threads; x++) { // for each thread
    int n_entries_calc;
    int i = 0;
    int j = 0;
    int e = 0;

    if (x < extra) {
      n_entries_calc = ept + 1; // add on extra entries to thread //TODO
    } else {
      n_entries_calc = ept; // otherwise just the number
    }

    thread_args[x] = (tAddArg *)malloc(
        sizeof(tAddArg) * (n_entries_calc)); // create memory for arguments
    stop = e + n_entries_calc;               // stop here
    index = 0;
    while (e < stop) { // we only went to 1131
      tAddArg ta;      // create argument struct
      ta.a = a;
      ta.b = b;
      ta.c = c;
      ta.i = i;
      ta.j = j;
      ta.n_calc = n_entries_calc;
      thread_args[x][index] = ta; // copy struct to memory

      i++;
      if (i == m) { // if we have reached last  row
        i = 0;      // wrap
        j++;
      }
      e = i + j * n; // calculate e
      index++;
    }
    pthread_create(&threads[x], NULL, calc_matrix_add_entry_thread,
                   thread_args[x]);
  }

  for (int x = 0; x < n_threads; x++) {
    pthread_join(threads[x], NULL);
    free(thread_args[x]);
  }
  free(threads);
  return 0;
}

void matrix_scalar(Matrix *a, double s) {
  for (int i = 0; i < a->m; i++) {
    for (int j = 0; j < a->n; j++) {
      *get_entry_address(a, i, j) = *get_entry_address(a, i, j) * s;
    }
  }
}

int calc_matrix_add_entry(Matrix *a, Matrix *b, Matrix *c, int i, int j) {

  double ent_a, ent_b;
  ent_a = *get_entry_address(a, i, j); // a_ik
  ent_b = *get_entry_address(b, i, j); // b_kj
  *get_entry_address(c, i, j) = ent_a + ent_b;
  return 0;
}
void *calc_matrix_add_entry_thread(void *v) {
  tAddArg *t = (tAddArg *)(v);
  int n_calc = (t[0]).n_calc;
  int i, j, x;
  double ent_a, ent_b;
  Matrix *a, *b, *c;

  for (x = 0; x < n_calc; x++) {
    a = t[x].a;
    b = t[x].b;
    c = t[x].c;
    i = t[x].i;
    j = t[x].j;
    ent_a = *get_entry_address(a, i, j); // a_ik
    ent_b = *get_entry_address(b, i, j); // b_kj

    *get_entry_address(c, i, j) = ent_a + ent_b;
  }
  return NULL;
}

double norm(Matrix *v) {
  int s = 0;
  for (int i = 0; i < v->m; i++) {
    s += pow(*get_entry_address(v, i, 0), 2);
  }
  return sqrt(s);
}

double dot(Matrix *a, Matrix *b) {
  double s;
  for (int i = 0; i < a->n; i++) {
    s += *get_entry_address(a, 0, i) * *get_entry_address(b, i, 0);
  }
  return s;
}
