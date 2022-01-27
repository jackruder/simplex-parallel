#include "matrix.h"
#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif

#ifndef STDIO_H
#define STDIO_H
#include <stdio.h>
#endif

#ifndef TIME_H
#define TIME_H
#include <time.h>
#endif /* ifndef TIME_H */

const double MAX_MATRIX_RAND_ENTRY = 99;

double *get_entry_address(Matrix *matrix, int i, int j) {
  return (matrix->data + j * matrix->m + i); // j * n_cols + i
}

double *get_column_address(Matrix *matrix, int j) {
  return (matrix->data + j * matrix->m);
}

void matrix_from_columns(Matrix *source, Matrix *sink, int *columns, int n) {

  for (int j = 0; j < n; j++) {
    double *col = get_column_address(source, columns[j]);
    for (int i = 0; i < source->m; i++) {
      *get_entry_address(sink, i, j) = col[i];
    }
  }
}

Matrix *get_column(Matrix *a, int col) {
  Matrix *matrix = matrix_empty(a->m, 1);

  matrix_from_columns(a, matrix, &col, 1);
  return matrix;
}

Matrix *matrix_from_array(double *array, int m, int n) { // row by row 2d array
  int i, j;
  Matrix *matrix;

  matrix = matrix_empty(m, n); // new matrix
  if (!matrix) {
    // fprintf(stderr,
    //       "matrix: Error: couldn't initialize matrix of size %dx%d\n", m,
    //      n);
    return NULL;
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      *(get_entry_address(matrix, i, j)) = *(array + i * n + j); // copy values
    }
  }

  return matrix;
}
Matrix *matrix_copy(Matrix *a) {
  Matrix *matrix = matrix_empty(a->m, a->n);
  if (!matrix) {
    // perror("matrix: Error: malloc() failed to allocate memory for matrix");
    return NULL;
  }

  for (int j = 0; j < matrix->n; j++) {
    for (int i = 0; i < matrix->m; i++) {
      *get_entry_address(matrix, i, j) = *get_entry_address(a, i, j);
    }
  }
  return matrix;
}

Matrix *matrix_empty(int m, int n) {
  Matrix *matrix =
      (Matrix *)malloc(sizeof(Matrix)); // create pointer to be returned
  if (!matrix) {
    // perror("matrix: Error: malloc() failed to allocate memory for matrix");
    return NULL;
  }

  matrix->data = (double *)malloc(sizeof(double) * m *
                                  n); // allocate enough space for a mxn matrix
  if (!(matrix->data)) {
    perror("matrix: Error: malloc() failed to allocate memory for matrix"
           "data");
    free(matrix);
    return NULL;
  }

  matrix->m = m;
  matrix->n = n;
  return matrix;
}

Matrix *matrix_rand(int m, int n) {
  int i, j;
  Matrix *matrix;
  matrix = matrix_empty(m, n); // new matrix
  if (!matrix) {
    // fprintf(stderr,
    //       "matrix: Error: couldn't initialize matrix of size %dx%d\n", m,
    //      n);
    return NULL;
  }
  for (i = 0; i < m; i++) {   // for each row
    for (j = 0; j < n; j++) { // for each entry in row
      *(get_entry_address(matrix, i, j)) =
          MAX_MATRIX_RAND_ENTRY / RAND_MAX *
          rand(); // fill matrix with random float from 0 to RAND_MAX
    }
  }
  return matrix;
}

void delete_matrix(Matrix *matrix) {
  free(matrix->data);
  free(matrix);
}

void print(Matrix *m) {
  for (int i = 0; i < m->m; i++) {
    for (int j = 0; j < m->n; j++) {
      printf("%lf\t", *get_entry_address(m, i, j));
    }
    printf("\n");
  }
}
