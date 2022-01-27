#ifndef MATRIX_H
#define MATRIX_H
#include "matrix.h"
#endif

typedef struct tMultArg { // thread argument structure
  Matrix *a;              // first input matrix
  Matrix *b;              // second input matrix
  Matrix *c;              // output matrix
  int i;                  // entry i
  int j;                  // entry j
  int n_calc;             // number of entries for this thread to calculate
} tMultArg;

typedef struct tAddArg { // thread argument structure
  Matrix *a;             // first input matrix
  Matrix *b;             // second input matrix
  Matrix *c;             // output matrix
  int i;                 // entry i
  int j;                 // entry j
  int n_calc;            // number of entries for this thread to calculate
} tAddArg;

// calculates (ab)_ij for matrices a and b for col i and row j, stores answer in
// c_ij
int calc_matrix_mult_entry(Matrix *a, Matrix *b, Matrix *c, int i, int j);

// callback function for threaded matrix multiplication
void *calc_matrix_mult_entry_thread(void *t);
// multiplies Matrices a and b and store the result in c
int matrix_mult(Matrix *a, Matrix *b, Matrix *c);

// multiply Matrices a and b and store the result in c, multithreaded
int matrix_mult_t(Matrix *a, Matrix *b, Matrix *c, int n_threads);
// calculates (ab)_ij for matrices a and b for col i and row j, stores answer in
// c_ij
int calc_matrix_add_entry(Matrix *a, Matrix *b, Matrix *c, int i, int j);

// callback function for threaded matrix addiplication
void *calc_matrix_add_entry_thread(void *t);
// adds Matrices a and b and store the result in c
int matrix_add(Matrix *a, Matrix *b, Matrix *c);

// add Matrices a and b and store the result in c, addithreaded
int matrix_add_t(Matrix *a, Matrix *b, Matrix *c, int n_threads);

// multiply a by scalar s
void matrix_scalar(Matrix *a, double s);

// returns l2 norm of a mx1 matrix (i.e. vector)
double norm(Matrix *m);

double dot(Matrix *a, Matrix *b);
