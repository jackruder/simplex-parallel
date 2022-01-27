#ifndef MATRIX_H
#define MATRIX_H
#include "matrix.h"
#endif

/* decomposes square matrix A (1:k, 1:k) where for i = 1:n-1,
 * A(i,i:k) is overwritten by U(i,i:k)
 * A(i+1:n, i) is overwritten by L(i+ 1:n, i)
 */
int lu(Matrix *A);

/* solve triangular Lx=b, where L is lower triangular and has ones on the
 * diagonal. b is overwritten by x
 * */
int forward(Matrix *l, Matrix *b);

/*Solve triangular Ux=b, where U is upper triangular.
 * b is overwritten by x
 */
int back(Matrix *u, Matrix *b);

/* Solve the system Ax = b, using the lu decomposition stored in A.
 * B is
 * overwritten with the solution
 */
int solveLU(Matrix *a, Matrix *b);
