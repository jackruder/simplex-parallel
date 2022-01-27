#ifndef MATRIX_H
#define MATRIX_H
#include "matrix.h"
#endif

#ifndef SIMPLEX_H
#define SIMPLEX_H
#include "simplex.h"
#endif /* ifndef SIMPLEX_H */

#ifndef MATRIXOPS_H
#define MATRIXOPS_H
#include "matrixops.h"
#endif /* ifndef MATRIXOPS_H */

#ifndef ALGEBRA_H
#define ALGEBRA_H
#include "algebra.h"
#endif /* ifndef ALGEBRA_H */

#ifndef STDIO_H
#define STDIO_H
#include <stdio.h>
#endif

#ifndef TIME_H
#define TIME_H
#include <time.h>
#endif

#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif /* ifndef STDLIB_H */

#ifndef STRING_H
#define STRING_H
#include <string.h>
#endif /* ifndef STRING_H */

#ifndef SYSTIME_H
#define SYSTIME_H
#include <sys/time.h>
#endif /* ifndef SYSTIME_H */
int main();
Matrix *inpMatrix();
Matrix *fromfile(char *str);
int main() {
  time_t seconds;
  suseconds_t microseconds;
  struct timeval start, stop;

  printf("A:\n");
  Matrix *A = fromfile("As.txt");
  print(A);
  printf("\nb:\n");
  Matrix *b = fromfile("bs.txt");
  print(b);
  printf("\nct:\n");
  Matrix *ct = fromfile("cts.txt");
  print(ct);

  printf("\nx:\n");
  Matrix *x = fromfile("xs.txt");
  print(x);

  gettimeofday(&start, NULL);
  simplex(A, b, ct, x, 2);
  gettimeofday(&stop, NULL);
  seconds = stop.tv_sec - start.tv_sec;
  microseconds = seconds * 1000000 + stop.tv_usec - start.tv_usec;
  printf("With %d threads: %ld microseconds.\n", 2, microseconds);

  gettimeofday(&start, NULL);
  simplex(A, b, ct, x, 0);
  gettimeofday(&stop, NULL);
  seconds = stop.tv_sec - start.tv_sec;
  microseconds = seconds * 1000000 + stop.tv_usec - start.tv_usec;
  printf("With no threads: %ld microseconds.\n", microseconds);

  return 0;
}
Matrix *fromfile(char *str) {
  FILE *f = fopen(str, "r");
  char *line = NULL;
  size_t len = 0;
  getline(&line, &len, f);
  int m = atoi(line);
  getline(&line, &len, f);
  int n = atoi(line);
  Matrix *r = matrix_empty(m, n);
  for (int i = 0; i < m; i++) {
    getline(&line, &len, f);
    for (int j = 0; j < n; j++) {
      *get_entry_address(r, i, j) = atof(line + 2 * j);
    }
  }
  return r;
}
Matrix *inpMatrix() {
  int i, j, m, n;

  printf("Enter number of Rows :");
  scanf("%d", &m);
  printf("Enter number of Cols :");
  scanf("%d", &n);

  printf("%d", m);
  printf("%d", n);

  Matrix *o = matrix_empty(m, n);

  printf("\nEnter matrix elements :\n");
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("Enter element [%d,%d] : ", i + 1, j + 1);
      scanf("%lf", get_entry_address(o, i, j));
    }
  }
  return o;
  // return NULL;
}
