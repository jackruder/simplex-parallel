#include "algebra.h"

int lu(Matrix *a) {
  int n = a->m;
  for (int k = 0; k < n - 1; k++) {
    for (int j = k + 1; j < n; j++) {
      *get_entry_address(a, j, k) =
          *get_entry_address(a, j, k) / *get_entry_address(a, k, k); // fix this
      for (int l = k + 1; l < n; l++) {
        *get_entry_address(a, j, l) =
            *get_entry_address(a, j, l) -
            *get_entry_address(a, j, k) * *get_entry_address(a, k, l);
      }
    }
  }
  return 0;
}

int forward(Matrix *l, Matrix *b) {
  for (int j = 0; j < l->n - 1; j++) {
    for (int k = j + 1; k < l->m; k++) {
      *get_entry_address(b, k, 0) =
          *get_entry_address(b, k, 0) -
          *get_entry_address(b, j, 0) * (*get_entry_address(l, k, j));
    }
  }
  return 0;
}

int back(Matrix *u, Matrix *b) {
  for (int j = u->n - 1; j > 0; --j) {
    *get_entry_address(b, j, 0) =
        *get_entry_address(b, j, 0) / *get_entry_address(u, j, j);
    for (int i = 0; i < j; i++) {
      *get_entry_address(b, i, 0) =
          *get_entry_address(b, i, 0) -
          *get_entry_address(b, j, 0) * (*get_entry_address(u, i, j));
    }
  }
  *get_entry_address(b, 0, 0) =
      *get_entry_address(b, 0, 0) / (*get_entry_address(u, 0, 0));
  return 0;
}

int solveLU(Matrix *a, Matrix *b) {
  forward(a, b);
  back(a, b);
  return 0;
}
