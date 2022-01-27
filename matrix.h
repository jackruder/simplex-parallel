typedef struct Matrix {
  double *data; // use a 1d array to store matrix, in column major order
  int m;        // rows
  int n;        // cols
} Matrix;

// generates a random matrix of size m x n, returns a Matrix pointer
Matrix *matrix_rand(int m, int n);

// allocates a new matrix from a source matrix with an array of column indices
// to include
void matrix_from_columns(Matrix *source, Matrix *sink, int *columns, int n);

// gets the column address of a matrix
double *get_column_address(Matrix *matrix, int j);

// handles freeing a dynamically allocated matrix
void delete_matrix(Matrix *matrix);

// gets the address of an element of a matrix
double *get_entry_address(Matrix *matrix, int i, int j);

// returns the column of the matrix in a vector form
Matrix *get_column(Matrix *a, int col);

// initializes an empty matrix given sizes, returns a matrix pointer
Matrix *matrix_empty(int m, int n);

// initializes a matrix from a 2d array given sizes, returns a Matrix pointer
Matrix *matrix_from_array(double *array, int m, int n);

// hard copy a matrix, returning pointer to the new matrix
Matrix *matrix_copy(Matrix *a);

// prints matrix to screen
void print(Matrix *m);
