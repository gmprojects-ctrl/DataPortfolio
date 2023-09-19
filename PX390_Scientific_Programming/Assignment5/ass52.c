#include<stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mkl_lapacke.h>

// Define M_PI if not define already
#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

// Custom array functions
////////////////////////////////////////////////////////////////////////////////////////////////////
// Scale an array
int scalar_arr(double *arr, int len, double scal)
{
  for (int i = 0; i < len; i++)
  {
    arr[i] *= scal;
  }
  return 0;
}
// Set array to a value
int set_arr(double *arr, int len, double value)
{
  for (int i = 0; i < len; i++)
  {
    arr[i] = value;
  }
  return 0;
}
// Copy the values of arr2 into arr1
int arr_copy(double *arr1, double *arr2, int len)
{
  for (int i = 0; i < len; i++)
  {
    arr1[i] = arr2[i];
  }
  return 0;
}
// Add the values of arr_2 to arr_1
int arr_add(double *arr1, double *arr2, int len)
{
  for (int i = 0; i < len; i++)
  {
    arr1[i] += arr2[i];
  }
  return 0;
}
// Print an array
int parr(double *arr, int len)
{
  printf("START\n");
  for (int i = 0; i < len; i++)
  {
    printf("%lf\n", arr[i]);
  }
  printf("END\n");
  return 0;
}

// Convergence check

int converge(double *arr1, double *arr2, int len)
{
  double tol = 1e-10;
  double val = 0;
  for (int i = 0; i < len; i++)
  {
    val += (arr1[i] - arr2[i]) * (arr1[i] - arr2[i]); // (x_i - y_i)^2 + ....
  }
  val = sqrt(val); //
  if (val < tol)
  {
    return 0;
  }
  return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Banded Matrix Functions
struct band_mat
{
  long ncol;       /* Number of columns in band matrix            */
  long nbrows;     /* Number of rows (bands in original matrix)   */
  long nbands_up;  /* Number of bands above diagonal           */
  long nbands_low; /* Number of bands below diagonal           */
  double *array;   /* Storage for the matrix in banded format  */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;   /* Number of rows of inverse matrix   */
  double *array_inv; /* Store the inverse if this is generated */
  int *ipiv;         /* Additional inverse information         */
};
typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns)
{
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low = nbands_lower;
  bmat->array = (double *)malloc(sizeof(double) * bmat->nbrows * bmat->ncol);
  bmat->nbrows_inv = bmat->nbands_up * 2 + bmat->nbands_low + 1;
  bmat->array_inv = (double *)malloc(sizeof(double) * (bmat->nbrows + bmat->nbands_low) * bmat->ncol);
  bmat->ipiv = (int *)malloc(sizeof(int) * bmat->ncol);
  if (bmat->array == NULL || bmat->array_inv == NULL)
  {
    return 0;
  }
  /* Initialise array to zero */
  long i;
  for (i = 0; i < bmat->nbrows * bmat->ncol; i++)
  {
    bmat->array[i] = 0.0;
  }
  return 1;
};

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat)
{
  free(bmat->array);
  free(bmat->array_inv);
  free(bmat->ipiv);
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column)
{
  int bandno = bmat->nbands_up + row - column;
  if (row < 0 || column < 0 || row >= bmat->ncol || column >= bmat->ncol)
  {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n", row, column, bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows * column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column)
{
  return *getp(bmat, row, column);
}

void setv(band_mat *bmat, long row, long column, double val)
{
  *getp(bmat, row, column) = val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b)
{
  /* Copy bmat array into the temporary store */
  int i, bandno;
  for (i = 0; i < bmat->ncol; i++)
  {
    for (bandno = 0; bandno < bmat->nbrows; bandno++)
    {
      bmat->array_inv[bmat->nbrows_inv * i + (bandno + bmat->nbands_low)] = bmat->array[bmat->nbrows * i + bandno];
    }
    x[i] = b[i];
  }

  long nrhs = 1;
  long ldab = bmat->nbands_low * 2 + bmat->nbands_up + 1;
  int info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  return info;
}

int printmat(band_mat *bmat)
{
  long i, j;
  for (i = 0; i < bmat->ncol; i++)
  {
    for (j = 0; j < bmat->nbrows; j++)
    {
      printf("%ld %ld %g \n", i, j, bmat->array[bmat->nbrows * i + j]);
    }
  }
  return 0;
}

/*Check that a grid point has valid coordinates */
int is_valid(long j, long p, long J, long P)
{
  return (j >= 0) && (j < J) && (p >= 0) && (p < P);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Folding & Indexing Functions

long fold(long i, long N)
{
  if (i < (N + 1) / 2)
  {
    return ((long)2 * i) % N;
  }
  else
  {
    return ((long)(2 * (N - i) - 1)) % N;
  }
}

/* Return the 1D element index corresponding to a particular grid point.
   We can rewrite this function without changing the rest of the code if
   we want to change the grid numbering scheme!
   Output: long integer with the index of the point
   Input:
   long j:  The X grid point index
   long k:  The Y grid point index
   long P:  The number of Y points.
*/
long indx(long j, long p, long P)
{
  return j * P + p;
}

/* Return the 2D point corresponding to a particular 1D grid index */
void gridp(long indx, long P, long *j, long *p)
{
  *j = (long)floor(indx / P); // If  L = j *P + p then j = floor(L/P) assuming p < j
  *p = indx - *j * P;         // p = L - j*P
}

// Return modulus fixed values (Hence no out of bounds errors)
long indxmod(long j, long p, long J, long P)
{
  long j_mod = (j + J) % J;
  long p_mod = (p + P) % P;
  return indx(j_mod, p_mod, P);
}

// Index fold
long indxfold(long i, long j, long N_t, long N_z)
{
  long fold_i = fold(i, N_t);
  long fold_j = fold(j, N_z);
  return indxmod(fold_i, fold_j, N_t, N_z);
}

int pmat(band_mat *bmat)
{
  FILE *fptr = fopen("mat.txt", "w");
  // Check whether we have succefully opened the file
  if (fptr == NULL)
  {
    return 1;
  }
  long i, j;
  for (i = 0; i < bmat->ncol; i++)
  {
    for (j = 0; j < bmat->ncol; j++)
    {
      fprintf(fptr, "%lf\t", getv(bmat, i, j));
    }
    fprintf(fptr, "\n");
  }
  fclose(fptr);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// IO Functions
// Reads input parameters, returns an error value (i.e. returns zero on success).
long read_input(long *N_t, long *N_z, double *f_t, long *I_min, char *fname)
{
  FILE *fptr = fopen(fname, "r");
  // Check whether we have successfully opened the file
  if (fptr == NULL)
    return 1;
  // Check whether we've read correct number of values.
  if (4 != fscanf(fptr, "%ld %ld %lf %ld", N_t, N_z, f_t, I_min))
  {
    return 1;
  }
  fclose(fptr);
  return 0;
}
// Store all the values in a neat struct
struct entry
{
  double Q11;
  double Q22;
  double Q12;
  double R;
};
typedef struct entry entry;

// Reading the Coefficients file

long read_coeff(entry *coeff, double *S, long nx, char *filename)
{
  FILE *fptr;
  // Open the file for reading
  fptr = fopen(filename, "r");
  if (fptr == NULL)
  {
    return 1;
  }
  // Read the values from the file
  for (int i = 0; i < nx; i++)
  {
    if (fscanf(fptr, "%lf %lf %lf %lf %lf", &coeff[i].Q11, &coeff[i].Q22, &coeff[i].Q12, &S[i], &coeff[i].R) != 5)
    {
      return 1;
    };
  }
  // Close the file
  fclose(fptr);
  return 0;
}

// Write to file
long write_output(double time, long N_t, long N_z, double *T, char *fname)
{
  FILE *fptr = fopen(fname, "w");
  // Check whether we have succefully opened the file
  if (fptr == NULL)
  {
    return 1;
  }
  // Write the arrays to the file in column form
  for (int i = 0; i < N_t; i++)
  {
    for (int j = 0; j < N_z; j++)
    {
      double theta_pi = 2 * M_PI * i / N_t;
      double zeta_pi = 2 * M_PI * j / N_z;
      double cTemp = T[indxmod(i, j, N_t, N_z)];
      if (fprintf(fptr, "%lf %lf %lf %lf\n", time, theta_pi, zeta_pi, cTemp) < 0)
      {
        return 1;
      }
    }
  }
  // Close the file
  fclose(fptr);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Main functions

int main()
{
  // Initialize all the necessary variables:
  long N_t = 0;
  long N_z = 0;
  double f_t = 0;
  long I_min = 0;

  // Read the input.txt
  if (read_input(&N_t, &N_z, &f_t, &I_min, "input.txt"))
  {
    printf("File read error\n");
    return 1;
  }

  // Calculating the number of columns
  long ncols = N_t * N_z;

  // Initialise S
  double *S = (double *)(malloc(sizeof(double) * ncols));

  // Initialise the Coeff
  entry *coeff = (struct entry *)(malloc(sizeof(struct entry) * ncols));
  // band_mat bmat;

  // Initialise T
  double *T = (double *)malloc(sizeof(double) * ncols);

  // Initialise T_prev
  double *T_prev = (double *)malloc(sizeof(double) * ncols);

  // Initialise T_t
  double *T_t = (double *)malloc(sizeof(double) * ncols);

  // Initialise diagonalarrr

  double *darr = (double *)malloc(sizeof(double) * ncols);

  // Initlaise fold_arr
  double *foldarr = (double *)malloc(sizeof(double) * ncols);
  
  if (S == NULL || T == NULL || T_prev == NULL || coeff == NULL || foldarr == NULL || darr == NULL || T_t == NULL)
  {
    printf("Allocation Error");
    return 1;
  }

  // Read the coefficients.txt
  if (read_coeff(coeff, S, ncols, "coefficients.txt"))
  {
    printf("File read error");
    return 1;
  }

  // Setting the previous Temperature to zero
  set_arr(T_prev, ncols, 0.0);
  set_arr(T, ncols, 0.0);

  // Calculating the time step
  double time_step = f_t / (I_min - 1);

  // Set maximum time length

  double I_MAX = 2 * I_min - 1;

  // Current time

  double ctime = 0;

  // Time Length counter

  double I_COUNTER = 0;

  // Changes in time
  double dt = 2 * M_PI / N_t;
  double dz = 2 * M_PI / N_z;
  double dtdz = dt * dz;
  double a = 1 / (4 * dt * dt);
  double b = 1 / (4 * dz * dz);
  double c = 1 / (4 * dtdz);

  // Creating the Matrix
  long nbands_low = ncols - 2;
  long nbands_up = nbands_low;
  band_mat bmat;
  init_band_mat(&bmat, nbands_low, nbands_up, ncols);

  // Different Values
  // a T_{i-2}{j}    a * coeff[index(i,j,N_z)].Q11
  // b T_{i-1}{j-1}  c * ( coeff[index(i-1,j,N_z)].Q12 + coeff[index(i,j-1,N_z)].Q12 )
  // c T_{i-1}{j+1}  -c * (coeff[index(i-1,j,N_z)].Q12 + coeff[index(i,j+1,N_z)].Q12)
  // d T_{i}{j-2}  b * (coeff[index(i,j-1,N_z)].Q22)
  // e T_{i}{j}   -a * () - b * () - coeff[index(i,j)].R - 1/time_step
  // f T_{i}{j+2}  b * coeff[index(i,j+1)].Q22
  // g T_{i+1}{j-1}  -c * coeff[index]
  // h T_{i+1}{j+1}
  // i T_{i+2}{j}

  // Values to store the current values of T_ij

  long l_i;
  long l_j;
  long m = 0;
  double current_value = 0;
  double past_value = 0;
  double part_a = 0;
  double part_b = 0;

  // Setting the values of the banded matrix
  for (int i = 0; i < ncols; i++)
  {
    gridp(i, N_z, &l_i, &l_j); // Get the Current values of i and j in the matrix

    long fold_row = indxfold(l_i, l_j, N_t, N_z);
    //printf("%ld %ld\n", l_i, l_j);

    // a T_{i-2}{j}
    m = indxfold(l_i - 2, l_j, N_t, N_z);
    current_value = a * coeff[indxmod(l_i - 1, l_j, N_t, N_z)].Q11;
    past_value = getv(&bmat, fold_row, m);
    setv(&bmat, fold_row, m, current_value + past_value);
    // setv(&bmat,i,m,1);

    // b T_{i-1}{j-1}
    m = indxfold(l_i - 1, l_j - 1, N_t, N_z);
    current_value = c * (coeff[indxmod(l_i - 1, l_j, N_t, N_z)].Q12 + coeff[indxmod(l_i, l_j - 1, N_t, N_z)].Q12);
    past_value = getv(&bmat, fold_row, m);
    setv(&bmat, fold_row, m, current_value + past_value);
    // setv(&bmat,i,m,2);
    // printf("\n%lf",c * ( coeff[indxmod(l_i-1,l_j,N_t,N_z)].Q12 + coeff[indxmod(l_i,l_j-1,N_t,N_z)].Q12 ));
    // printf("\n%lf",c);

    // c T_{i-1}{j+1}
    m = indxfold(l_i - 1, l_j + 1, N_t, N_z);
    current_value = -c * (coeff[indxmod(l_i - 1, l_j, N_t, N_z)].Q12 + coeff[indxmod(l_i, l_j + 1, N_t, N_z)].Q12);
    past_value = getv(&bmat, fold_row, m);
    setv(&bmat, fold_row, m, current_value + past_value);
    // setv(&bmat,i,m,3);

    // d T_{i}{j-2}
    m = indxfold(l_i, l_j - 2, N_t, N_z);
    current_value = b * coeff[indxmod(l_i, l_j - 1, N_t, N_z)].Q22;
    past_value = getv(&bmat, fold_row, m);
    setv(&bmat, fold_row, m, past_value + current_value);
    // setv(&bmat,i,m,4);

    // e T_{i}{j}
    m = indxfold(l_i, l_j, N_t, N_z);
    part_a = -a * (coeff[indxmod(l_i + 1, l_j, N_t, N_z)].Q11 - coeff[indxmod(l_i - 1, l_j, N_t, N_z)].Q11);
    part_b = -b * (coeff[indxmod(l_i, l_j + 1, N_t, N_z)].Q22 - coeff[indxmod(l_i, l_j - 1, N_t, N_z)].Q22);
    current_value = part_a + part_b - coeff[indxmod(l_i, l_j, N_t, N_z)].R - -1 / time_step;
    past_value = getv(&bmat, fold_row, m);
    // printf("%d %ld\n",i,m);
    setv(&bmat, fold_row, m, past_value + current_value);
    // setv(&bmat,i,m,555);

    // f T_{i}{j+2}
    m = indxfold(l_i, l_j + 2, N_t, N_z);
    current_value = b * coeff[indxmod(l_i, l_j + 1, N_t, N_z)].Q22;
    past_value = getv(&bmat, fold_row, m);
    setv(&bmat, fold_row, m, current_value + past_value);
    // setv(&bmat,i,m,6);

    // g T_{i+1}{j-1}
    m = indxfold(l_i + 1, l_j - 1, N_t, N_z);
    current_value = -c * (coeff[indxmod(l_i + 1, l_j, N_t, N_z)].Q12 + coeff[indxmod(l_i, l_j - 1, N_t, N_z)].Q12);
    past_value = getv(&bmat, fold_row, m);
    setv(&bmat, fold_row, m, current_value + past_value);
    // setv(&bmat,i,m,7);

    // h T_{i+1}{j+1}
    m = indxfold(l_i + 1, l_j + 1, N_t, N_z);
    current_value = c * (coeff[indxmod(l_i + 1, l_j, N_t, N_z)].Q12 + coeff[indxmod(l_i, l_j + 1, N_t, N_z)].Q12);
    past_value = getv(&bmat, fold_row, m);
    setv(&bmat, fold_row, m, current_value + past_value);
    // setv(&bmat,i,m,8);

    // i T_{i+2}{j}
    m = indxfold(l_i + 2, l_j, N_t, N_z);
    current_value = a * coeff[indxmod(l_i + 1, l_j, N_t, N_z)].Q11;
    past_value = getv(&bmat, fold_row, m);
    setv(&bmat, fold_row, m, current_value + past_value);
    // setv(&bmat,i,m,9);
  }

  // pmat(&bmat);
  //  Do a column swap of the first and second column
  // swapcolumn(&bmat,0,1);
  //  Do a swap of the last and first column;
  // swapcolumn(&bmat,ncols-2,ncols-1);
  // pmat(&bmat);
  //  Removes the previous version of output.txt

  // Get the current values at the diagonals

  for (int i = 0; i < ncols; i++)
  {
    darr[i] = getv(&bmat, i, i);
  }

  do
  {

    I_COUNTER += 1;
    time_step = f_t / I_COUNTER;
    ctime = 0;

    while (ctime < f_t)
    {
      // Changing the banded matrix  and folding the variables
      for (int i = 0; i < ncols; i++)
      {
        gridp(i, N_z, &l_i, &l_j);

        // e T_{i}{j}
        m = indxfold(l_i, l_j, N_t, N_z);
        part_a = -a * (coeff[indxmod(l_i + 1, l_j, N_t, N_z)].Q11 - coeff[indxmod(l_i - 1, l_j, N_t, N_z)].Q11);
        part_b = -b * (coeff[indxmod(l_i, l_j + 1, N_t, N_z)].Q22 - coeff[indxmod(l_i, l_j - 1, N_t, N_z)].Q22);
        current_value = part_a + part_b - coeff[indxmod(l_i, l_j, N_t, N_z)].R - -1 / time_step;
        past_value = getv(&bmat, m, m);
        // printf("%d %ld\n",i,m);
        setv(&bmat, m, m, past_value + current_value);
        // setv(&bmat,i,m,555);
        T_t[indxmod(l_i, l_j, N_t, N_z)] = T_prev[m];
        foldarr[m] = -1 * (S[indxmod(l_i, l_j, N_t, N_z)] + (T_t[indxmod(l_i, l_j, N_t, N_z)] / time_step));
      }
      // Add for Periodicity
      for (int i = 0; i < ncols; i++)
      {
        current_value = getv(&bmat, i, i);
        setv(&bmat, i, i, current_value + darr[i]);
      }

      // Solve the equation
      solve_Ax_eq_b(&bmat, T_prev, foldarr);

      // Unfold
      for (int i = 0; i < ncols; i++)
      {
        gridp(i, N_z, &l_i, &l_j);
        T[indxmod(l_i, l_j, N_z, N_t)] = T_prev[indxfold(l_i, l_j, N_z, N_t)];
      }

      ctime += time_step;
    }
  } while (converge(T, T_prev, ncols) != 0 && I_COUNTER < I_MAX);

  if (write_output(f_t, N_t, N_z, T, "output.txt"))
  {
    printf("Error in writing\n");
    return 1;
  };
  // Freeing all the data
  free(T);
  free(T_prev);
  free(T_t);
  free(darr);
  free(foldarr);
  free(S);
  free(coeff);
  finalise_band_mat(&bmat);
  return 0;
}