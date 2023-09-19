/*****************************************************************
 * This programme defines functions used to generate band
 * matrix for problems in which a solution at some grid point
 * k is coupled to a solution at a number of adjasont points.
 * This has been discussed for both, the implicit method
 * of finding solution to time dependent heat equation and
 * for a stationary solution of the heat equation with given
 * sources. The file also includes a main() function with an
 * example of how to solve a stationary heat equation with a
 * constant source.
 * 
 * Prior to compilation execute following lines on nenneke:
 * module purge
 * module load intel impi imkl
 * Then:
 * Compile:  gcc -o bandu band_utility.c -lm -lmkl -liomp5
 * Run: ./bandu
 * ****************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>

/* Define structure that holds band matrix information */
struct band_mat{
  long ncol;        /* Number of columns in band matrix */
  long nbrows;      /* Number of rows (bands in original matrix) */
  long nbands_up;   /* Number of bands above diagonal */
  long nbands_low;  /* Number of bands below diagonal */
  double *array;    /* Storage for the matrix in banded format */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix */
  double *array_inv;/* Store the inverse if this is generated */
  int *ipiv;        /* Additional inverse information */
};
/* Define a new type band_mat */
typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */ 
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
  bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
  bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
  if (bmat->array==NULL||bmat->array_inv==NULL) {
    return 0;
  }  
  /* Initialise array to zero */
  long i;
  for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    bmat->array[i] = 0.0;
  }
  return 1;
};

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat) {
  free(bmat->array);
  free(bmat->array_inv);
  free(bmat->ipiv);
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column);
}

/* Set an element of a band matrix to a desired value based on the pointer
   to a location in the band matrix, using the row and column indexes
   of the full matrix.           */
double setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
  return val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  /* Copy bmat array into the temporary store */
  int i,bandno;
  for(i=0;i<bmat->ncol;i++) { 
    for (bandno=0;bandno<bmat->nbrows;bandno++) {
      bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
    }
    x[i] = b[i];
  }

  long nrhs = 1;
  long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  return info;
}

int printmat(band_mat *bmat) {
  long i,j;
  for(i=0; i<bmat->ncol;i++) {
    for(j=0; j<bmat->nbrows; j++) {
      printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
    }
  }
  return 0;
}

/* An example of how to use the band matrix routines to solve a PDE:
   The equation solved is related to the steady state solution of the heat 
   diffusion equation.   
*/
int main() {
  band_mat bmat;
  long ncols = 1000;
  /* We have a three-point stencil (domain of numerical dependence) of
     our finite-difference equations:
     1 point to the left  -> nbands_low = 1
     1       to the right -> nbands_up  = 1
  */
  long nbands_low = 1;  
  long nbands_up  = 1;
  init_band_mat(&bmat, nbands_low, nbands_up, ncols);
  double *x = malloc(sizeof(double)*ncols);
  double *b = malloc(sizeof(double)*ncols);

    long i;
  /* Loop over the equation number and set the matrix
     values equal to the coefficients of the grid values. 
     Note that boundaries are treated with special cases */
  for(i=0; i<ncols; i++) {
    if(i>0)       {setv(&bmat,i,i-1,-1.0);};
    setv(               &bmat,i,i,   2.0);
    if(i<ncols-1) {setv(&bmat,i,i+1,-1.0);};
    /* Uniform source term in heat equation */
    b[i] = 1.0;   
  }
  
  /*  Print coefficient matrix for debugging: */ 
  /*  printmat(&bmat);           */

  solve_Ax_eq_b(&bmat, x, b);

  for(i=0; i<ncols; i++) {
    printf("%ld %g %g \n",i,x[i],b[i]);
  }

  finalise_band_mat(&bmat);
  return(0);
}
