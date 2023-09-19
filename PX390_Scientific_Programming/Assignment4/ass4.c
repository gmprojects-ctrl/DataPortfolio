#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<mkl_lapacke.h>

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
// Custom array function 

int scalar_arr(double* arr, int len,double scal){
  for(int i =0;i<len;i++){
    arr[i]*=scal;
  }
  return 0;
}

int parr(double* arr,int len){
  printf("START\n");
  for(int i =0; i<len;i++){
    printf("%lf\n",arr[i]);
  }
  printf("END\n");
  return 0;
}

// Reading/ Writing File functions in C
// Reads input parameters, returns an error value (i.e. returns zero on success).
long read_input(double* L, long* nx, double* freq, char *fname)
{
  FILE *fptr = fopen(fname, "r");
  // Check whether we have successfully opened the file
  if (fptr == NULL)
    return 1;
  // Check whether we've read correct number of values.
  if (3 != fscanf(fptr,"%lf %ld %lf",L,nx,freq))
  {
    return 1;
  }
  fclose(fptr);
  return 0;
}

// Reading the Coefficients file

long read_coeff(double* mu,double* K, double* q, long nx, char* filename) {
    FILE* fptr;
    // Open the file for reading
    fptr = fopen(filename, "r");
    if (fptr == NULL) {
        return 1;
    }
    // Read the values from the file
    for (int i = 0; i < nx; i++) {
        fscanf(fptr, "%lf %lf %lf",&mu[i],&K[i],&q[i]);
    }
    // Close the file
    fclose(fptr);
    return 0;
}



// Write to file

long write_output(double step,double* w, long nx,char *fname){
  double index =0;
  FILE *fptr = fopen(fname,"w");
  // Check whether we have succefully opened the file
  if (fptr == NULL){ 
    return 1;
  }
  // Write the arrays to the file in column form
  for (int i = 0; i < nx; i++) {
      fprintf(fptr, "%lf %lf\n", index, w[i]);
      index+=step;
  }
  // Close the file
  fclose(fptr);
  return 0;
}


int main(){
  double L=0.0;
  double freq=0.0;
  long nx =0;
  char *fname = "input.txt"; 
  if (read_input(&L, &nx,&freq,fname))
  {
      printf("File read error\n");
      return 1;
  }
  double step = L/(nx-1); //Calculate step
  double* mu = (double *)malloc(sizeof(double)*nx);
  double* K  = (double *)malloc(sizeof(double)*nx);
  double* q = (double *)malloc(sizeof(double)*nx);
  if (read_coeff(mu,K,q,nx,"coefficients.txt"))
  {
      printf("File read error\n");
      return 1;
  }
// Scale q by a factor of step^4
//parr(q,nx);
double scal = step * step * step * step;
scalar_arr(q,nx,scal);
// Solution matrix
double* w = (double *)malloc(sizeof(double)*nx);
// Solving the Band Matrix
band_mat bmat;
long ncols = nx;
/* We have a three-point stencil (domain of numerical dependence) of
   our finite-difference equations:
   1 point to the left  -> nbands_low = 1
   1       to the right -> nbands_up  = 1
*/
long nbands_low = 2;  
long nbands_up  = 2;
init_band_mat(&bmat, nbands_low, nbands_up, ncols);
for(int i=0; i<ncols;i++){
  double M = scal * mu[i] * freq * freq;
  if (i==0){
    setv(&bmat,i,i,K[i+1]+K[i+1]+4*K[i]+M);    //Set K[-1] = K[1]
    setv(&bmat,i,i+1,-2*K[i+1]-2*K[i]);
    setv(&bmat,i,i+2,K[i+1]);
  }
  else if(i==1){
    setv(&bmat,i,i-1,-2*K[i-1]-2*K[i]);
    setv(&bmat,i,i,K[i+1]+K[i-1]+4*K[i]+M);
    setv(&bmat,i,i+1,-2*K[i+1]-2*K[i]);
    setv(&bmat,i,i+2,K[i+1]);
  }
  else if(i==ncols-2){
    setv(&bmat,i,i-2,K[i-1]);
    setv(&bmat,i,i-1,-2*K[i-1]-2*K[i]);
    setv(&bmat,i,i,K[i+1]+K[i-1]+4*K[i]+M);
    setv(&bmat,i,i+1,-2*K[i+1]-2*K[i]);
  }
  else if(i==ncols-1){
    setv(&bmat,i,i-2,K[i-1]);
    setv(&bmat,i,i-1,-2*K[i-1]-2*K[i]);
    setv(&bmat,i,i,K[i-1]+K[i]+4*K[i]+M); //Set K[n+1] = K[n]   
  }
  else{
    setv(&bmat,i,i-2,K[i-1]);
    setv(&bmat,i,i-1,-2*K[i-1]-2*K[i]);
    setv(&bmat,i,i,K[i+1]+K[i-1]+4*K[i]+M);
    setv(&bmat,i,i+1,-2*K[i+1]-2*K[i]);
    setv(&bmat,i,i+2,K[i+1]);
  }
}
solve_Ax_eq_b(&bmat, w, q);
// Writing to output file
if(write_output(step,w,nx,"output.txt")){
  printf("Error writing to file");
  return 1;
}



//Freeing up memory
free(mu);
free(K);
free(q);
free(w);
finalise_band_mat(&bmat);
return 0;
}