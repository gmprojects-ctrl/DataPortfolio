//
// Insert list of bugs fixed here.
// No inclusion of math.h
// swap_mem func incorrectly defined (no sized and double pointers)
// removed nsteps as it wasn't used 
// inmput.txt -> input.txt
// Grid Spacing Error should be dx = lx /(nx-1)
// Incorect definition of malloc (no casting)
// y|| y_next == NULL should be y==NULL || y_next == NULL
// Need to use swap_memory not y=y_next
// Second forloop  should be nx-1 not nx

#include <stdlib.h>
#include <stdio.h>
#include <math.h> // No math.h included

void swap_mem(double **array1, double **array2);

long read_input(double *Lx_p, long *nx_p, double *c_p, double *tf_p, double *dt_p, long *ndt_p, double *S_p, char *fname);

double source(double S, double x);

int main(void)
{
  // **********
  // Parameters
  // **********

  long nx, ndt; //removed nsteps (unused variable)
  double lx, c, tf, dt, S;

  char *fname = "input.txt"; // Spelling Mistake
  if (read_input(&lx, &nx, &c, &tf, &dt, &ndt, &S, fname))
  {
    printf("File read error\n");
    return 1;
  }

  double dx = lx / (nx-1); // Grid Spacing Error
  if ((c * dt) / (dx * dx) > ((double)1 / 8)) //Casting Error
  {
    printf("Timestep too large\n");
    return 1;
  }

  // ************
  // Grid Storage
  // ************
  double *y, *y_next;

  /* Allocate memory according to size of nx */
  y = (double *)malloc(nx * sizeof(double));      // Missed casting
  y_next = (double *)malloc(nx * sizeof(double)); // Missed casting
  //printf("Orig Pointer Post: %p \t %p \n",y,y_next);

  if ((y == NULL) || (y_next == NULL))
  {
    printf("Allocation error\n");
    exit(1);
  }

  double x;

  // **************
  // initialisation
  // **************
  for (int j = 0; j < nx; j++)
  { // Syntax error
    y[j] = 0.;
  }

  double ctime = 0.;
  // Output at start of simulation.
  for (int j = 0; j < nx; j++)
  {
    x = j * dx - lx * 0.5;
    //printf("%d\n", j);
    printf("%g %g %g \n", ctime, x, y[j]);
  }
  long ntstep = 0;
  // loop over timesteps
  while (ctime < tf)
  {

    // loop over points
    for (int j = 1; j < nx-1 ; j++) //Should be nx-1
    {
      // Centered finite difference evaluation of 2nd deriv
      double lslope = (y[j + 1] + y[j - 1] - 2 * y[j]) / (dx * dx);
      y_next[j] = y[j] + dt * (source(S, x) + lslope * c);
    }

    // Copy next values at timestep to y array.
    // Need to swap memory arrays
    swap_mem(&y,&y_next);
    //printf("After Orig: %p \t %p\n",y,y_next);
    // Increment time.
    ctime += dt;
    ntstep++;
    if (ntstep % ndt == 0)
    { 
      for (int j = 0; j < nx; j++)
      {
        x = j * dx - lx * 0.5;
        printf("%g %g %g \n", ctime, x, y[j]);
      }
    }
  }
  //printf("%p \t %p",y,y_next);
  free(y_next);
  free(y); // y & y_next are memory addresses already, no need to free the address the address
}

// Reads input parameters, returns an error value (i.e. returns zero on success).
long read_input(double *Lx_p, long *nx_p, double *c_p, double *tf_p, double *dt_p, long *ndt_p, double *S_p, char *fname)
{
  FILE *fptr = fopen(fname, "r");
  // Check whether we have successfully opened the file
  if (fptr == NULL)
    return 1;
  // Check whether we've read correct number of values.
  if (7 != fscanf(fptr, "%lf %ld %lf %lf %lf %ld %lf", Lx_p, nx_p, c_p, tf_p, dt_p, ndt_p, S_p))
  {
    return 1;
  }
  fclose(fptr);
  return 0;
}

// Need to define the memory swap function
//void swap_mem(double **array1, double **array2, int size)
//{
//  double *temp = (double *)malloc(size * sizeof(double));
//  for (int i = 0; i < size; i++)
//  {
//    temp[i] = *array1[i];
//  }
//  for (int i = 0; i < size; i++)
//  {
//    *array1[i] = *array2[i];
//  }
//  for (int i = 0; i < size; i++)
//  {
//    *array2[i] = temp[i];
//  }
//  free(temp);
//}

void swap_mem(double** arr1, double** arr2){
  double* temp = *arr1;
  //printf("Memswap before:%p\t%p\t%p\n",arr1,arr2,temp);
  *arr1 = *arr2;
  *arr2 = temp;
  //printf("Memswap after:%p\t%p\t%p\n",arr1,arr2,temp);
  }


double source(double S, double x)
{
  return S * exp(-(pow(abs(x), 2.0)));
}
