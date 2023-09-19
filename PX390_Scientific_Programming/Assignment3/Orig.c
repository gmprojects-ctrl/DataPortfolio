//
//  Insert list of bugs fixed here.
//

#include <stdlib.h>
#include <stdio.h>


void swap_mem(double **array1, double **array2);

long read_input(double *Lx_p, long *nx_p, double *c_p, double *tf_p, double *dt_p, long *ndt_p, double *S_p, char *fname); 

double source(double S,double x);

int main(void) {
  // **********
  // Parameters
  // **********
  
  long nsteps, nx, ndt;
  double lx, c, tf, dt, S;

  char* fname = "imput.txt"; 
  if (read_input(&lx, &nx, &c, &tf, &dt, &ndt, &S, fname)) {
    printf("File read error\n");
    return 1;
  }

  double dx = lx/nx;

  if(c*dt/dx*dx>1/8) {
    printf("Timestep too large\n");
    return 1;
  }
  
  // ************
  // Grid Storage 
  // ************
  double *y, *y_next;

  /* Allocate memory according to size of nx */
  y      = malloc(sizeof(nx));
  y_next = malloc(sizeof(nx));
  if(y||y_next==NULL) {
    printf("Allocation error\n");
    exit(1);
  }
  
  int j;
  double x;

  // **************
  // initialisation 
  // **************
  for(j=0,j<=nx;j++;) {
    y[j] = 0;
  }
  
  double ctime;
  // Output at start of simulation.
  for (j=0; j<nx; j++) {
    x = j*dx-lx*0.5;
    printf("%g %g %g \n",ctime,x,y[j]);
  }

  long ntstep = 0;
  //loop over timesteps 
  while (ctime<tf){

    //loop over points 
    for (j=1; j<nx; j++) {
      // Centered finite difference evaluation of 2nd deriv
      double lslope = (y[j+1] + y[j-1] - 2*y[j])/(dx*dx);
      y_next[j] = y[j] + dt*(source(S,x) + lslope * c);
    }
    
    // Copy next values at timestep to y array.
    y=y_next;
     
    // Increment time.   
    ctime += dt;
    ntstep++;
    
    if (ntstep%ndt==0) {
      for (j=0; j<nx; j++ ) {
	x = j*dx-lx*0.5;
	printf("%g %g %g \n",ctime,x,y[j]);
      }
    }
  }

  free(&y);
  free(&y_next);
}

//Reads input parameters, returns an error value (i.e. returns zero on success).
long read_input(double *Lx_p, long *nx_p, double *c_p, double *tf_p, double *dt_p, long *ndt_p, double *S_p, char *fname) {
  FILE* fptr=fopen(fname,"r");
  //Check whether we have successfully opened the file
  if (fptr==NULL) return 1;
  //Check whether we've read correct number of values.
  if (7!=fscanf(fptr,"%lf %ld %lf %lf %lf %ld %lf", Lx_p, nx_p, c_p, tf_p, dt_p, ndt_p, S_p)) {
      return 1;
  }
  fclose(fptr);
  return 0;
}

double source(double S,double x) {
  return S*exp(-(pow(abs(x),2.0)));
}
