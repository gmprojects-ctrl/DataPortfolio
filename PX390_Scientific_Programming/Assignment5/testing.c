
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mkl_lapacke.h>

// Adapted from code in previous Assignments
// Reads input parameters, returns an error value 1, returns zero on success
long read_input(long *N_t, long *N_z, double *tf, long *I_m, char *fname) {
    FILE* fptr=fopen(fname,"r");
    // Check whether we have successfully opened the file
    if (fptr==NULL) {
		return 1;
	}
    // Check whether we've read correct number of values
    if (4!=fscanf(fptr,"%ld %ld %lf %ld", N_t, N_z, tf, I_m)) {
        return 1;
    }
    // Close file and return
    fclose(fptr);
    return 0;
}

// Adapted from code in previous Assignments
// Reads coeff parameters, returns an error value (returns zero on success)
long read_coeff(double *Q_11, double *Q_22, double *Q_12, double *S_input, double *R_input, char *fname, long n) {
    FILE* fptr=fopen(fname,"r");
    // Check whether we have successfully opened the file
    if (fptr==NULL) {
		return 1;
	}
	// Track the number of coefficients read
	long i = 0;
    // Read each line until EOF
    while (EOF != fscanf(fptr,"%lf %lf %lf %lf %lf", Q_11, Q_22, Q_12, S_input, R_input)) {
        // Increment the pointers to the next element to fill from file
		Q_11++;
        Q_12++;
        Q_22++;
		S_input++;
		R_input++;
		i++;
		// If there are too many lines in the coefficients.txt file
		if (i > n){
			return 2;
		}
    }
	// If there are too few lines in the coefficients.txt file
	if (i < n) {
		return 3;
	}
    // Close file and return number of lines read
    fclose(fptr);
    return 0;
}

// Adapted Assignment 4
// Write output to output.txt
long write_output(long theta, long zeta, double time_f, double *Temp) {
	// Inputs: 	theta 		- number of theta grid points
	//			zeta 		- number of zeta grid points
	//			time_f 		- time step
	//			Temp 		- array of temperature values
    FILE* output_file;
    // Name of output file name, set as a constant array of characters
    const char OUTPUT_FILE[] = "output.txt";
    // Open/Create output file
    if (!(output_file = fopen(OUTPUT_FILE, "w"))) {
        // Print error and return error value 1
        printf("Error creating or editing output file\n");
        return 1;
    }
	double pi = 3.14159265358979323846;
	// Theta interval size
	double theta_step = 2 * pi / theta;
	// Zeta interval size
	double zeta_step = 2 * pi / zeta;
	// Store current zeta_i and theta_i
	double theta_i, zeta_i;
	theta_i = 0;
	// For each theta_i
	for (long i = 0; i < theta; i++) {
		// Reset zeta_i to 0
		zeta_i = 0;
		// For each zeta_i
		for (long j = 0; j < zeta; j++) {
			// Print required output
			fprintf(output_file, "%lf %lf %lf %lf\n", time_f, theta_i, zeta_i, Temp[(i*zeta)+j]);
			// Increment zeta_i to zeta_i+1
			zeta_i += zeta_step;
		}
		// Increment theta_i to theta_i+1
		theta_i += theta_step;
	}
    // Close file and return success value 0
    fclose(output_file);
    return 0;
}

// ******************************************************************************************
// LAPACK matrix code adapted from 'band_utility_2d.c' on the Moodle page by Bogdan Hnat
// Most code is copied from my assignment 4
// ******************************************************************************************

// Band matrix
struct band_mat{
    // Number of columns in band matrix
    long ncol;
    // Number of rows in band matrix
    long nbrows;
    // Number of bands above the diagonal
    long nbands_up;
    // Number of bands below the diagonal
    long nbands_low;
    // Storage for the band matrix in banded format
    double *array;
    // Number of rows in inverse matrix
    long nbrows_inv;
    // Store the inverse matrix if generated
    double *array_inv;
    // Additional inverse information
    int *ipiv;
};

// Define a new type band_mat
typedef struct band_mat band_mat;

int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
    // Initialise a band matrix of a certain size, allocate memory and set the parameters
    bmat->nbrows = nbands_lower + nbands_upper + 1;
    bmat->ncol   = n_columns;
    bmat->nbands_up = nbands_upper;
    bmat->nbands_low= nbands_lower;
    bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
    bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
    bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
    bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
    // Check for memory allocation error
    if (bmat->array==NULL||bmat->array_inv==NULL) {
        return 0;
    }  
    // Initialise array to zero
    long i;
    for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
        bmat->array[i] = 0.0;
    }
    return 1;
};

void finalise_band_mat(band_mat *bmat) {
    // Finalise function to free memory
    free(bmat->array);
    free(bmat->array_inv);
    free(bmat->ipiv);
}

double *getp(band_mat *bmat, long row, long column) {
    // Get a pointer to a location in the band matrix using row and column indexes of full matrix
    int bandno = bmat->nbands_up + row - column;
    if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
        printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
        exit(1);
    }
    return &bmat->array[bmat->nbrows*column + bandno];
}

double getv(band_mat *bmat, long row, long column) {
    // Return the value of of a location in the band matrix using row and column indexes of full matrix
    return *getp(bmat,row,column);
}

void setv(band_mat *bmat, long row, long column, double val) {
    // Set an element of a band matrix to a desired value based on the pointer 
    // to a location in the band matrix, using the row and column indexes
    // of the full matrix.
    *getp(bmat,row,column) = val;
}

int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
    // Solve the equation Ax = b for a matrix A stored in band format and x and b real arrays
    // Copy bmat array into the temporary store
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

long indx(long i, long j, long Nt, long Nz) {
	// Find the line value L of i,j using L=i*N_zeta + j
	// Make sure i is in [0,...,N_theta - 1] 
	long a = (i+Nt)%Nt;
	// Make sure j is in [0,...,N_zeta - 1] 
	long b = (j+Nz)%Nz;
	// Return L
	return (a*Nz) + b;
}

long index_fold(long i, long j, long Nt, long Nz) {
	// Fold matrix of submatrices
	// Make sure k is in [0,...,N_theta - 1] 
	long k = (i+Nt)%Nt;
	// Folding algorithm
	if(k<(Nt+1)/2) {
		k = 2*i;
	} else {
		k = (2 * (Nt - i)) - 1;
	}
	// Fold elements in submatrices
	// Make sure l is in [0,...,N_zeta - 1]
	long l = (j+Nz)%Nz;
	// Folding algorithm
	if(l<(Nz+1)/2) {
		l = 2*j;
	} else {
		l = (2 * (Nz - j)) - 1;
	}
	// Return the matrix index of the new folded index k,l
	return indx(k, l, Nt, Nz);
}

int convergence_check(double *T,  double *T_temp, long n) {
	// Tolerance for convergence
	double tol = 1e-8;
	// For each element in the vector
	for (int i=0; i<n; i++) {
		// If percentage change is greater than chosen tolerance
		if (fabs((T[i]-T_temp[i])/T_temp[i])*100 > tol) {
			// Return 1 to indicate convergence has not been achieved
			return 1;
		}
	}
	// Convergence has been achieved if every element has converged
	return 0;
}

int main(void) {
	// Input parameters
	long N_theta, N_zeta, I_min;
	double t_f;
	
	// Input file name
    char* input = "input.txt"; 

    // Read input and make sure it has been read successfully
    if (read_input(&N_theta, &N_zeta, &t_f, &I_min, input)) {
        // Print error message and end code
        printf("File read error for 'input.txt'\n");
        return 1;
    }
	
	// Number of lines in coefficients file
	long n = N_theta * N_zeta;
	
	// Define coefficient arrays
    double *Q11, *Q22, *Q12, *S, *R;

    // Allocate memory according to size of n
    Q11		= (double *)malloc(n*sizeof(double));
    Q22		= (double *)malloc(n*sizeof(double));
    Q12		= (double *)malloc(n*sizeof(double));
	S		= (double *)malloc(n*sizeof(double));
	R		= (double *)malloc(n*sizeof(double));

    // Ensure none of the above are NULL
    if(Q11==NULL || Q22==NULL || Q12==NULL || S==NULL|| R==NULL) {
        // Print error message and end code
        printf("Allocation error\n");
        exit(1);
    }

    // Coefficients file name
    char* coeff = "coefficients.txt"; 
	
	int coeff_check = read_coeff(Q11, Q22, Q12, S, R, coeff, n);

    // Read coefficients and make sure they have been read successfully
    if (coeff_check == 1) {
        // Print error message and end code
        printf("File read error for 'coefficients.txt'\n");
        return 1;
    } else if (coeff_check == 2) {
        // Print error message and end code
        printf("Error in 'coefficients.txt', too many coefficients\n");
        return 1;
    } else if (coeff_check == 3) {
        // Print error message and end code
        printf("Error in 'coefficients.txt', too few coefficients\n");
        return 1;
    }
	
	// Vectors used in forming the matrix equation
	double *T_k, *T_temp, *y, *T, *curr_val_ij;
	T_k 			= (double *)malloc(n*sizeof(double));
	T_temp 			= (double *)malloc(n*sizeof(double));
	y 				= (double *)malloc(n*sizeof(double));
	T 				= (double *)malloc(n*sizeof(double));	
	curr_val_ij		= (double *)malloc(n*sizeof(double));	
	
	// Ensure none of the above are NULL
    if(T==NULL || T_k==NULL || T_temp==NULL || y==NULL || curr_val_ij==NULL) {
        // Print error message and end code
        printf("Allocation error\n");
        exit(1);
    }
	
	// Define for loop counters
    long i, j, k;
	
	// Initialise to an array of (double) zeros as this is the initial value
    for(j=0; j<n; j++) {
        T_k[j] = 0.0;
		T[j] = 0.0;
		T_temp[j] = 0.0;
		y[j] = 0.0;
		curr_val_ij[j] = 0.0;
    }
	
	// Solve PDE using band matrix
	// Create banded matrix object
	band_mat bmat;
	
	// Set number of bands
	long nbands_low = (N_theta*N_zeta)-1;
	
	// Set number of bands above equal to those below
	long nbands_up = nbands_low;

	// Initialise banded marix
	init_band_mat(&bmat, nbands_low, nbands_up, n);
	
	// Set value of pi
	double pi = 3.14159265358979323846;
	
	// Define a = 1/4*theta step squared, b = 1/4*zeta step squared
	// and c = 1/4*theta step*zeta step
	double a, b, c;
	//a = 1 / (4 * (2 * pi / N_theta) * (2 * pi / N_theta)); - rearranged to get:
	a = (N_theta * N_theta) / (16 * pi * pi);	
	//b = 1 / (4 * (2 * pi / N_theta) * (2 * pi / N_theta)); - rearranged to get:
	b = (N_zeta * N_zeta) / (16 * pi * pi);
	//c = 1 / (4 * (2 * pi / N_theta) * (2 * pi / N_zeta)); - rearranged to get:
	c = n / (16 * pi * pi);

	// Temp store of working matrix value
	double mat_val, curr_val;
	long f_idx;
	
	// Loop over all i,j combinations and set values in our matrix
	for(i=0; i<N_theta; i++) {
		for(j =0; j<N_zeta; j++) {
			// Set 8 values by finding their new indexes excluding T_ij
			// as this relies on time step which is not yet known
			
			// Find the folded row of the ij equation
			f_idx = index_fold(i, j, N_theta, N_zeta);
			
			// T_i-1, j-1 =
			mat_val =  c*(Q12[indx(i-1, j, N_theta, N_zeta)]+Q12[indx(i, j-1, N_theta, N_zeta)]);
			curr_val = getv(&bmat, f_idx, index_fold(i-1, j-1, N_theta, N_zeta));
			setv(&bmat, f_idx, index_fold(i-1, j-1, N_theta, N_zeta), mat_val + curr_val);
			
			// T_i-1, j =
			mat_val = 4*a*(Q11[indx(i-1, j, N_theta, N_zeta)] + Q11[indx(i, j, N_theta, N_zeta)])/2;
			curr_val = getv(&bmat, f_idx, index_fold(i-1, j, N_theta, N_zeta));
			setv(&bmat, f_idx, index_fold(i-1, j, N_theta, N_zeta), mat_val + curr_val);			
			
			// T_i-1, j+1 =
			mat_val = -1*c*(Q12[indx(i-1, j, N_theta, N_zeta)]+Q12[indx(i, j+1, N_theta, N_zeta)]);
			curr_val = getv(&bmat, f_idx, index_fold(i-1,j+1, N_theta, N_zeta));
			setv(&bmat, f_idx, index_fold(i-1,j+1, N_theta, N_zeta), mat_val + curr_val);
			
			// T_i, j-1 =
			mat_val = 4*b*(Q22[indx(i, j, N_theta, N_zeta)] + Q22[indx(i, j-1, N_theta, N_zeta)])/2;
			curr_val = getv(&bmat, f_idx, index_fold(i, j-1, N_theta, N_zeta));
			setv(&bmat, f_idx, index_fold(i, j-1, N_theta, N_zeta), mat_val + curr_val);
			
			// T_i, j is set in the do...while loop as time step is not yet known
			
			// T_i, j+1 =
			mat_val = 4*b*(Q22[indx(i, j, N_theta, N_zeta)] + Q22[indx(i, j+1, N_theta, N_zeta)])/2;
			curr_val = getv(&bmat, f_idx, index_fold(i, j+1, N_theta, N_zeta));
			setv(&bmat, f_idx, index_fold(i, j+1, N_theta, N_zeta), mat_val + curr_val);
			
			// T_i+1, j-1 =
			mat_val = -1*c*(Q12[indx(i+1,j, N_theta, N_zeta)] + Q12[indx(i,j-1, N_theta, N_zeta)]);
			curr_val = getv(&bmat, f_idx, index_fold(i+1, j-1, N_theta, N_zeta));
			setv(&bmat, f_idx, index_fold(i+1, j-1, N_theta, N_zeta), mat_val + curr_val);
			
			// T_i+1, j =
			mat_val = 4*a*(Q11[indx(i, j, N_theta, N_zeta)] + Q11[indx(i+1, j, N_theta, N_zeta)])/2;
			curr_val = getv(&bmat, f_idx, index_fold(i+1, j, N_theta, N_zeta));
			setv(&bmat, f_idx, index_fold(i+1, j, N_theta, N_zeta), mat_val + curr_val);			
			
			// T_i+1, j+1 =
			mat_val = c*(Q12[indx(i+1,j, N_theta, N_zeta)]+Q12[indx(i,j+1, N_theta, N_zeta)]);
			curr_val = getv(&bmat, f_idx, index_fold(i+1, j+1, N_theta, N_zeta));
			setv(&bmat, f_idx, index_fold(i+1, j+1, N_theta, N_zeta), mat_val + curr_val);
		}
	}

	// Current time step value
	double t_step;
	// Store the current L index and folded index
	long ij, f_ij;
	// Used to store part of the matrix value at ij
	double mat_val1, mat_val2;
	
	// Store the number of time intervals in use
	long N_int = I_min-1;
	// Set a maximum number of time intervals
	// +10 in case I_min is set to 0
	long I_max = (2 * I_min) + 10;
	
	// Find the current value of matrix at ij
	for (j=0; j<n; j++) {
		curr_val_ij[j] = getv(&bmat, j, j);
	}
	
	// Repeat method until convergent solution found, or no solution found
	do {
		// Increment number of intervals
		N_int += 1;
		// Set the time step
		t_step = t_f / N_int;
		// Repeat for each time interval k
		for(k=0; k<N_int; k++){
			// Set y = -S_i,j - T_{i,j,k-1}/(change in t)
			// For each value in the vector S and T_{k-1}
			for (i=0; i<N_theta; i++) {
				for(j=0; j<N_zeta; j++) {
					// Find the working index
					ij = indx(i, j, N_theta, N_zeta);
					// Find the index after folding
					f_ij = index_fold(i, j, N_theta, N_zeta);
					// Set T_temp to the previous time step values
					T_temp[ij] = T_k[f_ij];
					// Calculate y = -S-T_{k-1} / time step and fold it
					y[f_ij] = -1 * (S[ij]  + (T_temp[ij]/t_step));
					// Set ij value in banded matrix as we now know the time step
					// T_i, j =
					mat_val1 = -4*a*(((Q11[indx(i+1,j, N_theta, N_zeta)] + Q11[ij])/2) + ((Q11[indx(i-1,j, N_theta, N_zeta)] + Q11[ij])/2));
					mat_val2 = -4*b*(((Q22[indx(i,j+1, N_theta, N_zeta)] + Q22[ij])/2) + ((Q22[indx(i,j-1, N_theta, N_zeta)] + Q22[ij])/2));
					mat_val = mat_val1 + mat_val2 - R[ij]-(1/t_step);
					setv(&bmat, f_ij, f_ij, mat_val);
				}
			}
			
			// Add any previous values in each ij that may have occurred due to the periodic nature of ij
			for (j=0; j<n; j++) {
				curr_val = getv(&bmat, j, j);
				setv(&bmat, j, j, curr_val + curr_val_ij[j]);
			}
			
			// Solve A * T_{k+1} = b
			solve_Ax_eq_b(&bmat, T_k, y);
			
			// Commit 'unfolded' T_k to the next part of T
			for (i=0; i<N_theta; i++) {
				for(j=0; j<N_zeta; j++) {
					// 'Unfold' the result to save in T
					T[indx(i, j, N_theta, N_zeta)] = T_k[index_fold(i, j, N_theta, N_zeta)];
				}
			}			
		}	
	} while(convergence_check(T,  T_temp, n) == 1 && N_int < I_max);
	
	// Write output
    if (write_output(N_theta, N_zeta, t_f, T)) {
        // Print error message and end code
        printf("File write error\n");
        return 1;
    }
	
	// Free memory
	// Free vectors from input
    free(Q11);
    free(Q22);
    free(Q12);
	free(S);
	free(R);
	// Free vector that stored temp at final time
	free(T);
	// Free vectors used in matrix equation	
	free(T_k);
	free(T_temp);
	free(y);
	free(curr_val_ij);
	// Free banded matrix
	finalise_band_mat(&bmat);
	// Finish code
	return 0;
}
